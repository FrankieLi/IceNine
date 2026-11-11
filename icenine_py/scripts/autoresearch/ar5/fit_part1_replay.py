"""AR5 step 2, Part 1: replay of the stored seed steps of S8 with four refinement arms.

  scripts/autoresearch/ar5/fit_bg.sh <big|grad|mid> replay_<sample> fit_part1_replay.py \
      [--shard 0 --of 1] [--limit N]

For every seed step in the log of stored run S8 (voxel, reseed) the voxel is indexed again (the seed
step reads only the voxel's rays and uses ``default_rng(0)``; it does not read the BFS state),
the first 20 candidates are costed with the strict cost, and the best 5 by that cost are refined in
that order with each arm of ``fit_common.ARMS`` and S8's early stop (the first end with hit ratio
>= 0.95 is kept; if none, the lowest cost of up to 5 ends). Arm ``M`` is the stored refinement
and is checked against the log on EVERY step (candidates, costed, evaluations and hit ratio of each
refinement, the stop position); the plan asks for every 10th step.

One row per step: ``voxel, reseed, n_cand, costed, log_ok`` and, per arm, ``evals`` (costing plus
refinements), ``fin_evals``, ``fin_hit``, ``fin_cost``, ``fin_err`` (error of each refined end
against the truth, degrees; analysis only), ``taken`` (position), ``stop`` (position or -1),
``hit`` and ``err`` of the end taken. The hit ratio and the error are one evaluation each after the
refinement and are not counted (as in S8). Writes ``cache/replay_<sample>_<shard>of<n>.json.gz``.
"""

import argparse
import time
from typing import Any, Dict, List

import fit_common as FC
import fit_rules as FR
import numpy as np
import reseed_common as RC


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--shard", type=int, default=0)
    ap.add_argument("--of", type=int, default=1)
    ap.add_argument("--limit", type=int, default=0)
    a = ap.parse_args()
    from icenine import ray_indexing as RI

    log = FC.stored_e1(FC.SAMPLE)
    steps = [(i, e) for i, e in enumerate(log) if i % a.of == a.shard]
    if a.limit:
        steps = steps[: a.limit]
    env = RC.make_env()
    rec = env.hooks.rec
    T = env.T
    orig_lo = rec.local_optimization
    rows: List[Dict[str, Any]] = []
    t0 = time.time()
    n_bad = 0
    for n, (i, e) in enumerate(steps):
        v = int(e["voxel"])
        verts = env.hooks.vertices(v)
        ph = int(env.mic.voxels[v].phase)
        sc, crystal, params, strict_cf = rec._index_context(ph)
        res, _att = RI.index_voxel(sc, verts, crystal, params, RI.SEED_TOP_K)
        row: Dict[str, Any] = dict(
            step=i, voxel=v, reseed=bool(e["reseed"]), n_cand=len(res.R), arms={}
        )
        if len(res.R) == 0:
            row["costed"] = 0
            row["log_ok"] = bool(e["n_cand"] == 0)
            rows.append(row)
            continue
        R = res.R
        k = min(FC.N_COST, len(R))
        costs = np.array(
            [
                float(strict_cf.evaluate(np.asarray(R[j], np.float64), verts, ph).cost)
                for j in range(k)
            ]
        )
        order = [int(j) for j in np.argsort(costs, kind="stable")[: FC.N_REFINE]]
        row["costed"] = k
        for arm in FC.ARMS:
            fe: List[int] = []
            fh: List[float] = []
            fc: List[float] = []
            fr: List[float] = []
            ends: List[np.ndarray] = []
            for j in order:
                Rend, c_end, ev = FC.refine_arm(arm, rec, orig_lo, R[j], float(costs[j]), verts, ph)
                info = rec.evaluate_overlap(Rend, verts, ph)
                h = float(env.bfs._ratios(info)[0])
                fe.append(ev)
                fh.append(h)
                fc.append(c_end)
                fr.append(float(T.err(Rend[None], [v])[0]))
                ends.append(Rend)
                if h >= FR.TAU:
                    break
            sp = FR.stop_position(fh)
            taken, stopped = FR.early_stop(fh, fc)
            row["arms"][arm] = dict(
                evals=int(k + sum(fe)),
                fin_evals=fe,
                fin_hit=fh,
                fin_cost=fc,
                fin_err=fr,
                taken=int(taken),
                stop=-1 if sp is None else int(sp),
                hit=fh[taken],
                err=fr[taken],
            )
        m = row["arms"]["M"]
        row["log_ok"] = bool(
            e["n_cand"] == row["n_cand"]
            and e["costed"] == k
            and e["stop"] == m["stop"]
            and e["fin_evals"] == m["fin_evals"]
            and np.allclose(e["fin_hit"], m["fin_hit"], rtol=0, atol=1e-12)
            and int(e["evals"]) == m["evals"]
        )
        n_bad += int(not row["log_ok"])
        rows.append(row)
        if (n + 1) % 10 == 0 or n + 1 == len(steps):
            print(
                f"{FC.SAMPLE} shard {a.shard}: {n + 1}/{len(steps)} {time.time() - t0:.0f}s "
                f"log mismatches so far {n_bad}",
                flush=True,
            )
    name = f"replay_{FC.SAMPLE}_{a.shard}of{a.of}"
    if a.limit:
        name += f"_limit{a.limit}"
    RC.write_gz(rows, FC.CACHE / f"{name}.json.gz")
    (FC.CACHE / f"{name}.done").write_text(f"{time.time() - t0:.1f} mismatches={n_bad}\n")


if __name__ == "__main__":
    main()
