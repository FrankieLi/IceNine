"""AR5 step 3, before the runs: does the script layer's mixed refinement give what Part 1 computed?

  scripts/autoresearch/ar5/mix_bg.sh big check mix_check.py [--arm X200] [--n 20]

Takes 20 seed steps of ``big`` that have candidates (10 where Part 1's cache says the MC optimizer
ran on at least one candidate, 10 where it did not; evenly spaced in step order), runs the seed
step's candidate loop with the script layer's ``FC.refine_arm`` (the running version of the
mixed refinement, ``mix_common``) and S8's early stop, and compares with ``mix_rules.mixed_step``
on the replay cache: seed-step evaluations, the number of candidates on which MC ran, the kept
end's hit ratio (tolerance 1e-12) and its error (1e-9 degree), the stop position. Writes
``benchmarks/autoresearch/ar5/mix/check.json``.  Truth is used for the comparison only.
"""

import argparse
import json

import mix_common as MC  # noqa: E402  (first: sets the cache, installs the arms)

import fit_common as FC  # noqa: E402
import fit_part1_lib as L  # noqa: E402
import fit_rules as FR  # noqa: E402
import mix_rules as MR  # noqa: E402
import numpy as np  # noqa: E402
import reseed_common as RC  # noqa: E402


def pick_steps(rows, arm, n):
    cand = [r for r in rows if r["arms"]]
    mc = [
        r
        for r in cand
        if MR.mixed_step(r["costed"], r["arms"][MC.X_ARMS[arm]], r["arms"]["M"])["mc_ran"] > 0
    ]
    no = [r for r in cand if r not in mc]
    h = n // 2

    def spread(xs, k):
        idx = np.linspace(0, len(xs) - 1, k).round().astype(int)
        return [xs[i] for i in idx]

    return spread(mc, h) + spread(no, n - h)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--arm", default="X200")
    ap.add_argument("--n", type=int, default=20)
    a = ap.parse_args()
    assert FC.SAMPLE == "big"
    from icenine import ray_indexing as RI

    rows = L.load("big")
    chosen = pick_steps(rows, a.arm, a.n)
    env = RC.make_env()
    MC.HOLDER["bfs"] = env.bfs
    rec = env.hooks.rec
    T = env.T
    orig_lo = rec.local_optimization
    out_rows = []
    for row in chosen:
        v = int(row["voxel"])
        verts = env.hooks.vertices(v)
        ph = int(env.mic.voxels[v].phase)
        sc, crystal, params, strict_cf = rec._index_context(ph)
        res, _ = RI.index_voxel(sc, verts, crystal, params, RI.SEED_TOP_K)
        R = res.R
        k = min(FC.N_COST, len(R))
        costs = np.array(
            [
                float(strict_cf.evaluate(np.asarray(R[j], np.float64), verts, ph).cost)
                for j in range(k)
            ]
        )
        order = [int(j) for j in np.argsort(costs, kind="stable")[: FC.N_REFINE]]
        n0 = len(MC.LOG)
        hits, fcost, errs, evs = [], [], [], []
        for j in order:
            Rend, c_end, ev = FC.refine_arm(a.arm, rec, orig_lo, R[j], float(costs[j]), verts, ph)
            h = float(env.bfs._ratios(rec.evaluate_overlap(Rend, verts, ph))[0])
            hits.append(h)
            fcost.append(c_end)
            errs.append(float(T.err(Rend[None], [v])[0]))
            evs.append(ev)
            if h >= FR.TAU:
                break
        taken, stopped = FR.early_stop(hits, fcost)
        live = dict(
            evals=int(k + sum(evs)),
            mc_ran=int(sum(x["mc_ran"] for x in MC.LOG[n0:])),
            hit=hits[taken],
            err=errs[taken],
            stopped=bool(stopped),
            taken=int(taken),
        )
        cache = MR.mixed_step(row["costed"], row["arms"][MC.X_ARMS[a.arm]], row["arms"]["M"])
        ok = bool(
            live["evals"] == cache["evals"]
            and live["mc_ran"] == cache["mc_ran"]
            and live["taken"] == cache["taken"]
            and live["stopped"] == cache["stopped"]
            and abs(live["hit"] - cache["hit"]) <= 1e-12
            and abs(live["err"] - cache["err"]) <= 1e-9
        )
        out_rows.append(
            dict(step=row["step"], voxel=v, live=live, cache={f: cache[f] for f in live}, same=ok)
        )
        print(
            row["step"], v, "same" if ok else "DIFFERENT", live["evals"], cache["evals"], flush=True
        )
    res_ = dict(
        arm=a.arm,
        sample="big",
        n_steps=len(out_rows),
        n_with_mc=sum(1 for r in out_rows if r["cache"]["mc_ran"] > 0),
        n_same=sum(r["same"] for r in out_rows),
        evals_live=sum(r["live"]["evals"] for r in out_rows),
        evals_cache=sum(r["cache"]["evals"] for r in out_rows),
        all_same=all(r["same"] for r in out_rows),
        rows=out_rows,
    )
    (MC.MIX / "check.json").write_text(json.dumps(res_, indent=1) + "\n")
    print(json.dumps({k: v for k, v in res_.items() if k != "rows"}, indent=1))


if __name__ == "__main__":
    main()
