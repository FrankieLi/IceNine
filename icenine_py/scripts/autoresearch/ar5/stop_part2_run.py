"""AR5 step 1, part 2: a ladder before the grid search, on the grid-fallback cases of S8.

  RESEED_SAMPLE=grad|mid uv run python scripts/autoresearch/ar5/stop_part2_run.py [--limit N]

A case is a seed step of S8 (``e1_log.json.gz``) with ``stop`` = -1 and no refined hit ratio of
0.80 or more: its voxel went to a full grid search. For each case, in this order:

1. *What failed*: the indexer's candidate list (up to 1000), the rank of the first candidate
   within 2 degrees of the truth, whether it is among the 20 costed and the 5 refined (the 5
   cheapest of the 20 by strict cost, as S8 does), and a replay of S8's seed step whose refinement
   evaluations and hit ratios are compared with the stored log (reproduction check).
2. ``R1``: cost the first 100 candidates, refine the 5 cheapest with the seed step's local
   optimization (``refine_from_candidates``, one candidate, ``default_rng(0)``), stop at the first
   end with hit ratio >= 0.95; the end is that one, else the lowest final cost of the 5.
3. ``R3``: run level 0 of the grid search only (``AdaptiveVoxelReconstructor.reconstruct_voxel``
   with ``max_local_resolution`` set to 0 for the call and its final refinement replaced by a
   stub that keeps the candidates; the library is not changed), refine the 5 best orientations of
   its quick MC as above, same stop. The ladder's last rung is then the stored full grid search;
   the library cannot continue after level 0, so the full search pays level 0 again.
4. The stored grid fallback: evaluations (run.json reseeds, or the voxel's record minus the seed
   step's evaluations), hit ratio (record), error (the map of S8 before the pass).

Truth is used only to score. Evaluations are cost-function calls: strict costs of the costed
candidates, the refinements' counts, and the level-0 counts. Attaching rays, indexing and the one
overlap evaluation after a refinement are not counted (as in S8). Wall times are contended.
"""

import argparse
import json
import time
from typing import Any, Dict, List

import numpy as np
import stop_common as SC
import stop_rules as SR

import reseed_common as RC  # noqa: E402

TARGET = 0.95
N_REFINE = 5
DEG = 2.0


def cases() -> List[Dict[str, Any]]:
    log = RC.read_gz(SC.S8_DIR / "e1_log.json.gz")
    run = json.loads((SC.S8_DIR / "run.json").read_text())
    rec = RC.read_gz(SC.S8_DIR / "records.json.gz")
    out = []
    for e in log:
        if e["stop"] != -1 or (e["fin_hit"] and max(e["fin_hit"]) >= SR.SEED_TEST):
            continue
        v = int(e["voxel"])
        if e["reseed"]:
            rs = [x for x in run["stats"]["reseeds"] if x["voxel"] == v]
            assert len(rs) == 1, (v, len(rs))
            grid_evals = int(rs[0]["full_search_evaluations"])
        else:
            grid_evals = int(rec[str(v)]["n_evals"]) - int(e["evals"])
        out.append(
            dict(
                voxel=v,
                reseed=bool(e["reseed"]),
                stored_n_cand=e["n_cand"],
                stored_fin_evals=e["fin_evals"],
                stored_fin_hit=e["fin_hit"],
                stored_seed_step_evals=int(e["evals"]),
                grid_evals=grid_evals,
                grid_hit_record=float(rec[str(v)]["hit_ratio"]),
            )
        )
    return out


class Rungs:
    def __init__(self) -> None:
        self.env = RC.make_env()
        env = self.env
        self.bfs = env.bfs
        self.rec = env.bfs.reconstructor
        self.before = np.load(SC.PASS_DIR / "before_S8.npy")
        self.err_before = np.load(SC.errors_path(SC.SAMPLE))["S8_before"]
        from icenine import ray_indexing as RI
        from icenine.orientation_search import SearchCandidate

        self.RI, self.SC = RI, SearchCandidate

    def err(self, R: np.ndarray, v: int) -> float:
        return float(self.env.T.err(np.asarray(R)[None], np.array([v]))[0])

    def refine(self, vertices: Any, phase: int, R0: np.ndarray, cost: float) -> Dict[str, Any]:
        r = self.rec.refine_from_candidates(
            [self.SC(orientation=np.asarray(R0, np.float64).copy(), cost=float(cost))],
            vertices,
            phase,
            rng=np.random.default_rng(0),
        )
        ev = int(self.rec.last_eval_counts[1])
        info = self.rec.evaluate_overlap(np.asarray(r.orientation, np.float64), vertices, phase)
        hit = float(self.bfs._ratios(info)[0])
        return dict(R=np.asarray(r.orientation, np.float64), cost=float(r.cost), hit=hit, evals=ev)

    def refine_best(
        self, vertices: Any, phase: int, R: np.ndarray, costs: np.ndarray, v: int
    ) -> Dict[str, Any]:
        """Refine the 5 cheapest (by ``costs``), stop at the first with hit >= TARGET; the end is
        that one, else the lowest final cost of those refined."""
        order = [int(j) for j in np.argsort(costs, kind="stable")[:N_REFINE]]
        ends: List[Dict[str, Any]] = []
        stop = -1
        for pos, j in enumerate(order):
            r = self.refine(vertices, phase, R[j], costs[j])
            r["cand"] = j
            ends.append(r)
            if r["hit"] >= TARGET:
                stop = pos
                break
        best = ends[stop] if stop >= 0 else ends[int(np.argmin([x["cost"] for x in ends]))]
        return dict(
            evals=int(sum(x["evals"] for x in ends)),
            refined=[x["cand"] for x in ends],
            fin_evals=[x["evals"] for x in ends],
            fin_hit=[x["hit"] for x in ends],
            stop=stop,
            hit=best["hit"],
            cost=best["cost"],
            err=self.err(best["R"], v),
        )

    def case(self, c: Dict[str, Any]) -> Dict[str, Any]:
        env, rec = self.env, self.rec
        v = c["voxel"]
        vertices = env.hooks.vertices(v)
        phase = env.mic.voxels[v].phase
        t0 = time.time()
        sc, crystal, params, strict_cf = rec._index_context(phase)

        def strict(R: np.ndarray) -> float:
            return float(strict_cf.evaluate(np.asarray(R, np.float64), vertices, phase).cost)

        res, att = self.RI.index_voxel(sc, vertices, crystal, params, 1000)
        R = res.R
        n = len(R)
        Rt = np.repeat(env.Rtrue[v][None], n, 0)
        ang = env.mis(R, Rt) if n else np.zeros(0)
        near = np.where(ang <= DEG)[0]
        rank = int(near[0]) if len(near) else -1  # 0-based rank in the indexer's list
        # step 1: replay S8's seed step (20 costed, 5 refined) to compare with the stored log
        k20 = min(20, n)
        costs20 = np.array([strict(R[j]) for j in range(k20)])
        order20 = [int(j) for j in np.argsort(costs20, kind="stable")[:N_REFINE]]
        replay = [self.refine(vertices, phase, R[j], costs20[j]) for j in order20]
        replay_ok = [int(x["evals"]) for x in replay] == [int(x) for x in c["stored_fin_evals"]]
        replay_hit_ok = bool(
            np.allclose([x["hit"] for x in replay], c["stored_fin_hit"], atol=1e-9)
        )
        # R1: 100 costed (the first 20 are the same strict costs)
        k100 = min(100, n)
        costs100 = np.concatenate([costs20, [strict(R[j]) for j in range(k20, k100)]])
        r1 = self.refine_best(vertices, phase, R[:k100], costs100, v)
        r1["evals"] += int(k100)
        r1["costed"] = int(k100)
        # R3: level 0 of the grid search only
        r3 = self.level0(vertices, phase, v)
        # the stored grid fallback
        gR = self.before[v]
        g_cost, g_hit, _ = env.score3(gR, v)
        out = dict(
            c,
            n_cand=int(n),
            rank_near_truth=rank,
            near_truth_deg=float(ang[rank]) if rank >= 0 else None,
            n_within_2deg=int(len(near)),
            in_costed20=bool(0 <= rank < 20),
            in_first100=bool(0 <= rank < 100),
            in_refined5=bool(rank in order20),
            replay_evals_match=replay_ok,
            replay_hit_match=replay_hit_ok,
            replay_hit=[x["hit"] for x in replay],
            replay_evals=[int(x["evals"]) for x in replay],
            R1=r1,
            R3=r3,
            grid=dict(
                evals=c["grid_evals"],
                hit_record=c["grid_hit_record"],
                hit_map=float(g_hit),
                cost_map=float(g_cost),
                err=float(self.err_before[v]),
            ),
            wall_s=time.time() - t0,
        )
        out["ladder"] = SR.ladder_case(
            dict(evals=r1["evals"], hit=r1["hit"]),
            dict(evals=r3["evals"], hit=r3["hit"]),
            c["grid_evals"],
        )
        end_err = {"R1": r1["err"], "R3": r3["err"], "grid": out["grid"]["err"]}[
            out["ladder"]["rung"]
        ]
        out["ladder"]["err"] = float(end_err)
        out["ladder"]["right"] = bool(end_err <= SC.WRONG)
        return out

    def level0(self, vertices: Any, phase: int, v: int) -> Dict[str, Any]:
        rec = self.rec
        caught: Dict[str, Any] = {}
        orig = rec.refine_from_candidates
        old_max = rec.params.max_local_resolution

        def stub(cands: Any, *a: Any, **k: Any) -> Any:
            caught["c"] = list(cands)
            return cands[0]

        rec.refine_from_candidates = stub
        rec.params.max_local_resolution = 0
        try:
            rec.reconstruct_voxel(vertices, phase, rng=np.random.default_rng(0))
        finally:
            rec.refine_from_candidates = orig
            rec.params.max_local_resolution = old_max
        level0_evals = int(rec.last_eval_counts[2])
        cands = caught["c"]
        R = np.stack([x.orientation for x in cands])
        costs = np.array([x.cost for x in cands])
        top = np.argsort(costs, kind="stable")[:N_REFINE]
        r = self.refine_best(vertices, phase, R[top], costs[top], v)
        r["level0_evals"] = level0_evals
        r["evals"] += level0_evals
        r["n_level0_candidates"] = int(len(cands))
        return r


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--limit", type=int, default=0)
    a = ap.parse_args()
    cs = cases()
    if a.limit:
        cs = cs[: a.limit]
    print(SC.SAMPLE, "cases:", len(cs), flush=True)
    rr = Rungs()
    rows = []
    t0 = time.time()
    for c in cs:
        r = rr.case(c)
        rows.append(r)
        print(
            f"voxel {r['voxel']} rank {r['rank_near_truth']} n_cand {r['n_cand']} "
            f"replay {r['replay_evals_match']}/{r['replay_hit_match']} "
            f"R1 {r['R1']['evals']} hit {r['R1']['hit']:.3f} err {r['R1']['err']:.2f} | "
            f"R3 {r['R3']['evals']} hit {r['R3']['hit']:.3f} err {r['R3']['err']:.2f} | "
            f"grid {r['grid']['evals']} err {r['grid']['err']:.2f} | ladder {r['ladder']['rung']} "
            f"{r['ladder']['evals']:.0f} ({r['wall_s']:.0f} s)",
            flush=True,
        )
    name = f"part2_{SC.SAMPLE}" + (f"_limit{a.limit}" if a.limit else "")
    (SC.AR5 / "cache" / f"{name}.json").write_text(
        json.dumps(
            dict(
                sample=SC.SAMPLE,
                rows=rows,
                wall_s=time.time() - t0,
                timing_label="contended (another agent on the machine); not a timing result",
            ),
            default=float,
        )
        + "\n"
    )


if __name__ == "__main__":
    main()
