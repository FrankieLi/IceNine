"""AR1 end to end: the indexer's candidates checked with the strict cost, the best one finished with
the usual finisher (``refine_from_candidates``: FindOptimal and VarianceMinimizing, as the last
stage of a full search), against the stored full-search result of the same voxel.

Per voxel: indexing time (from recall.py's saved candidates), strict cost of the first K = 20
candidates, the choice among the first K (K = 1, 5, 20) by strict cost and its distance to the
truth, the finisher on the K = 5 choice (evaluations, wall time, final error), and the stored
error and evaluations.

  nohup uv run python scripts/autoresearch/ar1/e2e.py e0 --tag base > <log> 2>&1 &
  nohup uv run python scripts/autoresearch/ar1/e2e.py pilot --tag base > <log> 2>&1 &
Writes benchmarks/autoresearch/ar1/e2e_<kind>_<tag>.json (rows are appended as it goes).
"""

import argparse
import json
import sys
import time
from pathlib import Path
from typing import Any, Dict, List

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import ar1_common as AC  # noqa: E402

from icenine import ray_indexing as RI  # noqa: E402

K_FINISH = 5


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("kind", choices=["e0", "pilot"])
    ap.add_argument("--tag", default="base")
    ap.add_argument("--n", type=int, default=0)
    ap.add_argument("--part", default="0/1", help="i/n: this process takes voxels i::n")
    a = ap.parse_args()
    kind = "e0" if a.kind == "e0" else "full"
    stored = AC.e0_stored() if a.kind == "e0" else AC.pilot_stored()
    voxels = sorted(stored)
    if a.n:
        voxels = voxels[: a.n]
    i, n = (int(x) for x in a.part.split("/"))
    voxels = voxels[i::n]
    stem = f"{a.kind}_{a.tag}"
    rrows, cand = AC.load_recall(a.kind, a.tag)
    recall = {r["voxel"]: r for r in rrows}
    out_path = AC.OUT / f"e2e_{stem}_part{i}of{n}.json"
    AC.init(kind)
    from icenine.reconstructor import SearchCandidate

    rows: List[Dict[str, Any]] = []
    for v in voxels:
        e = AC.env(kind, v, need_sc=False)
        rec, lf, verts, phase, Rt = e["rec"], e["lf"], e["verts"], e["phase"], e["R_true"]
        sym = RI.proper_rotations(AC.symmetry_of(e))
        R = cand[f"v{v}"]
        if len(R) == 0:
            rows.append(dict(voxel=v, no_candidates=True, stored_err=stored[v]["err"]))
            continue
        t0 = time.perf_counter()
        costs = np.array(
            [float(lf.evaluate(R[j], verts, phase).cost) for j in range(len(R))], dtype=np.float64
        )
        t_strict = time.perf_counter() - t0
        dist = RI.misorientation_deg(R, Rt[None], sym)
        pick = {K: int(np.argmin(costs[: min(K, len(R))])) for K in (1, 5, 20)}
        best = pick[K_FINISH]
        # the finisher on the best three of the first five by strict cost, one run each (arm
        # "finish 1" is the first, "finish 3" keeps the lowest final cost of the three)
        order = np.argsort(costs[: min(K_FINISH, len(R))], kind="stable")[:3]
        fins: List[Dict[str, Any]] = []
        for j in order:
            t1 = time.perf_counter()
            res = rec.refine_from_candidates(
                [SearchCandidate(orientation=R[j].copy(), cost=float(costs[j]))],
                verts,
                phase,
                rng=np.random.default_rng(0),
            )
            fins.append(
                dict(
                    cand=int(j),
                    evals=int(rec.last_eval_counts[1]),
                    seconds=time.perf_counter() - t1,
                    err=float(
                        RI.misorientation_deg(np.asarray(res.orientation)[None], Rt[None], sym)[0]
                    ),
                    cost=float(res.cost),
                )
            )
        t_fin = fins[0]["seconds"]
        n_fin = fins[0]["evals"]
        err_fin = fins[0]["err"]
        rows.append(
            dict(
                voxel=v,
                stored_err=stored[v]["err"],
                stored_evals=stored[v]["evals"],
                n_cand=len(R),
                strict_costs=[float(c) for c in costs],
                dist=[float(d) for d in dist],
                pick=pick,
                err_pick={str(K): float(dist[j]) for K, j in pick.items()},
                finisher_evals=n_fin,
                strict_evals=len(R),
                err_final=err_fin,
                cost_final=fins[0]["cost"],
                finishers=fins,
                index_s=recall[v]["wall_s"],
                strict_s=t_strict,
                finisher_s=t_fin,
            )
        )
        print(
            v,
            f"stored={stored[v]['err']:.3f} pick5={dist[best]:.3f} final={err_fin:.3f}",
            f"fin_evals={n_fin} cost={fins[0]['cost']:.3f} fin_s={t_fin:.1f}",
            f"best3_err={min(fins, key=lambda f: f['cost'])['err']:.3f}",
            flush=True,
        )
        out_path.write_text(json.dumps(dict(rows=rows, label="contended")) + "\n")
    out_path.write_text(json.dumps(dict(rows=rows, label="contended")) + "\n")
    out_path.with_suffix(".done").write_text("done\n")


if __name__ == "__main__":
    main()
