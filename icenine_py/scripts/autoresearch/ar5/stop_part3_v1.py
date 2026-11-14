"""AR5 step 1, part 3d: check V1, a seed step at every flagged voxel, on ``mid``.

  RESEED_SAMPLE=mid uv run python scripts/autoresearch/ar5/stop_part3_v1.py [--limit N]
                   [--shard I --of K]

At each voxel flagged at the picked threshold (hit ratio of the map after the pass under it) run
S8's seed step: index, strict cost on the first 20 candidates, local optimization of the 5 cheapest
one by one (``refine_from_candidates``, ``default_rng(0)``), stop at the first end with hit ratio
0.95 or more (else the lowest final cost of the 5). Keep the end if its hit ratio is higher than the
voxel's after the pass. Evaluations: 20 strict costs plus the refinements' counts (attaching rays,
indexing and the overlap evaluation after a refinement are not counted, as in S8). Truth scores
fixed (wrong after the pass, right with V1) and broken (the reverse); no rule reads it.
Writes ``cache/v1_<sample>[_shard].json``; ``stop_part3_summary.py`` writes the slim result.
"""

import argparse
import json
import time
from typing import Any, Dict, List

import numpy as np
import stop_common as SC
import stop_part3_lib as L
import stop_rules as SR
from stop_part2_run import Rungs


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--shard", type=int, default=0)
    ap.add_argument("--of", type=int, default=1)
    ap.add_argument("--cap-min", type=float, default=0.0, help="stop after this many minutes")
    a = ap.parse_args()
    thr = json.loads((SC.AR5 / "part3_pick.json").read_text())["threshold"]
    M = L.load(SC.SAMPLE)
    flagged = np.where(M["hit"] < thr - SR.TOL)[0]
    if a.limit:
        flagged = flagged[: a.limit]
    flagged = flagged[a.shard :: a.of]
    # visit order: every 16th voxel first, then the next residue ..., so that a run stopped by its
    # cap has covered the layer evenly (a fixed order, no random numbers)
    flagged = flagged[np.concatenate([np.arange(s, len(flagged), 16) for s in range(16)])]
    print(SC.SAMPLE, "flagged in this shard", len(flagged), flush=True)
    rr = Rungs()
    env, rec = rr.env, rr.rec
    rows: List[Dict[str, Any]] = []
    t0 = time.time()
    capped = False
    for i, v in enumerate(flagged):
        if a.cap_min and time.time() - t0 > a.cap_min * 60.0:
            capped = True
            print("CAP reached after", i, "voxels", flush=True)
            break
        v = int(v)
        vertices = env.hooks.vertices(v)
        phase = env.mic.voxels[v].phase
        sc, crystal, params, strict_cf = rec._index_context(phase)
        res, _ = rr.RI.index_voxel(sc, vertices, crystal, params, 20)
        R = res.R
        row: Dict[str, Any] = dict(
            voxel=v,
            hit_after_pass=float(M["hit"][v]),
            wrong_after_pass=bool(M["wrong"][v]),
            n_cand=int(len(R)),
            evals=0,
            kept=False,
        )
        if len(R):
            k = min(20, len(R))
            costs = np.array(
                [
                    float(strict_cf.evaluate(np.asarray(R[j], np.float64), vertices, phase).cost)
                    for j in range(k)
                ]
            )
            r = rr.refine_best(vertices, phase, R[:k], costs, v)
            row.update(
                evals=int(k + r["evals"]),
                v1_hit=r["hit"],
                v1_err=r["err"],
                stop=r["stop"],
                kept=bool(r["hit"] > M["hit"][v] + SR.TOL),
            )
            row["wrong_with_v1"] = (
                bool(r["err"] > SC.WRONG) if row["kept"] else row["wrong_after_pass"]
            )
        else:
            row["wrong_with_v1"] = row["wrong_after_pass"]
        rows.append(row)
        if i % 100 == 0:
            print(i, len(flagged), f"{time.time() - t0:.0f} s", flush=True)
    name = f"v1_{SC.SAMPLE}" + (f"_limit{a.limit}" if a.limit else "")
    name += f"_shard{a.shard}of{a.of}" if a.of > 1 else ""
    (SC.AR5 / "cache" / f"{name}.json").write_text(
        json.dumps(
            dict(
                sample=SC.SAMPLE,
                threshold=thr,
                rows=rows,
                wall_s=time.time() - t0,
                capped=capped,
                n_flagged_in_shard=int(len(flagged)),
                cap_min=a.cap_min,
                timing_label="contended (another agent on the machine); not a timing result",
            ),
            default=float,
        )
        + "\n"
    )
    print("done", len(rows), f"{time.time() - t0:.0f} s", flush=True)


if __name__ == "__main__":
    main()
