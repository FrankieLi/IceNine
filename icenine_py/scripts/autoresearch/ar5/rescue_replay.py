"""AR5 step 9, phase A step 4: replay of rule ``G1`` on the 23 disclosed fallback cases.

  RESEED_SAMPLE=<grad|mid> uv run python scripts/autoresearch/ar5/rescue_replay.py

For each case of ``rescue/cases_<sample>.json`` (written by ``rescue_cost.py``: the failed seed
steps of arm ``A``, with the kept ends of the stored log; their orientations are rounded to 12
decimals there, so the replay is a cost estimate, not a bit-exact run) it applies rule ``G1``
exactly as the wrapper does (``rescue_rules.rescue_choose`` and ``rescue_common.cma_rescue``):
CMA from the kept ends in order of cost, first end with hit ratio >= 0.95. Records per case the
tries, evaluations, the hit ratio, cost and truth error of every try, and of the taken end;
the grid end it would replace (hit ratio, error). Writes ``rescue/replay_<sample>.json``. The gate
(``rescue_replay_gate.py``, plan "Decision rules") reads both files. Truth for scoring only.
"""

import json
import os
import time
from typing import Any, Dict, List

import rescue_common as RCM  # noqa: E402  (first)

import fit_common as FC  # noqa: E402
import numpy as np  # noqa: E402
import reseed_common as RC  # noqa: E402
import rescue_rules as RR  # noqa: E402


def main() -> None:
    s = FC.SAMPLE
    assert s in ("grad", "mid"), "the replay is for the disclosed samples"
    cases = json.loads((RCM.RESCUE / f"cases_{s}.json").read_text())
    env = RC.make_env()
    rec = env.hooks.rec
    T = env.T
    orig_lo = rec.local_optimization
    rows: List[Dict[str, Any]] = []
    t0 = time.time()
    for c in cases:
        v = int(c["voxel"])
        verts = env.hooks.vertices(v)
        ph = int(env.mic.voxels[v].phase)
        kept = c["kept"]

        def one(k: int) -> Dict[str, Any]:
            r = RCM.cma_rescue(
                rec,
                orig_lo,
                env.bfs._ratios,
                np.asarray(kept[k]["R"], np.float64).reshape(3, 3),
                verts,
                ph,
            )
            err = float(T.err(r["R"][None], [v])[0])
            return dict(hit=r["hit"], cost=r["cost"], evals=r["evals"], err=err)

        res = RR.rescue_choose([k["cost"] for k in kept], one)
        t = res["tries"][res["taken"]] if res["rescued"] else None
        rows.append(
            dict(
                sample=s,
                voxel=v,
                step=c["step"],
                n_kept=c["n_kept"],
                rescued=res["rescued"],
                n_tries=len(res["tries"]),
                evals=res["evals"],
                tries=res["tries"],
                taken_hit=None if t is None else t["hit"],
                taken_err=None if t is None else t["err"],
                taken_cost=None if t is None else t["cost"],
                wrong_rescued=bool(t is not None and t["err"] > 1.0),
                kept_best_hit=max(k["hit"] for k in kept),
                grid_end_hit=c["grid_end_hit"],
                grid_end_err=c["grid_end_err"],
                grid_end_cost=c["grid_end_cost"],
                grid_end_reached_095=bool(c["grid_end_hit"] >= 0.95),
            )
        )
        print(
            f"{s} voxel {v}: rescued={res['rescued']} tries={len(res['tries'])} "
            f"evals={res['evals']} {time.time() - t0:.0f}s",
            flush=True,
        )
    (RCM.RESCUE / f"replay_{s}.json").write_text(
        json.dumps(
            dict(
                sample=s,
                rule=RR.RULE,
                wall_s=time.time() - t0,
                wall_label="contended; not a timing result",
                rows=rows,
            ),
            indent=1,
        )
        + "\n"
    )
    os._exit(0)


if __name__ == "__main__":
    main()
