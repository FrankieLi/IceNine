"""AR5 step 10, phase A steps 1 and 4: run the grid search again on the disclosed fallback cases,
full (``--n 0``, the reproduction) or with the level-0 cut (``--n 1000``, the replay).

  RESEED_SAMPLE=<grad|mid> uv run python scripts/autoresearch/ar5/cut_replay.py --n 0|1000 \
      [--voxels 7167,5248] [--tag NAME]

For each failed seed step of arm ``A`` (``ar5/rescue/cases_<sample>.json``, 12 on ``mid``, 11 on
``grad``) it calls the reconstructor's grid search with the draws arm ``A`` gave it (generator of
(voxel, "grid", 0), as ``rescue_levels.py``), with the cut of ``cut_rules`` installed when
``--n`` > 0 (script layer: the discrete stage is wrapped; the library is unchanged). Records per
case: the evaluations (global + local), the level-0 list length and the number run, the end's hit
ratio, cost and truth error (symmetry-reduced, the library's angle), and its angle to the full
search's end stored by step 9 (``levels/<sample>_<voxel>.npz``, ``final_R``). Truth is for scoring
only. Writes ``cut/replay_<sample>_n<N><tag>.json`` (tracked, small). ``held`` is not accepted.
"""

import argparse
import contextlib
import io
import json
import os
import time
from typing import Any, Dict, List

import cut_common as CC  # noqa: E402  (first)

import cut_rules as CR  # noqa: E402
import fit_common as FC  # noqa: E402
import numpy as np  # noqa: E402
import order_draws as OD  # noqa: E402
import reseed_common as RC  # noqa: E402


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, required=True)
    ap.add_argument("--voxels", default="")
    ap.add_argument("--tag", default="")
    a = ap.parse_args()
    s = FC.SAMPLE
    assert s in ("grad", "mid"), "the replay is for the disclosed samples"
    cases = json.loads((CC.RCM.RESCUE / f"cases_{s}.json").read_text())
    if a.voxels:
        want = [int(x) for x in a.voxels.split(",")]
        cases = [c for c in cases if int(c["voxel"]) in want]
        assert len(cases) == len(want)
    env = RC.make_env()
    rec = env.hooks.rec
    T = env.T
    log: List[Dict[str, Any]] = []
    CC.install_cut(rec, a.n, log)
    rows: List[Dict[str, Any]] = []
    t0 = time.time()
    for c in cases:
        v = int(c["voxel"])
        verts = env.hooks.vertices(v)
        ph = int(env.mic.voxels[v].phase)
        log.clear()
        with contextlib.redirect_stdout(io.StringIO()):
            res = rec.reconstruct_voxel(
                voxel_vertices=verts,
                phase_index=ph,
                rng=OD.make_rng(OD.voxel_key(verts), "grid", 0),
            )
        g, loc, tot = rec.last_eval_counts
        R = np.asarray(res.orientation, np.float64)
        hit = float(env.bfs._ratios(rec.evaluate_overlap(R, verts, ph))[0])
        err = float(T.err(R[None], [v])[0])
        npz = CC.RCM.RCACHE / "levels" / f"{s}_{v}.npz"
        mis_full = None
        if npz.exists():
            Rf = np.load(npz)["final_R"]
            mis_full = float(env.mis(R[None], np.asarray(Rf, np.float64)[None])[0])
        l0 = [r for r in log if r["level0"]]
        rows.append(
            dict(
                sample=s,
                voxel=v,
                step=c["step"],
                n=a.n,
                evals=int(tot),
                global_evals=int(g),
                local_evals=int(loc),
                level0_listed=int(l0[0]["n_listed"]) if l0 else None,
                level0_run=int(l0[0]["n_run"]) if l0 else None,
                hit=hit,
                cost=float(res.cost),
                err=err,
                mis_to_full_end=mis_full,
                full_end_hit=float(c["grid_end_hit"]),
                full_end_err=float(c["grid_end_err"]),
                full_end_cost=float(c["grid_end_cost"]),
            )
        )
        print(
            f"{s} voxel {v}: n={a.n} evals={tot} err={err:.4f} hit={hit:.3f} "
            f"{time.time() - t0:.0f}s",
            flush=True,
        )
    out = CC.CUT / f"replay_{s}_n{a.n}{a.tag}.json"
    out.write_text(
        json.dumps(
            dict(
                sample=s,
                top_n=a.n,
                rule=CR.RULE if a.n else "no cut (full search)",
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
