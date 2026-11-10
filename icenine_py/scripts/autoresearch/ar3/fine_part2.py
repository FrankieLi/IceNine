"""AR3 step 4, part 2: step 3's stop readings again, on a fit that moves (the picked arm).

``part2(D, S, T, arm)`` builds, from the picked arm's fits (first fit, with a per-generation trace,
and the second fit on every fourth voxel), the same ``Data``-like object that step 3's summary
functions read, and calls them unchanged: ``pick_skip`` (skip the fit from what is known before
it), ``block_stop`` (stop inside the fit from its trace) and ``second_fit`` (a second fit on every
fourth voxel, and whether it helps). Their tables are renamed ``ar3c_*`` to ``ar3d_p2_*``. Rules
are picked on A and reported on B, as in step 3.

For the CMA fit a "block" is one generation (7 evaluations; the start's own evaluation is not a
block). The trace comes from ``run_fine.py --arm C --trace-cma``.
"""

from types import SimpleNamespace
from typing import Any, Dict

import gradient_common as G
import gradient_summary as GS
import numpy as np

import doc_tables as DT  # noqa: E402

NAMES = {
    "ar3c_skip": "ar3d_p2_skip",
    "ar3c_skip_level": "ar3d_p2_skip_level",
    "ar3c_converged": "ar3d_p2_converged",
    "ar3c_stop": "ar3d_p2_stop",
    "ar3c_stop_level": "ar3d_p2_stop_level",
    "ar3c_second": "ar3d_p2_second",
    "ar3c_second_wil": "ar3d_p2_second_wil",
    "ar3c_second_rule": "ar3d_p2_second_rule",
}


def part2(D: Any, S: Dict[str, Any], T: Dict[str, str], run: str) -> None:
    from fine_summary import FRUNS  # noqa: E402

    rows1 = G.load_gz(FRUNS / run / "refine.json.gz")
    rows2 = G.load_gz(FRUNS / run / "refine2.json.gz")
    n = D.n
    Pm = np.zeros((n, 3, 3))
    P0 = np.load(G.GRAD / "runs" / "P" / "maps.npz")["P"].astype(np.float64)
    Pm[:] = P0
    painted = np.array([r["voxel"] for r in rows1], dtype=int)
    Pm[painted] = D.starts[painted]
    Rm = Pm.copy()
    for r in rows1:
        Rm[r["voxel"]] = np.array(r["end"])
    R2m = Rm.copy()
    for r in rows2:
        R2m[r["voxel"]] = np.array(r["end"])
    d = SimpleNamespace(
        T=D.T,
        n=n,
        lev=D.lev,
        isA=D.isA,
        Pm=Pm,
        Rm=Rm,
        R2m=R2m,
        eP=D.T.err(Pm),
        eR=D.T.err(Rm),
        eR2=D.T.err(R2m),
        rR=rows1,
        rR2=rows2,
        vR=painted,
    )
    S2: Dict[str, Any] = dict(
        run=run,
        n_first=len(rows1),
        n_second=len(rows2),
        evals_first=int(sum(r["evals"] for r in rows1)),
        evals_second=int(sum(r["evals"] for r in rows2)),
    )
    T2: Dict[str, str] = {}
    GS.pick_skip(d, S2, T2)
    GS.block_stop(d, S2, T2)
    GS.second_fit(d, S2, T2)
    S["part2"] = S2
    for k, v in T2.items():
        T[NAMES[k]] = v
    _ = DT
