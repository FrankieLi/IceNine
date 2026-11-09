"""Shared loaders of the AR1 scripts: the E0 voxels (per-voxel clean images) and the pilot's
full-search voxels (clean full-sample images), each with its truth, SmoothCost and stored result.

The environments are those of Task F's ``outward_bench`` (E0: ``init_worker("e0")``; pilot:
``init_worker("full")``, the clean full-sample reconstruction with the pilot's config). Data that
is not in git (images, caches, stored runs) is read from the main checkout when the worktree has
no copy (``MAIN``).
"""

import json
import sys
from pathlib import Path
from typing import Any, Dict, List, Tuple

import numpy as np

HERE = Path(__file__).resolve().parent
ICE = HERE.parents[2]
MAIN = ICE.parents[3] / "icenine_py"  # read-only data of the main checkout (worktree parent)
for sub in (
    ICE / "scripts" / "eval_reduction",
    ICE / "scripts" / "phase_d",
    ICE / "scripts" / "common",
):
    sys.path.insert(0, str(sub))

OUT = ICE / "benchmarks" / "autoresearch" / "ar1"
E0_STORED = MAIN / "benchmarks" / "eval_reduction" / "E" / "cache" / "rank" / "e0_old_base"
PILOT_RUN = MAIN / "benchmarks" / "eval_reduction" / "G" / "g2" / "bfs" / "runs" / "keyoff_rerun"

_B: Any = None


def bench() -> Any:
    global _B
    if _B is None:
        import outward_bench as B  # noqa: E402

        _B = B
    return _B


def e0_stored() -> Dict[int, Dict[str, float]]:
    """Stored E0 full-search results (arm old_base): voxel -> err_deg, evals, cost_final."""
    out: Dict[int, Dict[str, float]] = {}
    for f in sorted(E0_STORED.glob("v*.npz")):
        z = np.load(f)
        out[int(z["voxel"])] = dict(
            err=float(z["err_deg"]), evals=float(z["evals"]), cost=float(z["cost_final"])
        )
    return out


def pilot_stored() -> Dict[int, Dict[str, Any]]:
    """Stored pilot (key-off rerun) full searches: full-sample voxel index -> row, err (deg),
    source, evaluations, hit ratio. Error from recon.mic against the truth (cubic symmetry)."""
    import g2_bfs

    e = g2_bfs.errors(PILOT_RUN)
    idx = np.load(ICE / "benchmarks" / "phase_d_pilot" / "region_index.npy")
    rec = json.loads((PILOT_RUN / "records.json").read_text())
    out: Dict[int, Dict[str, Any]] = {}
    for i, (src, nev) in enumerate(zip(e["source"], e["n_evals"])):
        # a full search costs about 330 000 evaluations, a local fit about 700; the 69 full searches
        # are the 62 voxels whose final source is "seed" plus 5 "revisit" and 2 "unresolved"
        # voxels that a reseed searched and a later revisit or a failed fit overwrote
        if nev > 100_000:
            out[int(idx[e["row"][i]])] = dict(
                row=int(e["row"][i]),
                err=float(e["err"][i]),
                source=str(src),
                evals=int(nev),
                hit=float(rec[str(i)]["hit_ratio"]),
                grain=int(e["grain"][i]),
            )
    return out


def init(kind: str) -> None:
    bench().init_worker(kind)


def env(kind: str, voxel: int) -> Dict[str, Any]:
    """Environment of one voxel with its SmoothCost (observed rays built)."""
    B = bench()
    e = B.voxel_env(voxel)
    e["sc"] = B.lazy(e, "sc")
    return e


def symmetry_of(e: Dict[str, Any]) -> Any:
    return e["rec"].setup.exp_setup.get_sample_symmetry()
