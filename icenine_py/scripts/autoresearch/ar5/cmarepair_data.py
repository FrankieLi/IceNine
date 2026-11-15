"""AR5 step 20: load the stored state after the boundary pass of one sample and arm (script layer;
no study imports, no image, no reconstructor). Applies step 19's frozen flag (``repair_rules``,
imported through ``cmarepair_rules``).

``load(sample, arm, truth=True)`` reads, from this step's cache (links, read only, made by
``cmarepair_link_cache.sh``; for ``held4`` the files of its own runs): ``before_<arm>.npy`` (the
map before the pass), ``D4_<arm>/tries.npz`` (the tries and the old ends' scores). It rebuilds the
state after the pass with ``after_pass_state`` and the flag with ``contested``. The neighbour
relation is the pass's own (``mic.get_neighbors``, radius 2 side lengths; D-20261115-N7), read from
the sample's grid (``gradient_grid.mic``: no orientation of the truth) or, with ``truth=True``,
from the truth map (the same voxels; scoring needs its orientations). ``truth=False`` never opens
the truth: the repair run uses it, so that ``held4``'s truth stays closed until its runs end.
"""

import json
from pathlib import Path
from typing import Any, Dict, List

import cmarepair_guard as CG
import cmarepair_rules as CR
import numpy as np

from icenine.mic_file import MicFile
from icenine.ray_indexing import proper_rotations
from icenine.symmetry import create_cubic_symmetry

ICE = CG.ICE
AR3 = ICE / "benchmarks" / "autoresearch" / "ar3"
STATS16 = ICE / "benchmarks" / "autoresearch" / "ar5" / "high500"
STATS18 = ICE / "benchmarks" / "autoresearch" / "ar5" / "learnskip"


def sym() -> np.ndarray:
    return proper_rotations(create_cubic_symmetry(3.61))


def sample_dir(s: str) -> Path:
    return AR3 / "gradient" / "sample" if s == "grad" else AR3 / s / "gradient" / "sample"


def _neighbours(mic: Any) -> List[List[int]]:
    return [
        sorted(
            int(j)
            for j in mic.get_neighbors(v, radius=2.0 * mic.voxels[v].side_length)
            if int(j) != v
        )
        for v in range(len(mic.voxels))
    ]


def read_grid(s: str) -> Dict[str, Any]:
    """The neighbour lists of the sample, from its grid (no truth is opened)."""
    mic = MicFile.read(str(CG.guarded(sample_dir(s) / "gradient_grid.mic")))
    return dict(nbrs=_neighbours(mic), n=len(mic.voxels))


def read_truth(s: str) -> Dict[str, Any]:
    """Truth orientations and the neighbour lists of the sample (the pass's own relation)."""
    mic = MicFile.read(str(CG.guarded(sample_dir(s) / "gradient_truth.mic")))
    R_t = np.array([v.orientation for v in mic.voxels], dtype=np.float64)
    return dict(R_t=R_t, nbrs=_neighbours(mic), n=len(mic.voxels))


def stored_stats(s: str, arm: str) -> Dict[str, Any]:
    """Step 16's (or, for J5L on mid, step 18's) stored totals: evaluations with the pass, wrong
    voxels after the pass and the median error after it."""
    if arm == "J5L":
        d = json.loads((STATS18 / "stats_mid.json").read_text())
        a = d["arms"]["J5L"]
    else:
        d = json.loads((STATS16 / f"stats_{s}.json").read_text())
        a = d["arms"]["J5"]
    out = dict(
        total=int(a["stages"]["total"]),
        wrong=int(a["after"]["wrong"]),
        median_err_deg=float(a["after"]["median_err_deg"]),
    )
    if arm != "J5L":
        s8 = d["arms"]["S8"]
        out.update(
            s8_total=int(s8["stages"]["total"]),
            s8_wrong=int(s8["after"]["wrong"]),
            s8_median_err_deg=float(s8["after"]["median_err_deg"]),
        )
    return out


def load(s: str, arm: str = "J5", truth: bool = True) -> Dict[str, Any]:
    """The state after the pass and the flag of one sample and arm. ``truth=False``: no truth map
    and no S8 map is opened (the repair run's form); the neighbours come from the grid."""
    cache = CG.cache_dir()
    p = cache / f"passes_{s}"
    Rb = np.load(CG.guarded(p / f"before_{arm}.npy")).astype(np.float64)
    Ra_stored = np.load(CG.guarded(p / f"after_{arm}.npy")).astype(np.float64)
    z = np.load(CG.guarded(p / f"D4_{arm}" / "tries.npz"))
    # the pass's tries are stored as float32; the pass's rule reads them orthonormalized
    U, _s, Vt = np.linalg.svd(z["end"].astype(np.float64))
    t_end = U @ Vt
    y = sym()
    R, cost, hit = CR.after_pass_state(Rb, z["old_score"], z["voxel"], t_end, z["cost"], z["hit"])
    other = CR.other_end_hit(R, Rb, z["old_score"][:, 1], z["voxel"], t_end, z["hit"], y)
    nb = read_truth(s) if truth else read_grid(s)
    flag, ratio, margin = CR.contested(hit, other, nb["nbrs"])
    flagged = [int(v) for v in np.where(flag)[0]]
    cands = CR.repair_candidates(R, nb["nbrs"], y, flagged)
    out: Dict[str, Any] = dict(
        sample=s,
        arm=arm,
        sym=y,
        R_before=Rb,
        R=R,
        cost=cost,
        hit=hit,
        other=other,
        ratio=ratio,
        margin=margin,
        flagged=flagged,
        cands=cands,
        fits=sum(len(cands[v]) + 1 for v in flagged),
        after_equals_stored=float(np.abs(R - Ra_stored).max()),
        R_after_stored=Ra_stored,
        **nb,
    )
    s8_path = p / "after_S8.npy"
    if truth and arm == "J5" and s8_path.exists():
        out["R_s8"] = np.load(CG.guarded(s8_path)).astype(np.float64)
    return out
