"""AR5 step 19: load the stored state after the boundary pass of one sample and arm (script layer;
no study imports, no image, no reconstructor). Applies the frozen ``repair_rules``.

``load(sample, arm)`` reads, from this step's cache (links, read only, made by
``repair_link_cache.sh``): ``before_<arm>.npy`` (the map before the pass), ``D4_<arm>/tries.npz``
(the tries and the old ends' scores) and the sample's truth map (scoring only). It rebuilds the
state after the pass with ``repair_rules.after_pass_state`` and the flag with
``repair_rules.contested``; the neighbour relation is the pass's own (``mic.get_neighbors`` with
radius 2 side lengths; D-20261115-N7).
"""

import json
from pathlib import Path
from typing import Any, Dict, List

import numpy as np
import repair_guard as RG
import repair_rules as RR

from icenine.mic_file import MicFile
from icenine.ray_indexing import proper_rotations
from icenine.symmetry import create_cubic_symmetry

ICE = RG.ICE
AR3 = ICE / "benchmarks" / "autoresearch" / "ar3"
STATS16 = ICE / "benchmarks" / "autoresearch" / "ar5" / "high500"
STATS18 = ICE / "benchmarks" / "autoresearch" / "ar5" / "learnskip"


def sym() -> np.ndarray:
    return proper_rotations(create_cubic_symmetry(3.61))


def sample_dir(s: str) -> Path:
    return AR3 / "gradient" / "sample" if s == "grad" else AR3 / s / "gradient" / "sample"


def read_truth(s: str) -> Dict[str, Any]:
    """Truth orientations and the neighbour lists of the sample (the pass's own relation)."""
    p = RG.guarded(sample_dir(s) / "gradient_truth.mic")
    mic = MicFile.read(str(p))
    n = len(mic.voxels)
    nbrs: List[List[int]] = [
        sorted(
            int(j)
            for j in mic.get_neighbors(v, radius=2.0 * mic.voxels[v].side_length)
            if int(j) != v
        )
        for v in range(n)
    ]
    R_t = np.array([v.orientation for v in mic.voxels], dtype=np.float64)
    return dict(R_t=R_t, nbrs=nbrs, n=n)


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


def load(s: str, arm: str = "J5") -> Dict[str, Any]:
    """The state after the pass and the flag of one sample and arm (truth read for scoring)."""
    cache = RG.cache_dir()
    p = cache / f"passes_{s}"
    Rb = np.load(RG.guarded(p / f"before_{arm}.npy")).astype(np.float64)
    Ra_stored = np.load(RG.guarded(p / f"after_{arm}.npy")).astype(np.float64)
    z = np.load(RG.guarded(p / f"D4_{arm}" / "tries.npz"))
    # the pass's tries are stored as float32; the pass's rule reads them orthonormalized
    U, _s, Vt = np.linalg.svd(z["end"].astype(np.float64))
    t_end = U @ Vt
    y = sym()
    R, cost, hit = RR.after_pass_state(Rb, z["old_score"], z["voxel"], t_end, z["cost"], z["hit"])
    other = RR.other_end_hit(R, Rb, z["old_score"][:, 1], z["voxel"], t_end, z["hit"], y)
    truth = read_truth(s)
    flag, ratio, margin = RR.contested(hit, other, truth["nbrs"])
    flagged = [int(v) for v in np.where(flag)[0]]
    cands = RR.repair_candidates(R, truth["nbrs"], y, flagged)
    s8_path = p / "after_S8.npy"
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
        **truth,
    )
    if arm == "J5" and s8_path.exists():
        out["R_s8"] = np.load(RG.guarded(s8_path)).astype(np.float64)
    return out
