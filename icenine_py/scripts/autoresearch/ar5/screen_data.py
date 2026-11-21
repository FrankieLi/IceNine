"""AR5 step 21: load the stored inputs of one sample for arm J5S (script layer; no image, no
reconstructor). Step 20's ``cmarepair_data`` with this step's guards (``screen_guard``).

``load_inputs(sample, truth=False)`` reads, from this step's cache (links, read only, made by
``screen_link_cache.sh``; for ``held5`` the files of its own runs): ``before_J5.npy`` (J5's map
before the pass) and ``D4_J5/tries.npz`` (J5's pass tries and the old ends' scores). The neighbour
relation is the pass's own (``mic.get_neighbors``, radius 2 side lengths), read from the sample's
grid (``gradient_grid.mic``: no orientation of the truth). ``truth=False`` never opens the truth:
the J5C and J5S runs use it, so that ``held5``'s truth stays closed until its runs end.

``stored_stats(sample)``: step 16's (held4: step 20's) stored totals and wrong counts of J5 and S8.
"""

import json
from pathlib import Path
from typing import Any, Dict, List

import numpy as np
import screen_guard as SG

from icenine.mic_file import MicFile
from icenine.ray_indexing import proper_rotations
from icenine.symmetry import create_cubic_symmetry

ICE = SG.ICE
AR3 = ICE / "benchmarks" / "autoresearch" / "ar3"
STATS = {
    "held4": ICE / "benchmarks" / "autoresearch" / "ar5" / "cmarepair" / "stats_held4.json",
}
STATS16 = ICE / "benchmarks" / "autoresearch" / "ar5" / "high500"


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
    mic = MicFile.read(str(SG.guarded(sample_dir(s) / "gradient_grid.mic")))
    return dict(nbrs=_neighbours(mic), n=len(mic.voxels))


def read_truth(s: str) -> Dict[str, Any]:
    """Truth orientations and the neighbour lists of the sample (scoring only)."""
    mic = MicFile.read(str(SG.guarded(sample_dir(s) / "gradient_truth.mic")))
    R_t = np.array([v.orientation for v in mic.voxels], dtype=np.float64)
    return dict(R_t=R_t, nbrs=_neighbours(mic), n=len(mic.voxels))


def stored_stats(s: str) -> Dict[str, Any]:
    """Totals (evaluations with the pass), wrong voxels and median error of J5 and S8."""
    p = STATS.get(s, STATS16 / f"stats_{s}.json")
    d = json.loads(p.read_text())
    if s == "held4":  # step 20's flat file
        return dict(
            total=int(d["j5_total"]),
            wrong=int(d["wrong_J5"]),
            median_err_deg=float(d["median_J5"]),
            s8_total=int(d["s8_total"]),
            s8_wrong=int(d["wrong_S8"]),
            s8_median_err_deg=float(d["median_S8"]),
        )
    j5, s8 = d["arms"]["J5"], d["arms"]["S8"]
    return dict(
        total=int(j5["stages"]["total"]),
        wrong=int(j5["after"]["wrong"]),
        median_err_deg=float(j5["after"]["median_err_deg"]),
        s8_total=int(s8["stages"]["total"]),
        s8_wrong=int(s8["after"]["wrong"]),
        s8_median_err_deg=float(s8["after"]["median_err_deg"]),
    )


def load_inputs(s: str, truth: bool = False) -> Dict[str, Any]:
    """J5's map before the pass, its stored pass tries and the old ends' scores; the neighbour
    lists from the grid (``truth=False``) or the truth map (``truth=True``, scoring only)."""
    p = SG.cache_dir() / f"passes_{s}"
    Rb = np.load(SG.guarded(p / "before_J5.npy")).astype(np.float64)
    z = np.load(SG.guarded(p / "D4_J5" / "tries.npz"))
    nb = read_truth(s) if truth else read_grid(s)
    return dict(
        sample=s,
        sym=sym(),
        R_before=Rb,
        old_score=np.asarray(z["old_score"], dtype=float),
        stored=dict(
            voxel=z["voxel"].astype(np.int64),
            cand=z["cand"].astype(np.int64),
            end=z["end"],
            cost=z["cost"],
            hit=z["hit"],
            conf=z["conf"],
            evals=z["evals"].astype(np.int64),
        ),
        **nb,
    )


def orthonormal(R: np.ndarray) -> np.ndarray:
    """Nearest rotation of each matrix (the stored tries' ends are float32)."""
    U, _s, Vt = np.linalg.svd(np.asarray(R, dtype=np.float64))
    return U @ Vt
