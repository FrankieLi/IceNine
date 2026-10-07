"""Shared paths, dataset loaders and metric helpers of the Q_max-8 cost proxy study.

Reuses scripts/findoptimal_robustness (common, features, e2_models) unchanged; this module only
adds the proxy's own cache / output directories and the aligned loading of the E2 dataset."""

import sys
from pathlib import Path
from typing import Any, Dict, List

import numpy as np

HERE = Path(__file__).resolve().parent
ICENINE_PY = HERE.parents[1]
FO = ICENINE_PY / "scripts" / "findoptimal_robustness"
sys.path.insert(0, str(FO))
sys.path.insert(0, str(HERE))

import common as C  # noqa: E402
import e2_models as M  # noqa: E402

CACHE = HERE / "cache"  # gitignored
OUT = ICENINE_PY / "benchmarks" / "coarse_proxy"
Q_LEVELS = (4.0, 5.0)  # |q| cut-offs of the low-Q extractors (Q3 contains no reflection)


def voxels() -> List[int]:
    info = dict(np.load(C.OUT_DIR / "voxels.npz"))
    return [int(v) for v in info["voxel_indices"]][: C.N_VOXELS]


def tasks() -> List[Any]:
    """(vidx, vpos, variant) of every E2 task file present, in `e2_models.load_dataset` order."""
    out = []
    for vpos, v in enumerate(voxels()):
        for var in C.VARIANTS:
            if (C.CACHE_DIR / "e2" / f"v{v}_{var}.npz").exists():
                out.append((v, vpos, var))
    return out


def load_per_task(directory: Path, keys: List[str]) -> Dict[str, np.ndarray]:
    """Concatenate the arrays `keys` of directory/v{vidx}_{variant}.npz over `tasks()`."""
    parts: Dict[str, List[np.ndarray]] = {k: [] for k in keys}
    for v, _, var in tasks():
        d = np.load(directory / f"v{v}_{var}.npz")
        for k in keys:
            parts[k].append(d[k])
    return {k: np.concatenate(p) for k, p in parts.items()}


def load_e2_aligned() -> Dict[str, np.ndarray]:
    """The E2 dataset (same order as e2_models.load_dataset) plus the candidate matrices R."""
    D = M.load_dataset()
    D["R"] = load_per_task(C.CACHE_DIR / "e2", ["R"])["R"]
    assert len(D["R"]) == len(D["err"])
    return D


def cost_col(D: Dict[str, np.ndarray]) -> np.ndarray:
    """The free Q_max-8 local cost of each candidate (cost_local, the second-to-last E2 column)."""
    return D["X"][:, -2]
