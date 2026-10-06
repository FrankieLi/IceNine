"""Shared statistics helpers for the experiment scripts (numpy/scipy only).

Import convention (same as scripts/findoptimal_robustness): the importing script puts
``icenine_py/scripts/common`` on ``sys.path`` and does ``import stats``.
"""

import importlib.util
import sys
from pathlib import Path
from typing import Any, Sequence, Tuple

import numpy as np
from scipy.stats import binomtest

_CSL_PATH = Path(__file__).resolve().parent.parent / "findoptimal_robustness" / "csl.py"


def wilson(k: int, n: int, z: float = 1.96) -> Tuple[float, float]:
    """Wilson score interval for k successes in n trials; (nan, nan) when n == 0."""
    if n <= 0:
        return (float("nan"), float("nan"))
    p = k / n
    denom = 1.0 + z * z / n
    centre = (p + z * z / (2 * n)) / denom
    half = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / denom
    return (float(max(0.0, centre - half)), float(min(1.0, centre + half)))


def mcnemar_exact(b: int, c: int) -> float:
    """Two-sided exact McNemar p-value from the discordant counts b and c."""
    n = b + c
    if n == 0:
        return 1.0
    return float(binomtest(min(b, c), n, 0.5, alternative="two-sided").pvalue)


def paired_discordant(wrong_a: Sequence[bool], wrong_b: Sequence[bool]) -> Tuple[int, int]:
    """(b, c): b = runs wrong under A only, c = runs wrong under B only."""
    a = np.asarray(wrong_a, dtype=bool)
    b = np.asarray(wrong_b, dtype=bool)
    if a.shape != b.shape:
        raise ValueError(f"shape mismatch: {a.shape} vs {b.shape}")
    return int(np.sum(a & ~b)), int(np.sum(~a & b))


def win_rate(
    err_a: Sequence[float], err_b: Sequence[float], tie: float = 0.002
) -> Tuple[float, float, int]:
    """(rate, tie_fraction, n): the fraction of pairs where A has the smaller error, with a tie
    (|a - b| < tie, strict) counting one half.

    Pairs with a non-finite value on either side are dropped; n is the number of pairs kept.
    With no pair left, returns (nan, nan, 0).
    """
    a = np.asarray(err_a, dtype=np.float64)
    b = np.asarray(err_b, dtype=np.float64)
    if a.shape != b.shape:
        raise ValueError(f"shape mismatch: {a.shape} vs {b.shape}")
    keep = np.isfinite(a) & np.isfinite(b)
    a, b = a[keep], b[keep]
    if a.size == 0:
        return (float("nan"), float("nan"), 0)
    d = a - b
    ties = np.abs(d) < tie
    wins = (d < 0) & ~ties
    return (float((wins.sum() + 0.5 * ties.sum()) / a.size), float(ties.mean()), int(a.size))


def reorder(
    values: np.ndarray, have_ids: Sequence[int], want_ids: Sequence[int], axis: int = 0
) -> np.ndarray:
    """Re-index ``values`` (axis ordered as ``have_ids``) into ``want_ids`` order.

    ``want_ids`` may be a subset of ``have_ids``. Raises ValueError on a missing or duplicated
    id, so a voxel-order mismatch cannot pass silently.
    """
    values = np.asarray(values)
    have = list(have_ids)
    if values.shape[axis] != len(have):
        raise ValueError(f"axis {axis} has {values.shape[axis]} entries but {len(have)} ids")
    if len(set(have)) != len(have):
        raise ValueError("have_ids contains duplicates")
    pos = {int(i): k for k, i in enumerate(have)}
    missing = [int(i) for i in want_ids if int(i) not in pos]
    if missing:
        raise ValueError(f"ids missing from have_ids: {missing[:10]}")
    idx = np.array([pos[int(i)] for i in want_ids], dtype=np.int64)
    return np.take(values, idx, axis=axis)


def _load_csl() -> Any:
    """Load csl.py under a private module name (no sys.path change)."""
    name = "_icenine_common_csl"
    if name in sys.modules:
        return sys.modules[name]
    spec = importlib.util.spec_from_file_location(name, _CSL_PATH)
    assert spec is not None and spec.loader is not None
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


def misorientation_deg_cubic(R_a: np.ndarray, R_b: np.ndarray) -> np.ndarray:
    """Cubic-symmetry-reduced misorientation angle in degrees of rotation matrices (batched).

    Thin wrapper over ``csl.reduced_misorientation_deg`` (scripts/findoptimal_robustness).
    """
    return np.asarray(_load_csl().reduced_misorientation_deg(R_a, R_b))
