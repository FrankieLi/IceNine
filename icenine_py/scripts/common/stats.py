"""Shared statistics helpers for the experiment scripts (numpy/scipy only).

Import convention (same as scripts/findoptimal_robustness): the importing script puts
``icenine_py/scripts/common`` on ``sys.path`` and does ``import stats``.
"""

import sys
from pathlib import Path
from typing import Sequence, Tuple

import numpy as np
from scipy.stats import binomtest

_CSL_DIR = Path(__file__).resolve().parent.parent / "findoptimal_robustness"


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


def win_rate(err_a: Sequence[float], err_b: Sequence[float], tie: float = 0.002) -> float:
    """Fraction of pairs where A has the smaller error; |a - b| <= tie counts one half.

    Pairs with a NaN on either side are dropped; returns NaN if no pair remains.
    """
    a = np.asarray(err_a, dtype=np.float64)
    b = np.asarray(err_b, dtype=np.float64)
    if a.shape != b.shape:
        raise ValueError(f"shape mismatch: {a.shape} vs {b.shape}")
    keep = ~(np.isnan(a) | np.isnan(b))
    a, b = a[keep], b[keep]
    if a.size == 0:
        return float("nan")
    wins = (a < b) & (np.abs(a - b) > tie)
    ties = np.abs(a - b) <= tie
    return float((wins.sum() + 0.5 * ties.sum()) / a.size)


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


def misorientation_deg_cubic(R_a: np.ndarray, R_b: np.ndarray) -> np.ndarray:
    """Cubic-symmetry-reduced misorientation angle in degrees of rotation matrices (batched).

    Thin wrapper over ``csl.reduced_misorientation_deg`` (scripts/findoptimal_robustness).
    """
    if str(_CSL_DIR) not in sys.path:
        sys.path.insert(0, str(_CSL_DIR))
    import csl  # noqa: PLC0415

    return np.asarray(csl.reduced_misorientation_deg(R_a, R_b))
