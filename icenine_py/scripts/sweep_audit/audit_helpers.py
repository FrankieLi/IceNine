"""Pure helpers for the Phase B2 audit of the April 2026 HP sweep and the hybrid benchmark."""

import re
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

LOG_PAT = re.compile(
    r"^\s+(\w+)\s+hp=\s*(\d+)\s+vox=\s*(\d+)\s+pert=(\d+)°\s+misori=\s*([\d.]+)°"
    r"\s+q=([\d.\-]+)\s+evals=\s*(\d+)\s+([\d.]+)s"
)
LOG_COLS = ["opt", "hp", "vox", "pert", "mis", "q", "ev", "t"]


def parse_log_line(line: str) -> Optional[Dict[str, object]]:
    """One per-run line of hp_sweep_<example>.log -> dict, or None if the line is not a run line.
    The misorientation is printed to 0.001 deg."""
    m = LOG_PAT.match(line)
    if m is None:
        return None
    g = m.groups()
    return dict(
        opt=g[0],
        hp=int(g[1]),
        vox=int(g[2]),
        pert=int(g[3]),
        mis=float(g[4]),
        q=float(g[5]),
        ev=int(g[6]),
        t=float(g[7]),
    )


def success(mis: Sequence[float], thr: float) -> np.ndarray:
    """Boolean success mask: final misorientation strictly below thr (degrees)."""
    return np.asarray(mis, dtype=float) < thr


def split_complete_blocks(
    keys: Sequence[Tuple[int, float]], block_len: int = 0
) -> List[Tuple[int, int]]:
    """Split a CSV's rows into the consecutive runs of the benchmark that wrote them. One run lists
    (voxel, perturbation) keys in increasing lexicographic order, so a new block starts wherever a
    key is not greater than the previous one. Returns [(start, stop)]; block_len is unused (kept
    for call compatibility)."""
    blocks: List[Tuple[int, int]] = []
    start = 0
    for i in range(1, len(keys)):
        if tuple(keys[i]) <= tuple(keys[i - 1]):
            blocks.append((start, i))
            start = i
    if len(keys):
        blocks.append((start, len(keys)))
    return blocks


def quantiles(x: Sequence[float], qs: Sequence[float] = (0.25, 0.5, 0.75, 0.9)) -> List[float]:
    """Quantiles of the finite values of x; NaNs if none."""
    a = np.asarray(x, dtype=float)
    a = a[np.isfinite(a)]
    return [float(v) for v in np.quantile(a, qs)] if a.size else [float("nan")] * len(qs)


def pick_best(rows: List[Dict[str, float]]) -> Dict[str, float]:
    """The row with the largest 'k' (successes); ties: larger 'n_evals' is worse, then smaller
    'hp'. Deterministic and independent of input order."""
    return sorted(rows, key=lambda r: (-r["k"], r["n_evals"], r["hp"]))[0]
