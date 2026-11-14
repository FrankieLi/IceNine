"""AR3 step 11: the pure parts of the sample builder for ``held`` (numpy and scipy only).

``held_make_sample.py`` uses these; ``tests/test_ar3_held_sample.py`` tests them on small inputs.
They are the steps of ``start_make_mid.py`` (the builder of ``mid``) written as functions, with the
same calls and the same order of random draws, so that ``mid`` is rebuilt from the same code path
when the arguments are 120 cells and seed 0.
"""

from typing import List, Sequence, Tuple

import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
from scipy.spatial import cKDTree


def draw_cell_seeds(n_voxels: int, n_cells: int, seed: int) -> np.ndarray:
    """Sorted indices of the voxel centres that are the Voronoi seeds (``default_rng(seed)``)."""
    rng = np.random.default_rng(seed)
    return np.sort(rng.choice(n_voxels, size=n_cells, replace=False))


def voronoi_cells(centres: np.ndarray, seed_idx: np.ndarray) -> np.ndarray:
    """Each voxel joins the nearest seed: the cell number (position in ``seed_idx``) per voxel."""
    return np.asarray(cKDTree(centres[seed_idx]).query(centres, k=1)[1])


def connected_grains(centres: np.ndarray, cell: np.ndarray, radius: float) -> np.ndarray:
    """Grain number per voxel: the connected pieces of one cell (centres within ``radius`` are
    neighbours). A cell that is split by the grid gives more than one grain."""
    n = len(centres)
    pairs = cKDTree(centres).query_pairs(radius, output_type="ndarray")
    same = cell[pairs[:, 0]] == cell[pairs[:, 1]]
    p = pairs[same]
    adj = coo_matrix((np.ones(len(p)), (p[:, 0], p[:, 1])), shape=(n, n))
    _ncomp, comp = connected_components(adj, directed=False)
    return comp.astype(np.int64)


def draw_levels(
    seed: int, n_grains: int, levels: Sequence[float]
) -> Tuple[np.random.Generator, List[float]]:
    """One gradient level per grain, drawn uniformly from ``levels`` with ``default_rng(seed)``.
    Returns the generator (its state goes on to ``gradient_orientations``) and the levels."""
    grng = np.random.default_rng(seed)
    idx = grng.integers(0, len(levels), size=n_grains)
    return grng, [float(levels[i]) for i in idx]


def grains_to_redraw(
    a: Sequence[int], b: Sequence[int], angle_deg: Sequence[float], min_sep_deg: float
) -> List[int]:
    """The later grain of every adjacent pair (a, b) whose misorientation is under ``min_sep_deg``,
    sorted and without repeats."""
    return sorted({int(max(i, j)) for i, j, d in zip(a, b, angle_deg) if d < min_sep_deg})
