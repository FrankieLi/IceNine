"""fine1: the pure parts of the generation-7 sample (numpy and scipy only, no file access).

``fine1_make_grid.py``, ``fine1_make_sample.py`` and ``fine1_gradient_check.py`` use these;
``tests/test_fine1.py`` tests them on hand-made triangles.

Grid rule: every triangle of a ``.mic`` grid is split into four, the way C++ splits it
(``ReconstructionStrategies.cpp`` ``SubdivideTriangles``, ``SetMinimumResolution``): three corner
triangles of the parent's orientation and the flipped middle one. A ``.mic`` row holds the
triangle's left vertex (x, y), a direction code (1 up, 2 down) and a generation; the side is the
header length over 2^generation. An up triangle has its vertices at (x, y), (x + s, y) and
(x + s/2, y + h); a down triangle at (x, y), (x + s/2, y - h) and (x + s, y); h = s sqrt(3)/2.
"""

from typing import Any, Dict, List, Sequence, Tuple

import numpy as np
from scipy.spatial import cKDTree

SQRT3_2 = float(np.sqrt(3.0) / 2.0)
UP, DOWN = 1, 2
GRADIENT_CAP_DEG = 0.45  # the gradient cap: largest neighbour-pair misorientation of a grain
GATE_PAIR_MAX_DEG = 0.5  # gate: the maximum over all neighbour pairs stays below this
GATE_SPREAD_MIN_DEG = 3.0  # gate: the largest grain spread reaches at least this
HIST_EDGES_DEG = (0.0, 0.01, 0.02, 0.05, 0.1, 0.2, 0.3, 0.4, 0.45, 0.5, 1.0e9)


# ---------------------------------------------------------------------------------------- grid


def split_once(xy: np.ndarray, up: np.ndarray, side: float) -> Tuple[np.ndarray, np.ndarray, float]:
    """Split every triangle (left vertex ``xy`` (n, 2), ``up`` (n,) bool, side ``side``) into four.

    Returns (xy (4n, 2), up (4n,), side / 2). Per parent the order is that of C++'s
    ``SubdivideTriangles``: the corner at the parent's left vertex, the corner one child side to
    its right, the corner at the apex (up parent: above; down parent: below), then the flipped
    middle triangle (it takes the apex corner's left vertex and the opposite direction).
    """
    xy = np.asarray(xy, dtype=np.float64).reshape(-1, 2)
    up = np.asarray(up, dtype=bool)
    c = side / 2.0
    h = c * SQRT3_2
    n = len(xy)
    out_xy = np.empty((n, 4, 2))
    out_up = np.empty((n, 4), dtype=bool)
    sign = np.where(up, 1.0, -1.0)  # the apex corner lies above (up) or below (down)
    out_xy[:, 0] = xy
    out_xy[:, 1] = xy + np.array([c, 0.0])
    out_xy[:, 2, 0] = xy[:, 0] + c / 2.0
    out_xy[:, 2, 1] = xy[:, 1] + sign * h
    out_xy[:, 3] = out_xy[:, 2]
    out_up[:, 0] = out_up[:, 1] = out_up[:, 2] = up
    out_up[:, 3] = ~up
    return out_xy.reshape(-1, 2), out_up.reshape(-1), c


def cpp_subdivide(
    x: float, y: float, up: bool, side: float, gen: int, gen_max: int
) -> List[Tuple[float, float, bool]]:
    """A literal port of C++ ``SubdivideTriangles`` for one parent triangle: the (left vertex x,
    left vertex y, points up) of its children at generation ``gen_max``, in C++'s order.
    ``side`` is the parent's side length (the header length over 2^gen)."""
    init = side * 2.0**gen  # the header side length (generation 0)
    up_diag = np.array([init / 2.0, init * np.sqrt(3.0) / 2.0])
    down_diag = np.array([init / 2.0, -init * np.sqrt(3.0) / 2.0])
    primary, secondary = (up_diag, down_diag) if up else (down_diag, up_diag)
    primary = primary / 2.0**gen_max
    secondary = secondary / 2.0**gen_max
    new_side = init / 2.0**gen_max
    left = np.array([x, y], dtype=np.float64)
    x_step = np.array([new_side, 0.0])
    n_points = int(2 ** (gen_max - gen))
    out: List[Tuple[float, float, bool]] = []
    for i in range(n_points):
        for j in range(n_points - i):
            cv = left + j * x_step
            out.append((float(cv[0]), float(cv[1]), up))
            if i > 0:
                out.append((float(cv[0]), float(cv[1]), not up))
        left = left + primary
    return out


def triangle_vertices(xy: np.ndarray, up: np.ndarray, side: Any) -> np.ndarray:
    """Vertices (n, 3, 2) of the triangles, as the readers build them. ``side``: scalar or (n,)."""
    xy = np.asarray(xy, dtype=np.float64).reshape(-1, 2)
    up = np.asarray(up, dtype=bool)
    s = np.broadcast_to(np.asarray(side, dtype=np.float64), (len(xy),))
    h = s * SQRT3_2
    out = np.empty((len(xy), 3, 2))
    out[:, 0] = xy
    out[:, 1, 0] = np.where(up, xy[:, 0] + s, xy[:, 0] + s / 2.0)
    out[:, 1, 1] = np.where(up, xy[:, 1], xy[:, 1] - h)
    out[:, 2, 0] = np.where(up, xy[:, 0] + s / 2.0, xy[:, 0] + s)
    out[:, 2, 1] = np.where(up, xy[:, 1] + h, xy[:, 1])
    return out


def triangle_area(verts: np.ndarray) -> np.ndarray:
    """Areas (n,) of triangles given as vertices (n, 3, 2)."""
    a, b, c = verts[:, 0], verts[:, 1], verts[:, 2]
    return (
        np.abs(
            (b[:, 0] - a[:, 0]) * (c[:, 1] - a[:, 1]) - (c[:, 0] - a[:, 0]) * (b[:, 1] - a[:, 1])
        )
        / 2.0
    )


def points_in_triangles(points: np.ndarray, verts: np.ndarray, tol: float) -> np.ndarray:
    """Bool (n,): point i lies inside (or within ``tol`` of) triangle i (``verts`` (n, 3, 2))."""
    a, b, c = verts[:, 0], verts[:, 1], verts[:, 2]

    def cross(p: np.ndarray, q: np.ndarray, r: np.ndarray) -> np.ndarray:
        return (q[:, 0] - p[:, 0]) * (r[:, 1] - p[:, 1]) - (q[:, 1] - p[:, 1]) * (r[:, 0] - p[:, 0])

    d1, d2, d3 = cross(a, b, points), cross(b, c, points), cross(c, a, points)
    # tol is a length; the cross product has the unit length^2, so scale it by an edge length
    e = np.linalg.norm(b - a, axis=1) * tol
    neg = (d1 < -e) | (d2 < -e) | (d3 < -e)
    pos = (d1 > e) | (d2 > e) | (d3 > e)
    return ~(neg & pos)


def grid_checks(
    parent_xy: np.ndarray,
    parent_up: np.ndarray,
    parent_side: float,
    child_xy: np.ndarray,
    child_up: np.ndarray,
    child_side: float,
    tol_mm: float = 1e-6,
) -> Dict[str, Any]:
    """The grid validity checks (plan): the count is four times the parent count, the area is the
    same, and every child lies inside its parent (children are in groups of four per parent, in
    parent order). Returns plain numbers and a bool ``ok``."""
    n = len(parent_xy)
    pv = triangle_vertices(parent_xy, parent_up, parent_side)
    cv = triangle_vertices(child_xy, child_up, child_side)
    parent_of = np.repeat(np.arange(n), 4) if len(child_xy) == 4 * n else None
    inside = False
    if parent_of is not None:
        inside = bool(np.all(points_in_triangles(cv[:, 0], pv[parent_of], tol_mm))) and bool(
            np.all(points_in_triangles(cv[:, 1], pv[parent_of], tol_mm))
        )
        inside = inside and bool(np.all(points_in_triangles(cv[:, 2], pv[parent_of], tol_mm)))
    pa = float(triangle_area(pv).sum())
    ca = float(triangle_area(cv).sum())
    out: Dict[str, Any] = dict(
        n_parent=n,
        n_child=len(child_xy),
        count_is_4x=bool(len(child_xy) == 4 * n),
        area_parent_mm2=pa,
        area_child_mm2=ca,
        area_rel_diff=float(abs(ca - pa) / pa) if pa > 0 else float("nan"),
        children_inside_parent=inside,
    )
    out["ok"] = bool(out["count_is_4x"] and inside and out["area_rel_diff"] < 1e-9)
    return out


# ----------------------------------------------------------------------------- neighbour pairs


def edge_pairs(centres: np.ndarray, side: float) -> np.ndarray:
    """Pairs (m, 2) of triangles of one side length that share an edge. Edge neighbours have
    their centres side/sqrt(3) apart (0.577 side); the next distances are 1.0 side (same
    direction) and 1.155 side (vertex only), so a cut at 0.7 side selects exactly the edge
    neighbours."""
    tree = cKDTree(np.asarray(centres, dtype=np.float64))
    p = tree.query_pairs(0.7 * side, output_type="ndarray")
    return np.asarray(p, dtype=np.int64).reshape(-1, 2)


def relative_angle_deg(Ma: np.ndarray, Mb: np.ndarray) -> np.ndarray:
    """Rotation angle (degrees) of Ma Mb^T, for stacks (n, 3, 3); accurate for tiny angles
    (arctan2 of the skew part and the trace part)."""
    R = np.einsum("nij,nkj->nik", Ma, Mb)
    skew = 0.5 * np.stack(
        [R[:, 2, 1] - R[:, 1, 2], R[:, 0, 2] - R[:, 2, 0], R[:, 1, 0] - R[:, 0, 1]], axis=1
    )
    sin_t = np.linalg.norm(skew, axis=1)
    cos_t = 0.5 * (np.trace(R, axis1=1, axis2=2) - 1.0)
    return np.degrees(np.arctan2(sin_t, cos_t))


# ----------------------------------------------------------------------------- the cap


def capped_gradient_orientations(
    centres: np.ndarray,
    grain_ids: np.ndarray,
    grain_orientations: Dict[int, np.ndarray],
    levels_deg: Dict[int, float],
    rng: np.random.Generator,
    pairs: np.ndarray,
    cap_deg: float = GRADIENT_CAP_DEG,
) -> Tuple[np.ndarray, List[Dict[str, Any]]]:
    """``icenine.gradient_sample.gradient_orientations`` (same draws, same model: one 3 by 2
    normal matrix per grain in increasing id order, rotation vector linear in position, scaled
    so that the grain's largest vector equals its level, left-multiplied onto the grain
    orientation) plus the gradient cap.

    The cap: after a grain's vectors are scaled to its level, the largest misorientation over
    the grain's neighbour pairs (``pairs``: edge-sharing voxel pairs; only pairs inside the grain
    count) is measured. If it is ``cap_deg`` or more, all the grain's vectors are scaled down by
    one factor, repeated until that maximum is ``cap_deg``. With no grain over the cap the result
    equals ``gradient_orientations`` bit for bit. Returns (orientations (n, 3, 3), one record per
    capped grain: id, drawn level, drawn largest step, kept largest step, kept spread, scale).
    """
    from icenine.gradient_sample import rotation_matrices

    xy = np.asarray(centres, dtype=np.float64)[:, :2]
    ids = np.asarray(grain_ids)
    out = np.empty((len(ids), 3, 3), dtype=np.float64)
    pairs = np.asarray(pairs, dtype=np.int64).reshape(-1, 2)
    same = ids[pairs[:, 0]] == ids[pairs[:, 1]] if len(pairs) else np.zeros(0, dtype=bool)
    pairs = pairs[same]
    pair_grain = ids[pairs[:, 0]] if len(pairs) else np.zeros(0, dtype=ids.dtype)
    capped: List[Dict[str, Any]] = []
    for g in sorted(int(x) for x in np.unique(ids)):
        M = rng.standard_normal((3, 2))
        sel = np.nonzero(ids == g)[0]
        Rg = np.asarray(grain_orientations[g], dtype=np.float64)
        out[sel] = Rg
        level = float(levels_deg[g])
        if level == 0.0 or len(sel) < 2:
            continue
        d = xy[sel] - xy[sel].mean(axis=0)
        rv = d @ M.T
        top = np.linalg.norm(rv, axis=1).max()
        if top <= 0.0:
            continue
        rv = rv * (np.radians(level) / top)
        gp = pairs[pair_grain == g]
        if len(gp):
            local = {int(v): i for i, v in enumerate(sel)}
            la = np.array([local[int(a)] for a in gp[:, 0]])
            lb = np.array([local[int(b)] for b in gp[:, 1]])
            Mrot = rotation_matrices(rv)
            drawn = float(relative_angle_deg(Mrot[la], Mrot[lb]).max())
            if drawn >= cap_deg:
                scale = 1.0
                cur = drawn
                for _ in range(50):
                    f = cap_deg / cur
                    scale *= f
                    Mrot = rotation_matrices(rv * scale)
                    cur = float(relative_angle_deg(Mrot[la], Mrot[lb]).max())
                    if abs(cur - cap_deg) <= 1e-12 * cap_deg:
                        break
                if cur > cap_deg:  # never leave it above the cap
                    scale *= cap_deg / cur * (1.0 - 1e-12)
                    Mrot = rotation_matrices(rv * scale)
                    cur = float(relative_angle_deg(Mrot[la], Mrot[lb]).max())
                rv = rv * scale
                Rv = rotation_matrices(rv)
                kept_spread = float(
                    np.degrees(np.linalg.norm(rv, axis=1).max())
                )  # rotation-vector length is the angle to the centre orientation
                capped.append(
                    dict(
                        id=g,
                        n_voxels=int(len(sel)),
                        level_drawn_deg=level,
                        max_step_drawn_deg=drawn,
                        max_step_kept_deg=cur,
                        level_kept_deg=kept_spread,
                        scale=float(scale),
                    )
                )
                out[sel] = Rv @ Rg
                continue
        out[sel] = rotation_matrices(rv) @ Rg
    return out, capped


# ----------------------------------------------------------------------------- the gate


def percentiles(values: np.ndarray, qs: Sequence[float] = (50, 90, 99, 99.9)) -> Dict[str, float]:
    """Percentiles of ``values`` as {"p50": ..}."""
    v = np.asarray(values, dtype=np.float64)
    return {f"p{q:g}": float(np.percentile(v, q)) for q in qs}


def histogram(values: np.ndarray, edges: Sequence[float] = HIST_EDGES_DEG) -> List[Dict[str, Any]]:
    """Counts per half-open bin [lo, hi) of the edges."""
    v = np.asarray(values, dtype=np.float64)
    counts, _ = np.histogram(v, bins=np.asarray(edges, dtype=np.float64))
    return [
        dict(lo=float(lo), hi=float(hi), count=int(c))
        for lo, hi, c in zip(edges[:-1], edges[1:], counts)
    ]


def gate(
    n_pairs: int,
    max_pair_deg: float,
    max_spread_deg: float,
    pair_limit: float = GATE_PAIR_MAX_DEG,
    spread_min: float = GATE_SPREAD_MIN_DEG,
) -> Dict[str, Any]:
    """The hard gate before the freeze: there are pairs, the maximum neighbour-pair misorientation
    is under ``pair_limit`` and the largest grain spread is at least ``spread_min``."""
    ok_pair = bool(n_pairs > 0 and max_pair_deg < pair_limit)
    ok_spread = bool(max_spread_deg >= spread_min)
    return dict(
        n_pairs=int(n_pairs),
        max_pair_deg=float(max_pair_deg),
        max_pair_limit_deg=pair_limit,
        max_pair_ok=ok_pair,
        max_spread_deg=float(max_spread_deg),
        max_spread_min_deg=spread_min,
        max_spread_ok=ok_spread,
        passed=bool(ok_pair and ok_spread),
    )


def nearest_centre_rows(left_xy: np.ndarray, n: int) -> np.ndarray:
    """Sorted row indices of the ``n`` voxels nearest the layer centre (0, 0), by the left
    vertex (x, y) of the .mic row: the expression of ``render_full.py --n-voxels``
    (``np.argsort`` of the distance, first ``n``)."""
    d = np.linalg.norm(np.asarray(left_xy, dtype=np.float64)[:, :2] - np.array([0.0, 0.0]), axis=1)
    return np.sort(np.argsort(d)[:n])
