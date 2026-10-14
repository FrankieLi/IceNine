"""Candidate overlap and de-duplication for the coarse levels of the no-start search.

Task C of the evaluation-reduction plan. Two parts:

1. A fast symmetry-reduced neighbour search for orientation matrices (``neighbour_lists``,
   ``greedy_nms``, ``component_labels``). It replaces the O(N^2 x |G|) reference
   ``scripts/findoptimal_robustness/csl.py::dedup_orientations`` and returns the same groups
   (tests/test_candidate_dedup.py). Method: every orientation becomes a unit quaternion; its
   representative in the fundamental zone is the member of its symmetry class with the largest
   |w|; a ``scipy.spatial.cKDTree`` holds the representatives plus, for the orientations near the
   zone boundary, the other class members that can lie within the tolerance of a point of the
   zone (the boundary copies, both quaternion signs); a radius query then finds every pair whose
   symmetry-reduced misorientation is at most the tolerance. The method does not assume a crystal
   system (any proper-rotation group given as quaternions).

2. ``DedupSpec`` (the value of the config key ``CandidateDedup``, off by default) and the two stage
   functions the level loop of ``AdaptiveVoxelReconstructor.reconstruct_voxel`` calls when the key
   is set: ``keep_before_quick_mc`` (clique filter C0, dedup C1) and ``distinct_order`` (C2).

Conventions: an orientation R maps crystal to sample coordinates and the crystal symmetry acts on
the right (R ~ R S), as in ``get_misorientation`` and ``csl.py``. Tolerances are angles in
radians; for unit quaternions the angle ``theta`` between two orientations has the chord
``2 sin(theta / 4)`` between q1 and the sign-aligned q2, which is the radius of the tree query.
"""

from dataclasses import dataclass
from typing import List, Optional, Tuple

import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
from scipy.spatial import cKDTree
from scipy.spatial.transform import Rotation

# ---------------------------------------------------------------------------
# Quaternion helpers (scalar first, [w, x, y, z])
# ---------------------------------------------------------------------------


def _matrix_to_wxyz(R: np.ndarray) -> np.ndarray:
    q = Rotation.from_matrix(np.asarray(R, dtype=np.float64)).as_quat()  # x, y, z, w
    return np.concatenate([q[..., 3:4], q[..., :3]], axis=-1)


def _qmul(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    """Hamilton product, broadcasting over leading axes."""
    w1, x1, y1, z1 = a[..., 0], a[..., 1], a[..., 2], a[..., 3]
    w2, x2, y2, z2 = b[..., 0], b[..., 1], b[..., 2], b[..., 3]
    return np.stack(
        [
            w1 * w2 - x1 * x2 - y1 * y2 - z1 * z2,
            w1 * x2 + x1 * w2 + y1 * z2 - z1 * y2,
            w1 * y2 - x1 * z2 + y1 * w2 + z1 * x2,
            w1 * z2 + x1 * y2 - y1 * x2 + z1 * w2,
        ],
        axis=-1,
    )


def _group_ops(symmetry_quats: Optional[np.ndarray]) -> np.ndarray:
    """The symmetry operators as unique quaternions (a quaternion and its negative are one
    rotation), the identity always included. Empty or None = no symmetry."""
    ident = np.array([[1.0, 0.0, 0.0, 0.0]])
    if symmetry_quats is None or len(symmetry_quats) == 0:
        return ident
    q = np.asarray(symmetry_quats, dtype=np.float64).reshape(-1, 4)
    q = np.concatenate([ident, q], axis=0)
    # canonical sign: the first component with magnitude above 1e-9 is positive
    lead = np.argmax(np.abs(q) > 1e-9, axis=1)
    sign = np.sign(q[np.arange(len(q)), lead])
    q = q * sign[:, None]
    _, first = np.unique(np.round(q, 9) + 0.0, axis=0, return_index=True)
    return q[np.sort(first)]


def chord_radius(tol_rad: float) -> float:
    """Euclidean distance between two unit quaternions (signs aligned) at angle ``tol_rad``."""
    return 2.0 * float(np.sin(max(float(tol_rad), 0.0) / 4.0))


# ---------------------------------------------------------------------------
# Symmetry-reduced neighbour search
# ---------------------------------------------------------------------------


def neighbour_lists(
    R: np.ndarray, tol_rad: float, symmetry_quats: Optional[np.ndarray] = None
) -> List[np.ndarray]:
    """For each orientation i, the indices j != i whose symmetry-reduced misorientation to i is at
    most ``tol_rad`` (sorted, unique; the relation is symmetric).

    R is (N, 3, 3). Exact up to floating-point rounding at the tolerance (the tree radius carries
    a 1e-9 relative margin): a boundary copy is stored when its |w| is within twice the chord
    radius of the class maximum, which holds every copy that can reach a point of the zone.
    """
    R = np.asarray(R, dtype=np.float64).reshape(-1, 3, 3)
    n = len(R)
    if n == 0:
        return []
    ops = _group_ops(symmetry_quats)
    q = _matrix_to_wxyz(R)
    eq = _qmul(q[:, None, :], ops[None, :, :])  # (N, K, 4): class members q S
    aw = np.abs(eq[..., 0])
    w_max = aw.max(axis=1)
    rep = eq[np.arange(n), aw.argmax(axis=1)]
    rep = np.where(rep[:, :1] < 0, -rep, rep)
    r = chord_radius(tol_rad)
    near = (w_max[:, None] - aw) <= 2.0 * r + 1e-12  # (N, K): copies kept in the tree
    idx_i, _ = np.nonzero(near)
    copies = eq[near]
    pts = np.concatenate([copies, -copies], axis=0)
    owner = np.concatenate([idx_i, idx_i])
    tree = cKDTree(pts)
    found = tree.query_ball_point(rep, r * (1.0 + 1e-9) + 1e-12)
    out: List[np.ndarray] = []
    for i, lst in enumerate(found):
        js = np.unique(owner[np.asarray(lst, dtype=np.int64)]) if len(lst) else np.empty(0, int)
        out.append(js[js != i])
    return out


def greedy_nms(
    R: np.ndarray,
    tol_rad: float,
    symmetry_quats: Optional[np.ndarray] = None,
    order: Optional[np.ndarray] = None,
    neighbours: Optional[List[np.ndarray]] = None,
) -> Tuple[np.ndarray, np.ndarray]:
    """Greedy non-maximum suppression in the given priority order (default: index order).

    An element is kept when no earlier kept element lies within ``tol_rad`` (the semantics of
    ``csl.dedup_orientations``); every other element is assigned to the first kept element that
    suppressed it. Returns (kept indices in priority order, leader[i] = the kept element
    that element i belongs to).
    """
    n = len(R)
    nb = neighbours if neighbours is not None else neighbour_lists(R, tol_rad, symmetry_quats)
    prio = np.arange(n) if order is None else np.asarray(order, dtype=np.int64)
    leader = np.full(n, -1, dtype=np.int64)
    kept: List[int] = []
    for i in prio:
        if leader[i] != -1:
            continue
        leader[i] = i
        kept.append(int(i))
        for j in nb[i]:
            if leader[j] == -1:
                leader[j] = i
    return np.array(kept, dtype=np.int64), leader


def component_labels(n: int, neighbours: List[np.ndarray]) -> np.ndarray:
    """Connected components of the neighbour graph (single-linkage groups), label per element."""
    rows = np.concatenate([np.full(len(nb), i, dtype=np.int64) for i, nb in enumerate(neighbours)])
    cols = np.concatenate([np.asarray(nb, dtype=np.int64) for nb in neighbours])
    g = coo_matrix((np.ones(len(rows)), (rows, cols)), shape=(n, n))
    return connected_components(g, directed=False)[1]


def merged_clique_labels(
    R: np.ndarray,
    clique: np.ndarray,
    link_rad: float,
    symmetry_quats: Optional[np.ndarray] = None,
    neighbours: Optional[List[np.ndarray]] = None,
) -> np.ndarray:
    """Label per candidate of its merged clique: the connected component of cliques, two cliques
    being linked when a candidate of one is within ``link_rad`` of a candidate of the other."""
    clique = np.asarray(clique, dtype=np.int64)
    nb = neighbours if neighbours is not None else neighbour_lists(R, link_rad, symmetry_quats)
    _, dense = np.unique(clique, return_inverse=True)
    n_cl = int(dense.max()) + 1 if len(dense) else 0
    rows, cols = [], []
    for i, js in enumerate(nb):
        for j in js:
            if dense[i] != dense[j]:
                rows.append(dense[i])
                cols.append(dense[j])
    g = coo_matrix((np.ones(len(rows)), (rows, cols)), shape=(n_cl, n_cl))
    comp = connected_components(g, directed=False)[1]
    return comp[dense]


# ---------------------------------------------------------------------------
# The opt-in specification (config key CandidateDedup) and the level-loop stages
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class DedupSpec:
    """What CandidateDedup does. Comma-separated, no spaces, e.g. ``clique=merged:1,pre=0.5``.

    clique=cpp:K|merged:K[@L]   before the quick MC keep the best K candidates (by discrete score)
                                of each clique: ``cpp`` = the candidates of one FZ point (C++
                                grouping), ``merged`` = connected components of cliques linked
                                when two candidates lie within L x diameter (default L = 1)
    pre=T                       before the quick MC drop candidates within T x diameter of a
                                better-scored one (symmetry-reduced), so the MC runs on distinct
                                inputs
    post=T                      after the quick MC, before the ranking cut: the kept quarter is
                                counted over candidates that are distinct within T x diameter
                                (a candidate within T x diameter of a better one is passed over
                                and ordered after the distinct ones)
    quarter=orig|filtered       the size of the kept quarter is counted from the discrete stage's
                                list (default; the slots freed by a filter go to the next
                                candidate) or from the list after the last filter that ran
    """

    clique: str = ""  # "", "cpp" or "merged"
    clique_k: int = 1
    clique_link: float = 1.0
    pre: float = 0.0
    post: float = 0.0
    quarter: str = "orig"

    @property
    def active(self) -> bool:
        return bool(self.clique) or self.pre > 0 or self.post > 0

    @classmethod
    def parse(cls, text: str) -> "DedupSpec":
        text = (text or "").strip()
        if text.lower() in ("", "off", "0", "none"):
            return cls()
        kw: dict = {}
        for part in text.split(","):
            if "=" not in part:
                raise ValueError(f"CandidateDedup: expected name=value, got {part!r}")
            name, val = (s.strip() for s in part.split("=", 1))
            if name == "clique":
                body, _, link = val.partition("@")
                kind, _, k = body.partition(":")
                if kind not in ("cpp", "merged"):
                    raise ValueError(f"CandidateDedup clique must be cpp or merged, got {kind!r}")
                kw["clique"] = kind
                kw["clique_k"] = int(k) if k else 1
                if kw["clique_k"] < 1:
                    raise ValueError("CandidateDedup clique K must be >= 1")
                if link:
                    kw["clique_link"] = float(link)
            elif name in ("pre", "post"):
                kw[name] = float(val)
                if kw[name] < 0:
                    raise ValueError(f"CandidateDedup {name} must be >= 0")
            elif name == "quarter":
                if val not in ("orig", "filtered"):
                    raise ValueError(
                        f"CandidateDedup quarter must be orig or filtered, got {val!r}"
                    )
                kw["quarter"] = val
            else:
                raise ValueError(f"CandidateDedup: unknown option {name!r}")
        return cls(**kw)

    def to_text(self) -> str:
        if not self.active:
            return "off"
        parts = []
        if self.clique:
            link = "" if self.clique_link == 1.0 else f"@{self.clique_link:g}"
            parts.append(f"clique={self.clique}:{self.clique_k}{link}")
        if self.pre > 0:
            parts.append(f"pre={self.pre:g}")
        if self.post > 0:
            parts.append(f"post={self.post:g}")
        if self.quarter != "orig":
            parts.append(f"quarter={self.quarter}")
        return ",".join(parts)


def best_k_per_group(group: np.ndarray, score: np.ndarray, k: int) -> np.ndarray:
    """Boolean mask: the ``k`` lowest-score members of each group (ties by index)."""
    order = np.argsort(np.asarray(score, dtype=np.float64), kind="stable")
    count: dict = {}
    keep = np.zeros(len(order), dtype=bool)
    for i in order:
        g = int(group[i])
        if count.get(g, 0) < k:
            keep[i] = True
            count[g] = count.get(g, 0) + 1
    return keep


def keep_before_quick_mc(
    spec: DedupSpec,
    R: np.ndarray,
    score: np.ndarray,
    clique: np.ndarray,
    diameter: float,
    symmetry_quats: Optional[np.ndarray],
) -> np.ndarray:
    """Boolean mask of the candidates the quick MC runs on (clique filter C0, then dedup C1).

    ``score`` is the discrete-stage cost (lower is better); ``clique`` the FZ point each
    candidate came from (-1 = unknown, e.g. an extra candidate: its own group)."""
    n = len(R)
    keep = np.ones(n, dtype=bool)
    if n == 0:
        return keep
    clique = np.asarray(clique, dtype=np.int64).copy()
    unknown = clique < 0
    clique[unknown] = clique.max() + 1 + np.arange(int(unknown.sum()))
    if spec.clique:
        group = clique
        if spec.clique == "merged":
            group = merged_clique_labels(R, clique, spec.clique_link * diameter, symmetry_quats)
        keep &= best_k_per_group(group, score, spec.clique_k)
    if spec.pre > 0:
        sel = np.nonzero(keep)[0]
        order = np.argsort(np.asarray(score)[sel], kind="stable")
        kept, _ = greedy_nms(R[sel], spec.pre * diameter, symmetry_quats, order=order)
        keep[:] = False
        keep[sel[kept]] = True
    return keep


def distinct_order(
    spec: DedupSpec,
    R_ranked: np.ndarray,
    diameter: float,
    symmetry_quats: Optional[np.ndarray],
) -> Tuple[np.ndarray, int]:
    """C2. ``R_ranked`` is in ranking order (best first). Returns (new order of the indices:
    the distinct candidates in rank order, then the passed-over ones in rank order, the number
    of distinct candidates)."""
    n = len(R_ranked)
    kept, _ = greedy_nms(R_ranked, spec.post * diameter, symmetry_quats)
    rest = np.setdiff1d(np.arange(n), kept, assume_unique=True)
    return np.concatenate([kept, rest]), int(len(kept))
