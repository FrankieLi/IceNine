"""Ray indexing: candidate orientations of a voxel from the observed rays (autoresearch AR1).

A *ray* is one diffracted beam, known from a spot seen on both detectors (``SmoothCost.observed_
rays``). Its unit lab direction ``d`` gives the scattering vector ``G = k (d - beam)`` (length
``2 k sin(theta)`` = the reciprocal-vector length ``q``). In the sample frame, ``g = Rz(-omega) G``
with omega the centre of the ray's frame (the forward model's own convention, ``lab = Rz(omega)
sample``): for a grain of orientation ``R`` the vector ``g`` equals ``R h`` for one reciprocal
vector ``h`` of the crystal.

Two rays of one grain fix its orientation: if ray 1 is of family F1 and ray 2 of family F2, and the
angle between ``g1`` and ``g2`` equals the angle between ``h1`` and some ``h2`` of F2, then the
rotation taking the triad (h1, h2) to (g1, g2) is a vote for the grain's orientation. A *family*
is an orbit of reciprocal vectors under the crystal's own proper rotations (nothing cubic-
specific: the operators come from the symmetry object, converted to the Cartesian frame of its
lattice; a hexagonal unit test covers the construction). Votes are reduced with the same operators
and clustered; each cluster is scored by the number of attached rays it explains and refined by a
least-squares (Procrustes) fit of the explained rays.

No cost function is evaluated here: a candidate is a proposal, to be checked with the strict cost.
"""

import time
from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
from scipy.spatial import cKDTree
from scipy.spatial.transform import Rotation

from .constants import KEV_OVER_HBAR_C_IN_ANG


# --------------------------------------------------------------------------------------------
# parameters and results
# --------------------------------------------------------------------------------------------
@dataclass
class IndexParams:
    """Tolerances (set from the measured errors of the sanity step, see the AR1 findings)."""

    q_max: float = 8.0  # reciprocal vectors up to this length (1/Angstrom)
    len_tol: float = 0.01  # relative |g| versus |h|
    ang_tol_deg: float = 0.7  # g versus R h, for explaining a ray
    pair_tol_deg: float = 1.0  # angle between two g versus between two h, for a vote
    min_sep_deg: float = 10.0  # skip ray pairs closer than this to (anti)parallel
    cluster_deg: float = 1.0  # votes closer than this (misorientation) are one cluster
    merge_deg: float = 1.0  # candidates closer than this are duplicates
    n_seeds: int = 150  # cluster seeds recounted with the symmetry operators
    max_pairs: int = 600_000  # random subsample of ray pairs above this
    refine_iters: int = 2
    seed: int = 0  # for the pair subsample only


@dataclass
class Candidates:
    """Result of one voxel: orientations sorted best first, with their scores."""

    R: np.ndarray  # (K, 3, 3) crystal-to-sample rotations, refined
    explained: np.ndarray  # (K,) number of attached rays explained after refinement
    votes: np.ndarray  # (K,) votes in the cluster
    n_rays: int = 0
    n_entries: int = 0
    n_pairs: int = 0
    n_votes: int = 0
    seconds: Dict[str, float] = field(default_factory=dict)
    vote_R: Optional[np.ndarray] = None  # all votes, kept only when asked (diagnostics)


# --------------------------------------------------------------------------------------------
# symmetry helpers (matrix based, any point group)
# --------------------------------------------------------------------------------------------
def cartesian_operators(symmetry: Any) -> List[np.ndarray]:
    """Point-group operators of a ``CrystalSymmetry`` in the Cartesian frame of its lattice.

    ``CrystalSymmetry.get_rotation_matrices()`` returns pymatgen's operators in the lattice
    (fractional) basis. They equal their Cartesian form only for orthogonal axes (cubic,
    tetragonal, orthorhombic); for a hexagonal lattice they are integer matrices that are not
    rotations. The Cartesian form is ``S_c = L^T S L^-T`` with ``L`` the lattice matrix (rows a, b,
    c), as in pymatgen's ``get_point_group_operations(cartesian=True)``. This is the frame of the
    reciprocal vectors of ``CrystalStructure`` (pymatgen Cartesian coordinates). Entries within
    1e-9 of an integer are set to that integer, so that orthogonal-axis operators stay exact.
    """
    mats = [np.asarray(m, dtype=np.float64) for m in symmetry.get_rotation_matrices()]
    structure = getattr(symmetry, "structure", None)
    if structure is None:
        return mats  # already Cartesian (no lattice to convert from)
    L = np.asarray(structure.lattice.matrix, dtype=np.float64)
    Lt, Lt_inv = L.T, np.linalg.inv(L.T)
    out: List[np.ndarray] = []
    for m in mats:
        c = Lt @ m @ Lt_inv
        near = np.abs(c - np.round(c)) < 1e-9
        c[near] = np.round(c[near])
        out.append(c)
    return out


def proper_rotations(symmetry: Any) -> np.ndarray:
    """(S, 3, 3) proper rotation matrices (det > 0) of a ``CrystalSymmetry`` in the Cartesian
    frame (``cartesian_operators``), duplicates removed. Raises ``ValueError`` if an operator is
    not orthonormal."""
    keep: List[np.ndarray] = []
    for m in cartesian_operators(symmetry):
        if np.abs(m @ m.T - np.eye(3)).max() > 1e-6:
            raise ValueError("a symmetry operator is not orthonormal in the Cartesian frame")
        if np.linalg.det(m) > 0 and not any(np.abs(m - k).max() < 1e-6 for k in keep):
            keep.append(m)
    return np.stack(keep)


def reduce_rotations(R: np.ndarray, sym: np.ndarray) -> np.ndarray:
    """Equivalent rotation ``R S`` of each ``R`` (m, 3, 3) closest to the identity (the largest
    trace over the crystal operators ``S``)."""
    m = R.reshape(-1, 9)
    tr = m @ sym.transpose(0, 2, 1).reshape(len(sym), 9).T  # tr(R S) = sum R_ij S_ji
    s = np.argmax(tr, axis=1)
    return np.einsum("mij,mjk->mik", R, sym[s])


def misorientation_deg(Ra: np.ndarray, Rb: np.ndarray, sym: np.ndarray) -> np.ndarray:
    """Symmetry-reduced angle (deg) between matching rotations ``Ra`` and ``Rb`` (broadcast
    (..., 3, 3)): the smallest rotation angle of ``Ra S Rb^T`` over the crystal operators."""
    M = np.einsum("...ji,...jk->...ik", Rb, Ra)  # Rb^T Ra
    flat = M.reshape(-1, 9)
    tr = flat @ sym.transpose(0, 2, 1).reshape(len(sym), 9).T  # tr(M S)
    c = np.clip((tr.max(axis=1) - 1.0) / 2.0, -1.0, 1.0)
    return np.degrees(np.arccos(c)).reshape(M.shape[:-2])


def _quat_pm(R: np.ndarray) -> np.ndarray:
    """Quaternions (x, y, z, w) with w >= 0 of rotation matrices (m, 3, 3)."""
    q = Rotation.from_matrix(R).as_quat()
    q[q[:, 3] < 0] *= -1.0
    return q


# --------------------------------------------------------------------------------------------
# the crystal's reciprocal vectors and their orbits
# --------------------------------------------------------------------------------------------
@dataclass
class CrystalModel:
    """Reciprocal vectors ``h`` (1/Angstrom, crystal frame) up to q_max and their orbits."""

    h: np.ndarray  # (K, 3)
    q: np.ndarray  # (K,)
    orbit: np.ndarray  # (K,) orbit label of each vector
    orbit_len: np.ndarray  # (O,) vector length of each orbit
    orbit_members: List[np.ndarray]  # indices into h, per orbit
    cos_table: Dict[Tuple[int, int], np.ndarray]  # (A, B) -> cos(angle(rep_A, h_b)), b in B
    sym: np.ndarray  # (S, 3, 3)

    @property
    def n_orbits(self) -> int:
        return len(self.orbit_len)


def build_crystal(h: np.ndarray, sym: np.ndarray, q_max: float = 8.0) -> CrystalModel:
    """Orbits of ``h`` under the proper rotations ``sym`` (a vector and its image are one orbit)."""
    h = np.asarray(h, dtype=np.float64)
    q = np.linalg.norm(h, axis=1)
    keep = (q <= q_max) & (q > 1e-9)
    h, q = h[keep], q[keep]
    tree = cKDTree(h)
    rows, cols = [], []
    for s in sym:
        d, j = tree.query(h @ s.T, k=1)
        if d.max() > 1e-4:
            raise ValueError("the reciprocal-vector list is not closed under the symmetry")
        rows.append(np.arange(len(h)))
        cols.append(j)
    n = len(h)
    adj = coo_matrix((np.ones(n * len(sym)), (np.concatenate(rows), np.concatenate(cols))), (n, n))
    n_orb, label = connected_components(adj, directed=False)
    members = [np.nonzero(label == o)[0] for o in range(n_orb)]
    orbit_len = np.array([q[m[0]] for m in members])
    hh = h / q[:, None]
    table: Dict[Tuple[int, int], np.ndarray] = {}
    for a in range(n_orb):
        rep = hh[members[a][0]]
        for b in range(n_orb):
            table[(a, b)] = hh[members[b]] @ rep
    return CrystalModel(h, q, label, orbit_len, members, table, sym)


# --------------------------------------------------------------------------------------------
# rays -> scattering vectors in the sample frame
# --------------------------------------------------------------------------------------------
def frame_centre_omegas(sc: Any) -> np.ndarray:
    """Omega (rad) at the centre of every frame of a ``SmoothCost`` (as ``attached_rays``)."""
    fi = sc._frame_index
    bins_c = np.array([np.nonzero(fi[:-1] == w)[0].mean() + 0.5 for w in range(sc.n_frames)])
    return np.asarray(sc.low + sc.width * bins_c, dtype=np.float64)


def scattering_vectors(sc: Any, idx: Optional[np.ndarray] = None) -> np.ndarray:
    """Sample-frame scattering vectors ``g = Rz(-omega) k (d - beam)`` (1/Angstrom) of the observed
    rays ``idx`` (default all); the beam is +x (SmoothCost's E4 frame check)."""
    rays = sc.observed_rays()
    if idx is None:
        idx = np.arange(len(rays.frame))
    k = KEV_OVER_HBAR_C_IN_ANG * sc.energy
    beam = np.asarray(sc.beam.numpy(), dtype=np.float64)
    G = k * (rays.dirn[idx] - beam)
    om = frame_centre_omegas(sc)[rays.frame[idx]]
    c, s = np.cos(om), np.sin(om)
    g = np.empty_like(G)
    g[:, 0] = c * G[:, 0] + s * G[:, 1]  # Rz(-omega)
    g[:, 1] = -s * G[:, 0] + c * G[:, 1]
    g[:, 2] = G[:, 2]
    return g


def crystal_from_smooth_cost(
    sc: Any, symmetry: Any, phase: int = 1, q_max: float = 8.0
) -> CrystalModel:
    """CrystalModel of a ``SmoothCost``'s phase (its reciprocal vectors) and a CrystalSymmetry."""
    h = sc._g[phase][0].numpy()
    return build_crystal(h, proper_rotations(symmetry), q_max)


# --------------------------------------------------------------------------------------------
# matching a rotation to rays
# --------------------------------------------------------------------------------------------
def match_rays(
    R: np.ndarray, g: np.ndarray, crystal: CrystalModel, ang_tol_deg: float, len_tol: float
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """For each ray, the reciprocal vector of ``R`` nearest in direction among those of the right
    length. Returns (index into crystal.h, angle in degrees, relative length error); the index is
    -1 and the angle is inf where no vector of the right length exists."""
    ph = crystal.h @ R.T  # (K, 3) R h
    gl = np.linalg.norm(g, axis=1)
    gh = g / gl[:, None]
    ok = np.abs(gl[:, None] / crystal.q[None, :] - 1.0) <= len_tol  # (n, K)
    cos = np.where(ok, gh @ (ph / crystal.q[:, None]).T, -2.0)
    j = np.argmax(cos, axis=1)
    cmax = cos[np.arange(len(g)), j]
    found = cmax > -1.5
    ang = np.where(found, np.degrees(np.arccos(np.clip(cmax, -1.0, 1.0))), np.inf)
    rel = np.where(found, np.abs(gl / crystal.q[j] - 1.0), np.inf)
    j = np.where(found, j, -1)
    return j, ang, rel


def explained_mask(
    R: np.ndarray, g: np.ndarray, crystal: CrystalModel, ang_tol_deg: float, len_tol: float
) -> np.ndarray:
    _, ang, rel = match_rays(R, g, crystal, ang_tol_deg, len_tol)
    return (ang <= ang_tol_deg) & (rel <= len_tol)


def refine_rotation(
    R: np.ndarray, g: np.ndarray, crystal: CrystalModel, p: IndexParams
) -> Tuple[np.ndarray, int]:
    """Least-squares (Procrustes) refinement of ``R`` on the rays it explains; returns the new
    rotation and the number of rays it explains."""
    for _ in range(p.refine_iters):
        j, ang, rel = match_rays(R, g, crystal, p.ang_tol_deg, p.len_tol)
        ok = (ang <= p.ang_tol_deg) & (rel <= p.len_tol)
        if ok.sum() < 2:
            break
        H = crystal.h[j[ok]]
        M = g[ok].T @ H  # sum g h^T
        U, _, Vt = np.linalg.svd(M)
        d = np.sign(np.linalg.det(U @ Vt))
        R = U @ np.diag([1.0, 1.0, d]) @ Vt
    return R, int(explained_mask(R, g, crystal, p.ang_tol_deg, p.len_tol).sum())


# --------------------------------------------------------------------------------------------
# voting
# --------------------------------------------------------------------------------------------
def _triad(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    """Orthonormal triads (m, 3, 3) (columns) from unit vectors a, b (m, 3)."""
    e1 = a
    e2 = b - (b * a).sum(1, keepdims=True) * a
    e2 /= np.linalg.norm(e2, axis=1, keepdims=True)
    e3 = np.cross(e1, e2)
    return np.stack([e1, e2, e3], axis=2)


def ring_entries(gl: np.ndarray, crystal: CrystalModel, len_tol: float) -> Tuple[np.ndarray, ...]:
    """(ray index, orbit index) of every ray-orbit pair whose lengths match within ``len_tol``."""
    rel = np.abs(gl[:, None] / crystal.orbit_len[None, :] - 1.0)
    r, o = np.nonzero(rel <= len_tol)
    return r, o


def cast_votes(g: np.ndarray, crystal: CrystalModel, p: IndexParams) -> Tuple[np.ndarray, int, int]:
    """Votes (m, 3, 3): rotations R with R h1 = g_a and R h2 = g_b for each ray pair (a, b) and
    family pair, h1 the representative of the first family and h2 each member of the second whose
    angle to h1 matches the angle between g_a and g_b. Returns (votes, n_entries, n_pairs)."""
    gl = np.linalg.norm(g, axis=1)
    gh = g / gl[:, None]
    er, eo = ring_entries(gl, crystal, p.len_tol)
    ne = len(er)
    if ne < 2:
        return np.zeros((0, 3, 3)), ne, 0
    ia, ib = np.triu_indices(ne, k=1)
    keep = er[ia] != er[ib]
    ia, ib = ia[keep], ib[keep]
    cosg = (gh[er[ia]] * gh[er[ib]]).sum(1)
    sep = np.cos(np.radians(p.min_sep_deg))
    keep = np.abs(cosg) < sep
    ia, ib, cosg = ia[keep], ib[keep], cosg[keep]
    if len(ia) > p.max_pairs:
        sel = np.random.default_rng(p.seed).choice(len(ia), p.max_pairs, replace=False)
        ia, ib, cosg = ia[sel], ib[sel], cosg[sel]
    n_pairs = len(ia)
    tol = np.radians(p.pair_tol_deg)
    oa, ob = eo[ia], eo[ib]
    hh = crystal.h / crystal.q[:, None]
    out: List[np.ndarray] = []
    key = oa * crystal.n_orbits + ob
    for kv in np.unique(key):
        sel = np.nonzero(key == kv)[0]
        A, B = int(kv // crystal.n_orbits), int(kv % crystal.n_orbits)
        tab = crystal.cos_table[(A, B)]  # (nb,)
        c = cosg[sel]
        width = tol * np.sqrt(np.maximum(1.0 - c * c, 1e-6))
        m = np.abs(tab[None, :] - c[:, None]) <= width[:, None]
        pi, mi = np.nonzero(m)
        if len(pi) == 0:
            continue
        pa, pb = ia[sel[pi]], ib[sel[pi]]
        h1 = np.broadcast_to(hh[crystal.orbit_members[A][0]], (len(pi), 3)).copy()
        h2 = hh[crystal.orbit_members[B][mi]]
        # a degenerate h2 (parallel to h1) cannot form a triad
        good = np.abs((h1 * h2).sum(1)) < 0.999
        if not good.any():
            continue
        pa, pb, h1, h2 = pa[good], pb[good], h1[good], h2[good]
        Fg = _triad(gh[er[pa]], gh[er[pb]])
        Fh = _triad(h1, h2)
        out.append(np.einsum("mij,mkj->mik", Fg, Fh))  # Fg Fh^T
    if not out:
        return np.zeros((0, 3, 3)), ne, n_pairs
    return np.concatenate(out), ne, n_pairs


def cluster_votes(
    votes: np.ndarray, sym: np.ndarray, p: IndexParams
) -> Tuple[np.ndarray, np.ndarray]:
    """Cluster rotations modulo the crystal symmetry. Returns (cluster rotations (C, 3, 3), vote
    counts (C,)) sorted by count, duplicates merged."""
    if len(votes) == 0:
        return np.zeros((0, 3, 3)), np.zeros(0, dtype=int)
    red = reduce_rotations(votes, sym)
    q = _quat_pm(red)
    tree = cKDTree(q)
    chord = 2.0 * np.sin(np.radians(p.cluster_deg) / 4.0)
    # coarse density at every vote (neighbours among the reduced representatives)
    cnt = np.asarray(tree.query_ball_point(q, chord, return_length=True), dtype=np.int64)
    seeds = np.argsort(-cnt, kind="stable")
    # symmetry-aware recount of the densest unassigned vote, whose members are then assigned
    assigned = np.zeros(len(q), dtype=bool)
    sel_R: List[np.ndarray] = []
    sel_n: List[int] = []
    for s in seeds:
        if len(sel_R) >= p.n_seeds:
            break
        if assigned[s]:
            continue
        Rs = red[s]
        eq = np.einsum("ij,sjk->sik", Rs, sym)  # all equivalents of the seed
        qe = _quat_pm(eq)
        pts = np.concatenate([qe, -qe])
        hits = set()
        for lst in tree.query_ball_point(pts, chord):
            hits.update(lst)
        idx = np.fromiter(hits, dtype=np.int64)
        assigned[idx] = True
        assigned[s] = True
        # mean of the members aligned to the seed (chordal L2 mean)
        mem = red[idx]
        tr = (
            np.einsum("ij,mik->mjk", Rs, mem).reshape(len(idx), 9)
            @ sym.transpose(0, 2, 1).reshape(len(sym), 9).T
        )
        al = np.einsum("mij,mjk->mik", mem, sym[np.argmax(tr, axis=1)])
        U, _, Vt = np.linalg.svd(al.sum(0))
        Rm = U @ np.diag([1.0, 1.0, np.sign(np.linalg.det(U @ Vt))]) @ Vt
        sel_R.append(Rm)
        sel_n.append(len(idx))
    order = np.argsort(-np.asarray(sel_n), kind="stable")
    return np.stack(sel_R)[order], np.asarray(sel_n)[order]


# --------------------------------------------------------------------------------------------
# the indexer
# --------------------------------------------------------------------------------------------
def index_rays(
    g: np.ndarray,
    crystal: CrystalModel,
    p: IndexParams,
    top_k: int = 20,
    keep_votes: bool = False,
) -> Candidates:
    """Candidate orientations from scattering vectors ``g`` (n, 3) in the sample frame."""
    t0 = time.perf_counter()
    n = len(g)
    if n < 2:
        return Candidates(np.zeros((0, 3, 3)), np.zeros(0, int), np.zeros(0, int), n_rays=n)
    votes, ne, npairs = cast_votes(g, crystal, p)
    t1 = time.perf_counter()
    cl_R, cl_n = cluster_votes(votes, crystal.sym, p)
    t2 = time.perf_counter()
    Rs: List[np.ndarray] = []
    ex: List[int] = []
    vn: List[int] = []
    for R, v in zip(cl_R, cl_n):
        Rr, e = refine_rotation(R, g, crystal, p)
        Rs.append(Rr)
        ex.append(e)
        vn.append(int(v))
    t3 = time.perf_counter()
    if not Rs:
        return Candidates(
            np.zeros((0, 3, 3)), np.zeros(0, int), np.zeros(0, int), n, ne, npairs, len(votes)
        )
    Rs_a, ex_a, vn_a = np.stack(Rs), np.asarray(ex), np.asarray(vn)
    order = np.lexsort((-vn_a, -ex_a))  # explained rays first, votes break ties
    keepR: List[int] = []
    for i in order:
        if len(keepR) >= top_k:
            break
        if (
            keepR
            and (
                misorientation_deg(Rs_a[i][None], Rs_a[keepR], crystal.sym) < 1.0 * p.merge_deg
            ).any()
        ):
            continue
        keepR.append(int(i))
    sel = np.asarray(keepR)
    t4 = time.perf_counter()
    return Candidates(
        Rs_a[sel],
        ex_a[sel],
        vn_a[sel],
        n,
        ne,
        npairs,
        len(votes),
        seconds=dict(vote=t1 - t0, cluster=t2 - t1, refine=t3 - t2, merge=t4 - t3, total=t4 - t0),
        vote_R=votes if keep_votes else None,
    )


def index_voxel(
    sc: Any,
    voxel_vertices: Any,
    crystal: CrystalModel,
    p: IndexParams,
    top_k: int = 20,
    keep_votes: bool = False,
) -> Tuple[Candidates, np.ndarray]:
    """Candidates of a voxel from the observed rays attached to it; also returns the attached ray
    indices. ``sc`` is a ``SmoothCost`` (E4 frame)."""
    att = sc.attached_rays(voxel_vertices)
    g = scattering_vectors(sc, att)
    res = index_rays(g, crystal, p, top_k, keep_votes)
    return res, att
