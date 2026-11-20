"""AR5 step 26 (rule IX): the faster indexer returns exactly what the reference returns.

``index_rays_fast`` / ``index_voxel_fast`` against ``index_rays`` / ``index_voxel`` on fixtures:
under 2 rays, no votes, tied votes and explained rays, more than 600,000 pairs, cubic and
hexagonal symmetry, a voxel through a stand-in ``SmoothCost``; and the neighbour count of the votes
against the k-d tree's own count, with and without repeated votes. Outputs must be bit-identical.
"""

from types import SimpleNamespace
from typing import Any, Tuple

import numpy as np
import pytest
from scipy.spatial import cKDTree
from scipy.spatial.transform import Rotation

from icenine import ray_indexing as RI
from icenine.symmetry import create_cubic_symmetry

A_LAT = 3.61
Q0 = 2.0 * np.pi / A_LAT


def fcc_vectors(max_n2: int = 20) -> np.ndarray:
    out = []
    r = range(-4, 5)
    for h in r:
        for k in r:
            for ll in r:
                n2 = h * h + k * k + ll * ll
                if n2 == 0 or n2 > max_n2:
                    continue
                if len({h % 2, k % 2, ll % 2}) == 1:
                    out.append((h, k, ll))
    return Q0 * np.array(out, dtype=np.float64)


@pytest.fixture(scope="module")
def crystal() -> RI.CrystalModel:
    sym = RI.proper_rotations(create_cubic_symmetry(A_LAT))
    return RI.build_crystal(fcc_vectors(), sym, q_max=8.0)


def rand_rot(seed: int) -> np.ndarray:
    return np.asarray(Rotation.random(random_state=seed).as_matrix())


def make_rays(R: np.ndarray, crystal: RI.CrystalModel, n: int, seed: int) -> np.ndarray:
    rng = np.random.default_rng(seed)
    idx = rng.choice(len(crystal.h), n, replace=False)
    return crystal.h[idx] @ R.T


def clutter(crystal: RI.CrystalModel, n: int, seed: int) -> np.ndarray:
    rng = np.random.default_rng(seed)
    d = rng.normal(size=(n, 3))
    d /= np.linalg.norm(d, axis=1, keepdims=True)
    return d * rng.choice(crystal.q, n)[:, None]


def noisy(g: np.ndarray, deg: float, seed: int) -> np.ndarray:
    rng = np.random.default_rng(seed)
    ax = rng.normal(size=g.shape)
    ax /= np.linalg.norm(ax, axis=1, keepdims=True)
    rot = Rotation.from_rotvec(ax * np.radians(rng.uniform(0, deg, (len(g), 1))))
    return np.asarray(rot.apply(g) * (1 + rng.uniform(-0.002, 0.002, (len(g), 1))))


def same(a: RI.Candidates, b: RI.Candidates) -> None:
    """Bit-identical rotations, counts and order, and the same bookkeeping numbers."""
    assert a.R.shape == b.R.shape
    assert np.array_equal(a.R, b.R)
    assert np.array_equal(a.explained, b.explained)
    assert np.array_equal(a.votes, b.votes)
    assert (a.n_rays, a.n_entries, a.n_pairs, a.n_votes) == (
        b.n_rays,
        b.n_entries,
        b.n_pairs,
        b.n_votes,
    )


def both(g: np.ndarray, crystal: RI.CrystalModel, top_k: int = 20) -> Tuple[Any, Any]:
    p = RI.IndexParams()
    a = RI.index_rays(g, crystal, p, top_k)
    b = RI.index_rays_fast(g, crystal, p, top_k)
    same(a, b)
    return a, b


def test_under_two_rays(crystal: RI.CrystalModel) -> None:
    for g in (np.zeros((0, 3)), np.zeros((1, 3)) + Q0 * 2):
        a, b = both(g, crystal)
        assert len(a.R) == 0 and len(b.R) == 0


def test_no_votes(crystal: RI.CrystalModel) -> None:
    """Two rays of ring lengths that are nearly parallel: no pair is far enough apart to vote."""
    q = crystal.orbit_len[0]
    g = np.array([[q, 0.0, 0.0], [q * np.cos(0.02), q * np.sin(0.02), 0.0]])
    a, b = both(g, crystal)
    assert len(a.R) == 0 and a.n_votes == 0
    assert len(b.R) == 0
    # two rays with no ring length at all
    a, b = both(np.array([[0.3, 0.1, 0.0], [0.0, 0.4, 0.1]]), crystal)
    assert len(a.R) == 0


def test_exact_vectors_tied_votes(crystal: RI.CrystalModel) -> None:
    """Noise-free rays: many votes are equal up to a rounding error and clusters tie in votes and
    explained rays; a repeated ray makes more ties."""
    for seed in range(3):
        R = rand_rot(100 + seed)
        g = make_rays(R, crystal, 25, seed)
        a, _ = both(g, crystal)
        assert a.explained[0] == 25
        both(np.concatenate([g, g[:6]]), crystal)  # the same ray twice
    both(make_rays(rand_rot(5), crystal, 30, 3), crystal, top_k=3)  # a small top_k


def test_noise_clutter_and_second_grain(crystal: RI.CrystalModel) -> None:
    R1, R2 = rand_rot(11), rand_rot(12)
    g = np.concatenate(
        [
            noisy(make_rays(R1, crystal, 35, 2), 0.4, 1),
            noisy(make_rays(R2, crystal, 22, 3), 0.4, 2),
            clutter(crystal, 120, 4),
        ]
    )
    g = g[np.random.default_rng(0).permutation(len(g))]
    both(g, crystal)
    both(g, crystal, top_k=5)


def test_more_than_max_pairs(crystal: RI.CrystalModel) -> None:
    """A random subsample of pairs (the generator is seeded), identical in both."""
    g = np.concatenate(
        [noisy(make_rays(rand_rot(3), crystal, 40, 1), 0.3, 5), clutter(crystal, 1500, 6)]
    )
    a, _ = both(g, crystal)
    assert a.n_pairs == RI.IndexParams().max_pairs


def test_hexagonal_symmetry() -> None:
    from pymatgen.core.lattice import Lattice
    from pymatgen.core.structure import Structure

    from icenine.symmetry import CrystalSymmetry

    st = Structure(
        Lattice.hexagonal(2.95, 4.68), ["Ti", "Ti"], [[1 / 3, 2 / 3, 0.25], [2 / 3, 1 / 3, 0.75]]
    )
    symmetry = CrystalSymmetry(st)
    rl = st.lattice.reciprocal_lattice
    r = range(-6, 7)
    hkl = np.array([(h, k, ll) for h in r for k in r for ll in r if (h, k, ll) != (0, 0, 0)])
    vec = rl.get_cartesian_coords(hkl.astype(np.float64))
    vec = vec[np.linalg.norm(vec, axis=1) <= 5.0]
    cr = RI.build_crystal(vec, RI.proper_rotations(symmetry), q_max=5.0)
    g = np.concatenate([noisy(make_rays(rand_rot(41), cr, 25, 6), 0.3, 7), clutter(cr, 40, 8)])
    a, _ = both(g, cr)
    assert len(a.R) >= 1
    both(make_rays(rand_rot(42), cr, 20, 9), cr)  # exact vectors


def test_vote_counts_match_tree(crystal: RI.CrystalModel) -> None:
    rng = np.random.default_rng(0)
    base = rng.normal(size=(300, 4))
    base /= np.linalg.norm(base, axis=1, keepdims=True)
    chord = 0.3
    mult = rng.integers(1, 5, 300)
    q = np.repeat(base, mult, axis=0)
    q = q + rng.normal(scale=1e-16, size=q.shape)  # a rounding error apart
    q = q[rng.permutation(len(q))]
    want = np.asarray(cKDTree(q).query_ball_point(q, chord, return_length=True))
    assert np.array_equal(RI._vote_counts_fast(q, chord), want)
    # no repeated vote: the plain path
    q1 = base
    want1 = np.asarray(cKDTree(q1).query_ball_point(q1, chord, return_length=True))
    assert np.array_equal(RI._vote_counts_fast(q1, chord), want1)


def test_vote_counts_many_multiplicities(monkeypatch: pytest.MonkeyPatch) -> None:
    """More distinct multiplicities than the cap: one tree of all votes, same counts."""
    rng = np.random.default_rng(1)
    base = rng.normal(size=(40, 4))
    base /= np.linalg.norm(base, axis=1, keepdims=True)
    q = np.repeat(base, np.arange(1, 41), axis=0)
    want = np.asarray(cKDTree(q).query_ball_point(q, 0.5, return_length=True))
    assert len(np.unique(np.arange(1, 41))) > RI._MAX_MULT_CLASSES
    assert np.array_equal(RI._vote_counts_fast(q, 0.5), want)
    monkeypatch.setattr(RI, "_MAX_MULT_CLASSES", 100)
    assert np.array_equal(RI._vote_counts_fast(q, 0.5), want)


def test_index_voxel_fast_same_as_index_voxel(crystal: RI.CrystalModel) -> None:
    """Through a stand-in for ``SmoothCost``: the attached ray indices and candidates agree."""
    from icenine.constants import KEV_OVER_HBAR_C_IN_ANG

    R = rand_rot(77)
    g = np.concatenate([noisy(make_rays(R, crystal, 30, 1), 0.3, 3), clutter(crystal, 30, 2)])
    energy = 64.351
    k = KEV_OVER_HBAR_C_IN_ANG * energy
    om = 0.37
    c, s = np.cos(om), np.sin(om)
    G = g @ np.array([[c, -s, 0.0], [s, c, 0.0], [0.0, 0.0, 1.0]]).T  # Rz(omega) g
    beam = np.array([1.0, 0.0, 0.0])
    rays = SimpleNamespace(frame=np.zeros(len(g), dtype=int), dirn=beam + G / k)
    width = om / 0.5
    sc = SimpleNamespace(
        observed_rays=lambda: rays,
        energy=energy,
        beam=SimpleNamespace(numpy=lambda: beam),
        _frame_index=np.array([0, 1]),
        n_frames=1,
        low=0.0,
        width=width,
        attached_rays=lambda verts: np.arange(len(g))[::-1],
    )
    p = RI.IndexParams()
    a, att_a = RI.index_voxel(sc, None, crystal, p, 20)
    b, att_b = RI.index_voxel_fast(sc, None, crystal, p, 20)
    assert np.array_equal(att_a, att_b)
    same(a, b)
    assert len(a.R) >= 1
