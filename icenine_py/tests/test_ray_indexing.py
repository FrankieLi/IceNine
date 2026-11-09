"""Tests of icenine.ray_indexing: symmetry helpers, orbits, scattering vectors, and the indexer on
synthetic exact vectors, with clutter and with a second grain."""

from types import SimpleNamespace
from typing import List, Tuple

import numpy as np
import pytest
from scipy.spatial.transform import Rotation

from icenine import ray_indexing as RI
from icenine.symmetry import create_cubic_symmetry

A_LAT = 3.61
Q0 = 2.0 * np.pi / A_LAT


@pytest.fixture(scope="module")
def sym() -> np.ndarray:
    return RI.proper_rotations(create_cubic_symmetry(A_LAT))


def fcc_vectors(max_n2: int = 20) -> np.ndarray:
    """FCC reciprocal vectors (all indices even or all odd), both signs, length up to |hkl|^2."""
    out = []
    r = range(-4, 5)
    for h in r:
        for k in r:
            for l in r:
                n2 = h * h + k * k + l * l
                if n2 == 0 or n2 > max_n2:
                    continue
                if len({h % 2, k % 2, l % 2}) == 1:
                    out.append((h, k, l))
    return Q0 * np.array(out, dtype=np.float64)


@pytest.fixture(scope="module")
def crystal(sym: np.ndarray) -> RI.CrystalModel:
    return RI.build_crystal(fcc_vectors(), sym, q_max=8.0)


def rand_rot(seed: int) -> np.ndarray:
    return np.asarray(Rotation.random(random_state=seed).as_matrix())


def test_proper_rotations_cubic(sym: np.ndarray) -> None:
    assert sym.shape == (24, 3, 3)
    for s in sym:
        assert np.allclose(s @ s.T, np.eye(3), atol=1e-9)
        assert np.linalg.det(s) > 0


def test_reduce_and_misorientation(sym: np.ndarray) -> None:
    R = np.stack([rand_rot(i) for i in range(20)])
    red = RI.reduce_rotations(R, sym)
    # same orientation class: zero misorientation to the original
    assert np.allclose(RI.misorientation_deg(red, R, sym), 0.0, atol=1e-5)
    # reduced rotations are at most the cubic maximum (62.8 deg) from the identity
    ang = np.degrees(np.arccos(np.clip((np.trace(red, axis1=1, axis2=2) - 1) / 2, -1, 1)))
    assert ang.max() < 62.9
    # an equivalent rotation R S has zero misorientation to R
    Rs = np.einsum("ij,jk->ik", R[0], sym[7])
    assert RI.misorientation_deg(Rs[None], R[:1], sym)[0] < 1e-5
    # a real 3 degree rotation is reported as 3 degrees
    d = Rotation.from_rotvec(np.radians(3.0) * np.array([0, 0, 1.0])).as_matrix()
    assert abs(RI.misorientation_deg((d @ R[1])[None], R[1:2], sym)[0] - 3.0) < 1e-6


def test_orbits(crystal: RI.CrystalModel) -> None:
    # every vector belongs to exactly one orbit, and an orbit has a single length
    seen = np.concatenate(crystal.orbit_members)
    assert sorted(seen.tolist()) == list(range(len(crystal.h)))
    for m, q in zip(crystal.orbit_members, crystal.orbit_len):
        assert np.allclose(crystal.q[m], q)
    # {111}: 8 vectors; {200}: 6; {220}: 12 (cubic, proper rotations)
    sizes = {
        round(float(q / Q0), 4): len(m) for q, m in zip(crystal.orbit_len, crystal.orbit_members)
    }
    assert sizes[round(np.sqrt(3), 4)] == 8
    assert sizes[2.0] == 6
    assert sizes[round(np.sqrt(8), 4)] == 12


def test_scattering_vector_convention() -> None:
    from icenine.constants import KEV_OVER_HBAR_C_IN_ANG

    omega = np.radians(np.array([10.0, 40.0]))
    G = np.array([[-3.0, 1.0, 2.0], [-2.0, -2.5, 1.0]])
    energy = 64.351
    k = KEV_OVER_HBAR_C_IN_ANG * energy
    d = np.array([1.0, 0.0, 0.0]) + G / k  # not unit, the convention test only needs d - beam
    rays = SimpleNamespace(frame=np.array([0, 1]), dirn=d)
    fi = np.array([0, 1, 2])
    sc = SimpleNamespace(
        observed_rays=lambda: rays,
        energy=energy,
        beam=SimpleNamespace(numpy=lambda: np.array([1.0, 0.0, 0.0])),
        _frame_index=fi,
        n_frames=2,
        low=0.0,
        width=omega[0] * 2 / 1.0,
    )
    # frame centres: bins 0 -> (0 + 0.5) * width = omega[0]; frame 1 -> 1.5 * width
    sc.width = omega[0] / 0.5
    om = RI.frame_centre_omegas(sc)
    assert np.allclose(om, [0.5 * sc.width, 1.5 * sc.width])
    g = RI.scattering_vectors(sc)
    for i in range(2):
        c, s = np.cos(om[i]), np.sin(om[i])
        Rz = np.array([[c, -s, 0], [s, c, 0], [0, 0, 1.0]])
        assert np.allclose(Rz @ g[i], G[i])  # lab = Rz(omega) sample


def make_rays(R: np.ndarray, crystal: RI.CrystalModel, n: int, seed: int) -> np.ndarray:
    rng = np.random.default_rng(seed)
    idx = rng.choice(len(crystal.h), n, replace=False)
    return crystal.h[idx] @ R.T


def test_index_exact_vectors(crystal: RI.CrystalModel, sym: np.ndarray) -> None:
    for seed in range(3):
        R = rand_rot(100 + seed)
        g = make_rays(R, crystal, 25, seed)
        res = RI.index_rays(g, crystal, RI.IndexParams())
        assert len(res.R) >= 1
        assert RI.misorientation_deg(res.R[0], R, sym) < 0.05
        assert res.explained[0] == 25
        assert res.n_rays == 25


def test_index_with_noise(crystal: RI.CrystalModel, sym: np.ndarray) -> None:
    rng = np.random.default_rng(5)
    R = rand_rot(7)
    g = make_rays(R, crystal, 30, 1)
    # rotate each vector by up to 0.4 degrees about a random axis, and scale it by up to 0.2%
    ax = rng.normal(size=g.shape)
    ax /= np.linalg.norm(ax, axis=1, keepdims=True)
    rot = Rotation.from_rotvec(ax * np.radians(rng.uniform(0, 0.4, (len(g), 1))))
    g = rot.apply(g) * (1 + rng.uniform(-0.002, 0.002, (len(g), 1)))
    res = RI.index_rays(g, crystal, RI.IndexParams())
    assert RI.misorientation_deg(res.R[0], R, sym) < 0.3
    assert res.explained[0] >= 28


def clutter_vectors(crystal: RI.CrystalModel, n: int, seed: int) -> np.ndarray:
    """Random directions with ring-matched lengths (a coincidence-hard clutter)."""
    rng = np.random.default_rng(seed)
    d = rng.normal(size=(n, 3))
    d /= np.linalg.norm(d, axis=1, keepdims=True)
    return d * rng.choice(crystal.q, n)[:, None]


def test_index_with_clutter_and_second_grain(crystal: RI.CrystalModel, sym: np.ndarray) -> None:
    R1, R2 = rand_rot(11), rand_rot(12)
    g = np.concatenate(
        [
            make_rays(R1, crystal, 35, 2),
            make_rays(R2, crystal, 22, 3),
            clutter_vectors(crystal, 120, 4),
        ]
    )
    rng = np.random.default_rng(0)
    g = g[rng.permutation(len(g))]
    res = RI.index_rays(g, crystal, RI.IndexParams(), top_k=5)
    d1 = RI.misorientation_deg(res.R, R1[None], sym)
    d2 = RI.misorientation_deg(res.R, R2[None], sym)
    assert d1.min() < 0.3 and d2.min() < 0.3  # both grains are in the list
    assert d1[0] < 0.3  # the grain with more rays comes first
    assert res.explained[0] >= 35
    assert res.explained[0] > res.explained[int(np.argmin(d2))]


def test_too_few_rays(crystal: RI.CrystalModel) -> None:
    res = RI.index_rays(np.zeros((1, 3)) + Q0 * 2, crystal, RI.IndexParams())
    assert len(res.R) == 0
    res = RI.index_rays(np.zeros((0, 3)), crystal, RI.IndexParams())
    assert len(res.R) == 0


def test_symmetry_independent_result(crystal: RI.CrystalModel, sym: np.ndarray) -> None:
    """A symmetry-equivalent orientation gives the same vectors, so the same candidate class."""
    R = rand_rot(21)
    g = make_rays(R, crystal, 20, 9)
    R_eq = R @ sym[11]
    g_eq = make_rays(R_eq, crystal, 20, 9)
    a = RI.index_rays(g, crystal, RI.IndexParams())
    b = RI.index_rays(g_eq, crystal, RI.IndexParams())
    assert RI.misorientation_deg(a.R[0], R, sym) < 0.05
    assert RI.misorientation_deg(b.R[0], R_eq, sym) < 0.05


def hcp_vectors(a: float, c: float, q_max: float) -> Tuple[np.ndarray, RI.CrystalModel]:
    """All reciprocal-lattice vectors of a hexagonal lattice up to ``q_max`` (pymatgen Cartesian
    frame, as ``CrystalStructure``), and the CrystalModel built with the hexagonal symmetry."""
    from pymatgen.core.lattice import Lattice
    from pymatgen.core.structure import Structure

    from icenine.symmetry import CrystalSymmetry

    st = Structure(
        Lattice.hexagonal(a, c), ["Ti", "Ti"], [[1 / 3, 2 / 3, 0.25], [2 / 3, 1 / 3, 0.75]]
    )
    symmetry = CrystalSymmetry(st)
    # the lattice-basis operators of a hexagonal crystal are not rotations
    raw = [np.asarray(m, dtype=np.float64) for m in symmetry.get_rotation_matrices()]
    assert max(np.abs(m @ m.T - np.eye(3)).max() for m in raw) > 0.1
    rl = st.lattice.reciprocal_lattice
    r = range(-6, 7)
    hkl = np.array([(h, k, l) for h in r for k in r for l in r if (h, k, l) != (0, 0, 0)])
    vec = rl.get_cartesian_coords(hkl.astype(np.float64))
    vec = vec[np.linalg.norm(vec, axis=1) <= q_max]
    return vec, RI.build_crystal(vec, RI.proper_rotations(symmetry), q_max=q_max)


def test_hexagonal_operators_and_indexing() -> None:
    """Hexagonal (hcp Ti): the operators are Cartesian rotations, the vector list is closed under
    them, the {10-10} family has 6 members, and exact vectors index to the right orientation."""
    vec, cr = hcp_vectors(2.95, 4.68, 5.0)
    sym = cr.sym
    assert sym.shape == (12, 3, 3)  # point group 622 (proper part of 6/mmm)
    for s in sym:
        assert np.allclose(s @ s.T, np.eye(3), atol=1e-9)
        assert np.linalg.det(s) > 0
    # the shortest in-plane family {10-10}: 6 vectors of one length, perpendicular to c (z)
    q100 = 4.0 * np.pi / (np.sqrt(3.0) * 2.95)
    o = int(np.argmin(np.abs(cr.orbit_len - q100)))
    assert abs(cr.orbit_len[o] - q100) < 1e-9
    assert len(cr.orbit_members[o]) == 6
    assert np.allclose(cr.h[cr.orbit_members[o]][:, 2], 0.0)
    R = rand_rot(41)
    g = make_rays(R, cr, 25, 6)
    res = RI.index_rays(g, cr, RI.IndexParams())
    assert RI.misorientation_deg(res.R[0], R, sym) < 0.05
    assert res.explained[0] == 25


def test_match_rays(crystal: RI.CrystalModel) -> None:
    R = rand_rot(31)
    g = make_rays(R, crystal, 10, 4)
    j, ang, rel = RI.match_rays(R, g, crystal, 0.5, 0.01)
    assert (ang < 1e-6).all() and (rel < 1e-9).all() and (j >= 0).all()
    # a vector of a length no ring has finds nothing
    j, ang, rel = RI.match_rays(R, np.array([[0.3, 0.1, 0.0]]), crystal, 0.5, 0.01)
    assert j[0] == -1 and np.isinf(ang[0])
