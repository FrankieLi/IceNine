"""AR5 step 27 (rule PX): ``max_pairs`` 50,000 in the indexer, passed through ``IndexParams``.

Checks: ``max_pairs`` reaches ``cast_votes``; a voxel under 50,000 pairs returns exactly what the
default (600,000) returns; a voxel over it keeps the true orientation among the candidates; the
hexagonal fixture runs; the arm wrapper passes the setting and logs the work; the pure judges.
"""

import sys
from pathlib import Path
from types import SimpleNamespace
from typing import Any, List

import numpy as np
import pytest
from scipy.spatial.transform import Rotation

from icenine import ray_indexing as RI
from icenine.symmetry import create_cubic_symmetry

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts" / "autoresearch" / "ar5"))
import px_rules as PR  # noqa: E402

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


def test_max_pairs_reaches_cast_votes(
    crystal: RI.CrystalModel, monkeypatch: pytest.MonkeyPatch
) -> None:
    seen: List[int] = []
    orig = RI.cast_votes

    def spy(g: Any, c: Any, p: Any) -> Any:
        seen.append(p.max_pairs)
        return orig(g, c, p)

    monkeypatch.setattr(RI, "cast_votes", spy)
    g = noisy(make_rays(rand_rot(1), crystal, 20, 1), 0.3, 2)
    RI.index_rays(g, crystal, RI.IndexParams(max_pairs=PR.NEW_MAX_PAIRS))
    assert seen == [50_000]


def test_under_50000_pairs_equals_base(crystal: RI.CrystalModel) -> None:
    """Few rays make few pairs: the subsample never triggers, so the output is bit-identical."""
    g = np.concatenate(
        [noisy(make_rays(rand_rot(2), crystal, 25, 3), 0.3, 4), clutter(crystal, 20, 5)]
    )
    a = RI.index_rays(g, crystal, RI.IndexParams())
    b = RI.index_rays(g, crystal, RI.IndexParams(max_pairs=PR.NEW_MAX_PAIRS))
    assert a.n_pairs < PR.NEW_MAX_PAIRS
    assert np.array_equal(a.R, b.R)
    assert np.array_equal(a.explained, b.explained)
    assert (a.n_pairs, a.n_votes) == (b.n_pairs, b.n_votes)


def test_over_50000_pairs_is_capped_and_keeps_truth(crystal: RI.CrystalModel) -> None:
    R = rand_rot(7)
    g = np.concatenate([noisy(make_rays(R, crystal, 45, 6), 0.3, 7), clutter(crystal, 400, 8)])
    a = RI.index_rays(g, crystal, RI.IndexParams())
    b = RI.index_rays(g, crystal, RI.IndexParams(max_pairs=PR.NEW_MAX_PAIRS))
    assert a.n_pairs > PR.NEW_MAX_PAIRS and b.n_pairs == PR.NEW_MAX_PAIRS
    assert b.n_votes < a.n_votes
    d = RI.misorientation_deg(b.R, np.broadcast_to(R, b.R.shape), crystal.sym)
    assert d.min() < 1.0  # the true orientation is among the candidates (clean fixture)


def test_hexagonal_runs_with_50000() -> None:
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
    R = rand_rot(41)
    g = np.concatenate([noisy(make_rays(R, cr, 25, 6), 0.3, 7), clutter(cr, 40, 8)])
    res = RI.index_rays(g, cr, RI.IndexParams(max_pairs=PR.NEW_MAX_PAIRS))
    assert len(res.R) >= 1
    big = np.concatenate([noisy(make_rays(R, cr, 40, 6), 0.3, 7), clutter(cr, 600, 8)])
    res2 = RI.index_rays(big, cr, RI.IndexParams(max_pairs=PR.NEW_MAX_PAIRS))
    assert res2.n_pairs <= PR.NEW_MAX_PAIRS and len(res2.R) >= 1


def test_wrapper_passes_setting_and_logs() -> None:
    seen: List[int] = []

    def fake(sc: Any, verts: Any, crystal: Any, p: Any, top_k: int = 20, keep: bool = False) -> Any:
        seen.append(p.max_pairs)
        return SimpleNamespace(n_rays=5, n_pairs=7, n_votes=11), np.arange(5)

    for arm, want in (("old", 600_000), ("new", 50_000)):
        log: List[Any] = []
        w = PR.make_index_voxel(fake, arm, log)
        res, att = w(None, None, None, RI.IndexParams(), 20)
        assert seen[-1] == want and len(att) == 5
        assert log[0]["n_pairs"] == 7 and log[0]["n_votes"] == 11 and log[0]["max_pairs"] == want
        assert PR.summed(log)["pairs"] == 7
    with pytest.raises(ValueError):
        PR.make_index_voxel(fake, "fast", [])


def test_judges() -> None:
    assert PR.judge_q1(0, 2000, 3)["ok"] and PR.judge_q1(10, 2000, 0)["ok"]
    assert not PR.judge_q1(11, 2000, 0)["ok"]
    assert PR.judge_c1(30, 30, 100, 100)["ok"] and not PR.judge_c1(50, 30, 100, 100)["ok"]
    assert PR.judge_c2(24_497_957)["ok"] and not PR.judge_c2(25_000_000)["ok"]
    assert PR.judge_c3(3, 0.001, 2)["ok"] and not PR.judge_c3(4, 0.0, 0)["ok"]
    assert PR.judge_t(0.5, 0.5)["ok"] and not PR.judge_t(0.49, 0.5)["ok"]
    assert PR.judge_h1(0.2, 5, 0.005)["ok"] and not PR.judge_h1(0.19, 0, 0.0)["ok"]
    assert PR.fewer(40, 100) == pytest.approx(0.6)
