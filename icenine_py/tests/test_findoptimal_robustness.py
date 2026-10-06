"""CSL tools and (later) feature / split helpers of scripts/findoptimal_robustness."""

import sys
from pathlib import Path

import numpy as np
from scipy.spatial.transform import Rotation

sys.path.insert(0, str(Path(__file__).parent.parent / "scripts" / "findoptimal_robustness"))

import csl  # noqa: E402

sys.path.insert(0, str(Path(__file__).parent))
sys.path.insert(0, str(Path(__file__).parent.parent / "scripts"))
sys.path.insert(0, str(Path(__file__).parent.parent / "benchmarks"))


def test_cubic_ops_group():
    assert csl.OPS.shape == (24, 3, 3)
    assert np.allclose(np.linalg.det(csl.OPS), 1.0)
    prod = csl.OPS[3] @ csl.OPS[7]
    assert any(np.allclose(prod, o, atol=1e-9) for o in csl.OPS)


def test_sigma3_relatives_are_60deg_about_111_and_four_distinct():
    R = Rotation.random(random_state=3).as_matrix()
    rel, labels = csl.csl_relatives(R, sigmas=[3])
    assert len(rel) == 4 and set(labels) == {"3"}
    for Rk in rel:
        M = R.T @ Rk  # crystal-frame misorientation, some symmetric variant
        # the reduced angle is 60 deg, and the axis (in the best symmetric variant) is <111>
        assert abs(csl.reduced_misorientation_deg(R, Rk) - 60.0) < 1e-6
        both = csl.OPS[:, None] @ M @ csl.OPS[None]
        rv = Rotation.from_matrix(both.reshape(-1, 3, 3)).as_rotvec()
        ang = np.degrees(np.linalg.norm(rv, axis=1))
        k = int(np.argmin(ang))
        ax = np.sort(np.abs(rv[k] / np.linalg.norm(rv[k])))
        assert np.allclose(ax, 1 / np.sqrt(3), atol=1e-6)
    # pairwise distinct modulo symmetry
    for i in range(4):
        for j in range(i):
            assert csl.reduced_misorientation_deg(rel[i], rel[j]) > 1.0


def test_relative_counts_and_classification():
    R = Rotation.random(random_state=5).as_matrix()
    counts = {}
    for e in csl.csl_table(11):
        rel, lab = csl.csl_relatives(R, sigmas=[e.sigma])
        counts[e.label] = sum(1 for x in lab if x == e.label)
        # every relative classifies as its own Sigma (or a lower one with the same rotation)
        c = csl.csl_classify(R, rel[0], max_sigma=11)
        assert c["sigma"] <= e.sigma and c["deviation"] < 0.05
    assert counts["3"] == 4
    assert counts["5"] == 6  # 36.87 deg <100>: 24 / |stabiliser| = 6 distinct orientations
    assert csl.csl_classify(R, R)["sigma"] == 1
    rand = Rotation.random(random_state=11).as_matrix()
    assert csl.csl_classify(R, rand)["angle"] > 1.0


def test_table_has_expected_sigmas():
    sig = {e.sigma for e in csl.csl_table(29)}
    assert sig == {3, 5, 7, 9, 11, 13, 15, 17, 19, 21, 23, 25, 27, 29}
    assert len(csl.csl_table(29)) == 21  # 3,5,7,9,11,13ab,15,17ab,19ab,21ab,23,25ab,27ab,29ab


def test_invariant_reflections_sigma3():
    R = np.eye(3)
    rel, _ = csl.csl_relatives(R, sigmas=[3])
    hk = np.array([[1, 1, 1], [1, 1, -1], [2, 0, 0], [1, 1, 0], [3, 1, 1]], dtype=float)
    inv = csl.invariant_reflection_mask(hk, R, rel)
    # twin about [111]: the (111) direction is shared with the relative whose axis is [111]
    assert inv[:, 0].any()
    assert not inv[:, 2].all()  # {200} is not preserved by every Sigma3 variant


def test_grain_disjoint_folds():
    import e2_models as M
    import common as C

    voxels = C.OUT_DIR / "voxels.npz"
    if not voxels.exists():
        import pytest

        pytest.skip("study voxel list not available")
    grain = np.load(voxels)["voxel_grain_id"]
    vpos = np.repeat(np.arange(len(grain)), 5)  # several candidates per voxel
    f = M.fold_of(vpos)
    assert set(f) == set(range(M.FOLDS))
    for v in range(len(grain)):  # a voxel is entirely in one fold
        assert len(set(f[vpos == v])) == 1
    for g in np.unique(grain):  # and a grain too
        assert len(set(f[np.isin(vpos, np.nonzero(grain == g)[0])])) == 1
    sizes = [len(set(vpos[f == k])) for k in range(M.FOLDS)]
    assert max(sizes) - min(sizes) <= 3
    assert np.array_equal(f, M.fold_of(vpos))  # deterministic


def test_feature_extractor_shape_and_determinism():
    import pytest

    from test_findoptimal_refactor import _build  # noqa: E402
    import features as F
    from optimizer_sweep import Geometry
    from icenine.reconstructor import _get_voxel_vertices

    rec, voxel, R_true = _build()
    lf = rec._make_local_cost_fn()
    phase = voxel.phase
    if phase not in lf._phase_recip_vecs:
        pytest.skip("no reflections for the voxel phase")
    geo = Geometry(2, 180, 2048, 2048)
    fe = F.FeatureExtractor(lf, geo, phase)
    rng = np.random.default_rng(0)
    keys = np.unique(rng.integers(0, 2 * 180 * 2048 * 2048, 200000))
    fe.set_image(keys)
    v = _get_voxel_vertices(voxel)
    f1 = fe.features(R_true, v, phase)
    f2 = fe.features(R_true, v, phase)
    assert f1.shape == (len(fe.feature_names()),)
    np.testing.assert_array_equal(f1, f2)
    assert np.isfinite(f1).all()
    # an empty image lights nothing: all hit features are zero
    fe.set_image(np.zeros(0, dtype=np.int64))
    f0 = fe.features(R_true, v, phase)
    names = fe.feature_names()
    assert f0[names.index("hit3")] == 0.0 and f0[names.index("log_n_pairs")] > 0
