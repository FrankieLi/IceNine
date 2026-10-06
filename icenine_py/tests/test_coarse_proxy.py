"""scripts/coarse_proxy and the optional q_max of findoptimal_robustness.features."""

import os
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).parent.parent
FO = ROOT / "scripts" / "findoptimal_robustness"
sys.path.insert(0, str(ROOT / "scripts" / "coarse_proxy"))
sys.path.insert(0, str(FO))
sys.path.insert(0, str(ROOT / "scripts"))
sys.path.insert(0, str(ROOT / "benchmarks"))

needs_cache = pytest.mark.skipif(
    not (FO / "cache" / "e2").exists()
    or not (ROOT / "benchmarks" / "findoptimal_robustness" / "voxels.npz").exists(),
    reason="findoptimal_robustness caches / voxel set not present",
)


def _import(name: str):
    """Import without leaking the thread settings the sweep scripts put into os.environ."""
    keep = dict(os.environ)
    try:
        return __import__(name)
    finally:
        for k in set(os.environ) - set(keep):
            del os.environ[k]


def env_voxel(C) -> int:
    return int(dict(np.load(C.OUT_DIR / "voxels.npz"))["voxel_indices"][0])


@pytest.fixture(scope="module")
def env():
    C = _import("common")
    F = _import("features")
    from optimizer_sweep import voxel_context

    C.init_worker(C.worker_args())
    W = C.get_worker()
    v = env_voxel(C)
    keys = np.load(C.CACHE_DIR / "images" / f"v{v}_clean.npz")["keys"]
    C.attach(keys)
    vctx = voxel_context(W.ctx, v)
    d = np.load(C.CACHE_DIR / "e0" / f"v{v}_clean.npz")
    rng = np.random.default_rng(0)
    from scipy.spatial.transform import Rotation

    # the truth, a perturbed truth and two random orientations
    R = [d["R_true"]]
    R.append(Rotation.from_rotvec(np.radians(0.3) * rng.normal(size=3)).as_matrix() @ d["R_true"])
    R += list(Rotation.random(2, random_state=1).as_matrix())
    return C, F, W, vctx, keys, R


@needs_cache
@pytest.mark.parametrize("q", [4.0, 5.0])
def test_low_q_peak_table_is_the_restricted_full_table(env, q):
    C, F, W, vctx, keys, Rs = env
    phase = vctx.voxel.phase
    full = F.FeatureExtractor(W.local_fn, W.ctx.geo, phase)
    low = F.FeatureExtractor(W.local_fn, W.ctx.geo, phase, q_max=q)
    full.set_image(keys)
    low.set_image(keys)
    assert (low.g_mag.numpy() <= q).all() and len(low.g_mag) < len(full.g_mag)
    # the plan's families: Q4 = {111}, {200}; Q5 adds {220} (Cu); nothing below Q = 3
    assert np.round(full.q_levels, 3)[: low.n_fam].tolist() == low.q_levels.tolist()
    assert low.n_fam == {4.0: 2, 5.0: 3}[q]
    with pytest.raises(ValueError):
        F.FeatureExtractor(W.local_fn, W.ctx.geo, phase, q_max=3.0)
    n_checked = 0
    for R in Rs:
        tf = full.peak_table(R, vctx.vertices)
        tl = low.peak_table(R, vctx.vertices)
        m = full.fam[tf["refl"]] < low.n_fam  # fam <= f(Q)
        assert np.array_equal(low.refl_full[tl["refl"]], tf["refl"][m])
        for k in ("det", "hit0", "hit_any", "hit1", "hit3"):
            assert np.array_equal(tl[k], tf[k][m]), k
        n_checked += int(m.sum())
        # the family aggregates of the shared families agree (hit rates; `frac` has another base)
        ff = full.features(R, vctx.vertices, phase, with_cost=False)
        fl = low.features(R, vctx.vertices, phase, with_cost=False)
        for j in range(low.n_fam):
            for off, k in enumerate(("hit0", "hit3")):
                assert fl[5 + 3 * j + off] == pytest.approx(ff[5 + 3 * j + off])
    assert n_checked > 0


@needs_cache
def test_q_max_none_leaves_features_unchanged(env):
    C, F, W, vctx, keys, Rs = env
    phase = vctx.voxel.phase
    a = F.FeatureExtractor(W.local_fn, W.ctx.geo, phase)
    b = F.FeatureExtractor(W.local_fn, W.ctx.geo, phase, q_max=None)
    a.set_image(keys)
    b.set_image(keys)
    assert a.refl_full.tolist() == list(range(len(a.g_mag))) and b.n_fam == a.n_fam == 8
    for R in Rs[:2]:
        fa, fb = a.features(R, vctx.vertices, phase), b.features(R, vctx.vertices, phase)
        assert np.array_equal(fa, fb)
    # and they equal the features stored in the E2 dataset (built before q_max existed)
    e2 = np.load(C.CACHE_DIR / "e2" / f"v{env_voxel(C)}_clean.npz")
    for i in (0, 5, 300, len(e2["R"]) - 1):
        assert np.allclose(a.features(e2["R"][i], vctx.vertices, phase), e2["X"][i], atol=1e-9)


@needs_cache
def test_low_q_costs_use_max_q_and_follow_the_images(env):
    C, F, W, vctx, keys, Rs = env
    phase = vctx.voxel.phase
    low = F.FeatureExtractor(W.local_fn, W.ctx.geo, phase, q_max=5.0)
    low.set_image(keys)
    f = low.features(Rs[0], vctx.vertices, phase)
    c0, c3 = low._low_q_cost_fns()
    assert len(c0._phase_recip_vecs[phase][1]) == len(low.g_mag) == 26
    assert f[-2] == pytest.approx(c0.evaluate(Rs[0].astype(np.float32), vctx.vertices, phase).cost)
    assert c0.pixel_radius == 0 and c3.pixel_radius == 3
    assert 0.0 <= f[-2] <= 1.0 and f[-1] <= f[-2] + 1.0


def test_folds_are_disjoint_by_grain_and_voxel():
    """Every grain, hence every voxel and both variants of it, lies in exactly one fold, and all
    four folds are used (the proxy's models are always applied to unseen grains)."""
    C = _import("common")
    import e2_models as M

    path = C.OUT_DIR / "voxels.npz"
    if not path.exists():
        pytest.skip("voxel set not present")
    grain = np.load(path)["voxel_grain_id"][: C.N_VOXELS]
    vpos = np.arange(len(grain))
    fold = M.fold_of(vpos)
    assert set(fold.tolist()) == set(range(M.FOLDS))
    for g in np.unique(grain):
        assert len(set(fold[grain == g].tolist())) == 1
    # a candidate array built like the dataset (two variants per voxel) inherits the voxel's fold
    vp2 = np.repeat(vpos, 2)
    assert np.array_equal(M.fold_of(vp2), np.repeat(fold, 2))


def test_reduced_angle_matrix_equals_csl_reference():
    import labels as L

    csl = _import("csl")
    from scipy.spatial.transform import Rotation

    Ra = Rotation.random(7, random_state=2).as_matrix()
    Rb = Rotation.random(5, random_state=3).as_matrix()
    A = L.reduced_angle_matrix(Ra, Rb)
    for i in range(len(Ra)):
        for j in range(len(Rb)):
            assert A[i, j] == pytest.approx(float(csl.reduced_misorientation_deg(Ra[i], Rb[j])))
    # symmetry-equivalent copies are at zero distance
    assert L.reduced_angle_matrix(Ra[:1] @ csl.OPS[5], Ra[:1])[0, 0] < 1e-3  # arccos near 0


def test_label_priority_truth_then_relative_then_own_cost():
    import labels as L

    csl = _import("csl")
    R_true = np.eye(3)
    rel, _ = csl.csl_relatives(R_true, sigmas=[3])
    from scipy.spatial.transform import Rotation

    tilt = Rotation.from_rotvec([0.0, 0.0, np.radians(1.0)]).as_matrix()
    far = Rotation.from_rotvec(np.radians(40.0) * np.array([1.0, 0, 0])).as_matrix()
    R = np.stack([tilt, rel[1] @ tilt, far, far])
    e2 = dict(
        R=R, err=np.array([1.0, 59.0, 40.0, 40.0]), source=np.array([0, 0, 0, 2]),
        cost=np.array([0.5, 0.6, 0.7, np.nan]),
    )  # fmt: skip
    lab = dict(rel_R0=rel, rel_post_cost=np.array([0.1, 0.2, 0.3, 0.4]), truth_full_cost=0.05)
    y, cat, near = L.build_labels(e2, lab)
    assert cat.tolist() == [0, 1, 2, 3]
    assert y[:3].tolist() == [0.05, 0.2, 0.7] and np.isnan(y[3])
    assert near[1] == pytest.approx(1.0, abs=1e-6)


def test_pruning_recall_levels_and_fractions():
    """recall(): the best basin candidate of a group must be among the max(1, int(n * frac)) best
    scored; groups without a basin candidate are skipped."""
    md = _import("models")
    err = np.array([10.0, 2.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0, 5.0, 6.0])
    D = dict(err=err)
    grp = [np.arange(8), np.arange(8, 10)]  # the second group has no candidate below 3 deg
    thr = np.array([3.0, 3.0])
    good = np.array([0, 9, 1, 2, 3, 4, 5, 6, 0, 0], float)  # basin candidate (index 1) scored best
    bad = np.array([9, 0, 8, 7, 6, 5, 4, 3, 0, 0], float)  # basin candidate ranked last
    assert md.recall(D, good, grp, 0.25, thr) == (1.0, 1)
    assert md.recall(D, bad, grp, 0.25, thr) == (0.0, 1)
    # frac 1/8 of 8 candidates keeps 1: the basin candidate at rank 1 (second best) is lost
    second = np.array([5, 4, 1, 2, 3, 0, 0, 0, 0, 0], float)
    assert md.recall(D, second, grp, 0.25, thr)[0] == 1.0
    assert md.recall(D, second, grp, 0.125, thr)[0] == 0.0
    # level-3 style threshold of 1 deg: a 2 deg candidate no longer counts
    assert md.recall(D, good, grp, 0.25, np.array([1.0, 1.0]))[1] == 0


def test_feature_sets_select_low_q_columns():
    md = _import("models")
    n = 5
    X = np.arange(n * 62, dtype=float).reshape(n, 62)
    D = dict(X=X, X4=np.zeros((n, 44)), X5=np.ones((n, 47)))
    assert md.feature_matrix("e2full", D).shape == (n, 62)
    c4 = md.feature_matrix("cache4", D)  # log_n_pairs + families 0-1 (3 columns each)
    assert c4.shape == (n, 7) and np.array_equal(c4[:, 0], X[:, 0])
    assert np.array_equal(c4[:, 1:], X[:, 5:11])
    assert md.feature_matrix("cache5", D).shape == (n, 10)
    l5c8 = md.feature_matrix("lowq5+c8", D)  # F-lowQ plus the free Q8 local cost (column -2)
    assert l5c8.shape == (n, 48) and np.array_equal(l5c8[:, -1], X[:, -2])
    assert md.feature_matrix("lowq4", D).shape == (n, 44)
    hand = md.untrained_score("hand", dict(X=X, X4=D["X4"], X5=D["X5"]))
    assert np.allclose(hand, X[:, 1:5].mean(axis=1))
    assert np.array_equal(md.untrained_score("cost8", D), -X[:, -2])
