"""scripts/coarse_proxy and the optional q_max of findoptimal_robustness.features."""

import os
import sys
from pathlib import Path
from types import ModuleType
from typing import Any, Tuple

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


def _import(name: str) -> ModuleType:
    """Import without leaking the thread settings the sweep scripts put into os.environ."""
    keep = dict(os.environ)
    try:
        return __import__(name)
    finally:
        for k in set(os.environ) - set(keep):
            del os.environ[k]


def env_voxel(C: ModuleType) -> int:
    return int(dict(np.load(C.OUT_DIR / "voxels.npz"))["voxel_indices"][0])


@pytest.fixture(scope="module")
def env() -> Tuple[Any, ...]:
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
def test_low_q_peak_table_is_the_restricted_full_table(env: Tuple[Any, ...], q: float) -> None:
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
    # the CSL-aware reflection masks are the full masks restricted to the kept reflections
    for sg in F.SIGMAS:
        assert np.array_equal(low.shared[sg], full.shared[sg][:, low.refl_full]), sg


@needs_cache
def test_q_max_none_leaves_features_unchanged(env: Tuple[Any, ...]) -> None:
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
        assert np.array_equal(a.features(e2["R"][i], vctx.vertices, phase), e2["X"][i])


@needs_cache
def test_low_q_costs_use_max_q_and_follow_the_images(env: Tuple[Any, ...]) -> None:
    """The two cost columns equal independently built cost functions with max_q = 5 (26
    reflections) at pixel radius 0 and 3, and follow the attached images."""
    from icenine.cost_functions import VoxelCostFunction

    C, F, W, vctx, keys, Rs = env
    phase, lf = vctx.voxel.phase, W.local_fn

    def reference(radius: int) -> Any:
        return VoxelCostFunction(
            simulator=lf.simulator, detector_list=lf.detector_list, range_map=lf.range_map,
            exp_data=lf.exp_data, sample=lf.sample, structure_list=lf.structure_list,
            mode="hard", eta_limit=lf.eta_limit, pixel_radius=radius, max_q=5.0,
            min_sin_eta=lf.min_sin_eta,
        )  # fmt: skip

    low = F.FeatureExtractor(W.local_fn, W.ctx.geo, phase, q_max=5.0)
    low.set_image(keys)
    ref0, ref3 = reference(0), reference(3)
    assert len(ref0._phase_recip_vecs[phase][1]) == len(low.g_mag) == 26
    R32 = Rs[1].astype(np.float32)
    f_clean = low.features(Rs[1], vctx.vertices, phase)
    assert f_clean[-2] == ref0.evaluate(R32, vctx.vertices, phase).cost
    assert f_clean[-1] == ref3.evaluate(R32, vctx.vertices, phase).cost
    # the realistic images of the same voxel: the cost changes and matches the reference again
    keys_all = np.load(C.CACHE_DIR / "images" / f"v{env_voxel(C)}_all.npz")["keys"]
    C.attach(keys_all)
    try:
        low.set_image(keys_all)
        ref0, ref3 = reference(0), reference(3)  # exp_data re-read at construction
        f_all = low.features(Rs[1], vctx.vertices, phase)
        assert f_all[-2] != f_clean[-2] or f_all[-1] != f_clean[-1]
        assert f_all[-2] == ref0.evaluate(R32, vctx.vertices, phase).cost
        assert f_all[-1] == ref3.evaluate(R32, vctx.vertices, phase).cost
    finally:
        C.attach(keys)


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


@needs_cache
@pytest.mark.parametrize("q_max", [5.0, None])
def test_features_batch_equals_per_candidate(env: Tuple[Any, ...], q_max: Any) -> None:
    """The batched pass returns exactly the per-candidate features (all columns, incl. the two
    costs) on the E2 candidates of the voxel, for several batch sizes, and counts the same cost
    evaluations."""
    C, F, W, vctx, keys, Rs = env
    phase, vert = vctx.voxel.phase, vctx.vertices
    fe = F.FeatureExtractor(W.local_fn, W.ctx.geo, phase, q_max=q_max)
    fe.set_image(keys)
    e2 = np.load(C.CACHE_DIR / "e2" / f"v{env_voxel(C)}_clean.npz")["R"]
    R = np.concatenate([np.stack(Rs), e2[:: max(1, len(e2) // 120)]])  # truth, random, E2 sample
    ref = np.stack([fe.features(r, vert, phase) for r in R])
    counters = [fn.eval_count for fn in fe._low_q_cost_fns()] if q_max else [W.local_fn.eval_count]
    for size in (1, 7, len(R)):
        got = np.concatenate(
            [fe.features_batch(R[i : i + size], vert, phase) for i in range(0, len(R), size)]
        )
        assert got.shape == ref.shape
        assert np.array_equal(got, ref), (size, np.abs(got - ref).max())
    # batch of 1 and the whole set add 2 evaluations per candidate to the counters
    after = [fn.eval_count for fn in fe._low_q_cost_fns()] if q_max else [W.local_fn.eval_count]
    n = 3 * len(R)  # three passes of the loop above
    assert [a - b for a, b in zip(after, counters)] == ([n, n] if q_max else [2 * n])
    # without costs: NaN cost columns, the same aggregates
    nc = fe.features_batch(R, vert, phase, with_cost=False)
    assert np.isnan(nc[:, -2:]).all() and np.array_equal(nc[:, :-2], ref[:, :-2])
    assert fe.features_batch(np.zeros((0, 3, 3)), vert, phase).shape == (0, ref.shape[1])


def test_batched_eq_decreases_with_batch_size() -> None:
    """keep_eighth.batched_eq: the proxy cost of a call is monotone in the batch size and the
    per-candidate cost falls with the batch size (T3 timing curve; skipped without it)."""
    if not (ROOT / "benchmarks" / "coarse_proxy" / "timing_batch.json").exists():
        pytest.skip("timing_batch.json not present")
    ke = _import("keep_eighth")
    one, fifty, big = (ke.batched_eq(np.array([n])) for n in (1, 50, 400))
    assert one < fifty < big
    assert one > fifty / 50 > 0 and fifty / 50 > big / 400 >= 0
    assert ke.batched_eq(np.array([], dtype=int)) == 0.0
