"""scripts/nn_hybrid: net stage vs the stored perturbation sweep, box sizing, and the finisher
reproducing the stored FindOptimal experiment B."""

import math
import os
import sys
import time
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(ROOT / "scripts" / "nn_hybrid"))
sys.path.insert(0, str(ROOT / "scripts"))
sys.path.insert(0, str(ROOT / "benchmarks"))

SWEEP = ROOT / "benchmarks" / "toy_orientation_sweep"
needs_data = pytest.mark.skipif(
    not (SWEEP / "findoptimal_b_raw.npz").exists()
    or not (ROOT / "scripts" / "toy_orientation_sweep_model_realistic_s0.pt").exists(),
    reason="sweep raw files / models not present",
)


def _import(name: str):
    """Import a scripts/nn_hybrid module without leaking the OMP / MKL thread settings that
    the sweep scripts put in os.environ when imported."""
    keep = dict(os.environ)
    try:
        if name == "summary":
            # Load by path: other study dirs (coarse_proxy, ...) also have a `summary` module, and
            # a plain import returns whichever one another test put in sys.modules first.
            import importlib.util

            spec = importlib.util.spec_from_file_location(
                "_nn_hybrid_summary", ROOT / "scripts" / "nn_hybrid" / "summary.py"
            )
            assert spec is not None and spec.loader is not None
            mod = importlib.util.module_from_spec(spec)
            spec.loader.exec_module(mod)
            return mod
        return __import__(name)
    finally:
        for k in set(os.environ) - set(keep):
            del os.environ[k]


def test_covariance_box_sizing():
    nh = _import("run")
    from icenine.config_file import ConfigFile
    from icenine.orientation_search import SearchParameters

    L = np.diag([0.1, 0.4, 0.2])  # sigma_max = 0.4
    assert nh.sigma_max_deg(L) == pytest.approx(0.4)
    Q = np.linalg.qr(np.random.default_rng(0).normal(size=(3, 3)))[0]
    assert nh.sigma_max_deg(Q @ L) == pytest.approx(0.4)  # rotation of the factor: same spectrum
    assert np.isnan(nh.sigma_max_deg(np.full((3, 3), np.nan)))
    default = 0.3292
    assert nh.covariance_box_deg(0.05, default) == pytest.approx(default)  # clipped below
    assert nh.covariance_box_deg(0.4, default) == pytest.approx(1.2)  # 3 sigma
    assert nh.covariance_box_deg(5.0, default) == pytest.approx(2.0)  # clipped above
    # the diameter handed to refine_from_candidates gives exactly that search radius
    d = nh.box_to_diameter_rad(1.2)
    assert max(d / 3.0, math.radians(0.2)) == pytest.approx(math.radians(1.2))
    with pytest.raises(AssertionError):
        nh.box_to_diameter_rad(0.1)
    # the default box of the ReconstructQ8 parameters is the documented 0.329 deg
    params = SearchParameters.from_config(ConfigFile.from_file(str(nh.fs.RECON_CONFIG)))
    assert nh.default_box_deg(params) == pytest.approx(default, abs=1e-3)


@pytest.fixture(scope="module")
def worker():
    """One in-process worker (physics, net, timer patches), torn down afterwards."""
    import torch

    nh = _import("run")
    threads = torch.get_num_threads()
    cfg, sweep = nh.worker_args(["realistic_s0"])
    nh.init_worker(cfg)
    try:
        yield nh, sweep
    finally:
        nh._T.uninstall()
        torch.set_num_threads(threads)


@needs_data
def test_net_stage_reproduces_sweep_err_angle(worker):
    """The net passes (x1, x2, x3) of realistic_s0 on one voxel and radius give the stored sweep's
    err_angle within 1e-3 deg, for the realistic and the clean variant."""
    import perturbation_sweep as ps

    nh, sweep = worker
    ctx = nh.fs._W.ctx
    mi = [str(m) for m in sweep["models"]].index("realistic_s0")
    vpos, ri = 0, 6
    vidx = int(sweep["voxel_indices"][vpos])
    n = 0
    for vb in nh.variant_batches(
        ctx, vidx, vpos, ri, sweep["n_roi"][vpos, ri], sweep["fail_pass1"][vpos, ri], nh.VARIANTS
    ):
        R, _L = nh.net_stage(
            ctx, nh._NETS["realistic_s0"], vidx, vb.prep1, vb.b1, vb.ok1, vb.delta0, vb.R_nom0,
            vb.vctx.R_true, vb.draws, vb.vctx.sources, vb.variant, vb.seed,
        )  # fmt: skip
        ref = sweep["err_angle"][vpos, ri, :, vb.vi, mi, :]  # (D, P)
        for p in range(ref.shape[1]):
            ok = np.isfinite(ref[:, p])
            assert ok.sum() > 0
            got = np.linalg.norm(ps.estimate_error_deg(R[:, p][ok], vb.vctx.R_true), axis=-1)
            assert np.abs(got - ref[:, p][ok]).max() < 1e-3
        assert (np.isnan(ref[:, 2]) == np.isnan(R[:, 2, 0, 0])).all()
        n += 1
    assert n == 2


@needs_data
def test_findoptimal_alone_reproduces_stored_experiment_b(worker):
    """Our case images, seed and refine_from_candidates reproduce findoptimal_b_raw R_final exactly
    for 2 cases (one per variant)."""
    nh, sweep = worker
    ctx = nh.fs._W.ctx
    b = np.load(SWEEP / "findoptimal_b_raw.npz")
    vpos, ri = 0, 6
    vidx = int(sweep["voxel_indices"][vpos])
    bpos = list(b["voxel_indices"]).index(vidx)  # findoptimal_b_raw's voxel axis is sorted by index
    jset = {0: 1, 1: 2}  # variant index -> direction
    checked = 0
    for vb in nh.variant_batches(
        ctx, vidx, vpos, ri, sweep["n_roi"][vpos, ri], sweep["fail_pass1"][vpos, ri], nh.VARIANTS
    ):
        j = jset[vb.vi]
        assert vb.have[j]
        nh.fs.attach_images(nh.case_keys(ctx, vb, j))
        seed = nh.b_seed(ctx.args, vpos, ri, j, vb.vi)
        R, cost, _ = nh.refine_fo(vb.R_nom0[j], vb.vctx, seed)
        assert np.array_equal(R, b["R_final"][bpos, ri, j, vb.vi])
        assert cost == b["cost_final"][bpos, ri, j, vb.vi]
        checked += 1
    assert checked == 2


# ---------------------------------------------------------------------------
# summary.py helpers
# ---------------------------------------------------------------------------


def test_wilson_interval():
    sm = _import("summary")  # nn_hybrid's summary uses the shared stats.wilson (T2)
    w = sm.shared_stats.wilson
    lo, hi = w(0, 1000)
    assert lo == pytest.approx(0.0, abs=1e-12) and hi == pytest.approx(0.0038, abs=2e-4)
    lo, hi = w(50, 100)
    assert (lo, hi) == pytest.approx((0.404, 0.596), abs=2e-3)
    assert all(np.isnan(w(0, 0)))
    assert w(1000, 1000)[1] == 1.0  # clamped (the old local copy gave 1.0000000000000002)


def test_win_rate_tie_handling():
    sm = _import("summary")
    a = np.array([0.10, 0.20, 0.301, 0.5, np.nan])
    b = np.array([0.20, 0.10, 0.300, 0.5, 0.1])
    win, tie, n = sm.win_rate(a, b)  # win, loss, tie (|d| < 0.002), tie, ignored
    assert n == 4 and tie == pytest.approx(0.5)
    assert win == pytest.approx((1 + 0.5 * 2) / 4)
    assert np.isnan(sm.win_rate(np.array([np.nan]), np.array([1.0]))[0])


def test_reorder_permuted_and_partial_voxel_lists():
    sm = _import("summary")
    arr = np.array([[10.0], [20.0], [30.0]])
    src = np.array([5, 1, 9])
    out = sm.reorder(arr, src, np.array([9, 5, 1]))
    assert out[:, 0].tolist() == [30.0, 10.0, 20.0]
    part = sm.reorder(arr, src, np.array([1, 7]))  # voxel 7 absent: NaN
    assert part[0, 0] == 20.0 and np.isnan(part[1, 0])
    ints = sm.reorder(np.array([3, 4]), np.array([1, 2]), np.array([2, 8]))
    assert ints.tolist() == [4, 0]


# ---------------------------------------------------------------------------
# stage_timer.py
# ---------------------------------------------------------------------------


class _Cls:
    def method(self, x):
        return x + 1

    @staticmethod
    def smeth(x):
        return x * 2

    def evaluate(self, x):
        return x


class _Sub(_Cls):
    pass  # inherits method / evaluate


def test_stage_timer_nesting_and_exclusive_time():
    st = _import("stage_timer").StageTimer()
    with st.stage("outer"):
        time.sleep(0.02)
        with st.stage("inner"):
            time.sleep(0.03)
    assert st.calls["outer"] == 1 and st.calls["inner"] == 1
    assert st.inclusive["outer"] >= st.inclusive["inner"] >= 0.03
    assert st.exclusive["outer"] == pytest.approx(
        st.inclusive["outer"] - st.inclusive["inner"], abs=1e-6
    )
    assert st.exclusive["inner"] == pytest.approx(st.inclusive["inner"])
    before = st.snapshot()
    with st.stage("inner"):
        pass
    assert set(st.delta_since(before, "inclusive")) == {"inner"}


def test_stage_timer_patch_uninstall_restores_identity():
    st = _import("stage_timer").StageTimer()
    orig_method = _Cls.__dict__["method"]
    orig_static = _Cls.__dict__["smeth"]
    st.patch(_Cls, "method", "m")
    st.patch(_Cls, "smeth", "s")
    st.patch(_Sub, "method", "inherited")  # attribute inherited from _Cls
    assert _Cls().method(1) == 2 and _Cls.smeth(3) == 6 and _Sub().method(1) == 2
    assert isinstance(_Cls.__dict__["smeth"], staticmethod)
    # _Sub's wrapper wraps _Cls's already-patched wrapper, so "m" is entered by both calls
    assert st.calls["m"] == 2 and st.calls["s"] == 1 and st.calls["inherited"] == 1
    st.uninstall()
    assert _Cls.__dict__["method"] is orig_method
    assert _Cls.__dict__["smeth"] is orig_static
    assert "method" not in _Sub.__dict__  # the inherited attribute is deleted again


def test_stage_timer_refuses_double_install_and_uses_callable_names():
    st = _import("stage_timer").StageTimer()
    with st.installed():
        st.patch(_Cls, "method", lambda a, k: "big" if a[1] > 5 else "small")
        with pytest.raises(RuntimeError):
            st.patch(_Cls, "method", "again")
        _Cls().method(10)
        _Cls().method(1)
    assert st.calls["big"] == 1 and st.calls["small"] == 1
    assert _Cls().method(1) == 2 and st.calls["big"] == 1  # uninstalled by the context manager


def test_stage_timer_evaluation_counts_by_role_without_strong_refs():
    import gc
    import weakref

    st = _import("stage_timer").StageTimer()
    with st.installed():
        st.count_evaluations(_Cls)
        a, b = _Cls(), _Cls()
        st.label(a, "global")
        st.label(b, "local")
        before = st.snapshot()
        a.evaluate(1), b.evaluate(1), b.evaluate(2)
        assert st.evaluations(b) == 2
        assert st.evaluations_by_role() == {"global": 1, "local": 2}
        assert st.evaluations_since(before) == {"global": 1, "local": 2}
        st.reset()
        assert st.evaluations_by_role() == {"global": 0, "local": 0}
        ref = weakref.ref(a)
        del a
        gc.collect()
        assert ref() is None  # the timer holds no strong reference
