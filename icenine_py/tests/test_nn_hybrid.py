"""scripts/nn_hybrid: net stage vs the stored perturbation sweep, box sizing, and the finisher
reproducing the stored FindOptimal experiment B."""

import math
import sys
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


def test_covariance_box_sizing():
    import run as nh
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
    import run as nh

    cfg, sweep = nh.worker_args(["realistic_s0"])
    nh.init_worker(cfg)
    return nh, sweep


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
