"""Helpers of scripts/perturbation_sweep.py that need no physics: voxel sampling and exclusion,
the perturbation convention (target angle exactly r), re-centring, error metrics, summary."""

import sys
from pathlib import Path

import numpy as np
import pytest
import torch
from scipy.spatial.transform import Rotation

sys.path.insert(0, str(Path(__file__).parent.parent / "scripts"))
import perturbation_sweep as ps  # noqa: E402

from icenine.orientation_eval import (  # noqa: E402
    error_summary,
    offsets_to_matrices,
)


def _rand_R(seed: int) -> np.ndarray:
    return Rotation.random(random_state=seed).as_matrix()


def test_excluded_voxel_indices_reads_dataset_files(tmp_path):
    a, b = tmp_path / "a.pt", tmp_path / "b.pt"
    torch.save({"voxel_indices": torch.tensor([5, 3, 9])}, a)
    torch.save({"voxel_indices": torch.tensor([3, 11])}, b)
    assert ps.excluded_voxel_indices([str(a), str(b)]).tolist() == [3, 5, 9, 11]


def test_eligible_candidates_apply_rmax_and_exclusion():
    r = np.array([0.0, 100.0, 499.0, 501.0, 300.0, 20.0])
    assert ps.eligible_candidates(r, [1, 4], 500.0).tolist() == [0, 2, 5]


def test_sample_voxels_is_seeded_and_never_returns_excluded():
    r = np.linspace(0, 400, 200)
    excluded = np.arange(0, 200, 7)
    cand = ps.eligible_candidates(r, excluded, 500.0)
    a = ps.sample_voxels(cand, 25, 3, lambda i: True)
    assert a == ps.sample_voxels(cand, 25, 3, lambda i: True)
    assert a != ps.sample_voxels(cand, 25, 4, lambda i: True)
    assert len(a) == len(set(a)) == 25
    assert not set(a) & set(excluded.tolist())


def test_sample_voxels_replaces_unusable_in_the_same_order():
    cand = np.arange(100)
    bad = {int(x) for x in np.random.default_rng(0).permutation(100)[:5]}
    plain = ps.sample_voxels(cand, 10, 1, lambda i: True)
    skip = ps.sample_voxels(cand, 10, 1, lambda i: i not in bad)
    assert not set(skip) & bad
    assert skip == [i for i in ps.sample_voxels(cand, 15, 1, lambda i: True) if i not in bad][:10]
    assert len(plain) == len(skip) == 10


@pytest.mark.parametrize("r", [0.05, 0.5, 1.0, 5.0])
def test_perturbation_angle_is_exactly_r(r):
    rng = np.random.default_rng(0)
    R_true = _rand_R(1)
    delta = ps.random_rotvecs(50, r, rng)
    assert np.allclose(np.linalg.norm(delta, axis=1), r, rtol=0, atol=1e-12)
    R_nom = ps.perturbed_nominal(R_true, delta)
    # the truth relative to the nominal is delta, so its angle is r ...
    rel = ps.relative_offset_deg(R_true, R_nom)
    assert np.allclose(rel, delta, atol=1e-10)
    assert np.allclose(np.linalg.norm(rel, axis=1), r, atol=1e-10)
    # ... in the dataset convention exp([delta]x) R_nom = R_true (what the renderer draws)
    assert np.allclose(offsets_to_matrices(delta, R_nom), R_true, atol=1e-12)
    # and the geodesic distance between the two orientations is r
    ang = np.degrees(
        np.arccos(
            np.clip((np.trace(R_true @ np.swapaxes(R_nom, 1, 2), axis1=1, axis2=2) - 1) / 2, -1, 1)
        )
    )
    assert np.allclose(ang, r, atol=1e-6)


def test_random_axes_are_isotropic():
    d = ps.random_rotvecs(20000, 1.0, np.random.default_rng(2))
    assert np.abs(d.mean(0)).max() < 0.02
    assert np.allclose((d**2).mean(0), 1 / 3, atol=0.02)


def test_recentring_keeps_the_same_truth():
    rng = np.random.default_rng(4)
    R_true = _rand_R(7)
    delta = ps.random_rotvecs(10, 2.0, rng)
    R_nom = ps.perturbed_nominal(R_true, delta)
    for _ in range(3):  # three passes with arbitrary (wrong) estimates
        delta_hat = delta * 0.3 + rng.normal(0, 0.2, size=delta.shape)
        R_nom, delta = ps.recentre(R_nom, delta_hat, R_true)
        # the re-rendered orientation exp([delta]x) R_nom is still exactly the true orientation
        assert np.allclose(offsets_to_matrices(delta, R_nom), R_true, atol=1e-12)
    # a perfect estimate leaves a zero offset and zero error
    R_new, d_new = ps.recentre(R_nom, delta, R_true)
    assert np.allclose(d_new, 0, atol=1e-9)
    assert np.allclose(ps.estimate_error_deg(R_new, R_true), 0, atol=1e-9)


def test_estimate_error_matches_error_summary():
    rng = np.random.default_rng(5)
    R_nom = _rand_R(8)
    truth = rng.normal(size=(30, 3))
    pred = truth + rng.normal(0, 0.05, size=(30, 3))
    R_est = offsets_to_matrices(pred, np.broadcast_to(R_nom, (30, 3, 3)))
    R_true = offsets_to_matrices(truth, np.broadcast_to(R_nom, (30, 3, 3)))
    ang = np.linalg.norm(ps.estimate_error_deg(R_est, R_true), axis=1)
    assert np.isclose(np.median(ang), error_summary(pred, truth)["median_angle"], atol=1e-9)


def test_mahalanobis_sq():
    chol = np.diag([0.1, 0.2, 0.4])[None]
    truth, pred = np.array([[0.1, 0.2, 0.4]]), np.zeros((1, 3))
    assert np.allclose(ps.mahalanobis_sq(truth, pred, chol), 3.0)


def test_grain_ids_and_boundary_flag():
    R1, R2 = _rand_R(1), _rand_R(2)
    g = ps.grain_ids(np.stack([R1, R1, R2, R2, R1]))
    assert g.tolist() == [g[0], g[0], g[2], g[2], g[0]] and g[0] != g[2]
    pos = np.array([[0, 0], [5, 0], [10, 0], [100, 0], [105, 0.0]])
    gid = np.array([0, 0, 1, 1, 1])
    assert ps.near_boundary(pos, gid, 6.0).tolist() == [False, True, True, False, False]


def test_load_model_roundtrip(tmp_path):
    from icenine.toy_orientation_model import GNLayerNet

    kw = dict(
        window_size=8,
        in_channels=2,
        context_dim=16,
        n_iter=2,
        frame_half_width=4,
        frame_width_rad=0.0175,
        pairing=False,
    )
    net = GNLayerNet(**kw)
    with torch.no_grad():
        for p in net.parameters():
            p.add_(0.01 * torch.randn_like(p))
    path = tmp_path / "m.pt"
    torch.save({"arch": "gn", "model_kwargs": kw, "state_dict": net.state_dict()}, path)
    net2, _ = ps.load_model(str(path))
    x = torch.rand(2, 5, 2, 8, 8)
    ctx, aux = torch.randn(2, 5, 16), {"nom_off": torch.zeros(2, 5, 3)}
    with torch.no_grad():
        a, b = net.eval()(x, ctx, aux), net2(x, ctx, aux)
    assert torch.allclose(a[0], b[0]) and torch.allclose(a[1], b[1])


def test_summarize_counts_failures_and_metrics():
    V, R, D, Va, Mo, P = 6, 2, 4, 2, 2, 3
    rng = np.random.default_rng(0)
    shape = (V, R, D, Va, Mo, P)
    fail = np.zeros((V, R, D, Va), dtype=np.int8)
    fail[0, 1, :2, :] = 4  # two failed cases at radius index 1 for voxel 0
    err = np.abs(rng.normal(0.05, 0.01, size=shape)).astype(np.float32)
    raw = dict(
        radii=np.array([0.5, 5.0]),
        models=np.array(["clean_s0", "clean_s1"]),
        voxel_r_perp_um=np.array([10.0, 50, 120, 200, 300, 400]),
        err_angle=err,
        err_x=err / 2,
        err_y=err / 2,
        err_z=err / 2,
        maha2=np.full(shape, 3.0, dtype=np.float32),
        fail_pass1=fail,
        stop_reason=np.zeros((V, R, D, Va, Mo), dtype=np.int8),
    )
    s = ps.summarize(raw)
    assert s["failures"]["clean"]["5"]["n_fail"] == 2
    assert s["failures"]["clean"]["5"]["few_present"] == 2
    assert s["failures"]["clean"]["0.5"]["n_fail"] == 0
    m = s["metrics"]["all"]["clean"]["clean_s0"]["5"][0]
    assert m["n_ok"] == V * D - 2 and m["frac_lt_0p1"] == 1.0 and m["frac_improved"] == 1.0
    assert abs(m["mean_maha2"] - 3.0) < 1e-6
    assert sum(s["tercile_n_voxels"]) == V
    assert "clean-trained" not in ps.format_summary(raw, s)  # runs on a tiny synthetic raw
