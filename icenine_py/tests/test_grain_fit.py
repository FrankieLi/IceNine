"""Full-grain fitting pieces (icenine/grain_fit.py).

Pure pieces (rotations, grouping, sub-triangles, Sinkhorn, fields) on synthetic inputs; the spot
model against ``SmoothCost.predict`` and the coverage points against the stock-rendered images of
the ThreeVoxels toy problem (real two-detector images)."""

import sys
from pathlib import Path
from typing import Any, List

import numpy as np
import pytest
import torch
from scipy.spatial.transform import Rotation

from icenine import grain_fit as GF

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))

DT = torch.float64


# ------------------------------------------------------------------------------ rotations
def test_rotvec_roundtrip_and_small_angle() -> None:
    rng = np.random.default_rng(0)
    w = rng.standard_normal((20, 3)) * 0.3
    w[0] = 0.0
    w[1] = 1e-9
    R = GF.rotvec_to_matrix(torch.as_tensor(w)).numpy()
    ref = Rotation.from_rotvec(w).as_matrix()
    assert np.abs(R - ref).max() < 1e-12
    assert np.abs(GF.matrix_to_rotvec(R) - w).max() < 1e-9
    # gradient finite at zero
    z = torch.zeros(3, dtype=DT, requires_grad=True)
    GF.rotvec_to_matrix(z).sum().backward()
    assert torch.isfinite(z.grad).all()


def _cubic_sym() -> np.ndarray:
    from icenine.ray_indexing import proper_rotations
    from icenine.symmetry import create_cubic_symmetry

    return proper_rotations(create_cubic_symmetry(3.6))


def test_align_and_mean_orientation_under_symmetry() -> None:
    sym = _cubic_sym()
    R0 = Rotation.random(random_state=3).as_matrix()
    rng = np.random.default_rng(1)
    small = Rotation.from_rotvec(rng.standard_normal((10, 3)) * 1e-3).as_matrix()
    R = small @ R0
    # replace by symmetry equivalents
    R = np.einsum("nij,njk->nik", R, sym[rng.integers(0, len(sym), 10)])
    A = GF.align_to(R, R0, sym)
    assert np.abs(A - small @ R0).max() < 1e-12
    M = GF.mean_orientation(R, sym, R0)
    from icenine.ray_indexing import misorientation_deg

    assert float(misorientation_deg(M[None], R0[None], sym)[0]) < 0.1


# ------------------------------------------------------------------------------ grouping
def _grid_neighbours(nx: int, ny: int) -> List[List[int]]:
    nb: List[List[int]] = []
    for j in range(ny):
        for i in range(nx):
            out = []
            for di, dj in ((1, 0), (-1, 0), (0, 1), (0, -1)):
                a, b = i + di, j + dj
                if 0 <= a < nx and 0 <= b < ny:
                    out.append(b * nx + a)
            nb.append(out)
    return nb


def test_group_grains_splits_by_misorientation_and_orders_by_size() -> None:
    sym = _cubic_sym()
    nx, ny = 6, 4
    nb = _grid_neighbours(nx, ny)
    Ra = Rotation.random(random_state=5).as_matrix()
    Rb = Rotation.random(random_state=6).as_matrix()
    R = np.empty((nx * ny, 3, 3))
    for j in range(ny):
        for i in range(nx):
            v = j * nx + i
            # a gradient of 0.5 deg per voxel inside grain A (still below the 2 deg threshold
            # between neighbours), grain B on the last two columns
            if i < 4:
                R[v] = Rotation.from_rotvec([0, 0, np.radians(0.5 * i)]).as_matrix() @ Ra
            else:
                R[v] = Rb
    lab = GF.group_grains(R, nb, sym, 2.0)
    assert set(lab[[j * nx + i for j in range(ny) for i in range(4)]]) == {0}
    assert set(lab[[j * nx + i for j in range(ny) for i in (4, 5)]]) == {1}
    members = [v for v in range(nx * ny) if lab[v] == 1]
    assert GF.shell_of(members, nb) == sorted(j * nx + 3 for j in range(ny))
    assert sorted(GF.boundary_members(members, nb)) == sorted(j * nx + 4 for j in range(ny))


# ------------------------------------------------------------------------------ sub-triangles
@pytest.mark.parametrize("n", [1, 2, 3, 5])
def test_subtriangle_centroids_cover_the_triangle(n: int) -> None:
    b = GF.subtriangle_bary(n)
    assert len(b) == n * n
    assert (b > 0).all() and (b.sum(1) < 1).all()
    # the mean of equal-area sub-triangle centroids is the triangle's centroid
    assert np.allclose(b.mean(0), [1 / 3, 1 / 3])
    col = torch.tensor([[0.0, 4.0, 0.0]], dtype=DT)
    row = torch.tensor([[0.0, 0.0, 2.0]], dtype=DT)
    pc, pr, m = GF.coverage_points(col, row, torch.as_tensor(b))
    assert abs(float(m.sum()) - 4.0) < 1e-12  # area of the triangle


# ------------------------------------------------------------------------------ Sinkhorn
def _cloud(n: int, seed: int) -> torch.Tensor:
    g = torch.Generator().manual_seed(seed)
    return torch.rand((1, n, 3), generator=g, dtype=DT) * torch.tensor([10.0, 10.0, 0.0], dtype=DT)


def test_sinkhorn_divergence_zero_at_identity_and_shift_squared() -> None:
    x = _cloud(60, 0)
    x[..., 2] = 5.5  # inside frame 5
    a = torch.ones((1, 60), dtype=DT)
    b = a / a.sum()
    d0 = GF.sinkhorn_divergence(x, a, x.clone(), b, blur_px=0.5, diameter_px=20.0)
    assert abs(float(d0)) < 1e-6
    y = x.clone()
    y[..., 0] += 2.0
    d = GF.sinkhorn_divergence(x, a, y, b, blur_px=0.2, diameter_px=20.0)
    assert abs(float(d) - 4.0) < 0.3  # W2^2 of a translation is the shift squared
    # frame hinge: an fpos inside the observed frame costs nothing, one frame off costs lam^2/4
    z = x.clone()
    z[..., 2] = 5.1
    assert abs(float(GF.sinkhorn_divergence(z, a, x, b, blur_px=0.5, diameter_px=20.0))) < 1e-6


def test_sinkhorn_gradient_points_toward_the_observed_spot() -> None:
    x = _cloud(40, 1)
    y = x.clone()
    y[..., 0] += 1.5
    s = torch.zeros(1, dtype=DT, requires_grad=True)
    a = torch.ones((1, 40), dtype=DT)
    b = a / a.sum()
    xs = x + torch.stack([s, torch.zeros_like(s), torch.zeros_like(s)], -1)
    for reach in (None, 3.0):
        s.grad = None
        d = GF.sinkhorn_divergence(xs, a, y, b, blur_px=0.5, diameter_px=20.0, reach_px=reach)
        d.sum().backward(retain_graph=True)
        assert float(s.grad) < 0  # moving +x lowers the distance


def test_sinkhorn_empty_prediction_gets_the_penalty() -> None:
    y = _cloud(10, 2)
    b = torch.full((1, 10), 0.1, dtype=DT)
    a = torch.zeros((1, 5), dtype=DT)
    x = torch.zeros((1, 5, 3), dtype=DT)
    d = GF.sinkhorn_divergence(x, a, y, b, blur_px=1.0, reach_px=3.0)
    assert float(d) == 9.0


def test_ls_objective_is_scale_free() -> None:
    P = torch.rand((2, 5, 5), dtype=DT)
    assert abs(float(GF.ls_objective([P], [3.0 * P])[0])) < 1e-12
    Q = torch.zeros_like(P)
    Q[0, 0, 0] = 1.0
    P2 = torch.zeros_like(P)
    P2[1, 1, 1] = 1.0
    assert abs(float(GF.ls_objective([P2], [Q])[0]) - 1.0) < 1e-12


# ------------------------------------------------------------------------------ fields
@pytest.mark.parametrize("order", [1, 2])
def test_field_fit_recovers_a_polynomial_field(order: int) -> None:
    sym = _cubic_sym()
    rng = np.random.default_rng(2)
    xy = rng.uniform(-1, 1, (40, 2))
    R_ref = Rotation.random(random_state=7).as_matrix()
    F = GF.FieldA.build(R_ref, xy, order=order)
    theta = rng.standard_normal(F.n_params) * np.radians(0.3)
    R = F.orientations(torch.as_tensor(theta)).numpy()
    # symmetry-equivalent inputs
    R = np.einsum("nij,njk->nik", R, sym[rng.integers(0, len(sym), len(R))])
    th = F.fit_theta(R, sym)
    # exp of a sum is not the sum of exps, but the rotation vectors here are below 1 deg, so
    # the fit of log(R R_ref^T) is linear in theta to second order
    assert np.abs(th - theta).max() < np.radians(0.01)


# ------------------------------------------------------------------------------ toy problem
def _toy() -> Any:
    from test_findoptimal_refactor import _build

    return _build()


def _sc(rec: Any) -> Any:
    from icenine.smooth_cost import SmoothCost

    s = rec.setup
    return SmoothCost(
        s.simulator,
        s.detector_list,
        s.range_map,
        s.exp_data,
        s.sample,
        s.structure_list,
        eta_limit=s.exp_setup.get_eta_limit(),
    )


def test_spot_model_equals_smooth_cost_predict() -> None:
    from icenine.reconstructor import _get_voxel_vertices

    rec, voxel, R_true = _toy()
    sc = _sc(rec)
    phase = int(voxel.phase)
    model = GF.SpotModel(sc, phase=phase)
    mic = rec.setup.sample.get_mic()
    Rs = np.stack([np.asarray(v.orientation, np.float64) for v in mic.voxels])
    Vs = np.stack([_get_voxel_vertices(v).numpy().astype(np.float64) for v in mic.voxels])
    rows = model.predict(torch.as_tensor(Rs), torch.as_tensor(Vs))
    for i in range(len(Rs)):
        ref = sc.predict(Rs[i : i + 1], torch.as_tensor(Vs[i]), phase)
        sel = rows.vox == i
        assert int(sel.sum()) == len(ref.refl)
        assert torch.equal(rows.refl[sel], ref.refl) and torch.equal(rows.branch[sel], ref.branch)
        assert torch.allclose(rows.fpos[sel], ref.fpos, atol=1e-12)
        assert torch.equal(rows.frame[sel], ref.frame)
        # SmoothCost samples the 3 vertices + the centroid
        assert torch.allclose(rows.col[sel], ref.col[..., :3], atol=1e-9)
        assert torch.allclose(rows.row[sel], ref.row[..., :3], atol=1e-9)


def test_coverage_points_land_on_stock_rendered_pixels() -> None:
    """At the true orientations, nearly all coverage mass falls on lit pixels of the stock-
    rendered ThreeVoxels images, in the predicted frame (pixel = floor of the coordinate)."""
    from icenine.reconstructor import _get_voxel_vertices

    rec, voxel, _R = _toy()
    sc = _sc(rec)
    model = GF.SpotModel(sc, phase=int(voxel.phase))
    mic = rec.setup.sample.get_mic()
    Rs = np.stack([np.asarray(v.orientation, np.float64) for v in mic.voxels])
    Vs = np.stack([_get_voxel_vertices(v).numpy().astype(np.float64) for v in mic.voxels])
    rows = model.predict(torch.as_tensor(Rs), torch.as_tensor(Vs))
    bary = torch.as_tensor(GF.subtriangle_bary(3))
    exp = rec.setup.exp_data
    on, tot = 0.0, 0.0
    for det in range(model.n_det):
        pc, pr, m, ok = GF._rows_points(rows, bary, det)
        for i in torch.nonzero(ok)[:, 0].tolist():
            f = int(rows.frame[i])
            img = exp.get_image(f, det)
            c = torch.floor(pc[i]).long()
            r = torch.floor(pr[i]).long()
            inb = (c >= 0) & (c < 2048) & (r >= 0) & (r < 2048)
            if not bool(inb.any()):
                continue
            lit = (img.get_pixel(c[inb], r[inb]) > 0).numpy()  # j = column, k = row
            mm = m[i].numpy()[inb.numpy()]
            tot += float(mm.sum())
            on += float((mm * lit).sum())
    assert tot > 0
    assert on / tot > 0.95


def test_touching_components_drop_far_spots() -> None:
    obs_grid = np.zeros((1, 20, 20), bool)
    obs_grid[0, 2:4, 2:4] = True  # near the prediction
    obs_grid[0, 15:17, 15:17] = True  # far away: another grain
    f, r, c = np.nonzero(obs_grid)
    obs = np.stack([c + 0.5, r + 0.5, f + 0.5], 1)
    pc = torch.tensor([[2.5, 3.1]], dtype=DT)
    pr = torch.tensor([[2.5, 2.9]], dtype=DT)
    keep = GF._touching_components(obs, (0, 1, 0, 20, 0, 20), pc, pr, torch.tensor([0]), 2)
    assert keep.sum() == 4 and (obs[keep][:, 0] < 5).all()
