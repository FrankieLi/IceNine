"""
Windowed local-refinement data generation for the toy orientation NN.

Real detector frames are 2048x2048 x (n_omega) per voxel -- far too large to
render per training sample. This module exploits the fact that local
refinement starts from a known nominal orientation: the fixed set of
observable peaks can be predicted once (`define_roi_set`), and only small
windows around those peaks' nominal pixel locations need to be rendered per
perturbed sample (`render_local_windows`), by reusing the same ray-tracing
primitives as `ForwardSimulation._simulate_peaks` (`diffraction_core.
get_reflected_ray_dir`/`build_reflected_ray`/`get_illuminated_pixel`) instead
of `Simulation.project_voxel_multi_detector`'s full-detector rasterization.

See icenine_py/MIGRATION_HISTORY.md for the approved plan this implements.
"""

import math
from dataclasses import dataclass
from typing import List, Optional, Tuple

import torch
from torch.utils.data import Dataset

from .constants import KEV_OVER_HBAR_C_IN_ANG
from .detector import Detector
from .diffraction_core import (
    build_reflected_ray,
    get_illuminated_pixel,
    get_reflected_ray_dir,
    get_scattering_omegas_torch,
)
from .image_data import ImageData
from .peak_filters import XDMEtaAcceptFn
from .sample import Sample
from .simulation import Simulation


@dataclass
class ROIPeak:
    """One fixed (reflection, omega-branch, detector) peak identity in an ROI set.

    Defined once at the nominal orientation; `render_local_windows` re-resolves
    the same reflection+branch under a perturbed orientation rather than
    running a fresh peak search, so peak identity stays constant across a
    dataset.
    """

    reflection_index: int
    omega_branch: int  # 1 or 2 -> get_scattering_omegas_torch's omega1/omega2
    detector_index: int
    g_hkl: torch.Tensor  # (3,) crystal-frame reciprocal vector, fixed
    form_intensity: float
    sin_2theta: float
    nominal_omega: float  # radians
    nominal_row: float  # window is centered here, fixed for the whole ROI set
    nominal_col: float


def _restore_and_rotate(sample: Sample, base_rotation: torch.Tensor, omega: float) -> None:
    """Reset sample's global rotation to base_rotation, then apply Rz(omega)."""
    sample.sample_to_lab_matrix[:3, :3] = base_rotation.clone()
    sample.rotate_z(omega)


def spot_overlaps_grid(pixels: List[Tuple[float, float]], num_cols: int, num_rows: int) -> bool:
    """True if a spot with the given (col, row) vertices overlaps the pixel grid.

    Uses the rasteriser's own vertex truncation (negative -> -1, else int()), so a
    spot is recorded exactly when the bounding box of its truncated vertices meets
    [0, num_cols-1] x [0, num_rows-1].
    """
    cols = [-1 if c < 0 else int(c) for c, _ in pixels]
    rows = [-1 if r < 0 else int(r) for _, r in pixels]
    return (
        max(cols) >= 0
        and min(cols) <= num_cols - 1
        and max(rows) >= 0
        and min(rows) <= num_rows - 1
    )


def _project_peak_on_detector(
    simulator: Simulation,
    sample: Sample,
    detector: Detector,
    voxel_vertices: torch.Tensor,
    scattering_dir: torch.Tensor,
    peak_filter: XDMEtaAcceptFn,
    require_grid: bool = True,
) -> Optional[Tuple[float, float, float, List[Tuple[float, float]]]]:
    """Project a voxel's 3 vertices onto one detector at the sample's current
    rotation state.

    Mirrors Simulation.project_voxel_multi_detector's inner loop, but returns
    raw (unrasterized) pixel coordinates instead of rasterizing into a
    full-size image, so callers can rasterize into an offset local window.

    Returns (row0, col0, intensity, [(col,row) x3]) -- row0/col0 are the
    centroid of the 3 projected vertices -- or None if the peak is filtered
    out, any vertex misses the detector plane, or (if require_grid) the spot does
    not overlap the detector's pixel grid (the rasteriser would clip it away, so
    it is never recorded).
    """
    reflected_dir = get_reflected_ray_dir(sample, scattering_dir, simulator.beam_direction)
    reflected_dir_n = reflected_dir / torch.norm(reflected_dir)
    accept, intensity = peak_filter(reflected_dir_n)
    if not accept:
        return None

    pixels = []
    for vi in range(3):
        ray = build_reflected_ray(sample, voxel_vertices[vi], reflected_dir)
        hit, col, row = get_illuminated_pixel(detector, ray)
        if not hit.item():
            return None
        pixels.append((float(col.item()), float(row.item())))

    if require_grid and not spot_overlaps_grid(pixels, detector.num_cols, detector.num_rows):
        return None

    row0 = sum(p[1] for p in pixels) / 3.0
    col0 = sum(p[0] for p in pixels) / 3.0
    return row0, col0, intensity, pixels


def _project_all_detectors(
    simulator: Simulation,
    sample: Sample,
    detector_list: List[Detector],
    voxel_vertices: torch.Tensor,
    scattering_dir: torch.Tensor,
    peak_filter: XDMEtaAcceptFn,
) -> Optional[List[Tuple[float, float, float, List[Tuple[float, float]]]]]:
    """Project a spot onto every detector, following ForwardSimulation._simulate_peaks:
    if the peak is filtered out or any vertex misses any detector *plane*, the
    peak is dropped on all detectors (returns None). Otherwise returns one
    (row0, col0, intensity, pixels) per detector, without the pixel-grid check.
    """
    results = []
    for detector in detector_list:
        r = _project_peak_on_detector(
            simulator,
            sample,
            detector,
            voxel_vertices,
            scattering_dir,
            peak_filter,
            require_grid=False,
        )
        if r is None:
            return None
        results.append(r)
    return results


def define_roi_set(
    nominal_orientation: torch.Tensor,
    voxel_vertices: torch.Tensor,
    sample: Sample,
    detector_list: List[Detector],
    range_map,
    exp_setup,
    structure_list,
    simulator: Simulation,
    phase_index: int = 0,
    detectors: str = "first",
) -> List[ROIPeak]:
    """Ray-trace the voxel's observable peaks at its nominal orientation once,
    fixing the peak identity (reflection, branch, detector) and window center
    (nominal_row/col) used for the whole dataset.

    detectors:
      "all"   -- the simulator's semantics: a peak is dropped everywhere if any spot
                 vertex misses any detector plane; otherwise it gets one entry per
                 detector whose pixel grid its spot overlaps (Stage 1 onwards).
      "first" -- one entry, on the first detector whose grid the spot overlaps
                 (Stage 0 behaviour, kept for reproducibility).
    """
    if detectors not in ("all", "first"):
        raise ValueError(f"detectors must be 'all' or 'first', got {detectors!r}")
    structure = structure_list[phase_index]
    reflections = structure.get_reflection_vectors()
    if not reflections:
        return []

    g_hkl_batch = torch.stack([torch.from_numpy(r.q_vec).float() for r in reflections])
    g_mag_batch = torch.tensor([r.q_mag for r in reflections], dtype=torch.float32)

    wavenumber = KEV_OVER_HBAR_C_IN_ANG * exp_setup.beam_energy
    sin_theta = (g_mag_batch / (2.0 * wavenumber)).clamp(-1.0, 1.0)
    sin_2theta = torch.sin(2.0 * torch.asin(sin_theta))

    g_lab_batch = (nominal_orientation @ g_hkl_batch.T).T
    g_mag_lab = torch.norm(g_lab_batch, dim=1)
    bragg = get_scattering_omegas_torch(
        g_lab_batch, g_mag_lab, exp_setup.beam_energy, exp_setup.get_beam_deflection_chi_laue()
    )

    eta_limit = exp_setup.get_eta_limit()
    base_rotation = sample.sample_to_lab_matrix[:3, :3].clone()

    roi_list: List[ROIPeak] = []
    try:
        for i in range(len(reflections)):
            if not bragg.observable[i]:
                continue
            scattering_dir = g_lab_batch[i] / g_mag_lab[i]
            for branch, omega_t in ((1, bragg.omega1[i]), (2, bragg.omega2[i])):
                omega = float(omega_t.item())
                if range_map.angle_to_wedge_index(omega) is None:
                    continue

                _restore_and_rotate(sample, base_rotation, omega)
                peak_filter = XDMEtaAcceptFn(
                    0.0, eta_limit, float(reflections[i].intensity), float(sin_2theta[i].item())
                )

                def entry(det_idx: int, row0: float, col0: float) -> ROIPeak:
                    return ROIPeak(
                        reflection_index=i,
                        omega_branch=branch,
                        detector_index=det_idx,
                        g_hkl=g_hkl_batch[i].clone(),
                        form_intensity=float(reflections[i].intensity),
                        sin_2theta=float(sin_2theta[i].item()),
                        nominal_omega=omega,
                        nominal_row=row0,
                        nominal_col=col0,
                    )

                if detectors == "all":
                    results = _project_all_detectors(
                        simulator,
                        sample,
                        detector_list,
                        voxel_vertices,
                        scattering_dir,
                        peak_filter,
                    )
                    if results is None:
                        continue
                    for det_idx, (row0, col0, _intensity, pixels) in enumerate(results):
                        det = detector_list[det_idx]
                        if spot_overlaps_grid(pixels, det.num_cols, det.num_rows):
                            roi_list.append(entry(det_idx, row0, col0))
                else:
                    for det_idx, detector in enumerate(detector_list):
                        result = _project_peak_on_detector(
                            simulator, sample, detector, voxel_vertices, scattering_dir, peak_filter
                        )
                        if result is not None:
                            roi_list.append(entry(det_idx, result[0], result[1]))
                            break
    finally:
        _restore_and_rotate(sample, base_rotation, 0.0)

    return roi_list


def render_local_windows(
    perturbed_orientation: torch.Tensor,
    roi_list: List[ROIPeak],
    voxel_vertices: torch.Tensor,
    sample: Sample,
    detector_list: List[Detector],
    range_map,
    exp_setup,
    simulator: Simulation,
    window_size: int = 32,
) -> Tuple[torch.Tensor, torch.Tensor]:
    """Render each ROI peak's fixed-size window under a perturbed orientation.

    Each window is centered on the peak's *nominal* pixel location (fixed at
    ROI-definition time), so the peak's actual position drifts within the
    window as the orientation perturbs -- this drift is exactly the signal
    the network learns to read. If a peak stops being observable (Bragg
    condition fails, falls outside the exposed omega range, a vertex misses any
    detector plane, or its spot leaves its detector's pixel grid) under this
    perturbation, its window is zero-filled and flagged in the returned mask.

    Returns:
        windows: (n_peaks, window_size, window_size) float32
        missing: (n_peaks,) bool -- True where the peak dropped out
    """
    n = len(roi_list)
    windows = torch.zeros(n, window_size, window_size, dtype=torch.float32)
    missing = torch.zeros(n, dtype=torch.bool)
    if n == 0:
        return windows, missing

    eta_limit = exp_setup.get_eta_limit()
    base_rotation = sample.sample_to_lab_matrix[:3, :3].clone()
    half = window_size / 2.0

    try:
        for i, roi in enumerate(roi_list):
            g_lab = perturbed_orientation @ roi.g_hkl
            g_mag = torch.norm(g_lab)
            bragg = get_scattering_omegas_torch(
                g_lab.unsqueeze(0),
                g_mag.unsqueeze(0),
                exp_setup.beam_energy,
                exp_setup.get_beam_deflection_chi_laue(),
            )
            if not bragg.observable[0]:
                missing[i] = True
                continue

            omega_t = bragg.omega1[0] if roi.omega_branch == 1 else bragg.omega2[0]
            omega = float(omega_t.item())
            if range_map.angle_to_wedge_index(omega) is None:
                missing[i] = True
                continue

            _restore_and_rotate(sample, base_rotation, omega)
            scattering_dir = g_lab / g_mag
            peak_filter = XDMEtaAcceptFn(0.0, eta_limit, roi.form_intensity, roi.sin_2theta)
            # Simulator semantics: a vertex missing any detector plane drops the peak
            # everywhere; then the spot must overlap its own detector's pixel grid.
            results = _project_all_detectors(
                simulator, sample, detector_list, voxel_vertices, scattering_dir, peak_filter
            )
            detector = detector_list[roi.detector_index]
            result = None if results is None else results[roi.detector_index]
            if result is None or not spot_overlaps_grid(
                result[3], detector.num_cols, detector.num_rows
            ):
                missing[i] = True
                continue

            _row0, _col0, intensity, pixels = result
            origin_col = roi.nominal_col - half
            origin_row = roi.nominal_row - half
            v0 = torch.tensor([pixels[0][0] - origin_col, pixels[0][1] - origin_row])
            v1 = torch.tensor([pixels[1][0] - origin_col, pixels[1][1] - origin_row])
            v2 = torch.tensor([pixels[2][0] - origin_col, pixels[2][1] - origin_row])

            canvas = ImageData(window_size, window_size)
            canvas.add_triangle_scanline(v0, v1, v2, intensity)
            windows[i] = torch.from_numpy(canvas.to_numpy())
    finally:
        _restore_and_rotate(sample, base_rotation, 0.0)

    return windows, missing


def sample_local_perturbations(
    nominal_orientation,
    n_samples: int,
    max_angle_deg: float,
    rng,
) -> Tuple[List, List]:
    """Random small-angle perturbations composed onto the nominal orientation.

    Uses the same near-identity quaternion parameterization as MCOptimizer's
    perturbation step (orientation_search.py): (x, y, z) sampled uniformly in
    a box, mapped to a near-identity quaternion via QuaternionGrid's
    barycentric parameterization, composed onto the nominal quaternion.

    Args:
        nominal_orientation: (3,3) rotation matrix (np.ndarray or torch.Tensor)
        n_samples: number of perturbations to draw
        max_angle_deg: approximate bound on the perturbation angle
        rng: numpy.random.Generator

    Returns:
        (matrices, quaternions): matrices is a list of (3,3) np.ndarray rotation
        matrices (physics-layer representation); quaternions is a list of
        (4,) np.ndarray [w,x,y,z] unit quaternions (network target representation).
    """
    import numpy as np

    from .sampling import QuaternionGrid, _quat_multiply, matrix_to_quaternion, quaternion_to_matrix

    grid_gen = QuaternionGrid()
    nominal_np = (
        nominal_orientation.numpy()
        if isinstance(nominal_orientation, torch.Tensor)
        else nominal_orientation
    )
    nominal_q = matrix_to_quaternion(nominal_np)

    max_angle_rad = math.radians(max_angle_deg)
    radius = math.tan(max_angle_rad) / math.sqrt(12.0)

    matrices = []
    quats = []
    for _ in range(n_samples):
        x = rng.uniform(-radius, radius)
        y = rng.uniform(-radius, radius)
        z = rng.uniform(-radius, radius)
        delta_q = grid_gen.get_near_identity_point(x, y, z)
        trial_q = _quat_multiply(delta_q, nominal_q)
        matrices.append(quaternion_to_matrix(trial_q))
        quats.append(trial_q)

    return matrices, quats


class OrientationDataset(Dataset):
    """Cached (windows, offset) pairs written by scripts/generate_toy_orientation_dataset.py.

    The file holds `windows` as a uint8 tensor (N, n_peaks, window, window) of
    thresholded (lit / not lit) pixels and `offsets_deg` (N, 3), the rotation-vector
    offsets from the nominal orientation in degrees. Windows are returned as float32.
    """

    def __init__(self, cache_path: str):
        data = torch.load(cache_path)
        self.windows: torch.Tensor = data["windows"]
        self.offsets_deg: torch.Tensor = data["offsets_deg"].float()
        self.meta = {k: v for k, v in data.items() if k not in ("windows", "offsets_deg")}

    def __len__(self) -> int:
        return self.windows.shape[0]

    def __getitem__(self, idx: int) -> Tuple[torch.Tensor, torch.Tensor]:
        return self.windows[idx].float(), self.offsets_deg[idx]


def quaternion_regression_loss(q_pred: torch.Tensor, q_true: torch.Tensor) -> torch.Tensor:
    """Training loss: 1 - |dot(q_pred, q_true)|, batched over (B, 4).

    Smooth near the optimum (dot -> 1), unlike a direct port of the
    arccos-based misorientation-degree metric, which has an unbounded
    gradient there.
    """
    dot = (q_pred * q_true).sum(dim=-1)
    return (1.0 - dot.abs()).mean()


def quat_misorientation_deg_batch(q1: torch.Tensor, q2: torch.Tensor) -> torch.Tensor:
    """Batched port of orientation_search._quat_misorientation_deg (2*arccos(|q1.q2|),
    in degrees). Eval-only: not used inside the training loss.
    """
    dot = (q1 * q2).sum(dim=-1).abs().clamp(max=1.0)
    return torch.rad2deg(2.0 * torch.acos(dot))


def cholesky_from_raw(raw: torch.Tensor, min_diag: float = 1e-4) -> torch.Tensor:
    """Lower-triangular Cholesky factors from 6 unconstrained numbers per sample.

    raw[..., 0:3] -> diagonal (softplus + min_diag), raw[..., 3:6] -> the strictly
    lower entries (1,0), (2,0), (2,1). Returns (..., 3, 3).
    """
    diag = torch.nn.functional.softplus(raw[..., 0:3]) + min_diag
    L10, L20, L21 = raw[..., 3], raw[..., 4], raw[..., 5]
    zero = torch.zeros_like(L10)
    lower = torch.stack([zero, zero, zero, L10, zero, zero, L20, L21, zero], dim=-1).reshape(
        *raw.shape[:-1], 3, 3
    )
    return torch.diag_embed(diag) + lower


def gaussian_nll_loss(
    mean: torch.Tensor, chol: torch.Tensor, target: torch.Tensor, beta: float = 0.0
) -> torch.Tensor:
    """Multivariate Gaussian negative log-likelihood with covariance L L^T.

    nll = 0.5 |L^-1 (target - mean)|^2 + sum log diag(L)   (constants dropped).
    beta > 0 gives the beta-NLL reweighting of Seitzer et al. (2022): each sample
    is weighted by stopgrad(det(cov)^(beta/3)).
    """
    r = (target - mean).unsqueeze(-1)
    z = torch.linalg.solve_triangular(chol, r, upper=False).squeeze(-1)
    log_diag = torch.log(torch.diagonal(chol, dim1=-2, dim2=-1))
    nll = 0.5 * (z**2).sum(-1) + log_diag.sum(-1)
    if beta > 0.0:
        weight = torch.exp(2.0 * log_diag.sum(-1) * beta / 3.0).detach()
        nll = nll * weight
    return nll.mean()


def mse_deg_loss(mean: torch.Tensor, target: torch.Tensor, scale_deg: float = 0.1) -> torch.Tensor:
    """0.5 * sum over axes of the squared error, in units of scale_deg**2, batch mean.

    Every axis has the same weight, so the stage axis is not down-weighted when the
    covariance says it is uncertain.
    """
    return (0.5 * ((target - mean) ** 2).sum(-1) / scale_deg**2).mean()


def decoupled_nll_loss(
    mean: torch.Tensor, chol: torch.Tensor, target: torch.Tensor, scale_deg: float = 0.1
) -> torch.Tensor:
    """Decoupled mean / covariance loss.

    MSE on the mean (equal weight per axis, see mse_deg_loss) plus the Gaussian NLL of
    the covariance evaluated at stopgrad(mean). The covariance term sends no gradient to
    the mean, so the mean's gradient is the plain MSE gradient regardless of sigma (the
    Seitzer et al. 2022 fix for the 1/sigma^2 scaling of the NLL mean gradient), while
    the covariance is still fitted to the residuals the mean actually makes.
    """
    return mse_deg_loss(mean, target, scale_deg) + gaussian_nll_loss(mean.detach(), chol, target)
