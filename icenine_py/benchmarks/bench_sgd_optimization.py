#!/usr/bin/env python3
"""
Riemannian SGD on SO(3) benchmark.

Compares five SGD-family orientation optimizers against the Riemannian Adam
baseline, all minimizing the same differentiable cost function:

  riemannian_sgd_plain      — pure steepest descent on SO(3): R ← R·exp(−lr·Ω)
                              No momentum, no adaptive scaling.
  riemannian_sgd_momentum   — heavy-ball momentum β=0.9 in so(3) coordinates.
                              m ← β·m + Ω;  R ← R·exp(−lr·m)
  riemannian_sgd_nesterov   — Nesterov look-ahead variant of momentum.
                              geoopt nesterov=True.
  riemannian_sgd_cosine     — SGD + momentum with cosine-annealing LR schedule.
                              lr_max=0.05 → lr_min=0.001 over n_steps.
  riemannian_sgld           — Stochastic Gradient Langevin Dynamics on SO(3).
                              Adds isotropic Brownian motion on T_R SO(3) at
                              each step, annealing temperature to zero.
                              Principled approach to escape flat cost landscape.
  riemannian_adam_manual    — Riemannian Adam (manual, from bench_riemannian_optimization.py)
                              Serves as Adam comparison baseline.

All six optimizers use the same voxels, perturbation axes (SEED=42), and sweep
parameters as bench_riemannian_optimization.py for direct comparison.

Why SGD variants can differ from Adam on SO(3):
  - Adam suppresses gradients via √m̂₂ denominator; in flat regions where
    ‖G‖≈0 early on, this can reduce the effective learning rate to lr/ε ≈ ∞
    (first step) or near-zero (when m2 accumulates small values). SGD has
    constant effective lr, which can be more aggressive in flat regions.
  - SGLD adds Langevin noise: noise_scale = √(2·lr·T), T annealing linearly.
    This allows probabilistic escape from flat/spurious basins — the same
    physics that the original adaptive MC search uses, but guided by the
    gradient.

Usage:
  cd icenine_py
  uv sync --extra riemannian
  uv run python benchmarks/bench_sgd_optimization.py --smoke-test --example threevoxels
  uv run python benchmarks/bench_sgd_optimization.py --example threevoxels
  uv run python benchmarks/bench_sgd_optimization.py --example manygrains

Outputs (icenine_py/benchmarks/):
  sgd_opt_threevoxels.csv / sgd_opt_manygrains.csv
  sgd_opt_convergence_{example}.png  — quality vs step, lines per optimizer
  sgd_opt_misori_{example}.png       — misorientation vs step
  sgd_opt_summary_{example}.png      — bar chart: final misori by optimizer × pert
"""

import argparse
import csv
import math
import os
import sys
import time
from pathlib import Path
from typing import Dict, List, Optional

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import torch
import torch.nn.functional as F  # noqa: F401 — imported for completeness

project_root = Path(__file__).parent.parent.parent
sys.path.insert(0, str(project_root / "icenine_py"))

benchmark_dir = Path(__file__).parent

# ---------------------------------------------------------------------------
# Sweep parameters (match bench_riemannian_optimization.py)
# ---------------------------------------------------------------------------

SCALES = [1, 2]
OMEGA_WINDOWS = [0, 1, 2]
PERTURBATIONS_DEG = [1.0, 2.0, 5.0]
N_STEPS = 100
LR = 0.01  # base lr for plain/momentum/sgld
LR_MAX = 0.05  # cosine schedule peak
LR_MIN = 0.001  # cosine schedule trough
MOMENTUM = 0.9  # for momentum / nesterov / cosine variants
SGLD_T_INIT = 0.01  # initial Langevin temperature (tunable via --sgld-temp)
SEED = 42
N_VOXELS_MANY = 20
QUALITY_THRESHOLD = 0.1
MAX_SCAN = 500

# Adam hyperparameters (for the baseline copy)
BETA1 = 0.9
BETA2 = 0.999
ADAM_EPS = 1e-8

OPTIMIZER_COLORS = {
    "riemannian_sgd_plain": "C0",
    "riemannian_sgd_momentum": "C1",
    "riemannian_sgd_nesterov": "C2",
    "riemannian_sgd_cosine": "C4",
    "riemannian_sgld": "C5",
    "riemannian_adam_manual": "C3",
}
OPTIMIZER_LINESTYLES = {
    "riemannian_sgd_plain": "-",
    "riemannian_sgd_momentum": "--",
    "riemannian_sgd_nesterov": "-.",
    "riemannian_sgd_cosine": ":",
    "riemannian_sgld": (0, (3, 1, 1, 1)),
    "riemannian_adam_manual": (0, (5, 2)),
}


# ---------------------------------------------------------------------------
# SO(3) helpers (identical to bench_riemannian_optimization.py)
# ---------------------------------------------------------------------------


def skew(theta: torch.Tensor) -> torch.Tensor:
    """3-vector → 3×3 skew-symmetric matrix (Lie algebra so(3))."""
    assert theta.shape == (3,), f"Expected shape (3,), got {theta.shape}"
    z = torch.zeros(1, dtype=theta.dtype, device=theta.device)
    row0 = torch.stack([z.squeeze(), -theta[2], theta[1]])
    row1 = torch.stack([theta[2], z.squeeze(), -theta[0]])
    row2 = torch.stack([-theta[1], theta[0], z.squeeze()])
    return torch.stack([row0, row1, row2])


def so3_log(R: np.ndarray) -> torch.Tensor:
    """Rotation matrix → axis-angle 3-vector (for initialization only)."""
    R_t = torch.from_numpy(R).float()
    trace = R_t.trace()
    cos_angle = ((trace - 1.0) / 2.0).clamp(-1.0, 1.0)
    angle = torch.acos(cos_angle)
    if angle.abs() < 1e-7:
        return torch.zeros(3)
    W = (R_t - R_t.T) / (2.0 * torch.sin(angle))
    axis = torch.stack([W[2, 1], W[0, 2], W[1, 0]])
    return angle * axis


def misorientation_deg(R1: np.ndarray, R2: np.ndarray) -> float:
    """Geodesic distance on SO(3) in degrees."""
    M = R1.T @ R2
    trace = np.trace(M)
    cos_angle = np.clip((trace - 1.0) / 2.0, -1.0, 1.0)
    return float(np.degrees(np.arccos(cos_angle)))


def rodrigues_np(axis: np.ndarray, angle_rad: float) -> np.ndarray:
    """Rotation matrix via Rodrigues formula (numpy, axis must be unit)."""
    K = np.array(
        [
            [0, -axis[2], axis[1]],
            [axis[2], 0, -axis[0]],
            [-axis[1], axis[0], 0],
        ],
        dtype=np.float64,
    )
    return np.eye(3) + math.sin(angle_rad) * K + (1 - math.cos(angle_rad)) * (K @ K)


def random_unit_axis(rng: np.random.Generator) -> np.ndarray:
    v = rng.standard_normal(3)
    return v / np.linalg.norm(v)


# ---------------------------------------------------------------------------
# Riemannian SO(3) utilities
# ---------------------------------------------------------------------------


def _project_to_tangent(R: torch.Tensor, G: torch.Tensor) -> torch.Tensor:
    """Project Euclidean gradient G (3×3) to tangent space at R ∈ SO(3).

    Returns Omega_skew ∈ so(3) (skew-symmetric 3×3).
    Formula: Omega_full = R^T @ G;  Omega_skew = (Omega_full − Omega_full^T) / 2
    """
    Omega_full = R.T @ G
    return (Omega_full - Omega_full.T) / 2.0


def _omega_to_vec(Omega: torch.Tensor) -> torch.Tensor:
    """Extract 3-vector from skew-symmetric matrix.

    Convention matches skew(): vec = [Omega[2,1], Omega[0,2], Omega[1,0]].
    """
    return torch.stack([Omega[2, 1], Omega[0, 2], Omega[1, 0]])


def _vec_to_skew(v: torch.Tensor) -> torch.Tensor:
    """3-vector → skew-symmetric matrix (same convention as skew())."""
    v1, v2, v3 = v[0], v[1], v[2]
    z = torch.zeros(1, dtype=v.dtype, device=v.device).squeeze()
    row0 = torch.stack([z, -v3, v2])
    row1 = torch.stack([v3, z, -v1])
    row2 = torch.stack([-v2, v1, z])
    return torch.stack([row0, row1, row2])


# ---------------------------------------------------------------------------
# geoopt import (shared by all geoopt-based optimizers)
# ---------------------------------------------------------------------------

_GEOOPT_AVAILABLE = False
try:
    import geoopt  # type: ignore

    _GEOOPT_AVAILABLE = True
except ImportError:
    pass


def _require_geoopt() -> None:
    if not _GEOOPT_AVAILABLE:
        raise RuntimeError("geoopt not installed. Run: uv sync --extra riemannian")


def _make_stiefel_param(R_init: np.ndarray) -> "geoopt.ManifoldParameter":
    """Create a geoopt ManifoldParameter on Stiefel(3,3) ≈ SO(3)."""
    manifold = geoopt.manifolds.Stiefel()
    return geoopt.ManifoldParameter(torch.from_numpy(R_init).float(), manifold=manifold)


# ---------------------------------------------------------------------------
# Optimizer 1: Riemannian SGD (plain — no momentum, no adaptive scaling)
# ---------------------------------------------------------------------------


def run_one_riemannian_sgd_plain(
    diff_fn,
    voxel,
    vertices,
    R_init: np.ndarray,
    scale: int,
    n_steps: int = N_STEPS,
    lr: float = LR,
) -> Dict:
    """Pure steepest descent on SO(3) via geoopt.RiemannianSGD (momentum=0).

    Update rule at each step:
        Ω = projection of ∂cost/∂R onto T_R SO(3)
        R ← R · exp(−lr · Ω)
    No moments, no adaptive learning rates. The learning rate is constant.
    """
    _require_geoopt()
    R_gt = voxel.orientation
    R = _make_stiefel_param(R_init)
    optimizer = geoopt.optim.RiemannianSGD([R], lr=lr, momentum=0.0)

    quality_hist: List[float] = []
    misori_hist: List[float] = []

    for step in range(n_steps + 1):
        if step > 0:
            optimizer.zero_grad()

        with torch.set_grad_enabled(step > 0):
            info = diff_fn.evaluate(R, vertices, phase_index=voxel.phase, scale=scale)

        quality_hist.append(float(info.quality.detach()))
        misori_hist.append(misorientation_deg(R_gt, R.detach().numpy()))

        if step > 0:
            info.cost.backward()
            optimizer.step()

    return {
        "quality_history": np.array(quality_hist),
        "misorientation_history": np.array(misori_hist),
        "R_final": R.detach().numpy(),
        "n_peaks": info.n_peaks,
    }


# ---------------------------------------------------------------------------
# Optimizer 2: Riemannian SGD with momentum
# ---------------------------------------------------------------------------


def run_one_riemannian_sgd_momentum(
    diff_fn,
    voxel,
    vertices,
    R_init: np.ndarray,
    scale: int,
    n_steps: int = N_STEPS,
    lr: float = LR,
    momentum: float = MOMENTUM,
) -> Dict:
    """SGD with heavy-ball momentum on SO(3).

    Update rule:
        m ← β·m + Ω                  (accumulate in so(3) coordinates)
        R ← R · exp(−lr · m)         (geodesic retraction)

    Under the left-trivialized flat connection on SO(3), parallel transport
    of the moment vector between steps is the identity — no rotation correction
    needed (same reasoning as in Riemannian Adam).
    """
    _require_geoopt()
    R_gt = voxel.orientation
    R = _make_stiefel_param(R_init)
    optimizer = geoopt.optim.RiemannianSGD([R], lr=lr, momentum=momentum)

    quality_hist: List[float] = []
    misori_hist: List[float] = []

    for step in range(n_steps + 1):
        if step > 0:
            optimizer.zero_grad()

        with torch.set_grad_enabled(step > 0):
            info = diff_fn.evaluate(R, vertices, phase_index=voxel.phase, scale=scale)

        quality_hist.append(float(info.quality.detach()))
        misori_hist.append(misorientation_deg(R_gt, R.detach().numpy()))

        if step > 0:
            info.cost.backward()
            optimizer.step()

    return {
        "quality_history": np.array(quality_hist),
        "misorientation_history": np.array(misori_hist),
        "R_final": R.detach().numpy(),
        "n_peaks": info.n_peaks,
    }


# ---------------------------------------------------------------------------
# Optimizer 3: Riemannian SGD with Nesterov momentum
# ---------------------------------------------------------------------------


def run_one_riemannian_sgd_nesterov(
    diff_fn,
    voxel,
    vertices,
    R_init: np.ndarray,
    scale: int,
    n_steps: int = N_STEPS,
    lr: float = LR,
    momentum: float = MOMENTUM,
) -> Dict:
    """Nesterov momentum on SO(3).

    Like SGD-momentum but evaluates gradient at the look-ahead point:
        R_la = R · exp(−lr·β·m)           (tentative step along current momentum)
        G_la = ∂cost/∂R at R_la
        Ω_la = project G_la to T_{R_la} SO(3)
        m ← β·m + Ω_la
        R ← R · exp(−lr · m)

    geoopt.RiemannianSGD implements Nesterov momentum when nesterov=True.
    """
    _require_geoopt()
    R_gt = voxel.orientation
    R = _make_stiefel_param(R_init)
    optimizer = geoopt.optim.RiemannianSGD([R], lr=lr, momentum=momentum, nesterov=True)

    quality_hist: List[float] = []
    misori_hist: List[float] = []

    for step in range(n_steps + 1):
        if step > 0:
            optimizer.zero_grad()

        with torch.set_grad_enabled(step > 0):
            info = diff_fn.evaluate(R, vertices, phase_index=voxel.phase, scale=scale)

        quality_hist.append(float(info.quality.detach()))
        misori_hist.append(misorientation_deg(R_gt, R.detach().numpy()))

        if step > 0:
            info.cost.backward()
            optimizer.step()

    return {
        "quality_history": np.array(quality_hist),
        "misorientation_history": np.array(misori_hist),
        "R_final": R.detach().numpy(),
        "n_peaks": info.n_peaks,
    }


# ---------------------------------------------------------------------------
# Optimizer 4: Riemannian SGD with cosine-annealing LR schedule
# ---------------------------------------------------------------------------


def run_one_riemannian_sgd_cosine(
    diff_fn,
    voxel,
    vertices,
    R_init: np.ndarray,
    scale: int,
    n_steps: int = N_STEPS,
    lr_max: float = LR_MAX,
    lr_min: float = LR_MIN,
    momentum: float = MOMENTUM,
) -> Dict:
    """SGD with momentum and cosine-annealing learning rate.

    lr_t = lr_min + 0.5·(lr_max − lr_min)·(1 + cos(π·t / T))

    Starts at lr_max (aggressive exploration), decays smoothly to lr_min
    (fine refinement). SGD is LR-sensitive; cosine annealing is the standard
    way to give it a fair comparison against Adam.
    """
    _require_geoopt()
    R_gt = voxel.orientation
    R = _make_stiefel_param(R_init)
    optimizer = geoopt.optim.RiemannianSGD([R], lr=lr_max, momentum=momentum)
    scheduler = torch.optim.lr_scheduler.CosineAnnealingLR(
        optimizer, T_max=max(n_steps, 1), eta_min=lr_min
    )

    quality_hist: List[float] = []
    misori_hist: List[float] = []

    for step in range(n_steps + 1):
        if step > 0:
            optimizer.zero_grad()

        with torch.set_grad_enabled(step > 0):
            info = diff_fn.evaluate(R, vertices, phase_index=voxel.phase, scale=scale)

        quality_hist.append(float(info.quality.detach()))
        misori_hist.append(misorientation_deg(R_gt, R.detach().numpy()))

        if step > 0:
            info.cost.backward()
            optimizer.step()
            scheduler.step()

    return {
        "quality_history": np.array(quality_hist),
        "misorientation_history": np.array(misori_hist),
        "R_final": R.detach().numpy(),
        "n_peaks": info.n_peaks,
    }


# ---------------------------------------------------------------------------
# Optimizer 5: Riemannian SGLD (Stochastic Gradient Langevin Dynamics)
# ---------------------------------------------------------------------------


def run_one_riemannian_sgld(
    diff_fn,
    voxel,
    vertices,
    R_init: np.ndarray,
    scale: int,
    n_steps: int = N_STEPS,
    lr: float = LR,
    t_init: float = SGLD_T_INIT,
) -> Dict:
    """Stochastic Gradient Langevin Dynamics on SO(3) — pure PyTorch.

    At each step t:
        Ω = project(R, ∂cost/∂R)              # Riemannian gradient (so(3))
        ω = _omega_to_vec(Ω)                   # 3-vector
        ε ~ N(0, I₃)                           # isotropic noise
        T_t = t_init · (1 − t/n_steps)         # linear temperature anneal
        noise_scale = √(2·lr·T_t)              # Langevin noise magnitude
        R ← R · exp(−lr·skew(ω) − noise_scale·skew(ε))

    Reference: Welling & Teh, "Bayesian Learning via Stochastic Gradient
    Langevin Dynamics", ICML 2011.

    The Langevin noise is isotropic in T_R SO(3) (3D), so it generates a
    proper Brownian motion on the manifold. As T_t → 0, SGLD reduces to
    standard Riemannian gradient descent.

    This is the most principled approach to escape the flat cost landscape:
    the same random exploration used by the original adaptive MC search, but
    guided by the gradient signal where it exists.
    """
    R_gt = voxel.orientation
    R = torch.from_numpy(R_init).float()

    quality_hist: List[float] = []
    misori_hist: List[float] = []

    for step in range(n_steps + 1):
        R_param = R.clone().requires_grad_(True)

        with torch.set_grad_enabled(step > 0):
            info = diff_fn.evaluate(R_param, vertices, phase_index=voxel.phase, scale=scale)

        quality_hist.append(float(info.quality.detach()))
        misori_hist.append(misorientation_deg(R_gt, R.numpy()))

        if step > 0:
            info.cost.backward()

            with torch.no_grad():
                G = R_param.grad  # (3,3) Euclidean gradient
                Omega = _project_to_tangent(R_param.detach(), G)
                omega_vec = _omega_to_vec(Omega)

                # Linear temperature anneal: T_t → 0 at step = n_steps
                T_t = t_init * max(0.0, 1.0 - step / n_steps)
                noise_scale = math.sqrt(2.0 * lr * T_t)

                noise_vec = torch.randn(3) * noise_scale
                step_vec = lr * omega_vec + noise_vec

                # Geodesic retraction with Langevin perturbation
                R = R_param.detach() @ torch.matrix_exp(-_vec_to_skew(step_vec))

    return {
        "quality_history": np.array(quality_hist),
        "misorientation_history": np.array(misori_hist),
        "R_final": R.numpy(),
        "n_peaks": info.n_peaks,
    }


# ---------------------------------------------------------------------------
# Optimizer 6: Riemannian Adam (manual) — Adam baseline for comparison
# ---------------------------------------------------------------------------


def run_one_riemannian_adam_manual(
    diff_fn,
    voxel,
    vertices,
    R_init: np.ndarray,
    scale: int,
    n_steps: int = N_STEPS,
    lr: float = LR,
    beta1: float = BETA1,
    beta2: float = BETA2,
    eps: float = ADAM_EPS,
) -> Dict:
    """Riemannian Adam baseline (copied from bench_riemannian_optimization.py).

    Included here so Adam can be directly compared to SGD variants in the same
    plots without merging CSV files.
    """
    R_gt = voxel.orientation
    R = torch.from_numpy(R_init).float()
    m1 = torch.zeros(3)
    m2 = torch.zeros(3)

    quality_hist: List[float] = []
    misori_hist: List[float] = []

    for step in range(n_steps + 1):
        R_param = R.clone().requires_grad_(True)

        with torch.set_grad_enabled(step > 0):
            info = diff_fn.evaluate(R_param, vertices, phase_index=voxel.phase, scale=scale)

        quality_hist.append(float(info.quality.detach()))
        misori_hist.append(misorientation_deg(R_gt, R.numpy()))

        if step > 0:
            info.cost.backward()

            with torch.no_grad():
                G = R_param.grad
                Omega_skew = _project_to_tangent(R_param.detach(), G)
                omega_vec = _omega_to_vec(Omega_skew)

                m1 = beta1 * m1 + (1 - beta1) * omega_vec
                m2 = beta2 * m2 + (1 - beta2) * omega_vec**2

                m1_hat = m1 / (1 - beta1**step)
                m2_hat = m2 / (1 - beta2**step)
                v_vec = m1_hat / (torch.sqrt(m2_hat) + eps)

                R = R_param.detach() @ torch.matrix_exp(-lr * _vec_to_skew(v_vec))

    return {
        "quality_history": np.array(quality_hist),
        "misorientation_history": np.array(misori_hist),
        "R_final": R.numpy(),
        "n_peaks": info.n_peaks,
    }


# ---------------------------------------------------------------------------
# Optimizer registry
# ---------------------------------------------------------------------------

OPTIMIZER_REGISTRY: Dict[str, object] = {
    "riemannian_adam_manual": run_one_riemannian_adam_manual,
    "riemannian_sgld": run_one_riemannian_sgld,
}
if _GEOOPT_AVAILABLE:
    OPTIMIZER_REGISTRY["riemannian_sgd_plain"] = run_one_riemannian_sgd_plain
    OPTIMIZER_REGISTRY["riemannian_sgd_momentum"] = run_one_riemannian_sgd_momentum
    OPTIMIZER_REGISTRY["riemannian_sgd_nesterov"] = run_one_riemannian_sgd_nesterov
    OPTIMIZER_REGISTRY["riemannian_sgd_cosine"] = run_one_riemannian_sgd_cosine


# ---------------------------------------------------------------------------
# Setup (identical to bench_riemannian_optimization.py)
# ---------------------------------------------------------------------------


def setup_example(example_dir: Path, basename: str, omega_windows: List[int]):
    """Load physics, cost functions, and mic for one example."""
    from icenine.config_file import ConfigFile
    from icenine.cost_functions import VoxelCostFunction
    from icenine.differentiable_cost import (
        DifferentiableCostFunction,
        MultiScaleImageStack,
        SparseImageStack,
    )
    from icenine.experiment_setup import XDMExperimentSetup
    from icenine.experimental_data import ExperimentalData
    from icenine.mic_file import MicFile
    from icenine.reconstructor import _get_voxel_vertices
    from icenine.sample import Sample
    from icenine.simulation import Simulation

    config_path = example_dir / "ConfigFiles" / "Example2.Simulation.config"
    data_dir = example_dir / "ScatteringData_Python"

    os.chdir(example_dir)

    config = ConfigFile.from_file(str(config_path))
    config.out_file_basename = basename

    exp_setup = XDMExperimentSetup(config)
    exp_setup.initialize_experiment()
    detector_list = exp_setup.get_detector_list()
    range_map = exp_setup.get_range_to_index_map()
    sample = Sample()
    exp_setup.initialize_sample(sample, detector_list[0])
    simulator = Simulation(exp_setup)
    structure_list = sample.get_structure_list()

    mic_path = Path(config.sample_filename)
    if not mic_path.is_absolute():
        mic_path = example_dir / mic_path
    mic = MicFile.read(str(mic_path))
    print(f"  Loaded .mic: {mic_path.name}  ({len(mic.voxels)} voxels)")

    print(f"  Loading SparseImageStack from {data_dir.name}/ ...")
    image_stack = SparseImageStack.from_image_directory(
        directory=str(data_dir),
        basename=basename,
        ext="d",
        serial_length=5,
        n_omega=180,
        n_detectors=2,
        num_rows=2048,
        num_cols=2048,
        binary=True,
    )
    print(f"  SparseImageStack: {image_stack.memory_bytes / 1024:.1f} KB")

    exp_data = ExperimentalData.from_image_directory(
        directory=str(data_dir),
        basename=basename,
        ext="d",
        serial_length=5,
        n_omega=180,
        n_detectors=2,
        num_rows=2048,
        num_cols=2048,
        mode="sparse",
    )

    hard_fn = VoxelCostFunction(
        simulator=simulator,
        detector_list=detector_list,
        range_map=range_map,
        exp_data=exp_data,
        sample=sample,
        structure_list=structure_list,
        mode="hard",
    )

    print("  Building shared downsampled base stacks [4x, 8x] ...")
    shared_ds = MultiScaleImageStack.build_shared_base(image_stack, [1, 4, 8])

    diff_fns: Dict[int, DifferentiableCostFunction] = {}
    for ow in omega_windows:
        print(f"  Building MultiScaleImageStack omega_window={ow} ...")
        ms = MultiScaleImageStack(
            image_stack,
            [1, 4, 8],
            omega_window=ow,
            _prebuilt_downsampled=shared_ds,
        )
        diff_fns[ow] = DifferentiableCostFunction(
            simulator=simulator,
            detector_list=detector_list,
            range_map=range_map,
            image_stack=ms,
            sample=sample,
            structure_list=structure_list,
        )

    return mic, hard_fn, diff_fns, _get_voxel_vertices


# ---------------------------------------------------------------------------
# Voxel selection (identical to bench_riemannian_optimization.py)
# ---------------------------------------------------------------------------


def select_voxels(
    mic,
    hard_fn,
    get_vertices,
    n: int,
    rng: np.random.Generator,
    threshold: float = QUALITY_THRESHOLD,
    max_scan: int = MAX_SCAN,
):
    from icenine.geometry import matrix_to_euler

    print(f"  Scanning up to {max_scan} voxels for hard quality > {threshold} ...")
    candidates = []
    for idx, voxel in enumerate(mic.voxels[:max_scan]):
        vertices = get_vertices(voxel)
        info = hard_fn.evaluate(
            orientation=voxel.orientation,
            voxel_vertices=vertices,
            phase_index=voxel.phase,
        )
        if info.quality > threshold:
            candidates.append((voxel, idx, info.quality))

    if len(candidates) < n:
        raise RuntimeError(f"Only {len(candidates)} candidates (need {n})")

    chosen_positions = rng.choice(len(candidates), size=n, replace=False)
    chosen = [candidates[i] for i in sorted(chosen_positions)]

    print(f"  Selected {n} voxels from {len(candidates)} candidates:")
    for voxel, vidx, q in chosen:
        euler = matrix_to_euler(voxel.orientation)
        print(
            f"    idx={vidx:5d}  φ1={euler[0]:7.2f}°  Φ={euler[1]:6.2f}°  "
            f"φ2={euler[2]:7.2f}°  hard_q={q:.4f}"
        )

    return [(v, vi) for v, vi, _ in chosen]


# ---------------------------------------------------------------------------
# Per-voxel sweep
# ---------------------------------------------------------------------------


def sweep_voxel_multi(
    voxel,
    vidx: int,
    diff_fns: Dict[int, object],
    get_vertices,
    perturbations_deg: List[float],
    rng: np.random.Generator,
    scales: List[int],
    omega_windows: List[int],
    optimizer_names: List[str],
    n_steps: int,
    optimizer_lrs: Optional[Dict[str, float]] = None,
) -> List[Dict]:
    """Run all (optimizer, scale, omega_window, perturbation) combos for one voxel.

    optimizer_lrs: optional per-optimizer lr override dict. Keys are optimizer names,
    values are learning rates. Optimizers not in the dict use their default lr.
    """
    vertices = get_vertices(voxel)
    rows = []

    # One random axis per perturbation — identical sequence to other benchmarks
    axes = [random_unit_axis(rng) for _ in perturbations_deg]

    for pi, (pert_deg, axis) in enumerate(zip(perturbations_deg, axes)):
        pert_rad = pert_deg * math.pi / 180.0
        R_pert = rodrigues_np(axis, pert_rad) @ voxel.orientation

        for scale in scales:
            for ow in omega_windows:
                for opt_name in optimizer_names:
                    run_fn = OPTIMIZER_REGISTRY[opt_name]
                    # Per-optimizer lr override (for fair comparison across optimizer families)
                    # cosine variant uses lr_max; all others use lr
                    extra_kwargs: Dict[str, object] = {}
                    if optimizer_lrs and opt_name in optimizer_lrs:
                        lr_val = optimizer_lrs[opt_name]
                        if opt_name == "riemannian_sgd_cosine":
                            extra_kwargs["lr_max"] = lr_val
                        else:
                            extra_kwargs["lr"] = lr_val
                    t0 = time.perf_counter()
                    result = run_fn(
                        diff_fns[ow],
                        voxel,
                        vertices,
                        R_pert,
                        scale=scale,
                        n_steps=n_steps,
                        **extra_kwargs,
                    )
                    elapsed = time.perf_counter() - t0

                    final_misori = result["misorientation_history"][-1]
                    final_quality = result["quality_history"][-1]
                    print(
                        f"    vox {vidx:5d}  pert={pert_deg:.0f}°  s={scale}  ω±{ow}"
                        f"  {opt_name:<30s}"
                        f"  misori={final_misori:.3f}°  q={final_quality:.4f}"
                        f"  ({elapsed:.0f}s)",
                        flush=True,
                    )

                    for step in range(n_steps + 1):
                        rows.append(
                            {
                                "optimizer": opt_name,
                                "voxel_idx": vidx,
                                "perturbation_deg": pert_deg,
                                "scale": scale,
                                "omega_window": ow,
                                "step": step,
                                "quality": float(result["quality_history"][step]),
                                "misorientation_deg": float(result["misorientation_history"][step]),
                            }
                        )

    return rows


# ---------------------------------------------------------------------------
# CSV save
# ---------------------------------------------------------------------------


def save_csv(rows: List[Dict], out_path: Path) -> None:
    fieldnames = [
        "optimizer",
        "voxel_idx",
        "perturbation_deg",
        "scale",
        "omega_window",
        "step",
        "quality",
        "misorientation_deg",
    ]
    with open(out_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    print(f"Saved CSV: {out_path}  ({len(rows)} rows)")


# ---------------------------------------------------------------------------
# Plots
# ---------------------------------------------------------------------------


def convergence_plot(
    csv_path: Path,
    out_quality: Path,
    out_misori: Path,
    title_prefix: str,
    perturbations_deg: List[float],
    optimizer_names: List[str],
    scale: int,
    ow: int,
    n_steps: int,
) -> None:
    """Quality and misorientation vs step — one panel per perturbation."""
    import pandas as pd

    df = pd.read_csv(csv_path)
    sub_conf = df[(df["scale"] == scale) & (df["omega_window"] == ow)]
    n_pert = len(perturbations_deg)

    for metric, ylabel, out_path in [
        ("quality", "Quality (0–1)", out_quality),
        ("misorientation_deg", "Misorientation (°)", out_misori),
    ]:
        fig, axes = plt.subplots(1, n_pert, figsize=(6 * n_pert, 5), sharey=(metric == "quality"))
        if n_pert == 1:
            axes = [axes]

        for ai, pert in enumerate(perturbations_deg):
            ax = axes[ai]
            sub = sub_conf[sub_conf["perturbation_deg"] == pert]

            for opt in optimizer_names:
                sel = sub[sub["optimizer"] == opt]
                if sel.empty:
                    continue
                mean_v = sel.groupby("step")[metric].mean()
                std_v = sel.groupby("step")[metric].std().fillna(0)
                color = OPTIMIZER_COLORS.get(opt, "gray")
                ls = OPTIMIZER_LINESTYLES.get(opt, "-")
                ax.plot(mean_v.index, mean_v.values, color=color, ls=ls, lw=2.0, label=opt)
                ax.fill_between(
                    mean_v.index, mean_v - std_v, mean_v + std_v, color=color, alpha=0.12
                )

            ax.set_xlabel("Step")
            ax.set_ylabel(ylabel if ai == 0 else "")
            ax.set_title(f"{title_prefix} — s={scale} ω±{ow}\nPerturbation = {pert:.0f}°")
            ax.legend(fontsize=7)
            ax.grid(True, alpha=0.3)
            ax.set_xlim(0, n_steps)
            if metric == "quality":
                ax.set_ylim(-0.02, 1.05)
            else:
                ax.set_ylim(bottom=0)
                ax.axhline(0, color="gray", lw=0.8, ls="--")

        fig.tight_layout()
        fig.savefig(out_path, dpi=150)
        plt.close(fig)
        print(f"Saved: {out_path}")


def summary_plot(
    csv_path: Path,
    out_path: Path,
    title_prefix: str,
    perturbations_deg: List[float],
    optimizer_names: List[str],
    scale: int,
    ow: int,
    n_steps: int,
) -> None:
    """Bar chart: final misorientation by optimizer × perturbation."""
    import pandas as pd

    df = pd.read_csv(csv_path)
    df_final = df[(df["step"] == n_steps) & (df["scale"] == scale) & (df["omega_window"] == ow)]

    n_pert = len(perturbations_deg)
    fig, axes = plt.subplots(1, n_pert, figsize=(5 * n_pert, 5), sharey=True)
    if n_pert == 1:
        axes = [axes]

    x = np.arange(len(optimizer_names))
    width = 0.6

    for ai, pert in enumerate(perturbations_deg):
        ax = axes[ai]
        sub = df_final[df_final["perturbation_deg"] == pert]

        means, stds, colors = [], [], []
        for opt in optimizer_names:
            vals = sub[sub["optimizer"] == opt]["misorientation_deg"]
            means.append(float(vals.mean()) if len(vals) > 0 else np.nan)
            stds.append(float(vals.std()) if len(vals) > 1 else 0.0)
            colors.append(OPTIMIZER_COLORS.get(opt, "gray"))

        bars = ax.bar(
            x, means, width, yerr=stds, color=colors, capsize=4, error_kw={"elinewidth": 1.2}
        )
        ax.set_xticks(x)
        ax.set_xticklabels(optimizer_names, rotation=30, ha="right", fontsize=7)
        ax.set_ylabel("Final misorientation (°)" if ai == 0 else "")
        ax.set_title(f"{title_prefix} — s={scale} ω±{ow}\nPerturbation = {pert:.0f}°")
        ax.axhline(0.1, color="green", lw=0.8, ls="--", label="0.1° target")
        ax.axhline(1.0, color="orange", lw=0.8, ls="--", label="1° threshold")
        ax.legend(fontsize=7)
        ax.grid(True, alpha=0.3, axis="y")

        for bar, mean in zip(bars, means):
            if not np.isnan(mean):
                ax.text(
                    bar.get_x() + bar.get_width() / 2,
                    mean + 0.05,
                    f"{mean:.2f}°",
                    ha="center",
                    va="bottom",
                    fontsize=7,
                )

    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    print(f"Saved: {out_path}")


# ---------------------------------------------------------------------------
# Per-example runner
# ---------------------------------------------------------------------------


def run_example(
    label: str,
    example_dir: Path,
    basename: str,
    n_voxels: Optional[int],
    scales: List[int],
    omega_windows: List[int],
    perturbations_deg: List[float],
    n_steps: int,
    optimizer_names: List[str],
    rng: np.random.Generator,
    smoke_test: bool = False,
    optimizer_lrs: Optional[Dict[str, float]] = None,
) -> None:
    print(f"\n{'='*60}")
    print(f"Example: {label}")
    print(f"geoopt available: {_GEOOPT_AVAILABLE}")
    print(f"Optimizers: {optimizer_names}")
    if optimizer_lrs:
        print(f"LR overrides: {optimizer_lrs}")
    print(f"{'='*60}")

    if smoke_test:
        n_voxels = min(n_voxels or 1, 1)
        n_steps = 5
        scales = [2]
        omega_windows = [1]
        perturbations_deg = [2.0]
        print("  [SMOKE TEST: 1 voxel, 5 steps, s=2 ω±1, pert=2°]")

    mic, hard_fn, diff_fns, get_vertices = setup_example(example_dir, basename, omega_windows)

    if n_voxels is None:
        from icenine.geometry import matrix_to_euler

        voxel_list = []
        for idx, voxel in enumerate(mic.voxels):
            vertices = get_vertices(voxel)
            info = hard_fn.evaluate(
                orientation=voxel.orientation,
                voxel_vertices=vertices,
                phase_index=voxel.phase,
            )
            if info.quality > QUALITY_THRESHOLD:
                euler = matrix_to_euler(voxel.orientation)
                voxel_list.append((voxel, idx))
                print(f"    idx={idx}  hard_q={info.quality:.4f}  [selected]")
        print(f"  Using {len(voxel_list)} voxels.")
    else:
        voxel_list = select_voxels(mic, hard_fn, get_vertices, n_voxels, rng)

    n_combos = (
        len(voxel_list)
        * len(perturbations_deg)
        * len(scales)
        * len(omega_windows)
        * len(optimizer_names)
    )
    print(
        f"\nRunning {len(voxel_list)} voxels × {len(perturbations_deg)} perts"
        f" × {len(scales)} scales × {len(omega_windows)} ow × {len(optimizer_names)} opts"
        f" = {n_combos} runs (n_steps={n_steps}) ..."
    )

    all_rows: List[Dict] = []
    t_start = time.perf_counter()

    for vi, (voxel, vidx) in enumerate(voxel_list):
        print(f"\n  Voxel {vi+1}/{len(voxel_list)}  (mic idx={vidx}) ...")
        rows = sweep_voxel_multi(
            voxel,
            vidx,
            diff_fns,
            get_vertices,
            perturbations_deg=perturbations_deg,
            rng=rng,
            scales=scales,
            omega_windows=omega_windows,
            optimizer_names=optimizer_names,
            n_steps=n_steps,
            optimizer_lrs=optimizer_lrs,
        )
        all_rows.extend(rows)

    elapsed = time.perf_counter() - t_start
    print(f"\nTotal time: {elapsed/60:.1f} min  ({len(all_rows)} rows)")

    base_tag = label.lower().replace(" ", "_").replace(".", "")
    # Embed sgd lr in tag when using per-optimizer lrs, so fair-run outputs are separate
    if optimizer_lrs:
        sgd_lr_vals = [v for k, v in optimizer_lrs.items() if k != "riemannian_adam_manual"]
        if sgd_lr_vals and len(set(sgd_lr_vals)) == 1:
            lr_str = f"_sgdlr{sgd_lr_vals[0]:.0e}".replace("-0", "").replace("+0", "")
        else:
            lr_str = ""
        tag = f"{base_tag}{lr_str}"
    else:
        tag = base_tag
    csv_path = benchmark_dir / f"sgd_opt_{tag}.csv"
    save_csv(all_rows, csv_path)

    # Plot for best omega_window (ow=1 or first) and scale=2 (or last)
    plot_scale = scales[-1]
    plot_ow = omega_windows[1] if len(omega_windows) > 1 else omega_windows[0]
    convergence_plot(
        csv_path,
        benchmark_dir / f"sgd_opt_convergence_{tag}.png",
        benchmark_dir / f"sgd_opt_misori_{tag}.png",
        title_prefix=label,
        perturbations_deg=perturbations_deg,
        optimizer_names=optimizer_names,
        scale=plot_scale,
        ow=plot_ow,
        n_steps=n_steps,
    )
    summary_plot(
        csv_path,
        benchmark_dir / f"sgd_opt_summary_{tag}.png",
        title_prefix=label,
        perturbations_deg=perturbations_deg,
        optimizer_names=optimizer_names,
        scale=plot_scale,
        ow=plot_ow,
        n_steps=n_steps,
    )


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Riemannian SGD on SO(3) orientation optimization benchmark"
    )
    parser.add_argument(
        "--smoke-test", action="store_true", help="Quick sanity: 1 voxel, 5 steps, s=2 ω±1, pert=2°"
    )
    parser.add_argument("--example", choices=["threevoxels", "manygrains", "both"], default="both")
    parser.add_argument(
        "--optimizer",
        choices=["all"]
        + list(
            {
                "riemannian_adam_manual",
                "riemannian_sgld",
                "riemannian_sgd_plain",
                "riemannian_sgd_momentum",
                "riemannian_sgd_nesterov",
                "riemannian_sgd_cosine",
            }
        ),
        default="all",
        help="Which optimizer(s) to run (default: all available)",
    )
    parser.add_argument("--n-steps", type=int, default=N_STEPS)
    parser.add_argument(
        "--sgld-temp",
        type=float,
        default=SGLD_T_INIT,
        help=f"SGLD initial temperature (default: {SGLD_T_INIT})",
    )
    parser.add_argument(
        "--sgd-lr",
        type=float,
        default=None,
        help="LR for SGD variants (plain/momentum/nesterov/sgld). "
        "Adam always uses --adam-lr. Use this for fair comparison "
        "(SGD needs ~10x smaller lr than Adam). Default: same as Adam (0.01).",
    )
    parser.add_argument(
        "--adam-lr",
        type=float,
        default=LR,
        help=f"LR for Adam baseline (default: {LR})",
    )
    args = parser.parse_args()

    # Allow sgld-temp override
    if args.sgld_temp != SGLD_T_INIT:
        SGLD_T_INIT = args.sgld_temp

    if args.optimizer == "all":
        optimizer_names = list(OPTIMIZER_REGISTRY.keys())
    else:
        if args.optimizer not in OPTIMIZER_REGISTRY:
            parser.error(
                f"Optimizer '{args.optimizer}' not available. "
                f"Available: {list(OPTIMIZER_REGISTRY.keys())}"
            )
        optimizer_names = [args.optimizer]

    print(f"Available optimizers: {list(OPTIMIZER_REGISTRY.keys())}")
    print(f"Running: {optimizer_names}")

    # Build per-optimizer lr dict for fair comparison
    sgd_variant_names = [
        "riemannian_sgd_plain",
        "riemannian_sgd_momentum",
        "riemannian_sgd_nesterov",
        "riemannian_sgd_cosine",
        "riemannian_sgld",
    ]
    optimizer_lrs: Dict[str, float] = {}
    optimizer_lrs["riemannian_adam_manual"] = args.adam_lr
    if args.sgd_lr is not None:
        for name in sgd_variant_names:
            optimizer_lrs[name] = args.sgd_lr
        print(f"LR: adam={args.adam_lr}, sgd_variants={args.sgd_lr}")
    else:
        print(f"LR: {args.adam_lr} (all optimizers, uniform)")

    rng = np.random.default_rng(seed=SEED)

    three_dir = project_root / "Examples" / "Example2.ThreeVoxels"
    many_dir = project_root / "Examples" / "Example2.ManyGrains"

    if args.example in ("threevoxels", "both"):
        run_example(
            label="ThreeVoxels",
            example_dir=three_dir,
            basename="3Grains.sim",
            n_voxels=None,
            scales=SCALES,
            omega_windows=OMEGA_WINDOWS,
            perturbations_deg=PERTURBATIONS_DEG,
            n_steps=args.n_steps,
            optimizer_names=optimizer_names,
            rng=rng,
            smoke_test=args.smoke_test,
            optimizer_lrs=optimizer_lrs if optimizer_lrs else None,
        )

    if args.example in ("manygrains", "both"):
        run_example(
            label="ManyGrains",
            example_dir=many_dir,
            basename="500Grains.sim",
            n_voxels=N_VOXELS_MANY,
            scales=SCALES,
            omega_windows=OMEGA_WINDOWS,
            perturbations_deg=PERTURBATIONS_DEG,
            n_steps=args.n_steps,
            optimizer_names=optimizer_names,
            rng=rng,
            smoke_test=args.smoke_test,
            optimizer_lrs=optimizer_lrs if optimizer_lrs else None,
        )
