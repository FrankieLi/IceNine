#!/usr/bin/env python3
"""
Comprehensive hyperparameter sweep — gradient methods vs. Monte Carlo optimizer.

Sweeps all meaningful hyperparameters for each optimizer family and benchmarks them
head-to-head against the existing MCOptimizer, recording timing, evaluation counts,
peak memory, starting orientation, ground truth, and final misorientation.

Also records subsampled optimization trajectories (angular step size and misorientation
from ground truth at every TRAJ_SUBSAMPLE gradient steps; every accepted MC move).

Optimizers and HP grids:
  riemannian_adam_geoopt   — lr × n_steps × beta1 = 42 configs
  riemannian_adam_manual   — lr × n_steps × beta1 = 42 configs (same grid)
  riemannian_sgd_plain     — lr × n_steps = 15 configs
  riemannian_sgd_momentum  — lr × momentum (n_steps=200) = 15 configs
  riemannian_sgld          — lr × T_init × n_steps = 45 configs
  mc_optimizer             — max_mc_steps × restarts × angular_step_frac = 36 configs

All gradient runs use scale=2, omega_window=1.
100 voxels for ManyGrains (5× the previous 20).

Usage:
  cd /Users/sfli/Research/IceNine/icenine_py
  uv sync --extra riemannian
  uv run python benchmarks/bench_hp_sweep.py --smoke-test --example threevoxels
  uv run python benchmarks/bench_hp_sweep.py --example threevoxels --optimizer gradient
  uv run python benchmarks/bench_hp_sweep.py --example threevoxels --optimizer mc
  uv run python benchmarks/bench_hp_sweep.py --example manygrains

Outputs (icenine_py/benchmarks/):
  hp_sweep_{example}.csv               — one row per run
  hp_sweep_trajectory_{example}.csv    — subsampled trajectory rows (run_id foreign key)
  hp_sweep_lr_sensitivity_{example}.png
  hp_sweep_nsteps_sensitivity_{example}.png
  hp_sweep_optimizer_comparison_{example}.png
  hp_sweep_trajectory_{example}.png
"""

import argparse
import csv
import math
import os
import sys
import time
import tracemalloc
from itertools import product
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
import torch
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

project_root = Path(__file__).parent.parent.parent
sys.path.insert(0, str(project_root / "icenine_py"))

benchmark_dir = Path(__file__).parent

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

FIXED_SCALE = 2
FIXED_OW = 1
PERTURBATIONS_DEG = [1.0, 2.0, 5.0]
SEED = 42
N_VOXELS_MANY = 100
QUALITY_THRESHOLD = 0.1
MAX_SCAN = 500
TRAJ_SUBSAMPLE = 10   # record gradient trajectory every N steps

# geoopt availability
_GEOOPT_AVAILABLE = False
try:
    import geoopt  # type: ignore
    _GEOOPT_AVAILABLE = True
except ImportError:
    pass


def _require_geoopt() -> None:
    if not _GEOOPT_AVAILABLE:
        raise RuntimeError("geoopt not installed. Run: uv sync --extra riemannian")


# ---------------------------------------------------------------------------
# Memory measurement
# ---------------------------------------------------------------------------

try:
    import psutil as _psutil
    _PSUTIL_AVAILABLE = True
except ImportError:
    _PSUTIL_AVAILABLE = False


def measure_run(fn, *args, **kwargs) -> Tuple[Any, float, float, float]:
    """Run fn(*args, **kwargs), returning (result, wall_sec, peak_heap_mb, rss_delta_mb)."""
    rss_before = 0.0
    if _PSUTIL_AVAILABLE:
        proc = _psutil.Process(os.getpid())
        rss_before = proc.memory_info().rss / 1024 ** 2

    tracemalloc.start()
    t0 = time.perf_counter()
    result = fn(*args, **kwargs)
    elapsed = time.perf_counter() - t0
    _, peak_bytes = tracemalloc.get_traced_memory()
    tracemalloc.stop()

    rss_delta = 0.0
    if _PSUTIL_AVAILABLE:
        rss_delta = proc.memory_info().rss / 1024 ** 2 - rss_before  # type: ignore[union-attr]

    return result, elapsed, peak_bytes / 1024 ** 2, rss_delta


# ---------------------------------------------------------------------------
# SO(3) helpers (identical to bench_riemannian_optimization.py)
# ---------------------------------------------------------------------------

def skew(theta: torch.Tensor) -> torch.Tensor:
    """3-vector → 3×3 skew-symmetric matrix."""
    assert theta.shape == (3,)
    z = torch.zeros(1, dtype=theta.dtype, device=theta.device)
    row0 = torch.stack([z.squeeze(), -theta[2],  theta[1]])
    row1 = torch.stack([theta[2],    z.squeeze(), -theta[0]])
    row2 = torch.stack([-theta[1],   theta[0],   z.squeeze()])
    return torch.stack([row0, row1, row2])


def misorientation_deg(R1: np.ndarray, R2: np.ndarray) -> float:
    """Geodesic distance on SO(3) in degrees."""
    M = R1.T @ R2
    trace = np.trace(M)
    cos_angle = np.clip((trace - 1.0) / 2.0, -1.0, 1.0)
    return float(np.degrees(np.arccos(cos_angle)))


def rodrigues_np(axis: np.ndarray, angle_rad: float) -> np.ndarray:
    """Rotation matrix via Rodrigues formula (numpy, axis must be unit)."""
    K = np.array([
        [0,        -axis[2],  axis[1]],
        [axis[2],   0,       -axis[0]],
        [-axis[1],  axis[0],  0      ],
    ], dtype=np.float64)
    return np.eye(3) + math.sin(angle_rad) * K + (1 - math.cos(angle_rad)) * (K @ K)


def random_unit_axis(rng: np.random.Generator) -> np.ndarray:
    v = rng.standard_normal(3)
    return v / np.linalg.norm(v)


def _project_to_tangent(R: torch.Tensor, G: torch.Tensor) -> torch.Tensor:
    Omega_full = R.T @ G
    return (Omega_full - Omega_full.T) / 2.0


def _omega_to_vec(Omega: torch.Tensor) -> torch.Tensor:
    return torch.stack([Omega[2, 1], Omega[0, 2], Omega[1, 0]])


def _vec_to_skew(v: torch.Tensor) -> torch.Tensor:
    v1, v2, v3 = v[0], v[1], v[2]
    z = torch.zeros(1, dtype=v.dtype, device=v.device).squeeze()
    row0 = torch.stack([z, -v3,  v2])
    row1 = torch.stack([v3,  z, -v1])
    row2 = torch.stack([-v2, v1,  z])
    return torch.stack([row0, row1, row2])


def _make_stiefel_param(R_init: np.ndarray) -> "geoopt.ManifoldParameter":
    manifold = geoopt.manifolds.Stiefel()
    return geoopt.ManifoldParameter(torch.from_numpy(R_init).float(), manifold=manifold)


def _make_traj_entry(
    event_idx: int, step: int, event_type: str,
    angular_step_deg: float, misori_gt_deg: float,
    quality: float, cur_step_rad: float,
) -> Dict:
    return {
        "event_idx": event_idx,
        "step": step,
        "event_type": event_type,
        "angular_step_deg": angular_step_deg,
        "misori_gt_deg": misori_gt_deg,
        "quality": quality,
        "cur_step_rad": cur_step_rad,
    }


# ---------------------------------------------------------------------------
# Gradient run_one_* functions — HP kwargs + trajectory (no full history)
# ---------------------------------------------------------------------------

def run_one_riemannian_adam_geoopt(
    diff_fn, voxel, vertices, R_init: np.ndarray, scale: int,
    n_steps: int, lr: float, beta1: float = 0.9, beta2: float = 0.999,
    traj_subsample: int = TRAJ_SUBSAMPLE,
) -> Dict:
    """Riemannian Adam via geoopt.RiemannianAdam — HP-sweep variant (no history arrays)."""
    _require_geoopt()
    R_gt = voxel.orientation
    R = _make_stiefel_param(R_init)
    optimizer = geoopt.optim.RiemannianAdam([R], lr=lr, betas=(beta1, beta2))

    R_prev_traj = R_init.copy()
    traj: List[Dict] = []
    event_idx = 0
    n_evals = 0
    info = None

    for step in range(n_steps + 1):
        if step > 0:
            optimizer.zero_grad()
        with torch.set_grad_enabled(step > 0):
            info = diff_fn.evaluate(R, vertices, phase_index=voxel.phase, scale=scale)
        n_evals += 1

        if traj_subsample > 0 and step % traj_subsample == 0:
            R_np = R.detach().numpy()
            traj.append(_make_traj_entry(
                event_idx, step, "grad_step",
                misorientation_deg(R_prev_traj, R_np),
                misorientation_deg(R_gt, R_np),
                float(info.quality.detach()),
                float("nan"),
            ))
            R_prev_traj = R_np.copy()
            event_idx += 1

        if step > 0 and info.cost.requires_grad:
            info.cost.backward()
            optimizer.step()

    R_final = R.detach().numpy()
    return {
        "R_final": R_final,
        "final_quality": float(info.quality.detach()),
        "final_misori_deg": misorientation_deg(R_gt, R_final),
        "n_evaluations": n_evals,
        "n_peaks": int(info.n_peaks) if info is not None else 0,
        "trajectory": traj,
    }


def run_one_riemannian_adam_manual(
    diff_fn, voxel, vertices, R_init: np.ndarray, scale: int,
    n_steps: int, lr: float, beta1: float = 0.9, beta2: float = 0.999,
    eps: float = 1e-8, traj_subsample: int = TRAJ_SUBSAMPLE,
) -> Dict:
    """Manual Riemannian Adam on SO(3) — HP-sweep variant."""
    R_gt = voxel.orientation
    R = torch.from_numpy(R_init).float()
    m1 = torch.zeros(3)
    m2 = torch.zeros(3)

    R_prev_traj = R_init.copy()
    traj: List[Dict] = []
    event_idx = 0
    n_evals = 0
    info = None

    for step in range(n_steps + 1):
        R_param = R.clone().requires_grad_(True)
        with torch.set_grad_enabled(step > 0):
            info = diff_fn.evaluate(R_param, vertices, phase_index=voxel.phase, scale=scale)
        n_evals += 1

        if traj_subsample > 0 and step % traj_subsample == 0:
            R_np = R.detach().numpy()
            traj.append(_make_traj_entry(
                event_idx, step, "grad_step",
                misorientation_deg(R_prev_traj, R_np),
                misorientation_deg(R_gt, R_np),
                float(info.quality.detach()),
                float("nan"),
            ))
            R_prev_traj = R_np.copy()
            event_idx += 1

        if step > 0 and info.cost.requires_grad:
            info.cost.backward()
            with torch.no_grad():
                G = R_param.grad
                Omega_skew = _project_to_tangent(R_param.detach(), G)
                omega_vec = _omega_to_vec(Omega_skew)
                m1 = beta1 * m1 + (1 - beta1) * omega_vec
                m2 = beta2 * m2 + (1 - beta2) * omega_vec ** 2
                m1_hat = m1 / (1 - beta1 ** step)
                m2_hat = m2 / (1 - beta2 ** step)
                v_vec = m1_hat / (torch.sqrt(m2_hat) + eps)
                R = R_param.detach() @ torch.matrix_exp(-lr * _vec_to_skew(v_vec))

    R_final = R.detach().numpy()
    return {
        "R_final": R_final,
        "final_quality": float(info.quality.detach()),
        "final_misori_deg": misorientation_deg(R_gt, R_final),
        "n_evaluations": n_evals,
        "n_peaks": int(info.n_peaks) if info is not None else 0,
        "trajectory": traj,
    }


def run_one_riemannian_sgd_plain(
    diff_fn, voxel, vertices, R_init: np.ndarray, scale: int,
    n_steps: int, lr: float, traj_subsample: int = TRAJ_SUBSAMPLE,
) -> Dict:
    """Plain Riemannian SGD on SO(3) — HP-sweep variant."""
    _require_geoopt()
    R_gt = voxel.orientation
    R = _make_stiefel_param(R_init)
    optimizer = geoopt.optim.RiemannianSGD([R], lr=lr, momentum=0.0)

    R_prev_traj = R_init.copy()
    traj: List[Dict] = []
    event_idx = 0
    n_evals = 0
    info = None

    for step in range(n_steps + 1):
        if step > 0:
            optimizer.zero_grad()
        with torch.set_grad_enabled(step > 0):
            info = diff_fn.evaluate(R, vertices, phase_index=voxel.phase, scale=scale)
        n_evals += 1

        if traj_subsample > 0 and step % traj_subsample == 0:
            R_np = R.detach().numpy()
            traj.append(_make_traj_entry(
                event_idx, step, "grad_step",
                misorientation_deg(R_prev_traj, R_np),
                misorientation_deg(R_gt, R_np),
                float(info.quality.detach()),
                float("nan"),
            ))
            R_prev_traj = R_np.copy()
            event_idx += 1

        if step > 0 and info.cost.requires_grad:
            info.cost.backward()
            optimizer.step()

    R_final = R.detach().numpy()
    return {
        "R_final": R_final,
        "final_quality": float(info.quality.detach()),
        "final_misori_deg": misorientation_deg(R_gt, R_final),
        "n_evaluations": n_evals,
        "n_peaks": int(info.n_peaks) if info is not None else 0,
        "trajectory": traj,
    }


def run_one_riemannian_sgd_momentum(
    diff_fn, voxel, vertices, R_init: np.ndarray, scale: int,
    n_steps: int, lr: float, momentum: float = 0.9,
    traj_subsample: int = TRAJ_SUBSAMPLE,
) -> Dict:
    """Riemannian SGD with momentum on SO(3) — HP-sweep variant."""
    _require_geoopt()
    R_gt = voxel.orientation
    R = _make_stiefel_param(R_init)
    optimizer = geoopt.optim.RiemannianSGD([R], lr=lr, momentum=momentum)

    R_prev_traj = R_init.copy()
    traj: List[Dict] = []
    event_idx = 0
    n_evals = 0
    info = None

    for step in range(n_steps + 1):
        if step > 0:
            optimizer.zero_grad()
        with torch.set_grad_enabled(step > 0):
            info = diff_fn.evaluate(R, vertices, phase_index=voxel.phase, scale=scale)
        n_evals += 1

        if traj_subsample > 0 and step % traj_subsample == 0:
            R_np = R.detach().numpy()
            traj.append(_make_traj_entry(
                event_idx, step, "grad_step",
                misorientation_deg(R_prev_traj, R_np),
                misorientation_deg(R_gt, R_np),
                float(info.quality.detach()),
                float("nan"),
            ))
            R_prev_traj = R_np.copy()
            event_idx += 1

        if step > 0 and info.cost.requires_grad:
            info.cost.backward()
            optimizer.step()

    R_final = R.detach().numpy()
    return {
        "R_final": R_final,
        "final_quality": float(info.quality.detach()),
        "final_misori_deg": misorientation_deg(R_gt, R_final),
        "n_evaluations": n_evals,
        "n_peaks": int(info.n_peaks) if info is not None else 0,
        "trajectory": traj,
    }


def run_one_riemannian_sgld(
    diff_fn, voxel, vertices, R_init: np.ndarray, scale: int,
    n_steps: int, lr: float, t_init: float = 0.01,
    traj_subsample: int = TRAJ_SUBSAMPLE,
) -> Dict:
    """Stochastic Gradient Langevin Dynamics on SO(3) — HP-sweep variant."""
    R_gt = voxel.orientation
    R = torch.from_numpy(R_init).float()

    R_prev_traj = R_init.copy()
    traj: List[Dict] = []
    event_idx = 0
    n_evals = 0
    info = None

    for step in range(n_steps + 1):
        R_param = R.clone().requires_grad_(True)
        with torch.set_grad_enabled(step > 0):
            info = diff_fn.evaluate(R_param, vertices, phase_index=voxel.phase, scale=scale)
        n_evals += 1

        if traj_subsample > 0 and step % traj_subsample == 0:
            R_np = R.detach().numpy()
            traj.append(_make_traj_entry(
                event_idx, step, "grad_step",
                misorientation_deg(R_prev_traj, R_np),
                misorientation_deg(R_gt, R_np),
                float(info.quality.detach()),
                float("nan"),
            ))
            R_prev_traj = R_np.copy()
            event_idx += 1

        if step > 0 and info.cost.requires_grad:
            info.cost.backward()
            with torch.no_grad():
                G = R_param.grad
                Omega = _project_to_tangent(R_param.detach(), G)
                omega_vec = _omega_to_vec(Omega)
                T_t = t_init * max(0.0, 1.0 - step / n_steps)
                noise_scale = math.sqrt(2.0 * lr * T_t)
                noise_vec = torch.randn(3) * noise_scale
                step_vec = lr * omega_vec + noise_vec
                R = R_param.detach() @ torch.matrix_exp(-_vec_to_skew(step_vec))

    R_final = R.detach().numpy()
    return {
        "R_final": R_final,
        "final_quality": float(info.quality.detach()),
        "final_misori_deg": misorientation_deg(R_gt, R_final),
        "n_evaluations": n_evals,
        "n_peaks": int(info.n_peaks) if info is not None else 0,
        "trajectory": traj,
    }


# ---------------------------------------------------------------------------
# MC run_one function
# ---------------------------------------------------------------------------

def run_one_mc(
    hard_fn,
    voxel,
    vertices,
    R_init: np.ndarray,
    R_gt: np.ndarray,
    max_mc_steps: int,
    successive_restarts: int,
    angular_step_frac: float,
    perturbation_rad: float,
    record_traj: bool = True,
) -> Dict:
    """Run MCOptimizer and return summary dict with trajectory."""
    from icenine.orientation_search import MCOptimizer

    mc = MCOptimizer(hard_fn, vertices, phase_index=voxel.phase)
    hard_fn.eval_count = 0
    angular_box = perturbation_rad * 1.5
    angular_step = angular_box * angular_step_frac
    traj: Optional[List[Dict]] = [] if record_traj else None

    result = mc.optimize(
        initial_orientation=R_init,
        angular_box_side=angular_box,
        angular_step=angular_step,
        max_mc_steps=max_mc_steps,
        max_restarts=successive_restarts,
        max_convergence_cost=0.0,
        trajectory=traj,
    )

    R_final = result.orientation
    n_peaks = 0
    if result.overlap_info is not None:
        n_peaks = int(getattr(result.overlap_info, "n_peaks", 0))

    # Enrich MC trajectory with event_idx, misori_gt_deg, and quality fields.
    # misori_gt_deg = NaN: we don't reconstruct R at each accepted move (only the final is known).
    # quality = NaN: not tracked per-accepted-move (only the final is known).
    if traj is not None:
        for i, rec in enumerate(traj):
            rec["event_idx"] = i
            rec["misori_gt_deg"] = float("nan")
            rec["quality"] = float("nan")

    return {
        "R_final": R_final,
        "final_quality": 1.0 - result.cost,
        "final_misori_deg": misorientation_deg(R_gt, R_final),
        "n_evaluations": hard_fn.eval_count,
        "n_peaks": n_peaks,
        "trajectory": traj if traj is not None else [],
    }


# ---------------------------------------------------------------------------
# HP grid builder
# ---------------------------------------------------------------------------

def build_hp_grids() -> Dict[str, List[Dict]]:
    """Return the full HP grid for each optimizer as a list of config dicts."""
    grids: Dict[str, List[Dict]] = {}

    # riemannian_adam_geoopt / riemannian_adam_manual: same grid (42 configs)
    adam_configs = [
        {"hp_id": i, "lr": lr, "n_steps": ns, "beta1": b1, "beta2": 0.999}
        for i, (lr, ns, b1) in enumerate(product(
            [1e-4, 5e-4, 1e-3, 5e-3, 0.01, 0.05, 0.1],
            [100, 200, 500],
            [0.9, 0.95],
        ))
    ]
    grids["riemannian_adam_geoopt"] = adam_configs
    grids["riemannian_adam_manual"] = [dict(c) for c in adam_configs]

    # riemannian_sgd_plain: 15 configs
    grids["riemannian_sgd_plain"] = [
        {"hp_id": i, "lr": lr, "n_steps": ns}
        for i, (lr, ns) in enumerate(product(
            [1e-4, 5e-4, 1e-3, 5e-3, 0.01],
            [100, 200, 500],
        ))
    ]

    # riemannian_sgd_momentum: 15 configs (n_steps=200 fixed)
    grids["riemannian_sgd_momentum"] = [
        {"hp_id": i, "lr": lr, "n_steps": 200, "momentum": m}
        for i, (lr, m) in enumerate(product(
            [1e-5, 5e-5, 1e-4, 5e-4, 1e-3],
            [0.5, 0.9, 0.99],
        ))
    ]

    # riemannian_sgld: 45 configs
    grids["riemannian_sgld"] = [
        {"hp_id": i, "lr": lr, "n_steps": ns, "t_init": T}
        for i, (lr, T, ns) in enumerate(product(
            [1e-4, 5e-4, 1e-3, 5e-3, 0.01],
            [0.001, 0.01, 0.1],
            [100, 200, 500],
        ))
    ]

    # mc_optimizer: 36 configs
    grids["mc_optimizer"] = [
        {
            "hp_id": i,
            "max_mc_steps": ms,
            "successive_restarts": sr,
            "angular_step_frac": asf,
        }
        for i, (ms, sr, asf) in enumerate(product(
            [100, 500, 1000, 3500],
            [0, 2, 5],
            [0.25, 0.5, 1.0],
        ))
    ]

    return grids


# ---------------------------------------------------------------------------
# HP grid → run function dispatch
# ---------------------------------------------------------------------------

_GRADIENT_OPTIMIZERS = {
    "riemannian_adam_geoopt": run_one_riemannian_adam_geoopt,
    "riemannian_adam_manual": run_one_riemannian_adam_manual,
    "riemannian_sgd_plain":   run_one_riemannian_sgd_plain,
    "riemannian_sgd_momentum": run_one_riemannian_sgd_momentum,
    "riemannian_sgld":        run_one_riemannian_sgld,
}


def _hp_kwargs_for_gradient(optimizer_name: str, hp_cfg: Dict) -> Dict:
    """Extract optimizer-specific kwargs from HP config dict."""
    base = {k: v for k, v in hp_cfg.items() if k not in ("hp_id",)}
    return base


# ---------------------------------------------------------------------------
# Setup (identical to bench_sgd_optimization.py)
# ---------------------------------------------------------------------------

def setup_example(example_dir: Path, basename: str):
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
        basename=basename, ext="d", serial_length=5,
        n_omega=180, n_detectors=2, num_rows=2048, num_cols=2048, binary=True,
    )
    print(f"  SparseImageStack: {image_stack.memory_bytes / 1024:.1f} KB")

    exp_data = ExperimentalData.from_image_directory(
        directory=str(data_dir),
        basename=basename, ext="d", serial_length=5,
        n_omega=180, n_detectors=2, num_rows=2048, num_cols=2048, mode="sparse",
    )

    hard_fn = VoxelCostFunction(
        simulator=simulator, detector_list=detector_list, range_map=range_map,
        exp_data=exp_data, sample=sample, structure_list=structure_list, mode="hard",
    )

    print(f"  Building MultiScaleImageStack omega_window={FIXED_OW} ...")
    shared_ds = MultiScaleImageStack.build_shared_base(image_stack, [1, 4, 8])
    ms = MultiScaleImageStack(
        image_stack, [1, 4, 8], omega_window=FIXED_OW,
        _prebuilt_downsampled=shared_ds,
    )
    diff_fn = DifferentiableCostFunction(
        simulator=simulator, detector_list=detector_list, range_map=range_map,
        image_stack=ms, sample=sample, structure_list=structure_list,
    )

    return mic, hard_fn, diff_fn, _get_voxel_vertices


# ---------------------------------------------------------------------------
# Voxel selection
# ---------------------------------------------------------------------------

def select_voxels(mic, hard_fn, get_vertices, n: int, rng: np.random.Generator,
                  threshold: float = QUALITY_THRESHOLD, max_scan: int = MAX_SCAN):
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
        print(f"    idx={vidx:5d}  φ1={euler[0]:7.2f}°  Φ={euler[1]:6.2f}°  "
              f"φ2={euler[2]:7.2f}°  hard_q={q:.4f}")

    return [(v, vi) for v, vi, _ in chosen]


# ---------------------------------------------------------------------------
# Incremental CSV writers
# ---------------------------------------------------------------------------

MAIN_FIELDNAMES = [
    "run_id", "optimizer", "hp_id", "lr", "n_steps", "n_evaluations",
    "beta1", "beta2", "momentum", "sgld_temp", "angular_step_frac",
    "successive_restarts", "scale", "omega_window", "perturbation_deg",
    "voxel_idx",
    "R_start_phi1", "R_start_Phi", "R_start_phi2",
    "R_gt_phi1", "R_gt_Phi", "R_gt_phi2",
    "final_misorientation_deg", "final_quality", "wall_time_sec",
    "n_peaks", "peak_memory_mb", "process_rss_delta_mb",
]

TRAJ_FIELDNAMES = [
    "run_id", "event_idx", "step", "event_type",
    "angular_step_deg", "misori_gt_deg", "quality", "cur_step_rad",
]


def open_incremental_csv(path: Path, fieldnames: List[str]):
    """Open a CSV for incremental writing; write header only if file is new."""
    is_new = not path.exists()
    f = open(path, "a", newline="")
    writer = csv.DictWriter(f, fieldnames=fieldnames, extrasaction="ignore")
    if is_new:
        writer.writeheader()
    return f, writer


# ---------------------------------------------------------------------------
# Per-run result-to-row helper
# ---------------------------------------------------------------------------

def _result_to_main_row(
    run_id: int,
    optimizer: str,
    hp_cfg: Dict,
    R_start: np.ndarray,
    R_gt: np.ndarray,
    result: Dict,
    wall_time: float,
    peak_mem: float,
    rss_delta: float,
    voxel_idx: int,
    perturbation_deg: float,
) -> Dict:
    from icenine.geometry import matrix_to_euler

    euler_start = matrix_to_euler(R_start)
    euler_gt = matrix_to_euler(R_gt)

    return {
        "run_id": run_id,
        "optimizer": optimizer,
        "hp_id": hp_cfg.get("hp_id", 0),
        "lr": hp_cfg.get("lr", float("nan")),
        "n_steps": hp_cfg.get("n_steps", hp_cfg.get("max_mc_steps", float("nan"))),
        "n_evaluations": result["n_evaluations"],
        "beta1": hp_cfg.get("beta1", float("nan")),
        "beta2": hp_cfg.get("beta2", float("nan")),
        "momentum": hp_cfg.get("momentum", float("nan")),
        "sgld_temp": hp_cfg.get("t_init", float("nan")),
        "angular_step_frac": hp_cfg.get("angular_step_frac", float("nan")),
        "successive_restarts": hp_cfg.get("successive_restarts", float("nan")),
        "scale": FIXED_SCALE if optimizer != "mc_optimizer" else float("nan"),
        "omega_window": FIXED_OW if optimizer != "mc_optimizer" else float("nan"),
        "perturbation_deg": perturbation_deg,
        "voxel_idx": voxel_idx,
        "R_start_phi1": float(euler_start[0]),
        "R_start_Phi": float(euler_start[1]),
        "R_start_phi2": float(euler_start[2]),
        "R_gt_phi1": float(euler_gt[0]),
        "R_gt_Phi": float(euler_gt[1]),
        "R_gt_phi2": float(euler_gt[2]),
        "final_misorientation_deg": result["final_misori_deg"],
        "final_quality": result["final_quality"],
        "wall_time_sec": wall_time,
        "n_peaks": result["n_peaks"],
        "peak_memory_mb": peak_mem,
        "process_rss_delta_mb": rss_delta,
    }


# ---------------------------------------------------------------------------
# Main sweep function
# ---------------------------------------------------------------------------

def run_hp_sweep(
    label: str,
    example_dir: Path,
    basename: str,
    voxel_list: List[Tuple],
    hard_fn,
    diff_fn,
    get_vertices,
    hp_grids: Dict[str, List[Dict]],
    optimizers_to_run: List[str],
    perturbations_deg: List[float],
    rng: np.random.Generator,
    main_csv_path: Path,
    traj_csv_path: Path,
    run_id_start: int = 0,
    smoke_test: bool = False,
) -> int:
    """Run the HP sweep and write results incrementally.

    Returns the final run_id (for continuation).
    """
    main_f, main_writer = open_incremental_csv(main_csv_path, MAIN_FIELDNAMES)
    traj_f, traj_writer = open_incremental_csv(traj_csv_path, TRAJ_FIELDNAMES)

    run_id = run_id_start
    t_start = time.perf_counter()
    total_runs = 0

    try:
        for optimizer_name in optimizers_to_run:
            hp_grid = hp_grids[optimizer_name]
            is_mc = optimizer_name == "mc_optimizer"

            print(f"\n{'='*60}")
            print(f"  Optimizer: {optimizer_name}  ({len(hp_grid)} HP configs)")
            print(f"{'='*60}")

            for hp_cfg in hp_grid:
                hp_id = hp_cfg["hp_id"]

                for vi, (voxel, vidx) in enumerate(voxel_list):
                    vertices = get_vertices(voxel)
                    R_gt = voxel.orientation

                    # One random axis per perturbation, fixed per (voxel, voxel_position)
                    # Use same SEED+vidx pattern for reproducibility
                    voxel_rng = np.random.default_rng(SEED + vidx)
                    axes = [random_unit_axis(voxel_rng) for _ in perturbations_deg]

                    for pi, (pert_deg, axis) in enumerate(zip(perturbations_deg, axes)):
                        pert_rad = pert_deg * math.pi / 180.0
                        R_start = rodrigues_np(axis, pert_rad) @ R_gt

                        if is_mc:
                            run_kwargs = dict(
                                hard_fn=hard_fn,
                                voxel=voxel,
                                vertices=vertices,
                                R_init=R_start,
                                R_gt=R_gt,
                                max_mc_steps=hp_cfg["max_mc_steps"],
                                successive_restarts=hp_cfg["successive_restarts"],
                                angular_step_frac=hp_cfg["angular_step_frac"],
                                perturbation_rad=pert_rad,
                            )
                            run_fn = run_one_mc
                        else:
                            run_fn_base = _GRADIENT_OPTIMIZERS[optimizer_name]
                            grad_kwargs = _hp_kwargs_for_gradient(optimizer_name, hp_cfg)
                            run_kwargs = dict(
                                diff_fn=diff_fn,
                                voxel=voxel,
                                vertices=vertices,
                                R_init=R_start,
                                scale=FIXED_SCALE,
                                **grad_kwargs,
                            )
                            run_fn = run_fn_base

                        result, wall_time, peak_mem, rss_delta = measure_run(
                            run_fn, **run_kwargs
                        )

                        # Build and write main summary row
                        main_row = _result_to_main_row(
                            run_id=run_id,
                            optimizer=optimizer_name,
                            hp_cfg=hp_cfg,
                            R_start=R_start,
                            R_gt=R_gt,
                            result=result,
                            wall_time=wall_time,
                            peak_mem=peak_mem,
                            rss_delta=rss_delta,
                            voxel_idx=vidx,
                            perturbation_deg=pert_deg,
                        )
                        main_writer.writerow(main_row)
                        main_f.flush()

                        # Write trajectory rows
                        traj = result.get("trajectory", [])
                        for rec in traj:
                            traj_row = {
                                "run_id": run_id,
                                "event_idx": rec["event_idx"],
                                "step": rec["step"],
                                "event_type": rec["event_type"],
                                "angular_step_deg": rec["angular_step_deg"],
                                "misori_gt_deg": rec.get("misori_gt_deg", float("nan")),
                                "quality": rec.get("quality", float("nan")),
                                "cur_step_rad": rec.get("cur_step_rad", float("nan")),
                            }
                            traj_writer.writerow(traj_row)
                        if traj:
                            traj_f.flush()

                        print(
                            f"  {optimizer_name:<28s} hp={hp_id:3d}"
                            f"  vox={vidx:5d}  pert={pert_deg:.0f}°"
                            f"  misori={main_row['final_misorientation_deg']:7.3f}°"
                            f"  q={main_row['final_quality']:.4f}"
                            f"  evals={main_row['n_evaluations']:5d}"
                            f"  {wall_time:.1f}s",
                            flush=True,
                        )

                        run_id += 1
                        total_runs += 1

                    if smoke_test and total_runs >= 3:
                        break
                if smoke_test and total_runs >= 3:
                    break
            if smoke_test and total_runs >= 3:
                break

    finally:
        main_f.close()
        traj_f.close()

    elapsed = time.perf_counter() - t_start
    print(f"\nTotal: {total_runs} runs in {elapsed/60:.1f} min")
    print(f"Saved main CSV: {main_csv_path}")
    print(f"Saved trajectory CSV: {traj_csv_path}")

    return run_id


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------

OPTIMIZER_COLORS = {
    "riemannian_adam_geoopt":   "C0",
    "riemannian_adam_manual":   "C1",
    "riemannian_sgd_plain":     "C2",
    "riemannian_sgd_momentum":  "C3",
    "riemannian_sgld":          "C4",
    "mc_optimizer":             "C5",
}


def plot_lr_sensitivity(csv_path: Path, out_path: Path, title_prefix: str) -> None:
    """Box plots of final misorientation distribution vs lr, one subplot per gradient optimizer."""
    try:
        import pandas as pd
    except ImportError:
        print("pandas not available, skipping lr sensitivity plot")
        return

    df = pd.read_csv(csv_path)
    grad_opts = [o for o in OPTIMIZER_COLORS if o != "mc_optimizer" and o in df["optimizer"].unique()]

    if not grad_opts:
        return

    ncols = 3
    nrows = math.ceil(len(grad_opts) / ncols)
    fig, axes = plt.subplots(nrows, ncols, figsize=(5 * ncols, 4 * nrows), sharey=True)
    axes_flat = axes.flat if hasattr(axes, "flat") else [axes]

    for ax, opt in zip(axes_flat, grad_opts):
        sub = df[df["optimizer"] == opt].dropna(subset=["lr"])
        if sub.empty:
            ax.set_visible(False)
            continue
        lr_vals = sorted(sub["lr"].unique())
        data = [sub[sub["lr"] == lr]["final_misorientation_deg"].values for lr in lr_vals]
        positions = list(range(len(lr_vals)))
        color = OPTIMIZER_COLORS.get(opt, "gray")
        bp = ax.boxplot(
            data,
            positions=positions,
            widths=0.5,
            patch_artist=True,
            medianprops={"color": "k", "lw": 1.5},
            flierprops={"marker": ".", "markersize": 3, "alpha": 0.4, "markeredgecolor": color},
            whiskerprops={"color": color},
            capprops={"color": color},
        )
        for patch in bp["boxes"]:
            patch.set_facecolor(color)
            patch.set_alpha(0.5)
        ax.set_xticks(positions)
        ax.set_xticklabels([f"{lr:.0e}" for lr in lr_vals], rotation=45, fontsize=7)
        ax.set_xlabel("Learning rate")
        ax.set_ylabel("Final misorientation (°)")
        ax.set_title(opt, fontsize=9)
        ax.grid(True, axis="y", alpha=0.3)

    # Hide unused subplots
    for ax in list(axes_flat)[len(grad_opts):]:
        ax.set_visible(False)

    fig.suptitle(f"{title_prefix} — LR sensitivity (misorientation distribution)", fontsize=11)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    print(f"Saved: {out_path}")


def plot_nsteps_sensitivity(csv_path: Path, out_path: Path, title_prefix: str) -> None:
    """Box plots of final misorientation distribution vs n_steps, one subplot per optimizer."""
    try:
        import pandas as pd
    except ImportError:
        return

    df = pd.read_csv(csv_path)
    opts = [o for o in OPTIMIZER_COLORS if o in df["optimizer"].unique()]
    if not opts:
        return

    ncols = 3
    nrows = math.ceil(len(opts) / ncols)
    fig, axes = plt.subplots(nrows, ncols, figsize=(5 * ncols, 4 * nrows), sharey=True)
    axes_flat = axes.flat if hasattr(axes, "flat") else [axes]

    for ax, opt in zip(axes_flat, opts):
        sub = df[df["optimizer"] == opt].dropna(subset=["n_steps"])
        if sub.empty:
            ax.set_visible(False)
            continue
        step_vals = sorted(sub["n_steps"].unique())
        data = [sub[sub["n_steps"] == s]["final_misorientation_deg"].values for s in step_vals]
        positions = list(range(len(step_vals)))
        color = OPTIMIZER_COLORS.get(opt, "gray")
        bp = ax.boxplot(
            data,
            positions=positions,
            widths=0.5,
            patch_artist=True,
            medianprops={"color": "k", "lw": 1.5},
            flierprops={"marker": ".", "markersize": 3, "alpha": 0.4, "markeredgecolor": color},
            whiskerprops={"color": color},
            capprops={"color": color},
        )
        for patch in bp["boxes"]:
            patch.set_facecolor(color)
            patch.set_alpha(0.5)
        ax.set_xticks(positions)
        ax.set_xticklabels([str(int(s)) for s in step_vals], fontsize=8)
        ax.set_xlabel("n_steps / max_mc_steps")
        ax.set_ylabel("Final misorientation (°)")
        ax.set_title(opt, fontsize=9)
        ax.grid(True, axis="y", alpha=0.3)

    # Hide unused subplots
    for ax in list(axes_flat)[len(opts):]:
        ax.set_visible(False)

    fig.suptitle(f"{title_prefix} — Steps sensitivity (misorientation distribution)", fontsize=11)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    print(f"Saved: {out_path}")


def plot_optimizer_comparison(csv_path: Path, out_path: Path, title_prefix: str) -> None:
    """Bar chart: best-HP misorientation per optimizer × perturbation."""
    try:
        import pandas as pd
    except ImportError:
        return

    df = pd.read_csv(csv_path)
    perts = sorted(df["perturbation_deg"].unique())
    opts = [o for o in OPTIMIZER_COLORS if o in df["optimizer"].unique()]
    if not opts or not perts:
        return

    fig, axes = plt.subplots(1, len(perts), figsize=(5 * len(perts), 5), sharey=True)
    if len(perts) == 1:
        axes = [axes]

    for ai, pert in enumerate(perts):
        ax = axes[ai]
        sub = df[df["perturbation_deg"] == pert]
        # For each optimizer, take the best-HP (minimum mean misorientation)
        best_per_opt = []
        for opt in opts:
            sub_o = sub[sub["optimizer"] == opt]
            if sub_o.empty:
                continue
            grp = sub_o.groupby("hp_id")["final_misorientation_deg"].mean()
            best_per_opt.append((opt, grp.min()))

        labels = [o for o, _ in best_per_opt]
        values = [v for _, v in best_per_opt]
        colors = [OPTIMIZER_COLORS.get(o, "gray") for o in labels]
        ax.bar(labels, values, color=colors, edgecolor="k", alpha=0.8)
        ax.set_title(f"{title_prefix}\nPerturbation = {pert:.0f}°")
        ax.set_ylabel("Best-HP mean misori (°)" if ai == 0 else "")
        ax.tick_params(axis="x", rotation=45)
        ax.grid(True, axis="y", alpha=0.3)

    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    print(f"Saved: {out_path}")


def plot_trajectory_step_sizes(traj_csv_path: Path, main_csv_path: Path,
                                out_path: Path, title_prefix: str) -> None:
    """Angular step size over event index, one line per optimizer (median across runs)."""
    try:
        import pandas as pd
    except ImportError:
        return

    traj = pd.read_csv(traj_csv_path)
    main = pd.read_csv(main_csv_path)[["run_id", "optimizer"]].drop_duplicates()
    traj = traj.merge(main, on="run_id", how="left")

    opts = [o for o in OPTIMIZER_COLORS if o in traj["optimizer"].unique()]
    if not opts:
        return

    fig, ax = plt.subplots(figsize=(10, 5))
    for opt in opts:
        sub = traj[traj["optimizer"] == opt]
        if sub.empty:
            continue
        grp = sub.groupby("event_idx")["angular_step_deg"]
        median_v = grp.median()
        p25 = grp.quantile(0.25)
        p75 = grp.quantile(0.75)
        color = OPTIMIZER_COLORS.get(opt, "gray")
        ax.plot(median_v.index, median_v.values, color=color, label=opt, lw=1.5)
        ax.fill_between(median_v.index, p25, p75, color=color, alpha=0.12)

    ax.set_xlabel("Event index (gradient step / MC accepted move)")
    ax.set_ylabel("Angular step size (°)")
    ax.set_title(f"{title_prefix} — Optimization trajectory step sizes")
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)
    ax.set_ylim(bottom=0)
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
    n_voxels_default: Optional[int],
    n_voxels_override: Optional[int],
    optimizers_to_run: List[str],
    perturbations_deg: List[float],
    smoke_test: bool,
) -> None:
    print(f"\n{'='*60}")
    print(f"Example: {label}")
    print(f"geoopt available: {_GEOOPT_AVAILABLE}")
    print(f"psutil available: {_PSUTIL_AVAILABLE}")
    print(f"Optimizers: {optimizers_to_run}")
    print(f"{'='*60}")

    if smoke_test:
        perturbations_deg = [2.0]
        n_voxels_override = 1
        print("  [SMOKE TEST: 1 voxel, 1 HP config, reduced steps]")

    mic, hard_fn, diff_fn, get_vertices = setup_example(example_dir, basename)

    n_voxels = n_voxels_override if n_voxels_override is not None else n_voxels_default
    rng = np.random.default_rng(SEED)

    if n_voxels is None:
        # Use all qualifying voxels (ThreeVoxels path)
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
                print(f"    idx={idx}  hard_q={info.quality:.4f}  "
                      f"φ1={euler[0]:.2f}°  Φ={euler[1]:.2f}°  φ2={euler[2]:.2f}°  [selected]")
        print(f"  Using all {len(voxel_list)} qualifying voxels.")
    else:
        voxel_list = select_voxels(mic, hard_fn, get_vertices, n_voxels, rng)

    hp_grids = build_hp_grids()

    if smoke_test:
        # Truncate HP grids to 1 config each; override n_steps to 5
        for opt_name in hp_grids:
            cfg = hp_grids[opt_name][0]
            if "n_steps" in cfg:
                cfg["n_steps"] = 5
            if "max_mc_steps" in cfg:
                cfg["max_mc_steps"] = 10
            hp_grids[opt_name] = [cfg]

    tag = label.lower().replace(" ", "_").replace(".", "")
    main_csv = benchmark_dir / f"hp_sweep_{tag}.csv"
    traj_csv = benchmark_dir / f"hp_sweep_trajectory_{tag}.csv"

    run_hp_sweep(
        label=label,
        example_dir=example_dir,
        basename=basename,
        voxel_list=voxel_list,
        hard_fn=hard_fn,
        diff_fn=diff_fn,
        get_vertices=get_vertices,
        hp_grids=hp_grids,
        optimizers_to_run=optimizers_to_run,
        perturbations_deg=perturbations_deg,
        rng=rng,
        main_csv_path=main_csv,
        traj_csv_path=traj_csv,
        smoke_test=smoke_test,
    )

    # Plots
    plot_lr_sensitivity(main_csv, benchmark_dir / f"hp_sweep_lr_sensitivity_{tag}.png", label)
    plot_nsteps_sensitivity(
        main_csv, benchmark_dir / f"hp_sweep_nsteps_sensitivity_{tag}.png", label
    )
    plot_optimizer_comparison(
        main_csv, benchmark_dir / f"hp_sweep_optimizer_comparison_{tag}.png", label
    )
    plot_trajectory_step_sizes(
        traj_csv, main_csv, benchmark_dir / f"hp_sweep_trajectory_{tag}.png", label
    )


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

ALL_GRADIENT_OPTS = [
    "riemannian_adam_geoopt",
    "riemannian_adam_manual",
    "riemannian_sgd_plain",
    "riemannian_sgd_momentum",
    "riemannian_sgld",
]
ALL_OPTS = ALL_GRADIENT_OPTS + ["mc_optimizer"]


def main() -> None:
    base_dir = project_root / "Examples"

    parser = argparse.ArgumentParser(description="HP sweep: gradient methods vs MC optimizer")
    parser.add_argument(
        "--example",
        choices=["threevoxels", "manygrains", "both"],
        default="threevoxels",
    )
    parser.add_argument(
        "--optimizer",
        choices=["all", "gradient", "mc"] + ALL_OPTS,
        default="all",
        help="Which optimizers to run. 'gradient' runs all gradient methods; 'mc' runs MC only.",
    )
    parser.add_argument("--smoke-test", action="store_true", help="Quick smoke test (1 voxel, reduced steps)")
    parser.add_argument("--n-voxels", type=int, default=None, help="Override number of voxels")
    parser.add_argument(
        "--plots-only",
        action="store_true",
        help="Skip running; just regenerate plots from existing CSVs",
    )
    args = parser.parse_args()

    if args.optimizer == "all":
        to_run = ALL_OPTS
    elif args.optimizer == "gradient":
        to_run = ALL_GRADIENT_OPTS
    elif args.optimizer == "mc":
        to_run = ["mc_optimizer"]
    else:
        to_run = [args.optimizer]

    # Filter out geoopt-dependent optimizers if not available
    if not _GEOOPT_AVAILABLE:
        geoopt_opts = {"riemannian_adam_geoopt", "riemannian_sgd_plain", "riemannian_sgd_momentum"}
        excluded = [o for o in to_run if o in geoopt_opts]
        if excluded:
            print(f"Warning: geoopt not available. Skipping: {excluded}")
        to_run = [o for o in to_run if o not in geoopt_opts]

    examples = []
    if args.example in ("threevoxels", "both"):
        examples.append((
            "ThreeVoxels",
            base_dir / "Example2.ThreeVoxels",
            "3Grains.sim",
            None,    # None → use all qualifying voxels (3 for ThreeVoxels)
        ))
    if args.example in ("manygrains", "both"):
        examples.append((
            "ManyGrains",
            base_dir / "Example2.ManyGrains",
            "500Grains.sim",
            N_VOXELS_MANY,
        ))

    for label, ex_dir, basename, n_vox_default in examples:
        if args.plots_only:
            tag = label.lower().replace(" ", "_").replace(".", "")
            benchmark_dir = Path(__file__).parent
            main_csv = benchmark_dir / f"hp_sweep_{tag}.csv"
            traj_csv = benchmark_dir / f"hp_sweep_trajectory_{tag}.csv"
            if not main_csv.exists():
                print(f"Skipping plots for {label}: {main_csv} not found")
                continue
            print(f"Regenerating plots for {label} from {main_csv}")
            plot_lr_sensitivity(main_csv, benchmark_dir / f"hp_sweep_lr_sensitivity_{tag}.png", label)
            plot_nsteps_sensitivity(
                main_csv, benchmark_dir / f"hp_sweep_nsteps_sensitivity_{tag}.png", label
            )
            plot_optimizer_comparison(
                main_csv, benchmark_dir / f"hp_sweep_optimizer_comparison_{tag}.png", label
            )
            if traj_csv.exists():
                plot_trajectory_step_sizes(
                    traj_csv, main_csv, benchmark_dir / f"hp_sweep_trajectory_{tag}.png", label
                )
        else:
            run_example(
                label=label,
                example_dir=ex_dir,
                basename=basename,
                n_voxels_default=n_vox_default,
                n_voxels_override=args.n_voxels,
                optimizers_to_run=to_run,
                perturbations_deg=PERTURBATIONS_DEG,
                smoke_test=args.smoke_test,
            )


if __name__ == "__main__":
    main()
