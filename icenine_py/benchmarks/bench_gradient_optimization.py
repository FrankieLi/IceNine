#!/usr/bin/env python3
"""
Gradient-based orientation optimization benchmark.

Tests whether Adam gradient descent can recover a ground-truth orientation
from a perturbed starting point, sweeping all 9 combinations of:
  scale         ∈ {0, 1, 2}  (1×, 4×, 8× downsampled)
  omega_window  ∈ {0, 1, 2}  (no blend, ±1 frame, ±2 frames)

at 3 starting perturbation distances: 1°, 2°, 5°.

Algorithm:
  theta = so3_log(R_perturbed)       # axis-angle init (3-vector)
  theta.requires_grad_(True)
  optimizer = Adam([theta], lr=0.01)
  for step in range(100):
      optimizer.zero_grad()
      R = torch.matrix_exp(skew(theta))
      info = diff_fn.evaluate(R, vertices, phase_index, scale=scale)
      info.cost.backward()
      optimizer.step()

Parameterization note:
  theta lives in the Lie algebra ℝ³ (unconstrained). matrix_exp(skew(theta))
  always returns a valid rotation matrix regardless of where Adam takes theta.
  Gradients flow through matrix_exp and skew without approximation.

Outputs (all in icenine_py/benchmarks/):
  grad_opt_threevoxels.csv     — per-step: voxel_idx, perturbation_deg, scale,
                                   omega_window, step, quality, misorientation_deg
  grad_opt_manygrains.csv      — same, 20 voxels
  grad_opt_convergence_*.png   — quality vs. step, one panel per perturbation
  grad_opt_summary_*.png       — 3×3 heatmap of final misorientation

Usage:
  cd /Users/sfli/Research/IceNine/icenine_py
  uv run python benchmarks/bench_gradient_optimization.py --smoke-test
  uv run python benchmarks/bench_gradient_optimization.py
"""

import argparse
import math
import os
import sys
import time
from pathlib import Path
from typing import List, Dict, Tuple, Optional

import numpy as np
import torch
import torch.nn.functional as F
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

project_root = Path(__file__).parent.parent.parent
sys.path.insert(0, str(project_root / "icenine_py"))

benchmark_dir = Path(__file__).parent

# ---------------------------------------------------------------------------
# Sweep parameters
# ---------------------------------------------------------------------------

# Scale 0 (full 2048² resolution) densifies images on-the-fly during evaluate(),
# which is ~2-3 min per 100-step run (vs ~5s for scale 1 or 2 with pre-densified
# stacks). Scale 0 also has no basin-widening — gradient is zero at >0.3° — so
# it cannot converge from the tested perturbations. Skip it in the sweep.
SCALES = [1, 2]
OMEGA_WINDOWS = [0, 1, 2]
PERTURBATIONS_DEG = [1.0, 2.0, 5.0]
N_STEPS = 100
LR = 0.01

SEED = 42
N_VOXELS_MANY = 20
QUALITY_THRESHOLD = 0.1
MAX_SCAN = 500

# Colors / linestyles for (scale, omega_window) combos
# scale: 0→C0, 1→C1, 2→C3
# omega_window: 0→':', 1→'-.', 2→'--'
SCALE_COLORS   = {0: "C0", 1: "C1", 2: "C3"}
OW_LINESTYLES  = {0: ":",   1: "-.", 2: "--"}


# ---------------------------------------------------------------------------
# SO(3) helpers
# ---------------------------------------------------------------------------

def skew(theta: torch.Tensor) -> torch.Tensor:
    """3-vector → 3×3 skew-symmetric matrix."""
    assert theta.shape == (3,), f"Expected shape (3,), got {theta.shape}"
    z = torch.zeros(1, dtype=theta.dtype, device=theta.device)
    row0 = torch.stack([z.squeeze(), -theta[2],  theta[1]])
    row1 = torch.stack([theta[2],    z.squeeze(), -theta[0]])
    row2 = torch.stack([-theta[1],   theta[0],   z.squeeze()])
    return torch.stack([row0, row1, row2])


def so3_log(R: np.ndarray) -> torch.Tensor:
    """Rotation matrix → axis-angle 3-vector (for initialization only).

    Uses the formula: theta = angle * axis, where
      angle = arccos((trace(R) - 1) / 2)
      axis  = (R - R^T) / (2 sin(angle))
    Returns zero vector for near-identity rotations.
    """
    R_t = torch.from_numpy(R).float()
    trace = R_t.trace()
    cos_angle = ((trace - 1.0) / 2.0).clamp(-1.0, 1.0)
    angle = torch.acos(cos_angle)

    if angle.abs() < 1e-7:
        return torch.zeros(3)

    # Skew-symmetric part
    W = (R_t - R_t.T) / (2.0 * torch.sin(angle))
    # Extract axis from skew-symmetric matrix
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
    K = np.array([
        [0,        -axis[2],  axis[1]],
        [axis[2],   0,       -axis[0]],
        [-axis[1],  axis[0],  0      ],
    ], dtype=np.float64)
    return np.eye(3) + math.sin(angle_rad) * K + (1 - math.cos(angle_rad)) * (K @ K)


def random_unit_axis(rng: np.random.Generator) -> np.ndarray:
    v = rng.standard_normal(3)
    return v / np.linalg.norm(v)


# ---------------------------------------------------------------------------
# Setup
# ---------------------------------------------------------------------------

def setup_example(example_dir: Path, basename: str, omega_windows: List[int]):
    """Load cost functions and mic for one example.

    Returns:
      mic, hard_fn, diff_fns, get_vertices
      diff_fns: dict mapping omega_window → DifferentiableCostFunction
                (each contains scales 0, 1, 2 internally)
    """
    import os as _os
    import psutil as _psutil

    def _rss():
        return _psutil.Process(_os.getpid()).memory_info().rss / 1024 / 1024

    def _mem(label: str, prev: float) -> float:
        r = _rss()
        print(f"  [MEM] {label:<50s} {r:6.0f} MB  Δ={r-prev:+.0f} MB", flush=True)
        return r

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

    os.chdir(example_dir)  # configs use relative paths

    _t = _rss()
    print(f"  [MEM] start setup  {_t:.0f} MB", flush=True)

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
    _t = _mem("physics setup done", _t)

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
    _t = _mem("SparseImageStack loaded", _t)

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
    _t = _mem("ExperimentalData loaded", _t)

    hard_fn = VoxelCostFunction(
        simulator=simulator,
        detector_list=detector_list,
        range_map=range_map,
        exp_data=exp_data,
        sample=sample,
        structure_list=structure_list,
        mode="hard",
    )
    _t = _mem("VoxelCostFunction built", _t)

    # Build downsampled stacks ONCE and share across all omega_window variants.
    # Without sharing, each MultiScaleImageStack re-densifies all 360 frames,
    # costing ~360 MB per call × n_omega_windows → ~1 GB unnecessary duplication.
    print("  Building shared downsampled base stacks [4x, 8x] ...")
    shared_ds = MultiScaleImageStack.build_shared_base(
        image_stack, downsample_factors=[1, 4, 8]
    )
    sz_mb = sum(s.images.numel() * 4 / 1024 / 1024 for s in shared_ds)
    print(f"  Shared base: {sz_mb:.0f} MB total")
    _t = _mem("build_shared_base done", _t)

    diff_fns: Dict[int, DifferentiableCostFunction] = {}
    for ow in omega_windows:
        print(f"  Building MultiScaleImageStack omega_window={ow} (reusing base) ...")
        ms = MultiScaleImageStack(
            image_stack,
            downsample_factors=[1, 4, 8],
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
        _t = _mem(f"DifferentiableCostFunction ow={ow}", _t)

    return mic, hard_fn, diff_fns, _get_voxel_vertices


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
        raise RuntimeError(
            f"Only {len(candidates)} candidates found (need {n}) — "
            f"lower threshold or increase max_scan"
        )

    chosen_positions = rng.choice(len(candidates), size=n, replace=False)
    chosen = [candidates[i] for i in sorted(chosen_positions)]

    print(f"  Selected {n} voxels from {len(candidates)} candidates:")
    for voxel, vidx, q in chosen:
        euler = matrix_to_euler(voxel.orientation)
        print(f"    idx={vidx:5d}  φ1={euler[0]:7.2f}°  Φ={euler[1]:6.2f}°  "
              f"φ2={euler[2]:7.2f}°  hard_q={q:.4f}")

    return [(v, vi) for v, vi, _ in chosen]


# ---------------------------------------------------------------------------
# Core optimization loop
# ---------------------------------------------------------------------------

def run_one(diff_fn, voxel, vertices, R_init: np.ndarray, scale: int,
            n_steps: int = N_STEPS, lr: float = LR) -> Dict:
    """Run Adam optimizer for n_steps starting from R_init.

    Returns:
      dict with keys:
        quality_history: (n_steps+1,) array — quality at step 0 (init), 1, ..., n_steps
        misorientation_history: (n_steps+1,) array — degrees from ground truth
        R_final: (3,3) ndarray — final rotation matrix
        n_peaks: int
    """
    R_gt = voxel.orientation  # ground truth

    theta = so3_log(R_init)
    theta = theta.detach().requires_grad_(True)
    optimizer = torch.optim.Adam([theta], lr=lr)

    quality_hist = []
    misori_hist = []

    for step in range(n_steps + 1):
        if step > 0:
            optimizer.zero_grad()

        with torch.set_grad_enabled(step > 0):
            R_t = torch.matrix_exp(skew(theta))
            info = diff_fn.evaluate(R_t, vertices, phase_index=voxel.phase, scale=scale)

        R_np = R_t.detach().numpy()
        quality_hist.append(float(info.quality.item()))
        misori_hist.append(misorientation_deg(R_gt, R_np))

        if step > 0:
            info.cost.backward()
            optimizer.step()

    R_final = torch.matrix_exp(skew(theta.detach())).numpy()
    return {
        "quality_history": np.array(quality_hist),
        "misorientation_history": np.array(misori_hist),
        "R_final": R_final,
        "n_peaks": info.n_peaks,
    }


# ---------------------------------------------------------------------------
# Full sweep
# ---------------------------------------------------------------------------

def sweep_voxel(voxel, vidx: int, diff_fns: Dict[int, object],
                get_vertices, perturbations_deg: List[float],
                rng: np.random.Generator, scales: List[int],
                omega_windows: List[int], n_steps: int) -> List[Dict]:
    """Run optimizer for all (scale, omega_window, perturbation) combos for one voxel.

    Returns a flat list of result dicts, each with metadata keys added.
    """
    vertices = get_vertices(voxel)
    rows = []

    # One random axis per perturbation (seeded via rng, consistent across configs)
    axes = [random_unit_axis(rng) for _ in perturbations_deg]

    for pi, (pert_deg, axis) in enumerate(zip(perturbations_deg, axes)):
        pert_rad = pert_deg * math.pi / 180.0
        R_pert = rodrigues_np(axis, pert_rad) @ voxel.orientation

        for scale in scales:
            for ow in omega_windows:
                t0 = time.perf_counter()
                result = run_one(
                    diff_fns[ow], voxel, vertices, R_pert,
                    scale=scale, n_steps=n_steps,
                )
                elapsed = time.perf_counter() - t0

                final_misori = result["misorientation_history"][-1]
                final_quality = result["quality_history"][-1]
                print(
                    f"    vox {vidx:5d}  pert={pert_deg:.0f}°  s={scale}  ω±{ow}"
                    f"  final_misori={final_misori:.3f}°  q={final_quality:.4f}"
                    f"  ({elapsed:.0f}s)",
                    flush=True,
                )

                for step in range(n_steps + 1):
                    rows.append({
                        "voxel_idx": vidx,
                        "perturbation_deg": pert_deg,
                        "scale": scale,
                        "omega_window": ow,
                        "step": step,
                        "quality": result["quality_history"][step],
                        "misorientation_deg": result["misorientation_history"][step],
                    })

    return rows


# ---------------------------------------------------------------------------
# CSV save
# ---------------------------------------------------------------------------

def save_csv(rows: List[Dict], out_path: Path):
    import csv
    fieldnames = ["voxel_idx", "perturbation_deg", "scale", "omega_window",
                  "step", "quality", "misorientation_deg"]
    with open(out_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    print(f"Saved CSV: {out_path}  ({len(rows)} rows)")


# ---------------------------------------------------------------------------
# Plots
# ---------------------------------------------------------------------------

def label_for(scale: int, ow: int) -> str:
    return f"s{scale} ω±{ow}"


def convergence_plot(csv_path: Path, out_path: Path, title_prefix: str,
                     scales: List[int], omega_windows: List[int],
                     perturbations_deg: List[float], n_steps: int):
    """Quality vs. step — one panel per perturbation, lines per (scale, ow) combo."""
    import pandas as pd

    df = pd.read_csv(csv_path)
    n_pert = len(perturbations_deg)
    fig, axes = plt.subplots(1, n_pert, figsize=(6 * n_pert, 5), sharey=True)
    if n_pert == 1:
        axes = [axes]

    for ai, pert in enumerate(perturbations_deg):
        ax = axes[ai]
        sub = df[df["perturbation_deg"] == pert]

        for scale in scales:
            for ow in omega_windows:
                sel = sub[(sub["scale"] == scale) & (sub["omega_window"] == ow)]
                if sel.empty:
                    continue
                # Mean over voxels at each step
                mean_q = sel.groupby("step")["quality"].mean()
                std_q  = sel.groupby("step")["quality"].std().fillna(0)
                steps = mean_q.index.values
                color = SCALE_COLORS[scale]
                ls    = OW_LINESTYLES[ow]
                ax.plot(steps, mean_q.values, color=color, ls=ls, lw=1.8,
                        label=label_for(scale, ow))
                ax.fill_between(steps, mean_q - std_q, mean_q + std_q,
                                color=color, alpha=0.12)

        ax.set_xlabel("Adam step")
        ax.set_ylabel("Quality (0–1)" if ai == 0 else "")
        ax.set_title(f"{title_prefix}\nPerturbation = {pert:.0f}°")
        ax.legend(fontsize=7, ncol=3)
        ax.grid(True, alpha=0.3)
        ax.set_xlim(0, n_steps)
        ax.set_ylim(-0.02, 1.05)

    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    print(f"Saved: {out_path}")


def misorientation_convergence_plot(csv_path: Path, out_path: Path, title_prefix: str,
                                    scales: List[int], omega_windows: List[int],
                                    perturbations_deg: List[float], n_steps: int):
    """Misorientation vs. step — one panel per perturbation."""
    import pandas as pd

    df = pd.read_csv(csv_path)
    n_pert = len(perturbations_deg)
    fig, axes = plt.subplots(1, n_pert, figsize=(6 * n_pert, 5), sharey=False)
    if n_pert == 1:
        axes = [axes]

    for ai, pert in enumerate(perturbations_deg):
        ax = axes[ai]
        sub = df[df["perturbation_deg"] == pert]

        for scale in scales:
            for ow in omega_windows:
                sel = sub[(sub["scale"] == scale) & (sub["omega_window"] == ow)]
                if sel.empty:
                    continue
                mean_m = sel.groupby("step")["misorientation_deg"].mean()
                std_m  = sel.groupby("step")["misorientation_deg"].std().fillna(0)
                steps = mean_m.index.values
                color = SCALE_COLORS[scale]
                ls    = OW_LINESTYLES[ow]
                ax.plot(steps, mean_m.values, color=color, ls=ls, lw=1.8,
                        label=label_for(scale, ow))
                ax.fill_between(steps, mean_m - std_m, mean_m + std_m,
                                color=color, alpha=0.12)

        ax.axhline(0, color="gray", lw=0.8, ls="--")
        ax.set_xlabel("Adam step")
        ax.set_ylabel("Misorientation (°)" if ai == 0 else "")
        ax.set_title(f"{title_prefix}\nPerturbation = {pert:.0f}°")
        ax.legend(fontsize=7, ncol=3)
        ax.grid(True, alpha=0.3)
        ax.set_xlim(0, n_steps)
        ax.set_ylim(bottom=0)

    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    print(f"Saved: {out_path}")


def summary_heatmap(csv_path: Path, out_path: Path, title_prefix: str,
                    scales: List[int], omega_windows: List[int],
                    perturbations_deg: List[float], n_steps: int):
    """3-panel heatmap: rows=scale, cols=omega_window, cell=final misorientation (°)."""
    import pandas as pd

    df = pd.read_csv(csv_path)
    df_final = df[df["step"] == n_steps]

    n_pert = len(perturbations_deg)
    fig, axes = plt.subplots(1, n_pert, figsize=(4 * n_pert, 4))
    if n_pert == 1:
        axes = [axes]

    for ai, pert in enumerate(perturbations_deg):
        ax = axes[ai]
        sub = df_final[df_final["perturbation_deg"] == pert]

        grid = np.full((len(scales), len(omega_windows)), np.nan)
        for si, scale in enumerate(scales):
            for oi, ow in enumerate(omega_windows):
                sel = sub[(sub["scale"] == scale) & (sub["omega_window"] == ow)]
                if not sel.empty:
                    grid[si, oi] = sel["misorientation_deg"].mean()

        # Color: green < 0.5°, yellow 0.5–2°, red > 2°
        im = ax.imshow(grid, vmin=0, vmax=5, cmap="RdYlGn_r", aspect="auto")
        plt.colorbar(im, ax=ax, label="Final misori (°)")

        ax.set_xticks(range(len(omega_windows)))
        ax.set_xticklabels([f"ω±{ow}" for ow in omega_windows])
        ax.set_yticks(range(len(scales)))
        ax.set_yticklabels([f"s{s}" for s in scales])
        ax.set_xlabel("omega_window")
        ax.set_ylabel("scale" if ai == 0 else "")
        ax.set_title(f"{title_prefix}\nPert={pert:.0f}°  (mean over voxels)")

        # Annotate cells with numeric value
        for si in range(len(scales)):
            for oi in range(len(omega_windows)):
                val = grid[si, oi]
                if not np.isnan(val):
                    ax.text(oi, si, f"{val:.2f}°", ha="center", va="center",
                            fontsize=9, color="black")

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
    n_voxels: Optional[int],   # None = use all voxels in mic
    scales: List[int],
    omega_windows: List[int],
    perturbations_deg: List[float],
    n_steps: int,
    rng: np.random.Generator,
    smoke_test: bool = False,
):
    print(f"\n{'='*60}")
    print(f"Example: {label}")
    print(f"{'='*60}")

    # Smoke-test overrides
    if smoke_test:
        n_voxels   = min(n_voxels or 1, 1)
        n_steps    = 5
        scales     = [2]
        omega_windows = [0, 1]
        perturbations_deg = [2.0]
        print("  [SMOKE TEST MODE: 1 voxel, 5 steps, scale=2 only]")

    mic, hard_fn, diff_fns, get_vertices = setup_example(
        example_dir, basename, omega_windows
    )

    if n_voxels is None:
        # Use ALL voxels (ThreeVoxels: just 3)
        all_voxels = [(v, i) for i, v in enumerate(mic.voxels)]
        # Filter by quality
        print(f"  Evaluating hard quality for {len(all_voxels)} voxels ...")
        voxel_list = []
        for voxel, vidx in all_voxels:
            vertices = get_vertices(voxel)
            info = hard_fn.evaluate(voxel.orientation, vertices, voxel.phase)
            if info.quality > QUALITY_THRESHOLD:
                voxel_list.append((voxel, vidx))
                print(f"    idx={vidx}  hard_q={info.quality:.4f}  [selected]")
            else:
                print(f"    idx={vidx}  hard_q={info.quality:.4f}  [skipped]")
        print(f"  Using {len(voxel_list)} voxels.")
    else:
        voxel_list = select_voxels(
            mic, hard_fn, get_vertices, n_voxels, rng
        )

    total = (len(voxel_list) * len(perturbations_deg)
             * len(scales) * len(omega_windows))
    print(f"\nSweeping {len(voxel_list)} voxels × {len(perturbations_deg)} perts"
          f" × {len(scales)} scales × {len(omega_windows)} ω_windows = {total} runs"
          f" of {n_steps} steps each ...")

    all_rows = []
    t_start = time.perf_counter()

    for vi, (voxel, vidx) in enumerate(voxel_list):
        print(f"\n  Voxel {vi+1}/{len(voxel_list)}  (mic idx={vidx}) ...")
        rows = sweep_voxel(
            voxel, vidx, diff_fns, get_vertices,
            perturbations_deg=perturbations_deg,
            rng=rng,
            scales=scales,
            omega_windows=omega_windows,
            n_steps=n_steps,
        )
        all_rows.extend(rows)

    elapsed = time.perf_counter() - t_start
    print(f"\nTotal sweep time: {elapsed/60:.1f} min  ({len(all_rows)} rows)")

    tag = label.lower().replace(" ", "_").replace(".", "")
    csv_path = benchmark_dir / f"grad_opt_{tag}.csv"
    save_csv(all_rows, csv_path)

    convergence_plot(
        csv_path,
        benchmark_dir / f"grad_opt_quality_convergence_{tag}.png",
        title_prefix=label,
        scales=scales, omega_windows=omega_windows,
        perturbations_deg=perturbations_deg, n_steps=n_steps,
    )
    misorientation_convergence_plot(
        csv_path,
        benchmark_dir / f"grad_opt_misori_convergence_{tag}.png",
        title_prefix=label,
        scales=scales, omega_windows=omega_windows,
        perturbations_deg=perturbations_deg, n_steps=n_steps,
    )
    summary_heatmap(
        csv_path,
        benchmark_dir / f"grad_opt_summary_{tag}.png",
        title_prefix=label,
        scales=scales, omega_windows=omega_windows,
        perturbations_deg=perturbations_deg, n_steps=n_steps,
    )


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Gradient optimization benchmark for differentiable cost function"
    )
    parser.add_argument(
        "--smoke-test",
        action="store_true",
        help="Quick sanity check: 1 voxel, 5 steps, scale=2, omega_window=[0,1], pert=2°",
    )
    parser.add_argument(
        "--example",
        choices=["threevoxels", "manygrains", "both"],
        default="both",
        help="Which example to run (default: both)",
    )
    args = parser.parse_args()

    rng = np.random.default_rng(seed=SEED)

    three_dir  = project_root / "Examples" / "Example2.ThreeVoxels"
    many_dir   = project_root / "Examples" / "Example2.ManyGrains"

    if args.example in ("threevoxels", "both"):
        run_example(
            label="ThreeVoxels",
            example_dir=three_dir,
            basename="3Grains.sim",
            n_voxels=None,              # use all voxels
            scales=SCALES,
            omega_windows=OMEGA_WINDOWS,
            perturbations_deg=PERTURBATIONS_DEG,
            n_steps=N_STEPS,
            rng=rng,
            smoke_test=args.smoke_test,
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
            n_steps=N_STEPS,
            rng=rng,
            smoke_test=args.smoke_test,
        )
