#!/usr/bin/env python3
"""
Hybrid Riemannian Adam + MC-restart optimizer head-to-head benchmark.

Compares the hybrid RiemannianAdamOptimizer (gradient descent with MC restarts)
against the pure MCOptimizer on the FindOptimal phase, using the actual
reconstruction pipeline.

Perturbations: {0.5°, 1°, 2°, 3°, 5°}, 20 voxels, RNG seed=42.
Both optimizers run from the same perturbed starting orientation per voxel.
Time profiling breaks down the hybrid run into Adam time vs. hard-eval time.

Usage:
  cd /Users/sfli/Research/IceNine/icenine_py
  uv sync --extra riemannian
  uv run python benchmarks/bench_hybrid_optimizer.py --smoke-test --example threevoxels
  uv run python benchmarks/bench_hybrid_optimizer.py --example threevoxels
  uv run python benchmarks/bench_hybrid_optimizer.py --example manygrains

Outputs (icenine_py/benchmarks/):
  bench_hybrid_{example}.csv              — one row per (voxel, perturbation, optimizer)
  bench_hybrid_success_rate_{example}.png — success rate vs. perturbation
  bench_hybrid_wall_time_{example}.png    — wall time vs. perturbation (stacked for hybrid)
  bench_hybrid_scatter_{example}.png      — MC misori vs. hybrid misori per run
"""

import argparse
import csv
import math
import os
import sys
import time
from pathlib import Path
from typing import Dict, List, Optional, Tuple

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

PERTURBATIONS_DEG = [0.5, 1.0, 2.0, 3.0, 5.0]
SEED = 42
N_VOXELS = 20
QUALITY_THRESHOLD = 0.1
MAX_SCAN = 500
FIXED_SCALE = 2       # scale index 2 = 8× downsample
FIXED_OW = 1          # omega_window for MultiScaleImageStack
ADAM_N_STEPS = 100
ADAM_LR = 1e-4
ADAM_BETA1 = 0.9
ADAM_BETA2 = 0.999
MC_MAX_STEPS = 3500
MC_RESTARTS = 2
SUCCESS_THRESHOLD_DEG = 0.5   # misorientation below this = success


# ---------------------------------------------------------------------------
# SO(3) helpers
# ---------------------------------------------------------------------------

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
    eta_limit = exp_setup.get_eta_limit()

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
        eta_limit=eta_limit,
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
        eta_limit=eta_limit,
    )

    return mic, hard_fn, diff_fn, _get_voxel_vertices


# ---------------------------------------------------------------------------
# Voxel selection
# ---------------------------------------------------------------------------

def select_voxels(mic, hard_fn, get_vertices, n: int, rng: np.random.Generator):
    from icenine.geometry import matrix_to_euler

    print(f"  Scanning up to {MAX_SCAN} voxels for hard quality > {QUALITY_THRESHOLD} ...")
    candidates = []
    for idx, voxel in enumerate(mic.voxels[:MAX_SCAN]):
        vertices = get_vertices(voxel)
        info = hard_fn.evaluate(
            orientation=voxel.orientation,
            voxel_vertices=vertices,
            phase_index=voxel.phase,
        )
        if info.quality > QUALITY_THRESHOLD:
            candidates.append((voxel, idx, info.quality, vertices))

    if len(candidates) == 0:
        raise RuntimeError("No qualifying voxels found")
    if len(candidates) < n:
        print(f"  Warning: only {len(candidates)} qualifying voxels; using all of them (requested {n})")
        n = len(candidates)

    chosen_positions = rng.choice(len(candidates), size=n, replace=False)
    chosen = [candidates[i] for i in sorted(chosen_positions)]

    print(f"  Selected {n} voxels from {len(candidates)} candidates:")
    for voxel, vidx, q, _ in chosen:
        euler = matrix_to_euler(voxel.orientation)
        print(f"    idx={vidx:5d}  φ1={euler[0]:7.2f}°  Φ={euler[1]:6.2f}°  "
              f"φ2={euler[2]:7.2f}°  hard_q={q:.4f}")

    return [(v, vi, verts) for v, vi, _, verts in chosen]


# ---------------------------------------------------------------------------
# Timed hybrid optimizer run
# ---------------------------------------------------------------------------

def run_hybrid_timed(
    hard_fn,
    diff_fn,
    voxel,
    vertices,
    R_start: np.ndarray,
    angular_box_side: float,
    n_steps: int,
    lr: float,
    scale: int,
    max_restarts: int,
    rng: np.random.Generator,
) -> Dict:
    """
    Run RiemannianAdamOptimizer and return result + per-phase timing.

    Returns dict with keys:
      R_final, final_misori_deg, final_quality,
      wall_time_sec, t_adam_sec, t_hard_eval_sec,
      n_hard_evals, n_diff_evals
    """
    import geoopt

    R_gt = voxel.orientation
    best_orientation = R_start.copy()
    t0_hard = time.perf_counter()
    best_info = hard_fn.evaluate(best_orientation, vertices, voxel.phase)
    t_hard_eval = time.perf_counter() - t0_hard
    best_cost = best_info.cost
    n_hard_evals = 1

    current_orientation = R_start.copy()
    t_adam_total = 0.0
    n_diff_evals = 0

    from icenine.sampling import QuaternionGrid, matrix_to_quaternion, quaternion_to_matrix, _quat_multiply
    grid_gen = QuaternionGrid()

    for restart in range(max_restarts + 1):
        manifold = geoopt.manifolds.Stiefel()
        R = geoopt.ManifoldParameter(
            torch.from_numpy(current_orientation).float(), manifold=manifold
        )
        optimizer = geoopt.optim.RiemannianAdam(
            [R], lr=lr, betas=(ADAM_BETA1, ADAM_BETA2)
        )

        t_adam_start = time.perf_counter()
        for step in range(n_steps + 1):
            if step > 0:
                optimizer.zero_grad()
            with torch.set_grad_enabled(step > 0):
                diff_info = diff_fn.evaluate(R, vertices, phase_index=voxel.phase, scale=scale)
            n_diff_evals += 1
            if step > 0 and diff_info.cost.requires_grad:
                diff_info.cost.backward()
                optimizer.step()
        t_adam_total += time.perf_counter() - t_adam_start

        # SVD re-orthogonalize
        candidate_np = R.detach().numpy()
        U, _, Vt = np.linalg.svd(candidate_np)
        candidate_np = U @ Vt

        t0_hard = time.perf_counter()
        hard_info = hard_fn.evaluate(candidate_np, vertices, voxel.phase)
        t_hard_eval += time.perf_counter() - t0_hard
        n_hard_evals += 1

        if hard_info.cost < best_cost:
            best_cost = hard_info.cost
            best_orientation = candidate_np.copy()
            best_info = hard_info

        if restart < max_restarts:
            half_box = angular_box_side / 2.0
            rx = rng.uniform(-half_box, half_box)
            ry = rng.uniform(-half_box, half_box)
            rz = rng.uniform(-half_box, half_box)
            perturb_q = grid_gen.get_near_identity_point(rx, ry, rz)
            best_q = matrix_to_quaternion(best_orientation)
            current_orientation = quaternion_to_matrix(_quat_multiply(perturb_q, best_q))

    final_misori = misorientation_deg(R_gt, best_orientation)
    final_quality = best_info.quality if hasattr(best_info, "quality") else float("nan")

    return {
        "R_final": best_orientation,
        "final_misori_deg": final_misori,
        "final_quality": float(final_quality),
        "t_adam_sec": t_adam_total,
        "t_hard_eval_sec": t_hard_eval,
        "n_hard_evals": n_hard_evals,
        "n_diff_evals": n_diff_evals,
    }


# ---------------------------------------------------------------------------
# MC run
# ---------------------------------------------------------------------------

def run_mc_timed(
    hard_fn,
    voxel,
    vertices,
    R_start: np.ndarray,
    angular_box_side: float,
    angular_step: float,
    max_mc_steps: int,
    max_restarts: int,
    rng: np.random.Generator,
) -> Dict:
    """
    Run MCOptimizer.optimize() and return result with timing + eval count.
    """
    from icenine.orientation_search import MCOptimizer

    mc_opt = MCOptimizer(
        cost_fn=hard_fn,
        voxel_vertices=vertices,
        phase_index=voxel.phase,
        rng=rng,
    )

    evals_before = hard_fn.eval_count
    t0 = time.perf_counter()
    result = mc_opt.optimize(
        initial_orientation=R_start,
        angular_box_side=angular_box_side,
        angular_step=angular_step,
        max_mc_steps=max_mc_steps,
        max_restarts=max_restarts,
        max_convergence_cost=0.0,
    )
    wall_time = time.perf_counter() - t0
    n_hard_evals = hard_fn.eval_count - evals_before

    R_gt = voxel.orientation
    final_misori = misorientation_deg(R_gt, result.orientation)
    final_quality = result.overlap_info.quality if result.overlap_info is not None else float("nan")

    return {
        "R_final": result.orientation,
        "final_misori_deg": final_misori,
        "final_quality": float(final_quality),
        "wall_time_sec": wall_time,
        "n_hard_evals": n_hard_evals,
    }


# ---------------------------------------------------------------------------
# CSV schema
# ---------------------------------------------------------------------------

FIELDNAMES = [
    "voxel_idx", "perturbation_deg",
    "optimizer",
    "final_misori_deg", "final_quality",
    "wall_time_sec",
    "t_adam_sec", "t_hard_eval_sec",
    "n_hard_evals", "n_diff_evals",
]


# ---------------------------------------------------------------------------
# Main benchmark loop
# ---------------------------------------------------------------------------

def run_benchmark(
    label: str,
    example_dir: Path,
    basename: str,
    smoke_test: bool = False,
) -> None:
    print(f"\n{'='*60}")
    print(f"  Hybrid Optimizer Benchmark — {label}")
    print(f"{'='*60}")

    rng = np.random.default_rng(SEED)
    mic, hard_fn, diff_fn, get_vertices = setup_example(example_dir, basename)

    n_voxels = 3 if smoke_test else N_VOXELS
    perturbations = [1.0, 2.0] if smoke_test else PERTURBATIONS_DEG

    voxel_list = select_voxels(mic, hard_fn, get_vertices, n_voxels, rng)

    out_csv = benchmark_dir / f"bench_hybrid_{label}.csv"
    is_new = not out_csv.exists()
    csv_file = open(out_csv, "a", newline="")
    writer = csv.DictWriter(csv_file, fieldnames=FIELDNAMES, extrasaction="ignore")
    if is_new:
        writer.writeheader()

    total_runs = len(voxel_list) * len(perturbations)
    run_count = 0

    for voxel, vidx, vertices in voxel_list:
        for pert_deg in perturbations:
            run_count += 1
            pert_rad = math.radians(pert_deg)
            axis = random_unit_axis(rng)
            R_pert = rodrigues_np(axis, pert_rad) @ voxel.orientation

            # Angular box matches HP sweep convention: perturbation * 1.5
            # MC angular step from HP sweep best: step_frac=0.5
            angular_box_side = pert_rad * 1.5
            mc_step_frac = 0.5
            angular_step = angular_box_side * mc_step_frac

            # Hybrid run
            rng_hybrid = np.random.default_rng(SEED + vidx * 1000 + int(pert_deg * 10))
            t0_total = time.perf_counter()
            hybrid_result = run_hybrid_timed(
                hard_fn=hard_fn,
                diff_fn=diff_fn,
                voxel=voxel,
                vertices=vertices,
                R_start=R_pert,
                angular_box_side=angular_box_side,
                n_steps=ADAM_N_STEPS,
                lr=ADAM_LR,
                scale=FIXED_SCALE,
                max_restarts=MC_RESTARTS,
                rng=rng_hybrid,
            )
            hybrid_wall = time.perf_counter() - t0_total

            # MC run (same starting orientation, fresh rng with same seed)
            rng_mc = np.random.default_rng(SEED + vidx * 1000 + int(pert_deg * 10))
            mc_result = run_mc_timed(
                hard_fn=hard_fn,
                voxel=voxel,
                vertices=vertices,
                R_start=R_pert,
                angular_box_side=angular_box_side,
                angular_step=angular_step,
                max_mc_steps=MC_MAX_STEPS,
                max_restarts=MC_RESTARTS,
                rng=rng_mc,
            )

            # Write hybrid row
            writer.writerow({
                "voxel_idx": vidx,
                "perturbation_deg": pert_deg,
                "optimizer": "hybrid_adam",
                "final_misori_deg": hybrid_result["final_misori_deg"],
                "final_quality": hybrid_result["final_quality"],
                "wall_time_sec": hybrid_wall,
                "t_adam_sec": hybrid_result["t_adam_sec"],
                "t_hard_eval_sec": hybrid_result["t_hard_eval_sec"],
                "n_hard_evals": hybrid_result["n_hard_evals"],
                "n_diff_evals": hybrid_result["n_diff_evals"],
            })

            # Write MC row
            writer.writerow({
                "voxel_idx": vidx,
                "perturbation_deg": pert_deg,
                "optimizer": "mc_optimizer",
                "final_misori_deg": mc_result["final_misori_deg"],
                "final_quality": mc_result["final_quality"],
                "wall_time_sec": mc_result["wall_time_sec"],
                "t_adam_sec": float("nan"),
                "t_hard_eval_sec": float("nan"),
                "n_hard_evals": mc_result["n_hard_evals"],
                "n_diff_evals": 0,
            })
            csv_file.flush()

            hybrid_ok = hybrid_result["final_misori_deg"] < SUCCESS_THRESHOLD_DEG
            mc_ok = mc_result["final_misori_deg"] < SUCCESS_THRESHOLD_DEG
            print(
                f"  [{run_count}/{total_runs}] vox={vidx} pert={pert_deg:.1f}°  "
                f"hybrid: {hybrid_result['final_misori_deg']:.3f}° "
                f"{'OK' if hybrid_ok else '--'} "
                f"({hybrid_wall:.2f}s adam={hybrid_result['t_adam_sec']:.2f}s "
                f"hard={hybrid_result['t_hard_eval_sec']:.2f}s "
                f"evals={hybrid_result['n_hard_evals']}H/{hybrid_result['n_diff_evals']}D)  "
                f"mc: {mc_result['final_misori_deg']:.3f}° "
                f"{'OK' if mc_ok else '--'} "
                f"({mc_result['wall_time_sec']:.2f}s evals={mc_result['n_hard_evals']})",
                flush=True,
            )

    csv_file.close()
    print(f"\n  Results written to {out_csv}")
    plot_results(label, out_csv)


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------

def plot_results(label: str, csv_path: Path) -> None:
    import csv as csv_mod

    rows = []
    with open(csv_path) as f:
        reader = csv_mod.DictReader(f)
        for row in reader:
            rows.append({
                "voxel_idx": int(row["voxel_idx"]),
                "perturbation_deg": float(row["perturbation_deg"]),
                "optimizer": row["optimizer"],
                "final_misori_deg": float(row["final_misori_deg"]),
                "final_quality": float(row["final_quality"]) if row["final_quality"] else float("nan"),
                "wall_time_sec": float(row["wall_time_sec"]),
                "t_adam_sec": float(row["t_adam_sec"]) if row["t_adam_sec"] else float("nan"),
                "t_hard_eval_sec": float(row["t_hard_eval_sec"]) if row["t_hard_eval_sec"] else float("nan"),
                "n_hard_evals": int(float(row["n_hard_evals"])),
            })

    perturbations = sorted(set(r["perturbation_deg"] for r in rows))
    optimizers = ["hybrid_adam", "mc_optimizer"]
    colors = {"hybrid_adam": "#2196F3", "mc_optimizer": "#FF5722"}
    labels = {"hybrid_adam": "Hybrid Adam", "mc_optimizer": "MC Optimizer"}

    def get_data(opt: str, pert: float, key: str) -> List[float]:
        return [r[key] for r in rows if r["optimizer"] == opt and r["perturbation_deg"] == pert]

    # Plot 1: Success rate vs perturbation
    fig, ax = plt.subplots(figsize=(8, 5))
    x = np.arange(len(perturbations))
    width = 0.35
    for i, opt in enumerate(optimizers):
        success_rates = []
        for pert in perturbations:
            misori = get_data(opt, pert, "final_misori_deg")
            rate = sum(m < SUCCESS_THRESHOLD_DEG for m in misori) / len(misori) if misori else 0.0
            success_rates.append(100.0 * rate)
        ax.bar(x + (i - 0.5) * width, success_rates, width, label=labels[opt],
               color=colors[opt], alpha=0.8)
    ax.set_xlabel("Perturbation (°)")
    ax.set_ylabel(f"Success rate (%) [misori < {SUCCESS_THRESHOLD_DEG}°]")
    ax.set_title(f"Hybrid Adam vs. MC — Success Rate ({label})")
    ax.set_xticks(x)
    ax.set_xticklabels([f"{p:.1f}°" for p in perturbations])
    ax.set_ylim(0, 105)
    ax.legend()
    ax.grid(axis="y", alpha=0.3)
    out_path = benchmark_dir / f"bench_hybrid_success_rate_{label}.png"
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {out_path.name}")

    # Plot 2: Wall time vs perturbation (stacked for hybrid: adam + hard_eval)
    fig, ax = plt.subplots(figsize=(8, 5))
    hybrid_adam_times = [
        np.mean([r["t_adam_sec"] for r in rows
                 if r["optimizer"] == "hybrid_adam" and r["perturbation_deg"] == pert
                 and not math.isnan(r["t_adam_sec"])])
        for pert in perturbations
    ]
    hybrid_hard_times = [
        np.mean([r["t_hard_eval_sec"] for r in rows
                 if r["optimizer"] == "hybrid_adam" and r["perturbation_deg"] == pert
                 and not math.isnan(r["t_hard_eval_sec"])])
        for pert in perturbations
    ]
    mc_times = [
        np.mean(get_data("mc_optimizer", pert, "wall_time_sec"))
        for pert in perturbations
    ]

    ax.bar(x - width / 2, hybrid_adam_times, width, label="Hybrid: Adam steps",
           color=colors["hybrid_adam"], alpha=0.8)
    ax.bar(x - width / 2, hybrid_hard_times, width, bottom=hybrid_adam_times,
           label="Hybrid: hard eval", color="#90CAF9", alpha=0.8)
    ax.bar(x + width / 2, mc_times, width, label=labels["mc_optimizer"],
           color=colors["mc_optimizer"], alpha=0.8)
    ax.set_xlabel("Perturbation (°)")
    ax.set_ylabel("Mean wall time (s)")
    ax.set_title(f"Hybrid Adam vs. MC — Wall Time ({label})")
    ax.set_xticks(x)
    ax.set_xticklabels([f"{p:.1f}°" for p in perturbations])
    ax.legend()
    ax.grid(axis="y", alpha=0.3)
    out_path = benchmark_dir / f"bench_hybrid_wall_time_{label}.png"
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {out_path.name}")

    # Plot 3: Scatter — MC misori vs hybrid misori per run
    fig, axes = plt.subplots(1, len(perturbations), figsize=(4 * len(perturbations), 4),
                             sharey=True, sharex=True)
    if len(perturbations) == 1:
        axes = [axes]
    max_misori = max(
        max((r["final_misori_deg"] for r in rows), default=5.0), 1.0
    )
    for ax, pert in zip(axes, perturbations):
        hybrid_misori = get_data("hybrid_adam", pert, "final_misori_deg")
        mc_misori = get_data("mc_optimizer", pert, "final_misori_deg")
        n = min(len(hybrid_misori), len(mc_misori))
        ax.scatter(mc_misori[:n], hybrid_misori[:n], alpha=0.6, s=20,
                   color="#5C6BC0", edgecolors="none")
        lim = max_misori * 1.05
        ax.plot([0, lim], [0, lim], "k--", lw=0.8, alpha=0.5, label="equal")
        ax.axhline(SUCCESS_THRESHOLD_DEG, color=colors["hybrid_adam"], lw=0.8,
                   alpha=0.5, linestyle=":")
        ax.axvline(SUCCESS_THRESHOLD_DEG, color=colors["mc_optimizer"], lw=0.8,
                   alpha=0.5, linestyle=":")
        ax.set_xlim(0, lim)
        ax.set_ylim(0, lim)
        ax.set_title(f"pert={pert:.1f}°")
        ax.set_xlabel("MC misori (°)")
        ax.grid(alpha=0.3)
    axes[0].set_ylabel("Hybrid misori (°)")
    fig.suptitle(f"MC vs. Hybrid misorientation per run ({label})", y=1.02)
    out_path = benchmark_dir / f"bench_hybrid_scatter_{label}.png"
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {out_path.name}")


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

EXAMPLES = {
    "threevoxels": (
        Path(__file__).parent.parent.parent / "Examples" / "Example2.ThreeVoxels",
        "3Grains.sim",
    ),
    "manygrains": (
        Path(__file__).parent.parent.parent / "Examples" / "Example2.ManyGrains",
        "500Grains.sim",
    ),
}


def main():
    parser = argparse.ArgumentParser(description="Hybrid Adam vs. MC optimizer benchmark")
    parser.add_argument("--example", choices=list(EXAMPLES.keys()),
                        default="threevoxels",
                        help="Which example dataset to use")
    parser.add_argument("--smoke-test", action="store_true",
                        help="Quick smoke test: 3 voxels, 2 perturbations")
    parser.add_argument("--plots-only", action="store_true",
                        help="Regenerate plots from existing CSV without running benchmark")
    args = parser.parse_args()

    example_dir, basename = EXAMPLES[args.example]
    label = args.example

    if args.plots_only:
        csv_path = benchmark_dir / f"bench_hybrid_{label}.csv"
        if not csv_path.exists():
            print(f"ERROR: {csv_path} does not exist. Run without --plots-only first.")
            sys.exit(1)
        plot_results(label, csv_path)
    else:
        run_benchmark(label, example_dir, basename, smoke_test=args.smoke_test)


if __name__ == "__main__":
    main()
