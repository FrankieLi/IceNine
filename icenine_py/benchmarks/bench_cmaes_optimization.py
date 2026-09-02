#!/usr/bin/env python3
"""
CMA-ES orientation optimization benchmark.

Tests whether CMA-ES (Covariance Matrix Adaptation Evolution Strategy) can recover
a ground-truth orientation from a perturbed starting point, comparing:
  - Hard cost (VoxelCostFunction, binary overlap) — strongest signal
  - Differentiable cost (DifferentiableCostFunction, scale=2 ω±1) — softer signal

All CMA-ES runs use the same perturbation axes as the Adam benchmark (seed=42),
allowing direct comparison of results.

Algorithm:
  theta0 = so3_log(R_perturbed).numpy()   # axis-angle init (3-vector, radians)
  sigma0 = perturbation_rad / 3           # initial step size
  es = cma.CMAEvolutionStrategy(theta0, sigma0, {'maxiter': 500, ...})
  while not es.stop():
      candidates = es.ask()
      costs = [cost_fn(R(t), ...) for t in candidates]
      es.tell(candidates, costs)

Parameterization: theta in Lie algebra R^3 → R = matrix_exp(skew(theta)) in SO(3).
CMA-ES receives plain Python floats (no gradients needed).

Outputs (icenine_py/benchmarks/):
  cmaes_opt_threevoxels.csv   — per-generation: voxel_idx, perturbation_deg, cost_fn,
                                  generation, n_evals, quality, misorientation_deg
  cmaes_opt_manygrains.csv    — same, 20 voxels
  cmaes_convergence_*.png     — quality vs. n_evals, hard vs diff, one panel per pert
  cmaes_misori_convergence_*.png — misorientation vs. n_evals
  cmaes_vs_adam_summary.png   — final misorientation: CMA-ES vs Adam side-by-side

Usage:
  cd icenine_py
  uv run python benchmarks/bench_cmaes_optimization.py --smoke-test --example threevoxels
  uv run python benchmarks/bench_cmaes_optimization.py --example threevoxels
  uv run python benchmarks/bench_cmaes_optimization.py --example manygrains
"""

import argparse
import csv
import math
import os
import sys
import time
from pathlib import Path
from typing import Dict, List, Optional

import cma
import numpy as np
import torch
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

project_root = Path(__file__).parent.parent.parent
sys.path.insert(0, str(project_root / "icenine_py"))

benchmark_dir = Path(__file__).parent

# ---------------------------------------------------------------------------
# Sweep parameters (match Adam benchmark for direct comparison)
# ---------------------------------------------------------------------------

PERTURBATIONS_DEG = [1.0, 2.0, 5.0]
CMA_MAXITER = 500
SEED = 42
N_VOXELS_MANY = 20
QUALITY_THRESHOLD = 0.1
MAX_SCAN = 500

COST_FN_LABELS = ["hard", "diff_s2_ow1"]
COST_FN_COLORS = {"hard": "k", "diff_s2_ow1": "C2"}
PERT_COLORS = {1.0: "C0", 2.0: "C1", 5.0: "C3"}


# ---------------------------------------------------------------------------
# SO(3) helpers (identical to bench_gradient_optimization.py)
# ---------------------------------------------------------------------------


def skew(theta: torch.Tensor) -> torch.Tensor:
    z = torch.zeros(1, dtype=theta.dtype, device=theta.device)
    row0 = torch.stack([z.squeeze(), -theta[2], theta[1]])
    row1 = torch.stack([theta[2], z.squeeze(), -theta[0]])
    row2 = torch.stack([-theta[1], theta[0], z.squeeze()])
    return torch.stack([row0, row1, row2])


def so3_log(R: np.ndarray) -> torch.Tensor:
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
    M = R1.T @ R2
    trace = np.trace(M)
    cos_angle = np.clip((trace - 1.0) / 2.0, -1.0, 1.0)
    return float(np.degrees(np.arccos(cos_angle)))


def rodrigues_np(axis: np.ndarray, angle_rad: float) -> np.ndarray:
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


def theta_to_R_np(theta_np: np.ndarray) -> np.ndarray:
    """Convert axis-angle 3-vector (numpy) → rotation matrix (numpy)."""
    theta_t = torch.tensor(theta_np, dtype=torch.float32)
    with torch.no_grad():
        R = torch.matrix_exp(skew(theta_t))
    return R.numpy()


# ---------------------------------------------------------------------------
# Setup (reuses pattern from bench_gradient_optimization.py)
# ---------------------------------------------------------------------------


def setup_example(example_dir: Path, basename: str):
    """Load hard cost fn, differentiable cost fn (s2 ω±1), and mic."""
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

    print("  Building MultiScaleImageStack scale=2 ω±1 ...")
    ms_s2_ow1 = MultiScaleImageStack(
        image_stack, [1, 4, 8], omega_window=1, _prebuilt_downsampled=shared_ds
    )
    diff_fn_s2_ow1 = DifferentiableCostFunction(
        simulator=simulator,
        detector_list=detector_list,
        range_map=range_map,
        image_stack=ms_s2_ow1,
        sample=sample,
        structure_list=structure_list,
    )

    cost_fns = {
        "hard": hard_fn,
        "diff_s2_ow1": diff_fn_s2_ow1,
    }
    return mic, cost_fns, _get_voxel_vertices


# ---------------------------------------------------------------------------
# Voxel selection (identical pattern to bench_gradient_optimization.py)
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
        raise RuntimeError(f"Only {len(candidates)} candidates found (need {n})")

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
# Core CMA-ES optimizer
# ---------------------------------------------------------------------------


def _eval_cost(cost_fn_label: str, cost_fn, voxel, vertices, theta_np: np.ndarray) -> float:
    """Evaluate cost for a single candidate theta. Returns float."""
    theta_t = torch.tensor(theta_np, dtype=torch.float32)
    with torch.no_grad():
        R = torch.matrix_exp(skew(theta_t))
        if cost_fn_label == "hard":
            R_np = R.numpy()
            info = cost_fn.evaluate(
                orientation=R_np,
                voxel_vertices=vertices,
                phase_index=voxel.phase,
            )
            return float(info.cost)
        else:
            info = cost_fn.evaluate(R, vertices, phase_index=voxel.phase, scale=2)
            return float(info.cost.item())


def run_cmaes(
    cost_fn_label: str, cost_fn, voxel, vertices, R_init: np.ndarray, maxiter: int = CMA_MAXITER
) -> Dict:
    """Run CMA-ES from R_init.

    Returns dict with keys:
      quality_history: (n_generations,) best quality per generation
      misorientation_history: (n_generations,) best misorientation per generation
      n_evals_history: (n_generations,) cumulative evaluations
      R_final: (3,3) ndarray — best rotation found
      n_evals_total: int
      n_generations: int
      converged: bool — True if CMA-ES stopped due to tolerance (not maxiter)
    """
    R_gt = voxel.orientation

    theta0 = so3_log(R_init).numpy().astype(np.float64)
    # sigma0: initial step size = |theta0| (the perturbation magnitude).
    # The basin around ground truth is ~0.5° wide; we start sigma0 at the
    # full perturbation distance so the population can spread back toward truth.
    # CMA-ES will shrink sigma as it converges.
    sigma0 = float(np.linalg.norm(theta0))
    if sigma0 < 1e-4:
        sigma0 = 1e-4  # minimum step size to avoid degenerate starts

    opts = {
        "maxiter": maxiter,
        "tolx": 1e-5,  # ~0.0006° — stop when step size converges
        "tolfun": 0,  # disable within-generation flat stop
        "tolfunhist": 0,  # disable cross-generation flat stop
        "tolflatfitness": maxiter,  # allow maxiter flat-fitness generations
        "tolstagnation": maxiter,  # don't stop on stagnation
        "verbose": -9,  # suppress all CMA-ES stdout
        "seed": SEED,
    }

    es = cma.CMAEvolutionStrategy(theta0, sigma0, opts)

    quality_history = []
    misori_history = []
    n_evals_history = []
    best_cost = np.inf
    best_theta = theta0.copy()

    while not es.stop():
        solutions = es.ask()
        costs = [_eval_cost(cost_fn_label, cost_fn, voxel, vertices, s) for s in solutions]
        es.tell(solutions, costs)

        # Track best so far this generation
        gen_best_idx = int(np.argmin(costs))
        gen_best_cost = costs[gen_best_idx]
        if gen_best_cost < best_cost:
            best_cost = gen_best_cost
            best_theta = solutions[gen_best_idx].copy()

        R_best = theta_to_R_np(best_theta)
        quality_history.append(1.0 - best_cost)
        misori_history.append(misorientation_deg(R_gt, R_best))
        n_evals_history.append(es.result.evaluations)

    R_final = theta_to_R_np(best_theta)
    stop_reason = es.stop()
    converged = not any(k in stop_reason for k in ("maxiter", "maxfevals"))

    return {
        "quality_history": np.array(quality_history),
        "misorientation_history": np.array(misori_history),
        "n_evals_history": np.array(n_evals_history),
        "R_final": R_final,
        "n_evals_total": int(es.result.evaluations),
        "n_generations": len(quality_history),
        "converged": converged,
        "stop_reason": str(stop_reason),
    }


# ---------------------------------------------------------------------------
# Per-voxel sweep
# ---------------------------------------------------------------------------


def sweep_voxel(
    voxel,
    vidx: int,
    cost_fns: Dict,
    get_vertices,
    perturbations_deg: List[float],
    rng: np.random.Generator,
    cost_fn_labels: List[str],
    maxiter: int,
    smoke_test: bool = False,
) -> List[Dict]:
    """Run CMA-ES for all (cost_fn, perturbation) combos for one voxel."""
    vertices = get_vertices(voxel)
    rows = []

    # One random axis per perturbation (same seed logic as Adam benchmark)
    axes = [random_unit_axis(rng) for _ in perturbations_deg]

    for pi, (pert_deg, axis) in enumerate(zip(perturbations_deg, axes)):
        pert_rad = pert_deg * math.pi / 180.0
        R_pert = rodrigues_np(axis, pert_rad) @ voxel.orientation

        for cf_label in cost_fn_labels:
            t0 = time.perf_counter()
            result = run_cmaes(cf_label, cost_fns[cf_label], voxel, vertices, R_pert, maxiter)
            elapsed = time.perf_counter() - t0

            final_misori = result["misorientation_history"][-1]
            final_quality = result["quality_history"][-1]
            print(
                f"    vox {vidx:5d}  pert={pert_deg:.0f}°  {cf_label:<14s}"
                f"  final_misori={final_misori:.3f}°  q={final_quality:.4f}"
                f"  gens={result['n_generations']}  evals={result['n_evals_total']}"
                f"  {'CONVERGED' if result['converged'] else 'maxiter'}"
                f"  ({elapsed:.0f}s)",
                flush=True,
            )

            for gi, (q, m, ne) in enumerate(
                zip(
                    result["quality_history"],
                    result["misorientation_history"],
                    result["n_evals_history"],
                )
            ):
                rows.append(
                    {
                        "voxel_idx": vidx,
                        "perturbation_deg": pert_deg,
                        "cost_fn": cf_label,
                        "generation": gi,
                        "n_evals": int(ne),
                        "quality": float(q),
                        "misorientation_deg": float(m),
                    }
                )

    return rows


# ---------------------------------------------------------------------------
# CSV save
# ---------------------------------------------------------------------------


def save_csv(rows: List[Dict], out_path: Path):
    fieldnames = [
        "voxel_idx",
        "perturbation_deg",
        "cost_fn",
        "generation",
        "n_evals",
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
    out_path_q: Path,
    out_path_m: Path,
    title_prefix: str,
    perturbations_deg: List[float],
    cost_fn_labels: List[str],
):
    """Quality and misorientation vs. n_evals — one panel per perturbation."""
    import pandas as pd

    df = pd.read_csv(csv_path)
    n_pert = len(perturbations_deg)

    for metric, ylabel, out_path in [
        ("quality", "Quality (0–1)", out_path_q),
        ("misorientation_deg", "Misorientation (°)", out_path_m),
    ]:
        fig, axes = plt.subplots(1, n_pert, figsize=(6 * n_pert, 5), sharey=(metric == "quality"))
        if n_pert == 1:
            axes = [axes]

        for ai, pert in enumerate(perturbations_deg):
            ax = axes[ai]
            sub = df[df["perturbation_deg"] == pert]

            for cf in cost_fn_labels:
                sel = sub[sub["cost_fn"] == cf]
                if sel.empty:
                    continue
                # Mean over voxels at each generation (use n_evals as x-axis)
                mean_v = sel.groupby("n_evals")[metric].mean()
                std_v = sel.groupby("n_evals")[metric].std().fillna(0)
                color = COST_FN_COLORS.get(cf, "gray")
                ax.plot(mean_v.index, mean_v.values, color=color, lw=2.0, label=cf)
                ax.fill_between(
                    mean_v.index, mean_v - std_v, mean_v + std_v, color=color, alpha=0.15
                )

            ax.set_xlabel("Function evaluations")
            ax.set_ylabel(ylabel if ai == 0 else "")
            ax.set_title(f"{title_prefix}\nPerturbation = {pert:.0f}°")
            ax.legend(fontsize=9)
            ax.grid(True, alpha=0.3)
            if metric == "quality":
                ax.set_ylim(-0.02, 1.05)
            else:
                ax.axhline(0, color="gray", lw=0.8, ls="--")

        fig.tight_layout()
        fig.savefig(out_path, dpi=150)
        plt.close(fig)
        print(f"Saved: {out_path}")


def vs_adam_summary(
    cmaes_csv: Path,
    adam_csv: Path,
    out_path: Path,
    title_prefix: str,
    perturbations_deg: List[float],
    cost_fn_labels: List[str],
):
    """Side-by-side bar chart: final misorientation for CMA-ES vs Adam."""
    import pandas as pd

    df_c = pd.read_csv(cmaes_csv)
    # Get final generation per (voxel, perturbation, cost_fn)
    df_c_final = (
        df_c.sort_values("generation")
        .groupby(["voxel_idx", "perturbation_deg", "cost_fn"])
        .last()
        .reset_index()
    )

    df_a = None
    adam_tags = []
    if adam_csv.exists():
        df_a = pd.read_csv(adam_csv)
        df_a_final = (
            df_a.sort_values("step")
            .groupby(["voxel_idx", "perturbation_deg", "scale", "omega_window"])
            .last()
            .reset_index()
        )
        # Use best Adam config (s=2, ow=1) for comparison
        df_a_best = df_a_final[(df_a_final["scale"] == 2) & (df_a_final["omega_window"] == 1)]
        adam_tags = [("adam_s2_ow1", df_a_best)]

    n_pert = len(perturbations_deg)
    fig, axes = plt.subplots(1, n_pert, figsize=(5 * n_pert, 5), sharey=True)
    if n_pert == 1:
        axes = [axes]

    all_labels = [cf for cf in cost_fn_labels] + [tag for tag, _ in adam_tags]
    x = np.arange(len(all_labels))
    width = 0.6

    for ai, pert in enumerate(perturbations_deg):
        ax = axes[ai]
        means = []
        stds = []

        for cf in cost_fn_labels:
            sub = df_c_final[
                (df_c_final["perturbation_deg"] == pert) & (df_c_final["cost_fn"] == cf)
            ]["misorientation_deg"]
            means.append(float(sub.mean()) if len(sub) > 0 else np.nan)
            stds.append(float(sub.std()) if len(sub) > 1 else 0.0)

        for tag, df_tag in adam_tags:
            sub = df_tag[df_tag["perturbation_deg"] == pert]["misorientation_deg"]
            means.append(float(sub.mean()) if len(sub) > 0 else np.nan)
            stds.append(float(sub.std()) if len(sub) > 1 else 0.0)

        colors = [COST_FN_COLORS.get(cf, "C4") for cf in cost_fn_labels]
        colors += ["C5"] * len(adam_tags)

        bars = ax.bar(
            x, means, width, yerr=stds, color=colors, capsize=4, error_kw={"elinewidth": 1.2}
        )
        ax.set_xticks(x)
        ax.set_xticklabels(all_labels, rotation=15, ha="right", fontsize=8)
        ax.set_ylabel("Final misorientation (°)" if ai == 0 else "")
        ax.set_title(f"{title_prefix}\nPerturbation = {pert:.0f}°")
        ax.axhline(0.1, color="green", lw=0.8, ls="--", label="0.1° target")
        ax.axhline(1.0, color="orange", lw=0.8, ls="--", label="1° threshold")
        ax.legend(fontsize=7)
        ax.grid(True, alpha=0.3, axis="y")

        # Annotate bars
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
    perturbations_deg: List[float],
    cost_fn_labels: List[str],
    rng: np.random.Generator,
    maxiter: int,
    smoke_test: bool = False,
):
    print(f"\n{'='*60}")
    print(f"Example: {label}")
    print(f"{'='*60}")

    if smoke_test:
        n_voxels = min(n_voxels or 1, 1)
        maxiter = 10
        perturbations_deg = [2.0]
        cost_fn_labels = cost_fn_labels[:1]  # just hard
        print("  [SMOKE TEST: 1 voxel, maxiter=10, pert=2°, hard cost only]")

    mic, cost_fns, get_vertices = setup_example(example_dir, basename)

    # Hard fn needed for voxel selection
    hard_fn = cost_fns["hard"]

    if n_voxels is None:
        # Use all voxels that meet threshold
        voxel_list = []
        for idx, voxel in enumerate(mic.voxels):
            vertices = get_vertices(voxel)
            info = hard_fn.evaluate(
                orientation=voxel.orientation,
                voxel_vertices=get_vertices(voxel),
                phase_index=voxel.phase,
            )
            if info.quality > QUALITY_THRESHOLD:
                voxel_list.append((voxel, idx))
                print(f"    idx={idx}  hard_q={info.quality:.4f}  [selected]")
        print(f"  Using {len(voxel_list)} voxels.")
    else:
        voxel_list = select_voxels(mic, hard_fn, get_vertices, n_voxels, rng)

    total = len(voxel_list) * len(perturbations_deg) * len(cost_fn_labels)
    print(
        f"\nRunning {len(voxel_list)} voxels × {len(perturbations_deg)} perts"
        f" × {len(cost_fn_labels)} cost fns = {total} CMA-ES runs (maxiter={maxiter}) ..."
    )

    all_rows = []
    t_start = time.perf_counter()

    for vi, (voxel, vidx) in enumerate(voxel_list):
        print(f"\n  Voxel {vi+1}/{len(voxel_list)}  (mic idx={vidx}) ...")
        rows = sweep_voxel(
            voxel,
            vidx,
            cost_fns,
            get_vertices,
            perturbations_deg=perturbations_deg,
            rng=rng,
            cost_fn_labels=cost_fn_labels,
            maxiter=maxiter,
            smoke_test=smoke_test,
        )
        all_rows.extend(rows)

    elapsed = time.perf_counter() - t_start
    print(f"\nTotal time: {elapsed/60:.1f} min  ({len(all_rows)} rows)")

    tag = label.lower().replace(" ", "_").replace(".", "")
    csv_path = benchmark_dir / f"cmaes_opt_{tag}.csv"
    save_csv(all_rows, csv_path)

    convergence_plot(
        csv_path,
        benchmark_dir / f"cmaes_quality_convergence_{tag}.png",
        benchmark_dir / f"cmaes_misori_convergence_{tag}.png",
        title_prefix=label,
        perturbations_deg=perturbations_deg,
        cost_fn_labels=cost_fn_labels,
    )

    adam_csv = benchmark_dir / f"grad_opt_{tag}.csv"
    vs_adam_summary(
        csv_path,
        adam_csv,
        benchmark_dir / f"cmaes_vs_adam_summary_{tag}.png",
        title_prefix=label,
        perturbations_deg=perturbations_deg,
        cost_fn_labels=cost_fn_labels,
    )


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="CMA-ES orientation benchmark")
    parser.add_argument(
        "--smoke-test",
        action="store_true",
        help="Quick sanity: 1 voxel, maxiter=10, pert=2°, hard only",
    )
    parser.add_argument("--example", choices=["threevoxels", "manygrains", "both"], default="both")
    parser.add_argument(
        "--maxiter",
        type=int,
        default=CMA_MAXITER,
        help=f"CMA-ES max generations (default {CMA_MAXITER})",
    )
    args = parser.parse_args()

    rng = np.random.default_rng(seed=SEED)

    three_dir = project_root / "Examples" / "Example2.ThreeVoxels"
    many_dir = project_root / "Examples" / "Example2.ManyGrains"

    if args.example in ("threevoxels", "both"):
        run_example(
            label="ThreeVoxels",
            example_dir=three_dir,
            basename="3Grains.sim",
            n_voxels=None,
            perturbations_deg=PERTURBATIONS_DEG,
            cost_fn_labels=COST_FN_LABELS,
            rng=rng,
            maxiter=args.maxiter,
            smoke_test=args.smoke_test,
        )

    if args.example in ("manygrains", "both"):
        run_example(
            label="ManyGrains",
            example_dir=many_dir,
            basename="500Grains.sim",
            n_voxels=N_VOXELS_MANY,
            perturbations_deg=PERTURBATIONS_DEG,
            cost_fn_labels=COST_FN_LABELS,
            rng=rng,
            maxiter=args.maxiter,
            smoke_test=args.smoke_test,
        )
