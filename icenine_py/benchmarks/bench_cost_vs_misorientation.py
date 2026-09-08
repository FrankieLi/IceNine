#!/usr/bin/env python3
"""
Cost vs. misorientation benchmark.

Evaluates hard cost (VoxelCostFunction) and differentiable cost
(DifferentiableCostFunction at 3 scales) over 31 misorientation angles
from 0° to 15°.

  1. Example2.ThreeVoxels  — voxel 0, 30 random axes per angle
     Variance = axis-to-axis spread at a single voxel

  2. Example2.ManyGrains   — 20 randomly selected voxels, 3 axes per voxel
     Variance = voxel-to-voxel spread; shows how different crystal
     orientations (Bragg peak configurations) change the landscape shape

Note on sampling: axes are drawn from np.random.randn(3), normalized.
Uniform on the sphere but NOT Haar-correct for SO(3) perturbations.
Adequate for landscape visualization.

Estimated runtime: ~2h (ThreeVoxels) + ~4h (ManyGrains) = ~6h total.

Usage:
  cd icenine_py
  uv run python benchmarks/bench_cost_vs_misorientation.py
"""

import math
import os
import time
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import torch

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------

project_root = Path(__file__).parent.parent.parent
benchmark_dir = Path(__file__).parent


# ---------------------------------------------------------------------------
# Physics helpers
# ---------------------------------------------------------------------------


def rodrigues(axis: np.ndarray, angle_rad: float) -> np.ndarray:
    """Rotation matrix via Rodrigues formula (axis must be unit vector)."""
    K = np.array(
        [
            [0, -axis[2], axis[1]],
            [axis[2], 0, -axis[0]],
            [-axis[1], axis[0], 0],
        ],
        dtype=np.float64,
    )
    return np.eye(3) + math.sin(angle_rad) * K + (1 - math.cos(angle_rad)) * (K @ K)


def random_unit_axes(n: int, rng: np.random.Generator) -> np.ndarray:
    """Draw n unit vectors uniformly on the sphere (not Haar-correct for SO(3))."""
    v = rng.standard_normal((n, 3))
    norms = np.linalg.norm(v, axis=1, keepdims=True)
    return v / norms


# ---------------------------------------------------------------------------
# Setup
# ---------------------------------------------------------------------------


def setup_example(example_dir: Path, basename: str):
    """Load cost functions, image stacks, and mic for one example."""
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

    print(f"  Loading images from {data_dir.name}/ ...")
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

    print("  Building MultiScaleImageStack [1x, 4x, 8x] ...")
    multi_stack = MultiScaleImageStack(image_stack, downsample_factors=[1, 4, 8])

    print("  Building MultiScaleImageStack [1x, 4x, 8x] with omega_window=1 ...")
    multi_stack_oblend = MultiScaleImageStack(
        image_stack, downsample_factors=[1, 4, 8], omega_window=1
    )

    hard_cost_fn = VoxelCostFunction(
        simulator=simulator,
        detector_list=detector_list,
        range_map=range_map,
        exp_data=exp_data,
        sample=sample,
        structure_list=structure_list,
        mode="hard",
    )

    diff_cost_fn = DifferentiableCostFunction(
        simulator=simulator,
        detector_list=detector_list,
        range_map=range_map,
        image_stack=multi_stack,
        sample=sample,
        structure_list=structure_list,
    )

    diff_cost_fn_oblend = DifferentiableCostFunction(
        simulator=simulator,
        detector_list=detector_list,
        range_map=range_map,
        image_stack=multi_stack_oblend,
        sample=sample,
        structure_list=structure_list,
    )

    return {
        "hard_cost_fn": hard_cost_fn,
        "diff_cost_fn": diff_cost_fn,
        "diff_cost_fn_oblend": diff_cost_fn_oblend,
        "mic": mic,
        "get_vertices": _get_voxel_vertices,
    }


# ---------------------------------------------------------------------------
# Voxel selection
# ---------------------------------------------------------------------------


def select_random_voxels(
    mic,
    hard_cost_fn,
    get_vertices,
    n_voxels: int,
    quality_threshold: float,
    max_scan: int,
    rng: np.random.Generator,
):
    """
    Scan up to max_scan voxels, collect those with hard quality > threshold,
    then randomly select n_voxels. Returns list of (voxel, idx) tuples.
    """
    from icenine.geometry import matrix_to_euler

    print(f"  Scanning up to {max_scan} voxels for quality > {quality_threshold} ...")
    candidates = []
    for idx, voxel in enumerate(mic.voxels[:max_scan]):
        vertices = get_vertices(voxel)
        info = hard_cost_fn.evaluate(
            orientation=voxel.orientation,
            voxel_vertices=vertices,
            phase_index=voxel.phase,
        )
        if info.quality > quality_threshold:
            candidates.append((voxel, idx, info.quality))
        if (idx + 1) % 50 == 0:
            print(
                f"    scanned {idx+1}/{max_scan}, {len(candidates)} qualifying so far", flush=True
            )

    if len(candidates) < n_voxels:
        raise RuntimeError(
            f"Only {len(candidates)} voxels with quality>{quality_threshold} "
            f"in first {max_scan}; need {n_voxels}"
        )

    chosen_idx = rng.choice(len(candidates), size=n_voxels, replace=False)
    chosen = [candidates[i] for i in sorted(chosen_idx)]

    print(f"  Selected {n_voxels} voxels (from {len(candidates)} candidates):")
    for voxel, vidx, q in chosen:
        euler = matrix_to_euler(voxel.orientation)
        print(
            f"    idx={vidx:5d}  phi1={euler[0]:7.2f}°  Phi={euler[1]:6.2f}°  "
            f"phi2={euler[2]:7.2f}°  hard_q={q:.4f}"
        )

    return [(v, vi) for v, vi, _ in chosen]


# ---------------------------------------------------------------------------
# Sweep
# ---------------------------------------------------------------------------


def sweep_one_voxel(
    voxel, hard_cost_fn, diff_cost_fn, get_vertices, angles_deg: np.ndarray, axes: np.ndarray
):
    """
    Evaluate all cost functions at each (angle, axis) pair for one voxel.
    axes: (n_axes, 3) — fixed set of rotation axes.
    Returns dict of arrays shape (n_angles, n_axes).
    """
    n_angles = len(angles_deg)
    n_axes = len(axes)
    hard = np.zeros((n_angles, n_axes))
    diff_s0 = np.zeros((n_angles, n_axes))
    diff_s1 = np.zeros((n_angles, n_axes))
    diff_s2 = np.zeros((n_angles, n_axes))

    base_orient = voxel.orientation
    vertices = get_vertices(voxel)

    for ai, angle_deg in enumerate(angles_deg):
        angle_rad = angle_deg * math.pi / 180.0
        for xi, axis in enumerate(axes):
            R = rodrigues(axis, angle_rad)
            perturbed = R @ base_orient
            perturbed_t = torch.from_numpy(perturbed).float()

            info_h = hard_cost_fn.evaluate(
                orientation=perturbed,
                voxel_vertices=vertices,
                phase_index=voxel.phase,
            )
            hard[ai, xi] = info_h.quality

            with torch.no_grad():
                for si, arr in enumerate([diff_s0, diff_s1, diff_s2]):
                    info_d = diff_cost_fn.evaluate(
                        perturbed_t, vertices, phase_index=voxel.phase, scale=si
                    )
                    arr[ai, xi] = info_d.quality.item()

    return {"hard": hard, "diff_s0": diff_s0, "diff_s1": diff_s1, "diff_s2": diff_s2}


# ---------------------------------------------------------------------------
# CSV / Plot — single voxel
# ---------------------------------------------------------------------------


def save_csv_single(path: Path, angles_deg: np.ndarray, results: dict):
    """Save per-sample CSV for single-voxel sweep."""
    n_angles, n_axes = results["hard"].shape
    rows = []
    for ai, angle in enumerate(angles_deg):
        for xi in range(n_axes):
            rows.append(
                f"{angle:.4f},{xi},"
                f"{results['hard'][ai, xi]:.6f},"
                f"{results['diff_s0'][ai, xi]:.6f},"
                f"{results['diff_s1'][ai, xi]:.6f},"
                f"{results['diff_s2'][ai, xi]:.6f},"
                f"{results['diff_s2_oblend'][ai, xi]:.6f}"
            )
    header = (
        "angle_deg,axis_idx,hard_quality,diff_quality_s0,"
        "diff_quality_s1,diff_quality_s2,diff_quality_s2_oblend"
    )
    with open(path, "w") as f:
        f.write(header + "\n")
        f.write("\n".join(rows) + "\n")
    print(f"  Saved CSV: {path}")


def save_plot_single(path: Path, angles_deg: np.ndarray, results: dict, title: str):
    """Plot quality vs misorientation for single voxel. Band = ±1 std over axes."""
    fig, ax = plt.subplots(figsize=(9, 5))

    styles = [
        ("hard", "Hard (binary overlap)", "k", "-", 2.0),
        ("diff_s0", "Diff s0 — 1× (2048²)", "C0", "-", 1.8),
        ("diff_s1", "Diff s1 — 4× (512²)", "C1", "--", 1.8),
        ("diff_s2", "Diff s2 — 8× (256², max-pool)", "C3", ":", 1.8),
        ("diff_s2_oblend", "Diff s2 — 8× + ω±1 blend", "C2", "-.", 1.8),
    ]

    for key, label, color, ls, lw in styles:
        arr = results[key]
        mean = arr.mean(axis=1)
        std = arr.std(axis=1)
        ax.plot(angles_deg, mean, color=color, ls=ls, lw=lw, label=label)
        ax.fill_between(angles_deg, mean - std, mean + std, color=color, alpha=0.15)

    ax.axvline(0, color="gray", lw=0.8, ls="--")
    ax.set_xlabel("Misorientation angle (degrees)")
    ax.set_ylabel("Quality (0–1)")
    ax.set_ylim(-0.02, 1.05)
    ax.set_xlim(angles_deg[0], angles_deg[-1])
    ax.set_title(title)
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)
    print(f"  Saved plot: {path}")


# ---------------------------------------------------------------------------
# CSV / Plot — multi-voxel
# ---------------------------------------------------------------------------


def save_csv_multi(path: Path, angles_deg: np.ndarray, voxel_indices: list, all_results: list):
    """
    Save multi-voxel CSV.
    Columns: angle_deg, voxel_idx, hard_quality, diff_quality_s0, s1, s2, s2_oblend
    Each row is the mean over axes for one (voxel, angle) pair.
    """
    rows = []
    for (voxel, vidx), res in zip(voxel_indices, all_results):
        for ai, angle in enumerate(angles_deg):
            rows.append(
                f"{angle:.4f},{vidx},"
                f"{res['hard'][ai].mean():.6f},"
                f"{res['diff_s0'][ai].mean():.6f},"
                f"{res['diff_s1'][ai].mean():.6f},"
                f"{res['diff_s2'][ai].mean():.6f},"
                f"{res['diff_s2_oblend'][ai].mean():.6f}"
            )
    header = (
        "angle_deg,voxel_idx,hard_quality,diff_quality_s0,"
        "diff_quality_s1,diff_quality_s2,diff_quality_s2_oblend"
    )
    with open(path, "w") as f:
        f.write(header + "\n")
        f.write("\n".join(rows) + "\n")
    print(f"  Saved CSV: {path}")


def save_plot_multi(path: Path, angles_deg: np.ndarray, all_results: list, title: str):
    """
    Plot multi-voxel cost landscape.
    For each cost function: thin semi-transparent lines per voxel (spaghetti),
    thick line = mean across voxels, shaded band = ±1 std across voxels.
    """
    styles = [
        ("hard", "Hard (binary overlap)", "k", "-", 2.2),
        ("diff_s0", "Diff s0 — 1× (2048²)", "C0", "-", 2.0),
        ("diff_s1", "Diff s1 — 4× (512²)", "C1", "--", 2.0),
        ("diff_s2", "Diff s2 — 8× (256², max-pool)", "C3", ":", 2.0),
        ("diff_s2_oblend", "Diff s2 — 8× + ω±1 blend", "C2", "-.", 2.0),
    ]

    # per_voxel_means[key] shape: (n_voxels, n_angles)
    per_voxel = {
        key: np.stack([res[key].mean(axis=1) for res in all_results], axis=0)
        for key in ("hard", "diff_s0", "diff_s1", "diff_s2", "diff_s2_oblend")
    }

    fig, ax = plt.subplots(figsize=(10, 6))

    for key, label, color, ls, lw in styles:
        arr = per_voxel[key]  # (n_voxels, n_angles)
        mean = arr.mean(axis=0)
        std = arr.std(axis=0)

        # Spaghetti: one thin line per voxel
        for vi in range(arr.shape[0]):
            ax.plot(angles_deg, arr[vi], color=color, ls=ls, lw=0.5, alpha=0.25)

        # Mean + std band
        ax.plot(angles_deg, mean, color=color, ls=ls, lw=lw, label=label)
        ax.fill_between(angles_deg, mean - std, mean + std, color=color, alpha=0.18)

    ax.axvline(0, color="gray", lw=0.8, ls="--")
    ax.set_xlabel("Misorientation angle (degrees)")
    ax.set_ylabel("Quality (0–1)")
    ax.set_ylim(-0.02, 1.05)
    ax.set_xlim(angles_deg[0], angles_deg[-1])
    ax.set_title(title)
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)
    print(f"  Saved plot: {path}")


# ---------------------------------------------------------------------------
# Run examples
# ---------------------------------------------------------------------------


def run_single_voxel(
    name: str,
    example_dir: Path,
    basename: str,
    voxel_idx: int,
    angles_deg: np.ndarray,
    n_axes: int,
    rng: np.random.Generator,
):
    from icenine.geometry import matrix_to_euler

    print(f"\n{'='*60}")
    print(f"Example: {name}  (single voxel, {n_axes} axes)")
    print(f"{'='*60}")

    setup = setup_example(example_dir, basename)
    mic, hard_fn, diff_fn, diff_fn_oblend, get_vertices = (
        setup["mic"],
        setup["hard_cost_fn"],
        setup["diff_cost_fn"],
        setup["diff_cost_fn_oblend"],
        setup["get_vertices"],
    )

    voxel = mic.voxels[voxel_idx]
    euler = matrix_to_euler(voxel.orientation)
    vertices = get_vertices(voxel)
    hard_info = hard_fn.evaluate(voxel.orientation, vertices, voxel.phase)
    print(
        f"  Voxel idx={voxel_idx}  "
        f"phi1={euler[0]:.2f}°  Phi={euler[1]:.2f}°  phi2={euler[2]:.2f}°  "
        f"hard_quality={hard_info.quality:.4f}"
    )

    axes = random_unit_axes(n_axes, rng)
    n_total = len(angles_deg) * n_axes
    print(f"  Sweeping {len(angles_deg)} angles × {n_axes} axes = {n_total} evals ...")
    t0 = time.perf_counter()

    hard_all = np.zeros((len(angles_deg), n_axes))
    diff_s0_all = np.zeros((len(angles_deg), n_axes))
    diff_s1_all = np.zeros((len(angles_deg), n_axes))
    diff_s2_all = np.zeros((len(angles_deg), n_axes))
    diff_s2ob_all = np.zeros((len(angles_deg), n_axes))

    base_orient = voxel.orientation
    for ai, angle_deg in enumerate(angles_deg):
        angle_rad = angle_deg * math.pi / 180.0
        for xi, axis in enumerate(axes):
            R = rodrigues(axis, angle_rad)
            perturbed = R @ base_orient
            perturbed_t = torch.from_numpy(perturbed).float()

            info_h = hard_fn.evaluate(perturbed, vertices, voxel.phase)
            hard_all[ai, xi] = info_h.quality

            with torch.no_grad():
                for si, arr in enumerate([diff_s0_all, diff_s1_all, diff_s2_all]):
                    info_d = diff_fn.evaluate(
                        perturbed_t, vertices, phase_index=voxel.phase, scale=si
                    )
                    arr[ai, xi] = info_d.quality.item()

                info_ob = diff_fn_oblend.evaluate(
                    perturbed_t, vertices, phase_index=voxel.phase, scale=2
                )
                diff_s2ob_all[ai, xi] = info_ob.quality.item()

        elapsed = time.perf_counter() - t0
        done = (ai + 1) * n_axes
        eta = elapsed / done * (n_total - done) if done > 0 else 0
        print(
            f"    angle={angle_deg:5.1f}°  "
            f"hard={hard_all[ai].mean():.4f}  "
            f"s0={diff_s0_all[ai].mean():.4f}  "
            f"s1={diff_s1_all[ai].mean():.4f}  "
            f"s2={diff_s2_all[ai].mean():.4f}  "
            f"s2ω={diff_s2ob_all[ai].mean():.4f}  "
            f"[{(ai+1)/len(angles_deg)*100:.0f}%  ETA {eta/60:.0f}m]",
            flush=True,
        )

    results = {
        "hard": hard_all,
        "diff_s0": diff_s0_all,
        "diff_s1": diff_s1_all,
        "diff_s2": diff_s2_all,
        "diff_s2_oblend": diff_s2ob_all,
    }

    tag = name.lower().replace(" ", "_")
    save_csv_single(benchmark_dir / f"misorientation_{tag}.csv", angles_deg, results)
    save_plot_single(
        benchmark_dir / f"misorientation_{tag}.png",
        angles_deg,
        results,
        f"{name}  —  voxel {voxel_idx}  ({n_axes} axes, ±1σ band)",
    )
    return results


def run_multi_voxel(
    name: str,
    example_dir: Path,
    basename: str,
    n_voxels: int,
    angles_deg: np.ndarray,
    n_axes: int,
    rng: np.random.Generator,
    quality_threshold: float = 0.1,
    max_scan: int = 500,
):
    print(f"\n{'='*60}")
    print(f"Example: {name}  ({n_voxels} voxels, {n_axes} axes each)")
    print(f"{'='*60}")

    setup = setup_example(example_dir, basename)
    mic, hard_fn, diff_fn, diff_fn_oblend, get_vertices = (
        setup["mic"],
        setup["hard_cost_fn"],
        setup["diff_cost_fn"],
        setup["diff_cost_fn_oblend"],
        setup["get_vertices"],
    )

    voxel_list = select_random_voxels(
        mic,
        hard_fn,
        get_vertices,
        n_voxels=n_voxels,
        quality_threshold=quality_threshold,
        max_scan=max_scan,
        rng=rng,
    )

    # Use the same set of axes for all voxels at each angle
    axes = random_unit_axes(n_axes, rng)
    n_total = n_voxels * len(angles_deg) * n_axes
    print(
        f"\n  Sweeping {n_voxels} voxels × {len(angles_deg)} angles × {n_axes} axes "
        f"= {n_total} evals ..."
    )

    all_results = []
    t0 = time.perf_counter()
    evals_done = 0

    for vi, (voxel, vidx) in enumerate(voxel_list):
        print(f"\n  --- Voxel {vi+1}/{n_voxels} (idx={vidx}) ---", flush=True)
        hard_all = np.zeros((len(angles_deg), n_axes))
        diff_s0_all = np.zeros((len(angles_deg), n_axes))
        diff_s1_all = np.zeros((len(angles_deg), n_axes))
        diff_s2_all = np.zeros((len(angles_deg), n_axes))
        diff_s2ob_all = np.zeros((len(angles_deg), n_axes))

        base_orient = voxel.orientation
        vertices = get_vertices(voxel)

        for ai, angle_deg in enumerate(angles_deg):
            angle_rad = angle_deg * math.pi / 180.0
            for xi, axis in enumerate(axes):
                R = rodrigues(axis, angle_rad)
                perturbed = R @ base_orient
                perturbed_t = torch.from_numpy(perturbed).float()

                info_h = hard_fn.evaluate(perturbed, vertices, voxel.phase)
                hard_all[ai, xi] = info_h.quality

                with torch.no_grad():
                    for si, arr in enumerate([diff_s0_all, diff_s1_all, diff_s2_all]):
                        info_d = diff_fn.evaluate(
                            perturbed_t, vertices, phase_index=voxel.phase, scale=si
                        )
                        arr[ai, xi] = info_d.quality.item()

                    info_ob = diff_fn_oblend.evaluate(
                        perturbed_t, vertices, phase_index=voxel.phase, scale=2
                    )
                    diff_s2ob_all[ai, xi] = info_ob.quality.item()

            evals_done += n_axes
            elapsed = time.perf_counter() - t0
            eta = elapsed / evals_done * (n_total - evals_done) if evals_done > 0 else 0
            print(
                f"    angle={angle_deg:5.1f}°  "
                f"hard={hard_all[ai].mean():.4f}  "
                f"s0={diff_s0_all[ai].mean():.4f}  "
                f"s1={diff_s1_all[ai].mean():.4f}  "
                f"s2={diff_s2_all[ai].mean():.4f}  "
                f"s2ω={diff_s2ob_all[ai].mean():.4f}  "
                f"[vox {vi+1}/{n_voxels}  ETA {eta/60:.0f}m]",
                flush=True,
            )

        all_results.append(
            {
                "hard": hard_all,
                "diff_s0": diff_s0_all,
                "diff_s1": diff_s1_all,
                "diff_s2": diff_s2_all,
                "diff_s2_oblend": diff_s2ob_all,
            }
        )

    tag = name.lower().replace(" ", "_")
    save_csv_multi(
        benchmark_dir / f"misorientation_{tag}.csv",
        angles_deg,
        voxel_list,
        all_results,
    )
    save_plot_multi(
        benchmark_dir / f"misorientation_{tag}.png",
        angles_deg,
        all_results,
        f"{name}  —  {n_voxels} random voxels  (spaghetti = per-voxel, band = ±1σ)",
    )
    return all_results


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    angles_deg = np.linspace(0, 15, 31)  # 0.5° steps
    rng = np.random.default_rng(seed=42)

    # Example 1: ThreeVoxels — single voxel, axis-to-axis variance
    run_single_voxel(
        name="ThreeVoxels",
        example_dir=project_root / "Examples" / "Example2.ThreeVoxels",
        basename="3Grains.sim",
        voxel_idx=0,
        angles_deg=angles_deg,
        n_axes=30,
        rng=rng,
    )

    # Example 2: ManyGrains — 20 random voxels, voxel-to-voxel variance
    run_multi_voxel(
        name="ManyGrains",
        example_dir=project_root / "Examples" / "Example2.ManyGrains",
        basename="500Grains.sim",
        n_voxels=20,
        angles_deg=angles_deg,
        n_axes=3,
        rng=rng,
        quality_threshold=0.1,
        max_scan=500,
    )

    print("\nDone.")
