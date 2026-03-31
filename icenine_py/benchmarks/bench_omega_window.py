"""
Compare omega_window=0, 1, 2 at scale 2 (8× downsampled).

Uses the same ManyGrains setup as bench_cost_vs_misorientation.py but
only evaluates DifferentiableCostFunction at scale=2 for three omega_window
values. Runs on 5 randomly-selected voxels with 5 axes each to keep runtime
short (~5–10 min).

Outputs:
  benchmarks/omega_window_comparison.png
  benchmarks/omega_window_comparison_gradient.png
"""

import math
import os
import sys
import time
from pathlib import Path

import numpy as np
import torch
import matplotlib.pyplot as plt

project_root = Path(__file__).parent.parent.parent
sys.path.insert(0, str(project_root / "icenine_py"))

from icenine.config_file import ConfigFile
from icenine.experiment_setup import XDMExperimentSetup
from icenine.simulation import Simulation
from icenine.sample import Sample
from icenine.cost_functions import VoxelCostFunction
from icenine.experimental_data import ExperimentalData
from icenine.mic_file import MicFile
from icenine.differentiable_cost import (
    SparseImageStack,
    MultiScaleImageStack,
    DifferentiableCostFunction,
)

benchmark_dir = Path(__file__).parent

OMEGA_WINDOWS = [0, 1, 2]
COLORS = {"0": "C3", "1": "C2", "2": "C5"}
LINESTYLES = {"0": ":", "1": "-.", "2": "--"}

EXAMPLE_DIR = project_root / "Examples" / "Example2.ManyGrains"
BASENAME = "500Grains.sim"
N_VOXELS = 5
N_AXES = 5
ANGLES_DEG = np.linspace(0, 15, 31)
SEED = 42


def rodrigues(axis: np.ndarray, angle_rad: float) -> np.ndarray:
    K = np.array([[0, -axis[2], axis[1]],
                  [axis[2], 0, -axis[0]],
                  [-axis[1], axis[0], 0]])
    return np.eye(3) + math.sin(angle_rad) * K + (1 - math.cos(angle_rad)) * (K @ K)


def random_unit_axes(n: int, rng: np.random.Generator) -> np.ndarray:
    v = rng.standard_normal((n, 3))
    v /= np.linalg.norm(v, axis=1, keepdims=True)
    return v


from icenine.reconstructor import _get_voxel_vertices


def central_diff(angles, values):
    grad = np.empty_like(values)
    grad[1:-1] = (values[2:] - values[:-2]) / (angles[2:] - angles[:-2])
    grad[0]    = (values[1]  - values[0])   / (angles[1]  - angles[0])
    grad[-1]   = (values[-1] - values[-2])  / (angles[-1] - angles[-2])
    return grad


def setup():
    config_path = EXAMPLE_DIR / "ConfigFiles" / "Example2.Simulation.config"
    data_dir    = EXAMPLE_DIR / "ScatteringData_Python"

    os.chdir(EXAMPLE_DIR)  # configs use relative paths

    config = ConfigFile.from_file(str(config_path))
    config.out_file_basename = BASENAME
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
        mic_path = EXAMPLE_DIR / mic_path
    mic = MicFile.read(str(mic_path))
    print(f"  Loaded .mic: {mic_path.name}  ({len(mic.voxels)} voxels)")

    print("  Loading SparseImageStack ...")
    image_stack = SparseImageStack.from_image_directory(
        directory=str(data_dir),
        basename=BASENAME,
        ext="d",
        serial_length=5,
        n_omega=180,
        n_detectors=2,
        num_rows=2048,
        num_cols=2048,
        binary=True,
    )
    print(f"  SparseImageStack: {image_stack.memory_bytes / 1024:.1f} KB")

    print("  Loading ExperimentalData for hard cost ...")
    exp_data = ExperimentalData.from_image_directory(
        directory=str(data_dir),
        basename=BASENAME,
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

    # Build downsampled base stacks ONCE and share across omega_window variants
    print("  Building shared downsampled base stacks [4x, 8x] ...")
    shared_ds = MultiScaleImageStack.build_shared_base(
        image_stack, downsample_factors=[1, 4, 8]
    )

    diff_fns = {}
    for ow in OMEGA_WINDOWS:
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

    return mic, hard_fn, diff_fns


def select_voxels(mic, hard_fn, rng, n=N_VOXELS, threshold=0.1, max_scan=500):
    from icenine.geometry import matrix_to_euler
    candidates = []
    for i, voxel in enumerate(mic.voxels[:max_scan]):
        vertices = _get_voxel_vertices(voxel)
        info = hard_fn.evaluate(voxel.orientation, vertices, voxel.phase)
        if info.quality > threshold:
            candidates.append((voxel, i, info.quality))
    chosen_idx = rng.choice(len(candidates), size=n, replace=False)
    chosen = [candidates[i] for i in sorted(chosen_idx)]
    print(f"  Selected {n} voxels from {len(candidates)} candidates:")
    for voxel, vidx, q in chosen:
        euler = matrix_to_euler(voxel.orientation)
        print(f"    idx={vidx:5d}  phi1={euler[0]:7.2f}°  Phi={euler[1]:6.2f}°  "
              f"phi2={euler[2]:7.2f}°  hard_q={q:.4f}")
    return [(v, vi) for v, vi, _ in chosen]


def sweep(voxel_list, diff_fns, rng):
    """Returns dict: omega_window -> (n_voxels, n_angles) array of mean-over-axes quality."""
    axes = random_unit_axes(N_AXES, rng)
    n_angles = len(ANGLES_DEG)
    results = {ow: np.zeros((len(voxel_list), n_angles)) for ow in OMEGA_WINDOWS}

    n_total = len(voxel_list) * n_angles * N_AXES
    t0 = time.perf_counter()
    done = 0

    for vi, (voxel, vidx) in enumerate(voxel_list):
        vertices = _get_voxel_vertices(voxel)
        buf = {ow: np.zeros((n_angles, N_AXES)) for ow in OMEGA_WINDOWS}

        for ai, angle_deg in enumerate(ANGLES_DEG):
            angle_rad = angle_deg * math.pi / 180.0
            for xi, axis in enumerate(axes):
                R = rodrigues(axis, angle_rad)
                perturbed_t = torch.from_numpy(R @ voxel.orientation).float()
                with torch.no_grad():
                    for ow in OMEGA_WINDOWS:
                        info = diff_fns[ow].evaluate(
                            perturbed_t, vertices, phase_index=voxel.phase, scale=2
                        )
                        buf[ow][ai, xi] = info.quality.item()
            done += N_AXES
            elapsed = time.perf_counter() - t0
            eta = elapsed / done * (n_total - done) if done else 0
            vals = "  ".join(f"ω±{ow}={buf[ow][ai].mean():.4f}" for ow in OMEGA_WINDOWS)
            print(f"    vox {vi+1}/{len(voxel_list)}  angle={angle_deg:5.1f}°  {vals}"
                  f"  [ETA {eta/60:.0f}m]", flush=True)

        for ow in OMEGA_WINDOWS:
            results[ow][vi] = buf[ow].mean(axis=1)

    return results


def plot_quality(results, out_path):
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    ax_q, ax_rel = axes

    for ow in OMEGA_WINDOWS:
        arr = results[ow]  # (n_voxels, n_angles)
        mean = arr.mean(axis=0)
        std  = arr.std(axis=0)
        label = f"Diff s2 — 8× + ω±{ow}" if ow > 0 else "Diff s2 — 8× (no blend)"
        color = COLORS[str(ow)]
        ls    = LINESTYLES[str(ow)]

        # Spaghetti
        for vi in range(arr.shape[0]):
            ax_q.plot(ANGLES_DEG, arr[vi], color=color, ls=ls, lw=0.5, alpha=0.25)

        ax_q.plot(ANGLES_DEG, mean, color=color, ls=ls, lw=2.0, label=label)
        ax_q.fill_between(ANGLES_DEG, mean - std, mean + std, color=color, alpha=0.15)

        # Relative quality (normalised to Q at 0°)
        q0 = mean[0] if mean[0] > 1e-9 else 1.0
        ax_rel.plot(ANGLES_DEG, mean / q0, color=color, ls=ls, lw=2.0, label=label)

    for ax, ylabel, title in [
        (ax_q,   "Quality (0–1)",          "Quality vs. misorientation — scale 2 (8×)"),
        (ax_rel, "Q / Q(0°)",              "Relative quality [normalised to Q at 0°]"),
    ]:
        ax.axvline(0, color="gray", lw=0.8, ls="--")
        ax.set_xlabel("Misorientation angle (degrees)")
        ax.set_ylabel(ylabel)
        ax.set_title(f"ManyGrains ({N_VOXELS} voxels) — {title}")
        ax.legend(fontsize=9)
        ax.grid(True, alpha=0.3)
        ax.set_xlim(ANGLES_DEG[0], ANGLES_DEG[-1])

    ax_q.set_ylim(-0.02, 1.05)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    print(f"Saved: {out_path}")


def plot_gradient(results, out_path):
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    ax_abs, ax_rel = axes

    for ow in OMEGA_WINDOWS:
        arr = results[ow]  # (n_voxels, n_angles)
        label = f"Diff s2 — 8× + ω±{ow}" if ow > 0 else "Diff s2 — 8× (no blend)"
        color = COLORS[str(ow)]
        ls    = LINESTYLES[str(ow)]

        grads = np.array([central_diff(ANGLES_DEG, arr[vi]) for vi in range(arr.shape[0])])
        mean_g = grads.mean(axis=0)
        std_g  = grads.std(axis=0)

        # Spaghetti
        for vi in range(grads.shape[0]):
            ax_abs.plot(ANGLES_DEG, grads[vi], color=color, ls=ls, lw=0.5, alpha=0.25)

        ax_abs.plot(ANGLES_DEG, mean_g, color=color, ls=ls, lw=2.0, label=label)
        ax_abs.fill_between(ANGLES_DEG, mean_g - std_g, mean_g + std_g,
                            color=color, alpha=0.15)

        # Relative: normalise each voxel's gradient by its own Q(0°)
        q0 = arr[:, 0:1]  # (n_voxels, 1)
        rel = grads / np.maximum(q0, 1e-9)
        mean_r = rel.mean(axis=0)
        std_r  = rel.std(axis=0)

        for vi in range(rel.shape[0]):
            ax_rel.plot(ANGLES_DEG, rel[vi], color=color, ls=ls, lw=0.5, alpha=0.25)

        ax_rel.plot(ANGLES_DEG, mean_r, color=color, ls=ls, lw=2.0, label=label)
        ax_rel.fill_between(ANGLES_DEG, mean_r - std_r, mean_r + std_r,
                            color=color, alpha=0.15)

    for ax, ylabel, title in [
        (ax_abs, "dQ / dθ  (per degree)",       "Absolute gradient"),
        (ax_rel, "(dQ/dθ) / Q(0°)  (per degree)", "Relative gradient"),
    ]:
        ax.axhline(0, color="gray", lw=0.8, ls="--")
        ax.set_xlabel("Misorientation angle (degrees)")
        ax.set_ylabel(ylabel)
        ax.set_title(f"ManyGrains ({N_VOXELS} voxels) — {title}")
        ax.legend(fontsize=9)
        ax.grid(True, alpha=0.3)
        ax.set_xlim(ANGLES_DEG[0], ANGLES_DEG[-1])

    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    print(f"Saved: {out_path}")


if __name__ == "__main__":
    rng = np.random.default_rng(seed=SEED)

    print("Setting up ManyGrains example ...")
    mic, hard_fn, diff_fns = setup()

    print("\nSelecting voxels ...")
    voxel_list = select_voxels(mic, hard_fn, rng)

    print(f"\nSweeping {len(voxel_list)} voxels × {len(ANGLES_DEG)} angles × {N_AXES} axes "
          f"× {len(OMEGA_WINDOWS)} omega_windows ...")
    results = sweep(voxel_list, diff_fns, rng)

    plot_quality(results, benchmark_dir / "omega_window_comparison.png")
    plot_gradient(results, benchmark_dir / "omega_window_comparison_gradient.png")
