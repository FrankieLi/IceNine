"""
Finite-difference gradient of quality vs. misorientation angle.

Reads the existing CSVs and plots dQ/d(angle) for each cost function variant.
The gradient is computed as the central finite difference of the mean quality
curve:  dQ/dθ ≈ (Q[i+1] - Q[i-1]) / (2 * Δθ),  with one-sided differences
at the endpoints.

Relative gradient (normalized by Q at θ=0) is also shown so curves with
different absolute quality levels can be compared.
"""

from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

benchmark_dir = Path(__file__).parent

COLS = [
    "hard_quality",
    "diff_quality_s0",
    "diff_quality_s1",
    "diff_quality_s2",
    "diff_quality_s2_oblend",
]

STYLES = [
    ("hard_quality", "Hard (binary)", "k", "-", 2.0),
    ("diff_quality_s0", "Diff s0 — 1×", "C0", "-", 1.6),
    ("diff_quality_s1", "Diff s1 — 4×", "C1", "--", 1.6),
    ("diff_quality_s2", "Diff s2 — 8×", "C3", ":", 1.6),
    ("diff_quality_s2_oblend", "Diff s2 — 8× + ω±1", "C2", "-.", 1.6),
]


def load_mean_curves(csv_path: Path, group_col: str) -> tuple:
    """Return (angles, mean_per_col) where mean_per_col is a dict col→array."""
    df = pd.read_csv(csv_path)
    grouped = df.groupby("angle_deg")[COLS].mean()
    angles = grouped.index.values
    means = {col: grouped[col].values for col in COLS}
    return angles, means


def central_diff(angles: np.ndarray, values: np.ndarray) -> np.ndarray:
    """Central finite difference dV/dθ, one-sided at boundaries."""
    grad = np.empty_like(values)
    grad[1:-1] = (values[2:] - values[:-2]) / (angles[2:] - angles[:-2])
    grad[0] = (values[1] - values[0]) / (angles[1] - angles[0])
    grad[-1] = (values[-1] - values[-2]) / (angles[-1] - angles[-2])
    return grad


def plot_gradients(csv_path: Path, out_path: Path, title_prefix: str, group_col: str = "voxel_idx"):
    angles, means = load_mean_curves(csv_path, group_col)

    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    ax_abs, ax_rel = axes

    q0 = {col: means[col][0] for col in COLS}  # quality at θ=0

    for col, label, color, ls, lw in STYLES:
        q = means[col]
        grad = central_diff(angles, q)

        # Absolute gradient (dQ/dθ, units: quality per degree)
        ax_abs.plot(angles, grad, color=color, ls=ls, lw=lw, label=label)

        # Relative gradient: (dQ/dθ) / Q(0)  — normalised by peak quality
        denom = q0[col] if q0[col] > 1e-9 else 1.0
        ax_rel.plot(angles, grad / denom, color=color, ls=ls, lw=lw, label=label)

    for ax, ylabel, subtitle in [
        (ax_abs, "dQ / dθ  (quality per degree)", "Absolute gradient"),
        (ax_rel, "(dQ/dθ) / Q(0)  (per degree)", "Relative gradient  [normalised by Q at 0°]"),
    ]:
        ax.axhline(0, color="gray", lw=0.8, ls="--")
        ax.set_xlabel("Misorientation angle (degrees)")
        ax.set_ylabel(ylabel)
        ax.set_title(f"{title_prefix} — {subtitle}")
        ax.legend(fontsize=9)
        ax.grid(True, alpha=0.3)
        ax.set_xlim(angles[0], angles[-1])

    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    print(f"Saved: {out_path}")


def plot_per_voxel_gradients(csv_path: Path, out_path: Path, title_prefix: str):
    """Spaghetti plot of per-voxel gradients for ManyGrains."""
    df = pd.read_csv(csv_path)
    voxels = df["voxel_idx"].unique()
    angles = np.sort(df["angle_deg"].unique())

    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    ax_abs, ax_rel = axes

    for col, label, color, ls, lw in STYLES:
        per_voxel = []
        for vid in voxels:
            sub = df[df["voxel_idx"] == vid].sort_values("angle_deg")
            q = sub[col].values
            g = central_diff(angles, q)
            per_voxel.append(g)
        per_voxel = np.array(per_voxel)  # (n_voxels, n_angles)
        mean_g = per_voxel.mean(axis=0)
        std_g = per_voxel.std(axis=0)

        # Spaghetti
        for vi in range(len(voxels)):
            ax_abs.plot(angles, per_voxel[vi], color=color, ls=ls, lw=0.4, alpha=0.2)
            ax_rel.plot(
                angles,
                per_voxel[vi] / (per_voxel[vi][0] if abs(per_voxel[vi][0]) > 1e-9 else 1),
                color=color,
                ls=ls,
                lw=0.4,
                alpha=0.2,
            )

        ax_abs.plot(angles, mean_g, color=color, ls=ls, lw=lw, label=label)
        ax_abs.fill_between(angles, mean_g - std_g, mean_g + std_g, color=color, alpha=0.15)

        # Relative (normalise each voxel by its own Q at θ=0, then take mean)
        q0_per_voxel = np.array(
            [df[df["voxel_idx"] == vid].sort_values("angle_deg")[col].iloc[0] for vid in voxels]
        )
        rel = per_voxel / np.maximum(q0_per_voxel[:, None], 1e-9)
        mean_rel = rel.mean(axis=0)
        std_rel = rel.std(axis=0)
        ax_rel.plot(angles, mean_rel, color=color, ls=ls, lw=lw, label=label)
        ax_rel.fill_between(angles, mean_rel - std_rel, mean_rel + std_rel, color=color, alpha=0.15)

    for ax, ylabel, subtitle in [
        (ax_abs, "dQ / dθ  (quality per degree)", "Absolute gradient"),
        (ax_rel, "(dQ/dθ) / Q(0)  (per degree)", "Relative gradient  [normalised by voxel Q(0)]"),
    ]:
        ax.axhline(0, color="gray", lw=0.8, ls="--")
        ax.set_xlabel("Misorientation angle (degrees)")
        ax.set_ylabel(ylabel)
        ax.set_title(f"{title_prefix} — {subtitle}")
        ax.legend(fontsize=9)
        ax.grid(True, alpha=0.3)
        ax.set_xlim(angles[0], angles[-1])

    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    print(f"Saved: {out_path}")


if __name__ == "__main__":
    # ThreeVoxels: single voxel, average over axes
    plot_gradients(
        csv_path=benchmark_dir / "misorientation_threevoxels.csv",
        out_path=benchmark_dir / "gradient_threevoxels.png",
        title_prefix="ThreeVoxels (voxel 0)",
        group_col="axis_idx",
    )

    # ManyGrains: spaghetti per voxel + mean ± std
    plot_per_voxel_gradients(
        csv_path=benchmark_dir / "misorientation_manygrains.csv",
        out_path=benchmark_dir / "gradient_manygrains.png",
        title_prefix="ManyGrains (20 voxels)",
    )
