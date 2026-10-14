"""Task B: figures from the ray walks, summary.json and the 2-D slices.

uv run python scripts/eval_reduction/cost_edges_figures.py --out benchmarks/eval_reduction/B/figures
"""

import argparse
import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from scipy.signal import find_peaks  # noqa: E402

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import cost_edges as CE  # noqa: E402
import cost_edges_run as RUN  # noqa: E402

CAUSES = RUN.CAUSES
COL = {
    "frame": "tab:red", "pixel": "tab:blue", "detector": "tab:green", "eta_cut": "tab:orange",
    "no_frame": "tab:purple", "bragg_limit": "tab:brown", "subpixel": "tab:cyan", "rule": "k",
    "denominator": "tab:gray", "unattributed": "magenta",
}  # fmt: skip
REG = ["fin", "mc", "grid"]
REGL = {"fin": "< 0.05 deg", "mc": "0.05-0.5 deg", "grid": "0.5-5 deg"}


def fig_ray(out: Path) -> None:
    z = np.load(RUN.CACHE / "full" / "v15901_d0.npz")
    r = z["radii|dense"].astype(float)
    c = z["full|cpp|cost|dense"].astype(float)
    parts = z["full|cpp|parts|dense"]
    evn = z["full|cpp|evn|dense"]
    sil = z["full|cpp|silent|dense"]
    fig, axs = plt.subplots(3, 2, figsize=(11, 8), gridspec_kw=dict(height_ratios=[3, 1.2, 1.2]))
    for col, hi in enumerate((0.1, 0.5)):
        m = r <= hi
        ax = axs[0, col]
        ax.plot(r[m], c[m], color="0.3", lw=0.8)
        d = np.nonzero((np.diff(c) != 0) & m[1:])[0]
        for t in d:
            cause = CAUSES[int(np.abs(parts[t]).argmax())]
            ax.plot(r[t + 1], c[t + 1], "o", ms=3, color=COL[cause])
        ax.set_ylabel("cost (1 - quality)")
        ax.set_title(f"complete cost, C++ rule, voxel 15901, one ray: 0 to {hi} deg", fontsize=9)
        for cause in ("frame", "pixel", "no_frame", "detector"):
            ax.plot([], [], "o", color=COL[cause], label=f"step dominated by: {cause}")
        ax.legend(fontsize=7, loc="lower right")
        # events per bin
        nb = 50
        edges = np.linspace(0, hi, nb + 1)
        mid = 0.5 * (r[:-1] + r[1:])
        for row, (arr, ttl) in enumerate(
            ((evn - sil, "events that changed the cost"), (sil, "silent events (cost unchanged)")),
            1,
        ):
            bot = np.zeros(nb)
            for j, t in enumerate(CE.EVENT_TYPES):
                if not arr[:, j].any():
                    continue
                h, _ = np.histogram(mid, edges, weights=arr[:, j])
                axs[row, col].bar(
                    edges[:-1],
                    h,
                    width=np.diff(edges),
                    bottom=bot,
                    align="edge",
                    color=COL[t],
                    label=t,
                )
                bot += h
            axs[row, col].set_ylabel("events per bin")
            axs[row, col].set_title(ttl, fontsize=8)
            if col == 1 and row == 1:
                axs[row, col].legend(fontsize=7, ncol=3)
        axs[2, col].set_xlabel("rotation away from the true orientation (deg)")
    fig.tight_layout()
    fig.savefig(out / "cost_vs_angle_ray.png", dpi=140)
    plt.close(fig)


def fig_rates(S: dict, out: Path) -> None:
    fig, axs = plt.subplots(2, 2, figsize=(11, 7.5), sharey=True)
    for i, kind in enumerate(("full", "e0")):
        for j, setting in enumerate(("full", "screen")):
            ax = axs[i, j]
            s = S[kind][f"{setting}|cpp"]
            types = [t for t in CE.EVENT_TYPES if any(s[r]["events_per_deg"][t] for r in REG)]
            x = np.arange(len(types))
            for k, r in enumerate(REG):
                cc = np.array([s[r]["cost_changing_events_per_deg"][t] for t in types])
                al = np.array([s[r]["events_per_deg"][t] for t in types])
                ax.bar(
                    x + (k - 1) * 0.27,
                    al,
                    0.25,
                    color=["lightgray"] * len(x),
                    edgecolor="k",
                    lw=0.4,
                )
                ax.bar(
                    x + (k - 1) * 0.27,
                    cc,
                    0.25,
                    color=[f"C{k}"] * len(x),
                    label=f"{REGL[r]} (cost-changing)" if True else None,
                )
            ax.set_yscale("log")
            ax.set_xticks(x)
            ax.set_xticklabels(types, rotation=20, fontsize=8)
            ax.set_title(
                f"{'full-sample voxels' if kind == 'full' else 'E0 voxels'}, "
                f"{'complete cost' if setting == 'full' else 'level-0 screen'}, C++ rule",
                fontsize=9,
            )
            if j == 0:
                ax.set_ylabel("events per degree of rotation")
            if i == 0 and j == 0:
                ax.legend(fontsize=7)
    fig.suptitle(
        "Events per degree by type: grey = all events, colour = those that changed the cost",
        fontsize=10,
    )
    fig.tight_layout()
    fig.savefig(out / "events_per_degree.png", dpi=140)
    plt.close(fig)


def fig_share(S: dict, out: Path) -> None:
    fig, axs = plt.subplots(1, 2, figsize=(12, 4.8), sharey=True)
    for ax, kind in zip(axs, ("full", "e0")):
        labels, data = [], []
        for setting in ("full", "screen"):
            for rule in ("old", "cpp", "r2"):
                for r in REG:
                    x = S[kind][f"{setting}|{rule}"][r]
                    tot = sum(x["gross"].values())
                    labels.append(f"{setting[:4]} {rule}\n{REGL[r]}")
                    data.append([x["gross"][c] / tot if tot else 0 for c in CAUSES])
        data = np.array(data)
        bot = np.zeros(len(data))
        for j, c in enumerate(CAUSES):
            if data[:, j].any():
                ax.bar(range(len(data)), data[:, j], bottom=bot, color=COL[c], label=c)
                bot += data[:, j]
        ax.set_xticks(range(len(data)))
        ax.set_xticklabels(labels, rotation=90, fontsize=6)
        ax.set_title(
            ("full-sample voxels" if kind == "full" else "E0 voxels")
            + ": share of the gross change of the cost",
            fontsize=9,
        )
    axs[0].legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(out / "variation_share.png", dpi=140)
    plt.close(fig)


def fig_minima(out: Path) -> None:
    z = np.load(RUN.CACHE / "full" / "v15901_d0.npz")
    fig, axs = plt.subplots(2, 1, figsize=(10, 6.5))
    for ax, (setting, rule) in zip(axs, (("full", "cpp"), ("screen", "cpp"))):
        r = z["radii|dense"].astype(float)
        c = z[f"{setting}|{rule}|cost|dense"].astype(float)
        ax.plot(r, c, color="0.3", lw=0.7)
        for p, mk in ((0.01, "o"), (0.05, "s")):
            pk, _ = find_peaks(-c, prominence=p)
            ax.plot(r[pk], c[pk], mk, mfc="none", ms=7, label=f"local minimum, prominence >= {p}")
        ax.set_title(f"voxel 15901, one ray, {setting}, C++ rule", fontsize=9)
        ax.set_ylabel("cost")
        ax.legend(fontsize=7)
    axs[1].set_xlabel("rotation away from the true orientation (deg)")
    fig.tight_layout()
    fig.savefig(out / "minima_example.png", dpi=140)
    plt.close(fig)


def fig_slices(out: Path) -> None:
    fig, axs = plt.subplots(2, 3, figsize=(13, 8))
    z = np.load(RUN.CACHE.parent / "slices" / "v15901.npz")
    for i, scale in enumerate(("fine", "wide")):
        ax_deg = z[f"{scale}|axis_deg"]
        ext = (ax_deg[0], ax_deg[-1], ax_deg[0], ax_deg[-1])
        for j, rule in enumerate(("old", "cpp", "r2")):
            ax = axs[i, j]
            im = ax.imshow(
                z[f"{scale}|{rule}|cost"],
                origin="lower",
                extent=ext,
                cmap="viridis",
                vmin=0,
                vmax=float(z[f"{scale}|{rule}|cost"].max()),
            )
            ax.figure.colorbar(im, ax=ax, shrink=0.8)
            fr = z[f"{scale}|{rule}|frame_changed"]
            cs = ax.contour(
                ax_deg, ax_deg, fr, levels=[0.5, 5.5, 10.5, 20.5], colors="w", linewidths=0.6
            )
            ax.set_title(
                {"old": "lenient rule", "cpp": "C++ rule", "r2": "ray rule R2"}[rule], fontsize=9
            )
            ax.set_xlabel("deg along axis 1")
            if j == 0:
                ax.set_ylabel(f"deg along axis 2 ({scale} slice)")
    fig.suptitle(
        "Voxel 15901, complete cost (colour) on a plane through the truth: +-0.05 deg (top), "
        "+-0.5 deg (bottom).\nWhite contours: 1, 6, 11 and 21 peaks compared with a different "
        "frame than at the truth",
        fontsize=10,
    )
    fig.savefig(out / "slice_maps.png", dpi=130)
    plt.close(fig)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    out = Path(a.out).resolve()
    S = json.loads((out.parent / "summary.json").read_text())
    out.mkdir(parents=True, exist_ok=True)
    fig_ray(out)
    fig_rates(S, out)
    fig_share(S, out)
    fig_minima(out)
    fig_slices(out)
    print("figures in", out)


if __name__ == "__main__":
    main()
