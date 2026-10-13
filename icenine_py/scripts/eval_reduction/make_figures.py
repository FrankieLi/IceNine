"""Task A figures (PNG; no data of their own: they read summary.json, attribution.json and
screen_timing.json).

uv run python scripts/eval_reduction/make_figures.py --a benchmarks/eval_reduction/A
"""

import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

STAGE_COLORS = {
    "screen": "#4c78a8",
    "recost": "#72b7b2",
    "rerecost": "#e45756",
    "quick": "#f58518",
    "find": "#54a24b",
    "variance": "#b279a2",
}


def fig_stages(root: Path, out: Path) -> None:
    A = json.loads((root / "step2" / "attribution.json").read_text())
    arms = [
        ("cpp", "Python, C++ rule"),
        ("cpp_edge", "Python, C++ rule + edge rules"),
        ("cpp_edge_interp", "Python, + omega interpolation"),
    ]
    cpp = {r["stage"]: r["cpp"] for r in A["stages_cpp"]["rows"]}
    bars = [("C++", cpp)] + [
        (lab, {r["stage"]: r["python"] for r in A[f"stages_{k}"]["rows"]}) for k, lab in arms
    ]
    fig, ax = plt.subplots(figsize=(8.5, 4.2))
    y = np.arange(len(bars))[::-1]
    for yi, (lab, d) in zip(y, bars):
        left = 0.0
        for s, c in STAGE_COLORS.items():
            v = d.get(s, 0.0)
            ax.barh(yi, v, left=left, color=c, label=s if yi == y[0] else None)
            left += v
        ax.text(left + 2000, yi, f"{left:,.0f}", va="center", fontsize=9)
    ax.set_yticks(y)
    ax.set_yticklabels([b[0] for b in bars])
    ax.set_xlabel("cost evaluations per no-start search (mean of six voxels)")
    ax.set_xlim(0, 290000)
    ax.legend(ncol=6, fontsize=8, loc="upper center", bbox_to_anchor=(0.5, -0.18))
    ax.set_title(
        "Where the evaluations go (red: C++'s second re-cost of each candidate)", fontsize=10
    )
    fig.tight_layout()
    fig.savefig(out / "stages_python_vs_cpp.png", dpi=130)
    plt.close(fig)


def fig_speed(root: Path, out: Path) -> None:
    r = json.loads((root / "step3" / "screen_timing.json").read_text())
    arms = [
        ("unbatched", "exact,\none per call\n(today)"),
        ("batched", "exact,\n32 FZ points\nper call"),
        ("interp1", "interpolated,\none FZ point\nper call"),
        ("interp32", "interpolated,\n32 FZ points\nper call"),
    ]
    fig, ax = plt.subplots(figsize=(7.5, 4.2))
    x = np.arange(len(arms))
    bottom = np.zeros(len(arms))
    for key, lab, col in (
        ("us_solve", "exact omega solve", "#4c78a8"),
        ("us_stage_d", "C overlap stage", "#54a24b"),
        ("us_other", "torch geometry and Python glue", "#bab0ac"),
    ):
        v = np.array([r[a][key] for a, _ in arms])
        ax.bar(x, v, bottom=bottom, label=lab, color=col)
        bottom += v
    for xi, tot, (a_, _) in zip(x, bottom, arms):
        untimed = r[a_]["us_per_orientation"]
        ax.text(xi, tot + 5, f"{tot:.0f} (untimed {untimed:.0f})", ha="center", fontsize=8)
    ax.set_xticks(x)
    ax.set_xticklabels([lab for _, lab in arms], fontsize=8)
    ax.set_ylabel("microseconds per orientation (single worker; bars from the timed run)")
    ax.legend(fontsize=8)
    ax.set_title("Level-0 screen, voxel 15901: the solve is 16% of today's time")
    fig.tight_layout()
    fig.savefig(out / "screen_time_split.png", dpi=130)
    plt.close(fig)


def fig_parity(root: Path, out: Path) -> None:
    S = json.loads((root / "step0" / "summary.json").read_text())
    labels = [
        "screen\n(peak counts)",
        "complete cost\n(confidence)",
        "complete cost\n(pixel counts)",
    ]
    default = [
        S["exact"]["EXS"]["conf_equal"]["pct"],
        S["exact"]["EXF"]["conf_equal"]["pct"],
        S["exact"]["EXF"]["counts_equal"]["pct"],
    ]
    edge = [
        S["exact"]["EXS_edge"]["conf_equal"]["pct"],
        S["exact"]["EXF_edge"]["conf_equal"]["pct"],
        S["exact"]["EXF_edge"]["counts_equal"]["pct"],
    ]
    fig, ax = plt.subplots(figsize=(7, 4))
    x = np.arange(3)
    ax.bar(x - 0.18, default, 0.36, label="Python default", color="#bab0ac")
    ax.bar(x + 0.18, edge, 0.36, label="with CppEdgeRules", color="#4c78a8")
    for xi, d, e in zip(x, default, edge):
        ax.text(xi - 0.18, d + 1, f"{d:.2f}", ha="center", fontsize=8)
        ax.text(xi + 0.18, e + 1, f"{e:.2f}", ha="center", fontsize=8)
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.set_ylim(0, 128)
    ax.set_ylabel(f"% of {S['n_orientations']:,} fixed orientations identical to C++")
    ax.legend(loc="upper left", fontsize=9)
    ax.set_title("Fixed-orientation parity with C++")
    fig.tight_layout()
    fig.savefig(out / "parity_with_cpp.png", dpi=130)
    plt.close(fig)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--a", required=True)
    a = ap.parse_args()
    root = Path(a.a).resolve()
    out = root / "figures"
    out.mkdir(exist_ok=True)
    fig_stages(root, out)
    fig_speed(root, out)
    fig_parity(root, out)
    print("figures in", out)


if __name__ == "__main__":
    main()
