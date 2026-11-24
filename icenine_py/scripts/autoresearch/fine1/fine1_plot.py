"""fine1 phase D: the two charts of the whole-layer result, and the pareto panel rows.

  uv run python scripts/autoresearch/fine1/fine1_plot.py

Reads ``benchmarks/autoresearch/fine1/score.json`` (``fine1_score.py``). Writes in
``benchmarks/autoresearch/fine1/``: ``fine7_cost_accuracy.png`` (evaluations against wrong voxels;
the main points count the 85,676 voxels every arm has an output for, with Wilson 95% bars; the
light points count all 98,280 voxels, where the C++ arms' voxels without output count as wrong),
``fine7_wall_time.png`` (wall hours and CPU hours per arm on a linear axis, the worker count and
the 13 ranks on 12 cores marked) and ``pareto_panel_fine7.json`` (the rows of the panel
``fine7_98280`` for the review page; the page itself is not edited here).
"""

import json
import sys
from pathlib import Path
from typing import Any, Dict

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
import fine1_paths as P  # noqa: E402

NAMES = {
    "CS": "stock C++",
    "CF": "corrected C++",
    "CSP": "stock C++",
    "CFP": "corrected C++",
    "J5S": "Python, latest recipe",
    "PS": "Python, stock + 3 Python-only keys",
}
COL = {
    "CS": "#c0392b",
    "CF": "#e67e22",
    "CSP": "#c0392b",
    "CFP": "#e67e22",
    "J5S": "#1f77b4",
    "PS": "#555555",
}
# label offsets (points) on the cost chart: (dx, dy, ha)
OFFS = {
    "CS": (12, 10, "left"),
    "CSP": (-12, 12, "right"),
    "CFP": (-12, -12, "right"),
    "CF": (12, -14, "left"),
    "J5S": (10, 8, "left"),
    "PS": (10, 8, "left"),
}


def tag(arm: str, d: Dict[str, Any]) -> str:
    w = "1 worker" if d["workers"] == 1 else "12 workers, 13 ranks on 12 cores"
    return f"{arm} ({NAMES[arm]}, {w}), {d['wall_s'] / 3600:.1f} h wall"


def main() -> int:
    s: Dict[str, Any] = json.loads((P.BENCH / "score.json").read_text())["arms"]
    n_c = next(iter(s.values()))["common_n"]
    fig, ax = plt.subplots(figsize=(10, 6))
    for arm, d in s.items():
        par = d["workers"] > 1
        # light points: all voxels (no output counts as wrong)
        lo, hi = d["wilson95"]
        n, y = d["n_voxels"], d["n_wrong"]
        ax.errorbar(
            d["evals_all"],
            y,
            yerr=[[y - lo * n], [hi * n - y]],
            fmt="s" if par else "o",
            color=COL[arm],
            mfc="white",
            alpha=0.45,
            ms=7,
            capsize=2,
        )
        # main points: the common voxels
        cy = d["common_wrong"]
        clo, chi = d["common_wilson95"]
        ax.errorbar(
            d["evals_all"],
            cy,
            yerr=[[cy - clo * n_c], [chi * n_c - cy]],
            fmt="s" if par else "o",
            color=COL[arm],
            mfc="white" if par else COL[arm],
            ms=9,
            capsize=3,
            mew=2,
        )
        dx, dy, ha = OFFS[arm]
        ax.annotate(
            tag(arm, d),
            (d["evals_all"], cy),
            textcoords="offset points",
            xytext=(dx, dy),
            ha=ha,
            fontsize=8,
        )
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("cost evaluations (all, log scale)")
    ax.set_ylabel(f"wrong voxels (over 1 degree, log scale), of {n_c:,} common voxels")
    ax.set_title(
        f"Whole fine layer, clean images. MAIN points: wrong voxels among the {n_c:,} voxels "
        "every arm has an output for.\nLight points: all 98,280 voxels, where the C++ arms' "
        "~10,000 voxels without output (our SampleRadius 0.50 set-up) count as wrong.\n"
        "Round = 1 worker; open square = 12 workers (13 MPI ranks on 12 cores, oversubscribed). "
        "All runs contended.",
        fontsize=8,
    )
    ax.grid(alpha=0.3, which="both")
    ax.set_xlim(2e7, 1.2e10)
    fig.tight_layout()
    fig.savefig(P.BENCH / "fine7_cost_accuracy.png", dpi=140)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(10, 4.8))
    arms = list(s)
    ypos = list(range(len(arms)))
    wall = [s[a]["wall_s"] / 3600 for a in arms]
    cpu = json.loads((P.BENCH / "score.json").read_text())["derived"]
    cpu_h = [cpu[a]["cpu_h"] for a in arms]
    ax.barh(
        [y + 0.2 for y in ypos],
        wall,
        height=0.38,
        color=[COL[a] for a in arms],
        hatch=["" if s[a]["workers"] == 1 else "xx" for a in arms],
        edgecolor="black",
        label="wall hours",
    )
    ax.barh(
        [y - 0.2 for y in ypos],
        cpu_h,
        height=0.38,
        color=[COL[a] for a in arms],
        alpha=0.35,
        label="CPU hours (user + system, all cores)",
    )
    for y, a, w, c in zip(ypos, arms, wall, cpu_h):
        extra = (
            "" if s[a]["workers"] == 1 else " (12 workers, 13 ranks on 12 cores: oversubscribed)"
        )
        ax.text(w, y + 0.2, f" {w:.1f} h wall{extra}", va="center", fontsize=8)
        ax.text(c, y - 0.2, f" {c:.1f} h CPU", va="center", fontsize=8)
    ax.set_yticks(ypos)
    ax.set_yticklabels([f"{a}: {NAMES[a]}" for a in arms], fontsize=8)
    ax.invert_yaxis()
    ax.set_xlabel("hours (linear axis); every run contended (a cpptools process on one core)")
    ax.set_title(
        "Wall-clock and CPU time per arm, whole fine layer (hatched = 12 fitting workers)",
        fontsize=9,
    )
    ax.set_xlim(0, 80)
    ax.legend(fontsize=8, loc="center right")
    fig.tight_layout()
    fig.savefig(P.BENCH / "fine7_wall_time.png", dpi=140)
    plt.close(fig)

    rows = [
        dict(
            method=f"{NAMES[a]} ({a})",
            arm=a,
            evaluations=d["evals_all"],
            wrong_on_common_85676=d["common_wrong"],
            wilson95_common=[d["common_wilson95"][0] * n_c, d["common_wilson95"][1] * n_c],
            wrong_all_98280=d["n_wrong"],
            wall_s=d["wall_s"],
            workers=d["workers"],
            label="contended" + ("; 13 ranks on 12 cores" if d["workers"] > 1 else ""),
        )
        for a, d in s.items()
    ]
    panel = dict(id="fine7_98280", title="Whole fine layer (98,280 voxels), clean", rows=rows)
    (P.BENCH / "pareto_panel_fine7.json").write_text(json.dumps(panel, indent=1) + "\n")
    print("wrote charts and panel rows")
    return 0


if __name__ == "__main__":
    sys.exit(main())
