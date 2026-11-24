"""fine1 phase D: the two charts of the whole-layer result, and the pareto panel rows.

  uv run python scripts/autoresearch/fine1/fine1_plot.py

Reads ``benchmarks/autoresearch/fine1/score.json`` (``fine1_score.py``) and the launcher's
``whole_<arm>.time`` files. Writes in ``benchmarks/autoresearch/fine1/``:
``fine7_cost_accuracy.png`` (evaluations against wrong voxels, Wilson 95% bars, wall time and
workers in each label), ``fine7_wall_time.png`` (wall seconds and CPU seconds per arm, the worker
count marked, so one worker and twelve fitting workers do not look alike) and
``pareto_panel_fine7.json`` (the rows of the panel ``fine7_98280`` for the review page's
cost-against-accuracy chart; the page itself is not edited here).
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
    "CS": "stock C++, 1 worker (CS)",
    "CF": "corrected C++, 1 worker (CF)",
    "CSP": "stock C++, 12 workers (CSP)",
    "CFP": "corrected C++, 12 workers (CFP)",
    "J5S": "Python, latest recipe (J5S)",
    "PS": "Python, stock (PS)",
}
OFFS = {
    "CS": (12, 8, "left"),
    "CSP": (-12, 16, "right"),
    "CFP": (-12, 2, "right"),
    "CF": (-12, -14, "right"),
    "J5S": (10, 6, "left"),
    "PS": (10, 6, "left"),
}
COL = {
    "CS": "#c0392b",
    "CF": "#e67e22",
    "CSP": "#c0392b",
    "CFP": "#e67e22",
    "J5S": "#1f77b4",
    "PS": "#555555",
}


def cpu_seconds(arm: str) -> float:
    t: Dict[str, float] = {}
    for ln in (P.LOGS / f"whole_{arm}.time").read_text().splitlines():
        p = ln.split()
        if len(p) == 2:
            t[p[0]] = float(p[1])
    return t["user"] + t["sys"]


def main() -> int:
    s: Dict[str, Any] = json.loads((P.BENCH / "score.json").read_text())["arms"]
    fig, ax = plt.subplots(figsize=(9, 5.6))
    for arm, d in s.items():
        lo, hi = d["wilson95"]
        n = d["n_voxels"]
        y = d["n_wrong"]
        ax.errorbar(
            d["evals_all"],
            y,
            yerr=[[y - lo * n], [hi * n - y]],
            fmt="o" if d["workers"] == 1 else "s",
            color=COL[arm],
            mfc=COL[arm] if d["workers"] == 1 else "white",
            ms=8,
            capsize=3,
        )
        plural = "s" if d["workers"] > 1 else ""
        lab = f"{arm}: {d['wall_s'] / 3600:.1f} h wall, {d['workers']} worker{plural}"
        dx, dy, ha = OFFS[arm]
        ax.annotate(
            lab, (d["evals_all"], y), textcoords="offset points", xytext=(dx, dy), ha=ha, fontsize=8
        )
        cy = d["common_wrong"]
        clo, chi = d["common_wilson95"]
        ax.errorbar(
            d["evals_all"],
            cy,
            yerr=[[cy - clo * d["common_n"]], [chi * d["common_n"] - cy]],
            fmt="x",
            color=COL[arm],
            ms=6,
            capsize=2,
            alpha=0.7,
        )
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("cost evaluations (all, log scale)")
    ax.set_ylabel("wrong voxels (over 1 degree; log scale), of 98,280")
    ax.set_title(
        "Whole fine layer (98,280 voxels, clean images): evaluations against wrong voxels\n"
        "filled = 1 worker, open square = 12 workers; round = all voxels (no output = wrong),\n"
        "cross = the 85,676 voxels every arm has an output for; all runs contended",
        fontsize=8,
    )
    ax.grid(alpha=0.3, which="both")
    ax.set_xlim(2e7, 1.2e10)
    fig.tight_layout()
    fig.savefig(P.BENCH / "fine7_cost_accuracy.png", dpi=140)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(9, 4.6))
    arms = list(s)
    ypos = range(len(arms))
    wall = [s[a]["wall_s"] / 3600 for a in arms]
    cpu = [cpu_seconds(a) / 3600 for a in arms]
    ax.barh([y + 0.2 for y in ypos], wall, height=0.38, color=[COL[a] for a in arms], label="wall")
    ax.barh(
        [y - 0.2 for y in ypos],
        cpu,
        height=0.38,
        color=[COL[a] for a in arms],
        alpha=0.35,
        hatch="//",
        label="CPU (user + system, all cores)",
    )
    for y, a, w, c in zip(ypos, arms, wall, cpu):
        ax.text(
            w,
            y + 0.2,
            f" {w:.1f} h wall, {s[a]['workers']} fitting worker"
            f"{'s' if s[a]['workers'] > 1 else ''}",
            va="center",
            fontsize=8,
        )
        ax.text(c, y - 0.2, f" {c:.1f} h CPU", va="center", fontsize=8)
    ax.set_yticks(list(ypos))
    ax.set_yticklabels([NAMES[a] for a in arms], fontsize=8)
    ax.invert_yaxis()
    ax.set_xscale("log")
    ax.set_xlabel("hours (log scale); every run contended (a cpptools process on one core)")
    ax.set_title("Wall-clock time and CPU time per arm, whole fine layer", fontsize=9)
    ax.legend(fontsize=8, loc="center right")
    ax.set_xlim(1, 1500)
    fig.tight_layout()
    fig.savefig(P.BENCH / "fine7_wall_time.png", dpi=140)
    plt.close(fig)

    rows = [
        dict(
            method=NAMES[a],
            arm=a,
            evaluations=d["evals_all"],
            wrong=d["n_wrong"],
            wilson95_wrong=[d["wilson95"][0] * d["n_voxels"], d["wilson95"][1] * d["n_voxels"]],
            wrong_on_common_85676=d["common_wrong"],
            wall_s=d["wall_s"],
            workers=d["workers"],
            label="contended",
        )
        for a, d in s.items()
    ]
    (P.BENCH / "pareto_panel_fine7.json").write_text(
        json.dumps(
            dict(id="fine7_98280", title="Whole fine layer (98,280 voxels), clean", rows=rows),
            indent=1,
        )
        + "\n"
    )
    print("wrote charts and panel rows")
    return 0


if __name__ == "__main__":
    sys.exit(main())
