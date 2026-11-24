"""tex1: the owner's plots of the whole-layer result (run after ``tex1_score.py``).

  uv run python scripts/autoresearch/tex1/tex1_plot.py [--arms CSP,RCP]

Reads ``benchmarks/autoresearch/tex1/score.json`` and the maps. Writes (PNG and SVG) in
``benchmarks/autoresearch/tex1/``, plots first:

  1. ``tex24_<arm>_maps``: the truth orientation map and the arm's orientation map (top row) and
     the misorientation map against each voxel's own truth, scale 0 to 2 degrees (over 2
     saturated), with the truth grain boundaries drawn on top (bottom row). ``plot_mic.py`` with
     its two opt-in flags ``--boundaries`` and ``--mis-max``. The map SVG (98,280 triangles) goes
     to the cache only.
  2. ``tex24_<arm>_distribution``: misorientation distribution, log x axis 0.001 to 63 degrees,
     the 1-degree line marked, all voxels and the four classes, the trivial answer dashed (errors
     under 0.001 sit in the first bin; voxels with no output are left out and counted in the
     title).
  3. ``tex24_cost_accuracy`` (the regroup page's Chart A: evaluations, log, against wrong voxels,
     the voxels both arms have an output for, Wilson 95%) and ``tex24_wall_time`` (Chart B: wall
     seconds, linear, hatched = 12 workers), the two ``tex24`` points; the ``fine7`` points are
     grey and labelled as another sample. Rows in ``charts_tex24.json`` (the page's format).
"""

import argparse
import gzip
import json
import shutil
import sys
import tempfile
from pathlib import Path
from typing import Any, Dict, List, Optional

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1]))  # scripts/ (plot_mic)
import plot_mic  # noqa: E402
import tex1_lib as F  # noqa: E402
import tex1_paths as P  # noqa: E402

NAMES = {"CSP": "stock C++, 12 workers", "RCP": "Python, recipe of 24 Nov"}
COL = {"CSP": "#c0392b", "RCP": "#1f77b4"}
GREY = "#9a9a9a"
FINE7_PANEL = P.ICE / "benchmarks" / "autoresearch" / "fine1" / "pareto_panel_fine7.json"
FINE7_ARMS = {"CSP": "stock C++ (CSP)", "J5S": "Python recipe (J5S)"}
BINS = np.logspace(-3, np.log10(63.0), 61)
MIS_MAX_DEG = 2.0


def histogram_fractions(err: np.ndarray, bins: np.ndarray = BINS) -> np.ndarray:
    """Share of the voxels per bin; errors under the first edge sit in the first bin and errors over
    the last in the last (nothing is dropped). Empty input gives zeros."""
    e = np.asarray(err, dtype=float)
    if len(e) == 0:
        return np.zeros(len(bins) - 1)
    c, _ = np.histogram(np.clip(e, bins[0], bins[-1]), bins=bins)
    return c / len(e)


def distribution_figure(
    err: np.ndarray,
    has: np.ndarray,
    triv: np.ndarray,
    vox_class: np.ndarray,
    arm: str,
    out_stem: Path,
) -> None:
    """Five panels (all voxels, sharp, noise, gradient, both): histogram of the arm's error over
    the voxels with output (solid) and of the trivial answer over the same voxels (dashed)."""
    fig, axes = plt.subplots(1, 5, figsize=(19, 3.8), sharey=True)
    centres = np.sqrt(BINS[:-1] * BINS[1:])
    n_miss = int((~has).sum())
    for ax, name in zip(axes, ("all",) + F.CLASS_NAMES):
        sel = has if name == "all" else has & (vox_class == F.CLASS_NAMES.index(name))
        ax.plot(centres, histogram_fractions(err[sel]), color=COL[arm], lw=1.6, label=arm)
        ax.plot(
            centres,
            histogram_fractions(triv[sel]),
            color="k",
            lw=1.2,
            ls="--",
            label="trivial answer (R_g)",
        )
        ax.axvline(1.0, color="#555555", lw=0.8, ls=":")
        ax.set_xscale("log")
        ax.set_xlim(BINS[0], BINS[-1])
        ax.set_title(f"{name}: {int(sel.sum()):,} voxels", fontsize=9)
        ax.set_xlabel("misorientation to the voxel's own truth (deg, log)")
        ax.grid(alpha=0.3, which="both")
    axes[0].set_ylabel("share of voxels per bin")
    axes[0].legend(fontsize=8)
    fig.suptitle(
        f"{arm} ({NAMES[arm]}): distribution of the error; dotted line = 1 degree; "
        f"{n_miss:,} voxels with no output left out",
        fontsize=10,
    )
    fig.tight_layout()
    for ext in ("png", "svg"):
        fig.savefig(f"{out_stem}.{ext}", dpi=140)
    plt.close(fig)


def fine7_points(path: Path = FINE7_PANEL) -> List[Dict[str, Any]]:
    """The fine7 CSP and J5S rows of fine1's panel (common-voxel counts): the grey points."""
    rows = json.loads(path.read_text())["rows"]
    return [r for r in rows if r["arm"] in FINE7_ARMS]


def chart_rows(score: Dict[str, Any]) -> Dict[str, Any]:
    """``charts_tex24.json``: the rows in the format of ``fine1``'s ``pareto_panel_fine7.json``
    (counts over the voxels both arms have an output for, and over all voxels)."""
    rows = []
    for arm, d in score["arms"].items():
        c = d["common_classes"]["all"]
        a = d["classes"]["all"]
        n_c = d["common_n"]
        rows.append(
            dict(
                method=f"{NAMES[arm]} ({arm})",
                arm=arm,
                evaluations=d["evals_all"],
                wrong_on_common=c["n_wrong"],
                common_n=n_c,
                wilson95_common=[c["wilson95"][0] * n_c, c["wilson95"][1] * n_c],
                wrong_all=a["n_wrong"],
                n_all=a["n"],
                wall_s=d["wall_s"],
                workers=d["workers"],
                label="contended"
                + ("; 13 ranks on 12 cores" if d["workers"] > 1 else "; single worker"),
            )
        )
    return dict(
        id="tex24_98280", title="Textured fine layer tex24 (98,280 voxels), clean", rows=rows
    )


def cost_accuracy_figure(panel: Dict[str, Any], grey: List[Dict[str, Any]], out_stem: Path) -> None:
    """Chart A: evaluations (log) against wrong voxels (log); tex24 points in colour, fine7 grey."""
    fig, ax = plt.subplots(figsize=(9, 5.6))
    for r in grey:
        ax.plot(r["evaluations"], r["wrong_on_common_85676"], "o", color=GREY, ms=9, mfc="white")
        ax.annotate(
            f"fine7: {FINE7_ARMS[r['arm']]} (another sample)",
            (r["evaluations"], r["wrong_on_common_85676"]),
            textcoords="offset points",
            xytext=(8, -12),
            fontsize=7,
            color=GREY,
        )
    for r in panel["rows"]:
        arm = r["arm"]
        y = max(r["wrong_on_common"], 0.5)
        lo, hi = r["wilson95_common"]
        ax.errorbar(
            r["evaluations"],
            y,
            yerr=[[max(y - lo, 0.0)], [max(hi - y, 0.0)]],
            fmt="s" if r["workers"] > 1 else "o",
            color=COL[arm],
            mfc="white" if r["workers"] > 1 else COL[arm],
            ms=10,
            capsize=3,
            mew=2,
        )
        ax.annotate(
            f"tex24 {arm}: {r['wrong_on_common']:,} of {r['common_n']:,} wrong, "
            f"{r['wall_s'] / 3600:.1f} h wall",
            (r["evaluations"], y),
            textcoords="offset points",
            xytext=(10, 10),
            fontsize=8,
        )
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("cost evaluations (all, log scale)")
    ax.set_ylabel("wrong voxels (over 1 degree, log scale)")
    ax.set_title(
        "Chart A, textured fine layer tex24 (coloured; Wilson 95%) and fine7 (grey, another "
        "sample).\nRound = 1 worker; open square = 12 workers (13 MPI ranks on 12 cores). "
        "Voxels both arms have an output for.",
        fontsize=8,
    )
    ax.grid(alpha=0.3, which="both")
    fig.tight_layout()
    for ext in ("png", "svg"):
        fig.savefig(f"{out_stem}.{ext}", dpi=140)
    plt.close(fig)


def wall_time_figure(panel: Dict[str, Any], grey: List[Dict[str, Any]], out_stem: Path) -> None:
    """Chart B: wall seconds (linear); hatched = 12 workers; fine7 bars grey."""
    fig, ax = plt.subplots(figsize=(9, 3.8))
    labels, vals, cols, hatch, edge = [], [], [], [], []
    for r in grey:
        labels.append(f"fine7 {FINE7_ARMS[r['arm']]} (another sample)")
        vals.append(r["wall_s"])
        cols.append(GREY)
        hatch.append("xx" if r["workers"] > 1 else "")
        edge.append("#666666")
    for r in panel["rows"]:
        labels.append(f"tex24 {r['arm']}: {NAMES[r['arm']]}")
        vals.append(r["wall_s"])
        cols.append(COL[r["arm"]])
        hatch.append("xx" if r["workers"] > 1 else "")
        edge.append("black")
    y = np.arange(len(vals))
    for yi, v, c, h, e in zip(y, vals, cols, hatch, edge):
        ax.barh(yi, v, color=c, hatch=h, edgecolor=e, height=0.6)
        ax.text(v, yi, f" {v:,.0f} s ({v / 3600:.2f} h)", va="center", fontsize=8)
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=8)
    ax.invert_yaxis()
    ax.set_xlabel("wall seconds (linear); hatched = 12 fitting workers; every run contended")
    ax.set_xlim(0, max(vals) * 1.3)
    ax.set_title("Chart B, wall-clock time", fontsize=9)
    fig.tight_layout()
    for ext in ("png", "svg"):
        fig.savefig(f"{out_stem}.{ext}", dpi=140)
    plt.close(fig)


def map_figure(
    arm: str,
    mic_path: str,
    out_png: Path,
    truth: Optional[str] = None,
    grainmap: Optional[str] = None,
) -> List[Dict[str, Any]]:
    """Truth and arm orientation maps (top) and the misorientation map with the truth grain
    boundaries, 0 to 2 degrees (bottom): ``plot_mic.plot_maps`` with its two opt-in arguments."""
    return plot_mic.plot_maps(
        truth or str(P.TRUTH),
        [mic_path],
        [f"{arm}: {NAMES[arm]}"],
        str(out_png),
        ref_label="tex24 truth",
        title=f"tex24, {arm}: orientation maps and misorientation to the truth "
        f"(0 to {MIS_MAX_DEG:g} deg, grain boundaries in black)",
        boundaries=grainmap or str(P.GRAINMAP),
        mis_max=MIS_MAX_DEG,
    )


def decompressed(arm: str, td: Path) -> str:
    run = P.RUNS / "whole" / arm
    if arm == "CSP":
        return str(run / "tex1_cpp.mic.LBFS")
    out = td / "recon.mic"
    with gzip.open(run / "recon.mic.gz", "rb") as fi, open(out, "wb") as fo:
        shutil.copyfileobj(fi, fo)
    return str(out)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--arms", default="CSP,RCP")
    a = ap.parse_args()
    arms = a.arms.split(",")
    score = json.loads((P.BENCH / "score.json").read_text())
    err_dir = P.SCORE
    triv = np.load(err_dir / "trivial_err.npy")
    vox_class = np.load(err_dir / "voxel_class.npy")
    for arm in arms:
        err = np.load(err_dir / f"{arm}_err.npy")
        has = np.isfinite(err)
        with tempfile.TemporaryDirectory() as td:
            mic = decompressed(arm, Path(td))
            st = map_figure(arm, mic, P.BENCH / f"tex24_{arm}_maps.png")
            (P.BENCH / f"tex24_{arm}_maps_stats.json").write_text(json.dumps(st, indent=1) + "\n")
        distribution_figure(err, has, triv, vox_class, arm, P.BENCH / f"tex24_{arm}_distribution")
        print("plotted", arm, flush=True)
    panel = chart_rows(score)
    (P.BENCH / "charts_tex24.json").write_text(json.dumps(panel, indent=1) + "\n")
    cost_accuracy_figure(panel, fine7_points(), P.BENCH / "tex24_cost_accuracy")
    wall_time_figure(panel, fine7_points(), P.BENCH / "tex24_wall_time")
    print("wrote charts and rows")
    return 0


if __name__ == "__main__":
    sys.exit(main())
