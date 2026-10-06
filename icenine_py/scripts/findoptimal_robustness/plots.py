#!/usr/bin/env python3
"""The four plots of the FindOptimal robustness study (benchmarks/findoptimal_robustness/*.png)."""

import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import common as C  # noqa: E402

# Okabe-Ito (colour-blind safe), one hue per entity
BLUE, ORANGE, GREEN, VERM, PURPLE, SKY, GREY = (
    "#0072B2",
    "#E69F00",
    "#009E73",
    "#D55E00",
    "#CC79A7",
    "#56B4E9",
    "#7f7f7f",
)
VLABEL = {"clean": "clean", "all": "realistic"}


def style(ax):
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.grid(alpha=0.25, lw=0.5)


def plot_where_lost():
    S = json.loads((C.OUT_DIR / "e0_summary.json").read_text())
    classes = ["right", "S1", "S2@0", "S2@1", "S2@2", "S2b", "S2c", "S3"]
    cols = [GREY, VERM, ORANGE, "#F0C060", "#F5DFA0", PURPLE, SKY, BLUE]
    names = [
        "right (< 1 deg)",
        "S1 truth never a level-0 candidate",
        "S2 pruned at level 0",
        "S2 pruned at level 1",
        "S2 pruned at level 2",
        "S2b kept, not regenerated",
        "S2c lost in last level",
        "S3 in hand-off, FindOptimal wrong",
    ]
    fig, ax = plt.subplots(figsize=(7.5, 4))
    for i, var in enumerate(C.VARIANTS):
        s = S[var]
        n = s["n_runs"]
        det = s["where_lost_detail"]
        vals = {"right": n - s["n_wrong"]}
        for c in classes[1:]:
            vals[c] = sum(
                v for k, v in det.items() if k == c or (c == "S2b" and k.startswith("S2b")) and True
            )
        vals["S2b"] = sum(v for k, v in det.items() if k.startswith("S2b"))
        for lvl in (0, 1, 2):
            vals[f"S2@{lvl}"] = det.get(f"S2@{lvl}", 0)
        left = 0
        for c, col, nm in zip(classes, cols, names):
            w = vals.get(c, 0) / n
            ax.barh(
                i,
                w,
                left=left,
                color=col,
                edgecolor="white",
                lw=1.5,
                label=nm if i == 0 else None,
                height=0.55,
            )
            if w > 0.04:
                ax.text(left + w / 2, i, f"{100 * w:.0f}%", ha="center", va="center", fontsize=8)
            left += w
    ax.set_yticks(range(2))
    ax.set_yticklabels([f"{VLABEL[v]}\n(n={S[v]['n_runs']} runs)" for v in C.VARIANTS])
    ax.set_xlabel("fraction of runs")
    ax.set_xlim(0, 1)
    style(ax)
    ax.legend(fontsize=7, loc="upper center", bbox_to_anchor=(0.5, -0.18), ncol=2, frameon=False)
    ax.set_title(
        "Where the truth is lost (full reconstructions, 3 seeds x 200 voxels)", fontsize=10
    )
    fig.tight_layout()
    fig.savefig(C.OUT_DIR / "where_lost.png", dpi=150)


def plot_fixes():
    F = json.loads((C.OUT_DIR / "fixes_summary.json").read_text())
    extra = {}
    p = C.OUT_DIR / "e2_endtoend.json"
    if p.exists():
        extra = json.loads(p.read_text())
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.6), sharey=True)
    for ax, var in zip(axes, C.VARIANTS):
        rows = list(F.get(var, [])) + list(extra.get(var, []))
        for r in rows:
            nm = r["name"]
            col = (
                GREY
                if nm.startswith("E0")
                else (
                    GREEN
                    if nm.startswith("F1") and "b" not in nm.split()[0]
                    else (
                        BLUE
                        if nm.startswith("F2")
                        else (
                            ORANGE
                            if nm.startswith("F3")
                            else PURPLE if nm.startswith("F4") else VERM
                        )
                    )
                )
            )
            ax.errorbar(
                r["evals_mean"] / 1e3,
                r["rate"],
                yerr=[[r["rate"] - r["lo"]], [r["hi"] - r["rate"]]],
                fmt="o",
                color=col,
                ms=5,
                capsize=2,
                lw=1,
            )
            ax.annotate(
                nm.replace(" (seed 0)", "").replace("CSL relatives ", ""),
                (r["evals_mean"] / 1e3, r["rate"]),
                fontsize=6.5,
                xytext=(4, 3),
                textcoords="offset points",
            )
        ax.set_xlabel("cost evaluations per answer (thousands)")
        ax.set_title(f"{VLABEL[var]}", fontsize=10)
        style(ax)
    axes[0].set_ylabel("wrong rate (> 1 deg), 95% Wilson CI")
    fig.suptitle("Wrong rate against cost for each fix and the classifier", fontsize=10)
    fig.tight_layout()
    fig.savefig(C.OUT_DIR / "fixes_wrong_rate_vs_cost.png", dpi=150)


def plot_roc():
    r = np.load(C.OUT_DIR / "e2_roc.npz")
    fig, ax = plt.subplots(figsize=(4.8, 4.6))
    for k, col, lab in (
        ("cost", GREY, "local cost"),
        ("LR", BLUE, "logistic regression"),
        ("GBT", GREEN, "gradient-boosted trees"),
        ("GBT-nocost", ORANGE, "GBT without cost features"),
    ):
        if f"{k}_fpr" in r.files:
            ax.plot(r[f"{k}_fpr"], r[f"{k}_tpr"], color=col, lw=1.6, label=lab)
    ax.plot([0, 1], [0, 1], color="#bbbbbb", lw=0.8, ls="--")
    ax.set_xlabel("false positive rate (trap scored as basin)")
    ax.set_ylabel("true positive rate (basin)")
    ax.set_title("ROC, basin (<1 deg) vs trap, contested candidates,\nheld-out voxels", fontsize=9)
    ax.legend(fontsize=8, frameon=False, loc="lower right")
    style(ax)
    fig.tight_layout()
    fig.savefig(C.OUT_DIR / "e2_roc.png", dpi=150)


def plot_learning():
    R = json.loads((C.OUT_DIR / "e2_results.json").read_text())["learning_curve"]
    sizes = [s for s in R if s != "cost_baseline"]
    x = [int(s) for s in sizes]
    fig, axes = plt.subplots(1, 2, figsize=(9, 3.8))
    for ax, (suffix, ttl, base) in zip(
        axes,
        (
            ("auc", "ROC-AUC, basin vs trap (set A)", R["cost_baseline"]["auc_A"]),
            (
                "prune",
                "pruning recall (truth basin kept at levels 0-2)",
                R["cost_baseline"]["prune_recall"],
            ),
        ),
    ):
        for kind, col in (("GBT", GREEN), ("LR", BLUE)):
            m = np.array([R[s][f"{kind}_{suffix}"][0] for s in sizes])
            sd = np.array([R[s][f"{kind}_{suffix}"][1] for s in sizes])
            ax.errorbar(x, m, yerr=sd, color=col, marker="o", ms=4, capsize=2, lw=1.4, label=kind)
        ax.axhline(base, color=GREY, ls="--", lw=1, label="local cost (no training)")
        ax.set_xlabel("training voxels (both variants each)")
        ax.set_title(ttl, fontsize=9)
        ax.set_xscale("log")
        ax.set_xticks(x)
        ax.set_xticklabels(x)
        style(ax)
    axes[0].legend(fontsize=8, frameon=False)
    fig.suptitle("Learning curve (mean +- sd over folds and subsets; held-out voxels)", fontsize=10)
    fig.tight_layout()
    fig.savefig(C.OUT_DIR / "e2_learning_curve.png", dpi=150)


if __name__ == "__main__":
    which = sys.argv[1:] or ["where", "fixes", "roc", "learning"]
    if "where" in which:
        plot_where_lost()
    if "fixes" in which:
        plot_fixes()
    if "roc" in which:
        plot_roc()
    if "learning" in which:
        plot_learning()
