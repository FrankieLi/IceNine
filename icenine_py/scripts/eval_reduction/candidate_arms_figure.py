"""Task C diagram: evaluations and wrong counts of every arm, side by side.

  uv run python scripts/eval_reduction/candidate_arms_figure.py

Reads ``benchmarks/eval_reduction/C/step3/summary.json`` (candidate_summary.py).
Left: evaluations of
the six full-sample voxels, as a share of the baseline's. Right: wrong voxels of the 200 E0 voxels
(bars; line = Wilson interval) with, per arm, the voxels that became wrong (broken) and that
became right (fixed) against the baseline. Writes ``figures/arms_evals_and_wrong.png``.
"""

import json
from pathlib import Path

import numpy as np

ICE = Path(__file__).resolve().parents[2]
BASE = ICE / "benchmarks" / "eval_reduction" / "C"


def main() -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    summ = json.loads((BASE / "step3" / "summary.json").read_text())
    arms = [a for a in summ if not a.startswith("_")]
    ev6 = [100 * summ[a]["six"]["evals_ratio"] for a in arms]
    wrong6 = [summ[a]["six"]["wrong"] for a in arms]
    e0 = [summ[a].get("e0") for a in arms]
    wr = [x["wrong"] if x else np.nan for x in e0]
    lo = [x["wilson_lo"] * x["n"] if x else np.nan for x in e0]
    hi = [x["wilson_hi"] * x["n"] if x else np.nan for x in e0]
    br = [x["broken"] if x else 0 for x in e0]
    fx = [x["fixed"] if x else 0 for x in e0]
    y = np.arange(len(arms))
    has30 = any("full30" in summ[a] for a in arms)
    ncol = 3 if has30 else 2
    fig, ax = plt.subplots(1, ncol, figsize=(13 + 5 * has30, 0.5 * len(arms) + 2.2), sharey=True)
    ax[0].barh(y, ev6, color="#4c78a8")
    for i, (e, w) in enumerate(zip(ev6, wrong6)):
        ax[0].text(e + 0.6, i, f"{e:.1f}%  wrong {w}/6", va="center", fontsize=8)
    ax[0].axvline(100, color="k", lw=0.8)
    ax[0].set_xlim(0, 125)
    ax[0].set_yticks(y)
    ax[0].set_yticklabels(
        [f"{a}\n{summ[a]['six']['spec']}" if len(summ[a]["six"]["spec"]) < 28 else a for a in arms],
        fontsize=7,
    )
    ax[0].invert_yaxis()
    ax[0].set_xlabel("evaluations, six full-sample voxels (% of the baseline's)")
    ax[0].set_title("Evaluations (six voxels, clean images)")
    ax[1].barh(y, wr, color="#e45756")
    ax[1].hlines(y, lo, hi, color="k", lw=1.2)
    for i, (w, b, f) in enumerate(zip(wr, br, fx)):
        ax[1].text(hi[i] + 1, i, f"{int(w)}   broken {b} / fixed {f}", va="center", fontsize=8)
    ax[1].set_xlim(0, 150)
    ax[1].set_xlabel("wrong voxels of 200 (E0, clean images; line = Wilson interval)")
    ax[1].set_title("Wrong voxels (E0, 200), paired against the baseline")
    if has30:
        e30 = [
            100 * summ[a]["full30"]["evals_ratio"] if "full30" in summ[a] else np.nan for a in arms
        ]
        ax[2].barh(y, e30, color="#54a24b")
        for i, a in enumerate(arms):
            f = summ[a].get("full30")
            if f:
                ax[2].text(
                    e30[i] + 0.6,
                    i,
                    f"{e30[i]:.1f}%  wrong {f['wrong']}/30 "
                    f"(broken {f['broken']}, fixed {f['fixed']})",
                    va="center",
                    fontsize=8,
                )
        ax[2].axvline(100, color="k", lw=0.8)
        ax[2].set_xlim(0, 150)
        ax[2].set_xlabel("evaluations, 30 full-sample voxels (% of the baseline's)")
        ax[2].set_title("30 full-sample voxels (clean images)")
    fig.tight_layout()
    (BASE / "figures").mkdir(exist_ok=True)
    fn = BASE / "figures" / "arms_evals_and_wrong.png"
    fig.savefig(fn, dpi=130)
    print(fn)


if __name__ == "__main__":
    main()
