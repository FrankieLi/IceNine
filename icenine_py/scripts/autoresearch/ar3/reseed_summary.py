"""AR3 step 8: the summary of parts 1 and 2 (reads only committed files under ``reseed/slim``).

  uv run python scripts/autoresearch/ar3/reseed_summary.py

Writes ``reseed/summary.json``, ``reseed/tables.md`` (marker blocks for ``sync_doc_tables.py``) and
the plots ``reseed/ar3_reseed_split.png`` (evaluations of one seed step by part, per level, on the
two samples) and ``reseed/ar3_reseed_classes.png`` (the classes of rejected fits per level).

Evaluations are raw counts of strict-cost calls. A seed step: costing of the candidates plus one
refinement (local optimization) per refined candidate; indexing is not counted.
"""

import gzip
import json
from pathlib import Path
from typing import Any, Dict, List, Optional

import ar3_common  # noqa: F401  (sets sys.path for doc_tables and stats)
import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from doc_tables import markdown_table, write_tables  # noqa: E402
from stats import wilson  # noqa: E402

HERE = Path(__file__).resolve().parent
OUT = HERE.parents[2] / "benchmarks" / "autoresearch" / "ar3" / "reseed"
SLIM = OUT / "slim"
TEST = 0.80
RIGHT = 1.0
CASES = [("big", "S5L"), ("big", "D4"), ("grad", "S5"), ("grad", "D4")]
NAME = {("big", "S5L"): "big S5", ("big", "D4"): "big D4", ("grad", "S5"): "grad S5"}
NAME[("grad", "D4")] = "grad D4"
START_BINS = [0.0, 0.1, 0.2, 0.5, np.inf]
START_LAB = ["start < 0.1", "0.1 to 0.2", "0.2 to 0.5", "> 0.5"]


def rj(name: str) -> Any:
    p = SLIM / name
    if name.endswith(".gz"):
        with gzip.open(p, "rt") as f:
            return json.load(f)
    return json.loads(p.read_text())


def step_stats(rows: List[Dict[str, Any]]) -> Dict[str, Any]:
    n = len(rows)
    if n == 0:
        return dict(n=0)
    nref = max(len(r["fin_evals"]) for r in rows)
    ref_means = [
        float(np.mean([r["fin_evals"][i] if i < len(r["fin_evals"]) else 0 for r in rows]))
        for i in range(nref)
    ]
    ev = np.array([r["evals"] for r in rows], float)
    costing = np.array([r["costed"] for r in rows], float)
    first_ok: List[Optional[int]] = []
    for r in rows:
        first_ok.append(
            next(
                (
                    i
                    for i in range(len(r["fin_hit"]))
                    if r["fin_hit"][i] >= TEST and r["fin_err"][i] <= RIGHT
                ),
                None,
            )
        )
    after = np.array(
        [0 if i is None else sum(r["fin_evals"][i + 1 :]) for r, i in zip(rows, first_ok)], float
    )
    kept_pos = [int(np.argmin(r["fin_costs"])) for r in rows]
    k0 = sum(1 for i in first_ok if i == 0)
    some = sum(1 for i in first_ok if i is not None)
    earlier = sum(
        1 for i, k in zip(first_ok, kept_pos) if i is not None and i < len(rows[0]["fin_evals"]) - 1
    )
    lo, hi = wilson(k0, n)
    return dict(
        n=n,
        evals=int(ev.sum()),
        mean_evals=float(ev.mean()),
        costing_mean=float(costing.mean()),
        costing_share=float(costing.sum() / ev.sum()),
        ref_means=ref_means,
        ref_share=float(1 - costing.sum() / ev.sum()),
        first_ok_at0=k0,
        first_ok_at0_ci=[lo, hi],
        some_ok=some,
        after_ok_evals=int(after.sum()),
        after_ok_share=float(after.sum() / ev.sum()),
        kept_hist=[int(sum(1 for k in kept_pos if k == j)) for j in range(nref)],
        right_kept=int(sum(1 for r, k in zip(rows, kept_pos) if r["fin_err"][k] <= RIGHT)),
        pass_kept=int(sum(1 for r, k in zip(rows, kept_pos) if r["fin_hit"][k] >= TEST)),
        earlier_ok=earlier,
    )


def load_steps() -> Dict[Any, List[Dict[str, Any]]]:
    d = {}
    for c in CASES:
        p = f"steps_{c[0]}_{c[1]}.json.gz"
        if (SLIM / p).exists():
            d[c] = rj(p)
    return d


def main() -> None:
    steps = load_steps()
    summ: Dict[str, Any] = {}
    tabs: Dict[str, str] = {}

    # ---- part 1.1: one seed step, split
    sp: Dict[str, Any] = {}
    rows_t: List[Dict[str, Any]] = []
    for c, rows in steps.items():
        s = step_stats(rows)
        sp[f"{c[0]}/{c[1]}"] = s
        d = dict(run=NAME[c], steps=s["n"], evals=s["evals"], mean=s["mean_evals"])
        d["costing"] = s["costing_mean"]
        for i in range(5):
            d[f"r{i + 1}"] = s["ref_means"][i] if i < len(s["ref_means"]) else None
        d["refine_share"] = 100.0 * s["ref_share"]
        d["first_ok"] = f"{s['first_ok_at0']} of {s['n']}"
        d["first_ok_ci"] = f"{s['first_ok_at0_ci'][0]:.2f} to {s['first_ok_at0_ci'][1]:.2f}"
        d["after_ok"] = 100.0 * s["after_ok_share"]
        rows_t.append(d)
    summ["steps"] = sp
    tabs["step_split"] = markdown_table(
        rows_t,
        [
            "run",
            "steps",
            "evals",
            "mean",
            "costing",
            "r1",
            "r2",
            "r3",
            "r4",
            "r5",
            "refine_share",
            "first_ok",
            "first_ok_ci",
            "after_ok",
        ],
        formats=dict(
            evals="d",
            mean=".0f",
            costing=".1f",
            r1=".0f",
            r2=".0f",
            r3=".0f",
            r4=".0f",
            r5=".0f",
            refine_share=".1f",
            after_ok=".1f",
        ),
    )

    # by level
    lvl: Dict[str, Any] = {}
    rows_t = []
    for c, rows in steps.items():
        for L in sorted({r["level"] for r in rows}):
            sel = [r for r in rows if r["level"] == L]
            s = step_stats(sel)
            lvl[f"{c[0]}/{c[1]}/{L}"] = s
            d = dict(
                run=NAME[c], level=L, steps=s["n"], mean=s["mean_evals"], costing=s["costing_mean"]
            )
            for i in range(5):
                d[f"r{i + 1}"] = s["ref_means"][i] if i < len(s["ref_means"]) else None
            d["first_ok"] = f"{s['first_ok_at0']} of {s['n']}"
            d["after_ok"] = 100.0 * s["after_ok_share"]
            rows_t.append(d)
    summ["steps_by_level"] = lvl
    tabs["step_by_level"] = markdown_table(
        rows_t,
        [
            "run",
            "level",
            "steps",
            "mean",
            "costing",
            "r1",
            "r2",
            "r3",
            "r4",
            "r5",
            "first_ok",
            "after_ok",
        ],
        formats=dict(
            level=".1f",
            mean=".0f",
            costing=".1f",
            r1=".0f",
            r2=".0f",
            r3=".0f",
            r4=".0f",
            r5=".0f",
            after_ok=".1f",
        ),
    )

    plot_split(steps, lvl)
    import reseed_summary2 as S2

    s2, t2 = S2.build()
    summ.update(s2)
    tabs.update(t2)
    (OUT / "summary.json").write_text(json.dumps(summ, indent=1, default=float) + "\n")
    write_tables(OUT / "tables.md", tabs)
    print(tabs["step_split"])


def plot_split(steps: Dict[Any, List[Dict[str, Any]]], lvl: Dict[str, Any]) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(11, 7), sharey=True)
    cols = ["#999999", "#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e"]
    labs = ["costing of candidates"] + [f"refinement {i}" for i in range(1, 6)]
    for ax, c in zip(axes.ravel(), CASES):
        if c not in steps:
            ax.set_visible(False)
            continue
        levels = sorted({r["level"] for r in steps[c]})
        bottom = np.zeros(len(levels))
        parts = [[lvl[f"{c[0]}/{c[1]}/{L}"]["costing_mean"] for L in levels]]
        for i in range(5):
            parts.append(
                [
                    (
                        lvl[f"{c[0]}/{c[1]}/{L}"]["ref_means"][i]
                        if i < len(lvl[f"{c[0]}/{c[1]}/{L}"]["ref_means"])
                        else 0.0
                    )
                    for L in levels
                ]
            )
        x = np.arange(len(levels))
        for p, col, lab in zip(parts, cols, labs):
            ax.bar(x, p, bottom=bottom, color=col, label=lab)
            bottom += np.array(p)
        ax.set_xticks(x)
        ax.set_xticklabels([f"{L:g}" for L in levels])
        ax.set_xlabel("gradient level of the true grain (degrees)")
        ax.set_title(NAME[c])
        for xi, L in zip(x, levels):
            ax.text(
                xi,
                0,
                f"n={lvl[f'{c[0]}/{c[1]}/{L}']['n']}",
                ha="center",
                va="bottom",
                fontsize=7,
                rotation=90,
                color="white",
            )
    axes[0, 0].set_ylabel("evaluations per seed step (mean)")
    axes[1, 0].set_ylabel("evaluations per seed step (mean)")
    axes[0, 0].legend(fontsize=7, loc="upper left")
    fig.tight_layout()
    fig.savefig(OUT / "ar3_reseed_split.png", dpi=110)
    plt.close(fig)


if __name__ == "__main__":
    main()
