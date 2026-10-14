"""Task B: doc tables from summary.json of cost_edges_summary.py (doc_tables blocks in tables.md).

uv run python scripts/eval_reduction/cost_edges_tables.py --out benchmarks/eval_reduction/B
"""

import argparse
import json
import sys
from collections import Counter
from pathlib import Path
from typing import Any, Dict, List

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[0] / "common"))
import cost_edges as CE  # noqa: E402
import doc_tables as DT  # noqa: E402

REG = ["fin", "mc", "grid"]
REGL_FULL = {"fin": "under 0.05 deg", "mc": "0.05 to 0.5 deg", "grid": "0.5 to 5 deg"}
# E0 samples 0.5 to 5 deg only in four 0.05-deg windows per ray (cost_edges_run.segments)
REGL_E0 = dict(REGL_FULL, grid="0.5 to 5 deg (4 windows of 0.05 deg)")
WALKL = {
    ("full", "dense"): "dense, 0 to 5 deg",
    ("full", "dense_to_0.5"): "dense, 0 to 0.5 deg only",
    ("e0", "dense"): "dense, 0 to 0.5 deg",
    ("e0", "profile"): "log-spaced, 121 points",
}
KINDL = {"full": "full-sample voxels (6)", "e0": "E0 voxels (30)"}
RULEL = {"old": "lenient", "cpp": "C++", "cppE": "C++ + CppEdgeRules", "r2": "ray rule R2"}
SETL = {"full": "complete cost", "screen": "level-0 screen"}
TYPES = CE.EVENT_TYPES + ["denominator"]


def f(x: float) -> str:
    if x == 0:
        return "0"
    if x >= 100:
        return f"{x:,.0f}"
    if x >= 10:
        return f"{x:.0f}"
    if x >= 1:
        return f"{x:.1f}"
    return f"{x:.2f}"


def wci(k: int, n: int, w: List[float]) -> str:
    """'k of n (lo-hi%)', 95% Wilson interval from summary.json."""
    return f"{k} of {n} ({100 * w[0]:.1f}-{100 * w[1]:.1f}%)"


def jump_bin(x: float) -> str:
    """The quarter-decade bin of a jump quantile (stored as the bin centre), as its range."""
    lo, hi = x / 10**0.125, x * 10**0.125
    return f"{lo:.2g}-{hi:.2g}"


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    out = Path(a.out).resolve()
    S = json.loads((out / "summary.json").read_text())
    T: Dict[str, str] = {}
    for kind in ("full", "e0"):
        if kind not in S:
            continue
        REGL = REGL_E0 if kind == "e0" else REGL_FULL
        # events per degree by type, C++ rule, complete cost and screen
        rows = []
        for t in CE.EVENT_TYPES:
            row: Dict[str, Any] = {"event type": t, "kind": CE.KIND[t]}
            for setting in ("full", "screen"):
                s = S[kind][f"{setting}|cpp"]
                for r in REG:
                    row[f"{SETL[setting]}, {REGL[r]}"] = (
                        f"{f(s[r]['cost_changing_events_per_deg'][t])} of "
                        f"{f(s[r]['events_per_deg'][t])}"
                    )
            rows.append(row)
        cols = ["event type", "kind"] + [
            f"{SETL[s]}, {REGL[r]}" for s in ("full", "screen") for r in REG
        ]
        T[f"evalB_edges_per_deg_{kind}"] = DT.markdown_table(rows, cols)
        # rules comparison
        rows = []
        for setting in ("full", "screen"):
            for rule in ("old", "cpp", "cppE", "r2"):
                s = S[kind][f"{setting}|{rule}"]
                row = {"cost": SETL[setting], "rule": RULEL[rule]}
                for r in REG:
                    x = s[r]
                    cc = sum(x["cost_changing_events_per_deg"].values())
                    row[f"{REGL[r]}: cost-changing events/deg"] = f(cc)
                    row[f"{REGL[r]}: steps with a cost change"] = (
                        f"{x['steps_cost_changed']} of {x['steps']}"
                    )
                    row[f"{REGL[r]}: steps on the cost = 1 plateau"] = f"{x['plateau_steps']}"
                rows.append(row)
        cols = ["cost", "rule"] + [
            f"{REGL[r]}: {k}"
            for r in REG
            for k in (
                "cost-changing events/deg",
                "steps with a cost change",
                "steps on the cost = 1 plateau",
            )
        ]
        T[f"evalB_rules_{kind}"] = DT.markdown_table(rows, cols)
        # gross change share by cause (C++ rule, complete cost and screen)
        rows = []
        for setting in ("full", "screen"):
            for r in REG:
                x = S[kind][f"{setting}|cpp"][r]
                tot = sum(x["gross"].values())
                row = {
                    "cost": SETL[setting],
                    "regime": REGL[r],
                    "net variation": f(x["net_variation"]),
                    "gross change": f(tot),
                }
                for c in TYPES + ["unattributed"]:
                    row[c] = f"{100 * x['gross'][c] / tot:.1f}%" if tot else "n/a"
                rows.append(row)
        T[f"evalB_gross_share_{kind}"] = DT.markdown_table(
            rows, ["cost", "regime", "net variation", "gross change"] + TYPES + ["unattributed"]
        )
        # jump sizes
        rows = []
        for c in ["frame", "pixel", "detector", "eta_cut", "no_frame", "bragg_limit"]:
            row = {"type": c}
            for r in REG:
                q = S[kind]["full|cpp"][r]["jump_quantiles"][c]
                row[REGL[r]] = (
                    f"{jump_bin(q['p50'])} / {jump_bin(q['p90'])} (n={q['n']})"
                    if q["n"]
                    else "none"
                )
            rows.append(row)
        T[f"evalB_jumps_{kind}"] = DT.markdown_table(rows, ["type"] + [REGL[r] for r in REG])
        # denominator
        rows = []
        for setting in ("full", "screen"):
            for rule in ("old", "cpp", "r2"):
                row = {"cost": SETL[setting], "rule": RULEL[rule]}
                for r in REG:
                    x = S[kind][f"{setting}|{rule}"][r]
                    row[REGL[r]] = (
                        f"{x['steps_denominator_only']} of {x['steps_denominator_changed']}"
                    )
                rows.append(row)
        T[f"evalB_denominator_{kind}"] = DT.markdown_table(
            rows, ["cost", "rule"] + [REGL[r] for r in REG]
        )
        # qualification flips by cause
        rows = []
        for rule in ("old", "cpp", "r2"):
            for setting in ("full", "screen"):
                row = {"cost": SETL[setting], "rule": RULEL[rule]}
                tot = Counter()
                deg = 0.0
                for r in REG:
                    x = S[kind][f"{setting}|{rule}"][r]
                    deg += x["degrees"]
                    tot.update(x["flips"])
                for t in [
                    "frame",
                    "pixel",
                    "detector",
                    "eta_cut",
                    "no_frame",
                    "rule",
                    "unattributed",
                ]:
                    row["cause: " + t] = f"{tot[t]}"
                row["per degree"] = f(sum(tot.values()) / deg)
                rows.append(row)
        T[f"evalB_flips_{kind}"] = DT.markdown_table(
            rows,
            ["cost", "rule"]
            + [
                "cause: " + t
                for t in [
                    "frame",
                    "pixel",
                    "detector",
                    "eta_cut",
                    "no_frame",
                    "rule",
                    "unattributed",
                ]
            ]
            + ["per degree"],
        )
        # minima
        mk = kind + "_minima"
        if mk in S:
            rows = []
            for setting in ("full", "screen"):
                for rule in ("old", "cpp", "r2"):
                    segs = ("dense", "profile") if kind == "e0" else ("dense", "dense_to_0.5")
                    for seg in segs:
                        m = S[mk][f"{setting}|{rule}|{seg}"]
                        for p, v in m["by_prominence"].items():
                            if not v["n_rays"]:
                                continue
                            ci = Counter(v["cause_in"])
                            top = ", ".join(f"{k} {n}" for k, n in ci.most_common(3)) or "none"
                            rows.append(
                                {
                                    "cost": SETL[setting],
                                    "rule": RULEL[rule],
                                    "walk": WALKL[(kind, seg)],
                                    "prominence": p,
                                    "rays with a minimum (95% Wilson)": wci(
                                        v["with_min"], v["n_rays"], v["with_min_wilson95"]
                                    ),
                                    "minima": v["n_min"],
                                    "median depth": (
                                        f"{np.median(v['depth']):.4f}" if v["depth"] else "n/a"
                                    ),
                                    "median width (deg)": (
                                        f"{np.median(v['width']):.3f}"
                                        if v["width"]
                                        else ("n/a (log grid)" if seg == "profile" else "n/a")
                                    ),
                                    "median position (deg)": (
                                        f"{np.median(v['rpos']):.2f}" if v["rpos"] else "n/a"
                                    ),
                                    "cause of the dip (top 3)": top,
                                }
                            )
            T[f"evalB_minima_{kind}"] = DT.markdown_table(
                rows,
                [
                    "cost",
                    "rule",
                    "walk",
                    "prominence",
                    "rays with a minimum (95% Wilson)",
                    "minima",
                    "median depth",
                    "median width (deg)",
                    "median position (deg)",
                    "cause of the dip (top 3)",
                ],
            )
            rows = []
            for setting in ("full", "screen"):
                for rule in ("old", "cpp", "r2"):
                    m = S[mk][f"{setting}|{rule}|dense"]
                    below = "rays whose lowest cost is below the truth's (95% Wilson)"
                    rows.append(
                        {
                            "cost": SETL[setting],
                            "rule": RULEL[rule],
                            "walk": WALKL[(kind, "dense")],
                            below: wci(
                                m["rays_min_below_truth"],
                                m["n_rays"],
                                m["rays_min_below_truth_wilson95"],
                            ),
                        }
                    )
            T[f"evalB_truth_lowest_{kind}"] = DT.markdown_table(
                rows,
                [
                    "cost",
                    "rule",
                    "walk",
                    "rays whose lowest cost is below the truth's (95% Wilson)",
                ],
            )
    DT.write_tables(out / "tables.md", T)
    print(f"wrote {len(T)} tables to", out / "tables.md")


if __name__ == "__main__":
    main()
