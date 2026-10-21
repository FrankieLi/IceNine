"""Summary of the Task F carry-option comparison (D-20261021-N1).

Reads ``benchmarks/eval_reduction/Fcompare/cache/<set>/`` (runs, contended) and
``cache/timing_<set>/`` (the single-worker pass on a subset, or absent), writes
``summary_<set>.json``, ``tables_<set>.md`` (marker-delimited blocks), the figure and
``doc_tables_all.md``.

Eight options (rows): F4 (1-frame window), F4 (5-frame window), F5, F5 then F4, policy 1 (F4,
then F5+F4 where F4's own acceptance test rejects it; 1-frame window in both F4 parts), policy 5
(the same with the 5-frame window in both F4 parts), CMA-ES sigma0 1 deg, CMA-ES sigma0 2 deg. A
policy's cost is F4's plus, where F4 was rejected, F5+F4's (F5+F4 is run from the original start);
its end point is F4's, or F5+F4's where F4 was rejected.

Captured = final orientation within the threshold of the truth (tolerance 1e-3 deg, as Task F).
Strict-equivalent evaluations use Task F's current weights (``F/timing/<set>.json`` and
``<set>_f5.json``, from Task F's timing voxels, not this comparison's).

  uv run python scripts/eval_reduction/outward_compare_summary.py e0
  uv run python scripts/eval_reduction/outward_compare_summary.py full
  uv run python scripts/eval_reduction/outward_compare_summary.py merge
  uv run python scripts/eval_reduction/outward_compare_summary.py figure
"""

import argparse
import json
import math
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

import numpy as np

HERE = Path(__file__).resolve().parent
ICE = HERE.parents[1]
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(ICE / "scripts" / "common"))
import doc_tables  # noqa: E402
import outward_summary as OS  # noqa: E402
import stats as S  # noqa: E402

OUT = ICE / "benchmarks" / "eval_reduction" / "Fcompare"
F_OUT = ICE / "benchmarks" / "eval_reduction" / "F"
OFFSETS = [1.0, 2.0, 3.0, 5.0]
OFFSET_IDX = [2, 3, 4, 5]
TOL = 1e-3
# option key -> (label, kind, parts); kind "run" = one method's record; "policy" = (first, fallback)
OPTIONS: Dict[str, Tuple[str, str, Tuple[str, ...]]] = {
    "F4": ("F4, 1-frame window", "run", ("f4",)),
    "F4w": ("F4, 5-frame window", "run", ("f4w",)),
    "F5": ("F5 expanding shells", "run", ("f5",)),
    "F54": ("F5 then F4", "run", ("f54",)),
    "P1": ("Policy, 1-frame window", "policy", ("f4", "f54")),
    "P5": ("Policy, 5-frame window", "policy", ("f4w", "f54w")),
    "C1": ("CMA-ES, sigma0 1 deg", "run", ("f2_1",)),
    "C2": ("CMA-ES, sigma0 2 deg", "run", ("f2_2",)),
}
COLOR = {"F4": "#d62728", "F4w": "#e377c2", "F5": "#ff7f0e", "F54": "#9467bd", "P1": "#2ca02c",
         "P5": "#17becf", "C1": "#1f77b4", "C2": "#08306b"}  # fmt: skip
MARK = {"F4": "o", "F4w": "o", "F5": "s", "F54": "s", "P1": "^", "P5": "^", "C1": "D", "C2": "D"}
PAIRS = [  # (name, base option, other option): fixed = base wrong and other right
    ("Policy 1 against F4 (1-frame)", "F4", "P1"),
    ("Policy 5 against F4 (5-frame)", "F4w", "P5"),
    ("F4 5-frame against F4 1-frame", "F4", "F4w"),
    ("Policy 5 against policy 1", "P1", "P5"),
]
SET_LABEL = {"e0": "E0 (per-voxel clean images)", "full": "full sample (clean images)"}


# --------------------------------------------------------------------------------------------
# loading and per-option arrays
# --------------------------------------------------------------------------------------------
def load(tag: str) -> Dict[str, np.ndarray]:
    cols: Dict[str, List[np.ndarray]] = {}
    for f in sorted((OUT / "cache" / tag).glob("v*_o*.npz")):
        z = np.load(f)
        v, oi = f.stem[1:].split("_o")
        n = len(z["err_start"])
        cols.setdefault("voxel", []).append(np.full(n, int(v)))
        cols.setdefault("oi", []).append(np.full(n, int(oi)))
        cols.setdefault("dir", []).append(np.arange(n))
        for k in z.files:
            if "|" in k or k == "err_start":
                cols.setdefault(k, []).append(z[k])
    return {k: np.concatenate(v) for k, v in cols.items()}


def weights_for(kind: str) -> Dict[str, float]:
    tm = json.loads((F_OUT / "timing" / f"{kind}.json").read_text())
    f5t = OS.f5_timing(kind, tm["voxels"])
    return OS.weights(tm["ms_per_evaluation"], f5t)


def method_arrays(D: Dict[str, np.ndarray], m: str, w: Dict[str, float]) -> Dict[str, np.ndarray]:
    """Per-run counts by kind, strict-equivalent count, seconds and end point error of method m."""
    n_jac = D[f"{m}|n_jac"]
    a = dict(
        err=D[f"{m}|err"],
        strict=D[f"{m}|n_strict"].astype(float),
        screen=D[f"{m}|n_screen"].astype(float),
        value=(D[f"{m}|n_smooth"] - n_jac).astype(float),
        jac=n_jac.astype(float),
        secs=D[f"{m}|secs"],
        flag=D[f"{m}|flag"],
    )
    if m in ("f5", "f54", "f54w"):
        se5 = D["f5|n_strict"] + w["screen_f5"] * D["f5|n_screen"] + w["val_f5"] * D["f5|n_smooth"]
        if m == "f5":
            se = se5
        else:
            assert np.array_equal(D[f"{m}|n_screen"], D["f5|n_screen"])
            n_val = D[f"{m}|n_smooth"] - D["f5|n_smooth"] - n_jac
            se = se5 + (D[f"{m}|n_strict"] - D["f5|n_strict"]) + w["val"] * n_val + w["jac"] * n_jac
    else:
        se = a["strict"] + w["screen"] * a["screen"] + w["val"] * a["value"] + w["jac"] * a["jac"]
    a["se"] = se
    return a


def option_arrays(
    D: Dict[str, np.ndarray], w: Dict[str, float]
) -> Dict[str, Dict[str, np.ndarray]]:
    base = {m: method_arrays(D, m, w) for m in ("f2_1", "f2_2", "f4", "f4w", "f5", "f54", "f54w")}
    out: Dict[str, Dict[str, np.ndarray]] = {}
    for key, (_, kind, parts) in OPTIONS.items():
        if kind == "run":
            out[key] = dict(base[parts[0]], fallback=np.zeros(len(D["voxel"]), dtype=bool))
            continue
        first, fb = base[parts[0]], base[parts[1]]
        rej = first["flag"] < 0.5
        o = {k: np.where(rej, fb["err"], first["err"]) if k == "err" else first[k] + rej * fb[k]
             for k in ("err", "strict", "screen", "value", "jac", "secs", "se")}  # fmt: skip
        o["fallback"] = rej
        o["flag"] = np.where(rej, fb["flag"], first["flag"])
        out[key] = o
    return out


# --------------------------------------------------------------------------------------------
# formatting helpers
# --------------------------------------------------------------------------------------------
def share(k: int, n: int) -> str:
    lo, hi = S.wilson(k, n)
    return f"{k} of {n}, {100 * k / n:.1f}% ({100 * lo:.1f}-{100 * hi:.1f})"


def q(x: np.ndarray, p: float) -> float:
    x = x[np.isfinite(x)]
    return float(np.percentile(x, p)) if len(x) else math.nan


def fnum(x: float) -> str:
    if not np.isfinite(x):
        return "n/a"
    return f"{x:.0f}" if abs(x) >= 100 else (f"{x:.1f}" if abs(x) >= 10 else f"{x:.2f}")


def captured(err: np.ndarray, thr: float) -> np.ndarray:
    return err < thr - TOL


# --------------------------------------------------------------------------------------------
# summary of one set
# --------------------------------------------------------------------------------------------
def summarize(kind: str) -> None:
    D = load(kind)
    w = weights_for(kind)
    O = option_arrays(D, w)
    # cross-check the strict-equivalent counts against Task F's own function (F4, F4w, F5, F54)
    for m, key in (("f4", "F4"), ("f4w", "F4w"), ("f5", "F5"), ("f54", "F54")):
        ref = OS.strict_equiv(D, m, w, math.nan)
        assert np.allclose(ref, O[key]["se"]), m
    T: Optional[Dict[str, np.ndarray]] = None
    pf: Optional[Dict[str, Any]] = None
    if (OUT / "cache" / f"timing_{kind}").exists():
        T = load(f"timing_{kind}")
        OT = option_arrays(T, w)
        pf = json.loads((OUT / "timing" / f"preflight_{kind}.json").read_text())
    nvox = len(set(int(v) for v in D["voxel"]))
    summ: Dict[str, Any] = dict(
        set=kind, n_voxels=nvox, weights=w, timing_label=("single-worker" if T else "none"),
        timing_voxels=sorted(set(int(v) for v in T["voxel"])) if T else [],
    )  # fmt: skip
    tables: Dict[str, str] = {}
    pre = f"evalFc_{kind}_"

    # time of each option per run: single-worker where the timing pass exists, else contended
    def time_cell(key: str, oi: int) -> Tuple[float, float, int]:
        if T is None:
            return math.nan, math.nan, 0
        s = T["oi"] == oi
        return q(OT[key]["secs"][s], 50), q(OT[key]["secs"][s], 90), int(s.sum())

    main_rows: Dict[float, List[Dict[str, Any]]] = {}
    raw_rows: Dict[float, List[Dict[str, Any]]] = {}
    for off, oi in zip(OFFSETS, OFFSET_IDX):
        sel = D["oi"] == oi
        n = int(sel.sum())
        mr, rr = [], []
        for key, (label, _, _) in OPTIONS.items():
            o = O[key]
            k1 = int(captured(o["err"][sel], 1.0).sum())
            k005 = int(captured(o["err"][sel], 0.05).sum())
            med_t, p90_t, nt = time_cell(key, oi)
            fb = int(o["fallback"][sel].sum()) if OPTIONS[key][1] == "policy" else None
            mr.append({
                "option": label, "within 1 deg": share(k1, n), "within 0.05 deg": share(k005, n),
                "fell back": "n/a" if fb is None else f"{fb} of {n}",
                "strict-equivalent, median": fnum(q(o["se"][sel], 50)),
                "strict-equivalent, 90th pct": fnum(q(o["se"][sel], 90)),
                "seconds per run, median": fnum(med_t), "seconds per run, 90th pct": fnum(p90_t),
            })  # fmt: skip
            rr.append({
                "option": label,
                "strict": fnum(q(o["strict"][sel], 50)), "screen": fnum(q(o["screen"][sel], 50)),
                "smooth value": fnum(q(o["value"][sel], 50)),
                "smooth Jacobian": fnum(q(o["jac"][sel], 50)),
                "strict-equivalent, median": fnum(q(o["se"][sel], 50)),
                "strict-equivalent, 90th pct": fnum(q(o["se"][sel], 90)),
            })  # fmt: skip
            summ.setdefault("options", {}).setdefault(key, {})[f"{off:g}"] = dict(
                n=n, within_1=k1, within_005=k005, fallback=fb,
                se_median=q(o["se"][sel], 50), se_p90=q(o["se"][sel], 90),
                raw_median={c: q(o[c][sel], 50) for c in ("strict", "screen", "value", "jac")},
                secs_single_median=med_t, secs_single_p90=p90_t, n_timing_runs=nt,
                secs_contended_median=q(o["secs"][sel], 50),
            )  # fmt: skip
        main_rows[off], raw_rows[off] = mr, rr
        tables[f"{pre}main_{off:g}deg"] = doc_tables.markdown_table(mr, list(mr[0]))
        tables[f"{pre}raw_{off:g}deg"] = doc_tables.markdown_table(rr, list(rr[0]))

    # ---- the option with the lowest median cost among those whose share within 1 deg (point
    # estimate, no interval) reaches a level; ties in cost: the one with the larger share
    crows = []
    for off in OFFSETS:
        od = {key: summ["options"][key][f"{off:g}"] for key in OPTIONS}
        for level in (100.0, 99.0, 95.0, 90.0):
            row: Dict[str, Any] = {
                "offset (deg)": f"{off:g}",
                "share within 1 deg at least": f"{level:g}%",
            }
            for col, field in (("fewest strict-equivalent evaluations (median)", "se_median"),
                               ("fewest seconds per run (median)", "secs_single_median"),
            ):  # fmt: skip
                ok = [k for k, d in od.items()
                      if 100 * d["within_1"] / d["n"] >= level - 1e-9
                      and np.isfinite(d[field])]  # fmt: skip
                if not ok:
                    row[col] = "none reaches it" if field == "se_median" else "n/a"
                    continue
                best = min(ok, key=lambda k: (od[k][field], -od[k]["within_1"]))
                row[col] = f"{OPTIONS[best][0]} ({fnum(od[best][field])})"
            crows.append(row)
    summ["cheapest"] = crows
    tables[pre + "cheapest"] = doc_tables.markdown_table(crows, list(crows[0]))

    # ---- per-voxel counts within 1 deg (runs cluster by voxel), voxels below 20 of 20 for any
    vox = sorted(set(int(v) for v in D["voxel"]))
    for off, oi in zip(OFFSETS, OFFSET_IDX):
        rows = []
        for v in vox:
            sel = (D["oi"] == oi) & (D["voxel"] == v)
            ks = {key: int(captured(O[key]["err"][sel], 1.0).sum()) for key in OPTIONS}
            nv = int(sel.sum())
            if any(k < nv for k in ks.values()):
                rows.append({"voxel": v, **{key: f"{k} of {nv}" for key, k in ks.items()}})
        summ.setdefault("per_voxel_below_full", {})[f"{off:g}"] = rows
        if rows:
            tables[f"{pre}pervoxel_{off:g}deg"] = doc_tables.markdown_table(rows, list(rows[0]))
        else:
            tables[f"{pre}pervoxel_{off:g}deg"] = (
                f"Every option captured all runs of every voxel at {off:g} deg."
            )
        # how many voxels have any run missed, per option
        summ.setdefault("voxels_with_a_miss", {})[f"{off:g}"] = {
            key: sum(
                1 for v in vox
                if not captured(O[key]["err"][(D["oi"] == oi) & (D["voxel"] == v)], 1.0).all()
            )
            for key in OPTIONS
        }  # fmt: skip

    # ---- paired, run by run, within 1 deg
    prow = []
    for off, oi in zip(OFFSETS, OFFSET_IDX):
        sel = D["oi"] == oi
        for name, a, b in PAIRS:
            ca, cb = captured(O[a]["err"], 1.0), captured(O[b]["err"], 1.0)
            fixed, broken = ~ca & cb & sel, ca & ~cb & sel
            both_wrong, both_right = ~ca & ~cb & sel, ca & cb & sel
            vf = len(set(D["voxel"][fixed])), len(set(D["voxel"][broken]))
            prow.append({
                "offset (deg)": f"{off:g}", "comparison": name, "fixed": int(fixed.sum()),
                "broken": int(broken.sum()), "both wrong": int(both_wrong.sum()),
                "both right": int(both_right.sum()),
                "voxels with a fixed run": vf[0], "voxels with a broken run": vf[1],
            })  # fmt: skip
            summ.setdefault("paired", {}).setdefault(f"{off:g}", {})[name] = prow[-1]
    tables[pre + "paired"] = doc_tables.markdown_table(prow, list(prow[0]))

    # ---- contended against single-worker seconds on the timing runs
    if T is not None:
        crow = []
        for off, oi in zip(OFFSETS, OFFSET_IDX):
            s = T["oi"] == oi
            for key, (label, _, _) in OPTIONS.items():
                sw, ct = q(OT[key]["secs"][s], 50), None
                # the same runs in the main (contended) pass
                ct_arr = []
                for v, d in set(zip(T["voxel"][s].tolist(), T["dir"][s].tolist())):
                    m_ = (D["oi"] == oi) & (D["voxel"] == v) & (D["dir"] == d)
                    ct_arr.append(float(O[key]["secs"][m_][0]))
                ct = float(np.median(ct_arr))
                crow.append({"offset (deg)": f"{off:g}", "option": label,
                             "single-worker, median (s)": fnum(sw),
                             "contended, same runs, median (s)": fnum(ct),
                             "ratio of medians": f"{ct / sw:.2f}" if sw > 0 else "n/a",
                             })  # fmt: skip
        tables[pre + "contention"] = doc_tables.markdown_table(crow, list(crow[0]))
        # does the timing pass reproduce the main run's results?
        same = {}
        for key in OPTIONS:
            eq = []
            for i in range(len(T["voxel"])):
                m_ = (
                    (D["oi"] == T["oi"][i])
                    & (D["voxel"] == T["voxel"][i])
                    & (D["dir"] == T["dir"][i])
                )
                eq.append(bool(np.isclose(O[key]["err"][m_][0], OT[key]["err"][i], atol=1e-6)
                               and O[key]["se"][m_][0] == OT[key]["se"][i]))  # fmt: skip
            same[key] = f"{sum(eq)} of {len(eq)}"
        summ["timing_reproduces_main"] = same
        summ["timing_preflight"] = dict(
            loadavg_at_start=pf.get("loadavg") if pf else None,
            busy_processes=pf.get("busy_processes") if pf else None,
        )  # fmt: skip

    (OUT / f"summary_{kind}.json").write_text(json.dumps(summ, indent=1, default=float) + "\n")
    doc_tables.write_tables(OUT / f"tables_{kind}.md", tables)
    print(f"{kind}: {len(tables)} tables")


def figure_all() -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.ticker

    fdir = OUT / "figures"
    fdir.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(2, 4, figsize=(17, 8), sharey=True, sharex="row")
    for r, kind in enumerate(("e0", "full")):
        s = json.loads((OUT / f"summary_{kind}.json").read_text())
        for c, off in enumerate(OFFSETS):
            ax = axes[r][c]
            for key in OPTIONS:
                d = s["options"][key][f"{off:g}"]
                x = (
                    d["secs_single_median"]
                    if s["timing_label"] == "single-worker"
                    else d["secs_contended_median"]
                )
                lo, hi = S.wilson(d["within_1"], d["n"])
                y = 100 * d["within_1"] / d["n"]
                ax.errorbar(x, y, yerr=[[y - 100 * lo], [100 * hi - y]], color=COLOR[key],
                            fmt=MARK[key],
                            mfc=COLOR[key] if key not in ("F4w", "P5", "F54", "C2") else "white",
                            ms=8, capsize=2, label=OPTIONS[key][0])  # fmt: skip
            ax.set_xscale("log")
            ax.set_ylim(25, 102)
            ax.xaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, _: f"{v:g}"))
            ax.xaxis.set_minor_formatter(
                matplotlib.ticker.FuncFormatter(
                    lambda v, _: f"{v:g}" if v in (0.03, 0.05, 0.2, 0.3, 0.5, 2, 3) else ""
                )
            )
            ax.tick_params(axis="x", which="both", labelsize=7)
            ax.set_title(f"{SET_LABEL[kind]}, {off:g} deg", fontsize=9)
            ax.set_xlabel(
                "median seconds per run"
                + (" (single-worker)" if s["timing_label"] == "single-worker" else " (contended)"),
                fontsize=8,
            )
            ax.set_ylabel("share within 1 deg (%)", fontsize=8)
            ax.grid(alpha=0.3)
    axes[0][0].legend(fontsize=7, loc="lower right")
    fig.tight_layout()
    fig.savefig(fdir / "capture_vs_seconds.png", dpi=130)
    plt.close(fig)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("kind", choices=["e0", "full", "merge", "figure"])
    a = ap.parse_args()
    if a.kind == "merge":
        merged: Dict[str, str] = {}
        for f in sorted(OUT.glob("tables_*.md")):
            merged.update(doc_tables.read_tables(f))
        doc_tables.write_tables(OUT / "doc_tables_all.md", merged)
        print(f"{len(merged)} tables -> doc_tables_all.md")
    elif a.kind == "figure":
        figure_all()
    else:
        summarize(a.kind)


if __name__ == "__main__":
    main()
