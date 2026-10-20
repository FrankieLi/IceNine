"""Task F summary: capture rate, error and evaluations per method, start offset and voxel set.

Reads ``benchmarks/eval_reduction/F/cache/<tag>/v<voxel>_o<offset>.npz`` (outward_bench.py), writes
``summary_<tag>.json``, ``tables_<tag>.md`` (marker-delimited blocks for the findings doc) and the
figures. A run is *captured* at 1 deg (0.1, 0.05 deg) when the method's final orientation is within
that angle of the truth (angle of R R_true^T, no symmetry reduction, offsets are at most 5 deg).

Evaluation counts. Raw: strict + screen + smooth evaluations, each counted once (the convention of
the evaluation ledger). Strict-equivalent: each kind weighted by its single-worker time relative to
a strict evaluation, at the batch size the method uses it. F4 and the F4 part of F5+F4 use single
calls (``timing/<kind>.json``, first timing voxel; E0: strict 0.422 ms, E4 value 0.66 ms, E4
value-and-Jacobian 10.8 ms). F5 and the F5 part of F5+F4 use F5's own call pattern
(``timing/<kind>_f5.json``, same voxel: the screen a shell at a time in chunks of up to 256, E4 on
32 rotations at a time through ``evaluate_many``; E0: 0.0175 and 0.367 ms per orientation). The
first weights (screen timed in batches of 64, every E4 value as a single call) overstated F5's time
by a factor 1.2 to 1.6 against its single-worker wall time; they are kept as the ``*_v1`` tables.
The windows Gauss-Newton (f3) makes no cost evaluations; its entry is its measured single-worker
time per run divided by the strict evaluation time.
Choice rule (Task F decision point): among the methods that can run on the full sample, the capture
rate (1 deg) per strict-equivalent evaluation to stop, at 2 and at 5 deg (score = 1000 x share /
median); at each offset a method's rank is 1 + the number of methods with a strictly higher score
(ties share the better rank); the ranks at the two offsets are added, the two lowest rank sums are
carried, and ties in the sum go to the larger mean of the two scores.

  uv run python scripts/eval_reduction/outward_summary.py e0
  uv run python scripts/eval_reduction/outward_summary.py full
  uv run python scripts/eval_reduction/outward_summary.py merge   # doc_tables_all.md
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
import stats as S  # noqa: E402

OUT = ICE / "benchmarks" / "eval_reduction" / "F"
OFFSETS = [0.1, 0.5, 1.0, 2.0, 3.0, 5.0]
CKPT = [3, 10, 30, 100, 300, 1000, 3000, 10000, 30000]
METHODS = ["f1", "f2_05", "f2_1", "f2_2", "f3", "f4", "f4w", "f5", "f54"]
LABEL = {
    "f1": "F1 FindOptimal, box = offset",
    "f2_05": "F2 CMA-ES, sigma0 0.5 deg",
    "f2_1": "F2 CMA-ES, sigma0 1 deg",
    "f2_2": "F2 CMA-ES, sigma0 2 deg",
    "f3": "F3 centroid Gauss-Newton (windows)",
    "f4": "F4 Gauss-Newton on E4 rays",
    "f4w": "F4 with a 5-frame matching window",
    "f5": "F5 expanding shells",
    "f54": "F5 then F4",
}
SHORT = {"f1": "F1", "f2_05": "F2 0.5", "f2_1": "F2 1", "f2_2": "F2 2", "f3": "F3", "f4": "F4",
         "f4w": "F4 w5", "f5": "F5", "f54": "F5+F4"}  # fmt: skip
COLOR = {"f1": "#7f7f7f", "f2_05": "#aec7e8", "f2_1": "#1f77b4", "f2_2": "#08306b", "f3": "#2ca02c",
         "f4": "#d62728", "f4w": "#e377c2", "f5": "#ff7f0e", "f54": "#9467bd"}  # fmt: skip
TOL = 1e-3  # deg, as outward_bench.CAPTURE_TOL (a start exactly thr out measures thr - 1e-4)
FULL_OK = ["f1", "f2_05", "f2_1", "f2_2", "f4", "f4w", "f5", "f54"]  # run on the full sample
# the methods at the decision point (F4 with the 5-frame window was added after it)
CARRY_POOL = ["f1", "f2_05", "f2_1", "f2_2", "f4", "f5", "f54"]


def load(tag: str) -> Dict[str, Any]:
    """Concatenate the cache of one run: arrays per ``method|field`` plus voxel/offset/direction."""
    files = sorted((OUT / "cache" / tag).glob("v*_o*.npz"))
    cols: Dict[str, List[np.ndarray]] = {}
    for f in files:
        z = np.load(f)
        v, oi = f.stem[1:].split("_o")
        n = len(z["err_start"])
        cols.setdefault("voxel", []).append(np.full(n, int(v)))
        cols.setdefault("oi", []).append(np.full(n, int(oi)))
        cols.setdefault("dir", []).append(np.arange(n))
        zs = [z]
        for extra in sorted((OUT / "cache").glob(f"{tag}_f4w")):  # the f4w run of the same cases
            if (extra / f.name).exists():
                zs.append(np.load(extra / f.name))
        for zz in zs:
            for k in zz.files:
                if (zz is z or "|" in k) and zz[k].ndim >= 1 and len(zz[k]) == n:
                    cols.setdefault(k, []).append(zz[k])
    return {k: np.concatenate(v) for k, v in cols.items()}


def weights(timing: Dict[str, float], f5: Optional[Dict[str, float]] = None) -> Dict[str, float]:
    """Weights relative to a strict evaluation. ``screen_f5`` and ``val_f5`` are F5's screen and E4
    at its own batch sizes (``timing/<kind>_f5.json``); without it (the v1 weights) they equal the
    single-call ``screen`` and ``val``."""
    s = timing["strict"]
    w = dict(screen=timing["screen"] / s, val=timing["e4_value"] / s,
             jac=timing["e4_value_and_jacobian"] / s)  # fmt: skip
    w["screen_f5"] = f5["screen_ms_per_orientation"] / s if f5 else w["screen"]
    w["val_f5"] = f5["smooth_ms_per_orientation"] / s if f5 else w["val"]
    return w


def f5_timing(kind: str, timing_voxels: Sequence[int]) -> Optional[Dict[str, float]]:
    """F5's per-orientation times on the first timing voxel (the voxel of the other weights)."""
    path = OUT / "timing" / f"{kind}_f5.json"
    if not path.exists():
        return None
    pv = json.loads(path.read_text())["per_voxel"]
    return dict(pv[str(timing_voxels[0])])


def strict_equiv(D: Dict[str, Any], m: str, w: Dict[str, float], f3_se: float) -> np.ndarray:
    if m == "f3":
        return np.full(len(D["voxel"]), f3_se)
    if m in ("f5", "f54"):
        se5 = (D["f5|n_strict"] + w["screen_f5"] * D["f5|n_screen"]
               + w["val_f5"] * D["f5|n_smooth"])  # fmt: skip
        if m == "f5":
            return se5
        # F5+F4 shares the F5 run of the same case (outward_bench.run_f54); the rest is F4's
        assert np.array_equal(D["f54|n_screen"], D["f5|n_screen"])
        n_str = D["f54|n_strict"] - D["f5|n_strict"]
        n_jac = D["f54|n_jac"]
        n_val = D["f54|n_smooth"] - D["f5|n_smooth"] - n_jac
        return se5 + n_str + w["val"] * n_val + w["jac"] * n_jac
    n_val = D[f"{m}|n_smooth"] - D[f"{m}|n_jac"]
    return (D[f"{m}|n_strict"] + w["screen"] * D[f"{m}|n_screen"] + w["val"] * n_val
            + w["jac"] * D[f"{m}|n_jac"])  # fmt: skip


def pct(k: int, n: int) -> str:
    lo, hi = S.wilson(k, n)
    return f"{100 * k / n:.0f}% ({100 * lo:.0f}-{100 * hi:.0f})" if n else "n/a"


def count_ci(k: int, n: int) -> str:
    """``k of n`` with the 95% Wilson interval in percent (one decimal, for small shares)."""
    if not n:
        return "none"
    lo, hi = S.wilson(k, n)
    return f"{k} of {n} ({100 * lo:.1f}-{100 * hi:.1f}%)"


def q(x: np.ndarray, p: float) -> float:
    x = x[np.isfinite(x)]
    return float(np.percentile(x, p)) if len(x) else math.nan


def fmt_n(x: float) -> str:
    return "n/a" if not np.isfinite(x) else (f"{x:.0f}" if x >= 100 else f"{x:.3g}")


def summarize(tag: str, kind: str) -> None:
    D = load(tag)
    methods = [m for m in METHODS if f"{m}|err" in D]
    tm_path = OUT / "timing" / f"{kind}.json"
    if not tm_path.exists():
        tm_path = OUT / "timing" / "e0.json"
    tm_json = json.loads(tm_path.read_text())
    timing = tm_json["ms_per_evaluation"]
    f5t = f5_timing(tm_path.stem, tm_json["voxels"])
    w = weights(timing, f5t)
    w_v1 = weights(timing)
    # f3 single-worker time per run from the timing run
    f3_se = math.nan
    tz = OUT / "cache" / "timing_e0"
    if tz.exists():
        T = load("timing_e0")
        if "f3|secs" in T:
            f3_se = float(np.median(T["f3|secs"])) * 1e3 / timing["strict"]
    se = {m: strict_equiv(D, m, w, f3_se) for m in methods}
    se_v1 = {m: strict_equiv(D, m, w_v1, f3_se) for m in methods}
    raw = {m: (D[f"{m}|n_strict"] + D[f"{m}|n_screen"] + D[f"{m}|n_smooth"]) if m != "f3"
           else np.full(len(D["voxel"]), math.nan) for m in methods}  # fmt: skip
    sets: Dict[str, np.ndarray] = {}
    trows = [
        ("strict cost", "all methods", timing["strict"]),
        ("screen (pixel radius 3), per orientation, batches of 64", "v1 weights only",
         timing["screen"]),
        ("E4 value, one rotation per call", "F4", timing["e4_value"]),
        ("E4 value and Jacobian, one rotation per call", "F4", timing["e4_value_and_jacobian"]),
    ]  # fmt: skip
    if f5t:
        trows += [
            ("screen, per orientation, a shell at a time (chunks of up to 256)", "F5",
             f5t["screen_ms_per_orientation"]),
            ("E4 value, per orientation, 32 rotations per call", "F5",
             f5t["smooth_ms_per_orientation"]),
        ]  # fmt: skip
    tables_timing = doc_tables.markdown_table(
        [
            {
                "kind of evaluation": k,
                "used for": u,
                "milliseconds": v,
                "times a strict evaluation": v / timing["strict"],
            }
            for k, u, v in trows
        ],
        ["kind of evaluation", "used for", "milliseconds", "times a strict evaluation"],
        {"milliseconds": ".4f", "times a strict evaluation": ".3f"},
    )
    tables_timing_voxels = ""
    if (OUT / "timing" / f"{tm_path.stem}_f5.json").exists():
        pv = json.loads((OUT / "timing" / f"{tm_path.stem}_f5.json").read_text())["per_voxel"]
        tables_timing_voxels = doc_tables.markdown_table(
            [
                {
                    "voxel": int(v),
                    "weights from it": "yes" if int(v) == tm_json["voxels"][0] else "",
                    "F4 E4 value (ms)": r["f4_e4_value_ms"],
                    "F4 E4 value and Jacobian (ms)": r["f4_e4_value_and_jacobian_ms"],
                    "F5 screen (ms per orientation)": r["screen_ms_per_orientation"],
                    "F5 E4 (ms per orientation)": r["smooth_ms_per_orientation"],
                    "strict inside F5 (ms)": r["strict_ms_per_orientation"],
                }
                for v, r in pv.items()
            ],
            ["voxel", "weights from it", "F4 E4 value (ms)", "F4 E4 value and Jacobian (ms)",
             "F5 screen (ms per orientation)", "F5 E4 (ms per orientation)",
             "strict inside F5 (ms)"],
            {"F4 E4 value (ms)": ".3f", "F4 E4 value and Jacobian (ms)": ".2f",
             "F5 screen (ms per orientation)": ".4f", "F5 E4 (ms per orientation)": ".3f",
             "strict inside F5 (ms)": ".3f"},
        )  # fmt: skip
    if kind == "e0":
        vs = json.loads((OUT / "voxel_sets.json").read_text())
        sets = {
            "near-axis": np.isin(D["voxel"], vs["near"]),
            "far-from-axis": np.isin(D["voxel"], vs["far"]),
        }
    allm = np.ones(len(D["voxel"]), dtype=bool)
    tables: Dict[str, str] = {}
    summ: Dict[str, Any] = {"tag": tag, "kind": kind, "n_runs": int(len(allm)),
                            "weights_vs_strict": w, "weights_vs_strict_v1": w_v1,
                            "f3_strict_equivalent": f3_se}  # fmt: skip

    def cell(m: str, oi: int, mask: np.ndarray, thr: float) -> Tuple[int, int]:
        sel = mask & (D["oi"] == oi)
        return int((D[f"{m}|err"][sel] < thr - TOL).sum()), int(sel.sum())

    prefix = "evalF_" + kind + "_"
    tables[prefix + "timing"] = tables_timing
    if tables_timing_voxels:
        tables[prefix + "timing_voxels"] = tables_timing_voxels
    present = [oi for oi in range(len(OFFSETS)) if (D["oi"] == oi).any()]
    ocols = [f"{OFFSETS[oi]:g} deg" for oi in present]
    # ---- capture matrices
    for name, thr in (("cap1", 1.0), ("cap01", 0.1), ("cap005", 0.05)):
        rows = []
        for m in methods:
            r: Dict[str, Any] = {"method": LABEL[m]}
            for oi, c in zip(present, ocols):
                k, n = cell(m, oi, allm, thr)
                r[c] = pct(k, n)
                summ.setdefault(name, {}).setdefault(m, {})[f"{OFFSETS[oi]:g}"] = [k, n]
            rows.append(r)
        tables[prefix + name] = doc_tables.markdown_table(rows, ["method"] + ocols)
    # ---- median final error
    rows = []
    for m in methods:
        r = {"method": LABEL[m]}
        for oi, c in zip(present, ocols):
            r[c] = float(np.median(D[f"{m}|err"][D["oi"] == oi]))
            summ.setdefault("median_err", {}).setdefault(m, {})[f"{OFFSETS[oi]:g}"] = r[c]
        rows.append(r)
    tables[prefix + "err"] = doc_tables.markdown_table(
        rows, ["method"] + ocols, {c: ".4f" for c in ocols}
    )
    # ---- evaluations to stop (median and 90th percentile) and to capture
    for name, arr, label in (("stop_raw", raw, "raw"), ("stop_se", se, "strict-equivalent")):
        rows = []
        for m in methods:
            r = {"method": LABEL[m]}
            for oi, c in zip(present, ocols):
                x = arr[m][D["oi"] == oi]
                r[c] = f"{fmt_n(q(x, 50))} ({fmt_n(q(x, 90))})"
                summ.setdefault(name, {}).setdefault(m, {})[f"{OFFSETS[oi]:g}"] = [
                    q(x, 50),
                    q(x, 90),
                ]
            rows.append(r)
        tables[prefix + name] = doc_tables.markdown_table(rows, ["method"] + ocols)
    for name, key, thr in (("tocap1", "cap1", 1.0), ("tocap005", "cap005", 0.05)):
        rows = []
        for m in methods:
            r = {"method": LABEL[m]}
            for oi, c in zip(present, ocols):
                x = D[f"{m}|{key}"][D["oi"] == oi]
                got = np.isfinite(x)
                r[c] = (
                    f"{fmt_n(q(x[got], 50))} ({fmt_n(q(x[got], 90))}) n={int(got.sum())}"
                    if got.any()
                    else "n/a"
                )
                summ.setdefault("to_" + key, {}).setdefault(m, {})[f"{OFFSETS[oi]:g}"] = [
                    q(x[got], 50),
                    q(x[got], 90),
                    int(got.sum()),
                ]
            rows.append(r)
        tables[prefix + name] = doc_tables.markdown_table(rows, ["method"] + ocols)
    # ---- wall time per run, contended (main run) and single-worker (timing run)
    rows = []
    for m in methods:
        r = {"method": LABEL[m]}
        for oi, c in zip(present, ocols):
            r[c] = float(np.median(D[f"{m}|secs"][D["oi"] == oi]))
        rows.append(r)
    tables[prefix + "secs"] = doc_tables.markdown_table(
        rows, ["method"] + ocols, {c: ".3f" for c in ocols}
    )
    if (OUT / "cache" / f"timing_{kind}").exists():
        T = load(f"timing_{kind}")
        rows = []
        for m in [x for x in methods if f"{x}|secs" in T]:
            r = {"method": LABEL[m]}
            for oi, c in zip(present, ocols):
                sel = T["oi"] == oi
                r[c] = float(np.median(T[f"{m}|secs"][sel])) if sel.any() else None
            rows.append(r)
        tables[prefix + "secs_single"] = doc_tables.markdown_table(
            rows, ["method"] + ocols, {c: ".3f" for c in ocols}
        )
        summ["secs_single"] = {
            m: {
                f"{OFFSETS[oi]:g}": float(np.median(T[f"{m}|secs"][T["oi"] == oi]))
                for oi in present
                if (T["oi"] == oi).any()
            }
            for m in methods
            if f"{m}|secs" in T
        }
        # strict-equivalent against single-worker wall time (both in strict evaluations) on the
        # timing run's cases: the median per-run ratio, with the v1 and the current weights
        wall = {m: T[f"{m}|secs"] * 1e3 / timing["strict"] for m in methods if f"{m}|secs" in T}
        rows = []
        for m in [x for x in methods if x in wall and x != "f3"]:
            r = {"method": LABEL[m]}
            for oi, c in zip(present, ocols):
                sel = T["oi"] == oi
                if not sel.any():
                    continue
                r1 = float(np.median(strict_equiv(T, m, w_v1, f3_se)[sel] / wall[m][sel]))
                r2 = float(np.median(strict_equiv(T, m, w, f3_se)[sel] / wall[m][sel]))
                r[c] = f"{r1:.2f} -> {r2:.2f}"
                summ.setdefault("se_over_wall", {}).setdefault(m, {})[f"{OFFSETS[oi]:g}"] = dict(
                    v1=r1, current=r2, n=int(sel.sum()))  # fmt: skip
            rows.append(r)
        cols = [c for oi, c in zip(present, ocols) if (T["oi"] == oi).any()]
        tables[prefix + "se_vs_wall"] = doc_tables.markdown_table(rows, ["method"] + cols)
    # ---- the method's own acceptance test (strict cost <= 0.8) against the truth
    rows = []
    for m in [x for x in ("f4", "f5", "f54") if x in methods]:
        r = {"method": LABEL[m]}
        for oi, c in zip(present, ocols):
            sel = D["oi"] == oi
            acc = D[f"{m}|flag"][sel] > 0.5
            ok = D[f"{m}|err"][sel] < 1.0 - TOL
            r[c] = (
                f"{int(acc.sum())} accepted, {int((acc & ~ok).sum())} of them wrong; "
                f"{int((~acc & ok).sum())} of {int((~acc).sum())} rejected were right"
            )
            summ.setdefault("acceptance", {}).setdefault(m, {})[f"{OFFSETS[oi]:g}"] = dict(
                n=int(sel.sum()), accepted=int(acc.sum()), accepted_wrong=int((acc & ~ok).sum()),
                rejected=int((~acc).sum()), rejected_right=int((~acc & ok).sum()))  # fmt: skip
        rows.append(r)
    tables[prefix + "acc"] = doc_tables.markdown_table(rows, ["method"] + ocols)
    # ---- near / far axis
    if sets:
        for thr, tn, tname in (
            (1.0, "1 deg", "nearfar_cap1"),
            (0.05, "0.05 deg", "nearfar_cap005"),
        ):
            rows = []
            for m in methods:
                r = {"method": LABEL[m]}
                for sname, mask in sets.items():
                    for oi in (3, 4, 5):
                        k, n = cell(m, oi, mask, thr)
                        r[f"{sname}, {OFFSETS[oi]:g} deg"] = pct(k, n)
                        summ.setdefault("nearfar", {}).setdefault(m, {})[
                            f"{sname}|{OFFSETS[oi]:g}|{thr:g}"
                        ] = [k, n]
                rows.append(r)
            tables[prefix + tname] = doc_tables.markdown_table(
                rows, ["method"] + [c for c in rows[0] if c != "method"]
            )

    # ---- choice (the carry rule is applied to the methods that existed at the decision point;
    # a second table adds F4 with the 5-frame window, which came after)
    def rank_table(
        cand: List[str], cost: Dict[str, np.ndarray], scale: Optional[Dict[str, Dict[float, float]]]
    ) -> Tuple[Dict[str, Dict[float, float]], Dict[str, int], List[str]]:
        score: Dict[str, Dict[float, float]] = {m: {} for m in cand}
        rank_sum = {m: 0 for m in cand}
        for oi in (3, 5):
            sc_ = {}
            sel = D["oi"] == oi
            for m in cand:
                rate = float((D[f"{m}|err"][sel] < 1.0 - TOL).mean())
                f = scale[m][OFFSETS[oi]] if scale else 1.0
                sc_[m] = 1000.0 * rate / (float(np.median(cost[m][sel])) * f)
                score[m][OFFSETS[oi]] = sc_[m]
            for m in cand:
                rank_sum[m] += 1 + sum(sc_[o] > sc_[m] for o in cand)  # ties share the better rank
        order = sorted(cand, key=lambda m: (rank_sum[m], -np.mean(list(score[m].values()))))
        return score, rank_sum, order

    # wall-time scaling: the median ratio of single-worker wall time to strict-equivalent on the
    # timing run's cases (current weights), per method and offset
    wall_scale: Optional[Dict[str, Dict[float, float]]] = None
    if "se_over_wall" in summ and all(
        m in summ["se_over_wall"] for m in CARRY_POOL if m in methods
    ):
        wall_scale = {
            m: {float(k): 1.0 / v["current"] for k, v in d.items()}
            for m, d in summ["se_over_wall"].items()
        }
    variants = [("choice", CARRY_POOL, se, None), ("choice_with_f4w", FULL_OK, se, None),
                ("choice_v1", CARRY_POOL, se_v1, None)]  # fmt: skip
    if wall_scale is not None:
        variants.append(("choice_wall", CARRY_POOL, se, wall_scale))
    for tname, pool, cost, scale in variants:
        cand = [m for m in methods if m in pool]
        if "f4w" in pool and "f4w" not in methods:
            continue
        score, rank_sum, order = rank_table(cand, cost, scale)
        carry = order[:2] if tname != "choice_with_f4w" else []
        rows = []
        for m in order:
            rows.append(
                {
                    "method": LABEL[m],
                    "score 2 deg": score[m][2.0],
                    "score 5 deg": score[m][5.0],
                    "rank sum": rank_sum[m],
                    "mean score (tie-break)": float(np.mean(list(score[m].values()))),
                    "carried": "yes" if m in carry else "",
                }
            )
        tables[prefix + tname] = doc_tables.markdown_table(
            rows, ["method", "score 2 deg", "score 5 deg", "rank sum", "mean score (tie-break)",
                   "carried"],
            {"score 2 deg": ".3g", "score 5 deg": ".3g",
             "mean score (tie-break)": ".3g"})  # fmt: skip
        summ[tname] = {"score": {m: {f"{k:g}": v for k, v in score[m].items()} for m in cand},
                       "rank_sum": rank_sum, "carry": carry, "candidates": cand,
                       "order": order}  # fmt: skip
    carry = summ["choice"]["carry"]
    # ---- policy: F4 first, F5 then F4 only where F4's own acceptance test rejects it
    for f4m, pname in (("f4", "policy"), ("f4w", "policy_f4w")):
        if f4m not in methods or "f54" not in methods:
            continue
        rej = D[f"{f4m}|flag"] < 0.5
        err_p = np.where(rej, D["f54|err"], D[f"{f4m}|err"])
        se_p = se[f4m] + rej * se["f54"]
        rows = []
        for oi, c in zip(present, ocols):
            sel = D["oi"] == oi
            r = {
                "offset": c,
                "runs": int(sel.sum()),
                f"{SHORT[f4m]} rejected, F5+F4 run": int((rej & sel).sum()),
            }
            for thr, nm in ((1.0, "1 deg"), (0.1, "0.1 deg"), (0.05, "0.05 deg")):
                k = int((err_p[sel] < thr - TOL).sum())
                r["within " + nm] = pct(k, int(sel.sum()))
            r["median evals (strict-equivalent)"] = fmt_n(q(se_p[sel], 50))
            r["90th percentile"] = fmt_n(q(se_p[sel], 90))
            summ.setdefault("policy_" + f4m, {})[f"{OFFSETS[oi]:g}"] = dict(
                n=int(sel.sum()), fallback=int((rej & sel).sum()),
                within_1=int((err_p[sel] < 1.0 - TOL).sum()),
                within_01=int((err_p[sel] < 0.1 - TOL).sum()),
                within_005=int((err_p[sel] < 0.05 - TOL).sum()),
                median_se=q(se_p[sel], 50), p90_se=q(se_p[sel], 90))  # fmt: skip
            rows.append(r)
        tables[prefix + pname] = doc_tables.markdown_table(rows, list(rows[0]))
    # ---- per voxel: runs within 1 deg out of 20, at each offset of 2 deg or more (runs of one
    # voxel share its images, so the Wilson intervals above, which treat runs as independent,
    # are too narrow when the voxels differ)
    vox = sorted(set(int(v) for v in D["voxel"]))
    pv_m = [m for m in ("f1", "f2_1", "f2_2", "f4", "f4w", "f5", "f54") if m in methods]
    for oi in [x for x in (3, 4, 5) if x in present]:
        for m in pv_m:
            ks = []
            for v in vox:
                sel = (D["oi"] == oi) & (D["voxel"] == v)
                ks.append([int((D[f"{m}|err"][sel] < 1.0 - TOL).sum()), int(sel.sum())])
            summ.setdefault("per_voxel_cap1", {}).setdefault(m, {})[f"{OFFSETS[oi]:g}"] = dict(
                zip([str(v) for v in vox], ks))  # fmt: skip
    if kind == "full":
        rows = []
        for m in pv_m:
            r = {"method": LABEL[m]}
            for v in vox:
                k, n = summ["per_voxel_cap1"][m]["5"][str(v)]
                r[f"voxel {v}"] = f"{k} of {n}"
            rows.append(r)
        tables[prefix + "pervoxel_cap1_5"] = doc_tables.markdown_table(rows, list(rows[0]))
    else:
        rows = []
        for m in pv_m:
            r = {"method": LABEL[m]}
            for oi in [x for x in (3, 4, 5) if x in present]:
                d = summ["per_voxel_cap1"][m][f"{OFFSETS[oi]:g}"]
                lost = [v for v, (k, n) in d.items() if k < n]
                r[f"{OFFSETS[oi]:g} deg: voxels with a lost run"] = f"{len(lost)} of {len(d)}"
                r[f"{OFFSETS[oi]:g} deg: fewest within 1 deg"] = (
                    f"{min(k for k, n in d.values())} of {list(d.values())[0][1]}")  # fmt: skip
            rows.append(r)
        tables[prefix + "pervoxel_cap1"] = doc_tables.markdown_table(rows, list(rows[0]))
    # ---- acceptance margins, all offsets pooled (and 5 deg), with Wilson intervals (runs counted
    # as independent; a run is wrong when it ends 1 deg or more out)
    rows = []
    for m in [x for x in ("f4", "f4w") if x in methods]:
        for lab, sel in (("all offsets", np.ones(len(D["voxel"]), dtype=bool)),
                         ("5 deg", D["oi"] == 5)):  # fmt: skip
            acc = (D[f"{m}|flag"] > 0.5) & sel
            rej = (D[f"{m}|flag"] < 0.5) & sel
            ok = D[f"{m}|err"] < 1.0 - TOL
            c, e = D[f"{m}|cost"], D[f"{m}|err"]
            na, nr = int(acc.sum()), int(rej.sum())
            r = {
                "method": LABEL[m],
                "runs": lab,
                "accepted but wrong": count_ci(int((acc & ~ok).sum()), na),
                "rejected but right": count_ci(int((rej & ok).sum()), nr),
                "largest accepted strict cost": float(c[acc].max()) if na else math.nan,
                "smallest rejected strict cost": float(c[rej].min()) if nr else math.nan,
                "largest error accepted (deg)": float(e[acc].max()) if na else math.nan,
                "smallest error rejected (deg)": float(e[rej].min()) if nr else math.nan,
            }  # fmt: skip
            summ.setdefault("acceptance_margin", {}).setdefault(m, {})[lab] = dict(
                n_accepted=na, accepted_wrong=int((acc & ~ok).sum()), n_rejected=nr,
                rejected_right=int((rej & ok).sum()),
                max_cost_accepted=r["largest accepted strict cost"],
                min_cost_rejected=r["smallest rejected strict cost"],
                max_err_accepted=r["largest error accepted (deg)"],
                min_err_rejected=r["smallest error rejected (deg)"])  # fmt: skip
            rows.append(r)
    if rows:
        tables[prefix + "acc_margin"] = doc_tables.markdown_table(
            rows, list(rows[0]),
            {"largest accepted strict cost": ".3f", "smallest rejected strict cost": ".3f",
             "largest error accepted (deg)": ".3f",
             "smallest error rejected (deg)": ".2f"})  # fmt: skip
    # ---- what F4's evaluations are made of: Jacobians, E4 values, strict (median per run)
    for m in [x for x in ("f4", "f4w") if x in methods]:
        for oi in present:
            sel = D["oi"] == oi
            nj = D[f"{m}|n_jac"][sel]
            nv = D[f"{m}|n_smooth"][sel] - nj
            share = w["jac"] * nj / se[m][sel]
            summ.setdefault("composition", {}).setdefault(m, {})[f"{OFFSETS[oi]:g}"] = dict(
                median_jac=float(np.median(nj)), median_val=float(np.median(nv)),
                median_strict=float(np.median(D[f"{m}|n_strict"][sel])),
                median_raw=float(np.median(raw[m][sel])),
                median_share_se_from_jac=float(np.median(share)))  # fmt: skip
    if "f4" in methods and "f54" in methods:  # the policy's cost with the v1 weights
        rej = D["f4|flag"] < 0.5
        se_p1 = se_v1["f4"] + rej * se_v1["f54"]
        summ["policy_f4_v1"] = {
            f"{OFFSETS[oi]:g}": dict(
                median_se=q(se_p1[D["oi"] == oi], 50), p90_se=q(se_p1[D["oi"] == oi], 90)
            )  # fmt: skip
            for oi in present
        }
    summ["stop_se_v1"] = {
        m: {f"{OFFSETS[oi]:g}": [q(se_v1[m][D["oi"] == oi], p) for p in (50, 90)] for oi in present}
        for m in ("f5", "f54")
        if m in methods
    }
    doc_tables.write_tables(OUT / f"tables_{tag}.md", tables)
    (OUT / f"summary_{tag}.json").write_text(json.dumps(summ, indent=1, default=float) + "\n")
    figures(D, methods, se, kind, tag, sets)
    print("carry:", carry)


def figures(
    D: Dict[str, Any],
    methods: Sequence[str],
    se: Dict[str, np.ndarray],
    kind: str,
    tag: str,
    sets: Dict[str, np.ndarray],
) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fdir = OUT / "figures"
    fdir.mkdir(exist_ok=True)
    present = [oi for oi in range(len(OFFSETS)) if (D["oi"] == oi).any()]
    nc = 3
    nr = int(math.ceil(len(present) / nc))
    err = lambda m, sel: D[f"{m}|err"][sel]  # noqa: E731
    rate = lambda m, sel, thr: float((err(m, sel) < thr - TOL).mean())  # noqa: E731

    # ---- 1. capture rate against evaluations, one panel per offset
    fig, axes = plt.subplots(nr, nc, figsize=(14, 3.6 * nr + 0.6), sharey=True, squeeze=False)
    for ax in axes.ravel()[len(present) :]:
        ax.axis("off")
    for k, (oi, ax) in enumerate(zip(present, axes.ravel())):
        sel = D["oi"] == oi
        for m in methods:
            if m == "f3":
                continue
            y = (D[f"{m}|curve"][sel] < 1.0 - TOL).mean(axis=0)
            ax.plot(
                CKPT,
                y,
                marker="o",
                ms=3,
                color=COLOR[m],
                label=SHORT[m],
                ls="--" if m == "f5" else "-",
                lw=2.2 if m in ("f4", "f54") else 1.4,
            )
        ax.set_xscale("log")
        ax.set_title(f"start {OFFSETS[oi]:g} deg out (n = {int(sel.sum())})")
        ax.set_ylim(-0.02, 1.02)
        ax.grid(alpha=0.3)
        if k % nc == 0:
            ax.set_ylabel("share of runs within 1 deg of the truth")
    axes[0, 0].legend(fontsize=8, loc="lower right")
    fig.suptitle(f"Capture rate against evaluations ({kind}). F3 is not drawn: it makes no cost "
                 "evaluations. F5 lies under F5+F4.")  # fmt: skip
    fig.supxlabel(
        "cost evaluations so far, raw count (strict + screen + smooth, each counted once; "
        "not time: one F4 Jacobian costs about 26 strict evaluations, one screen about 0.04)"
    )
    fig.tight_layout()
    fig.savefig(fdir / f"capture_vs_evals_{tag}.png", dpi=130)
    plt.close(fig)

    # ---- 2. capture against the start offset (1 deg and 0.05 deg)
    xs = [OFFSETS[oi] for oi in present]
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2), sharey=True)
    for ax, thr in zip(axes, (1.0, 0.05)):
        for m in methods:
            y = [rate(m, D["oi"] == oi, thr) for oi in present]
            ax.plot(xs, y, marker="o", color=COLOR[m], label=SHORT[m])
        ax.set_xscale("log")
        ax.set_xlabel("start offset from the truth (deg)")
        ax.set_title(f"final error under {thr:g} deg")
        ax.grid(alpha=0.3)
    axes[0].set_ylabel("share of runs")
    axes[0].legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(fdir / f"capture_vs_offset_{tag}.png", dpi=130)
    plt.close(fig)

    # ---- 3. efficiency: capture rate against median strict-equivalent evaluations to stop
    eff = [oi for oi in (2, 3, 5) if oi in present]
    fig, axes = plt.subplots(1, len(eff), figsize=(5 * len(eff), 4.4), sharey=True, squeeze=False)
    for ax, oi in zip(axes[0], eff):
        sel = D["oi"] == oi
        for m in methods:
            x = float(np.median(se[m][sel]))
            ax.scatter([x], [rate(m, sel, 1.0)], s=70, color=COLOR[m], label=SHORT[m], marker="o")
            ax.scatter([x], [rate(m, sel, 0.05)], s=70, facecolors="none", edgecolors=COLOR[m],
                       marker="s")  # fmt: skip
        ax.set_xscale("log")
        ax.set_title(f"start {OFFSETS[oi]:g} deg out")
        ax.set_xlabel("median evaluations to stop (strict-equivalent, time-weighted)")
        ax.grid(alpha=0.3)
    axes[0, 0].set_ylabel("share of runs (filled: within 1 deg; open square: within 0.05 deg)")
    axes[0, 0].legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(fdir / f"efficiency_{tag}.png", dpi=130)
    plt.close(fig)

    # ---- 4. near / far axis
    if sets:
        fig, axes = plt.subplots(1, 2, figsize=(11, 4.2), sharey=True)
        for ax, (sname, mask) in zip(axes, sets.items()):
            for m in methods:
                y = [rate(m, mask & (D["oi"] == oi), 1.0) for oi in present]
                ax.plot(xs, y, marker="o", color=COLOR[m], label=SHORT[m])
            ax.set_xscale("log")
            ax.set_title(f"{sname} voxels: final error under 1 deg")
            ax.set_xlabel("start offset (deg)")
            ax.grid(alpha=0.3)
        axes[0].set_ylabel("share of runs")
        axes[0].legend(fontsize=8)
        fig.tight_layout()
        fig.savefig(fdir / f"nearfar_{tag}.png", dpi=130)
        plt.close(fig)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("kind", choices=["e0", "full", "merge"])
    ap.add_argument("--tag", default="")
    a = ap.parse_args()
    if a.kind == "merge":  # every tables_<tag>.md into doc_tables_all.md for sync_doc_tables.py
        merged: Dict[str, str] = {}
        for f in sorted(OUT.glob("tables_*.md")):
            merged.update(doc_tables.read_tables(f))
        doc_tables.write_tables(OUT / "doc_tables_all.md", merged)
        print(f"{len(merged)} tables -> doc_tables_all.md")
        return
    summarize(a.tag or a.kind, a.kind)


if __name__ == "__main__":
    main()
