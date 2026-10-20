"""Task F summary: capture rate, error and evaluations per method, start offset and voxel set.

Reads ``benchmarks/eval_reduction/F/cache/<tag>/v<voxel>_o<offset>.npz`` (outward_bench.py), writes
``summary_<tag>.json``, ``tables_<tag>.md`` (marker-delimited blocks for the findings doc) and the
figures. A run is *captured* at 1 deg (0.1, 0.05 deg) when the method's final orientation is within
that angle of the truth (angle of R R_true^T, no symmetry reduction, offsets are at most 5 deg).

Evaluation counts. Raw: strict + screen + smooth evaluations, each counted once (the convention of
the evaluation ledger). Strict-equivalent: each kind weighted by its single-worker time relative to
a strict evaluation (``timing/e0.json``: strict 0.41 ms, screen 0.023 ms, E4 value 0.66 ms, E4
value-and-Jacobian 10.8 ms). The windows Gauss-Newton (f3) makes no cost evaluations; its entry is
its measured single-worker time per run divided by the strict evaluation time.
Choice rule (Task F decision point): among the methods that can run on the full sample, the capture
rate (1 deg) per strict-equivalent evaluation to stop, at 2 and at 5 deg; ranks at the two offsets
added, the two lowest rank sums are carried, ties by the larger mean score.

  uv run python scripts/eval_reduction/outward_summary.py e0
  uv run python scripts/eval_reduction/outward_summary.py full
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


def weights(timing: Dict[str, float]) -> Dict[str, float]:
    s = timing["strict"]
    return dict(screen=timing["screen"] / s, val=timing["e4_value"] / s,
                jac=timing["e4_value_and_jacobian"] / s)  # fmt: skip


def strict_equiv(D: Dict[str, Any], m: str, w: Dict[str, float], f3_se: float) -> np.ndarray:
    if m == "f3":
        return np.full(len(D["voxel"]), f3_se)
    n_val = D[f"{m}|n_smooth"] - D[f"{m}|n_jac"]
    return (D[f"{m}|n_strict"] + w["screen"] * D[f"{m}|n_screen"] + w["val"] * n_val
            + w["jac"] * D[f"{m}|n_jac"])  # fmt: skip


def pct(k: int, n: int) -> str:
    lo, hi = S.wilson(k, n)
    return f"{100 * k / n:.0f}% ({100 * lo:.0f}-{100 * hi:.0f})" if n else "n/a"


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
    timing = json.loads(tm_path.read_text())["ms_per_evaluation"]
    w = weights(timing)
    # f3 single-worker time per run from the timing run
    f3_se = math.nan
    tz = OUT / "cache" / "timing_e0"
    if tz.exists():
        T = load("timing_e0")
        if "f3|secs" in T:
            f3_se = float(np.median(T["f3|secs"])) * 1e3 / timing["strict"]
    se = {m: strict_equiv(D, m, w, f3_se) for m in methods}
    raw = {m: (D[f"{m}|n_strict"] + D[f"{m}|n_screen"] + D[f"{m}|n_smooth"]) if m != "f3"
           else np.full(len(D["voxel"]), math.nan) for m in methods}  # fmt: skip
    sets: Dict[str, np.ndarray] = {}
    if kind == "e0":
        vs = json.loads((OUT / "voxel_sets.json").read_text())
        sets = {
            "near-axis": np.isin(D["voxel"], vs["near"]),
            "far-from-axis": np.isin(D["voxel"], vs["far"]),
        }
    allm = np.ones(len(D["voxel"]), dtype=bool)
    tables: Dict[str, str] = {}
    summ: Dict[str, Any] = {"tag": tag, "kind": kind, "n_runs": int(len(allm)),
                            "weights_vs_strict": w, "f3_strict_equivalent": f3_se}  # fmt: skip

    def cell(m: str, oi: int, mask: np.ndarray, thr: float) -> Tuple[int, int]:
        sel = mask & (D["oi"] == oi)
        return int((D[f"{m}|err"][sel] < thr - TOL).sum()), int(sel.sum())

    prefix = "evalF_" + kind + "_"
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
        rows = []
        for m in methods:
            r = {"method": LABEL[m]}
            for sname, mask in sets.items():
                for thr, tn in ((1.0, "1 deg"), (0.05, "0.05 deg")):
                    for oi in (3, 4, 5):
                        k, n = cell(m, oi, mask, thr)
                        r[f"{sname}, {OFFSETS[oi]:g} deg, <{tn}"] = pct(k, n)
                        summ.setdefault("nearfar", {}).setdefault(m, {})[
                            f"{sname}|{OFFSETS[oi]:g}|{thr:g}"
                        ] = [k, n]
            rows.append(r)
        cols = ["method"] + [c for c in rows[0] if c != "method"]
        tables[prefix + "nearfar"] = doc_tables.markdown_table(rows, cols)

    # ---- choice (the carry rule is applied to the methods that existed at the decision point;
    # a second table adds F4 with the 5-frame window, which came after)
    def rank_table(
        cand: List[str],
    ) -> Tuple[Dict[str, Dict[float, float]], Dict[str, int], List[str]]:
        score: Dict[str, Dict[float, float]] = {m: {} for m in cand}
        rank_sum = {m: 0 for m in cand}
        for oi in (3, 5):
            sc_ = {}
            sel = D["oi"] == oi
            for m in cand:
                rate = float((D[f"{m}|err"][sel] < 1.0 - TOL).mean())
                sc_[m] = 1000.0 * rate / float(np.median(se[m][sel]))
                score[m][OFFSETS[oi]] = sc_[m]
            for m in cand:
                rank_sum[m] += 1 + sum(sc_[o] > sc_[m] for o in cand)  # ties share the better rank
        order = sorted(cand, key=lambda m: (rank_sum[m], -np.mean(list(score[m].values()))))
        return score, rank_sum, order

    for tname, pool in (("choice", CARRY_POOL), ("choice_with_f4w", FULL_OK)):
        cand = [m for m in methods if m in pool]
        if "f4w" in pool and "f4w" not in methods:
            continue
        score, rank_sum, order = rank_table(cand)
        carry = order[:2] if tname == "choice" else []
        rows = []
        for m in order:
            rows.append(
                {
                    "method": LABEL[m],
                    "score 2 deg": score[m][2.0],
                    "score 5 deg": score[m][5.0],
                    "rank sum": rank_sum[m],
                    "carried": "yes" if m in carry else "",
                }
            )
        tables[prefix + tname] = doc_tables.markdown_table(
            rows, ["method", "score 2 deg", "score 5 deg", "rank sum", "carried"],
            {"score 2 deg": ".3g", "score 5 deg": ".3g"})  # fmt: skip
        summ[tname] = {"score": {m: {f"{k:g}": v for k, v in score[m].items()} for m in cand},
                       "rank_sum": rank_sum, "carry": carry, "candidates": cand}  # fmt: skip
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
    fig.supxlabel("cost evaluations so far (strict + screen + smooth, each counted once)")
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
    ap.add_argument("kind", choices=["e0", "full"])
    ap.add_argument("--tag", default="")
    a = ap.parse_args()
    summarize(a.tag or a.kind, a.kind)


if __name__ == "__main__":
    main()
