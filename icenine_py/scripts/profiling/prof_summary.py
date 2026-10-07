#!/usr/bin/env python3
"""Summary of the run-time profiling (Task 3): timing tables of U0 / U1 / U2, stage breakdowns,
paired accuracy, the contention factor, profiler highlights and the "helps" verdicts, all from the
stored per-run records (scripts/profiling/cache) plus the committed accuracy of Tasks 1-2
(benchmarks/coarse_proxy/summary.json, benchmarks/nn_hybrid/*_raw.npz).

  uv run python scripts/profiling/prof_summary.py
  (writes benchmarks/profiling/summary.{txt,json}, u0_runs.json.gz and seeded_runs.json.gz)
"""

import gzip
import io
import json
import sys
from collections import defaultdict
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import prof_common as PC  # noqa: E402

CACHE, OUT = PC.CACHE, PC.OUT
VAR_NAME = {"clean": "clean", "all": "realistic"}
VARIANTS = ["clean", "all"]
RADII = [0.05, 0.1, 0.25, 0.5, 0.75, 1.0, 1.5, 2.0, 3.0, 5.0]
BANDS = {  # use case -> radius indices
    "U1 (r 0.05, 0.1)": [0, 1],
    "U2 (r 0.5 .. 3)": [3, 4, 5, 6, 7, 8],
    "U2 stress (r 5)": [9],
}
U0_PIPES = ["baseline", "F1", "F1b", "e2", "proxy", "proxy+F1"]
U0_ROW = {  # Task 2 end-to-end row names
    "baseline": "E0 baseline",
    "F1": "F1 (Sigma<=29)",
    "F1b": "F1b",
    "e2": "E2 rerank (GBT, full pass)",
    "proxy": "proxy p_i",
    "proxy+F1": "proxy p_i + F1",
}
U0_ROW3 = {
    "baseline": "E0 baseline (3 seeds)",
    "F1": "F1 (Sigma<=29, 3 seeds)",
    "proxy": "proxy p_i (3 seeds)",
    "proxy+F1": "proxy p_i + F1 (3 seeds)",
}
TOL_MED = 0.005  # "median error within +0.005 deg"
FALLBACK = "H0+fallback"
HG_RI = {6, 7, 8, 9}  # radii at which HG is run (1.5, 2, 3, 5 deg)


def f3(x: float) -> str:
    return "nan" if x is None or not np.isfinite(x) else f"{x:.3g}"


def fr(x: float) -> str:
    return "nan" if x is None or not np.isfinite(x) else f"{100 * x:.3g}%"


def load_dir(d: Path) -> List[Dict[str, Any]]:
    out: List[Dict[str, Any]] = []
    for f in sorted(d.glob("v*.json")):
        out += json.loads(f.read_text())["records"]
    return out


def case_median(recs: List[Dict[str, Any]], key: Tuple[str, ...], value: Any) -> Dict[Tuple, float]:
    """Median over repetitions of value(r), per case key."""
    g: Dict[Tuple, List[float]] = defaultdict(list)
    for r in recs:
        g[tuple(r[k] for k in key)].append(value(r))
    return {k: float(np.median(v)) for k, v in g.items()}


def rep_spread(recs: List[Dict[str, Any]], key: Tuple[str, ...], value: Any) -> Dict[str, float]:
    """Run-to-run spread on the repeated cases: (max - min) / median of the repeats."""
    g: Dict[Tuple, List[float]] = defaultdict(list)
    for r in recs:
        g[tuple(r[k] for k in key)].append(value(r))
    s = [(max(v) - min(v)) / np.median(v) for v in g.values() if len(v) >= 2]
    return dict(
        n_cases=len(s),
        median=float(np.median(s)) if s else float("nan"),
        p90=float(np.percentile(s, 90)) if s else float("nan"),
    )


def tstat(x: Sequence[float]) -> Dict[str, float]:
    return (
        PC.summarize_times(list(x))
        if len(x)
        else dict(n=0, median=np.nan, p10=np.nan, p90=np.nan, mean=np.nan)
    )


def wilson(k: int, n: int) -> Tuple[float, float, float]:
    return PC.wilson(k, n)


def tfmt(s: Dict[str, float]) -> str:
    return f"{f3(s['median'])} [{f3(s['p10'])}-{f3(s['p90'])}]"


# ---------------------------------------------------------------------------
# U0
# ---------------------------------------------------------------------------

STAGE_GROUPS = [
    ("discrete search (excl. evaluate)", lambda k: k.startswith("discrete_")),
    ("quick MC (excl. evaluate)", lambda k: k.startswith("quick_mc_")),
    ("FindOptimal (excl. evaluate)", lambda k: k == "find_optimal"),
    ("VarianceMinimizing (excl. evaluate)", lambda k: k == "variance"),
    ("cost evaluations, global", lambda k: k == "evaluate_global"),
    ("cost evaluations, local", lambda k: k == "evaluate_local"),
    ("rank_key proxy: set_image", lambda k: k == "proxy_setimage"),
    ("rank_key proxy: features", lambda k: k == "proxy_features"),
    ("rank_key proxy: prediction", lambda k: k == "proxy_predict"),
]


def u0_tables(recs: List[Dict[str, Any]], t2: Dict[str, Any]) -> Tuple[Dict[str, Any], List[str]]:
    lines: List[str] = []
    res: Dict[str, Any] = {}
    tval = lambda r: r.get("wall_total", r["wall"])  # noqa: E731
    by_var: Dict[str, Dict[str, Dict[str, Any]]] = {}
    for var in VARIANTS:
        rv = [r for r in recs if r["variant"] == var]
        res[var] = {}
        by_var[var] = {}
        lines.append(
            f"-- U0, {VAR_NAME[var]}: {len({r['vidx'] for r in rv})} voxels, seed 0, "
            f"single worker; wall seconds per voxel (median [10-90%]) --"
        )
        lines.append(
            f"  {'pipeline':10s} {'n':>3s} {'time s':>20s} {'mean s':>7s} {'glob ev':>8s} "
            f"{'loc ev':>7s} {'all ev':>8s} {'wrong':>9s} {'med err(right)':>14s} {'=stored':>8s}"
        )
        for p in U0_PIPES:
            rp = [r for r in rv if r["pipe"] == p]
            tm = case_median(rp, ("vidx",), tval)
            st = tstat(list(tm.values()))
            r0 = [r for r in rp if r["rep"] == 0]
            wrong = sum(r["err"] > 1.0 for r in r0)
            right = [r["err"] for r in r0 if r["err"] <= 1.0]
            ge, le, ea = [], [], []
            for r in r0:
                if "evals_rec" in r:
                    ge.append(r["evals_rec"][0])
                    le.append(r["evals_rec"][1])
                    ea.append(r["evals_all"]["global"] + r["evals_all"]["local"])
            if p in ("F1", "proxy+F1"):
                par = "baseline" if p == "F1" else "proxy"
                pe = {(r["vidx"]): r for r in rv if r["pipe"] == par and r["rep"] == 0}
                ge = [pe[r["vidx"]]["evals_rec"][0] for r in r0]
                le = [pe[r["vidx"]]["evals_rec"][1] + r["evals_all"].get("local", 0) for r in r0]
                ea = [
                    sum(pe[r["vidx"]]["evals_all"].values()) + sum(r["evals_all"].values())
                    for r in r0
                ]
            ident = sum(r["identical_to_stored"] for r in rp)
            d = dict(
                time=st,
                wrong=wrong,
                n=len(r0),
                med_err_right=float(np.median(right)) if right else np.nan,
                evals_global=float(np.mean(ge)),
                evals_local=float(np.mean(le)),
                evals_all=float(np.mean(ea)),
                identical=ident,
                n_runs=len(rp),
                per_voxel=tm,
            )
            res[var][p] = d
            by_var[var][p] = d
            lines.append(
                f"  {p:10s} {st['n']:3d} {tfmt(st):>20s} {f3(st['mean']):>7s} "
                f"{f3(np.mean(ge)):>8s} "
                f"{f3(np.mean(le)):>7s} {f3(np.mean(ea)):>8s} {wrong:>4d}/{len(r0):<4d} "
                f"{f3(d['med_err_right']):>14s} "
                f"{ident:>3d}/{len(rp):<4d}"
            )
        sp = rep_spread([r for r in rv], ("vidx", "pipe"), tval)
        if sp["n_cases"]:
            lines.append(
                f"  repeats (3 x on {sp['n_cases'] // len(U0_PIPES)} voxels): per-case "
                f"(max-min)/median of the 3 timings: median {fr(sp['median'])}, "
                f"90th pct {fr(sp['p90'])}"
            )
        else:
            lines.append("  no repeats (realistic)")
        res[var]["repeat_spread"] = sp
        # accuracy of this subset vs the 200-voxel run of Task 2
        lines.append("  accuracy on these voxels vs the 200-voxel runs of Task 2 (wrong rate):")
        for p in U0_PIPES:
            row = _t2_row(t2, var, p)
            d = res[var][p]
            lines.append(
                f"    {p:10s} subset {d['wrong']}/{d['n']} = {fr(d['wrong'] / d['n'])}; "
                f"200-voxel run {fr(row['rate'])} [{fr(row['lo'])}, {fr(row['hi'])}]"
            )
        lines.append("")
    # stage breakdown
    lines.append(
        "-- U0 stage breakdown: mean seconds per run (exclusive times; their sum is the "
        "wall time of the pipeline's runs) --"
    )
    res["stages"] = {}
    for var in VARIANTS:
        for p in ("baseline", "F1", "F1b", "e2", "proxy", "proxy+F1"):
            rp = [r for r in recs if r["variant"] == var and r["pipe"] == p and r["rep"] == 0]
            acc: Dict[str, float] = defaultdict(float)
            for r in rp:
                ex = r["stages"]["exclusive"]
                for name, sel in STAGE_GROUPS:
                    acc[name] += sum(v for k, v in ex.items() if sel(k)) / len(rp)
                if p in ("F1", "proxy+F1"):
                    acc["F1 post hoc: total"] += r["stages"]["inclusive"].get("f1", 0.0) / len(rp)
                acc["wall"] += (
                    tval(r) / len(rp) if p not in ("F1", "proxy+F1") else (r["wall"] / len(rp))
                )  # noqa: E501
            res["stages"][f"{var}/{p}"] = dict(acc)
    for var in VARIANTS:
        lines.append(f"  {VAR_NAME[var]}")
        names = [n for n, _ in STAGE_GROUPS]
        lines.append(
            f"  {'pipeline':10s} "
            + " ".join(f"{n[:22]:>22s}" for n in names[:6])
            + f" {'proxy set/feat/pred':>20s} {'wall':>7s}"
        )
        for p in ("baseline", "F1b", "e2", "proxy", "F1", "proxy+F1"):
            a = res["stages"][f"{var}/{p}"]
            wall = a["wall"]
            cells = " ".join(
                f"{f3(a.get(n, 0.0)):>7s} ({100 * a.get(n, 0.0) / wall:4.1f}%)".rjust(22)
                for n in names[:6]
            )
            px = "/".join(f3(a.get(n, 0.0)) for n in names[6:9])
            lines.append(f"  {p:10s} {cells} {px:>20s} {f3(wall):>7s}")
    lines.append(
        "  (F1 and proxy+F1 rows: the F1 post hoc step alone; its wall = quick MC of "
        "the relatives + FindOptimal of the top relatives)"
    )
    return res, lines


def _t2_row(t2: Dict[str, Any], var: str, pipe: str, three: bool = False) -> Dict[str, Any]:
    if three:
        rows = t2["multi_seed"][f"p_i/{var}"]["rows"]
        name = U0_ROW3[pipe]
    else:
        rows = t2["end_to_end"][var]
        name = U0_ROW[pipe]
    for r in rows:
        if r["name"] == name:
            return r
    raise KeyError(name)


def u0_verdicts(res: Dict[str, Any], t2: Dict[str, Any]) -> Tuple[Dict[str, Any], List[str]]:
    """The plan's U0 rule: helps vs a comparator if the median wall time per voxel is >= 20% lower
    at matched accuracy (wrong rate <= comparator's Wilson upper bound and median error within
    +0.005 deg), or the wrong rate is lower with non-overlapping intervals at <= 10% extra time."""
    out: Dict[str, Any] = {}
    lines = ["-- U0 verdicts (per variant; 'helps' overall only if it holds in both) --"]
    for p in U0_PIPES[1:]:
        for c in ("baseline", "F1", "F1b"):
            if p == c:
                continue
            ok_v: List[bool] = []
            parts = []
            for var in VARIANTS:
                three = p == "proxy+F1" and c in ("baseline", "F1")
                rp, rc = _t2_row(t2, var, p, three), _t2_row(t2, var, c, three)
                tp = res[var][p]["time"]["median"]
                tc = res[var][c]["time"]["median"]
                matched = rp["rate"] <= rc["hi"] and (
                    rp["median_err_right"] <= rc["median_err_right"] + TOL_MED
                )
                saving = 1 - tp / tc
                a = matched and saving >= 0.20
                b = rp["hi"] < rc["lo"] and tp <= 1.10 * tc
                ok = bool(a or b)
                ok_v.append(ok)
                parts.append(
                    f"{VAR_NAME[var]}: time {f3(tp)} vs {f3(tc)} s ({100 * (tp / tc - 1):+.1f}%), "
                    f"wrong {fr(rp['rate'])} [{fr(rp['lo'])}, {fr(rp['hi'])}] vs {fr(rc['rate'])} "
                    f"[{fr(rc['lo'])}, {fr(rc['hi'])}] ({'3-seed' if three else 'seed 0'}), "
                    f"med err {f3(rp['median_err_right'])} vs {f3(rc['median_err_right'])}, "
                    f"matched {matched}, (a) {a}, (b) {b} -> {ok}"
                )
            out[f"{p} vs {c}"] = dict(helps=all(ok_v), per_variant=dict(zip(VARIANTS, ok_v)))
            lines.append(f"  {p} vs {c}: helps = {all(ok_v)}")
            lines += [f"      {s}" for s in parts]
    return out, lines


# ---------------------------------------------------------------------------
# U1 / U2
# ---------------------------------------------------------------------------


def seeded_tables(
    recs: List[Dict[str, Any]], full: Any, u0: Dict[str, Any]
) -> Tuple[Dict[str, Any], List[str]]:
    lines: List[str] = []
    res: Dict[str, Any] = {}
    pipes = ["H0", "H1", "H3", "MCr", "HG"]
    key = ("vidx", "ri", "j", "variant", "pipe")
    tprod = case_median(recs, key, lambda r: r["total_production"])
    tmeas = case_median(recs, key, lambda r: r["total_measured"])
    nprod = case_median(recs, key, lambda r: r["net_production"])
    nmeas = case_median(recs, key, lambda r: r["net_measured"])
    evals = case_median(recs, key, lambda r: r["evals"])
    rep0 = {tuple(r[k] for k in key): r for r in recs if r["rep"] == 0}
    errs = {k: r["err"] for k, r in rep0.items()}
    fb = {k: r["fallback"] for k, r in rep0.items()}
    _STORE.update(tprod=tprod, full=full)

    def sel(pipe: str, var: str, ris: Sequence[int], d: Dict[Tuple, float]) -> List[float]:
        return [v for k, v in d.items() if k[4] == pipe and k[3] == var and k[1] in ris]

    # per-radius table
    lines.append(
        "-- U1 / U2 per radius: production-equivalent wall seconds per case, median "
        "[10-90%] (both variants pooled), 400 cases per pipeline, single worker --"
    )
    lines.append(f"  {'r':>5s} " + " ".join(f"{p:>20s}" for p in pipes))
    res["per_radius"] = {}
    for ri, r in enumerate(RADII):
        row, cells = {}, []
        for p in pipes:
            vals = [v for var in VARIANTS for v in sel(p, var, [ri], tprod)]
            row[p] = tstat(vals)
            cells.append(f"{tfmt(row[p]) if vals else '-':>20s}")
        res["per_radius"][str(r)] = row
        lines.append(f"  {r:5.2f} " + " ".join(cells))
    lines.append("")
    # band tables
    res["bands"] = {}
    for band, ris in BANDS.items():
        res["bands"][band] = {}
        for var in VARIANTS:
            lines.append(f"-- {band}, {VAR_NAME[var]} --")
            lines.append(
                f"  {'pipe':5s} {'n':>4s} {'prod. time s':>20s} {'measured s':>11s} "
                f"{'net prod/meas s':>16s} {'evals':>7s} {'wrong (subset)':>15s} "
                f"{'med err':>8s} {'wrong full run [95%]':>24s} {'med err full':>12s} "
                f"{'fallb':>5s}"
            )
            for p in pipes:
                keys = [k for k in tprod if k[4] == p and k[3] == var and k[1] in ris]
                if not keys:
                    continue
                t = tstat([tprod[k] for k in keys])
                tm = tstat([tmeas[k] for k in keys])
                npd = float(np.median([nprod[k] for k in keys]))
                nmd = float(np.median([nmeas[k] for k in keys]))
                ev = float(np.mean([evals[k] for k in keys]))
                e = np.array([errs[k] for k in keys])
                kw = int((e > 1.0).sum())
                med = float(np.median(e))
                fs_ = full.stats_band(p, ris, VARIANTS.index(var))
                d = dict(
                    time=t,
                    measured=tm,
                    net_production=npd,
                    net_measured=nmd,
                    evals=ev,
                    n=len(keys),
                    wrong=kw,
                    med_err=med,
                    full=fs_,
                    fallbacks=int(sum(fb[k] for k in keys)),
                )
                res["bands"][band].setdefault(var, {})[p] = d
                lines.append(
                    f"  {p:5s} {len(keys):4d} {tfmt(t):>20s} {f3(tm['median']):>11s} "
                    f"{f3(npd) + '/' + f3(nmd):>16s} {f3(ev):>7s} "
                    f"{kw:>4d}/{len(keys):<4d} {fr(kw / len(keys)):>5s} {f3(med):>8s} "
                    f"{fr(fs_['wrong'])} [{fr(fs_['wrong_lo'])},{fr(fs_['wrong_hi'])}]".ljust(0)
                    + f" {f3(fs_['median']):>12s} {d['fallbacks']:>5d}"
                )
            lines.append("")
    # agreement of the subset with the full runs on the very same cases
    lines.append(
        "-- agreement with the full runs of Task 1 on the same cases (voxels 0-4, "
        "directions 0-3): median |err_timing - err_full| (deg) and fraction of cases whose "
        "wrong flag is the same --"
    )
    res["agreement"] = {}
    for p in ("H0", "H1", "H3", "HG", "MCr"):
        d_all, same_flag, n = [], 0, 0
        for k, e in errs.items():
            vidx, ri, j, var, pp = k
            if pp != p:
                continue
            ef = full.case_err(p, vidx, ri, j, VARIANTS.index(var))
            if not np.isfinite(ef):
                continue
            d_all.append(abs(e - ef))
            same_flag += int((e > 1.0) == (ef > 1.0))
            n += 1
        if n:
            res["agreement"][p] = dict(
                n=n,
                median_abs_diff=float(np.median(d_all)),
                frac_exact=float(np.mean(np.array(d_all) < 1e-9)),
                same_wrong_flag=same_flag / n,
            )
            lines.append(
                f"  {p:4s} n={n}: median |diff| {f3(np.median(d_all))}, exactly equal in "
                f"{fr(np.mean(np.array(d_all) < 1e-9))}, same wrong flag "
                f"{fr(same_flag / n)}"
            )
    lines.append("")
    # net time: measured vs production-equivalent
    lines.append(
        "-- network stage time per case (s), median over cases at all radii: measured "
        "(harness rendering included) vs production-equivalent (rendering excluded) --"
    )
    res["net_time"] = {}
    for p in ("H1", "H3", "HG"):
        a = [v for k, v in nprod.items() if k[4] == p]
        b = [v for k, v in nmeas.items() if k[4] == p]
        res["net_time"][p] = dict(production=float(np.median(a)), measured=float(np.median(b)))
        lines.append(
            f"  {p}: measured {f3(np.median(b))}, production-equivalent {f3(np.median(a))}"
        )
    # breakdown by stage for H3 and HG (means of exclusive seconds)
    lines.append("")
    lines.append(
        "-- seeded stage breakdown: mean exclusive seconds per case (both variants, all "
        "radii where the pipeline runs) --"
    )
    res["stages"] = {}
    stage_names = [
        "prepare",
        "decode",
        "forward",
        "render",
        "render_physics",
        "render_distractors",
        "render_realism",
        "gn",
        "gn_extract",
        "gn_solve",
        "net",
        "finisher",
        "find_optimal",
        "variance",
        "evaluate_local",
        "evaluate_mc",
        "mc_other",
    ]
    lines.append(f"  {'pipe':5s} " + " ".join(f"{n[:9]:>9s}" for n in stage_names))
    for p in pipes:
        rp = [r for r in recs if r["pipe"] == p and r["rep"] == 0]
        acc = {
            n: float(np.mean([r["stages"]["exclusive"].get(n, 0.0) for r in rp]))
            for n in stage_names
        }
        res["stages"][p] = acc
        lines.append(f"  {p:5s} " + " ".join(f"{acc[n]:9.4f}" for n in stage_names))
    lines.append(
        "  (prepare = prepare_nominal: ROI / observer / window spec; decode = "
        "decode_windows; forward = the network forward pass; render* = harness rendering; "
        "gn* = Gauss-Newton; finisher = FindOptimal or the MC, exclusive of its nested "
        "stages; evaluate_* = cost-function calls)"
    )
    return res, lines


_STORE: Dict[str, Any] = {}


def band_info(pipe: str, var: str, ris: Sequence[int]) -> Dict[str, Any]:
    """Timing statistics (this study) and full-run accuracy (Task 1) of a pipeline over radii."""
    vals = [v for k, v in _STORE["tprod"].items() if k[4] == pipe and k[3] == var and k[1] in ris]
    full = _STORE["full"].stats_band(pipe, ris, VARIANTS.index(var))
    return dict(time=tstat(vals), full=full)


def u2_eval(
    p: str, var: str, ris: Sequence[int], u0: Dict[str, Any], t2: Dict[str, Any]
) -> Dict[str, Any]:
    """The plan's U2 rule for one NN pipeline, variant and radius set. Times are MEANS per case
    for the pipeline, H0 and the fallbacks (H0 + fallback = mean H0 time + H0's wrong rate x the
    mean time of the fallback pipeline, with an ideal trigger). Fallbacks: the baseline
    reconstruction (symmetry-agnostic) and F1b (cubic-specific, Task 2 seed-0 accuracy)."""
    b, h = band_info(p, var, ris), band_info("H0", var, ris)
    fp, f0 = b["full"], h["full"]
    p0, n0 = f0["wrong"], f0["n"]
    t_p, t_h0 = b["time"]["mean"], h["time"]["mean"]
    t_full = u0[var]["baseline"]["time"]["mean"]
    t_f1b = u0[var]["F1b"]["time"]["mean"]
    full_row, f1b_row = _t2_row(t2, var, "baseline"), _t2_row(t2, var, "F1b")
    w_fb, w_fb2 = p0 * full_row["rate"], p0 * f1b_row["rate"]
    hi_fb = wilson(int(round(w_fb * n0)), n0)[2]
    hi_fb2 = wilson(int(round(w_fb2 * n0)), n0)[2]
    opts = {
        "H0": (f0["wrong"], f0["wrong_hi"], f0["median"]),
        "full reconstruction": (full_row["rate"], full_row["hi"], full_row["median_err_right"]),
        FALLBACK: (w_fb, hi_fb, f0["median"]),
        "H0 + F1b fallback (cubic)": (w_fb2, hi_fb2, f0["median"]),
    }
    cubic = "H0 + F1b fallback (cubic)"
    agn = {k: v for k, v in opts.items() if k != cubic}

    def judge(o: Dict[str, Any], sp: float) -> Tuple[str, bool, bool, bool]:
        best_ = min(o, key=lambda k: o[k][0])
        m_ = fp["wrong"] <= o[best_][1] and fp["median"] <= o[best_][2] + TOL_MED
        a_ = bool(m_ and sp >= 1.5)
        only_ = bool(fp["wrong"] <= 0.01 and all(v[0] > 0.01 for v in o.values()))
        return best_, bool(m_), a_, only_

    t_fb, t_fb2 = t_h0 + p0 * t_full, t_h0 + p0 * t_f1b
    sp1, sp2 = t_fb / t_p, t_fb2 / t_p
    # symmetry-agnostic comparators (H0, full reconstruction, H0 + baseline fallback) ...
    best, matched, a, only = judge(agn, sp1)
    # ... and with the cubic-specific H0 + F1b fallback added (faster must hold against both)
    best_c, matched_c, a_c, only_c = judge(opts, min(sp1, sp2))
    return dict(
        helps=bool(a or only), a=a, only=only, matched=bool(matched), best=best, t_p=t_p,
        helps_cubic=bool(a_c or only_c), a_cubic=a_c, only_cubic=only_c,
        matched_cubic=matched_c, best_cubic=best_c,
        t_h0=t_h0, t_fb=t_fb, t_fb2=t_fb2, speed_fb=sp1, speed_fb_f1b=sp2,
        speedup_full=t_full / t_p, wrong=fp["wrong"], wrong_lo=fp["wrong_lo"],
        wrong_hi=fp["wrong_hi"], median=fp["median"], h0_wrong=p0, h0_median=f0["median"],
        w_fb=w_fb, w_fb2=w_fb2, t_full=t_full, full_wrong=full_row["rate"],
    )  # fmt: skip


def u2_line(p: str, var: str, e: Dict[str, Any]) -> str:
    return (
        f"{VAR_NAME[var]}: {p} mean {f3(e['t_p'])} s, wrong {fr(e['wrong'])} "
        f"[{fr(e['wrong_lo'])}, {fr(e['wrong_hi'])}], med err {f3(e['median'])}; H0 mean "
        f"{f3(e['t_h0'])} s wrong {fr(e['h0_wrong'])}; H0+fallback (ideal trigger, baseline "
        f"reconstruction {f3(e['t_full'])} s) {f3(e['t_fb'])} s wrong {fr(e['w_fb'])}; "
        f"H0+F1b fallback {f3(e['t_fb2'])} s wrong {fr(e['w_fb2'])}; full reconstruction wrong "
        f"{fr(e['full_wrong'])}; best non-NN option {e['best']}; matched {e['matched']}; "
        f"{f3(e['speed_fb'])}x faster than H0+fallback, {f3(e['speed_fb_f1b'])}x than H0+F1b "
        f"fallback; symmetry-agnostic comparators: best {e['best']}, matched {e['matched']}, (a) "
        f"{e['a']}, only <=1% {e['only']} -> {e['helps']}; with the cubic F1b fallback: best "
        f"{e['best_cubic']}, matched {e['matched_cubic']}, (a) {e['a_cubic']}, only <=1% "
        f"{e['only_cubic']} -> {e['helps_cubic']}; speed-up vs full "
        f"reconstruction {f3(e['speedup_full'])}x"
    )


def seeded_verdicts(
    res: Dict[str, Any], u0: Dict[str, Any], t2: Dict[str, Any]
) -> Tuple[Dict[str, Any], List[str]]:
    out: Dict[str, Any] = {}
    lines: List[str] = ["-- U1 / U2 verdicts --"]
    # U1: helps only if it beats H0 on time AND accuracy. The plan does not define "beats on
    # accuracy"; the rule used here has a 0.002 deg tie (as the paired win rates of Task 1) or
    # non-overlapping wrong-rate intervals; the literal reading (any lower median) is also given.
    for p in ("H1", "H3"):
        ok_v, lit_v, parts = [], [], []
        for var in VARIANTS:
            b = res["bands"]["U1 (r 0.05, 0.1)"][var]
            tp, t0 = b[p]["time"]["median"], b["H0"]["time"]["median"]
            fp, f0 = b[p]["full"], b["H0"]["full"]
            acc = (fp["median"] < f0["median"] - 0.002) or (fp["wrong_hi"] < f0["wrong_lo"])
            lit = (fp["median"] < f0["median"]) or (fp["wrong_hi"] < f0["wrong_lo"])
            ok = bool(tp < t0 and acc)
            ok_v.append(ok)
            lit_v.append(bool(tp < t0 and lit))
            parts.append(
                f"{VAR_NAME[var]}: time {f3(tp)} vs {f3(t0)} s ({100 * (tp / t0 - 1):+.1f}%), "
                f"full-run median err {f3(fp['median'])} vs {f3(f0['median'])}, wrong "
                f"{fr(fp['wrong'])} vs {fr(f0['wrong'])}; beats on accuracy (0.002 deg tie) "
                f"{acc}, literal {lit}, on time {tp < t0} -> {ok}"
            )
        out[f"U1 {p} vs H0"] = dict(helps=all(ok_v), helps_literal=all(lit_v))
        lines.append(
            f"  U1 {p} vs H0: helps = {all(ok_v)} (literal rule without the tie: {all(lit_v)})"
        )
        lines += [f"      {s}" for s in parts]
    for var in VARIANTS:
        b = res["bands"]["U1 (r 0.05, 0.1)"][var]
        lines.append(
            f"  (reference, U1 {VAR_NAME[var]}: MC told r {f3(b['MCr']['time']['median'])} s, "
            f"full-run median err {f3(b['MCr']['full']['median'])}; "
            f"H0 {f3(b['H0']['time']['median'])} s)"
        )
    # U2: pooled bands (means) and every radius
    out["U2"] = {}
    sets = [
        ("U2 (r 0.5 .. 3)", ("H1", "H3"), [3, 4, 5, 6, 7, 8]),
        ("U2 HG band (r 1.5 .. 3)", ("HG",), [6, 7, 8]),
        ("U2 stress (r 5)", ("H1", "H3", "HG"), [9]),
    ]
    for label, pp, ris in sets:
        for p in pp:
            ev = {var: u2_eval(p, var, ris, u0, t2) for var in VARIANTS}
            out["U2"].setdefault(label, {})[p] = dict(
                ev,
                helps_overall=all(e["helps"] for e in ev.values()),
                helps_overall_cubic=all(e["helps_cubic"] for e in ev.values()),
            )
            lines.append(
                f"  {label} {p}: helps = {out['U2'][label][p]['helps_overall']} "
                f"(with the cubic F1b fallback: {out['U2'][label][p]['helps_overall_cubic']})"
            )
            lines += [f"      {u2_line(p, v, e)}" for v, e in ev.items()]
    lines.append("  per radius (mean times; helps needs both variants):")
    for ri in range(3, 10):
        for p in ("H1", "H3", "HG"):
            if p == "HG" and ri not in HG_RI:
                continue
            ev = {var: u2_eval(p, var, [ri], u0, t2) for var in VARIANTS}
            out["U2"].setdefault(f"r={RADII[ri]:g}", {})[p] = dict(
                ev,
                helps_overall=all(e["helps"] for e in ev.values()),
                helps_overall_cubic=all(e["helps_cubic"] for e in ev.values()),
            )
            c, r_ = ev["clean"], ev["all"]
            lines.append(
                f"    r={RADII[ri]:g} {p}: helps "
                f"{out['U2'][f'r={RADII[ri]:g}'][p]['helps_overall']} "
                f"(cubic fallback {out['U2'][f'r={RADII[ri]:g}'][p]['helps_overall_cubic']})"
                f"; mean time clean/realistic {f3(c['t_p'])}/{f3(r_['t_p'])} s vs H0+fallback "
                f"{f3(c['t_fb'])}/{f3(r_['t_fb'])} s ({f3(c['speed_fb'])}x/{f3(r_['speed_fb'])}x), "
                f"vs H0+F1b fallback {f3(c['t_fb2'])}/{f3(r_['t_fb2'])} s "
                f"({f3(c['speed_fb_f1b'])}x/{f3(r_['speed_fb_f1b'])}x); wrong {fr(c['wrong'])}/"
                f"{fr(r_['wrong'])} vs H0 {fr(c['h0_wrong'])}/{fr(r_['h0_wrong'])}, "
                f"H0+fallback {fr(c['w_fb'])}/{fr(r_['w_fb'])}; speed-up vs full reconstruction "
                f"{f3(c['speedup_full'])}x/{f3(r_['speedup_full'])}x"
            )
    lines.append(
        "  per radius, mean production seconds per case [speed-up vs the full reconstruction]:"
    )
    for var in VARIANTS:
        t_full = u0[var]["baseline"]["time"]["mean"]
        lines.append(f"    {VAR_NAME[var]} (full reconstruction mean {f3(t_full)} s)")
        for ri in range(3, 10):
            cells = []
            for p in ("H0", "H1", "H3", "HG", "MCr"):
                t = band_info(p, var, [ri])["time"]
                if t["n"]:
                    cells.append(f"{p} {f3(t['mean'])} s [{f3(t_full / t['mean'])}x]")
            lines.append(f"      r={RADII[ri]:g}: " + ", ".join(cells))
    return out, lines


# ---------------------------------------------------------------------------
# Full-run accuracy of Task 1 (reused from scripts/nn_hybrid/summary.py)
# ---------------------------------------------------------------------------


class FullRuns:
    """Error arrays (voxel, radius, direction, variant) of the full runs of Task 1."""

    def __init__(self) -> None:
        import importlib.util

        p = PC.ICENINE_PY / "scripts" / "nn_hybrid" / "summary.py"
        spec = importlib.util.spec_from_file_location("nn_hybrid_summary", p)
        assert spec is not None and spec.loader is not None
        self.NS = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(self.NS)
        self.st = self.NS.Study(PC.ICENINE_PY / "benchmarks" / "nn_hybrid")
        for pipe in ("H0", "H1", "H3", "HG"):
            assert self.st.load(pipe), f"missing raw of {pipe}"
        self.err = {p_: self.st.pipe_err(p_) for p_ in ("H0", "H1", "H3", "HG")}
        self.err["MCr"] = self.st.baseline_err("MC")  # optimizer sweep, unreduced
        self.vpos = {int(v): i for i, v in enumerate(self.st.vox)}

    def stats_band(self, pipe: str, ris: Sequence[int], vi: int) -> Dict[str, float]:
        e = np.concatenate([self.NS.pick(self.err[pipe], ri, vi, self.st.mask) for ri in ris])
        return self.NS.stats(e, None)

    def case_err(self, pipe: str, vidx: int, ri: int, j: int, vi: int) -> float:
        return float(self.err[pipe][self.vpos[vidx], ri, j, vi])


# ---------------------------------------------------------------------------
# Contention, isolation, profilers
# ---------------------------------------------------------------------------


def contention(
    u1: List[Dict[str, Any]],
    u10: List[Dict[str, Any]],
    s1: List[Dict[str, Any]],
    s10: List[Dict[str, Any]],
    w10: Dict[str, Dict[str, float]],
) -> Tuple[Dict[str, Any], List[str]]:
    lines = ["-- contention: 10 workers vs 1 worker, same cases (median of per-run time ratios) --"]
    out: Dict[str, Any] = {}
    tv = lambda r: r.get("wall_total", r["wall"])  # noqa: E731
    a = case_median(u1, ("vidx", "variant", "pipe"), tv)
    b = case_median([r for r in u10 if r["rep"] == 0], ("vidx", "variant", "pipe"), tv)
    ratios = {}
    for p in U0_PIPES:
        r = [b[k] / a[k] for k in b if k in a and k[2] == p]
        if r:
            ratios[p] = dict(
                n=len(r), median=float(np.median(r)), p10=PC.pct(r, 10), p90=PC.pct(r, 90)
            )
    out["u0"] = ratios
    allr = [b[k] / a[k] for k in b if k in a]
    out["u0_all"] = dict(n=len(allr), median=float(np.median(allr)) if allr else float("nan"))
    lines.append(
        f"  U0 (matched voxel x variant x pipeline, n={len(allr)}): median ratio "
        f"{f3(out['u0_all']['median'])}"
    )
    for p, d in ratios.items():
        lines.append(
            f"    {p:9s} n={d['n']:3d} median {f3(d['median'])} [{f3(d['p10'])}-{f3(d['p90'])}]"
        )
    key = ("vidx", "ri", "j", "variant", "pipe")
    a2 = case_median(s1, key, lambda r: r["total_production"])
    b2 = case_median([r for r in s10 if r["rep"] == 0], key, lambda r: r["total_production"])
    ratios2 = {}
    for p in ("H0", "H1", "H3", "MCr", "HG"):
        r = [b2[k] / a2[k] for k in b2 if k in a2 and k[4] == p]
        if r:
            ratios2[p] = dict(
                n=len(r), median=float(np.median(r)), p10=PC.pct(r, 10), p90=PC.pct(r, 90)
            )
    out["seeded"] = ratios2
    allr2 = [b2[k] / a2[k] for k in b2 if k in a2]
    out["seeded_all"] = dict(
        n=len(allr2), median=float(np.median(allr2)) if allr2 else float("nan")
    )
    lines.append(
        f"  U1/U2 (matched cases, n={len(allr2)}): median ratio "
        f"{f3(out['seeded_all']['median'])}"
    )
    for p, d in ratios2.items():
        lines.append(
            f"    {p:9s} n={d['n']:3d} median {f3(d['median'])} [{f3(d['p10'])}-{f3(d['p90'])}]"
        )
    for name, w in w10.items():
        lines.append(f"  {name}: 10-worker wall {f3(w['wall_s'])} s for {int(w['tasks'])} tasks")
    return out, lines


def isolation_lines(tag_dirs: Dict[str, Path]) -> Tuple[Dict[str, Any], List[str]]:
    out: Dict[str, Any] = {}
    lines = [
        "-- isolation conditions --",
        "  note: start / end / wall below are the original runs; flagged tasks were re-timed"
        " later on 2026-10-06 (MIGRATION_HISTORY, Re-timing) and the task-start load"
        " statistics include those re-runs",
    ]
    for name, d in tag_dirs.items():
        if not d.exists():
            continue
        s = (
            json.loads((d / "isolation_start.json").read_text())
            if (d / "isolation_start.json").exists()
            else {}
        )
        e = (
            json.loads((d / "isolation_end.json").read_text())
            if (d / "isolation_end.json").exists()
            else {}
        )
        wall = json.loads((d / "wall.json").read_text()) if (d / "wall.json").exists() else {}
        la, busy = [], set()
        for f in sorted(d.glob("v*.json")):
            iso = json.loads(f.read_text())["isolation"]
            la.append(iso["loadavg"][0])
            busy.update(iso["processes_over_5pct_cpu"])
        out[name] = dict(start=s, end=e, wall=wall, task_loadavg1=la, busy=sorted(busy))
        lines.append(
            f"  {name}: start {s.get('time')} load {[round(x, 2) for x in s.get('loadavg', [])]}, "
            f"{s.get('power')}, battery {s.get('battery', '')[:40]}; end {e.get('time')} load "
            f"{[round(x, 2) for x in e.get('loadavg', [])]}; OMP/MKL threads "
            f"{s.get('OMP_NUM_THREADS')}/{s.get('MKL_NUM_THREADS')}, torch threads "
            f"{s.get('torch_threads')}; wall {f3(wall.get('wall_s', float('nan')))} s for "
            f"{wall.get('tasks')} tasks; task-start 1-min load: "
            f"median {f3(np.median(la)) if la else 'nan'}, "
            f"max {f3(max(la)) if la else 'nan'}; other processes above 5% CPU at task starts: "
            f"{sorted(busy) or 'none'}; thermal: {s.get('thermal', '')[:90]}"
        )
    return out, lines


def profiler_lines() -> Tuple[Dict[str, Any], List[str]]:
    out: Dict[str, Any] = {}
    lines = ["-- profilers (cProfile top 30 per pipeline: benchmarks/profiling/cprofile_*.txt) --"]
    for name in ("u0", "seeded"):
        p = OUT / f"cprofile_{name}_shares.json"
        if not p.exists():
            continue
        d = json.loads(p.read_text())
        out[f"cprofile_{name}"] = d
        lines.append(f"  cProfile {name}: share of the total self time (tottime) by code class")
        for pipe, sh in d.items():
            tot = sh.get("_total_tottime_s", float("nan"))
            parts = ", ".join(
                f"{k.split(' (')[0]} {fr(v)}"
                for k, v in sh.items()
                if not k.startswith("_") and k != "wall_under_profiler_s"
            )
            lines.append(f"    {pipe:9s} total {f3(tot)} s: {parts}")
    p = OUT / "mps_vs_cpu.json"
    if p.exists():
        d = json.loads(p.read_text())
        out["mps"] = d
        lines.append(
            "  forward pass, median ms per call [10-90%] (net realistic_s0, torch "
            f"{d['torch']}, MPS available {d['mps_available']}):"
        )
        for b in ("batch1", "batch20"):
            if b not in d:
                continue
            row = d[b]
            cells = [
                f"{k} {f3(v['median'])} [{f3(v['p10'])}-{f3(v['p90'])}]"
                for k, v in row.items()
                if isinstance(v, dict)
            ]
            lines.append(f"    batch {row['batch']}: " + "; ".join(cells))
            if "max_abs_diff_mean_deg" in row:
                lines.append(
                    f"      max |MPS - CPU| mean {row['max_abs_diff_mean_deg']:.2e} deg, "
                    f"chol {row['max_abs_diff_chol']:.2e}"
                )
            if "mps_resident" in row:
                c1 = row["cpu_1thread"]["median"]
                lines.append(
                    f"      MPS speed-up over CPU 1 thread: resident "
                    f"{f3(c1 / row['mps_resident']['median'])}x, with transfer "
                    f"{f3(c1 / row['mps_with_transfer']['median'])}x"
                )
    p = OUT / "overhead.json"
    if p.exists():
        d = json.loads(p.read_text())
        out["overhead"] = d
        lines.append(
            f"  stage-wrapper overhead: "
            f"{d['per_call']['overhead_s'] * 1e6:.2f} us per evaluate call "
            f"({d['paired']['evaluate_calls']} calls in a baseline run = "
            f"{100 * d['share_of_run_from_per_call']:.2f}% of it); paired plain vs wrapped run "
            f"of voxel {d['paired']['voxel']}: {d['paired']['median_plain_s']:.2f} vs "
            f"{d['paired']['median_wrapped_s']:.2f} s (x{d['paired']['wrapped_over_plain']:.4f})"
        )
    p = OUT / "torch_profiler_net.txt"
    if p.exists():
        txt = p.read_text().splitlines()
        lines.append("  torch.profiler (CPU), net x3 stage, record_function ranges:")
        lines += ["    " + t.strip() for t in txt[2:12] if t.strip().startswith("net::")]
    return out, lines


def write_gz(path: Path, obj: Any) -> None:
    """Deterministic gzip (no mtime, no file name in the header): regenerating is byte-stable."""
    with open(path, "wb") as raw:
        with gzip.GzipFile(filename="", fileobj=raw, mode="wb", mtime=0) as gz:
            with io.TextIOWrapper(gz, encoding="utf-8") as f:
                json.dump(obj, f, default=float)


def compact(recs: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
    keep = sorted({k for r in recs for k in r} - {"R_final", "ref_R"})  # sorted: stable key order
    out = []
    for r in recs:
        d = {k: r[k] for k in keep if k in r and k != "stages"}
        d["stages_exclusive"] = {k: round(v, 5) for k, v in r["stages"]["exclusive"].items()}
        out.append(d)
    return out


def main() -> None:
    t2 = json.loads((PC.ICENINE_PY / "benchmarks" / "coarse_proxy" / "summary.json").read_text())
    u1 = load_dir(CACHE / "u0" / "w1")
    s1 = load_dir(CACHE / "seeded" / "w1")
    u10 = load_dir(CACHE / "u0" / "w10") if (CACHE / "u0" / "w10").exists() else []
    s10 = load_dir(CACHE / "seeded" / "w10") if (CACHE / "seeded" / "w10").exists() else []
    full = FullRuns()
    summary: Dict[str, Any] = {}
    L: List[str] = [
        "Run-time profiling (Task 3). Caveats: per-voxel images (<= 3 distractor sources), not "
        "full-sample renders; the rendering of the detector windows is harness work (excluded "
        "from 'production' times); Apple M2 Max (8 P + 4 E cores), thermals not controlled.",
        "",
    ]
    iso, il = isolation_lines(
        {
            "U0 1 worker": CACHE / "u0" / "w1",
            "U1/U2 1 worker": CACHE / "seeded" / "w1",
            "U0 10 workers": CACHE / "u0" / "w10",
            "U1/U2 10 workers": CACHE / "seeded" / "w10",
        }
    )
    summary["isolation"] = iso
    L += il + [""]
    if u1:
        u0res, ul = u0_tables(u1, t2)
        summary["u0"] = {
            v: {
                p: {k: x for k, x in d.items() if k != "per_voxel"}
                for p, d in vv.items()
                if isinstance(d, dict) and "time" in d
            }
            for v, vv in u0res.items()
            if v in VARIANTS
        }
        summary["u0_stages"] = u0res["stages"]
        L += ul + [""]
        uv, uvl = u0_verdicts(u0res, t2)
        summary["u0_verdicts"] = uv
        L += uvl + [""]
        write_gz(OUT / "u0_runs.json.gz", compact(u1 + u10))
    if s1 and u1:
        sres, sl = seeded_tables(s1, full, u0res)
        summary["seeded"] = {k: v for k, v in sres.items()}
        L += sl + [""]
        sv, svl = seeded_verdicts(sres, u0res, t2)
        summary["seeded_verdicts"] = sv
        L += svl + [""]
        write_gz(OUT / "seeded_runs.json.gz", compact(s1 + s10))
    if u10 and s10 and u1 and s1:
        w10 = {}
        for name, d in (("U0", CACHE / "u0" / "w10"), ("U1/U2", CACHE / "seeded" / "w10")):
            if (d / "wall.json").exists():
                w10[name] = json.loads((d / "wall.json").read_text())
        cres, cl = contention(u1, u10, s1, s10, w10)
        summary["contention"] = cres
        L += cl + [""]
    pr, pl = profiler_lines()
    summary["profilers"] = pr
    L += pl
    PC.write_json(OUT / "summary.json", summary)
    (OUT / "summary.txt").write_text("\n".join(L) + "\n")
    print("\n".join(L))


if __name__ == "__main__":
    main()
