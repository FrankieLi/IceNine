#!/usr/bin/env python3
"""
Phase B3 summary: head-to-head of local finishers (bench.py caches) -> summary.json, tables.md.

Usage (from icenine_py/):
  uv run python scripts/finisher_bench/finisher_summary.py [--cache run] [--out DIR]

Error = angle between the result and the truth (symmetry-agnostic: local searches; the sweep and
T5 starts are within 3 deg). Wrong = error above 1 deg. Cost gap = cost(result) - cost(truth).
Pairing is against (i) MC as deployed (its own natural budget: a median 201 evaluations), case by
case, with a win rate (tie = |difference| < 0.002 deg, counted one half), an exact sign test over
the non-tied cases, McNemar on wrong and on < 0.02 deg, and a voxel-clustered sign test (the
per-voxel median difference, voxels with a non-tied median).
"""

import argparse
import glob
import json
import math
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

import numpy as np
from scipy.spatial.transform import Rotation
from scipy.stats import binomtest

HERE = Path(__file__).resolve().parent
ICENINE_PY = HERE.parents[1]
sys.path.insert(0, str(HERE.parent / "common"))

from doc_tables import markdown_table, write_tables  # noqa: E402
from stats import mcnemar_exact, paired_discordant, wilson, win_rate  # noqa: E402

CACHE_DIR = HERE / "cache"
OUT_DIR = ICENINE_PY / "benchmarks" / "finisher_bench"
VARIANTS = ["clean", "realistic"]
TIE_DEG = 0.002
WRONG_DEG = 1.0
WRONG_CUT = WRONG_DEG + 1e-6  # starts at exactly 1 deg (r = 1) must not flip on rounding
BUDGETS = [250, 1000, 2600, 10000]
LABEL = {
    "mc_deployed": "(i) MC deployed (200 steps)",
    "mc_april": "(i) MC April sweep (3500 steps)",
    "mc_local": "(ii) MC local restarts",
    "es_box": "(iii) (1+1)-ES, step0 0.1317 deg",
    "es_002": "(iii) (1+1)-ES, step0 0.02 deg",
    "nm": "(iv) Nelder-Mead",
    "cma_005": "(v) CMA-ES, sigma0 0.05 deg",
    "cma_02": "(v) CMA-ES, sigma0 0.2 deg",
    "vm_small": "(vi) VarianceMinimizing, quarter box",
    "gn": "(vii) centroid Huber GN",
    "adam": "(viii) hybrid Riemannian Adam",
    "finisher": "default finisher (MC + VM)",
    "start": "start (no move)",
}
ANYTIME = ["mc_april", "mc_local", "es_box", "es_002", "nm", "cma_005", "cma_02", "vm_small"]
SETS = {
    "T5": ("T5 (300 cases)", ("H3", "H0")),
    "T5-H3": ("T5 H3 starts (200)", ("H3",)),
    "T5-H0": ("T5 H0 starts (100)", ("H0",)),
    "SW": ("sweep r 0.05-1 deg", ("SW",)),
}
RADII = {0: 0.05, 1: 0.1, 2: 0.25, 3: 0.5, 5: 1.0}


# ---------------------------------------------------------------------------
# Loading
# ---------------------------------------------------------------------------


def angles(R: np.ndarray, Rt: np.ndarray) -> np.ndarray:
    """Angle (deg) of R Rt^T, broadcast over leading axes."""
    R, Rt = np.broadcast_arrays(R, Rt)
    M = (R @ np.swapaxes(Rt, -1, -2)).reshape(-1, 3, 3)
    return np.degrees(Rotation.from_matrix(M).magnitude()).reshape(R.shape[:-2])


def load(cache: Path) -> Dict[str, Any]:
    files = sorted(glob.glob(str(cache / "*.npz")))
    assert files, f"no caches in {cache}"
    parts = [dict(np.load(f)) for f in files]
    d: Dict[str, Any] = {}
    meta_kind, meta_vox, meta_ri, meta_j = [], [], [], []
    for p in parts:
        n = len(p["dirs"])
        meta_kind += [str(p["kind"])] * n
        meta_vox += [int(p["vidx"])] * n
        meta_ri += [int(p["ri"])] * n
        meta_j += [int(j) for j in p["dirs"]]
    d["kind"], d["vox"] = np.array(meta_kind), np.array(meta_vox)
    d["ri"], d["j"] = np.array(meta_ri), np.array(meta_j)
    d["methods"] = [str(m) for m in parts[0]["methods"]]
    d["ckpts"] = [int(c) for c in parts[0]["ckpts"]]
    for k in parts[0]:
        if k in ("dirs", "vidx", "ri", "kind", "methods", "ckpts"):
            continue
        d[k] = np.concatenate([p[k] for p in parts], axis=1)  # (variant, case, ...)
    d["gn_secs"] = np.concatenate([p["gn_secs"] for p in parts], axis=1)
    d["n_files"] = len(files)
    return d


# ---------------------------------------------------------------------------
# Rows: one (method, budget) with per-case error / cost / evaluations
# ---------------------------------------------------------------------------


def build_rows(d: Dict[str, Any]) -> List[Dict[str, Any]]:
    Rt = d["R_true"]
    ct = d["cost_true"]
    rows: List[Dict[str, Any]] = []

    def add(key: str, budget: str, err: np.ndarray, cost: np.ndarray, used: np.ndarray) -> None:
        rows.append(dict(key=key, budget=budget, err=err, gap=cost - ct, used=used.astype(float)))

    add("start", "0", angles(d["R_start"], Rt), d["cost_start"], np.zeros_like(ct))
    nck = len(d["ckpts"])
    for m in ["mc_deployed"] + ANYTIME:
        mi = d["methods"].index(m)
        if m == "mc_deployed":  # natural stop (about 201 evaluations)
            k = nck - 1
            err = angles(d["res_R"][:, :, mi, k], Rt)
            add(m, "natural", err, d["res_cost"][:, :, mi, k], d["res_used"][:, :, mi, k])
            continue
        for b in BUDGETS:
            k = d["ckpts"].index(b)
            err = angles(d["res_R"][:, :, mi, k], Rt)
            add(m, str(b), err, d["res_cost"][:, :, mi, k], d["res_used"][:, :, mi, k])
    add("finisher", "natural", angles(d["fin_R"], Rt), d["fin_cost"], d["fin_evals"])
    add("gn", "natural", angles(d["gn_R"], Rt), d["gn_cost"], np.zeros_like(ct))
    add("adam", "natural", angles(d["adam_R"], Rt), d["adam_cost"], d["adam_hard"])
    return rows


# ---------------------------------------------------------------------------
# Statistics
# ---------------------------------------------------------------------------


def frac(k: int, n: int) -> Dict[str, Any]:
    lo, hi = wilson(int(k), int(n))
    return dict(k=int(k), n=int(n), p=k / max(n, 1), lo=lo, hi=hi)


def fs(f: Dict[str, Any]) -> str:
    return f"{100 * f['p']:.1f}% ({100 * f['lo']:.1f}-{100 * f['hi']:.1f})"


def sign_test(diff: np.ndarray, tie: float = TIE_DEG) -> Dict[str, Any]:
    """diff = err_method - err_reference; negative = the method is better. Exact two-sided."""
    better, worse = int((diff < -tie).sum()), int((diff > tie).sum())
    n = better + worse
    p = float(binomtest(better, n, 0.5).pvalue) if n else 1.0
    return dict(better=better, worse=worse, ties=int(len(diff) - n), p=p)


def voxel_sign_test(diff: np.ndarray, vox: np.ndarray, tie: float = TIE_DEG) -> Dict[str, Any]:
    meds = np.array([np.median(diff[vox == v]) for v in np.unique(vox)])
    out = sign_test(meds, tie)
    out["n_vox"] = int(len(meds))
    return out


def row_stats(
    row: Dict[str, Any], dep: Dict[str, Any], vi: int, mask: np.ndarray, vox: np.ndarray
) -> Dict[str, Any]:
    e, g, u = row["err"][vi][mask], row["gap"][vi][mask], row["used"][vi][mask]
    n = int(mask.sum())
    ed = dep["err"][vi][mask]
    out: Dict[str, Any] = dict(
        n=n,
        evals_median=float(np.median(u)),
        err_median=float(np.median(e)),
        err_q25=float(np.quantile(e, 0.25)),
        err_q75=float(np.quantile(e, 0.75)),
        err_p90=float(np.quantile(e, 0.9)),
        lt001=frac((e < 0.01).sum(), n),
        lt002=frac((e < 0.02).sum(), n),
        wrong=frac((e > WRONG_CUT).sum(), n),
        gap_median=float(np.median(g)),
        gap_q75=float(np.quantile(g, 0.75)),
        below_truth_cost=frac((g < 0).sum(), n),
    )
    wr, tie_frac, _ = win_rate(e, ed, TIE_DEG)
    out["vs_deployed"] = dict(
        win_rate=wr,
        tie_frac=tie_frac,
        sign=sign_test(e - ed),
        voxel_sign=voxel_sign_test(e - ed, vox[mask]),
    )
    b, c = paired_discordant(e > WRONG_CUT, ed > WRONG_CUT)
    out["vs_deployed"]["mcnemar_wrong"] = dict(
        only_method_wrong=b, only_deployed_wrong=c, p=mcnemar_exact(b, c)
    )
    b, c = paired_discordant(e >= 0.02, ed >= 0.02)
    out["vs_deployed"]["mcnemar_lt002"] = dict(
        only_method_fails=b, only_deployed_fails=c, p=mcnemar_exact(b, c)
    )
    return out


def fmt_p(p: float) -> str:
    return f"{p:.2g}" if p >= 1e-3 else f"{p:.1e}"


# ---------------------------------------------------------------------------
# Tables
# ---------------------------------------------------------------------------


def rowname(key: str, budget: str) -> str:
    return LABEL[key] + ("" if budget in ("natural", "0") else f" @ {budget}")


def build_summary(d: Dict[str, Any]) -> Dict[str, Any]:
    rows = build_rows(d)
    dep = next(r for r in rows if r["key"] == "mc_deployed")
    S: Dict[str, Any] = dict(
        n_files=d["n_files"],
        n_cases=int(len(d["kind"])),
        budgets=BUDGETS,
        tie_deg=TIE_DEG,
        wrong_deg=WRONG_DEG,
        sets={},
    )
    for sname, (_lab, kinds) in SETS.items():
        mask = np.isin(d["kind"], kinds)
        if not mask.any():  # partial caches only
            continue
        S["sets"][sname] = dict(
            n_cases=int(mask.sum()), n_voxels=int(len(np.unique(d["vox"][mask])))
        )
        for vi, var in enumerate(VARIANTS):
            S["sets"][sname][var] = {
                f"{r['key']}@{r['budget']}": row_stats(r, dep, vi, mask, d["vox"]) for r in rows
            }
    # sweep by radius (realistic and clean), at natural / 2600
    S["sw_by_r"] = {}
    for vi, var in enumerate(VARIANTS):
        S["sw_by_r"][var] = {}
        for ri, r in RADII.items():
            mask = (d["kind"] == "SW") & (d["ri"] == ri)
            if mask.sum() == 0:
                continue
            S["sw_by_r"][var][str(r)] = {
                f"{x['key']}@{x['budget']}": row_stats(x, dep, vi, mask, d["vox"])
                for x in rows
                if x["budget"] in ("natural", "2600", "0")
            }
            S["sw_by_r"][var][str(r)]["n"] = int(mask.sum())
    S["strata"] = build_strata(d, rows)
    S["deployed_mc"] = build_deployed(d)
    S["extras"] = build_extras(d)
    S["units"] = build_units(d)
    S["checks"] = build_checks(d)
    return S


def build_strata(d: Dict[str, Any], rows: List[Dict[str, Any]]) -> Dict[str, Any]:
    """Cases split by what the deployed MC did: no improvement over the start (never accepted a
    lower cost) or an improvement. Median error of a few methods at 1000 evaluations."""
    mi = d["methods"].index("mc_deployed")
    k = len(d["ckpts"]) - 1
    improved = d["res_cost"][:, :, mi, k] < d["cost_start"] - 1e-12  # (V, N)
    out: Dict[str, Any] = {}
    pick = [("mc_deployed", "natural"), ("start", "0")] + [
        (m, "1000") for m in ("es_box", "es_002", "nm", "cma_005", "cma_02", "vm_small")
    ]
    for sname, kinds in (("T5-H3", ("H3",)), ("SW", ("SW",))):
        base = np.isin(d["kind"], kinds)
        for vi, var in enumerate(VARIANTS):
            for lab, sel in (("no_improvement", ~improved[vi]), ("improved", improved[vi])):
                mask = base & sel
                cell: Dict[str, Any] = dict(n=int(mask.sum()))
                for key, b in pick:
                    r = next(x for x in rows if x["key"] == key and x["budget"] == b)
                    if mask.sum():
                        e = r["err"][vi][mask]
                        cell[f"{key}@{b}"] = dict(
                            err_median=float(np.median(e)),
                            lt002=frac((e < 0.02).sum(), int(mask.sum())),
                        )
                out[f"{sname}|{var}|{lab}"] = cell
    return out


def build_deployed(d: Dict[str, Any]) -> Dict[str, Any]:
    """How the deployed MC ended: no improvement over the start, restarts exhausted (ended before
    its 200 steps; the max-convergence stop needs a cost below 1e-4, which does not occur here),
    and the overlap of the two."""
    mi = d["methods"].index("mc_deployed")
    k = len(d["ckpts"]) - 1
    improved = d["res_cost"][:, :, mi, k] < d["cost_start"] - 1e-12
    early = d["total_evals"][:, :, mi] < 200
    out: Dict[str, Any] = {}
    for sname, (_lab, kinds) in SETS.items():
        base = np.isin(d["kind"], kinds)
        for vi, var in enumerate(VARIANTS):
            n = int(base.sum())
            ni, ee = ~improved[vi] & base, early[vi] & base
            out[f"{sname}|{var}"] = dict(
                n=n,
                no_improvement=frac(int(ni.sum()), n),
                ended_early=frac(int(ee.sum()), n),
                both=int((ni & ee).sum()),
                early_not_noimp=int((ee & ~ni).sum()),
                noimp_not_early=int((ni & ~ee).sum()),
            )
    return out


def build_extras(d: Dict[str, Any]) -> Dict[str, Any]:
    """Overlap of the T5 H0 and sweep case sets, and how often a method's result at 250
    evaluations is the same orientation as at 10000 (it found nothing better afterwards)."""
    key = list(zip(d["kind"], d["vox"].tolist(), d["ri"].tolist(), d["j"].tolist()))
    h0 = {k[1:] for k in key if k[0] == "H0"}
    sw = {k[1:] for k in key if k[0] == "SW"}
    out: Dict[str, Any] = dict(t5_h0_cases=len(h0), sw_cases=len(sw), overlap=len(h0 & sw))
    k250, kmax = d["ckpts"].index(250), d["ckpts"].index(10000)
    same: Dict[str, Any] = {}
    for m in ["mc_local", "es_box", "es_002", "nm", "cma_005", "cma_02", "vm_small"]:
        mi = d["methods"].index(m)
        eq = np.all(d["res_R"][:, :, mi, k250] == d["res_R"][:, :, mi, kmax], axis=(-1, -2))
        for vi, var in enumerate(VARIANTS):
            for sname in ("T5-H3", "SW"):
                mask = np.isin(d["kind"], SETS[sname][1])
                same[f"{m}|{sname}|{var}"] = frac(int(eq[vi][mask].sum()), int(mask.sum()))
    out["same_at_250_and_10000"] = same
    return out


def build_units(d: Dict[str, Any]) -> Dict[str, Any]:
    def q(x: np.ndarray) -> List[float]:
        return [float(v) for v in np.quantile(x, [0.25, 0.5, 0.75])]

    out: Dict[str, Any] = {}
    for vi, var in enumerate(VARIANTS):
        out[var] = dict(
            gn=dict(
                passes=3,
                ok=frac(int((d["gn_ok"][vi] > 0).sum()), int(d["gn_ok"][vi].size)),
                n_spots_used_q=q(d["gn_n_used"][vi]),
                wall_s_per_case_q=q(d["gn_secs"][vi]),
            ),
            adam=dict(
                hard_evals_q=q(d["adam_hard"][vi]),
                diff_steps_q=q(d["adam_diff"][vi]),
                wall_s_per_case_q=q(d["adam_secs"][vi]),
            ),
            finisher=dict(evals_q=q(d["fin_evals"][vi]), wall_s_per_case_q=q(d["fin_secs"][vi])),
        )
        for m in ["mc_deployed"] + ANYTIME:
            mi = d["methods"].index(m)
            out[var][m] = dict(
                total_evals_q=q(d["total_evals"][vi][:, mi].astype(float)),
                wall_s_to_10000_q=q(d["secs"][vi][:, mi]),
            )
    return out


def build_checks(d: Dict[str, Any]) -> Dict[str, Any]:
    raw = d["fin_vs_raw_maxabs"]
    fin = raw[np.isfinite(raw)]
    mi = d["methods"].index("mc_deployed")
    over = int((d["total_evals"] > 10000).sum())
    return dict(
        finisher_vs_stored_cases=int(fin.size),
        finisher_vs_stored_max_abs=float(fin.max()) if fin.size else math.nan,
        finisher_vs_stored_equal=int((fin == 0).sum()),
        max_total_evals=int(d["total_evals"].max()),
        runs_over_budget=over,
        deployed_total_evals_q=[
            float(v) for v in np.quantile(d["total_evals"][:, :, mi].astype(float), [0, 0.5, 1])
        ],
    )


# ---------------------------------------------------------------------------
# Table text
# ---------------------------------------------------------------------------


def tables(S: Dict[str, Any], timing: Optional[Dict[str, Any]]) -> Dict[str, str]:
    T: Dict[str, str] = {}
    keys_order = ["start@0", "mc_deployed@natural"]
    for m in ["mc_april", "mc_local", "es_box", "es_002", "nm", "cma_005", "cma_02", "vm_small"]:
        keys_order += [f"{m}@{b}" for b in BUDGETS]
    keys_order += ["finisher@natural", "gn@natural", "adam@natural"]

    def name(k: str) -> str:
        m, b = k.split("@")
        return rowname(m, b)

    long_cols = [
        "method", "evals", "median err", "p90 err", "<0.01 deg", "<0.02 deg", "wrong >1 deg",
        "cost gap median", "win vs (i)", "sign p", "voxel-sign p", "McNemar p (wrong)",
        "McNemar p (<0.02)",
    ]  # fmt: skip
    for sname in ("T5", "SW"):
        for var in VARIANTS:
            if sname not in S["sets"]:
                continue
            rows = []
            for k in keys_order:
                s = S["sets"][sname][var][k]
                v = s["vs_deployed"]
                rows.append(
                    {
                        "method": name(k),
                        "evals": s["evals_median"],
                        "median err": s["err_median"],
                        "p90 err": s["err_p90"],
                        "<0.01 deg": fs(s["lt001"]),
                        "<0.02 deg": fs(s["lt002"]),
                        "wrong >1 deg": f"{s['wrong']['k']}/{s['wrong']['n']} " + fs(s["wrong"]),
                        "cost gap median": s["gap_median"],
                        "win vs (i)": v["win_rate"],
                        "sign p": fmt_p(v["sign"]["p"]),
                        "voxel-sign p": fmt_p(v["voxel_sign"]["p"]),
                        "McNemar p (wrong)": fmt_p(v["mcnemar_wrong"]["p"]),
                        "McNemar p (<0.02)": fmt_p(v["mcnemar_lt002"]["p"]),
                    }
                )
            T[f"long_{sname}_{var}"] = markdown_table(
                rows,
                long_cols,
                formats={"evals": ".0f", "median err": ".4f", "p90 err": ".4f",
                         "cost gap median": ".4f", "win vs (i)": ".2f"},
            )  # fmt: skip
    # compact headline: method x budget -> median err / <0.02 / wrong
    for sname in SETS:
        for var in VARIANTS:
            if sname not in S["sets"]:
                continue
            rows = []
            lines = [("start (no move)", "start@0", None), ("(i) MC deployed (200 steps)",
                     "mc_deployed@natural", None)]  # fmt: skip
            for m in ANYTIME:
                lines.append((LABEL[m], m, BUDGETS))
            for m in ("finisher", "gn", "adam"):
                lines.append((LABEL[m], f"{m}@natural", None))
            for lab, key, buds in lines:
                row: Dict[str, Any] = {"method": lab}
                if buds is None:
                    s = S["sets"][sname][var][key]
                    row["natural"] = cell(s)
                    row["evals (median)"] = f"{s['evals_median']:.0f}"
                else:
                    for b in BUDGETS:
                        row[str(b)] = cell(S["sets"][sname][var][f"{key}@{b}"])
                rows.append(row)
            T[f"head_{sname}_{var}"] = markdown_table(
                rows, ["method", "natural", "evals (median)"] + [str(b) for b in BUDGETS]
            )
    # sweep by radius at 2600 and natural
    for var in VARIANTS:
        if not S["sw_by_r"][var]:
            continue
        rows = []
        for key in ("start@0", "mc_deployed@natural", "mc_april@2600", "mc_local@2600",
                    "es_box@2600", "es_002@2600", "nm@2600", "cma_005@2600", "cma_02@2600",
                    "vm_small@2600", "finisher@natural", "gn@natural", "adam@natural"):  # fmt: skip
            row = {"method": name(key)}
            for r in sorted(S["sw_by_r"][var], key=float):
                row[f"r={r}"] = cell(S["sw_by_r"][var][r][key])
            rows.append(row)
        T[f"sw_by_r_{var}"] = markdown_table(
            rows, ["method"] + [f"r={r}" for r in sorted(S["sw_by_r"][var], key=float)]
        )
    # strata
    rows = []
    for k, c in S["strata"].items():
        for mk, v in c.items():
            if mk == "n":
                continue
            rows.append(
                dict(
                    stratum=k.replace("|", " "),
                    n=c["n"],
                    method=name(mk),
                    med=v["err_median"],
                    lt=fs(v["lt002"]),
                )
            )
    T["strata"] = markdown_table(
        rows, ["stratum", "n", "method", "med", "lt"], formats={"med": ".4f"}
    )
    rows = []
    for k, v in S["deployed_mc"].items():
        rows.append(
            dict(
                set=k.replace("|", " "), n=v["n"], no_improvement=fs(v["no_improvement"]),
                ended_early=fs(v["ended_early"]), both=v["both"],
                early_only=v["early_not_noimp"], noimp_only=v["noimp_not_early"],
            )
        )  # fmt: skip
    T["deployed_mc"] = markdown_table(
        rows, ["set", "n", "no_improvement", "ended_early", "both", "early_only", "noimp_only"]
    )
    # units
    rows = []
    for var in VARIANTS:
        u = S["units"][var]
        g, a, f = u["gn"], u["adam"], u["finisher"]
        sq = g["n_spots_used_q"]
        rows += [
            dict(
                variant=var,
                method="(vii) centroid Huber GN",
                unit="3 window renders + solves",
                count=f"3 (spots used q25/50/75 {sq[0]:.0f}/{sq[1]:.0f}/{sq[2]:.0f}; "
                f"estimate for {fs(g['ok'])})",
                wall=f"{g['wall_s_per_case_q'][1]:.2f}",
            ),
            dict(
                variant=var,
                method="(viii) hybrid Riemannian Adam",
                unit="hard evals + Adam steps",
                count=f"{a['hard_evals_q'][1]:.0f} hard + {a['diff_steps_q'][1]:.0f} "
                "differentiable fwd+bwd (medians)",
                wall=f"{a['wall_s_per_case_q'][1]:.2f}",
            ),
            dict(
                variant=var,
                method="default finisher (MC + VM)",
                unit="hard evals",
                count=f"{f['evals_q'][1]:.0f} (q25-q75 {f['evals_q'][0]:.0f}"
                f"-{f['evals_q'][2]:.0f})",
                wall=f"{f['wall_s_per_case_q'][1]:.2f}",
            ),
        ]
    T["units"] = markdown_table(
        rows, ["variant", "method", "unit", "count", "wall"],
        labels=None,
    )  # fmt: skip
    if timing:
        T["timing"] = timing_table(timing)
    return T


def cell(s: Dict[str, Any]) -> str:
    return (
        f"{s['err_median']:.4f} / {100 * s['lt002']['p']:.0f}% / " f"{100 * s['wrong']['p']:.1f}%"
    )


def timing_table(t: Dict[str, Any]) -> str:
    rows = []
    for k, v in t["methods"].items():
        base, _, bud = k.partition("@")
        rows.append(
            dict(
                method=LABEL.get(base, base),
                budget=bud or v["budget"],
                evals=v["evals_median"],
                wall=v["wall_median_s"],
                per_eval=v["ms_per_eval"],
                q=f"{v['wall_q25_s']:.3f}-{v['wall_q75_s']:.3f}",
            )
        )
    return markdown_table(
        rows,
        ["method", "budget", "evals", "wall", "per_eval", "q"],
        formats={"evals": ".0f", "wall": ".3f", "per_eval": ".3f"},
    )


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--cache", default="run")
    ap.add_argument("--out", default=str(OUT_DIR))
    args = ap.parse_args()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    d = load(CACHE_DIR / args.cache)
    S = build_summary(d)
    tpath = out / "timing.json"
    timing = json.loads(tpath.read_text()) if tpath.exists() else None
    if timing:
        S["timing"] = timing
    meta = out / "run_meta.json"
    if meta.exists():
        S["run"] = json.loads(meta.read_text())
    (out / "summary.json").write_text(json.dumps(S, indent=1, default=float))
    write_tables(out / "tables.md", tables(S, timing))
    print(f"{S['n_cases']} cases x 2 variants from {S['n_files']} files; wrote {out}")


if __name__ == "__main__":
    main()
