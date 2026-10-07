#!/usr/bin/env python3
"""
Summarise the T5 finisher-diagnosis cache into benchmarks/finisher_diagnosis/
(summary.json, summary.txt, tables.md written with doc_tables.write_tables).

Usage (from icenine_py/):
  uv run python scripts/finisher_diagnosis/diag_summary.py [--cache DIR]
"""

import argparse
import glob
import json
import sys
from pathlib import Path
from typing import Any, Dict, List

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent / "common"))

import diagnose as D  # noqa: E402
from doc_tables import markdown_table, write_tables  # noqa: E402
from stats import wilson  # noqa: E402

TOL = 1e-9  # costs are float64 of float32-evaluated quantities; "equal" within TOL
RADII = [0.05, 0.1, 0.25, 0.5, 0.75, 1.0, 1.5, 2.0, 3.0]


def load(cache: Path, pipe: str) -> Dict[str, np.ndarray]:
    files = sorted(glob.glob(str(cache / f"{pipe}_v*_r*.npz")))
    parts = [np.load(f) for f in files]
    out: Dict[str, np.ndarray] = {}
    for k in parts[0].files:
        if k in ("vidx", "ri", "pipe"):
            continue
        out[k] = np.concatenate([p[k] for p in parts])
    out["ri_case"] = np.concatenate([np.full(len(p["dirs"]), int(p["ri"])) for p in parts])
    out["vidx_case"] = np.concatenate([np.full(len(p["dirs"]), int(p["vidx"])) for p in parts])
    return out


def frac(mask: np.ndarray) -> Dict[str, Any]:
    n, k = int(mask.size), int(mask.sum())
    lo, hi = wilson(k, n)
    return dict(k=k, n=n, p=k / max(n, 1), lo=lo, hi=hi)


def fmt(f: Dict[str, Any]) -> str:
    return f"{f['k']}/{f['n']} ({100 * f['p']:.0f}%, {100 * f['lo']:.0f}-{100 * f['hi']:.0f})"


def q(x: np.ndarray, qs=(0.25, 0.5, 0.75)) -> List[float]:
    return [float(v) for v in np.nanquantile(x, qs)]


def summarise(d: Dict[str, np.ndarray]) -> Dict[str, Any]:
    n = len(d["cost_res"])
    gap = d["cost_res"] - d["cost_true"]
    above = gap > TOL
    s: Dict[str, Any] = dict(n=n)
    s["rerun_maxabs_max"] = float(d["rerun_maxabs"].max())
    s["result_above_truth"] = frac(above)
    s["gap_quartiles"] = q(gap)
    s["ang_start_true_med"] = float(np.median(d["ang_start_true"]))
    s["ang_res_true_quartiles"] = q(d["ang_res_true"])
    # (a) geodesic
    ts, gc = d["geo_t"][0], d["geo_cost"]
    c0 = gc[:, 0]
    inter = (ts > 0) & (ts < 1)
    barrier = gc[:, inter].max(axis=1) - c0
    s["geo"] = dict(
        barrier_gt0=frac(barrier > TOL),
        barrier_quartiles=q(barrier),
        barrier_gt_gap=frac(barrier > np.maximum(gap, 0) + TOL),
        some_point_below_result=frac((gc[:, inter].min(axis=1) < c0 - TOL)),
        min_point_below_truth_plus_tol=frac(gc[:, inter].min(axis=1) <= d["cost_true"] + TOL),
        end_equals_truth=float(np.abs(gc[:, -1] - d["cost_true"]).max()),
        # first t where the path cost rises above the result's: how close to the result
        first_rise_t_median=float(
            np.median(
                [
                    ts[1:][g[1:] > g[0] + TOL][0] if (g[1:] > g[0] + TOL).any() else np.nan
                    for g in gc
                ]
            )
        ),
        monotone_nonincreasing=frac((np.diff(gc, axis=1) <= TOL).all(axis=1)),
    )
    # (c) granularity
    gt, gr = d["gran_truth"], d["gran_res"]  # (n, angles, dirs)
    ang = D.GRAN_ANGLES_DEG
    rows = []
    for ai, a in enumerate(ang):
        dt = gt[:, ai] - d["cost_true"][:, None]
        dr = gr[:, ai] - d["cost_res"][:, None]
        rows.append(
            dict(
                angle_deg=float(a),
                truth_same=float(np.mean(np.abs(dt) <= TOL)),
                truth_lower=float(np.mean(dt < -TOL)),
                truth_med_abs=float(np.median(np.abs(dt))),
                res_same=float(np.mean(np.abs(dr) <= TOL)),
                res_lower=float(np.mean(dr < -TOL)),
                res_med_abs=float(np.median(np.abs(dr))),
            )
        )
    s["gran"] = rows
    # cost quantum: smallest nonzero |difference| between distinct costs among a case's evals
    quanta = []
    for i in range(n):
        v = np.unique(np.round(np.concatenate([gt[i].ravel(), [d["cost_true"][i]]]), 12))
        if len(v) > 1:
            quanta.append(np.diff(v).min())
    s["quantum_median"] = float(np.median(quanta)) if quanta else float("nan")
    # share of the gap explained by granularity: best cost within 0.05 deg of the truth
    j05 = int(np.argmin(np.abs(ang - 0.05)))
    best05 = gt[:, : j05 + 1].min(axis=(1, 2)) - d["cost_true"]
    s["truth_best_within_0.05deg"] = dict(
        lower_than_truth=frac(best05 < -TOL),
        drop_quartiles=q(-best05),
        gap_over_drop_median=(
            float(np.median(gap[best05 < -TOL] / -best05[best05 < -TOL]))
            if (best05 < -TOL).any()
            else float("nan")
        ),
    )
    s["res_local_min_0.05deg"] = frac(gr[:, : j05 + 1].min(axis=(1, 2)) >= d["cost_res"] - TOL)
    s["counts_truth_median"] = [float(v) for v in np.median(d["counts_truth"], axis=0)]
    s["counts_res_median"] = [float(v) for v in np.median(d["counts_res"], axis=0)]
    # (b) continuations
    cont = {}
    for kind in D.CONT_NAMES:
        c, a, mv = d[f"c_{kind}_cost"], d[f"c_{kind}_ang"], d[f"c_{kind}_moved"]
        ref = d["cost_true"] if kind != "from_truth" else d["cost_true"]
        start_cost = d["cost_res"] if kind != "from_truth" else d["cost_true"]
        start_ang = d["ang_res_true"] if kind != "from_truth" else np.zeros(n)
        cont[kind] = dict(
            improves=frac(c < start_cost - TOL),
            reaches_truth_cost=frac(c <= ref + TOL),
            cost_drop_quartiles=q(start_cost - c),
            remaining_gap_quartiles=q(c - ref),
            gap_closed_median=(
                float(np.median(np.clip((start_cost - c)[gap > TOL] / gap[gap > TOL], 0, None)))
                if kind != "from_truth"
                else float("nan")
            ),
            toward_truth_quartiles=q(start_ang - a),
            closer_to_truth=frac(a < start_ang - 1e-6),
            ang_end_quartiles=q(a),
            moved_quartiles=q(mv),
        )
        vm = d[f"c_{kind}_vm"]
        if np.isfinite(vm[:, 0]).any():
            cont[kind]["vm_capped"] = frac(vm[:, D.VM_KEYS.index("capped")] > 0)
    s["cont"] = cont
    five = ["mc_long", "mc_smallstep", "mc_smallbox", "vm_long", "vm_smallbox"]
    best = np.minimum.reduce([d[f"c_{k}_cost"] for k in five])
    pos = gap > TOL
    s["best_of_five"] = dict(
        reaches_truth_cost=frac(best <= d["cost_true"] + TOL),
        gap_closed_median=float(np.median(((d["cost_res"] - best)[pos] / gap[pos]))),
    )
    # (d) stopping
    fo, vm = d["fo_log"], d["vm_log"]
    L, V = D.LOG_KEYS, D.VM_KEYS
    stop = fo[:, L.index("stop")].astype(int)
    s["stop"] = dict(
        fo_codes={nm: int((stop == i).sum()) for i, nm in enumerate(D.STOP_NAMES[:3])},
        fo_steps_run_quartiles=q(fo[:, L.index("steps_run")]),
        fo_n_accept_quartiles=q(fo[:, L.index("n_accept")]),
        fo_last_accept_quartiles=q(fo[:, L.index("last_accept")]),
        fo_restarts_max=int(fo[:, L.index("n_restarts")].max()),
        fo_final_step_deg_quartiles=q(fo[:, L.index("final_step_deg")]),
        fo_min_ergodic_quartiles=q(fo[:, L.index("min_ergodic")]),
        fo_since_improve_quartiles=q(fo[:, L.index("since_improve")]),
        fo_cost_drop_quartiles=q(fo[:, L.index("cost_start")] - fo[:, L.index("cost_end")]),
        fo_converged_flag=frac(d["fo_converged_flag"]),
        vm_steps_taken_quartiles=q(vm[:, V.index("steps_taken")]),
        vm_budget_extension_quartiles=q(
            vm[:, V.index("max_steps_end")] - float(D.fs._W.rec.params.max_mc_steps)
            if D.fs._W is not None
            else vm[:, V.index("max_steps_end")] - 200.0
        ),
        vm_n_improve_quartiles=q(vm[:, V.index("n_improve")]),
        vm_last_improve_frac_quartiles=q(vm[:, V.index("last_improve")] / vm[:, V.index("n_sub")]),
        vm_final_radius_deg_quartiles=q(vm[:, V.index("final_radius_deg")]),
        vm_final_variance_quartiles=q(vm[:, V.index("final_variance")]),
        vm_cost_drop_quartiles=q(vm[:, V.index("cost_start")] - vm[:, V.index("cost_end")]),
        vm_improves=frac(vm[:, V.index("cost_end")] < vm[:, V.index("cost_start")] - TOL),
        vm_gives_final=frac(vm[:, V.index("cost_end")] < fo[:, L.index("cost_end")] - TOL),
        vm_capped_in_default=frac(vm[:, V.index("capped")] > 0),
    )
    # by radius
    by_r = []
    for ri in range(D.N_RADII):
        m = d["ri_case"] == ri
        if not m.any():
            continue
        by_r.append(
            dict(
                radius_deg=RADII[ri],
                n=int(m.sum()),
                above=frac(above[m]),
                gap_med=float(np.median(gap[m])),
                ang_res_med=float(np.median(d["ang_res_true"][m])),
                mc_long_reach=frac(d["c_mc_long_cost"][m] <= d["cost_true"][m] + TOL),
                vm_long_reach=frac(d["c_vm_long_cost"][m] <= d["cost_true"][m] + TOL),
            )
        )
    s["by_radius"] = by_r
    return s


def mt(hdr: List[str], rows: List[List[Any]]) -> str:
    return markdown_table([dict(zip(hdr, r)) for r in rows], hdr)


def tables(S: Dict[str, Dict[str, Any]]) -> Dict[str, str]:
    t: Dict[str, str] = {}
    hdr = ["set", "n", "result above truth", "gap q25/50/75", "angle res-truth q25/50/75 (deg)"]
    rows = [
        [
            k, s["n"], fmt(s["result_above_truth"]),
            "/".join(f"{x:.3f}" for x in s["gap_quartiles"]),
            "/".join(f"{x:.3f}" for x in s["ang_res_true_quartiles"]),
        ]
        for k, s in S.items()
    ]  # fmt: skip
    t["t5_overview"] = mt(hdr, rows)
    hdr = ["set", "path cost rises above result", "barrier q25/50/75", "a point below result",
           "a point <= truth cost", "monotone"]  # fmt: skip
    rows = [
        [
            k, fmt(s["geo"]["barrier_gt0"]),
            "/".join(f"{x:.3f}" for x in s["geo"]["barrier_quartiles"]),
            fmt(s["geo"]["some_point_below_result"]),
            fmt(s["geo"]["min_point_below_truth_plus_tol"]),
            fmt(s["geo"]["monotone_nonincreasing"]),
        ]
        for k, s in S.items()
    ]  # fmt: skip
    t["t5_geodesic"] = mt(hdr, rows)
    hdr = ["set", "continuation", "improves", "reaches truth cost", "gap closed (median)",
           "toward truth q25/50/75 (deg)", "moved q50 (deg)"]  # fmt: skip
    rows = []
    for k, s in S.items():
        for kind, c in s["cont"].items():
            rows.append(
                [
                    k, kind, fmt(c["improves"]), fmt(c["reaches_truth_cost"]),
                    f"{c['gap_closed_median']:.2f}",
                    "/".join(f"{x:.3f}" for x in c["toward_truth_quartiles"]),
                    f"{c['moved_quartiles'][1]:.3f}",
                ]
            )  # fmt: skip
    t["t5_continuations"] = mt(hdr, rows)
    hdr = ["set", "angle (deg)", "truth: same cost", "truth: lower", "truth: median |d|",
           "result: same cost", "result: lower", "result: median |d|"]  # fmt: skip
    rows = []
    for k, s in S.items():
        for g in s["gran"]:
            rows.append(
                [
                    k, g["angle_deg"], f"{g['truth_same']:.2f}", f"{g['truth_lower']:.2f}",
                    f"{g['truth_med_abs']:.4f}", f"{g['res_same']:.2f}", f"{g['res_lower']:.2f}",
                    f"{g['res_med_abs']:.4f}",
                ]
            )  # fmt: skip
    t["t5_granularity"] = mt(hdr, rows)
    hdr = ["set", "FO stop codes (budget/restarts/cost)", "FO steps q50", "FO accepts q50",
           "FO last accept q50", "FO final step q50 (deg)", "VM steps q50", "VM improves",
           "VM supplies final", "VM final radius q50 (deg)"]  # fmt: skip
    rows = []
    for k, s in S.items():
        st = s["stop"]
        rows.append(
            [
                k, "/".join(str(v) for v in st["fo_codes"].values()),
                f"{st['fo_steps_run_quartiles'][1]:.0f}", f"{st['fo_n_accept_quartiles'][1]:.0f}",
                f"{st['fo_last_accept_quartiles'][1]:.0f}",
                f"{st['fo_final_step_deg_quartiles'][1]:.4f}",
                f"{st['vm_steps_taken_quartiles'][1]:.0f}", fmt(st["vm_improves"]),
                fmt(st["vm_gives_final"]), f"{st['vm_final_radius_deg_quartiles'][1]:.3f}",
            ]
        )  # fmt: skip
    t["t5_stopping"] = mt(hdr, rows)
    hdr = ["set", "radius (deg)", "n", "result above truth", "gap q50", "angle res-truth q50 (deg)",
           "mc_long reaches truth cost", "vm_long reaches truth cost"]  # fmt: skip
    rows = []
    for k, s in S.items():
        for r in s["by_radius"]:
            rows.append(
                [
                    k, r["radius_deg"], r["n"], f"{r['above']['k']}/{r['above']['n']}",
                    f"{r['gap_med']:.3f}", f"{r['ang_res_med']:.3f}",
                    f"{r['mc_long_reach']['k']}/{r['n']}", f"{r['vm_long_reach']['k']}/{r['n']}",
                ]
            )  # fmt: skip
    t["t5_by_radius"] = mt(hdr, rows)
    return t


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--cache", default=str(D.CACHE_DIR / "full"))
    ap.add_argument("--out", default=str(D.OUT_DIR))
    args = ap.parse_args()
    S = {}
    for pipe in ("H3", "H0"):
        S[pipe] = summarise(load(Path(args.cache), pipe))
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    (out / "summary.json").write_text(json.dumps(S, indent=1))
    write_tables(out / "tables.md", tables(S))
    print((out / "tables.md").read_text())


if __name__ == "__main__":
    main()
