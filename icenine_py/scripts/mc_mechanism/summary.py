#!/usr/bin/env python3
"""
Phase B1 summary: the MC restart / step-collapse mechanism and the improvement-probability curves.
Reads scripts/mc_mechanism/cache/<run|pilot> (mc_trace.py) and writes
benchmarks/mc_mechanism/{summary.json,tables.md}.

Usage (from icenine_py/):
  uv run python scripts/mc_mechanism/summary.py [--cache run]
"""

import argparse
import glob
import json
import math
import sys
from pathlib import Path
from typing import Any, Dict, List, Tuple

import numpy as np

HERE = Path(__file__).resolve().parent
ICENINE_PY = HERE.parents[1]
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent / "common"))

import mc_helpers as H  # noqa: E402
from doc_tables import markdown_table, write_tables  # noqa: E402
from stats import wilson  # noqa: E402

OUT_DIR = ICENINE_PY / "benchmarks" / "mc_mechanism"
T5_CACHE = ICENINE_PY / "scripts" / "finisher_diagnosis" / "cache" / "full"
VARIANTS = ["clean", "realistic"]
LOG_KEYS = [
    "stop",
    "steps_run",
    "n_accept",
    "last_accept",
    "n_restarts",
    "final_step_deg",
    "min_ergodic",
    "since_improve",
    "cost_start",
    "cost_end",
]
STOP_NAMES = {0: "step budget", 1: "restarts exhausted", 2: "cost converged"}
PT_LABEL = ["result", "truth", "0.005", "0.01", "0.02", "0.05"]
R50_DEG = 0.0005  # Phase A: plateau radius r50 (median, all four sets)
QUANT_DEG = 0.0125  # Phase A2: centroid-quantisation scale (not a bound)


def frac(k: int, n: int) -> Dict[str, Any]:
    lo, hi = wilson(int(k), int(n))
    return dict(k=int(k), n=int(n), p=k / max(n, 1), lo=lo, hi=hi)


def ff(f: Dict[str, Any]) -> str:
    return f"{f['k']}/{f['n']} ({100 * f['p']:.0f}%, {100 * f['lo']:.0f}-{100 * f['hi']:.0f})"


def q3(x: Any) -> List[float]:
    a = np.asarray(x, dtype=float)
    a = a[np.isfinite(a)]
    return [float(v) for v in np.quantile(a, [0.25, 0.5, 0.75])] if a.size else [math.nan] * 3


def qt(q: List[float], nd: int) -> str:
    return "/".join(f"{v:.{nd}f}" for v in q)


def load(cache: Path) -> Tuple[Dict[str, Dict[str, np.ndarray]], Dict[str, np.ndarray]]:
    """{pipe: {key: array (variant, N, ...)}} concatenated over tasks in sorted file order, plus
    the T5 fo_log of the same cases (realistic reference, shape (N, 10))."""
    files = sorted(glob.glob(str(cache / "*.npz")))
    acc: Dict[str, Dict[str, List[np.ndarray]]] = {}
    t5: Dict[str, List[np.ndarray]] = {"H3": [], "H0": []}
    t5s: Dict[str, List[np.ndarray]] = {"H3": [], "H0": []}
    t5v: Dict[str, List[np.ndarray]] = {"H3": [], "H0": []}
    meta: Dict[str, np.ndarray] = {}
    for f in files:
        d = np.load(f)
        pipe = str(d["pipe"])
        a = acc.setdefault(pipe, {})
        for k in d.files:
            if k in ("dirs", "vidx", "ri", "pipe", "steps_deg", "offsets_deg"):
                continue
            if k in ("box_deg", "step0_deg", "n_prop", "max_mc_steps", "max_restarts"):
                meta[k] = d[k]
                continue
            a.setdefault(k, []).append(d[k])
        meta["steps_deg"], meta["offsets_deg"] = d["steps_deg"], d["offsets_deg"]
        ref = T5_CACHE / Path(f).name
        if ref.exists():
            r5 = np.load(ref)
            t5[pipe].append(r5["fo_log"])
            t5s[pipe].append(r5["ang_start_true"])
            t5v[pipe].append(r5["vm_log"])
        a.setdefault("vidx_case", []).append(np.full(d["result"].shape[1], int(d["vidx"])))
    out = {
        p: {k: np.concatenate(v, axis=(0 if k == "vidx_case" else 1)) for k, v in a.items()}
        for p, a in acc.items()
    }
    meta2 = {f"t5_{p}": (np.concatenate(v) if v else np.zeros((0, 10))) for p, v in t5.items()}
    meta2.update(
        {f"t5vm_{p}": (np.concatenate(v) if v else np.zeros((0, 11))) for p, v in t5v.items()}
    )
    meta2.update(
        {f"t5start_{p}": (np.concatenate(v) if v else np.zeros(0)) for p, v in t5s.items()}
    )
    return out, {**meta, **meta2}


def load_end(cache: Path) -> Dict[str, Dict[str, np.ndarray]]:
    """mc_end_curve.py output: {pipe: {key: (variant, N, ...)}}, same file order as load()."""
    acc: Dict[str, Dict[str, List[np.ndarray]]] = {}
    for f in sorted(glob.glob(str(cache / "*.npz"))):
        d = np.load(f)
        a = acc.setdefault(str(d["pipe"]), {})
        for k in d.files:
            if k in ("dirs", "vidx", "ri", "pipe", "steps_deg", "n_prop"):
                continue
            a.setdefault(k, []).append(d[k])
        a.setdefault("steps_deg", [d["steps_deg"]])
        a.setdefault("n_prop", [d["n_prop"]])
    out: Dict[str, Dict[str, np.ndarray]] = {}
    for p, a in acc.items():
        out[p] = {
            k: (v[0] if k in ("steps_deg", "n_prop") else np.concatenate(v, axis=1))
            for k, v in a.items()
        }
    return out


def end_stats(vi: int, E: Dict[str, np.ndarray]) -> Dict[str, Any]:
    """Improvement probability at the MC stage's own output, for the runs that accepted at least
    once (their step has collapsed) and for those that never accepted."""
    s = E["steps_deg"]
    P, C = E["c_p_improve"][vi], E["c_cost_prog"][vi]  # (N, s)
    own, ownstep, dist, acc = (
        E["own_p_improve"][vi],
        E["own_step"][vi],
        E["dist_mc"][vi],
        E["n_accept"][vi],
    )
    n_prop = int(E["n_prop"])
    out: Dict[str, Any] = dict(steps_deg=s.tolist(), n_prop=n_prop)
    for name, m in (("improving", acc > 0), ("noaccept", acc == 0)):
        n = int(m.sum())
        if n == 0:
            continue
        Pm = P[m]
        best = Pm.max(axis=1)
        sP = np.array([s[int(np.argmax(r))] if r.max() > 0 else math.nan for r in Pm])
        cm = C[m].mean(axis=0)
        sC = H.argmax_step(s, cm)
        stuck = (own[m] == 0) & (best >= 0.05)
        out[name] = dict(
            n=n,
            own_step_q=q3(ownstep[m]),
            P_own_q=q3(own[m]),
            P_own_zero=frac(int((own[m] == 0).sum()), n),
            P_best_q=q3(best),
            s_best_P_q=q3(sP),
            s_best_P_over_dist_q=q3(sP / dist[m]),
            stuck_zero_own_but_best_ge_0p05=frac(int(stuck.sum()), n),
            dist_mc_q=q3(dist[m]),
            s_star_cost_mean=sC,
            P_med=np.median(Pm, axis=0).tolist(),
            P_q25=np.quantile(Pm, 0.25, axis=0).tolist(),
            P_q75=np.quantile(Pm, 0.75, axis=0).tolist(),
            cost_prog_mean=cm.tolist(),
            P_zero_everywhere=frac(int((best == 0).sum()), n),
        )
    return out


def mechanism(
    pipe: str, vi: int, D: Dict[str, np.ndarray], meta: Dict[str, np.ndarray]
) -> Dict[str, Any]:
    """Mechanism checks and run-end statistics for one set and variant."""
    log = D["log"][vi]  # (N, 10)
    N = len(log)
    L = {k: log[:, i] for i, k in enumerate(LOG_KEYS)}
    ev, step, merg, nsince = (
        D["tr_event"][vi],
        D["tr_step_deg"][vi],
        D["tr_min_erg"][vi],
        D["tr_n_since"][vi],
    )
    bang = D["tr_best_ang"][vi]
    max_steps = int(meta["max_mc_steps"])
    acc0 = L["n_accept"] == 0
    out: Dict[str, Any] = dict(
        n=N,
        box_deg=float(meta["box_deg"]),
        step0_deg=float(meta["step0_deg"]),
        max_mc_steps=max_steps,
        max_restarts=int(meta["max_restarts"]),
    )
    # --- (1) restart threshold: every restart / exhaustion fires with n_since == min_ergodic
    rs = (ev == H.EV_RESTART) | (ev == H.EV_EXHAUSTED)
    out["restart_events"] = int(rs.sum())
    out["restart_nsince_eq_minerg"] = frac(int((nsince[rs] == merg[rs]).sum()), int(rs.sum()))
    out["min_erg_initial"] = float(np.median(merg[:, 0]))
    out["min_erg_formula_initial"] = H.min_ergodic(
        math.radians(float(meta["box_deg"])), math.radians(float(meta["step0_deg"])), max_steps
    )
    # --- (2) each global improvement halves the step and multiplies min_ergodic by ~8
    g = ev == H.EV_GLOBAL
    g[:, -1] = False
    idx = np.argwhere(g)
    nxt = np.isfinite(step[idx[:, 0], idx[:, 1] + 1])
    ratio = step[idx[nxt, 0], idx[nxt, 1] + 1] / step[idx[nxt, 0], idx[nxt, 1]]
    merat = merg[idx[nxt, 0], idx[nxt, 1] + 1] / merg[idx[nxt, 0], idx[nxt, 1]]
    out["improvements_checked"] = int(nxt.sum())
    out["halving_exact"] = frac(int((ratio == 0.5).sum()), int(len(ratio)))
    out["minerg_ratio_q"] = q3(merat)
    # --- (3) stop reason vs accepts
    out["stop_counts"] = {nm: int((L["stop"] == c).sum()) for c, nm in STOP_NAMES.items()}
    out["noaccept"] = frac(int(acc0.sum()), N)
    out["noaccept_exhaust"] = frac(int((acc0 & (L["stop"] == 1)).sum()), int(acc0.sum()))
    out["exhaust_noaccept"] = frac(
        int((acc0 & (L["stop"] == 1)).sum()), int((L["stop"] == 1).sum())
    )
    out["accept_budget"] = frac(int((~acc0 & (L["stop"] == 0)).sum()), int((~acc0).sum()))
    out["noaccept_steps_run_q"] = q3(L["steps_run"][acc0])
    out["noaccept_restarts_q"] = q3(L["n_restarts"][acc0])
    # interval between restart events in no-accept runs (steps), and local accepts inside them
    iv: List[float] = []
    loc: List[float] = []
    for i in np.nonzero(acc0)[0]:
        e = np.nonzero(rs[i])[0]
        prev = -1
        for t in e:
            iv.append(t - prev)
            loc.append(int((ev[i, prev + 1 : t + 1] == H.EV_LOCAL).sum()))
            prev = t
    out["noaccept_interval_q"] = q3(iv)
    out["noaccept_local_accepts_per_interval_q"] = q3(loc)
    # --- any restart after the first global improvement?
    first_g = np.where(g.any(axis=1), np.argmax(g, axis=1), max_steps)
    after = [
        (ev[i, first_g[i] + 1 :] == H.EV_RESTART).any()
        | (ev[i, first_g[i] + 1 :] == H.EV_EXHAUSTED).any()
        for i in range(N)
        if first_g[i] < max_steps
    ]
    out["restart_after_first_accept"] = frac(int(np.sum(after)), len(after))
    out["minerg_after_first_accept_q"] = q3(
        [
            merg[i, first_g[i] + 1]
            for i in range(N)
            if first_g[i] + 1 < max_steps and np.isfinite(merg[i, first_g[i] + 1])
        ]
    )
    # --- (4) improving runs: where the run ends
    imp = ~acc0
    Li = {k: v[imp] for k, v in L.items()}
    ii = np.nonzero(imp)[0]
    ang_end = np.array(
        [bang[i, int(L["steps_run"][i]) - 1] for i in ii]
    )  # distance of the MC output
    res_dist = D["ang_res_true"][vi][
        imp
    ]  # distance of the finisher result (after VarianceMinimizing)
    fs = Li["final_step_deg"]
    tail = Li["steps_run"] - Li["last_accept"] - 1
    out["n_improving"] = int(imp.sum())
    out["n_accept_q"] = q3(Li["n_accept"])
    out["last_accept_q"] = q3(Li["last_accept"])
    out["tail_q"] = q3(tail)
    out["tail_frac_q"] = q3(tail / Li["steps_run"])
    out["final_step_q"] = q3(fs)
    out["final_step_over_dist_q"] = q3(fs / ang_end)
    out["dist_end_q"] = q3(ang_end)
    out["result_dist_q"] = q3(res_dist)
    out["final_step_lt_r50"] = frac(int((fs < R50_DEG).sum()), int(imp.sum()))
    out["final_step_lt_quant"] = frac(int((fs < QUANT_DEG).sum()), int(imp.sum()))
    out["final_step_lt_dist"] = frac(int((fs < ang_end).sum()), int(imp.sum()))
    out["final_step_lt_dist_over_10"] = frac(int((fs < ang_end / 10).sum()), int(imp.sum()))
    # --- travel ceiling: the accepted moves of a run that never restarts again form a geometric
    # series (each trial rotation is at most the current step, which halves at each accept)
    two_s0 = 2.0 * H.corner_angle_ratio() * float(meta["step0_deg"])  # ceiling on the path length
    ta, ev_i = D["tr_trial_ang"][vi][imp], ev[imp]
    path = np.array([np.nansum(ta[i][ev_i[i] == H.EV_GLOBAL]) for i in range(int(imp.sum()))])
    out["two_step0_deg"] = two_s0
    out["corner_ratio"] = H.corner_angle_ratio()
    out["path_q"] = q3(path)
    out["path_max"] = float(path.max())
    out["path_lt_two_step0"] = frac(int((path < two_s0).sum()), int(imp.sum()))
    out["improving_with_restart"] = frac(int((Li["n_restarts"] > 0).sum()), int(imp.sum()))
    sa = meta.get(f"t5start_{pipe}")
    if sa is not None and len(sa) == N:
        far = sa > two_s0
        out["start_err_q"] = q3(sa)
        out["start_beyond_two_step0"] = frac(int(far.sum()), N)
        out["final_err_q_start_beyond"] = q3(D["ang_res_true"][vi][far])
        out["final_err_q_start_within"] = q3(D["ang_res_true"][vi][~far])
        out["n_start_within"] = int((~far).sum())
    # --- lock-in: after the k-th accept (step s after halving) the run can still move at most
    # 2 * corner * s in total, so once the distance to the truth exceeds that, the truth is out
    # of reach of the rest of the run (triangle inequality; no restart can fire after an accept)
    reach_f = 2.0 * H.corner_angle_ratio()
    lk: List[Dict[str, float]] = []
    for i in np.nonzero(imp)[0]:
        for k, t in enumerate(np.nonzero(ev[i] == H.EV_GLOBAL)[0], start=1):
            if bang[i, t] > reach_f * step[i, t] * 0.5:
                lk.append(
                    dict(
                        i=int(i),
                        k=k,
                        at=float(t),
                        dist=float(bang[i, t]),
                        step=float(step[i, t] * 0.5),
                        final=float(bang[i, int(L["steps_run"][i]) - 1]),
                    )
                )
                break
    out["lockin"] = dict(
        n_locked=frac(len(lk), int(imp.sum())),
        k_q=q3([x["k"] for x in lk]),
        at_q=q3([x["at"] for x in lk]),
        step_q=q3([x["step"] for x in lk]),
        dist_q=q3([x["dist"] for x in lk]),
        final_ge_floor=frac(
            int(sum(x["final"] >= x["dist"] - reach_f * x["step"] - 1e-9 for x in lk)), len(lk)
        ),
        reach_factor=reach_f,
    )
    # --- by accept number: step and distance to the truth after the k-th global improvement
    byk = []
    for k in range(1, 11):
        st, an, at = [], [], []
        for i in np.nonzero(imp)[0]:
            e = np.nonzero(g[i] | (ev[i] == H.EV_GLOBAL))[0]
            if len(e) >= k:
                t = e[k - 1]
                at.append(t)
                an.append(bang[i, t])
                st.append(step[i, t] * 0.5)  # step after halving
        if st:
            byk.append(
                dict(
                    k=k,
                    n=len(st),
                    step=float(np.median(st)),
                    dist=float(np.median(an)),
                    at=float(np.median(at)),
                )
            )
    out["by_accept"] = byk
    return out


def curve_stats(
    pipe: str, vi: int, D: Dict[str, np.ndarray], meta: Dict[str, np.ndarray]
) -> Dict[str, Any]:
    s = np.asarray(meta["steps_deg"], float)
    P = D["pc_p_improve"][vi]  # (N, pts, s)
    C = D["pc_cost_prog"][vi]
    G = D["pc_dist_prog"][vi]
    d0 = D["pc_dist0"][vi]  # (N, pts)
    N = P.shape[0]
    log = D["log"][vi]
    imp = log[:, 2] > 0
    fstep = log[:, 5]
    groups = {
        "result": [0],
        "truth": [1],
        "0.005": [2, 3],
        "0.01": [4, 5],
        "0.02": [6, 7],
        "0.05": [8, 9],
    }
    out: Dict[str, Any] = dict(steps_deg=s.tolist(), points={})
    for name, ix in groups.items():
        Pm = P[:, ix].mean(axis=1)  # (N, s): per case, mean over the point's directions
        Cm = C[:, ix].mean(axis=1)
        Gm = G[:, ix].mean(axis=1)
        dm = d0[:, ix].mean(axis=1)
        pm = np.median(Pm, axis=0)
        cm = Cm.mean(axis=0)
        gm = Gm.mean(axis=0)
        sP, sC, sG = H.argmax_step(s, pm), H.argmax_step(s, cm), H.argmax_step(s, gm)
        per_case_sC = np.array([H.argmax_step(s, Cm[i]) or math.nan for i in range(N)])
        ent = dict(
            d0_med=float(np.median(dm)),
            P_med=pm.tolist(),
            P_q25=np.quantile(Pm, 0.25, axis=0).tolist(),
            P_q75=np.quantile(Pm, 0.75, axis=0).tolist(),
            cost_prog_mean=cm.tolist(),
            dist_prog_mean=gm.tolist(),
            s_star_P=sP,
            s_star_cost=sC,
            s_star_dist=sG,
            s_star_cost_percase_q=q3(per_case_sC),
            P_at_s_star_cost=float(pm[list(s).index(sC)]) if sC is not None else math.nan,
            s_star_over_d0=(
                (sC / float(np.median(dm))) if (sC is not None and np.median(dm) > 0) else math.nan
            ),
        )
        if name == "truth":
            ni = D["pc_n_improve"][vi][:, 1, :]  # (N, s) improving proposals out of N_PROP
            ent.update(
                any_improvement=frac(int((ni.sum(axis=1) > 0).sum()), N),
                P_pooled_max_over_s=float(ni.sum(axis=0).max() / (N * float(meta["n_prop"]))),
            )
        if name == "result":
            # the finisher's own final step against the curve at its result (improving runs)
            Pf = np.array(
                [np.interp(np.log(max(fstep[i], 1e-12)), np.log(s), Pm[i]) for i in range(N)]
            )
            Cf = np.array(
                [np.interp(np.log(max(fstep[i], 1e-12)), np.log(s), Cm[i]) for i in range(N)]
            )
            clamp = fstep < s[0]
            vm = meta.get(f"t5vm_{pipe}")
            if vm is not None and len(vm) == N and vi == 1:
                vr = vm[:, 3]  # final VarianceMinimizing subregion radius (deg), a step-equivalent
                Pv = np.array(
                    [np.interp(np.log(max(vr[i], 1e-12)), np.log(s), Pm[i]) for i in range(N)]
                )
                ent.update(
                    vm_final_radius_q=q3(vr),
                    P_at_vm_radius_q=q3(Pv),
                    vm_n_improve_q=q3(vm[:, 1]),
                    vm_steps_q=q3(vm[:, 6]),
                    vm_radius_clamped=frac(int(((vr < s[0]) | (vr > s[-1])).sum()), N),
                )
            ent.update(
                final_step_q_improving=q3(fstep[imp]),
                P_at_final_step_q_improving=q3(Pf[imp]),
                clamped_improving=frac(int(clamp[imp].sum()), int(imp.sum())),
                P_at_initial_step_q_noaccept=q3(
                    [
                        np.interp(np.log(float(meta["step0_deg"])), np.log(s), Pm[i])
                        for i in np.nonzero(~imp)[0]
                    ]
                ),
                cost_prog_ratio_best_over_final=q3(
                    [Cm[i].max() / Cf[i] if Cf[i] > 0 else math.nan for i in np.nonzero(imp)[0]]
                ),
                evals_per_improvement_at_final_q=q3(1.0 / np.maximum(Pf[imp], 1e-3)),
                s_star_cost_over_final_step=(sC / float(np.median(fstep[imp]))) if sC else math.nan,
                P_best_q=q3(Pm.max(axis=1)),
                P_zero_everywhere=frac(int((Pm.max(axis=1) == 0).sum()), N),
            )
        out["points"][name] = ent
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--cache", default="run")
    args = ap.parse_args()
    data, meta = load(HERE / "cache" / args.cache)
    S: Dict[str, Any] = dict(cache=args.cache)
    # --- reproduction of T5
    rep = {}
    for pipe, D in data.items():
        real = D["rerun_maxabs"][1]
        t5 = meta[f"t5_{pipe}"]
        same = (
            bool(len(t5) == len(real) and np.array_equal(D["log"][1], t5, equal_nan=True))
            if len(t5)
            else None
        )
        rep[pipe] = dict(
            n=int(real.size),
            realistic_result_identical=frac(int((real == 0).sum()), int(real.size)),
            log_equals_t5=same,
            clean_result_identical_to_realistic_result=frac(
                int((D["rerun_maxabs"][0] == 0).sum()), int(real.size)
            ),
        )
    S["reproduction"] = rep
    assert all(r["realistic_result_identical"]["k"] == r["n"] for r in rep.values())
    S["setup"] = {
        k: float(v)
        for k, v in meta.items()
        if k in ("box_deg", "step0_deg", "n_prop", "max_mc_steps", "max_restarts")
    }
    S["finisher_evals_q"] = {
        f"{p}|{vn}": q3(D["rerun_evals"][vi])
        for p, D in data.items()
        for vi, vn in enumerate(VARIANTS)
    }  # evaluations of the whole finisher (MC + VarianceMinimizing) per case
    S["mechanism"], S["curves"] = {}, {}
    for pipe, D in data.items():
        for vi, vn in enumerate(VARIANTS):
            S["mechanism"][f"{pipe}|{vn}"] = mechanism(pipe, vi, D, meta)
            S["curves"][f"{pipe}|{vn}"] = curve_stats(pipe, vi, D, meta)
    import re

    comp: Dict[str, Any] = dict(
        n_prop=int(meta["n_prop"]), n_points=10, n_steps=len(meta["steps_deg"])
    )
    comp["curve_evals_per_case_variant"] = comp["n_prop"] * comp["n_points"] * comp["n_steps"]
    for tag, fn in (("trace", "trace.log"), ("end", "end.log")):
        lp = HERE / "cache" / "logs" / fn
        m = re.search(r"finished; wall (\d+)s", lp.read_text()) if lp.exists() else None
        comp[f"wall_{tag}_s_contended_10_workers"] = int(m.group(1)) if m else None
    comp["finisher_evals_median_all_cases"] = float(
        np.median(np.concatenate([D["rerun_evals"][1] for D in data.values()]))
    )
    S["compute"] = comp
    end = load_end(HERE / "cache" / ("end_" + args.cache))
    S["end_curve"] = {
        f"{p}|{vn}": end_stats(vi, E) for p, E in end.items() for vi, vn in enumerate(VARIANTS)
    }
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    (OUT_DIR / "summary.json").write_text(json.dumps(S, indent=1, default=float))
    write_tables(OUT_DIR / "tables.md", tables(S))
    print((OUT_DIR / "tables.md").read_text())


def tables(S: Dict[str, Any]) -> Dict[str, str]:
    T: Dict[str, str] = {}
    keys = list(S["mechanism"])
    ec = S.get("end_curve", {})

    def sv(k: str) -> Tuple[str, str]:
        pipe, vn = k.split("|")
        return pipe, vn

    # --- the rules
    rows = []
    for k in keys:
        m, (pipe, vn) = S["mechanism"][k], sv(k)
        rows.append(
            {
                "set": pipe,
                "variant": vn,
                "n": m["n"],
                "box_deg": m["box_deg"],
                "step0_deg": m["step0_deg"],
                "min_erg_0": m["min_erg_initial"],
                "formula_0": m["min_erg_formula_initial"],
                "restarts at n_since = min_erg": ff(m["restart_nsince_eq_minerg"]),
                "step halved at each improvement": ff(m["halving_exact"]),
                "min_erg ratio q25/50/75": qt(m["minerg_ratio_q"], 2),
                "restart after first accept": ff(m["restart_after_first_accept"]),
                "min_erg after first accept": qt(m["minerg_after_first_accept_q"], 0),
            }
        )
    T["b1_rules"] = markdown_table(
        rows,
        list(rows[0]),
        {"box_deg": ".4f", "step0_deg": ".4f", "min_erg_0": ".0f", "formula_0": "d"},
    )
    # --- stop reasons
    rows = []
    for k in keys:
        m, (pipe, vn) = S["mechanism"][k], sv(k)
        sc = m["stop_counts"]
        rows.append(
            {
                "set": pipe,
                "variant": vn,
                "n": m["n"],
                "step budget": sc["step budget"],
                "restarts exhausted": sc["restarts exhausted"],
                "converged": sc["cost converged"],
                "no accept": ff(m["noaccept"]),
                "exhausted given no accept": ff(m["noaccept_exhaust"]),
                "no accept given exhausted": ff(m["exhaust_noaccept"]),
                "step budget given accepted": ff(m["accept_budget"]),
                "steps run (no accept)": qt(m["noaccept_steps_run_q"], 0),
                "restarts (no accept)": qt(m["noaccept_restarts_q"], 0),
                "steps between restarts": qt(m["noaccept_interval_q"], 0),
                "local accepts in between": qt(m["noaccept_local_accepts_per_interval_q"], 0),
            }
        )
    T["b1_stop"] = markdown_table(rows, list(rows[0]))
    # --- where improving runs end
    rows = []
    for k in keys:
        m, (pipe, vn) = S["mechanism"][k], sv(k)
        rows.append(
            {
                "set": pipe,
                "variant": vn,
                "n": m["n_improving"],
                "accepts": qt(m["n_accept_q"], 0),
                "last accept step": qt(m["last_accept_q"], 0),
                "steps after last accept": qt(m["tail_q"], 0),
                "final step (deg)": qt(m["final_step_q"], 5),
                "MC output to truth (deg)": qt(m["dist_end_q"], 4),
                "finisher result to truth (deg)": qt(m["result_dist_q"], 4),
                "step / MC distance": qt(m["final_step_over_dist_q"], 4),
                "step < 0.0005": ff(m["final_step_lt_r50"]),
                "step < 0.0125": ff(m["final_step_lt_quant"]),
                "step < distance": ff(m["final_step_lt_dist"]),
                "step < distance/10": ff(m["final_step_lt_dist_over_10"]),
            }
        )
    T["b1_end"] = markdown_table(rows, list(rows[0]))
    # --- lock-in
    rows = []
    for k in keys:
        m, (pipe, vn) = S["mechanism"][k], sv(k)
        lk = m["lockin"]
        rows.append(
            {
                "set": pipe,
                "variant": vn,
                "n": m["n_improving"],
                "locked in": ff(lk["n_locked"]),
                "at accept k": qt(lk["k_q"], 0),
                "at run step": qt(lk["at_q"], 0),
                "step then (deg)": qt(lk["step_q"], 5),
                "distance then (deg)": qt(lk["dist_q"], 4),
                "MC output >= floor": ff(lk["final_ge_floor"]),
            }
        )
    T["b1_lockin"] = markdown_table(rows, list(rows[0]))
    # --- travel
    rows = []
    for k in keys:
        m, (pipe, vn) = S["mechanism"][k], sv(k)
        nan3 = [math.nan] * 3
        rows.append(
            {
                "set": pipe,
                "variant": vn,
                "n": m["n_improving"],
                "ceiling (deg)": m["two_step0_deg"],
                "path length q25/50/75": qt(m["path_q"], 3),
                "path max": m["path_max"],
                "path < ceiling": ff(m["path_lt_two_step0"]),
                "restart before first accept": ff(m["improving_with_restart"]),
                "start error": qt(m.get("start_err_q", nan3), 3),
                "start > ceiling": (
                    ff(m["start_beyond_two_step0"]) if "start_beyond_two_step0" in m else "n/a"
                ),
                "result error given start > ceiling": qt(
                    m.get("final_err_q_start_beyond", nan3), 3
                ),
                "result error given start <= ceiling": qt(
                    m.get("final_err_q_start_within", nan3), 3
                ),
            }
        )
    T["b1_travel"] = markdown_table(
        rows, list(rows[0]), {"ceiling (deg)": ".4f", "path max": ".4f"}
    )
    # --- by accept number
    rows = []
    for k in keys:
        pipe, vn = sv(k)
        for b in S["mechanism"][k]["by_accept"]:
            rows.append(
                {
                    "set": pipe,
                    "variant": vn,
                    "k": b["k"],
                    "n": b["n"],
                    "run step": b["at"],
                    "step after (deg)": b["step"],
                    "distance to truth (deg)": b["dist"],
                }
            )
    fm_ba = {"run step": ".0f", "step after (deg)": ".5f", "distance to truth (deg)": ".4f"}
    T["b1_by_accept"] = markdown_table(rows, list(rows[0]), fm_ba)
    T["b1_by_accept_h3"] = markdown_table(
        [
            {k: v for k, v in r.items() if k not in ("set", "variant")}
            for r in rows
            if r["set"] == "H3" and r["variant"] == "realistic"
        ],
        ["k", "n", "run step", "step after (deg)", "distance to truth (deg)"],
        fm_ba,
    )
    # --- curves
    s = S["curves"][keys[0]]["steps_deg"]
    cols = ["set", "variant", "point", "d0"] + [f"{x:g}" for x in s]
    for kind, nm, nd in (("P_med", "b1_P", 2), ("cost_prog_mean", "b1_cost_prog", 4)):
        rows = []
        for k in keys:
            pipe, vn = sv(k)
            for name, e in S["curves"][k]["points"].items():
                r = dict(set=pipe, variant=vn, point=name, d0=e["d0_med"])
                for x, v in zip(s, e[kind]):
                    r[f"{x:g}"] = v
                rows.append(r)
        fm = {"d0": ".4f", **{f"{x:g}": f".{nd}f" for x in s}}
        T[nm] = markdown_table(rows, cols, fm)
        T[nm + "_h3"] = markdown_table([r for r in rows if r["set"] == "H3"], cols, fm)
    rows = []
    for k in keys:
        pipe, vn = sv(k)
        for name, e in S["curves"][k]["points"].items():
            rows.append(
                {
                    "set": pipe,
                    "variant": vn,
                    "point": name,
                    "d0 (deg)": e["d0_med"],
                    "s max P": e["s_star_P"],
                    "s max cost progress": e["s_star_cost"],
                    "s max distance progress": e["s_star_dist"],
                    "per-case s max cost progress q25/50/75": qt(e["s_star_cost_percase_q"], 4),
                    "P at s max cost progress": e["P_at_s_star_cost"],
                    "s / d0": e["s_star_over_d0"],
                }
            )
    T["b1_sstar"] = markdown_table(
        rows,
        list(rows[0]),
        {
            "d0 (deg)": ".4f",
            "s max P": ".4f",
            "s max cost progress": ".4f",
            "s max distance progress": ".4f",
            "P at s max cost progress": ".2f",
            "s / d0": ".2f",
        },
    )
    rows = []
    for k in keys:
        pipe, vn = sv(k)
        e = S["curves"][k]["points"]["result"]
        rows.append(
            {
                "set": pipe,
                "variant": vn,
                "d0 (deg)": e["d0_med"],
                "s max cost progress": e["s_star_cost"],
                "s max cost progress / MC final step": e["s_star_cost_over_final_step"],
                "MC final step (improving runs)": qt(e["final_step_q_improving"], 5),
                "P at MC final step": qt(e["P_at_final_step_q_improving"], 3),
                "step below grid": ff(e["clamped_improving"]),
                "evals per improvement at final step": qt(e["evals_per_improvement_at_final_q"], 0),
                "P best on grid": qt(e["P_best_q"], 2),
                "cost progress best / final": (
                    qt(e["cost_progress_ratio_best_over_final"], 1)
                    if "cost_progress_ratio_best_over_final" in e
                    else qt(e["cost_prog_ratio_best_over_final"], 1)
                ),
                "VarianceMinimizing final radius": (
                    qt(e["vm_final_radius_q"], 4) if "vm_final_radius_q" in e else "n/a"
                ),
                "P at that radius": (
                    qt(e["P_at_vm_radius_q"], 2) if "P_at_vm_radius_q" in e else "n/a"
                ),
                "P at initial step (no-accept runs)": qt(e["P_at_initial_step_q_noaccept"], 2),
                "P = 0 at every step": ff(e["P_zero_everywhere"]),
            }
        )
    T["b1_final_vs_best"] = markdown_table(
        rows,
        list(rows[0]),
        {
            "d0 (deg)": ".4f",
            "s max cost progress": ".4f",
            "s max cost progress / MC final step": ".1f",
        },
    )
    # --- truth
    rows = []
    for k in keys:
        pipe, vn = sv(k)
        e = S["curves"][k]["points"]["truth"]
        rows.append(
            {
                "set": pipe,
                "variant": vn,
                "cases with an improving proposal at the truth": ff(e["any_improvement"]),
                "pooled P at the best step": e["P_pooled_max_over_s"],
            }
        )
    T["b1_truth"] = markdown_table(rows, list(rows[0]), {"pooled P at the best step": ".3f"})
    # --- at the MC output
    rows = []
    for k, e in ec.items():
        pipe, vn = sv(k)
        for grp in ("improving", "noaccept"):
            if grp not in e:
                continue
            g = e[grp]
            rows.append(
                {
                    "set": pipe,
                    "variant": vn,
                    "runs": grp,
                    "n": g["n"],
                    "own step (deg)": qt(g["own_step_q"], 5),
                    "output to truth (deg)": qt(g["dist_mc_q"], 4),
                    "P at own step": qt(g["P_own_q"], 3),
                    "P = 0 at own step": ff(g["P_own_zero"]),
                    "best P on grid": qt(g["P_best_q"], 2),
                    "step of best P": qt(g["s_best_P_q"], 4),
                    "P = 0 at own step, best P >= 0.05": ff(g["stuck_zero_own_but_best_ge_0p05"]),
                    "s max cost progress (mean)": g["s_star_cost_mean"],
                }
            )
    if rows:
        T["b1_end_curve_stats"] = markdown_table(
            rows, list(rows[0]), {"s max cost progress (mean)": ".4f"}
        )
        s_end = next(iter(ec.values()))["steps_deg"]
        rows = []
        for k, e in ec.items():
            pipe, vn = sv(k)
            if "improving" not in e:
                continue
            r = {"set": pipe, "variant": vn, "n": e["improving"]["n"]}
            for x, v in zip(s_end, e["improving"]["P_med"]):
                r[f"{x:g}"] = v
            rows.append(r)
        T["b1_end_curve"] = markdown_table(rows, list(rows[0]), {f"{x:g}": ".3f" for x in s_end})
    return T


if __name__ == "__main__":
    main()
