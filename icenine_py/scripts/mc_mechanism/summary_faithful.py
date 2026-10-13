#!/usr/bin/env python3
"""B1 rerun with the C++-faithful MC and VarianceMinimizing: the headline mechanism numbers,
before (cache/run, cache/end_run) and after (cache/run_mcfaithful, cache/end_run_mcfaithful).

  uv run python scripts/mc_mechanism/mc_trace.py run --tag mcfaithful --workers 10
  uv run python scripts/mc_mechanism/mc_end_curve.py run --tag mcfaithful --workers 10
  uv run python scripts/mc_mechanism/summary_faithful.py

Writes benchmarks/mc_mechanism/mc_faithful/{summary.json,tables.md}. Trace semantics after the
port (mc_trace.TracedMC): an event is recorded per step; the block outcome sits on the block's
last step (EV_GLOBAL success, EV_RESTART / EV_EXHAUSTED failure).
"""

import glob
import json
import sys
from pathlib import Path
from typing import Any, Dict, List, Tuple

import numpy as np

HERE = Path(__file__).resolve().parent
ICE = HERE.parents[1]
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent / "common"))
import mc_helpers as H  # noqa: E402
from doc_tables import markdown_table, write_tables  # noqa: E402
from stats import wilson  # noqa: E402

OUT = ICE / "benchmarks" / "mc_mechanism" / "mc_faithful"
CACHE = HERE / "cache"
VARIANTS = ["clean", "realistic"]
LOG = ["stop", "steps_run", "n_accept", "last_accept", "n_restarts", "final_step_deg",
       "min_ergodic", "since_improve", "cost_start", "cost_end"]  # fmt: skip


def load(cache: Path) -> Dict[str, Dict[str, np.ndarray]]:
    out: Dict[str, Dict[str, np.ndarray]] = {}
    for f in sorted(glob.glob(str(cache / "*.npz"))):
        z = np.load(f)
        pipe = str(z["pipe"])
        d = out.setdefault(pipe, {})
        for k in z.files:
            if z[k].ndim >= 2 and k not in ("steps_deg", "offsets_deg"):
                d.setdefault(k, []).append(z[k])  # type: ignore[arg-type]
    return {
        p: {k: np.concatenate(v, axis=1) for k, v in d.items()} for p, d in out.items()
    }  # (variant, N, ...)


def frac(k: int, n: int) -> str:
    lo, hi = wilson(int(k), int(n))
    return f"{k}/{n} ({100 * k / max(n, 1):.0f}%, {100 * lo:.0f}-{100 * hi:.0f})"


def q3(x: np.ndarray, f: str = ".4f") -> str:
    q = np.nanquantile(x, [0.25, 0.5, 0.75])
    return "/".join(format(v, f) for v in q)


def mech(d: Dict[str, np.ndarray], vi: int) -> Dict[str, Any]:
    log = d["log"][vi]
    L = {k: log[:, i] for i, k in enumerate(LOG)}
    n = len(log)
    ev = d["tr_event"][vi]
    nr = ev != H.EV_NOT_RUN
    acc = L["n_accept"] > 0
    stop = L["stop"]
    # restart events after the run's first global improvement
    after_first = 0
    restarts_total = 0
    for i in range(n):
        g = np.nonzero(ev[i] == H.EV_GLOBAL)[0]
        r = np.nonzero((ev[i] == H.EV_RESTART) | (ev[i] == H.EV_EXHAUSTED))[0]
        restarts_total += len(r)
        if len(g):
            after_first += int((r > g[0]).sum())
    last = nr.sum(axis=1) - 1
    best_ang = d["tr_best_ang"][vi][np.arange(n), np.maximum(last, 0)]
    res_err = d["ang_res_true"][vi]
    out: Dict[str, Any] = dict(
        n=n,
        stop_budget=int((stop == 0).sum()),
        stop_restarts=int((stop == 1).sum()),
        stop_converged=int((stop == 2).sum()),
        any_improve=frac(int(acc.sum()), n),
        exhausted_given_noimp=frac(int((~acc & (stop == 1)).sum()), int((~acc).sum())),
        noimp_given_exhausted=frac(int((~acc & (stop == 1)).sum()), int((stop == 1).sum())),
        restart_after_first_improve=f"{after_first}/{restarts_total} restarts",
        runs_restart_after_first=frac(
            int(
                sum(
                    bool(
                        (
                            ((ev[i] == H.EV_RESTART) | (ev[i] == H.EV_EXHAUSTED))[
                                (np.nonzero(ev[i] == H.EV_GLOBAL)[0][0] + 1 if acc[i] else 0) :
                            ]
                        ).any()
                    )
                    for i in range(n)
                    if acc[i]
                )
            ),
            int(acc.sum()),
        ),
        steps_run=q3(L["steps_run"], ".0f"),
        n_accept=q3(L["n_accept"], ".0f"),
        n_restarts=q3(L["n_restarts"], ".0f"),
        final_step_deg=q3(L["final_step_deg"], ".5f"),
        mc_out_to_truth=q3(best_ang, ".4f"),
        mc_out_to_truth_median=float(np.nanmedian(best_ang)),
        finisher_err=q3(res_err, ".4f"),
        finisher_err_median=float(np.median(res_err)),
        finisher_lt002=frac(int((res_err < 0.02).sum()), n),
        finisher_wrong=frac(int((res_err > 1.000001).sum()), n),
        finisher_evals=q3(d["rerun_evals"][vi], ".0f"),
        finisher_evals_median=float(np.median(d["rerun_evals"][vi])),
    )  # fmt: skip
    return out


def curves(d: Dict[str, np.ndarray], steps: np.ndarray, vi: int) -> Dict[str, Any]:
    P = d["pc_p_improve"][vi]  # (N, 10, 15): points result, truth, 8 offsets
    out: Dict[str, Any] = {}
    for pi, name in ((0, "result"), (1, "truth")):
        med = np.nanmedian(P[:, pi, :], axis=0)
        out[name] = dict(P_med=[float(x) for x in med],
                         best_step_med=float(np.nanmedian([H.argmax_step(steps, p) or np.nan
                                                           for p in P[:, pi, :]])))  # fmt: skip
    return out


def main() -> None:
    old, new = load(CACHE / "run"), load(CACHE / "run_mcfaithful")
    steps = np.load(sorted(glob.glob(str(CACHE / "run" / "*.npz")))[0])["steps_deg"]
    S: Dict[str, Any] = dict(steps_deg=[float(x) for x in steps])
    rows: List[Dict[str, Any]] = []
    for pipe in ("H3", "H0"):
        for vi, var in enumerate(VARIANTS):
            for tag, dd in (("before", old), ("after", new)):
                m = mech(dd[pipe], vi)
                S[f"{pipe}|{var}|{tag}"] = m
                S[f"{pipe}|{var}|{tag}|curves"] = curves(dd[pipe], steps, vi)
                rows.append(dict(set=pipe, variant=var, mc=tag, **m))
    cols_stop = ["set", "variant", "mc", "n", "stop_budget", "stop_restarts", "stop_converged",
                 "any_improve", "exhausted_given_noimp", "noimp_given_exhausted",
                 "runs_restart_after_first"]  # fmt: skip
    cols_end = ["set", "variant", "mc", "steps_run", "n_accept", "n_restarts", "final_step_deg",
                "mc_out_to_truth", "finisher_err", "finisher_lt002", "finisher_wrong",
                "finisher_evals"]  # fmt: skip
    tables = {
        "mcf_b1_stop": markdown_table(rows, cols_stop),
        "mcf_b1_end": markdown_table(rows, cols_end),
    }
    # improvement-probability curve at the finisher's result and at the truth, H3 realistic
    crow = []
    for pipe in ("H3", "H0"):
        for tag in ("before", "after"):
            c = S[f"{pipe}|realistic|{tag}|curves"]
            for pt in ("result", "truth"):
                pm = c[pt]["P_med"]
                crow.append(dict(set=pipe, mc=tag, point=pt,
                                 **{f"s{steps[i]:g}": pm[i] for i in (1, 2, 4, 7, 9, 11)},
                                 best_step=c[pt]["best_step_med"]))  # fmt: skip
    tables["mcf_b1_curve"] = markdown_table(
        crow,
        ["set", "mc", "point"] + [f"s{steps[i]:g}" for i in (1, 2, 4, 7, 9, 11)] + ["best_step"],
        formats={f"s{steps[i]:g}": ".2f" for i in (1, 2, 4, 7, 9, 11)} | {"best_step": ".4f"},
    )  # fmt: skip
    # the same curve at the MC stage's own output
    eo = glob.glob(str(CACHE / "end_run_mcfaithful" / "*.npz"))
    if eo:
        e_old, e_new = load(CACHE / "end_run"), load(CACHE / "end_run_mcfaithful")
        t_old, t_new = old, new  # traces: the step the last block ran at
        gsteps = np.load(sorted(glob.glob(str(CACHE / "end_run" / "*.npz")))[0])["steps_deg"]
        erow = []
        for pipe in ("H3", "H0"):
            for tag, dd, tt in (("before", e_old, t_old), ("after", e_new, t_new)):
                d = dd[pipe]
                vi = 1
                imp = d["n_accept"][vi] > 0
                own = d["own_p_improve"][vi][imp]
                P = d["c_p_improve"][vi][imp]
                tr = tt[pipe]
                nr = tr["tr_event"][vi] != H.EV_NOT_RUN
                last_i = nr.sum(axis=1) - 1
                last_step = tr["tr_step_deg"][vi][np.arange(len(last_i)), np.maximum(last_i, 0)]
                gl = np.log(np.asarray(gsteps))
                nearest = np.argmin(np.abs(gl[None, :] - np.log(last_step)[:, None]), axis=1)
                P_last = d["c_p_improve"][vi][np.arange(len(nearest)), nearest][imp]
                erow.append(dict(
                    set=pipe, mc=tag, n_improving=int(imp.sum()),
                    mc_out_to_truth=float(np.median(d["dist_mc"][vi][imp])),
                    P_at_own_step=float(np.median(own)),
                    last_block_step=float(np.median(last_step[imp])),
                    P_at_last_block_step=float(np.median(P_last)),
                    best_step=float(np.nanmedian([H.argmax_step(gsteps, p) or np.nan for p in P])),
                    P_at_0p0005=float(np.median(P[:, list(gsteps).index(0.0005)])),
                    P_at_0p005=float(np.median(P[:, list(gsteps).index(0.005)])),
                ))  # fmt: skip
        S["end_curve"] = erow
        tables["mcf_b1_endcurve"] = markdown_table(
            erow,
            ["set", "mc", "n_improving", "mc_out_to_truth", "last_block_step",
             "P_at_last_block_step", "P_at_own_step", "best_step",
             "P_at_0p0005", "P_at_0p005"],
            formats={"mc_out_to_truth": ".4f", "P_at_own_step": ".2f", "best_step": ".4f",
                     "last_block_step": ".5f", "P_at_last_block_step": ".2f",
                     "P_at_0p0005": ".2f", "P_at_0p005": ".2f"},
        )  # fmt: skip
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "summary.json").write_text(json.dumps(S, indent=1, default=float))
    write_tables(OUT / "tables.md", tables)
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
