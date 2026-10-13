#!/usr/bin/env python3
"""
Phase B2 (finisher/MC study): did we deploy too soon? Re-analysis of the April 2026 HP sweep and the
hybrid-optimizer benchmark. No optimizer is run; inputs are files already in the repository:

  benchmarks/hp_sweep_manygrains.log          per-run lines of the 58,500-run ManyGrains sweep (the
                                              per-run CSV was never committed; the log is on disk,
                                              gitignored). Misorientation is printed to 0.001 deg.
  benchmarks/hp_sweep_trajectory_manygrains.csv  step records, linked to the log by run order
                                              (checked below); gitignored, on disk.
  benchmarks/hp_sweep_threevoxels.csv          1,755 runs, 3 near-axis voxels
  benchmarks/bench_hybrid_{manygrains,threevoxels}.csv   hybrid Adam vs MC
  Examples/Example2.ManyGrains/SimInput/rand_500grains_1mm_inFZ.mic   voxel positions -> r_perp

Two success criteria appear: the recorded headline table (96% etc.) is reproduced by "final
misorientation < 1.0 deg" whatever the start offset (THR_SWEEP), and the hybrid benchmark and later
studies use "< 0.5 deg" (THR_FINE). Both are reported.

Output: benchmarks/sweep_audit/{summary.json,tables.md}.

Usage (from icenine_py/):
  uv run python scripts/sweep_audit/audit.py
"""

import json
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional

import numpy as np
import pandas as pd
from scipy import stats as sps

HERE = Path(__file__).resolve().parent
ICENINE_PY = HERE.parents[1]
BENCH = ICENINE_PY / "benchmarks"
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent / "common"))
sys.path.insert(0, str(BENCH))

import audit_helpers as A  # noqa: E402
from doc_tables import markdown_table, write_tables  # noqa: E402
from stats import mcnemar_exact, wilson  # noqa: E402

OUT_DIR = BENCH / "sweep_audit"
B1_SUMMARY = BENCH / "mc_mechanism" / "summary.json"
MIC = (
    ICENINE_PY.parent
    / "Examples"
    / "Example2.ManyGrains"
    / "SimInput"
    / "rand_500grains_1mm_inFZ.mic"
)
THR_SWEEP = 1.0  # reproduces the recorded headline table
THR_FINE = 0.5  # bench_hybrid_optimizer.py SUCCESS_THRESHOLD_DEG, later studies
PREC_DEG = [0.5, 0.1, 0.05, 0.02]
OPTS = [
    "riemannian_adam_geoopt",
    "riemannian_adam_manual",
    "riemannian_sgd_plain",
    "riemannian_sgd_momentum",
    "riemannian_sgld",
    "mc_optimizer",
]
# The headline table recorded in MIGRATION_HISTORY (success % at start offsets 1 / 2 / 5 deg) and
# the hp id the table's text names (None: not identifiable from the text).
RECORDED = {
    "riemannian_adam_geoopt": (96, 41, 0),
    "riemannian_adam_manual": (94, 49, 0),
    "riemannian_sgd_plain": (88, 52, 0),
    "riemannian_sgd_momentum": (94, 51, 0),
    "riemannian_sgld": (92, 50, 0),
    "mc_optimizer": (92, 40, 6),
}
DOCUMENTED_HP = {  # hp ids of the configs the table names (lr/steps/beta1, T, MC 3500/2/0.5)
    "riemannian_adam_geoopt": 0,
    "riemannian_adam_manual": 2,
    "riemannian_sgd_plain": 0,
    "riemannian_sgd_momentum": 0,
    "riemannian_sgld": 3,
    "mc_optimizer": 31,
}


def frac(k: int, n: int) -> Dict[str, Any]:
    lo, hi = wilson(int(k), int(n))
    return dict(k=int(k), n=int(n), p=k / max(n, 1), lo=lo, hi=hi)


def ff(f: Dict[str, Any]) -> str:
    return f"{f['k']}/{f['n']} ({100 * f['p']:.0f}%, {100 * f['lo']:.0f}-{100 * f['hi']:.0f})"


def qtxt(q: List[float], nd: int = 3) -> str:
    return "/".join(f"{x:.{nd}f}" for x in q)


def hp_table() -> Dict[Any, Dict[str, Any]]:
    import bench_hp_sweep as b  # the sweep's own grids

    out: Dict[Any, Dict[str, Any]] = {}
    for opt, cfgs in b.build_hp_grids().items():
        for c in cfgs:
            d = dict(c)
            d["budget"] = int(c["max_mc_steps"] if opt == "mc_optimizer" else c["n_steps"])
            out[(opt, int(c["hp_id"]))] = d
    return out


def hp_text(opt: str, c: Dict[str, Any]) -> str:
    if opt == "mc_optimizer":
        return (
            f"{c['max_mc_steps']} steps, {c['successive_restarts']} restarts, "
            f"step frac {c['angular_step_frac']}"
        )
    s = f"lr {c['lr']:g}, {c['n_steps']} steps"
    if "beta1" in c:
        s += f", b1 {c['beta1']}"
    if "momentum" in c:
        s += f", mom {c['momentum']}"
    if "t_init" in c:
        s += f", T {c['t_init']}"
    return s


def load_manygrains(hp: Dict[Any, Dict[str, Any]]) -> pd.DataFrame:
    rows = []
    for ln in open(BENCH / "hp_sweep_manygrains.log", errors="replace"):
        r = A.parse_log_line(ln)
        if r is not None:
            rows.append(r)
    df = pd.DataFrame(rows)
    df["run_id"] = np.arange(len(df))
    df["budget"] = [hp[(o, h)]["budget"] for o, h in zip(df.opt, df.hp)]
    from icenine.mic_file import MicFile

    mic = MicFile.read(str(MIC))
    pos = np.array([v.position for v in mic.voxels[:500]], dtype=float)
    r_all = np.hypot(pos[:, 0], pos[:, 1]) * 1e3  # um, as generate_toy_orientation_dataset
    df["r_perp"] = r_all[df.vox.values]
    return df


def load_threevoxels(hp: Dict[Any, Dict[str, Any]]) -> pd.DataFrame:
    t = pd.read_csv(BENCH / "hp_sweep_threevoxels.csv")
    df = pd.DataFrame(
        dict(
            opt=t.optimizer,
            hp=t.hp_id,
            vox=t.voxel_idx,
            pert=t.perturbation_deg.astype(int),
            mis=t.final_misorientation_deg,
            ev=t.n_evaluations,
        )
    )
    df["budget"] = [hp[(o, h)]["budget"] for o, h in zip(df.opt, df.hp)]
    return df


def config_table(df: pd.DataFrame, pert: int, thr: float) -> pd.DataFrame:
    """Per (opt, hp) success counts at the given start offset and threshold."""
    d = df[df.pert == pert]
    g = d.groupby(["opt", "hp"])
    return pd.DataFrame(
        dict(
            k=g.mis.apply(lambda m: int((m < thr).sum())),
            n=g.mis.size(),
            n_evals=g.ev.median(),
            budget=g.budget.first(),
        )
    ).reset_index()


def sel(df: pd.DataFrame, opt: str, hp: int, pert: Optional[int] = None) -> pd.DataFrame:
    d = df[(df.opt == opt) & (df.hp == hp)]
    return d if pert is None else d[d.pert == pert]


def pct_triple(mg: pd.DataFrame, o: str, hpid: int) -> tuple:
    return tuple(
        round(100 * float((sel(mg, o, hpid, p).mis < THR_SWEEP).mean())) for p in (1, 2, 5)
    )


def reference_config(mg: pd.DataFrame, hp: Dict[Any, Dict[str, Any]], o: str) -> Dict[str, Any]:
    """The config of the recorded table: among all hps whose 1/2/5 success percentages (< 1 deg)
    are within 1 point of the recorded ones, the documented hp if it is among them, else the
    smallest hp id. Returns hp id and whether the documented config reproduces the record."""
    rec = RECORDED[o]
    ok = [
        h
        for (oo, h) in hp
        if oo == o and all(abs(a - b) <= 1 for a, b in zip(pct_triple(mg, o, h), rec))
    ]
    doc = DOCUMENTED_HP[o]
    use = doc if doc in ok else (min(ok) if ok else doc)
    return dict(hp=use, documented=doc, documented_reproduces=doc in ok, matched=bool(ok))


def prec_stats(mis: np.ndarray) -> Dict[str, Any]:
    """Precision of all runs in `mis` and of those with final error < THR_SWEEP."""
    ok = mis[mis < THR_SWEEP]
    return dict(
        n=int(mis.size),
        n_ok=int(ok.size),
        q=A.quantiles(ok),
        lt={str(t): frac(int((mis < t).sum()), int(mis.size)) for t in PREC_DEG},
    )


HEADERS = {
    "hybrid_precision": {
        "ex": "example",
        "opt": "optimizer",
        "n": "runs (offsets <= 1 deg)",
        "q": "final error q25/50/75",
        "mn": "minimum",
        "lt05": "< 0.5 deg",
        "lt01": "< 0.1 deg",
        "lt002": "< 0.02 deg",
    },
    "mg_headline": {
        "opt": "optimizer",
        "cfg": "config",
        "evals": "evaluations",
        "s1": "< 1 deg, 1 deg start",
        "s2": "< 1 deg, 2 deg start",
        "s5": "< 1 deg, 5 deg start",
        "f1": "< 0.5 deg, 1 deg start",
        "f2": "< 0.5 deg, 2 deg start",
        "rec": "recorded %",
    },
    "mg_precision": {
        "opt": "optimizer",
        "set": "runs",
        "n": "n",
        "ok": "n < 1 deg",
        "q": "error of the < 1 deg runs q25/50/75/90",
        "lt05": "< 0.5 deg",
        "lt01": "< 0.1 deg",
        "lt005": "< 0.05 deg",
        "lt002": "< 0.02 deg",
    },
    "mg_by_rperp": {
        "set": "runs",
        "rbin": "r_perp",
        "sw": "< 1 deg",
        "fine": "< 0.5 deg",
        "med_ok": "median error of the < 1 deg runs",
    },
    "mg_budget": {
        "opt": "optimizer",
        "budget": "budget class",
        "cfg": "config",
        "evals": "evaluations",
        "s1": "< 1 deg, 1 deg start",
        "s2": "< 1 deg, 2 deg start",
        "f1": "< 0.5 deg, 1 deg start",
        "med_ok": "median error of the < 1 deg runs",
        "lt002": "< 0.02 deg, 1 deg start",
        "all1": "all configs in class, < 1 deg",
        "allf": "all configs in class, < 0.5 deg",
    },
    "mc_stopping": {
        "steps": "max steps",
        "restarts": "max restarts",
        "n": "n",
        "early": "ended before the budget",
        "noacc": "never accepted",
        "e_noacc": "ended early given never accepted",
        "e_acc": "ended early given accepted",
        "ev": "median evaluations",
    },
    "mc_headline_timing": {
        "pert": "start offset (deg)",
        "n": "n",
        "early": "ended before the budget",
        "nacc": "median accepts",
        "last": "step of last accept q25/50/75",
        "before10": "last accept before 10% of budget",
        "cur": "step at last accept (deg) q25/50/75",
        "mis": "median final error (deg)",
    },
    "hybrid": {
        "ex": "example",
        "pert": "start offset (deg)",
        "n": "n",
        "hyb": "hybrid < 0.5 deg",
        "mc": "MC < 0.5 deg",
        "hyb1": "hybrid < 1 deg",
        "mc1": "MC < 1 deg",
        "hmed": "hybrid median error of < 0.5 deg runs",
        "mmed": "MC median error of < 0.5 deg runs",
        "hall": "hybrid median error, all",
        "mall": "MC median error, all",
        "hev": "hybrid evaluations",
        "mev": "MC evaluations",
        "ht": "hybrid time (s)",
        "mt": "MC time (s)",
    },
}


def retitle(md: str, mapping: Dict[str, str]) -> str:
    """Replace the header cells of a markdown table (first line) by readable names."""
    head, rest = md.split("\n", 1)
    cells = [c.strip() for c in head.strip().strip("|").split("|")]
    return "| " + " | ".join(mapping.get(c, c) for c in cells) + " |\n" + rest


def main() -> None:
    hp = hp_table()
    mg = load_manygrains(hp)
    tv = load_threevoxels(hp)
    summ: Dict[str, Any] = dict(
        thr_sweep=THR_SWEEP,
        thr_fine=THR_FINE,
        n_manygrains_runs=int(len(mg)),
        n_threevoxels_runs=int(len(tv)),
    )
    tables: Dict[str, str] = {}

    vox_r = mg.drop_duplicates("vox").r_perp
    summ["mg_voxels"] = dict(
        n=int(mg.vox.nunique()),
        r_perp_q=[float(x) for x in np.quantile(vox_r, [0, 0.05, 0.5, 0.95, 1])],
        frac_beyond_120=float((vox_r > 120).mean()),
    )

    # ---- headline (reproduction from the log) -----------------------------------------------
    ref: Dict[str, Dict[str, Any]] = {o: reference_config(mg, hp, o) for o in OPTS}
    rows, head = [], {}
    for o in OPTS:
        hpid = ref[o]["hp"]
        d1, d2, d5 = (sel(mg, o, hpid, p) for p in (1, 2, 5))
        f = {
            "s1": frac(int((d1.mis < THR_SWEEP).sum()), len(d1)),
            "s2": frac(int((d2.mis < THR_SWEEP).sum()), len(d2)),
            "s5": frac(int((d5.mis < THR_SWEEP).sum()), len(d5)),
            "f1": frac(int((d1.mis < THR_FINE).sum()), len(d1)),
            "f2": frac(int((d2.mis < THR_FINE).sum()), len(d2)),
        }
        head[o] = dict(
            cfg=hp_text(o, hp[(o, hpid)]),
            recorded=list(RECORDED[o]),
            **ref[o],
            evals=float(d1.ev.median()),
            time_med=float(sel(mg, o, hpid).t.median()),
            **f,
        )
        rows.append(
            dict(
                opt=o,
                cfg=head[o]["cfg"],
                evals=head[o]["evals"],
                s1=ff(f["s1"]),
                s2=ff(f["s2"]),
                s5=ff(f["s5"]),
                f1=ff(f["f1"]),
                f2=ff(f["f2"]),
                rec="/".join(map(str, RECORDED[o])),
            )
        )
    summ["mg_headline"] = head
    tables["mg_headline"] = markdown_table(
        rows, ["opt", "cfg", "evals", "s1", "s2", "s5", "f1", "f2", "rec"], {"evals": ".0f"}
    )

    # ---- precision --------------------------------------------------------------------------
    rows, prec = [], {}
    for o in OPTS:
        best = prec_stats(sel(mg, o, head[o]["hp"], 1).mis.values)
        allc = prec_stats(mg[(mg.opt == o) & (mg.pert == 1)].mis.values)
        prec[o] = dict(ref=best, all_configs=allc)
        allst = mg[mg.opt == o].mis.values  # all start offsets (1, 2 and 5 deg)
        prec[o]["all_configs_all_starts_lt002"] = frac(int((allst < 0.02).sum()), len(allst))
        for lab, p in (("recorded config", best), ("all configs", allc)):
            rows.append(
                dict(
                    opt=o,
                    set=lab,
                    n=p["n"],
                    ok=p["n_ok"],
                    q=qtxt(p["q"]),
                    lt05=ff(p["lt"]["0.5"]),
                    lt01=ff(p["lt"]["0.1"]),
                    lt005=ff(p["lt"]["0.05"]),
                    lt002=ff(p["lt"]["0.02"]),
                )
            )
    summ["mg_precision_1deg"] = prec
    tables["mg_precision"] = markdown_table(
        rows, ["opt", "set", "n", "ok", "q", "lt05", "lt01", "lt005", "lt002"]
    )

    # ---- by distance from the rotation axis ---------------------------------------------------
    tert = np.quantile(vox_r, [1 / 3, 2 / 3])
    mg["rbin"] = np.digitize(mg.r_perp.values, tert)
    names = [f"<{tert[0]:.0f} um", f"{tert[0]:.0f}-{tert[1]:.0f} um", f">={tert[1]:.0f} um"]
    rows, rp = [], {}

    def add_rows(lab: str, d: pd.DataFrame) -> None:
        for bi, nm in enumerate(names):
            x = d[d.rbin == bi].mis.values
            f1 = frac(int((x < THR_SWEEP).sum()), len(x))
            f2 = frac(int((x < THR_FINE).sum()), len(x))
            ok = x[x < THR_SWEEP]
            rp[f"{lab}|{nm}"] = dict(
                s_sweep=f1, s_fine=f2, med_ok=float(np.median(ok)) if ok.size else float("nan")
            )
            rows.append(
                dict(set=lab, rbin=nm, sw=ff(f1), fine=ff(f2), med_ok=rp[f"{lab}|{nm}"]["med_ok"])
            )

    for lab, o in (
        ("Adam (geoopt), recorded config", "riemannian_adam_geoopt"),
        ("MC, recorded config", "mc_optimizer"),
    ):
        d = sel(mg, o, head[o]["hp"], 1)
        add_rows(lab, d)
        rho, p = sps.spearmanr(d.r_perp.values, d.mis.values)
        rp[f"{lab}|spearman"] = dict(rho=float(rho), p=float(p), n=int(len(d)))
    for lab, o in (
        ("Adam (geoopt), all configs", "riemannian_adam_geoopt"),
        ("MC, all configs", "mc_optimizer"),
    ):
        add_rows(lab, mg[(mg.opt == o) & (mg.pert == 1)])
        d = mg[(mg.opt == o) & (mg.pert == 1)]
        rho, p = sps.spearmanr(d.r_perp.values, d.mis.values)  # pooled runs share voxels
        rp[f"{lab}|spearman"] = dict(rho=float(rho), p=float(p), n=int(len(d)))
    for lab, o in (
        ("Adam (geoopt), recorded config", "riemannian_adam_geoopt"),
        ("MC, recorded config", "mc_optimizer"),
    ):
        d = sel(tv, o, head[o]["hp"], 1).sort_values("vox")
        rp[f"{lab}|ThreeVoxels"] = dict(mis=[float(x) for x in d.mis.values])
        rows.append(
            dict(
                set=lab + " (same hp id)",
                rbin="12 um (ThreeVoxels)",
                sw=f"{int((d.mis < THR_SWEEP).sum())}/{len(d)}",
                fine=f"{int((d.mis < THR_FINE).sum())}/{len(d)}",
                med_ok=float(d.mis.median()),
            )
        )
    for lab, o in (
        ("Adam (geoopt), all configs", "riemannian_adam_geoopt"),
        ("MC, all configs", "mc_optimizer"),
    ):
        d = tv[(tv.opt == o) & (tv.pert == 1)]
        f1 = frac(int((d.mis < THR_SWEEP).sum()), len(d))
        f2 = frac(int((d.mis < THR_FINE).sum()), len(d))
        ok = d.mis[d.mis < THR_SWEEP]
        rp[f"{lab}|ThreeVoxels"] = dict(s_sweep=f1, s_fine=f2)
        rows.append(
            dict(
                set=lab + " (3 voxels x configs)",
                rbin="12 um (ThreeVoxels)",
                sw=ff(f1),
                fine=ff(f2),
                med_ok=float(ok.median()) if len(ok) else float("nan"),
            )
        )
    summ["mg_by_rperp"] = dict(tertile_edges_um=[float(x) for x in tert], cells=rp)
    tables["mg_by_rperp"] = markdown_table(
        rows, ["set", "rbin", "sw", "fine", "med_ok"], {"med_ok": ".3f"}
    )

    # ---- matched budget -----------------------------------------------------------------------
    classes = [("<=200 steps", [100, 200]), ("500-1000 steps", [500, 1000]), ("3500 steps", [3500])]
    ct1 = config_table(
        mg, 1, THR_FINE
    )  # best config per budget class: most runs < 0.5 deg at 1 deg
    rows, bud = [], {}
    for o in OPTS:
        for cn, bl in classes:
            sub = ct1[(ct1.opt == o) & ct1.budget.isin(bl)]
            if sub.empty:
                continue
            b = A.pick_best(sub.to_dict("records"))  # type: ignore[arg-type]
            hpid = int(b["hp"])
            r1, r2 = sel(mg, o, hpid, 1).mis.values, sel(mg, o, hpid, 2).mis.values
            ok1 = r1[r1 < THR_SWEEP]
            allc = mg[(mg.opt == o) & mg.budget.isin(bl) & (mg.pert == 1)].mis.values
            x = dict(
                hp=hpid,
                cfg=hp_text(o, hp[(o, hpid)]),
                evals=float(b["n_evals"]),
                s1=frac(int(len(ok1)), 100),
                s2=frac(int((r2 < THR_SWEEP).sum()), 100),
                f1=frac(int((r1 < THR_FINE).sum()), 100),
                med_ok=float(np.median(ok1)) if ok1.size else float("nan"),
                lt002=frac(int((r1 < 0.02).sum()), 100),
                allcfg_s1=frac(int((allc < THR_SWEEP).sum()), len(allc)),
                allcfg_f1=frac(int((allc < THR_FINE).sum()), len(allc)),
            )
            bud[f"{o}|{cn}"] = x
            rows.append(
                dict(
                    opt=o,
                    budget=cn,
                    cfg=x["cfg"],
                    evals=x["evals"],
                    s1=ff(x["s1"]),
                    s2=ff(x["s2"]),
                    f1=ff(x["f1"]),
                    med_ok=x["med_ok"],
                    lt002=ff(x["lt002"]),
                    all1=ff(x["allcfg_s1"]),
                    allf=ff(x["allcfg_f1"]),
                )
            )
    summ["mg_budget"] = bud
    tables["mg_budget"] = markdown_table(
        rows,
        ["opt", "budget", "cfg", "evals", "s1", "s2", "f1", "med_ok", "lt002", "all1", "allf"],
        {"evals": ".0f", "med_ok": ".3f"},
    )

    # ---- MC stopping in the sweep (log evals + trajectory events) ------------------------------
    tr = pd.read_csv(BENCH / "hp_sweep_trajectory_manygrains.csv")
    g = tr[tr.event_type == "grad_step"].groupby("run_id").tail(1).set_index("run_id")
    j = mg.set_index("run_id").join(g, how="inner")
    jj = j[j.mis < 1.5]  # non-diverged runs: the last grad_step record should match the final value
    dd = np.abs(jj.mis - jj.misori_gt_deg)
    summ["trajectory_link_check"] = dict(
        n=int(len(jj)),
        median_abs_diff=float(dd.median()),
        frac_within_0p05=float((dd < 0.05).mean()),
        frac_within_0p05_all_gradient_runs=float((np.abs(j.mis - j.misori_gt_deg) < 0.05).mean()),
    )
    assert summ["trajectory_link_check"]["frac_within_0p05"] > 0.9
    acc = tr[tr.event_type == "mc_accept"]
    la = acc.groupby("run_id").agg(n_acc=("step", "size"), last_step=("step", "max"))
    last_cur = (
        acc.sort_values(["run_id", "step"])
        .groupby("run_id")
        .tail(1)
        .set_index("run_id")
        .cur_step_rad
    )
    mc = mg[mg.opt == "mc_optimizer"].set_index("run_id").join(la)
    mc["n_acc"] = mc.n_acc.fillna(0)
    mc["last_cur_deg"] = np.degrees(last_cur.reindex(mc.index))
    mc["early"] = mc.ev < mc.budget + 1
    rows, mcs = [], {}
    for steps in (100, 500, 1000, 3500):
        for rs in (0, 2, 5):
            ids = [
                h
                for (o, h), c in hp.items()
                if o == "mc_optimizer" and c["budget"] == steps and c["successive_restarts"] == rs
            ]
            d = mc[mc.hp.isin(ids)]
            noacc = d.n_acc == 0
            x = dict(
                n=int(len(d)),
                early=frac(int(d.early.sum()), len(d)),
                noacc=frac(int(noacc.sum()), len(d)),
                early_given_noacc=frac(int((d.early & noacc).sum()), int(noacc.sum())),
                early_given_acc=frac(int((d.early & ~noacc).sum()), int((~noacc).sum())),
                med_evals=float(d.ev.median()),
            )
            mcs[f"{steps}|{rs}"] = x
            rows.append(
                dict(
                    steps=steps,
                    restarts=rs,
                    n=x["n"],
                    early=ff(x["early"]),
                    noacc=ff(x["noacc"]),
                    e_noacc=ff(x["early_given_noacc"]),
                    e_acc=ff(x["early_given_acc"]),
                    ev=x["med_evals"],
                )
            )
    summ["mc_stopping"] = mcs
    tables["mc_stopping"] = markdown_table(
        rows, ["steps", "restarts", "n", "early", "noacc", "e_noacc", "e_acc", "ev"], {"ev": ".0f"}
    )
    hpid = head["mc_optimizer"]["hp"]
    rows, hl = [], {}
    for p in (1, 2, 5):
        d = mc[(mc.hp == hpid) & (mc.pert == p)]
        a = d[d.n_acc > 0]
        x = dict(
            n=int(len(d)),
            early=frac(int(d.early.sum()), len(d)),
            n_acc_med=float(a.n_acc.median()),
            last_step_q=A.quantiles(a.last_step, (0.25, 0.5, 0.75)),
            last_before_10pct=frac(int((a.last_step < 0.1 * d.budget.iloc[0]).sum()), len(a)),
            last_cur_deg_q=A.quantiles(a.last_cur_deg, (0.25, 0.5, 0.75)),
            mis_med=float(d.mis.median()),
        )
        hl[str(p)] = x
        rows.append(
            dict(
                pert=p,
                n=x["n"],
                early=ff(x["early"]),
                nacc=x["n_acc_med"],
                last=qtxt(x["last_step_q"], 0),
                before10=ff(x["last_before_10pct"]),
                cur=qtxt(x["last_cur_deg_q"], 4),
                mis=x["mis_med"],
            )
        )
    summ["mc_headline_timing"] = hl
    tables["mc_headline_timing"] = markdown_table(
        rows,
        ["pert", "n", "early", "nacc", "last", "before10", "cur", "mis"],
        {"nacc": ".0f", "mis": ".3f"},
    )

    # ---- hybrid benchmark ---------------------------------------------------------------------
    hyb: Dict[str, Any] = {}
    prec_rows: List[Dict[str, Any]] = []
    rows = []
    for ex in ("manygrains", "threevoxels"):
        h = pd.read_csv(BENCH / f"bench_hybrid_{ex}.csv")
        pair_keys = list(zip(h.voxel_idx, h.perturbation_deg))[::2]  # rows come as (hybrid, mc)
        blocks = A.split_complete_blocks(pair_keys, 0)
        full = max(blocks, key=lambda b: b[1] - b[0])
        used = h.iloc[2 * full[0] : 2 * full[1]]
        hyb[ex] = dict(
            rows_in_file=int(len(h)),
            blocks=[int(b[1] - b[0]) for b in blocks],
            rows_used=int(len(used)),
        )
        for p in sorted(used.perturbation_deg.unique()):
            sub = used[used.perturbation_deg == p]
            hy, m = sub[sub.optimizer == "hybrid_adam"], sub[sub.optimizer == "mc_optimizer"]
            n = len(hy)
            okh = hy.final_misori_deg[hy.final_misori_deg < THR_FINE]
            okm = m.final_misori_deg[m.final_misori_deg < THR_FINE]
            cell = dict(
                n=n,
                hybrid=frac(len(okh), n),
                mc=frac(len(okm), n),
                hybrid_1=frac(int((hy.final_misori_deg < THR_SWEEP).sum()), n),
                mc_1=frac(int((m.final_misori_deg < THR_SWEEP).sum()), n),
                hybrid_med_ok=float(okh.median()) if len(okh) else float("nan"),
                mc_med_ok=float(okm.median()) if len(okm) else float("nan"),
                hybrid_med_all=float(hy.final_misori_deg.median()),
                mc_med_all=float(m.final_misori_deg.median()),
                hybrid_evals=[float(hy.n_hard_evals.median()), float(hy.n_diff_evals.median())],
                mc_evals=float(m.n_hard_evals.median()),
                hybrid_time=float(hy.wall_time_sec.median()),
                mc_time=float(m.wall_time_sec.median()),
            )
            hyb[ex][f"pert{p:g}"] = cell
            rows.append(
                dict(
                    ex=ex,
                    pert=p,
                    n=n,
                    hyb=ff(cell["hybrid"]),
                    mc=ff(cell["mc"]),
                    hyb1=ff(cell["hybrid_1"]),
                    mc1=ff(cell["mc_1"]),
                    hmed=cell["hybrid_med_ok"],
                    mmed=cell["mc_med_ok"],
                    hall=cell["hybrid_med_all"],
                    mall=cell["mc_med_all"],
                    hev=f"{cell['hybrid_evals'][0]:.0f}+{cell['hybrid_evals'][1]:.0f}d",
                    mev=cell["mc_evals"],
                    ht=cell["hybrid_time"],
                    mt=cell["mc_time"],
                )
            )
        low = used[used.perturbation_deg <= 1.0]
        for lab, o in (("hybrid", "hybrid_adam"), ("mc", "mc_optimizer")):
            e = low[low.optimizer == o].final_misori_deg.values
            hyb[ex][f"{lab}_low_offsets"] = dict(
                n=int(len(e)),
                min=float(e.min()),
                q=A.quantiles(e, (0.25, 0.5, 0.75)),
                **{f"lt{t}": frac(int((e < t).sum()), len(e)) for t in (0.5, 0.1, 0.05, 0.02)},
            )
            prec_rows.append(
                dict(
                    ex=ex,
                    opt=lab,
                    n=len(e),
                    q=qtxt(A.quantiles(e, (0.25, 0.5, 0.75))),
                    mn=float(e.min()),
                    lt05=ff(hyb[ex][f"{lab}_low_offsets"]["lt0.5"]),
                    lt01=ff(hyb[ex][f"{lab}_low_offsets"]["lt0.1"]),
                    lt002=ff(hyb[ex][f"{lab}_low_offsets"]["lt0.02"]),
                )
            )
        hy = used[used.optimizer == "hybrid_adam"].set_index(["voxel_idx", "perturbation_deg"])
        mm = used[used.optimizer == "mc_optimizer"].set_index(["voxel_idx", "perturbation_deg"])
        mm = mm.reindex(hy.index)
        sh, sm = (hy.final_misori_deg < THR_FINE).values, (mm.final_misori_deg < THR_FINE).values
        b_, c_ = int((sh & ~sm).sum()), int((~sh & sm).sum())
        vdiff = pd.Series(sh.astype(int) - sm.astype(int), index=hy.index).groupby(level=0).sum()
        pos, neg = int((vdiff > 0).sum()), int((vdiff < 0).sum())
        hyb[ex]["paired"] = dict(
            pairs=int(len(hy)),
            hybrid_only=b_,
            mc_only=c_,
            mcnemar_p=mcnemar_exact(b_, c_),
            voxels=int(len(vdiff)),
            voxels_hybrid_ahead=pos,
            voxels_mc_ahead=neg,
            voxel_sign_p=mcnemar_exact(pos, neg),
        )
    summ["hybrid"] = hyb
    tables["hybrid"] = markdown_table(
        rows,
        [
            "ex",
            "pert",
            "n",
            "hyb",
            "mc",
            "hyb1",
            "mc1",
            "hmed",
            "mmed",
            "hall",
            "mall",
            "hev",
            "mev",
            "ht",
            "mt",
        ],
        {
            "hmed": ".3f",
            "mmed": ".3f",
            "hall": ".3f",
            "mall": ".3f",
            "mev": ".0f",
            "ht": ".2f",
            "mt": ".2f",
        },
    )

    tables["hybrid_precision"] = markdown_table(
        prec_rows, ["ex", "opt", "n", "q", "mn", "lt05", "lt01", "lt002"], {"mn": ".3f"}
    )
    # ---- measured vs deployed ------------------------------------------------------------------
    b1_path = B1_SUMMARY
    if b1_path.exists():
        b1 = json.loads(b1_path.read_text())
        setup = b1["setup"]
        fin = b1["curves"]["H3|realistic"]["points"]["result"]["d0_med"]
        st = b1["mechanism"]["H3|realistic"].get("start_err_q", [float("nan")] * 3)
        evq = b1["finisher_evals_q"]["H3|realistic"]
        m_ref = head["mc_optimizer"]
        pa, pm = prec["riemannian_adam_geoopt"]["ref"], prec["mc_optimizer"]["ref"]
        hl = summ["mc_headline_timing"]
        bm = bud["mc_optimizer|<=200 steps"]
        b3 = bud["mc_optimizer|3500 steps"]
        rp_cell = summ["mg_by_rperp"]["cells"]["Adam (geoopt), recorded config|spearman"]
        vq = summ["mg_voxels"]["r_perp_q"]
        l1 = lambda f: f"{f['k']}/{f['n']}"  # noqa: E731
        gap = [
            dict(
                aspect="success criterion",
                measured=(
                    f"final error < {THR_SWEEP:g} deg from a 1/2/5 deg start; recorded 96/41/0% "
                    "(Adam), 92/40/6% (MC). At a 1 deg start this means 'ended below the start'"
                ),
                deployed="a finisher meant to end within about 0.01-0.03 deg of the truth",
                gap="the criterion does not measure the deployed target",
            ),
            dict(
                aspect="precision of the 'successes'",
                measured=(
                    f"median final error {pa['q'][1]:.2f} deg (Adam) and "
                    f"{pm['q'][1]:.2f} deg (MC); "
                    f"under 0.1 deg: {l1(pa['lt']['0.1'])} and {l1(pm['lt']['0.1'])}; "
                    f"under 0.02 deg: {l1(pa['lt']['0.02'])} and {l1(pm['lt']['0.02'])}"
                ),
                deployed=f"finisher result median {fin:.4f} deg from the truth (T5, H3 realistic)",
                gap=(
                    f"the sweep's successes are about {pm['q'][1] / fin:.0f}x (MC) coarser than "
                    "the "
                    "finisher's result; the 0.01-0.03 deg scale was not measured"
                ),
            ),
            dict(
                aspect="MC settings",
                measured=(
                    f"{m_ref['cfg']}; box = 1.5 x the start offset, step = 0.5 x box "
                    "(0.75 deg at a 1 deg start)"
                ),
                deployed=(
                    f"{int(setup['max_mc_steps'])} steps, {int(setup['max_restarts'])} restarts, "
                    f"box {setup['box_deg']:.3f} deg, step {setup['step0_deg']:.4f} deg"
                ),
                gap="different step scale, box and budget",
            ),
            dict(
                aspect="starting point",
                measured="the truth rotated by exactly 1, 2 or 5 deg about a random axis",
                deployed=(
                    "a coarse-search or network start; T5 H3 start error "
                    f"{st[1]:.3f} deg (q25-q75 {st[0]:.3f}-{st[2]:.3f})"
                ),
                gap=(
                    f"the sweep's starts are {1.0 / st[1]:.0f}x (1 deg) to {5.0 / st[1]:.0f}x "
                    "(5 deg) the T5 H3 median start error"
                ),
            ),
            dict(
                aspect="budget",
                measured=(
                    "Adam/SGD 101-501 evaluations (gradient evaluations); MC 101-3501; "
                    "hybrid 4 hard + 303 differentiable"
                ),
                deployed=(
                    f"MC {int(setup['max_mc_steps'])} steps, then VarianceMinimizing; whole "
                    f"finisher median {evq[1]:.0f} evaluations (T5 H3)"
                ),
                gap="MC at 101 and at 3501 evaluations give the same result (next row)",
            ),
            dict(
                aspect="matched budget (MC)",
                measured=(
                    f"100 steps (best config): {ff(bm['f1'])} under 0.5 deg at a 1 deg start, "
                    f"median error of the < 1 deg runs {bm['med_ok']:.3f}; 3500 steps (best "
                    f"config, 5 restarts): {ff(b3['f1'])}, {b3['med_ok']:.3f}; pooled over each "
                    f"class {ff(bm['allcfg_f1'])} vs {ff(b3['allcfg_f1'])}; the recorded 3500-step "
                    f"2-restart config: {ff(m_ref['f1'])}, {pm['q'][1]:.3f}"
                ),
                deployed="200 steps (not run in the sweep)",
                gap=(
                    "more MC steps bought no accuracy in the sweep between 100 and 3500 steps; "
                    "200 steps were not run"
                ),
            ),
            dict(
                aspect="data and voxels",
                measured=(
                    f"ManyGrains: 100 voxels, r_perp {vq[0]:.0f}-{vq[4]:.0f} um (median "
                    f"{vq[2]:.0f}); ThreeVoxels: 3 voxels at 12 um"
                ),
                deployed="per-voxel finisher on the reconstruction's voxels",
                gap=(
                    "no detectable r_perp dependence inside the sweep's range "
                    "(Spearman of error vs r_perp, "
                    f"recorded Adam config: rho {rp_cell['rho']:.2f}, p {rp_cell['p']:.2f}, "
                    "n 100 voxels)"
                ),
            ),
            dict(
                aspect="MC stopping in the sweep",
                measured=(
                    "recorded MC config at a 1 deg start: last accepted move before 10% of the "
                    f"budget in {ff(hl['1']['last_before_10pct'])}"
                ),
                deployed="same rule, 200 steps",
                gap="consistent with the step collapse measured in B1; never examined in the sweep",
            ),
        ]
        summ["gap"] = gap
        tables["gap"] = markdown_table(gap, ["aspect", "measured", "deployed", "gap"])
    tables = {k: retitle(v, HEADERS.get(k, {})) for k, v in tables.items()}
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    (OUT_DIR / "summary.json").write_text(json.dumps(summ, indent=1, default=float))
    write_tables(OUT_DIR / "tables.md", tables)
    print("wrote", OUT_DIR)
    for o in OPTS:
        print(
            o,
            head[o]["recorded"],
            pct_triple(mg, o, head[o]["hp"]),
            head[o]["documented_reproduces"],
        )


if __name__ == "__main__":
    main()
