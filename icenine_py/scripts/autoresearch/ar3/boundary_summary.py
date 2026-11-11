"""AR3 step 5: summary of part 0, part 1 (the pass on the stored maps) and part 2 (runs K, KE).

  uv run python scripts/autoresearch/ar3/boundary_summary.py
  uv run python scripts/autoresearch/ar3/boundary_ratio.py      # optional, after this script

Run order: this script first, then ``boundary_ratio.py``. This script rewrites ``tables.md``
completely but keeps an ``ar3e_ratio`` block that is already there (``boundary_ratio.py`` writes
that block), so running it alone does not drop the ratio table; on a fresh ``tables.md`` run
``boundary_ratio.py`` after it.

Reads ``boundary/part0.json``, ``boundary/runs/<map>/tries.npz`` and ``pass.json`` (and
``<map>_second`` when present), the run.json files of the compared runs, ``boundary/runs/K/
progress.json`` (run K stopped at its cap and has no run.json), and the stored ``recon.mic`` of
D4, E4, KE, D and D2 (the old-end score check). Writes ``boundary/summary.json``,
``boundary/tables.md``, ``boundary/ar3_boundary.png``, ``boundary/ar3_boundary_pipeline.png``
and ``boundary/cache/after_<map>.npy`` (the map after the picked rule, input of the second pass).
Truth is used to score and to pick the rule on set A (even true-grain id), reported on B (odd).
Wilson 95% intervals, no p-values; a gain or a loss only when intervals do not overlap, and only
on one footing (the pass on both sides). One run, one voxel order. Every wall time here was
measured under contention and is not a timing result.
"""

import gzip
import json
from typing import Any, Dict, List

import boundary_common as B
import boundary_lib as BL
import numpy as np
from fine_summary import overlap, wil, wil_str

import doc_tables as DT  # noqa: E402
import stats as ST  # noqa: E402

RUN_EVALS = {  # evaluations of the run that made the map (refine stage included)
    "D4": ("fine/runs/D4", None),
    "P4R4": ("fine/runs/P4", "evals_r4"),
    "E4": ("fine/runs/E4", None),
    "D": ("gradient/runs/D", None),
    "PR": ("gradient/runs/P", "evals_R"),
    "D2": ("layer/runs/D", None),
    "K": ("boundary/runs/K", None),
    "KE": ("boundary/runs/KE", None),
}
PAIRS = [
    ("D4", "D"),
    ("P4R4", "PR"),
    ("E4", "D2"),
    ("K", "D4"),
    ("K", "D"),
    ("KE", "E4"),
    ("KE", "D2"),
]
RH_GRID = [0.90, 0.92, 0.94, 0.96, 0.98, 1.00]
CQ_GRID = [25, 50, 75, 90, 100]  # cost quantiles (pooled step 4 maps); 100 = cost not used
# stored maps whose recon.mic carries the run's own scores (columns 9, 10, 11: confidence, cost,
# hit ratio); PR is stored as float32 orientations only, P4R4 is pieced together from two files
RECON_MIC = {
    "D4": "fine/runs/D4/recon.mic.gz",
    "E4": "fine/runs/E4/recon.mic.gz",
    "KE": "boundary/runs/KE/recon.mic.gz",
    "D": "gradient/runs/D/recon.mic.gz",
    "D2": "layer/runs/D/recon.mic",
}
OLD_END_TOL = 1e-3
# the pipelines of the lead plot: (map, sample panel)
PLOT_PIPES = [
    ("D4", "grad"),
    ("D", "grad"),
    ("P4R4", "grad"),
    ("PR", "grad"),
    ("E4", "layer"),
    ("KE", "layer"),
    ("D2", "layer"),
]
# label placement in the lead plot: (x offset, y offset in points, horizontal alignment)
LABEL_POS = {
    "D4": (-10, 0, "right"),
    "P4R4": (10, 10, "left"),
    "PR": (-10, -26, "right"),
    "D": (10, 0, "left"),
    "KE": (0, 34, "center"),
    "E4": (10, 0, "left"),
    "D2": (10, 0, "left"),
}
K_STAGES = ("n_evals_seed", "n_evals_neighbor", "n_evals_revisit", "n_evals_reseed")


def run_evals(m: str) -> int:
    d, extra = RUN_EVALS[m]
    j = json.loads((B.C.OUT / d / "run.json").read_text())
    e = int(j["evals_total"])
    return e + (int(j[extra]) if extra else 0)


def four(err: np.ndarray) -> Dict[str, Any]:
    n = len(err)
    return dict(
        n=n,
        median=float(np.median(err)),
        p95=float(np.percentile(err, 95)),
        gt01=wil(int((err > 0.1).sum()), n),
        wrong=wil(int((err > B.WRONG).sum()), n),
    )


def main() -> None:
    S: Dict[str, Any] = {}
    T: Dict[str, str] = {}
    p0 = json.loads((B.BND / "part0.json").read_text())
    have = [m for m in B.ALL_MAPS if (B.RUNS / m / "tries.npz").exists()]
    envs: Dict[str, B.Env] = {}
    for m in have:
        k = B.ALL_MAPS[m][0]
        if k not in envs:
            envs[k] = B.Env(k)
    R0: Dict[str, np.ndarray] = {}
    E0: Dict[str, np.ndarray] = {}
    TR: Dict[str, BL.Tries] = {}
    for m in have:
        env = envs[B.ALL_MAPS[m][0]]
        R0[m] = B.load_map(m, env.n)
        E0[m] = env.err(R0[m])
        TR[m] = BL.Tries(B.RUNS / m)
    # ------------------------------------------------------------ part 0 table
    rows = []
    for m, st in p0["maps"].items():
        rs = p0["rows"][m]
        nw = [r for r in rs if r["newly"]]
        both = [
            r
            for r in nw
            if r["cand_within_1"] and (r["better_cost"] or r["better_hit"] or r["better_conf"])
        ]
        rows.append(
            dict(
                map=m,
                boundary=st["n_boundary"],
                cands=st["n_candidates"],
                wrong=st["n_wrong"],
                wrong_cand1=sum(r["cand_within_1"] for r in rs),
                tr_cost=sum(r["better_cost"] for r in rs),
                tr_hit=sum(r["better_hit"] for r in rs),
                tr_conf=sum(r["better_conf"] for r in rs),
                newly=len(nw),
                newly_both=len(both),
            )
        )
    S["part0"] = dict(rows=rows, gate=p0["gate"])
    # new keys (review fix 7): shares of boundary voxels, candidates per boundary voxel, totals
    shares: Dict[str, Dict[str, Any]] = {}
    p0_table = []
    for r in rows:
        n_m = envs[B.ALL_MAPS[r["map"]][0]].n
        shares[r["map"]] = dict(
            n=n_m,
            boundary_share=r["boundary"] / n_m,
            cands_per_boundary=r["cands"] / r["boundary"],
        )
        p0_table.append(
            dict(
                r,
                bnd_share=f"{100 * r['boundary'] / n_m:.1f}%",
                cands_per_bnd=r["cands"] / r["boundary"],
            )
        )
    S["part0_totals"] = dict(
        wrong=sum(r["wrong"] for r in rows),
        tr_hit=sum(r["tr_hit"] for r in rows),
        tr_cost=sum(r["tr_cost"] for r in rows),
        per_map=shares,
    )
    T["ar3e_part0"] = DT.markdown_table(
        p0_table,
        [
            "map",
            "boundary",
            "bnd_share",
            "cands",
            "cands_per_bnd",
            "wrong",
            "wrong_cand1",
            "tr_cost",
            "tr_hit",
            "tr_conf",
            "newly",
            "newly_both",
        ],
        {"boundary": ",", "cands": ",", "cands_per_bnd": ".2f"},
    )
    # ------------------------------------------------------------ rules on the step 4 maps
    s4 = [m for m in B.STEP4 if m in TR]
    after: Dict[str, Dict[str, np.ndarray]] = {m: {} for m in have}
    took: Dict[str, Dict[str, np.ndarray]] = {m: {} for m in have}
    errs: Dict[str, Dict[str, np.ndarray]] = {m: {} for m in have}
    for m in have:
        env = envs[B.ALL_MAPS[m][0]]
        for r in BL.RULES:
            after[m][r], took[m][r] = BL.apply_rule(TR[m], R0[m], r)
            errs[m][r] = env.err(after[m][r])
    selA = {m: envs[B.ALL_MAPS[m][0]].isA for m in have}
    rule_rows: List[Dict[str, Any]] = []
    pickA: Dict[str, Any] = {}
    for part, name in ((True, "A"), (False, "B")):
        for r in BL.RULES:
            k = n = ch = fx = bk = 0
            for m in s4:
                sel = selA[m] if part else ~selA[m]
                c = BL.counts_after(E0[m], errs[m][r], sel)
                n += c["n"]
                k += c["wrong_after"]
                fx += c["fixed"]
                bk += c["broken"]
                ch += int((took[m][r] & sel).sum())
            w = wil(k, n)
            pickA[f"{name}/{r}"] = dict(wrong=w, changed=ch, fixed=fx, broken=bk)
            rule_rows.append(
                dict(set=name, rule=r, wrong=wil_str(w), changed=ch, fixed=fx, broken=bk)
            )
    fewest = min(BL.RULES, key=lambda r: pickA[f"A/{r}"]["wrong"]["k"])
    close = [
        r for r in BL.RULES if overlap(pickA[f"A/{r}"]["wrong"], pickA[f"A/{fewest}"]["wrong"])
    ]
    pick = min(close, key=lambda r: pickA[f"A/{r}"]["changed"])
    S["pick"] = dict(fewest_wrong_on_A=fewest, overlapping=close, picked=pick, table=pickA)
    T["ar3e_rules"] = DT.markdown_table(
        rule_rows, ["set", "rule", "wrong", "changed", "fixed", "broken"]
    )
    # new keys (review fix 6): how close the pick is, and hit against cost paired by voxel
    margin: Dict[str, Any] = {}
    for other in BL.RULES:
        if other == pick:
            continue
        wp, wo = pickA[f"A/{pick}"]["wrong"], pickA[f"A/{other}"]["wrong"]
        margin[other] = dict(
            pick_hi_A=wp["hi"],
            other_lo_A=wo["lo"],
            gap_A=wo["lo"] - wp["hi"],
            gap_A_voxels=(wo["lo"] - wp["hi"]) * wp["n"],
            overlap_A=overlap(wp, wo),
            overlap_B=overlap(pickA[f"B/{pick}"]["wrong"], pickA[f"B/{other}"]["wrong"]),
            changed_A_pick=pickA[f"A/{pick}"]["changed"],
            changed_A_other=pickA[f"A/{other}"]["changed"],
            tie_break_if_overlap=min(
                (pick, other), key=lambda r: pickA[f"A/{r}"]["changed"]  # type: ignore
            ),
        )
    S["pick_margin"] = margin
    pp_rows = []
    S["pick_paired"] = {}
    for part, name in ((True, "A"), (False, "B")):
        for other in BL.RULES:
            if other == pick:
                continue
            wo_p = wp_o = both = n = 0
            for m in s4:
                sel = selA[m] if part else ~selA[m]
                w_pick = errs[m][pick][sel] > B.WRONG
                w_oth = errs[m][other][sel] > B.WRONG
                a_only, b_only = ST.paired_discordant(w_pick, w_oth)
                wp_o += a_only
                wo_p += b_only
                both += int((w_pick & w_oth).sum())
                n += int(sel.sum())
            d = dict(
                set=name,
                pick=pick,
                other=other,
                n=n,
                pick_wrong_only=wp_o,
                other_wrong_only=wo_p,
                both_wrong=both,
            )
            S["pick_paired"][f"{name}/{other}"] = d
            pp_rows.append(d)
    T["ar3e_pick_paired"] = DT.markdown_table(
        pp_rows,
        ["set", "pick", "other", "n", "pick_wrong_only", "other_wrong_only", "both_wrong"],
        {"n": ","},
    )
    # ------------------------------------------------------------ per map, picked rule
    pm_rows, ot_rows, er_rows, lv_rows, pass_rows = [], [], [], [], []
    S["maps"] = {}
    for m in have:
        env = envs[B.ALL_MAPS[m][0]]
        t = TR[m]
        ra = after[m][pick]
        ea = errs[m][pick]
        np.save(B.CACHE / f"after_{m}.npy", ra)
        c = BL.counts_after(E0[m], ea, np.ones(env.n, bool))
        moved = env.mis(ra, R0[m]) > 0.1
        d = dict(counts=c, changed=int(took[m][pick].sum()), moved_gt01=int(moved.sum()))
        ev_run = run_evals(m)
        ev_pass = int(t.info["evals"])
        d.update(
            boundary=int(t.info["n_boundary_all"]),
            tries=int(t.info["n_tries"]),
            evals_pass=ev_pass,
            evals_run=ev_run,
            share=ev_pass / ev_run,
            wall_tries_s=t.info["wall_tries_s"],
            timing_label=t.info["timing_label"],
        )
        S["maps"][m] = d
        pass_rows.append(
            dict(
                map=m,
                boundary=d["boundary"],
                tries=d["tries"],
                evals_pass=ev_pass,
                evals_run=ev_run,
                share=f"{100 * ev_pass / ev_run:.1f}%",
            )
        )
        pm_rows.append(
            dict(
                map=m,
                wrong_before=c["wrong_before"],
                wrong_after=c["wrong_after"],
                changed=d["changed"],
                moved_gt01=d["moved_gt01"],
                fixed=c["fixed"],
                broken=c["broken"],
                both_wrong=c["both_wrong"],
            )
        )
        for r in BL.RULES:
            if r == pick:
                continue
            cc = BL.counts_after(E0[m], errs[m][r], np.ones(env.n, bool))
            ot_rows.append(
                dict(
                    map=m,
                    rule=r,
                    wrong_after=cc["wrong_after"],
                    changed=int(took[m][r].sum()),
                    fixed=cc["fixed"],
                    broken=cc["broken"],
                )
            )
        for g, sel in BL.groups(env):
            for lab, e in (("before", E0[m]), ("after", ea)):
                f = four(e[sel])
                S["maps"][m].setdefault("err", {})[f"{g}/{lab}"] = f
                if g == "all" or env.level is not None:
                    row = dict(
                        map=m,
                        group=g,
                        when=lab,
                        n=f["n"],
                        wrong=wil_str(f["wrong"]),
                        gt01=wil_str(f["gt01"]),
                        median=f["median"],
                        p95=f["p95"],
                    )
                    (er_rows if g == "all" else lv_rows).append(row)
    cols = ["map", "group", "when", "n", "wrong", "gt01", "median", "p95"]
    T["ar3e_pass"] = DT.markdown_table(
        pass_rows,
        ["map", "boundary", "tries", "evals_pass", "evals_run", "share"],
        {k: "," for k in ("boundary", "tries", "evals_pass", "evals_run")},
    )
    T["ar3e_paired"] = DT.markdown_table(pm_rows, list(pm_rows[0]))
    T["ar3e_other_rules"] = DT.markdown_table(
        ot_rows, ["map", "rule", "wrong_after", "changed", "fixed", "broken"]
    )
    T["ar3e_error"] = DT.markdown_table(er_rows, cols, {"n": ",", "median": ".4f", "p95": ".3f"})
    T["ar3e_level"] = DT.markdown_table(lv_rows, cols, {"n": ",", "median": ".4f", "p95": ".3f"})
    # ------------------------------------------------------------ whole pipeline
    pipe = []
    for a, b in PAIRS:
        if a not in have or b not in have:
            continue
        env = envs[B.ALL_MAPS[a][0]]
        ea_, eb_ = errs[a][pick], errs[b][pick]
        for base_lab, eb, evb in (
            ("run only", E0[b], run_evals(b)),
            ("run plus pass", eb_, run_evals(b) + int(TR[b].info["evals"])),
        ):
            c = BL.counts_after(eb, ea_, np.ones(env.n, bool))
            wa, wb = wil(c["wrong_before"], env.n), wil(c["wrong_after"], env.n)
            verdict = (
                "no gain or loss (intervals overlap)"
                if overlap(wa, wb)
                else ("gain" if c["wrong_after"] < c["wrong_before"] else "loss")
            )
            pipe.append(
                dict(
                    this=a + " plus pass",
                    against=f"{b} ({base_lab})",
                    evals_this=run_evals(a) + int(TR[a].info["evals"]),
                    evals_against=evb,
                    wrong_this=wil_str(wb),
                    wrong_against=wil_str(wa),
                    fixed=c["fixed"],
                    broken=c["broken"],
                    both_wrong=c["both_wrong"],
                    verdict=verdict,
                )
            )
    S["pipeline"] = pipe
    # review fix 3: one footing (the pass on both sides) in its own table with a verdict; the
    # one-sided rows (pass on one side only) in a second table with no verdict
    pipe_cols = [
        "this",
        "against",
        "evals_this",
        "evals_against",
        "wrong_this",
        "wrong_against",
        "fixed",
        "broken",
        "both_wrong",
    ]
    fmt_ev = {"evals_this": ",", "evals_against": ","}
    T["ar3e_pipeline"] = DT.markdown_table(
        [p for p in pipe if p["against"].endswith("(run plus pass)")],
        pipe_cols + ["verdict"],
        fmt_ev,
    )
    T["ar3e_pipeline_onesided"] = DT.markdown_table(
        [p for p in pipe if p["against"].endswith("(run only)")], pipe_cols, fmt_ev
    )
    # new key (addition A): the points of the lead plot, run plus pass on every side
    pts = []
    for m, panel in PLOT_PIPES:
        if m not in have:
            continue
        c = S["maps"][m]["counts"]
        w = wil(c["wrong_after"], c["n"])
        pts.append(
            dict(
                map=m,
                panel=panel,
                evals=run_evals(m) + int(TR[m].info["evals"]),
                wrong=c["wrong_after"],
                n=c["n"],
                lo=w["lo"],
                hi=w["hi"],
            )
        )
    S["pipeline_plot"] = pts
    # ------------------------------------------------------------ evaluations by stage (BFS runs)
    kr = []
    for m in ("K", "D4", "D", "KE", "E4", "D2"):
        rp = B.C.OUT / RUN_EVALS[m][0] / "run.json"
        if not rp.exists():
            continue
        st = json.loads(rp.read_text())["stats"]
        g = lambda k: int(st.get(k, 0))  # noqa: E731
        kr.append(
            dict(
                run=m,
                seed=g("n_evals_seed"),
                neighbor=g("n_evals_neighbor"),
                revisit=g("n_evals_revisit"),
                reseed=g("n_evals_reseed"),
                total=int(json.loads(rp.read_text())["evals_total"]),
                entries=g("n_seeds") + g("n_reseeds"),
                index_used=g("index_used"),
                grid_fallbacks=g("index_fallback_seed_test"),
                grid_evals=g("index_fallback_grid_evals"),
                first_accepted=g("n_neighbor_accepted_first"),
                neighbor_fits=g("n_neighbor_fits"),
                retries=g("n_retry_attempted_neighbor") + g("n_retry_attempted_revisit"),
                retries_accepted=g("n_retry_accepted_neighbor") + g("n_retry_accepted_revisit"),
                wall_s=float(json.loads(rp.read_text())["wall_total_s"]),
            )
        )
    S["stages"] = kr
    # new keys (review fix 5): run K stopped at its cap, so its numbers come from its last
    # progress record; the derived sums quoted on the page are computed here
    kp = B.RUNS / "K" / "progress.json"
    k_row = None
    if kp.exists():
        pj = json.loads(kp.read_text())
        st = pj["stats"]
        g = lambda k: int(st.get(k, 0))  # noqa: E731
        kd = dict(
            source="boundary/runs/K/progress.json",
            status="unfinished, stopped at its cap",
            elapsed_s=float(pj["elapsed_s"]),
            seeds=g("n_seeds"),
            reseeds=g("n_reseeds"),
            entries=g("n_seeds") + g("n_reseeds"),
            index_used=g("index_used"),
            grid_fallbacks=g("index_fallback_seed_test"),
            grid_evals=g("index_fallback_grid_evals"),
            evals_seed=g("n_evals_seed"),
            evals_neighbor=g("n_evals_neighbor"),
            evals_revisit=g("n_evals_revisit"),
            evals_reseed=g("n_evals_reseed"),
            evals_total=sum(g(k) for k in K_STAGES),
            evals_neighbor_revisit=g("n_evals_neighbor") + g("n_evals_revisit"),
            neighbor_fits=g("n_neighbor_fits"),
            first_accepted=g("n_neighbor_accepted_first"),
            retries=g("n_retry_attempted_neighbor") + g("n_retry_attempted_revisit"),
            retries_accepted=g("n_retry_accepted_neighbor") + g("n_retry_accepted_revisit"),
        )
        S["K_progress"] = kd
        k_row = dict(
            run=f"K (unfinished, at {kd['elapsed_s']:,.0f} s)",
            seed=kd["evals_seed"],
            neighbor=kd["evals_neighbor"],
            revisit=kd["evals_revisit"],
            reseed=kd["evals_reseed"],
            total=kd["evals_total"],
            entries=kd["entries"],
            index_used=kd["index_used"],
            grid_fallbacks=kd["grid_fallbacks"],
            grid_evals=kd["grid_evals"],
            neighbor_fits=kd["neighbor_fits"],
            first_accepted=kd["first_accepted"],
            retries=kd["retries"],
            retries_accepted=kd["retries_accepted"],
            wall_s=kd["elapsed_s"],
        )
    by_run = {r["run"]: r for r in kr}
    derived: Dict[str, Any] = {}
    for m in ("KE", "E4", "D4", "D2"):
        if m in by_run:
            derived[f"{m}_evals_neighbor_revisit"] = by_run[m]["neighbor"] + by_run[m]["revisit"]
    if "KE" in S["maps"]:
        derived["KE_pass_share_of_run"] = S["maps"]["KE"]["share"]
    S["derived"] = derived
    if kr:
        stage_rows = ([k_row] if k_row else []) + kr
        stage_rows = [
            {**{k: v for k, v in r.items() if k != "wall_s"}, "wall_s_contended": r["wall_s"]}
            for r in stage_rows
        ]
        T["ar3e_stages"] = DT.markdown_table(
            stage_rows,
            [
                "run",
                "seed",
                "neighbor",
                "revisit",
                "reseed",
                "total",
                "entries",
                "index_used",
                "grid_fallbacks",
                "grid_evals",
                "neighbor_fits",
                "first_accepted",
                "retries",
                "retries_accepted",
                "wall_s_contended",
            ],
            {
                **{
                    k: "," for k in ("seed", "neighbor", "revisit", "reseed", "total", "grid_evals")
                },
                "wall_s_contended": ".0f",
            },
        )
    # ------------------------------------------------------------ second pass
    sec = []
    for m in have:
        d2 = B.RUNS / f"{m}_second"
        if not (d2 / "tries.npz").exists():
            continue
        env = envs[B.ALL_MAPS[m][0]]
        t2 = BL.Tries(d2)
        ra = after[m][pick]
        r2, tk2 = BL.apply_rule(t2, ra, pick)
        e2 = env.err(r2)
        c = BL.counts_after(errs[m][pick], e2, np.ones(env.n, bool))
        sec.append(
            dict(
                map=m,
                boundary=int(t2.info["n_boundary_all"]),
                tries=int(t2.info["n_tries"]),
                changed_again=int(tk2.sum()),
                changed_first=int(took[m][pick].sum()),
                wrong_after_first=c["wrong_before"],
                wrong_after_second=c["wrong_after"],
                fixed=c["fixed"],
                broken=c["broken"],
            )
        )
    S["second"] = sec
    # new key (review fix 11): the evaluations of each second pass (its pass.json)
    S["second_evals"] = {
        r["map"]: int(
            json.loads((B.RUNS / f"{r['map']}_second" / "pass.json").read_text())["evals"]
        )
        for r in sec
    }
    if sec:
        sec_t = [dict(r, evals_second=S["second_evals"][r["map"]]) for r in sec]
        T["ar3e_second"] = DT.markdown_table(
            sec_t,
            list(sec[0]) + ["evals_second"],
            {"boundary": ",", "tries": ",", "evals_second": ","},
        )
    # new key (addition B): the run's stored scores of the old ends against the pass's fresh ones
    S["old_end_check"] = old_end_check(TR)
    T["ar3e_old_end"] = DT.markdown_table(
        [dict(map=m, **v) for m, v in S["old_end_check"].items()],
        ["map", "n", "within_tol", "share_within_tol", "max_abs_diff"],
        {"n": ",", "within_tol": ",", "share_within_tol": ".2%", "max_abs_diff": ".4f"},
    )
    # ------------------------------------------------------------ skip rule
    if len(s4) == 3:
        S["skip"], sk = skip_study(envs, R0, E0, TR, took, errs, pick, s4)
        T.update(sk)
        S["skip_evals"], T["ar3e_skip_evals"] = skip_evals(envs, TR, S["skip"], s4)
        S["skip_cost_axis"] = dict(
            cost_quantile_picked=S["skip"]["cost_quantile"],
            cost_grid_edge=bool(S["skip"]["cost_quantile"] in (CQ_GRID[0], CQ_GRID[-1])),
            uninformative=bool(S["skip"]["cost_quantile"] in (CQ_GRID[0], CQ_GRID[-1])),
            rule="plan: a pick at the edge of its grid is called uninformative",
        )
    # keep the ratio table that boundary_ratio.py writes (run order: this script, then that one)
    tab_path = B.BND / "tables.md"
    if tab_path.exists():
        old_tabs = DT.read_tables(tab_path)
        if "ar3e_ratio" in old_tabs and "ar3e_ratio" not in T:
            T["ar3e_ratio"] = old_tabs["ar3e_ratio"]
    DT.write_tables(tab_path, T)
    (B.BND / "summary.json").write_text(json.dumps(S, indent=1, default=float) + "\n")
    figure(S, have)
    figure_pipeline(S)
    print("pick", pick, json.dumps(S["pick"]["table"])[:300])


def old_end_check(TR: Dict[str, BL.Tries]) -> Dict[str, Dict[str, Any]]:
    """Per stored .mic map: on how many voxels the run's own (confidence, cost, hit ratio) match
    the pass's fresh strict score of the same end within OLD_END_TOL (all three at once)."""
    out: Dict[str, Dict[str, Any]] = {}
    for m, rel in RECON_MIC.items():
        if m not in TR:
            continue
        p = B.C.OUT / rel
        opener: Any = gzip.open if p.suffix == ".gz" else open
        with opener(p, "rt") as f:
            lines = f.read().split("\n")[1:]
        a = np.array([[float(x) for x in ln.split()] for ln in lines if ln.strip()])
        stored = a[:, [10, 11, 9]]  # cost, hit ratio, confidence (the order of old_score)
        d = np.abs(stored - TR[m].old_score)
        ok = d.max(axis=1) <= OLD_END_TOL
        out[m] = dict(
            n=int(len(ok)),
            within_tol=int(ok.sum()),
            share_within_tol=float(ok.mean()),
            max_abs_diff=float(d.max()),
            tol=OLD_END_TOL,
        )
    return out


def skip_evals(envs, TR, info, s4):  # type: ignore
    """Evaluations of the pass with the picked skip rule: the stored evaluations of the tries
    that the rule keeps (review fix 12), per set and per map (all voxels)."""
    per: Dict[str, Any] = {}
    agg = {"A": [0, 0], "B": [0, 0]}
    rows = []
    for m in s4:
        env = envs[B.ALL_MAPS[m][0]]
        sc = TR[m].old_score
        rel = BL.rel_hit(env, sc[:, 1])
        skip = (rel >= info["rel_hit"]) & (sc[:, 0] <= info["cost_threshold"])
        ev = TR[m].evals.astype(np.int64)
        kept = ~skip[TR[m].voxel]
        for name, sel in (("A", env.isA), ("B", ~env.isA)):
            in_set = sel[TR[m].voxel]
            agg[name][0] += int(ev[in_set].sum())
            agg[name][1] += int(ev[in_set & kept].sum())
        per[m] = dict(evals_pass=int(ev.sum()), evals_kept=int(ev[kept].sum()))
        per[m]["kept_share"] = per[m]["evals_kept"] / per[m]["evals_pass"]
        rows.append(dict(scope=m, **per[m]))
    for name, (a, b) in agg.items():
        per[f"set {name}"] = dict(evals_pass=a, evals_kept=b, kept_share=b / a)
        rows.append(dict(scope=f"set {name} (3 maps)", **per[f"set {name}"]))
    tab = DT.markdown_table(
        rows,
        ["scope", "evals_pass", "evals_kept", "kept_share"],
        {"evals_pass": ",", "evals_kept": ",", "kept_share": ".1%"},
    )
    return per, tab


def skip_study(envs, R0, E0, TR, took, errs, pick, s4):  # type: ignore
    feats: Dict[str, Dict[str, np.ndarray]] = {}
    for m in s4:
        env = envs[B.ALL_MAPS[m][0]]
        sc = TR[m].old_score
        feats[m] = dict(rel=BL.rel_hit(env, sc[:, 1]), cost=sc[:, 0])
    pooled = np.concatenate([feats[m]["cost"] for m in s4])
    cq = {q: float(np.percentile(pooled, q)) for q in CQ_GRID}
    res = {}
    for part in ("A", "B"):
        for rh in RH_GRID:
            for q in CQ_GRID:
                tot_t = sk_t = fx = lost = bk_av = 0
                wrong_noskip = wrong_skip = n = 0
                for m in s4:
                    env = envs[B.ALL_MAPS[m][0]]
                    sel = env.isA if part == "A" else ~env.isA
                    f = feats[m]
                    skip = (f["rel"] >= rh) & (f["cost"] <= cq[q])
                    tv = BL.tries_per_voxel(TR[m], env.n)
                    wb, wa = E0[m] > B.WRONG, errs[m][pick] > B.WRONG
                    fixed, broken = wb & ~wa, ~wb & wa
                    tot_t += int(tv[sel].sum())
                    sk_t += int(tv[sel & skip].sum())
                    fx += int((fixed & sel).sum())
                    lost += int((fixed & sel & skip).sum())
                    bk_av += int((broken & sel & skip).sum())
                    n += int(sel.sum())
                    wrong_noskip += int((wa & sel).sum())
                    wrong_skip += (
                        int((wa & sel).sum())
                        + int((fixed & sel & skip).sum())
                        - int((broken & sel & skip).sum())
                    )
                res[(part, rh, q)] = dict(
                    tries=tot_t,
                    skipped=sk_t,
                    fixed=fx,
                    lost=lost,
                    broken_avoided=bk_av,
                    n=n,
                    wrong_noskip=wrong_noskip,
                    wrong_skip=wrong_skip,
                )
    ok = [(k, v) for k, v in res.items() if k[0] == "A" and v["lost"] <= 0.01 * v["fixed"]]
    kbest, vbest = max(ok, key=lambda kv: (kv[1]["skipped"], -kv[1]["lost"], kv[0][1]))
    _, rh, q = kbest
    edge = rh in (RH_GRID[0], RH_GRID[-1])
    info = dict(
        rel_hit=rh,
        cost_quantile=q,
        cost_threshold=cq[q],
        edge_of_grid=bool(edge),
        cost_not_used=bool(q == 100),
        cost_quantiles=cq,
    )
    rows = []
    for part in ("A", "B"):
        v = res[(part, rh, q)]
        rows.append(
            dict(
                set=part,
                tries=v["tries"],
                skipped=f"{v['skipped']} ({100 * v['skipped'] / max(v['tries'], 1):.1f}%)",
                fixed=v["fixed"],
                lost=v["lost"],
                lost_ci=wil_str(wil(v["lost"], max(v["fixed"], 1))),
                broken_avoided=v["broken_avoided"],
                wrong_noskip=v["wrong_noskip"],
                wrong_skip=v["wrong_skip"],
            )
        )
    # grid on A, rel_hit alone (cost not used)
    grid = []
    for rh_ in RH_GRID:
        v = res[("A", rh_, 100)]
        grid.append(
            dict(
                rel_hit=rh_,
                skipped=f"{100 * v['skipped'] / max(v['tries'], 1):.1f}%",
                lost=v["lost"],
                fixed=v["fixed"],
                broken_avoided=v["broken_avoided"],
            )
        )
    return info, {
        "ar3e_skip": DT.markdown_table(
            rows,
            [
                "set",
                "tries",
                "skipped",
                "fixed",
                "lost",
                "lost_ci",
                "broken_avoided",
                "wrong_noskip",
                "wrong_skip",
            ],
            {"tries": ","},
        ),
        "ar3e_skip_grid": DT.markdown_table(
            grid, ["rel_hit", "skipped", "lost", "fixed", "broken_avoided"]
        ),
    }


def figure(S: Dict[str, Any], have: List[str]) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(1, 2, figsize=(12, 4.2))
    for i, m in enumerate(have):
        c = S["maps"][m]["counts"]
        for j, (k, col) in enumerate((("wrong_before", "gray"), ("wrong_after", "tab:red"))):
            w = wil(c[k], c["n"])
            x = i + (j - 0.5) * 0.3
            p = c[k] / c["n"]
            ax[0].errorbar(
                x, 100 * p, yerr=[[100 * (p - w["lo"])], [100 * (w["hi"] - p)]], fmt="o", color=col
            )
        ax[1].bar(i, S["maps"][m]["share"] * 100, color="tab:blue")
    for a in ax:
        a.set_xticks(range(len(have)))
        a.set_xticklabels(have)
    ax[0].set_ylabel("wrong voxels (%), Wilson 95%")
    ax[0].set_title(f"Before (gray) and after (red) the pass, rule {S['pick']['picked']}")
    ax[1].set_ylabel("pass evaluations, % of the run's")
    ax[1].set_title("Cost of the pass")
    fig.tight_layout()
    fig.savefig(B.BND / "ar3_boundary.png", dpi=110)


def figure_pipeline(S: Dict[str, Any]) -> None:
    """Lead plot (addition A): wrong voxels after the pass against the evaluations of run plus
    pass, one panel per sample, every pipeline with the pass, Wilson 95% intervals as bars."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    titles = {
        "grad": "Gradient sample (step 3)",
        "layer": "Step 2's sample (one orientation per grain)",
    }
    step4 = set(B.STEP4) | {"KE"}
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.4))
    for ax, panel in zip(axes, ("grad", "layer")):
        pts = [p for p in S["pipeline_plot"] if p["panel"] == panel]
        for p in pts:
            x = p["evals"] / 1e6
            col = "tab:red" if p["map"] in step4 else "tab:gray"
            ax.errorbar(
                x,
                p["wrong"],
                yerr=[[p["wrong"] - p["lo"] * p["n"]], [p["hi"] * p["n"] - p["wrong"]]],
                fmt="o",
                color=col,
                capsize=3,
            )
            dx, dy, ha = LABEL_POS.get(p["map"], (8, 6, "left"))
            ax.annotate(
                f"{p['map']} + pass\n{p['wrong']} wrong, {x:.1f}M",
                (x, p["wrong"]),
                textcoords="offset points",
                xytext=(dx, dy),
                ha=ha,
                va="center",
                fontsize=9,
                color=col,
                bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="none", alpha=0.85),
            )
        xmax = max(p["evals"] for p in pts) / 1e6
        ax.set_xlim(0, xmax * 1.35)
        ymax = max(p["hi"] * p["n"] for p in pts)
        ax.set_ylim(-0.05 * ymax - 1, ymax * 1.3 + 2)
        ax.set_xlabel("evaluations of run plus pass (millions)")
        ax.set_ylabel("wrong voxels after the pass (of 24,570), Wilson 95%")
        ax.set_title(titles[panel])
    fig.suptitle(
        "Every pipeline with the same pass: red = step 4 runs and KE, gray = runs they are "
        "compared with",
        fontsize=10,
    )
    fig.tight_layout()
    fig.savefig(B.BND / "ar3_boundary_pipeline.png", dpi=110)


if __name__ == "__main__":
    main()
