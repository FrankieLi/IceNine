"""AR3 step 5: summary of part 0, part 1 (the pass on the stored maps) and part 2 (runs K, KE).

  uv run python scripts/autoresearch/ar3/boundary_summary.py

Reads ``boundary/part0.json``, ``boundary/runs/<map>/tries.npz`` and ``pass.json`` (and
``<map>_second`` when present), the run.json files of the compared runs. Writes
``boundary/summary.json``, ``boundary/tables.md``, ``boundary/ar3_boundary.png`` and
``boundary/cache/after_<map>.npy`` (the map after the picked rule, input of the second pass).
Truth is used to score and to pick the rule on set A (even true-grain id), reported on B (odd).
Wilson 95% intervals, no p-values; a gain or a loss only when intervals do not overlap.
One run, one voxel order.
"""

import json
from typing import Any, Dict, List, Optional

import boundary_common as B
import boundary_lib as BL
import numpy as np
from fine_summary import overlap, wil, wil_str

import doc_tables as DT  # noqa: E402

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
    T["ar3e_part0"] = DT.markdown_table(
        rows,
        [
            "map",
            "boundary",
            "cands",
            "wrong",
            "wrong_cand1",
            "tr_cost",
            "tr_hit",
            "tr_conf",
            "newly",
            "newly_both",
        ],
        {k: "," for k in ("boundary", "cands")},
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
    T["ar3e_pipeline"] = DT.markdown_table(
        pipe,
        [
            "this",
            "against",
            "evals_this",
            "evals_against",
            "wrong_this",
            "wrong_against",
            "fixed",
            "broken",
            "both_wrong",
            "verdict",
        ],
        {"evals_this": ",", "evals_against": ","},
    )
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
    if kr:
        T["ar3e_stages"] = DT.markdown_table(
            kr,
            ["run", "seed", "neighbor", "revisit", "reseed", "total", "entries", "index_used",
             "grid_fallbacks", "grid_evals", "neighbor_fits", "first_accepted", "retries",
             "retries_accepted", "wall_s"],
            {k: "," for k in ("seed", "neighbor", "revisit", "reseed", "total", "grid_evals")}
            | {"wall_s": ".0f"},
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
    if sec:
        T["ar3e_second"] = DT.markdown_table(sec, list(sec[0]), {"boundary": ",", "tries": ","})
    # ------------------------------------------------------------ skip rule
    if len(s4) == 3:
        S["skip"], sk = skip_study(envs, R0, E0, TR, took, errs, pick, s4)
        T.update(sk)
    DT.write_tables(B.BND / "tables.md", T)
    (B.BND / "summary.json").write_text(json.dumps(S, indent=1, default=float) + "\n")
    figure(S, have)
    print("pick", pick, json.dumps(S["pick"]["table"])[:300])


def skip_study(envs, R0, E0, TR, took, errs, pick, s4):  # type: ignore
    feats: Dict[str, Dict[str, np.ndarray]] = {}
    for m in s4:
        env = envs[B.ALL_MAPS[m][0]]
        sc = TR[m].old_score
        feats[m] = dict(rel=BL.rel_hit(env, sc[:, 1]), cost=sc[:, 0])
    pooled = np.concatenate([feats[m]["cost"] for m in s4])
    cq = {q: float(np.percentile(pooled, q)) for q in CQ_GRID}
    res = {}
    out_rows: Dict[str, List[Dict[str, Any]]] = {"A": [], "B": []}
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


if __name__ == "__main__":
    main()
