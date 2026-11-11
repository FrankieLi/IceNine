"""AR5 step 2: summary.json, tables.md and the plots of the fit study.

  uv run python scripts/autoresearch/ar5/fit_summary.py

Reads only files that are small enough to commit: ``part0_<sample>.json``, ``pick.json``,
``part1_mid.json`` (written here from the cached replay of ``mid``, when it exists),
``base_check.json`` and ``part2_<name>_<sample>.json``. Anything not run yet is written as
"not evaluated". The Part 1 numbers of ``big`` and ``grad`` are recomputed here from the cached
replay shards (``cache/replay_*``) so the tables follow the same arrays as the pick.
Writes ``summary.json``, ``tables.md`` (blocks for ``scripts/dev/sync_doc_tables.py``) and the
plots ``ar5_fit_stages.png`` and ``ar5_fit_arms.png`` under ``benchmarks/autoresearch/ar5/fit/``.
"""

import json
import os
from typing import Any, Dict, List, Optional

os.environ.setdefault("RESEED_SAMPLE", "big")

import fit_common as FC  # noqa: E402
import fit_part1_lib as L  # noqa: E402
import fit_rules as FR  # noqa: E402
import matplotlib  # noqa: E402
import numpy as np  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from doc_tables import markdown_table, write_tables  # noqa: E402
from stats import wilson  # noqa: E402

FIT = FC.FIT
SAMPLES = ("big", "grad", "mid")
NAMES = {"big": "big", "grad": "first gradient sample", "mid": "mid"}
RUN_NAMES = ("S10", "S10a", "S10b")
STAGES = ("seed_step", "grid_fallback", "neighbour", "revisit", "pass_")
STAGE_LABEL = {
    "seed_step": "seed step",
    "grid_fallback": "grid fallback",
    "neighbour": "neighbour fits",
    "revisit": "revisits",
    "pass_": "boundary pass",
}


def rj(name: str) -> Optional[Dict[str, Any]]:
    p = FIT / name
    return json.loads(p.read_text()) if p.exists() else None


def ci(k: int, n: int) -> str:
    lo, hi = wilson(k, n)
    return f"{100 * lo:.1f} to {100 * hi:.1f}%"


HEADERS: Dict[str, Dict[str, str]] = {
    "part0": dict(
        steps="seed steps",
        refs="refinements",
        ref_evals="refinement evaluations",
        ref_mean="mean per refinement",
        fb_runs="MC fallbacks",
        fb_share="share of fits (%)",
        fb_evals="fallback evaluations",
        fb_mean="mean per fallback",
        fb_won="ended lower than CMA",
        fb_won_ci="Wilson 95% interval",
        total="evaluations of S8 with the pass",
        share="refinements + fallbacks (% of total)",
    ),
    "part1": dict(
        passed="pass the seed test",
        passed_ci="Wilson 95% interval",
        m_only_p="passes: M only",
        arm_only_p="passes: arm only",
        right="within 1 degree",
        m_only_r="within 1 deg: M only",
        arm_only_r="within 1 deg: arm only",
        evals="seed-step evaluations",
        fall="fewer than M (%)",
    ),
    "err": {
        "steps": "steps with candidates",
        "median": "median error (deg)",
        "p90": "90th percentile (deg)",
        "within": "within 0.2 deg (%)",
    },
    "pick": dict(
        passed="pass the seed test",
        right="within 1 degree",
        evals="seed-step evaluations",
        a="rule (a)",
        b="rule (b)",
        chosen="picked",
    ),
    "mid_pass": dict(a="rule (a)", b="rule (b)", fall="evaluations fall (%)", c="rule (c)"),
    "part2": dict(
        seed="seed step",
        grid="grid fallback",
        neigh="neighbour fits",
        rev="revisits",
        pas="boundary pass",
        total="total",
        gfb="grid fallbacks (count)",
        wrong_b="wrong before the pass",
        wrong_a="wrong after the pass",
        fix_b="fixed / broken before",
        fix_a="fixed / broken after",
        med="median error after (deg)",
    ),
    "fits": dict(
        s8_fits="S8 fits",
        s8_fb="S8 MC fallbacks",
        s8_won="S8 MC ended lower",
        s10_fits="S10 fits",
        s10_nofall="S10 CMA cost did not fall",
        s8_neigh_revisit="S8 neighbour + revisit evaluations",
        s10_neigh_revisit="S10 neighbour + revisit evaluations",
        seeds8="S8 seed steps",
        seeds10="S10 seed steps",
    ),
    "rule": dict(
        fewer="fewer than S8 (%)",
        a="rule (a)",
        net="broken minus fixed",
        b="rule (b)",
        meddiff="median error minus S8's (deg)",
        c="rule (c)",
        passes="passes",
    ),
    "pred": dict(n="prediction", obs="observed", verdict="verdict"),
    "shares": dict(
        t8="S8 total",
        t10="S10 total",
        saved="saved",
        seed="S8 seed step (% of total)",
        grid="S8 grid fallback (%)",
        nr="S8 neighbour + revisit (%)",
        pas="S8 boundary pass (%)",
        nr_saved="saved in neighbour + revisit (% of S8 total)",
    ),
}


def pretty(name: str, md: str) -> str:
    """Rename the header cells of a markdown table (the keys used while building it)."""
    lines = md.split("\n")
    m = HEADERS.get(name, {})
    cells = [c.strip() for c in lines[0].strip("|").split("|")]
    lines[0] = "| " + " | ".join(m.get(c, c) for c in cells) + " |"
    return "\n".join(lines)


def wrong_text(side: Dict[str, Any]) -> str:
    """'count (Wilson 95% interval)' of a scored side."""
    lo, hi = 100 * side["wilson_lo"], 100 * side["wilson_hi"]
    return f"{side['wrong']} ({lo:.2f} to {hi:.2f}%)"


def part1_counts(sample: str) -> Optional[Dict[str, Any]]:
    """Per-sample arm counts. mid is read from ``part1_mid.json`` (written by fit_part1_mid.py)."""
    if sample == "mid":
        return rj("part1_mid.json")
    rows = L.load(sample)
    if not rows:
        return None
    return dict(
        n_steps=len(rows),
        counts=L.counts(rows),
        paired_passed={a: L.paired(rows, a, "passed") for a in FC.ARMS if a != "M"},
        paired_right={a: L.paired(rows, a, "right") for a in FC.ARMS if a != "M"},
        log_mismatches_M=sum(1 for r in rows if not r["log_ok"]),
        err_q=L.err_quantiles(rows),
    )


def main() -> None:
    S: Dict[str, Any] = dict(samples=list(SAMPLES))
    S["part0"] = {s: rj(f"part0_{s}.json") for s in SAMPLES}
    S["pick"] = rj("pick.json")
    S["part1"] = {s: part1_counts(s) for s in SAMPLES}
    S["part1_mid_pass"] = rj("part1_mid_pass.json")
    S["base_check"] = rj("base_check.json")
    S["part2"] = {n: {s: rj(f"part2_{n}_{s}.json") for s in SAMPLES} for n in RUN_NAMES}
    T: Dict[str, str] = {}
    pk0 = S["pick"]
    if pk0:
        c = pk0["counts_pooled"]
        S["derived"] = dict(
            net_fewer_passes={a: c["M"]["n_pass"] - c[a]["n_pass"] for a in FC.ARMS if a != "M"},
            net_fewer_right={a: c["M"]["n_right"] - c[a]["n_right"] for a in FC.ARMS if a != "M"},
            replayed_steps_all_samples=sum(
                S["part1"][s]["n_steps"] for s in SAMPLES if S["part1"][s]
            ),
            sample_grains=dict(big=32, grad=497, mid=120),
            voxels_per_sample=24570,
        )

    # ---- Part 0
    rows0: List[Optional[Dict[str, Any]]] = []
    labels = []
    for s in SAMPLES:
        d = S["part0"][s]
        labels.append(NAMES[s])
        if d is None:
            rows0.append(None)
            continue
        lo, hi = wilson(d["mc_fallback_won"], d["mc_fallback_runs"])
        rows0.append(
            dict(
                sample=NAMES[s],
                steps=d["seed_steps"],
                refs=d["seed_refinements"],
                ref_evals=d["seed_refinement_evals"],
                ref_mean=d["seed_refinement_evals_mean"],
                fb_runs=d["mc_fallback_runs"],
                fb_share=100 * d["mc_fallback_share"],
                fb_evals=d["mc_fallback_evals"],
                fb_mean=d["mc_fallback_evals_mean"],
                fb_won=d["mc_fallback_won"],
                fb_won_ci=f"{100 * lo:.2f} to {100 * hi:.2f}%",
                total=d["evals_with_pass"],
                share=100 * d["share_both"],
            )
        )
    T["part0"] = markdown_table(
        rows0,
        [
            "sample",
            "steps",
            "refs",
            "ref_evals",
            "ref_mean",
            "fb_runs",
            "fb_share",
            "fb_evals",
            "fb_mean",
            "fb_won",
            "fb_won_ci",
            "total",
            "share",
        ],
        formats=dict(ref_mean=".0f", fb_share=".1f", fb_mean=".0f", share=".1f", total="d"),
        labels=[NAMES[s] for s in SAMPLES],
    )

    # ---- Part 1 counts
    rows1: List[Optional[Dict[str, Any]]] = []
    lab1: List[str] = []
    for s in SAMPLES:
        d = S["part1"][s]
        for arm in FC.ARMS:
            lab1.append(f"{NAMES[s]} {arm}")
            if d is None:
                rows1.append(None)
                continue
            c, m, n = d["counts"][arm], d["counts"]["M"], d["n_steps"]
            pp = d["paired_passed"].get(arm)
            pr = d["paired_right"].get(arm)
            rows1.append(
                dict(
                    sample=NAMES[s],
                    arm=arm,
                    steps=n,
                    passed=c["n_pass"],
                    passed_ci=ci(c["n_pass"], n),
                    m_only_p=pp["m_only"] if pp else 0,
                    arm_only_p=pp["arm_only"] if pp else 0,
                    right=c["n_right"],
                    m_only_r=pr["m_only"] if pr else 0,
                    arm_only_r=pr["arm_only"] if pr else 0,
                    evals=c["evals"],
                    fall=100 * FR.fall(c["evals"], m["evals"]),
                )
            )
    T["part1"] = markdown_table(
        rows1,
        [
            "sample",
            "arm",
            "steps",
            "passed",
            "passed_ci",
            "m_only_p",
            "arm_only_p",
            "right",
            "m_only_r",
            "arm_only_r",
            "evals",
            "fall",
        ],
        formats=dict(fall=".1f", evals="d"),
        labels=lab1,
    )

    # ---- end errors of the arms
    rows_e: List[Optional[Dict[str, Any]]] = []
    lab_e: List[str] = []
    for s in SAMPLES:
        d = S["part1"][s]
        for arm in FC.ARMS:
            lab_e.append(f"{NAMES[s]} {arm}")
            if d is None or "err_q" not in d:
                rows_e.append(None)
                continue
            q = d["err_q"][arm]
            rows_e.append(
                dict(
                    sample=NAMES[s],
                    arm=arm,
                    steps=q["n"],
                    median=q["median"],
                    p90=q["p90"],
                    within=100 * q["within_0_2"],
                )
            )
    T["err"] = markdown_table(
        rows_e,
        ["sample", "arm", "steps", "median", "p90", "within"],
        formats=dict(median=".3f", p90=".3f", within=".1f"),
        labels=lab_e,
    )
    # ---- pick
    pk = S["pick"]
    if pk:
        rp = []
        for arm in FC.ARMS:
            c = pk["counts_pooled"][arm]
            rp.append(
                dict(
                    arm=arm,
                    passed=c["n_pass"],
                    right=c["n_right"],
                    evals=c["evals"],
                    a="reference" if arm == "M" else str(pk["rule_a"][arm]),
                    b="reference" if arm == "M" else str(pk["rule_b"][arm]),
                    chosen="yes" if arm == pk["pick"] else "no",
                )
            )
        T["pick"] = markdown_table(
            rp, ["arm", "passed", "right", "evals", "a", "b", "chosen"], formats=dict(evals="d")
        )
    # ---- mid pass rule
    mp = S["part1_mid_pass"]
    if mp:
        rm = []
        for arm, d in mp["arms"].items():
            rm.append(
                dict(
                    arm=arm,
                    a=str(d["a"]),
                    b=str(d["b"]),
                    fall=100 * d["evals_fall"],
                    c=str(d["c"]),
                    passes=str(d["passes"]),
                )
            )
        T["mid_pass"] = markdown_table(
            rm, ["arm", "a", "b", "fall", "c", "passes"], formats=dict(fall=".1f")
        )

    # ---- Part 2
    rows2: List[Optional[Dict[str, Any]]] = []
    lab2: List[str] = []
    for n in ("S8",) + RUN_NAMES:
        if n != "S8" and not any(S["part2"][n][s] for s in SAMPLES):
            continue
        for s in SAMPLES:
            d = S["part2"][n][s] if n != "S8" else None
            if n == "S8":
                # stored S8 and its pass, from any finished run's file (the stored side)
                src = next((S["part2"][k][s] for k in RUN_NAMES if S["part2"][k][s]), None)
                if not src or src.get("unfinished"):
                    continue
                st, b, a = src["stages_s8"], src["s8_before"], src["s8_after"]
                gf = src["grid_fallbacks_s8"]
                pb = pa = None
            else:
                if d is None or d.get("unfinished"):
                    lab2.append(f"{n} {NAMES[s]}")
                    rows2.append(None)
                    continue
                st, b, a, gf = d["stages"], d["before"], d["after"], d["grid_fallbacks"]
                pb, pa = d["paired_before"], d["paired_after"]
            lab2.append(f"{n} {NAMES[s]}")
            rows2.append(
                dict(
                    run=n,
                    sample=NAMES[s],
                    seed=st["seed_step"],
                    grid=st["grid_fallback"],
                    neigh=st["neighbour"],
                    rev=st["revisit"],
                    pas=st["pass_"],
                    total=st["total"],
                    gfb=gf,
                    wrong_b=wrong_text(b),
                    wrong_a=wrong_text(a),
                    fix_b="" if pb is None else f"{pb['fixed']} / {pb['broken']}",
                    fix_a="" if pa is None else f"{pa['fixed']} / {pa['broken']}",
                    med=a["median_err_deg"],
                )
            )
    T["part2"] = markdown_table(
        rows2,
        [
            "run",
            "sample",
            "seed",
            "grid",
            "neigh",
            "rev",
            "pas",
            "total",
            "gfb",
            "wrong_b",
            "wrong_a",
            "fix_b",
            "fix_a",
            "med",
        ],
        formats=dict(seed="d", grid="d", neigh="d", rev="d", pas="d", total="d", med=".4f"),
        labels=lab2,
    )
    # ---- what the fallback did in S8 and what S10 saved
    rf: List[Optional[Dict[str, Any]]] = []
    lf: List[str] = []
    for s in SAMPLES:
        d, p0 = S["part2"]["S10"][s], S["part0"][s]
        lf.append(NAMES[s])
        if d is None or d.get("unfinished") or p0 is None:
            rf.append(None)
            continue
        w = d["wrapper"]
        rf.append(
            dict(
                sample=NAMES[s],
                s8_fits=p0["neighbour_fit_calls"],
                s8_fb=p0["mc_fallback_runs"],
                s8_won=p0["mc_fallback_won"],
                s10_fits=w["wrapper_calls"],
                s10_nofall=w["cma_not_fell"],
                s8_neigh_revisit=d["stages_s8"]["neighbour"] + d["stages_s8"]["revisit"],
                s10_neigh_revisit=d["stages"]["neighbour"] + d["stages"]["revisit"],
                seeds8=p0["seed_steps"],
                seeds10=d["seed_steps"],
            )
        )
    T["fits"] = markdown_table(
        rf,
        [
            "sample",
            "s8_fits",
            "s8_fb",
            "s8_won",
            "s10_fits",
            "s10_nofall",
            "s8_neigh_revisit",
            "s10_neigh_revisit",
            "seeds8",
            "seeds10",
        ],
        formats=dict(s8_neigh_revisit="d", s10_neigh_revisit="d"),
        labels=lf,
    )
    # ---- where S8's evaluations are, and what S10 saved
    rs: List[Optional[Dict[str, Any]]] = []
    ls: List[str] = []
    for s in SAMPLES:
        d = S["part2"]["S10"][s]
        ls.append(NAMES[s])
        if d is None or d.get("unfinished"):
            rs.append(None)
            continue
        a8, a10 = d["stages_s8"], d["stages"]
        t8 = a8["total"]
        rs.append(
            dict(
                sample=NAMES[s],
                t8=t8,
                t10=a10["total"],
                saved=t8 - a10["total"],
                seed=100 * a8["seed_step"] / t8,
                grid=100 * a8["grid_fallback"] / t8,
                nr=100 * (a8["neighbour"] + a8["revisit"]) / t8,
                pas=100 * a8["pass_"] / t8,
                nr_saved=100
                * (a8["neighbour"] + a8["revisit"] - a10["neighbour"] - a10["revisit"])
                / t8,
            )
        )
    T["shares"] = markdown_table(
        rs,
        ["sample", "t8", "t10", "saved", "seed", "grid", "nr", "pas", "nr_saved"],
        formats=dict(
            t8="d", t10="d", saved="d", seed=".1f", grid=".1f", nr=".1f", pas=".1f", nr_saved=".1f"
        ),
        labels=ls,
    )
    # ---- rule verdicts
    rr = []
    rl = []
    for n in RUN_NAMES:
        for s in SAMPLES:
            d = S["part2"][n][s]
            if d is None or d.get("unfinished"):
                continue
            r = d["rule"]
            rl.append(f"{n} {NAMES[s]}")
            rr.append(
                dict(
                    run=n,
                    sample=NAMES[s],
                    fewer=100 * r["fewer"],
                    a=str(r["a"]),
                    net=r["net_broken"],
                    b=str(r["b"]),
                    meddiff=r["med_diff"],
                    c=str(r["c"]),
                    passes=str(r["passes"]),
                )
            )
    if rr:
        T["rule"] = markdown_table(
            rr,
            ["run", "sample", "fewer", "a", "net", "b", "meddiff", "c", "passes"],
            formats=dict(fewer=".1f", meddiff="+.4f"),
            labels=rl,
        )
    bc = S["base_check"]
    if bc:
        T["base"] = markdown_table(
            [
                dict(
                    item="evaluations before the pass",
                    new=bc["evals_total"][0],
                    stored=bc["evals_total"][1],
                ),
                dict(
                    item="wrong voxels before the pass",
                    new=bc["wrong_new"],
                    stored=bc["wrong_stored"],
                ),
                dict(
                    item="seed steps",
                    new=bc["seed_steps"][0],
                    stored=bc["seed_steps"][1],
                ),
            ],
            ["item", "new", "stored"],
            formats=dict(new="d", stored="d"),
        )
    T["pred"] = predictions(S)
    T = {k: pretty(k, v) for k, v in T.items()}
    write_tables(FIT / "tables.md", T)
    (FIT / "summary.json").write_text(json.dumps(S, indent=1, default=float) + "\n")
    plots(S)
    print("summary written;", sorted(T))


def predictions(S: Dict[str, Any]) -> str:
    """The five predictions of the plan, each with the observed numbers and a verdict."""
    rows: List[Dict[str, Any]] = []
    NE = "not evaluated"

    def done(n: str, s: str) -> Optional[Dict[str, Any]]:
        d = S["part2"][n][s]
        return d if d and not d.get("unfinished") else None

    # 1: the MC fallback ended lower than the CMA end in under 1% of its runs, every sample
    shares = []
    ok1 = True
    for s in SAMPLES:
        p = S["part0"][s]
        if p is None:
            ok1 = None
            continue
        sh = 100 * p["mc_fallback_won"] / p["mc_fallback_runs"]
        shares.append(f"{NAMES[s]} {p['mc_fallback_won']} of {p['mc_fallback_runs']} ({sh:.2f}%)")
        ok1 = ok1 and sh < 1.0 if ok1 is not None else ok1
    rows.append(
        dict(n=1, obs="; ".join(shares), verdict=NE if ok1 is None else "agrees" if ok1 else "no")
    )
    # 2: C60 is the pick, passes on mid, and cuts the seed-step evaluations by 90% or more
    pk = S["pick"]
    cuts = []
    ok_cut = True
    for s in SAMPLES:
        d = S["part1"][s]
        if d is None:
            ok_cut = None
            continue
        f = 100 * FR.fall(d["counts"]["C60"]["evals"], d["counts"]["M"]["evals"])
        cuts.append(f"{NAMES[s]} {f:.1f}%")
        ok_cut = ok_cut and f >= 90.0 if ok_cut is not None else ok_cut
    v2 = "no (the pick is " + (pk["pick"] if pk else "?") + ")"
    rows.append(
        dict(
            n=2,
            obs=f"pick {pk['pick'] if pk else NE}; C60 cuts the seed step by " + ", ".join(cuts),
            verdict=v2
            + ("; the cut holds" if ok_cut else "; the cut fails" if ok_cut is False else ""),
        )
    )
    # 3: C0 fails rule (a)
    if pk:
        c0 = pk["counts_pooled"]["C0"]["n_pass"]
        m0 = pk["counts_pooled"]["M"]["n_pass"]
        rows.append(
            dict(
                n=3,
                obs=f"C0 passes the seed test in {c0} of {pk['n_steps']} steps, M in {m0}",
                verdict="agrees" if not pk["rule_a"]["C0"] else "no",
            )
        )
    # 4: S10 with the pass under 4.0 / 14 / 13.5 million and passes its rule on all three
    lim = dict(big=4.0e6, grad=14.0e6, mid=13.5e6)
    obs4 = []
    ok4: Optional[bool] = True
    for s in SAMPLES:
        d = done("S10", s)
        if d is None:
            ok4 = None
            obs4.append(f"{NAMES[s]} {NE}")
            continue
        t = d["stages"]["total"]
        obs4.append(f"{NAMES[s]} {t} (limit {int(lim[s])}; rule passes {d['rule']['passes']})")
        if ok4 is not None:
            ok4 = ok4 and t < lim[s] and bool(d["rule"]["passes"])
    rows.append(
        dict(n=4, obs="; ".join(obs4), verdict=NE if ok4 is None else "agrees" if ok4 else "no")
    )
    # 5: grid fallbacks of S10 within 2 of S8's
    obs5 = []
    ok5: Optional[bool] = True
    for s in SAMPLES:
        d = done("S10", s)
        if d is None:
            ok5 = None
            obs5.append(f"{NAMES[s]} {NE}")
            continue
        obs5.append(f"{NAMES[s]} {d['grid_fallbacks']} against {d['grid_fallbacks_s8']}")
        if ok5 is not None:
            ok5 = ok5 and abs(d["grid_fallbacks"] - d["grid_fallbacks_s8"]) <= 2
    rows.append(
        dict(n=5, obs="; ".join(obs5), verdict=NE if ok5 is None else "agrees" if ok5 else "no")
    )
    return markdown_table(rows, ["n", "obs", "verdict"])


def plots(S: Dict[str, Any]) -> None:
    # plot 1: evaluations by stage, S8 and each finished S10 run, three samples
    runs = [
        n
        for n in RUN_NAMES
        if any(S["part2"][n][s] and not S["part2"][n][s].get("unfinished") for s in SAMPLES)
    ]
    if runs:
        fig, axes = plt.subplots(1, 3, figsize=(13, 4.2), sharey=False)
        cols = dict(
            seed_step="#4c72b0",
            grid_fallback="#c44e52",
            neighbour="#55a868",
            revisit="#dd8452",
            pass_="#8172b3",
        )
        for ax, s in zip(axes, SAMPLES):
            labels, data = [], []
            ref = next(
                (
                    S["part2"][n][s]
                    for n in runs
                    if S["part2"][n][s] and not S["part2"][n][s].get("unfinished")
                ),
                None,
            )
            if ref:
                labels.append("S8")
                data.append(ref["stages_s8"])
            for n in runs:
                d = S["part2"][n][s]
                if d and not d.get("unfinished"):
                    labels.append(n)
                    data.append(d["stages"])
            bottom = np.zeros(len(data))
            for k in STAGES:
                v = np.array([d[k] for d in data], float) / 1e6
                ax.bar(labels, v, bottom=bottom, color=cols[k], label=STAGE_LABEL[k])
                bottom += v
            ax.set_title(NAMES[s])
            ax.set_ylabel("evaluations (millions)")
        axes[0].legend(fontsize=8)
        fig.suptitle("Evaluations by stage, with the boundary pass")
        fig.tight_layout()
        fig.savefig(FIT / "ar5_fit_stages.png", dpi=120)
        plt.close(fig)
    # plot 2: seed refinement arms, evaluations against the share of steps within 1 degree
    fig, axes = plt.subplots(1, 3, figsize=(13, 4.0))
    for ax, s in zip(axes, SAMPLES):
        d = S["part1"][s]
        ax.set_title(NAMES[s])
        if not d:
            ax.text(0.5, 0.5, "not evaluated", ha="center", transform=ax.transAxes)
            continue
        n = d["n_steps"]
        for arm in FC.ARMS:
            c = d["counts"][arm]
            ax.scatter(max(c["evals"], 1) / 1e6, 100 * c["n_right"] / n, s=60, label=arm)
            ax.scatter(max(c["evals"], 1) / 1e6, 100 * c["n_pass"] / n, s=60, marker="x", color="k")
        ax.set_xscale("log")
        ax.set_xlabel("seed-step evaluations (millions, log)")
        ax.set_ylabel("share of steps (%)")
    axes[0].legend(fontsize=8, title="dot: within 1 degree; cross: seed test")
    fig.tight_layout()
    fig.savefig(FIT / "ar5_fit_arms.png", dpi=120)
    plt.close(fig)


if __name__ == "__main__":
    main()
