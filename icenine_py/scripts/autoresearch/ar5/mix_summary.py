"""AR5 step 3: tables, plots and ``summary.json`` of the mixed seed refinement.

  uv run python scripts/autoresearch/ar5/mix_summary.py

Reads ``pick.json``, ``part1.json``, ``check.json`` and ``part2_<run>_<sample>.json`` in
``benchmarks/autoresearch/ar5/mix/`` and the numbers of runs S8 and S10 in
``benchmarks/autoresearch/ar5/fit/``; writes ``tables.md`` (tables in the marker format of
``doc_tables``; ``scripts/dev/sync_doc_tables.py`` copies them into the findings page),
``summary.json`` and two plots. Every number of the page comes from these files.
"""

import json
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

HERE = Path(__file__).resolve().parent
ICE = HERE.parents[2]
sys.path.insert(0, str(ICE / "scripts" / "common"))
MIX = ICE / "benchmarks" / "autoresearch" / "ar5" / "mix"
FIT = ICE / "benchmarks" / "autoresearch" / "ar5" / "fit"

from doc_tables import markdown_table, write_tables  # noqa: E402

SAMPLES = ("big", "grad", "mid")
NAMES = {"big": "big", "grad": "first gradient sample", "mid": "mid"}
ARMS = ("M", "X60", "X200")
STAGES = ("seed_step", "grid_fallback", "neighbour", "revisit", "pass_")
STAGE_LABEL = dict(
    seed_step="seed step",
    grid_fallback="grid fallback",
    neighbour="neighbour fits",
    revisit="revisits",
    pass_="boundary pass",
)
PRED2 = dict(big=4.4e6, grad=14.2e6, mid=17.5e6)  # prediction 2: S11 totals below these


def load(p: Path) -> Optional[Dict[str, Any]]:
    return json.loads(p.read_text()) if p.exists() else None


def pct(x: float) -> float:
    return 100.0 * x


def ci(side: Dict[str, Any]) -> str:
    return f"{pct(side['wilson_lo']):.2f} to {pct(side['wilson_hi']):.2f}%"


def fmt_int(n: int) -> str:
    return f"{n:,}"


def derived() -> Dict[str, Any]:
    """Numbers derived for the review fixes (D-20261111-N16): the share of S11's saving by stage
    against S8 and S10, the cost of one grid fallback, the uncounted hit-ratio checks of the
    mixed seed refinement, and the plan's what-if step (from part1.json)."""
    shares: Dict[str, Any] = {}
    fallback: Dict[str, Any] = {}
    unc: Dict[str, Any] = {}
    for s in SAMPLES:
        d10 = load(FIT / f"part2_S10_{s}.json")
        d11 = load(MIX / f"part2_S11_{s}.json")
        if not d11 or d11.get("unfinished"):
            continue
        st11 = d11["stages"]
        refs = dict(S8=(d10["stages_s8"], d10["grid_fallbacks_s8"]))
        refs["S10"] = (d10["stages"], d10["grid_fallbacks"])
        shares[s] = {}
        for ref, (st, _) in refs.items():
            save = st["total"] - st11["total"]
            shares[s][ref] = dict(
                total_saving=int(save),
                **{k: float((st[k] - st11[k]) / save) for k in STAGES},
            )
        per = {}
        for name, (st, g) in (*refs.items(), ("S11", (st11, d11["grid_fallbacks"]))):
            per[name] = dict(
                grid_evals=int(st["grid_fallback"]),
                fallbacks=int(g),
                per_fallback=(float(st["grid_fallback"] / g) if g else None),
            )
        fallback[s] = per
        n = int(d11["mc_share"]["candidates"])
        unc[s] = dict(
            mixed_candidates=n, total_evals=int(st11["total"]), share=float(n / st11["total"])
        )
    p1 = load(MIX / "part1.json")
    return dict(
        stage_shares=shares,
        grid_fallback_cost=fallback,
        uncounted_hit_checks=unc,
        plan_whatif=p1.get("plan_whatif") if p1 else None,
    )


def tables() -> Dict[str, str]:
    p1 = load(MIX / "part1.json")
    pk = load(MIX / "pick.json")
    ck = load(MIX / "check.json")
    T: Dict[str, str] = {}
    # ---- pick (pooled over the pick samples)
    rows = []
    pooled = pk["pooled"]["arms"]
    for a in ARMS:
        c = pooled[a]
        rows.append(
            dict(
                arm=a,
                passes=c["n_pass"],
                right=c["n_right"],
                evals=c["evals"],
                fewer=pct(1 - c["evals"] / pooled["M"]["evals"]),
                mc=pct(c["share_mc_ran"]),
                missing=c["steps_missing"],
                ra=("n/a" if a == "M" else str(pk["rule_a"][a])),
                rb=("n/a" if a == "M" else str(pk["rule_b"][a])),
            )
        )
    T["pick"] = markdown_table(
        rows,
        ["arm", "passes", "right", "evals", "fewer", "mc", "missing", "ra", "rb"],
        {"fewer": ".1f", "mc": ".1f"},
    ).replace(
        "| arm | passes | right | evals | fewer | mc | missing | ra | rb |",
        "| arm | seed-test passes | ends within 1 degree | seed-step evaluations | fewer than M (%)"
        " | candidates where MC ran (%) | steps with a missing MC end | rule (a) | rule (b) |",
    )
    # ---- part 1 per sample
    rows = []
    for s in SAMPLES:
        d = p1["by_sample"][s]
        n = d["n_steps"]
        for a in ARMS:
            c = d["arms"][a]
            ref = d["arms"]["M"]
            row = dict(
                sample=NAMES[s],
                arm=a,
                steps=n,
                passes=c["n_pass"],
                p_m=("n/a" if a == "M" else c["paired_passed"]["m_only"]),
                p_a=("n/a" if a == "M" else c["paired_passed"]["arm_only"]),
                right=c["n_right"],
                r_m=("n/a" if a == "M" else c["paired_right"]["m_only"]),
                r_a=("n/a" if a == "M" else c["paired_right"]["arm_only"]),
                evals=c["evals"],
                fewer=pct(1 - c["evals"] / ref["evals"]),
                mc=pct(c["share_mc_ran"]),
                missing=c["steps_missing"],
            )
            rows.append(row)
    T["part1"] = markdown_table(
        rows,
        [
            "sample",
            "arm",
            "steps",
            "passes",
            "p_m",
            "p_a",
            "right",
            "r_m",
            "r_a",
            "evals",
            "fewer",
            "mc",
            "missing",
        ],  # fmt: skip
        {"fewer": ".1f", "mc": ".1f"},
    ).replace(
        "| sample | arm | steps | passes | p_m | p_a | right | r_m | r_a | evals | fewer | mc |"
        " missing |",
        "| sample | arm | steps | seed-test passes | passes with M only | passes with the arm only"
        " | ends within 1 degree | M only | arm only | seed-step evaluations | fewer than M (%)"
        " | candidates where MC ran (%) | steps with a missing MC end |",
    )
    # ---- check
    T["check"] = markdown_table(
        [
            dict(
                steps=ck["n_steps"],
                withmc=ck["n_with_mc"],
                same=ck["n_same"],
                live=ck["evals_live"],
                cache=ck["evals_cache"],
            )
        ],
        ["steps", "withmc", "same", "live", "cache"],
    ).replace(
        "| steps | withmc | same | live | cache |",
        "| seed steps checked | with MC run on a candidate | identical to Part 1 | evaluations,"
        " script layer | evaluations, Part 1 |",
    )
    # ---- whole-layer runs
    runs: Dict[str, Dict[str, Any]] = {s: {} for s in SAMPLES}
    for s in SAMPLES:
        for r in ("S11", "S11a"):
            runs[s][r] = load(MIX / f"part2_{r}_{s}.json")
        runs[s]["S10"] = load(FIT / f"part2_S10_{s}.json")
    rows = []
    for s in SAMPLES:
        d10 = runs[s]["S10"]
        base = d10["stages_s8"]
        for name, st, g in (
            ("S8", base, d10["grid_fallbacks_s8"]),
            ("S10", d10["stages"], d10["grid_fallbacks"]),
            *(
                (r, runs[s][r]["stages"], runs[s][r]["grid_fallbacks"])
                for r in ("S11", "S11a")
                if runs[s][r] and not runs[s][r].get("unfinished")
            ),
        ):
            rows.append(
                dict(
                    sample=NAMES[s],
                    run=name,
                    seed=st["seed_step"],
                    grid=st["grid_fallback"],
                    nb=st["neighbour"],
                    rv=st["revisit"],
                    pre=st["before_pass"],
                    ps=st["pass_"],
                    tot=st["total"],
                    f8=pct(1 - st["total"] / base["total"]),
                    f10=(None if name == "S8" else pct(1 - st["total"] / d10["stages"]["total"])),
                    gf=g,
                )
            )
        for r in ("S11", "S11a"):
            if runs[s][r] and runs[s][r].get("unfinished"):
                rows.append(dict(sample=NAMES[s], run=f"{r} (unfinished: 90 minute cap)"))
    cols = ["sample", "run", "seed", "grid", "nb", "rv", "pre", "ps", "tot", "f8", "f10", "gf"]
    T["stages"] = markdown_table(rows, cols, {"f8": ".1f", "f10": ".1f"}).replace(
        "| sample | run | seed | grid | nb | rv | pre | ps | tot | f8 | f10 | gf |",
        "| sample | run | seed step | grid fallback | neighbour fits | revisits | before the pass"
        " | boundary pass | total with the pass | fewer than S8 (%) | fewer than S10 (%)"
        " | grid fallbacks (count) |",
    )
    # ---- accuracy and rule
    rows = []
    for s in SAMPLES:
        for r in ("S11", "S11a"):
            d = runs[s][r]
            if not d or d.get("unfinished"):
                continue
            ru = d["rule"]
            v = d.get("vs_s10")
            rows.append(
                dict(
                    sample=NAMES[s],
                    run=r,
                    wb=d["before"]["wrong"],
                    wbi=ci(d["before"]),
                    wa=d["after"]["wrong"],
                    wai=ci(d["after"]),
                    w8b=d["s8_before"]["wrong"],
                    w8=d["s8_after"]["wrong"],
                    w10b=v["before_s10"]["wrong"] if v else None,
                    w10=v["after_s10"]["wrong"] if v else None,
                    fxb=d["paired_before"]["fixed"],
                    bkb=d["paired_before"]["broken"],
                    fx=d["paired_after"]["fixed"],
                    bk=d["paired_after"]["broken"],
                    bw=d["paired_after"]["both_wrong"],
                    fxb10=v["paired_before"]["fixed"] if v else None,
                    bkb10=v["paired_before"]["broken"] if v else None,
                    fx10=v["paired_after"]["fixed"] if v else None,
                    bk10=v["paired_after"]["broken"] if v else None,
                    net=ru["net_broken"],
                    med=d["after"]["median_err_deg"],
                    med8=d["s8_after"]["median_err_deg"],
                    dm=ru["med_diff"],
                    ra=str(ru["a"]),
                    rb=str(ru["b"]),
                    rc=str(ru["c"]),
                    ps=str(ru["passes"]),
                    fw=pct(ru["fewer"]),
                )
            )
    cols = ["sample", "run", "wb", "wbi", "wa", "wai", "w8b", "w8", "w10b", "w10"]
    T["accuracy"] = markdown_table(rows, cols).replace(
        "| sample | run | wb | wbi | wa | wai | w8b | w8 | w10b | w10 |",
        "| sample | run | wrong before the pass | 95% interval | wrong after the pass"
        " | 95% interval | S8 before the pass | S8 after the pass | S10 before the pass"
        " | S10 after the pass |",
    )
    cols = ["sample", "run", "fxb", "bkb", "fx", "bk", "bw", "fxb10", "bkb10", "fx10", "bk10"]
    T["paired"] = markdown_table(rows, cols).replace(
        "| sample | run | fxb | bkb | fx | bk | bw | fxb10 | bkb10 | fx10 | bk10 |",
        "| sample | run | fixed vs S8, before the pass | broken vs S8, before the pass"
        " | fixed vs S8, after the pass | broken vs S8, after the pass"
        " | wrong in both, after the pass | fixed vs S10, before the pass"
        " | broken vs S10, before the pass | fixed vs S10, after the pass"
        " | broken vs S10, after the pass |",
    )
    cols = ["sample", "run", "net", "med", "med8", "dm", "ra", "rb", "rc", "ps", "fw"]
    fmts = {"med": ".5f", "med8": ".5f", "dm": ".5f", "fw": ".1f"}
    T["rule"] = markdown_table(rows, cols, fmts).replace(
        "| sample | run | net | med | med8 | dm | ra | rb | rc | ps | fw |",
        "| sample | run | broken minus fixed vs S8 (after the pass) | median error, S11 (deg)"
        " | median error, S8 (deg) | difference (deg) | rule: at least 20% fewer evaluations"
        " | rule: broken minus fixed at most 5 | rule: median at most 0.005 degree above"
        " | passes | fewer evaluations than S8 (%) |",
    )
    # evaluations without the grid-fallback stage (it varies with the path; see the page)
    rows = []
    for s in SAMPLES:
        d10 = runs[s]["S10"]
        d11 = runs[s]["S11"]
        if not d11 or d11.get("unfinished"):
            continue
        b8 = d10["stages_s8"]["total"] - d10["stages_s8"]["grid_fallback"]
        b10 = d10["stages"]["total"] - d10["stages"]["grid_fallback"]
        b11 = d11["stages"]["total"] - d11["stages"]["grid_fallback"]
        rows.append(
            dict(
                sample=NAMES[s],
                s8=b8,
                s10=b10,
                s11=b11,
                f8=pct(1 - b11 / b8),
                f10=pct(1 - b11 / b10),
                gsave=d10["stages_s8"]["grid_fallback"] - d11["stages"]["grid_fallback"],
                tsave=d10["stages_s8"]["total"] - d11["stages"]["total"],
                gsave10=d10["stages"]["grid_fallback"] - d11["stages"]["grid_fallback"],
                tsave10=d10["stages"]["total"] - d11["stages"]["total"],
            )
        )
    T["nogrid"] = markdown_table(
        rows,
        ["sample", "s8", "s10", "s11", "f8", "f10", "gsave", "tsave", "gsave10", "tsave10"],
        {"f8": ".1f", "f10": ".1f"},
    ).replace(
        "| sample | s8 | s10 | s11 | f8 | f10 | gsave | tsave | gsave10 | tsave10 |",
        "| sample | S8 without grid fallbacks | S10 without grid fallbacks"
        " | S11 without grid fallbacks | S11 fewer than S8 (%) | S11 fewer than S10 (%)"
        " | grid-fallback evaluations S11 saves against S8"
        " | all evaluations S11 saves against S8"
        " | grid-fallback evaluations S11 saves against S10"
        " | all evaluations S11 saves against S10 |",
    )
    # ---- predictions
    rows = []
    for s in SAMPLES:
        d = runs[s]["S11"]
        if not d or d.get("unfinished"):
            continue
        v = d["vs_s10"]
        rows.append(
            dict(
                sample=NAMES[s],
                tot=d["stages"]["total"],
                lim=PRED2[s],
                gf=d["grid_fallbacks"],
                gf8=d["grid_fallbacks_s8"],
                sd=d["stages"]["seed_step"],
                p1=v["seed_step_evals_part1"],
                dev=pct(v["seed_step_run_over_part1"]),
                mc=pct(d["mc_share"]["share"]) if "mc_share" in d else None,
                mcn=d["mc_share"]["candidates"] if "mc_share" in d else None,
            )
        )
    T["predictions"] = markdown_table(
        rows,
        ["sample", "tot", "lim", "gf", "gf8", "sd", "p1", "dev", "mc", "mcn"],
        {"lim": ".0f", "dev": ".1f", "mc": ".1f"},
    ).replace(
        "| sample | tot | lim | gf | gf8 | sd | p1 | dev | mc | mcn |",
        "| sample | S11 total with the pass | prediction 2: below | S11 grid fallbacks"
        " | S8 grid fallbacks | S11 seed-step evaluations | Part 1 count, same sample"
        " | difference (%) | candidates where MC ran (%) | candidates refined |",
    )
    # ---- review fixes (D-20261111-N16): shares by stage, cost of a fallback, uncounted checks
    D = derived()
    rows = []
    for ref in ("S8", "S10"):
        for s in SAMPLES:
            if s not in D["stage_shares"]:
                continue
            sh = D["stage_shares"][s][ref]
            rows.append(
                dict(
                    sample=NAMES[s],
                    vs=ref,
                    save=sh["total_saving"],
                    **{k: pct(sh[k]) for k in STAGES},
                )
            )
    T["shares"] = markdown_table(
        rows, ["sample", "vs", "save", *STAGES], {k: ".1f" for k in STAGES}
    ).replace(
        "| sample | vs | save | seed_step | grid_fallback | neighbour | revisit | pass_ |",
        "| sample | against | evaluations S11 saves | seed step (%) | grid fallback (%)"
        " | neighbour fits (%) | revisits (%) | boundary pass (%) |",
    )
    rows = []
    for s in SAMPLES:
        if s not in D["grid_fallback_cost"]:
            continue
        fb = D["grid_fallback_cost"][s]
        row: Dict[str, Any] = dict(sample=NAMES[s])
        for r in ("S8", "S10", "S11"):
            row[f"g{r}"] = fb[r]["grid_evals"]
            row[f"n{r}"] = fb[r]["fallbacks"]
            pf = fb[r]["per_fallback"]
            row[f"p{r}"] = "n/a" if pf is None else f"{pf:,.0f}".replace(",", "")
        rows.append(row)
    gcols = [f"{x}{r}" for r in ("S8", "S10", "S11") for x in ("g", "n", "p")]
    T["fallback_cost"] = markdown_table(rows, ["sample", *gcols]).replace(
        "| sample | " + " | ".join(gcols) + " |",
        "| sample | S8 grid-fallback evaluations | S8 fallbacks | S8 per fallback"
        " | S10 grid-fallback evaluations | S10 fallbacks | S10 per fallback"
        " | S11 grid-fallback evaluations | S11 fallbacks | S11 per fallback |",
    )
    rows = [
        dict(sample=NAMES[s], n=u["mixed_candidates"], tot=u["total_evals"], sh=pct(u["share"]))
        for s, u in D["uncounted_hit_checks"].items()
    ]
    T["uncounted"] = markdown_table(rows, ["sample", "n", "tot", "sh"], {"sh": ".4f"}).replace(
        "| sample | n | tot | sh |",
        "| sample | uncounted hit-ratio checks (one per mixed candidate) | counted evaluations"
        " of S11 with the pass | uncounted as a share of counted (%) |",
    )
    wi = D["plan_whatif"]
    if wi:
        rows = []
        for s in SAMPLES:
            w = wi["by_sample"][s]
            for x in w["steps_differ"]:
                rows.append(
                    dict(
                        sample=NAMES[s],
                        step=x["step"],
                        cand=x["candidate"] + 1,
                        ch=x["cma_hit"],
                        cc=x["cma_cost"],
                        ce=x["cma_evals"],
                        mh=x["mc_hit"],
                        mc=x["mc_cost"],
                        me=x["mc_evals"],
                        ne=x["next_cma_evals"],
                        nh=x["next_cma_hit"],
                        er=w["evals_rule"],
                        ev=w["evals_variant"],
                    )
                )
        fm = {"ch": ".5f", "cc": ".5f", "mh": ".5f", "mc": ".5f", "nh": ".3f"}
        cols = [
            "sample",
            "step",
            "cand",
            "ch",
            "cc",
            "ce",
            "mh",
            "mc",
            "me",
            "ne",
            "nh",
            "er",
            "ev",
        ]
        T["whatif"] = markdown_table(rows, cols, fm).replace(
            "| " + " | ".join(cols) + " |",
            "| sample | replay step | candidate | CMA end hit ratio | CMA end cost"
            " | CMA evaluations | MC end hit ratio | MC end cost | MC evaluations"
            " | next candidate: CMA evaluations | next candidate: CMA end hit ratio"
            " | X200 seed-step evaluations, this page | the plan's what-if |",
        )
    return T


def plots(T: Dict[str, str]) -> None:
    cols = dict(
        seed_step="#4c72b0", grid_fallback="#c44e52", neighbour="#55a868", revisit="#dd8452",
        pass_="#8172b3",
    )  # fmt: skip
    fig, axes = plt.subplots(1, 3, figsize=(13, 4.2))
    for ax, s in zip(axes, SAMPLES):
        d10 = load(FIT / f"part2_S10_{s}.json")
        labels, data = ["S8", "S10"], [d10["stages_s8"], d10["stages"]]
        for r in ("S11", "S11a"):
            d = load(MIX / f"part2_{r}_{s}.json")
            if d and not d.get("unfinished"):
                labels.append(r)
                data.append(d["stages"])
        bottom = np.zeros(len(data))
        for k in STAGES:
            v = np.array([x[k] for x in data], float) / 1e6
            ax.bar(labels, v, bottom=bottom, color=cols[k], label=STAGE_LABEL[k])
            bottom += v
        ax.set_title(NAMES[s])
        ax.set_ylabel("evaluations (millions)")
    axes[0].legend(fontsize=8)
    fig.suptitle("Evaluations by stage, with the boundary pass")
    fig.tight_layout()
    fig.savefig(MIX / "ar5_mix_stages.png", dpi=120)
    plt.close(fig)
    # plot 2: Part 1 arms (evaluations against the share of steps within 1 degree)
    p1 = load(MIX / "part1.json")
    fig, axes = plt.subplots(1, 3, figsize=(13, 4.0))
    for ax, s in zip(axes, SAMPLES):
        d = p1["by_sample"][s]
        n = d["n_steps"]
        for a in ARMS:
            c = d["arms"][a]
            ax.scatter(c["evals"] / 1e6, 100 * c["n_right"] / n, s=60, label=a)
            ax.scatter(c["evals"] / 1e6, 100 * c["n_pass"] / n, s=60, marker="x", color="k")
        ax.set_title(NAMES[s])
        ax.set_xlabel("seed-step evaluations (millions)")
        ax.set_ylabel("share of steps (%)")
    axes[0].legend(fontsize=8, title="dot: within 1 degree; cross: seed test")
    fig.tight_layout()
    fig.savefig(MIX / "ar5_mix_arms.png", dpi=120)
    plt.close(fig)


def main() -> None:
    T = tables()
    write_tables(MIX / "tables.md", T)
    plots(T)
    summary: Dict[str, Any] = dict(
        pick=load(MIX / "pick.json")["pick"],
        part2={r: {s: load(MIX / f"part2_{r}_{s}.json") for s in SAMPLES} for r in ("S11", "S11a")},
        **derived(),
    )
    (MIX / "summary.json").write_text(json.dumps(summary, indent=1) + "\n")
    print("wrote tables.md, summary.json and two plots")


if __name__ == "__main__":
    main()
