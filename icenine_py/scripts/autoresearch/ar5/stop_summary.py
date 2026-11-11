"""AR5 step 1: summary of all parts -> ``summary.json``, ``tables.md``, the plots.

  uv run python scripts/autoresearch/ar5/stop_summary.py

Reads the result files of ``benchmarks/autoresearch/ar5/`` (parts 1 to 3); writes ``summary.json``
(every number the findings page quotes), ``tables.md`` (marker blocks for
``scripts/dev/sync_doc_tables.py``) and three plots. Part 4 (whole-layer run T1) was not run because
no rule passed; the page says so.
"""

import json
import sys
from pathlib import Path
from typing import Any, Dict, List

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

HERE = Path(__file__).resolve().parent
AR5 = HERE.parents[2] / "benchmarks" / "autoresearch" / "ar5"
sys.path.insert(0, str(HERE.parents[1] / "common"))
sys.path.insert(0, str(HERE))
import doc_tables as DT  # noqa: E402
import stop_rules as SR  # noqa: E402

SAMPLES = ("big", "grad", "mid")
NAME = {"big": "big", "grad": "first gradient sample", "mid": "mid"}


def J(name: str) -> Any:
    return json.loads((AR5 / name).read_text())


def main() -> None:
    p1, p1p = J("part1_results.json"), J("part1_pick.json")
    p2, cases = J("part2_summary.json"), J("part2_cases.json")
    p3f, p3p = J("part3_flags.json"), J("part3_pick.json")
    vt, v1 = J("part3_votes.json"), J("part3_v1.json")["summary"]
    tables: Dict[str, str] = {}
    summary: Dict[str, Any] = {}

    # ---- Part 1
    rows = []
    for name, per in p1["table"].items():
        for s in SAMPLES:
            r = per[s]
            full_gain = p1["table"]["full"][s]["net_gain"]
            rows.append(
                dict(
                    setting=name,
                    sample=NAME[s],
                    kept=r["tries_kept"],
                    of=r["tries"],
                    evals=r["evals"],
                    wrong=r["wrong_after"],
                    fixed=r["fixed"],
                    broken=r["broken"],
                    gain=r["net_gain"],
                    lost=100.0 * (1 - r["net_gain"] / full_gain),
                )
            )
    tables["p1_rules"] = DT.markdown_table(
        rows,
        ["setting", "sample", "kept", "of", "evals", "wrong", "fixed", "broken", "gain", "lost"],
        {"kept": ",d", "of": ",d", "evals": ",d", "lost": ".1f"},
        labels=None,
    )
    pooled_rows = []
    full_ev = p1p["full_pass"]["evals"]
    full_gain = p1p["full_pass"]["net_gain"]
    for r in [dict(name="full pass", evals=full_ev, net_gain=full_gain, tries_kept=None)] + p1p[
        "rows"
    ]:
        pooled_rows.append(
            dict(
                setting=r["name"],
                evals=r["evals"],
                cut=100.0 * (1 - r["evals"] / full_ev),
                gain=r["net_gain"],
                lost=100.0 * (1 - r["net_gain"] / full_gain),
                ok="yes" if r["net_gain"] >= 0.95 * full_gain - SR.TOL else "no",
            )
        )
    tables["p1_pick"] = DT.markdown_table(
        pooled_rows,
        ["setting", "evals", "cut", "gain", "lost", "ok"],
        {"evals": ",d", "cut": ".1f", "lost": ".1f"},
    )
    pr = p1["pass_rule_mid"]
    summary["part1"] = dict(
        pick=p1p["pick"],
        pooled_full_evals=full_ev,
        pooled_full_net_gain=full_gain,
        min_net_gain=p1p["min_net_gain"],
        p1_pick_d=p1p["p1_pick_d"],
        p1_has_qualifying_setting=p1p["p1_has_qualifying_setting"],
        pass_rule_mid=pr,
        prediction1=p1["prediction1"]["holds"],
        prediction1_pick_samples=p1["prediction1"]["holds_on_pick_samples"],
        max_cut_with_loss_within_5pc=max(
            [
                100 * (1 - r["evals"] / full_ev)
                for r in p1p["rows"]
                if r["net_gain"] >= 0.95 * full_gain - SR.TOL
            ]
            or [0.0]
        ),
        passes=pr["passes"],
    )

    # ---- Part 2
    fail_rows, rung_rows = [], []
    for s in ("grad", "mid"):
        for c in cases[s]:
            fail_rows.append(
                dict(
                    sample=NAME[s],
                    voxel=c["voxel"],
                    kind="reseed" if c["reseed"] else "first seed",
                    cand=c["n_cand"],
                    rank=c["rank_near_truth"] + 1 if c["rank_near_truth"] >= 0 else None,
                    c20="yes" if c["in_costed20"] else "no",
                    r5="yes" if c["in_refined5"] else "no",
                    s8best=max(c["stored_fin_hit"]),
                    gridhit=c["grid"]["hit_map"],
                    griderr=c["grid"]["err"],
                )
            )
            rung_rows.append(
                dict(
                    sample=NAME[s],
                    voxel=c["voxel"],
                    r1hit=c["R1"]["hit"],
                    r1ev=c["R1"]["evals"],
                    r1err=c["R1"]["err"],
                    r3hit=c["R3"]["hit"],
                    r3ev=c["R3"]["evals"],
                    r3err=c["R3"]["err"],
                    ends=c["ladder"]["rung"],
                    ladev=int(c["ladder"]["evals"]),
                    gridev=c["grid_evals"],
                )
            )
    tables["p2_failed"] = DT.markdown_table(
        fail_rows,
        ["sample", "voxel", "kind", "cand", "rank", "c20", "r5", "s8best", "gridhit", "griderr"],
        {"s8best": ".3f", "gridhit": ".3f", "griderr": ".2f"},
    )
    tables["p2_rungs"] = DT.markdown_table(
        rung_rows,
        [
            "sample",
            "voxel",
            "r1hit",
            "r1ev",
            "r1err",
            "r3hit",
            "r3ev",
            "r3err",
            "ends",
            "ladev",
            "gridev",
        ],
        {
            "r1hit": ".3f",
            "r3hit": ".3f",
            "r1err": ".2f",
            "r3err": ".2f",
            "r1ev": ",d",
            "r3ev": ",d",
            "ladev": ",d",
            "gridev": ",d",
        },
    )
    tot_rows = []
    for k, label in (("grad", NAME["grad"]), ("mid", "mid"), ("pooled", "both")):
        d = p2["pooled"] if k == "pooled" else p2["per_sample"][k]
        tot_rows.append(
            dict(
                group=label,
                n=d["n"],
                r1pass=d["r1_pass_seed_test"],
                r1right=d["r1_right"],
                r3pass=d["r3_pass_seed_test"],
                r3right=d["r3_right"],
                ends=f"{d['ladder_rung']['R1']} / {d['ladder_rung']['R3']} / "
                f"{d['ladder_rung']['grid']}",
                ladev=int(d["ladder_evals"]),
                gridev=int(d["grid_evals"]),
                share=100.0 * d["ladder_share_of_grid"],
                rl=d["right_ladder"],
                rg=d["right_grid"],
            )
        )
    tables["p2_totals"] = DT.markdown_table(
        tot_rows,
        [
            "group",
            "n",
            "r1pass",
            "r1right",
            "r3pass",
            "r3right",
            "ends",
            "ladev",
            "gridev",
            "share",
            "rl",
            "rg",
        ],
        {"ladev": ",d", "gridev": ",d", "share": ".1f"},
    )
    summary["part2"] = dict(
        pooled=p2["pooled"],
        per_sample=p2["per_sample"],
        prediction2=p2["prediction2"],
        prediction3=p2["prediction3"],
    )

    # ---- Part 3
    thr = p3p["threshold"]
    flag_rows = []
    for s in SAMPLES:
        for k, row in p3f["tables"][s].items():
            a = row["all"]
            flag_rows.append(
                dict(
                    sample=NAME[s],
                    thr=float(k),
                    flagged=a["flagged"],
                    share=100.0 * a["flagged"] / a["n"],
                    wf=a["wrong_flagged"],
                    wu=a["wrong_unflagged"],
                    nw=a["n_wrong"],
                )
            )
    tables["p3_flag"] = DT.markdown_table(
        flag_rows,
        ["sample", "thr", "flagged", "share", "wf", "wu", "nw"],
        {"thr": ".2f", "flagged": ",d", "share": ".1f"},
    )
    src_rows = []
    key = f"{thr:.2f}"
    for s in SAMPLES:
        for g in ("seed", "neighbor", "revisit", "reseed", "boundary", "all"):
            a = p3f["tables"][s][key][g]
            src_rows.append(
                dict(
                    sample=NAME[s],
                    group=g,
                    n=a["n"],
                    flagged=a["flagged"],
                    wrong=a["n_wrong"],
                    wf=a["wrong_flagged"],
                )
            )
    tables["p3_source"] = DT.markdown_table(
        src_rows, ["sample", "group", "n", "flagged", "wrong", "wf"], {"n": ",d", "flagged": ",d"}
    )
    left_rows = []
    for s in SAMPLES:
        a = p3f["at_threshold"][s]
        lo, hi = a["wrong_unflagged_interval"][2:]
        left_rows.append(
            dict(
                sample=NAME[s],
                n=a["n"],
                flagged=a["flagged"],
                share=100 * a["share_flagged"],
                wrong=a["n_wrong"],
                caught=a["wrong_flagged"],
                left=a["wrong_unflagged"],
                unflagged=a["unflagged"],
                lo=100 * lo,
                hi=100 * hi,
            )
        )
    tables["p3_left"] = DT.markdown_table(
        left_rows,
        ["sample", "n", "flagged", "share", "wrong", "caught", "left", "unflagged", "lo", "hi"],
        {"n": ",d", "flagged": ",d", "unflagged": ",d", "share": ".1f", "lo": ".3f", "hi": ".3f"},
    )
    vrows = []
    for g in (*SAMPLES, "all"):
        d = vt["groups"][g]
        c = d["best_cut"]
        vrows.append(
            dict(
                group=NAME.get(g, "all three"),
                nw=d["n_wrong_flagged"],
                nr=d["n_right_flagged"],
                mw=d["median_wrong"],
                mr=d["median_right"],
                cut=c["cut"],
                w=f"{c['wrong_ge']} of {c['n_wrong']}",
                r=f"{c['right_ge']} of {c['n_right']}",
                sep="yes" if d["separates"] else "no",
            )
        )
    tables["p3_votes"] = DT.markdown_table(
        vrows,
        ["group", "nw", "nr", "mw", "mr", "cut", "w", "r", "sep"],
        {"mw": ".0f", "mr": ".0f", "cut": ".0f"},
    )
    v1rows = [
        dict(
            thr=v1["threshold"],
            flagged=v1["flagged_total"],
            done=v1["voxels_processed"],
            ev=v1["evals_processed"],
            mean=v1["mean_evals_per_voxel"],
            est=v1["evals_whole_flag_estimate"],
            share=100 * v1["share_of_run_estimate"],
            kept=v1["kept"],
            wrong=v1["wrong_among_processed"],
            fixed=v1["fixed"],
            broken=v1["broken"],
            after=v1["wrong_after_v1"],
        )
    ]
    tables["p3_v1"] = DT.markdown_table(
        v1rows,
        [
            "thr",
            "flagged",
            "done",
            "ev",
            "mean",
            "share",
            "kept",
            "wrong",
            "fixed",
            "broken",
            "after",
        ],
        {
            "thr": ".2f",
            "flagged": ",d",
            "done": ",d",
            "ev": ",d",
            "mean": ",.0f",
            "share": ".1f",
            "kept": ",d",
        },
    )
    summary["part3"] = dict(
        threshold=thr,
        flags=p3f["at_threshold"],
        pooled_pick=p3f["pooled_pick"],
        prediction4=p3f["prediction4"],
        votes=vt["groups"],
        v1=v1,
        prediction5_holds=bool(not vt["groups"]["all"]["separates"]),
    )

    # ---- predictions and rules
    p1s, p2s, p3s = summary["part1"], summary["part2"], summary["part3"]
    pred = [
        (
            "1",
            "a try-test setting cuts 70% of the pass and loses at most 5% of its gain",
            "does not agree",
            f"best cut within the 5% loss: {p1s['max_cut_with_loss_within_5pc']:.1f}%",
        ),
        (
            "2",
            "a candidate near the truth is in the first 100 and not the 20 costed, in half the "
            "cases",
            "does not agree",
            f"{p2s['prediction2']['cases_with_candidate_in_first100_not_costed20']} of "
            f"{p2s['prediction2']['n']}",
        ),
        (
            "3",
            "the ladder passes its rule",
            "does not agree",
            f"{100 * p2s['pooled']['ladder_share_of_grid']:.1f}% of the grid evaluations",
        ),
        (
            "4",
            "a hit-ratio threshold flags 80% of wrong voxels and at most 5% of voxels",
            "agrees" if p3s["prediction4"]["holds"] else "does not agree",
            f"first part holds ({100 * p3s['pooled_pick']['recall']:.0f}% of wrong voxels caught); "
            f"second part fails ({100 * p3s['pooled_pick']['share_flagged']:.1f}% of voxels "
            f"flagged on the pick samples, {100 * p3s['flags']['mid']['share_flagged']:.1f}% on "
            "mid)",
        ),
        (
            "5",
            "the vote count does not separate wrong from right among flagged voxels",
            "agrees" if p3s["prediction5_holds"] else "does not agree",
            "a weak separation exists; "
            f"cut {vt['groups']['all']['best_cut']['cut']:.0f}: "
            f"{vt['groups']['all']['best_cut']['wrong_ge']} of "
            f"{vt['groups']['all']['best_cut']['n_wrong']} wrong, "
            f"{vt['groups']['all']['best_cut']['right_ge']} of "
            f"{vt['groups']['all']['best_cut']['n_right']} right at or above",
        ),
        (
            "6",
            "T1 needs under 14 M (gradient) and 22 M (mid) evaluations",
            "not evaluated",
            "no rule passed, so T1 was not run",
        ),
    ]
    tables["predictions"] = DT.markdown_table(
        [dict(n=a, what=b, result=c, numbers=d) for a, b, c, d in pred],
        ["n", "what", "result", "numbers"],
    )
    rules = [
        (
            "try test (Part 1)",
            "on mid: at least half of the pass cut, at most 5% of the gain lost",
            "not met",
            f"{100 * pr['cut']:.1f}% cut, "
            f"{100 * (1 - pr['net_gain'] / pr['net_gain_full']):.1f}% of the gain lost "
            f"({pr['net_gain']} of {pr['net_gain_full']})",
        ),
        (
            "ladder (Part 2)",
            "at most 2/3 of the grid evaluations, at most one case fewer right",
            "not met",
            f"{100 * p2s['pooled']['ladder_share_of_grid']:.1f}% of the evaluations; right "
            f"{p2s['pooled']['right_ladder']} against {p2s['pooled']['right_grid']}",
        ),
        (
            "V1 (Part 3)",
            "fixes more than it breaks, fixes a quarter, costs under 5% of the run",
            "met" if v1["passes"] else "not met",
            f"fixed {v1['fixed']}, broken {v1['broken']}, {v1['wrong_among_processed']} wrong "
            f"processed; about {100 * v1['share_of_run_estimate']:.0f}% of the run",
        ),
    ]
    tables["rules"] = DT.markdown_table(
        [dict(rule=a, what=b, verdict=c, numbers=d) for a, b, c, d in rules],
        ["rule", "what", "verdict", "numbers"],
    )
    summary["rules"] = [dict(rule=a, what=b, verdict=c, numbers=d) for a, b, c, d in rules]
    summary["predictions"] = [dict(n=a, what=b, result=c, numbers=d) for a, b, c, d in pred]
    (AR5 / "summary.json").write_text(json.dumps(summary, indent=1, default=float) + "\n")
    DT.write_tables(AR5 / "tables.md", tables)

    # ---- plots
    plot_try_test(p1)
    plot_flag(p3f)
    plot_votes()


def plot_try_test(p1: Dict[str, Any]) -> None:
    fig, axes = plt.subplots(1, 3, figsize=(11, 3.6), sharey=True)
    for ax, s in zip(axes, SAMPLES):
        full = p1["table"]["full"][s]
        for name, per in p1["table"].items():
            if name == "full":
                continue
            r = per[s]
            x = 100.0 * r["tries_kept"] / r["tries"]
            y = 100.0 * (1 - r["net_gain"] / full["net_gain"])
            ax.scatter(
                x,
                y,
                color="C0" if name.startswith("P1") else "C1" if name.startswith("P2") else "C2",
            )
            ax.annotate(
                name.replace("P3 d=0.2 ", "P3 "),
                (x, y),
                fontsize=6,
                xytext=(3, 3),
                textcoords="offset points",
            )
        ax.axhline(5, color="k", ls="--", lw=0.8)
        ax.set_title(NAME[s], fontsize=9)
        ax.set_xlabel("tries kept (% of the pass)")
    axes[0].set_ylabel("net gain lost (% of the full pass)")
    fig.suptitle("Try test: tries kept against net gain lost (dashed: 5% loss)", fontsize=10)
    fig.tight_layout()
    fig.savefig(AR5 / "ar5_try_test.png", dpi=110)
    plt.close(fig)


def plot_flag(p3f: Dict[str, Any]) -> None:
    fig, ax = plt.subplots(figsize=(5.5, 3.6))
    for s in SAMPLES:
        xs, ys = [], []
        for k, row in p3f["tables"][s].items():
            a = row["all"]
            xs.append(100.0 * a["flagged"] / a["n"])
            ys.append(100.0 * a["wrong_flagged"] / a["n_wrong"])
        ax.plot(xs, ys, marker="o", label=NAME[s])
    ax.axvline(5, color="k", ls="--", lw=0.8)
    ax.axhline(80, color="k", ls=":", lw=0.8)
    ax.set_xlabel("voxels flagged (% of all voxels), thresholds 0.50 to 0.95")
    ax.set_ylabel("wrong voxels flagged (%)")
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(AR5 / "ar5_flag.png", dpi=110)
    plt.close(fig)


def plot_votes() -> None:
    z = np.load(AR5 / "part3_votes.npz")
    wv, rv = [], []
    for s in SAMPLES:
        fl, wr, v = z[f"{s}_flagged"], z[f"{s}_wrong"], z[f"{s}_votes"]
        wv.append(v[fl & wr])
        rv.append(v[fl & ~wr])
    wv, rv = np.concatenate(wv), np.concatenate(rv)
    fig, ax = plt.subplots(figsize=(6, 3.2))
    for y, d, c in ((1, wv, "C3"), (0, rv, "C0")):
        jit = ((np.arange(len(d)) * 0.6180339887) % 1.0 - 0.5) * 0.5  # fixed spread, no random
        ax.scatter(d, y + jit, s=3, alpha=0.4, color=c)
    ax.set_yticks([0, 1])
    ax.set_yticklabels([f"right (n={len(rv)})", f"wrong (n={len(wv)})"])
    ax.set_xlabel("highest vote count, unexplained rays (flagged voxels)")
    fig.tight_layout()
    fig.savefig(AR5 / "ar5_votes.png", dpi=110)
    plt.close(fig)


if __name__ == "__main__":
    main()
