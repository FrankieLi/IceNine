"""AR5 step 12, phase A step 6: validity, plot, tables, pareto rows, predictions, verdict.

  uv run python scripts/autoresearch/ar5/stallhit_summary.py

Reads ``stallhit/stats_<sample>.json`` and ``stallhit/unmatched_<sample>.json``
(``stallhit_stats.py``) and ``stallhit/gate.json``. Writes ``ar5_stallhit_evals_wrong.png`` (made
first: run evaluations against wrong voxels, S8, G and H, before and after the pass),
``tables.md`` (blocks for ``sync_doc_tables.py``), ``summary.json`` (every number of the tables, the
predictions, the rules, the verdict) and ``pareto_rows.json`` (the format of
``ar5/patience/pareto_rows.json``: arm H on the five samples, arms S8 and G on ``held2``, with a
status field). No p-values; paired counts are by voxel.
"""

import json
import sys
from pathlib import Path
from typing import Any, Dict, List

import numpy as np

HERE = Path(__file__).resolve().parent
SH = HERE.parents[2] / "benchmarks" / "autoresearch" / "ar5" / "stallhit"
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent.parent / "common"))

import doc_tables as DT  # noqa: E402
import stallhit_judge as SJ  # noqa: E402
import stats as ST  # noqa: E402

SAMPLES = SJ.SAMPLES
N = 24570
STATUS_PENDING = "study result of step 12 (rule PH(3000)); no key yet; not a candidate arm"


def f(x: float, d: int = 5) -> str:
    return f"{x:.{d}f}"


def n(x: int) -> str:
    return f"{int(x):,}"


def plot(st: Dict[str, Any]) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.ticker import LogFormatter

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.8))
    col = dict(
        held2="tab:purple", held="tab:red", mid="tab:blue", grad="tab:green", big="tab:orange"
    )
    for ax, tag, title in zip(axes, ("before", "after"), ("before the pass", "after the pass")):
        for s in SAMPLES:
            a = st[s]["arms"]
            key = "before_pass" if tag == "before" else "total"
            pts = [
                ("S8", "^", a["H"]["stages_s8"][key], a["H"][f"s8_{tag}"]),
                ("G", "o", a["G"]["stages"][key], a["G"][tag]),
                ("H", "s", a["H"]["stages"][key], a["H"][tag]),
            ]
            for arm, mk, x, side in pts:
                k = side["wrong"]
                lo, hi = N * side["wilson_lo"], N * side["wilson_hi"]
                ax.errorbar(
                    x / 1e6,
                    k,
                    yerr=[[k - lo], [hi - k]],
                    fmt=mk,
                    color=col[s],
                    mfc="none" if arm == "S8" else col[s],
                    ms=6,
                    capsize=2,
                    lw=0.8,
                    label=f"{s} {arm}" if tag == "before" else None,
                )
        ax.set_xscale("log")
        ax.xaxis.set_minor_formatter(LogFormatter(labelOnlyBase=False, minor_thresholds=(2, 0.4)))
        ax.tick_params(axis="x", which="minor", labelsize=7)
        ax.set_xlabel(
            "run evaluations (millions, log%s)" % (", with the pass" if tag == "after" else "")
        )
        ax.set_ylabel("wrong voxels (more than 1 degree), Wilson 95%")
        ax.set_title(title)
        ax.grid(alpha=0.3)
    fig.legend(loc="lower center", fontsize=7, ncol=8, frameon=False)
    fig.suptitle(
        "S8 (triangle), arm G = step 9 (circle), arm H = G with stall check PH(3000) (square)",
        fontsize=10,
    )
    fig.tight_layout(rect=(0, 0.14, 1, 1))
    fig.savefig(SH / "ar5_stallhit_evals_wrong.png", dpi=110)


def pareto(st: Dict[str, Any], status: str) -> None:
    rows: List[Dict[str, Any]] = []

    def add(s: str, arm: str, order: str, d: Dict[str, Any], side_key: str, after: bool) -> None:
        side = d[side_key]
        total = d["stages_key"]["total" if after else "before_pass"]
        rows.append(
            dict(
                sample=s,
                arm=arm,
                order=order,
                boundary_pass=after,
                evals_total=total,
                evals_per_voxel=round(total / N, 2),
                wrong=side["wrong"],
                n=N,
                interval_pct=[round(100 * side["wilson_lo"], 2), round(100 * side["wilson_hi"], 2)],
                source="icenine_py/benchmarks/autoresearch/ar5/stallhit/stats_%s.json" % s,
                status=status,
            )
        )

    for s in SAMPLES:
        h = dict(st[s]["arms"]["H"])
        h["stages_key"] = h["stages"]
        for after in (False, True):
            add(
                s,
                "H",
                "interleaved + rescue G1 + stall check PH(3000)",
                h,
                "after" if after else "before",
                after,
            )
        if s == "held2":
            g = dict(st[s]["arms"]["G"])
            g["stages_key"] = g["stages"]
            s8 = dict(st[s]["arms"]["H"])
            s8["stages_key"] = s8["stages_s8"]
            for after in (False, True):
                add(s, "G", "interleaved + rescue G1", g, "after" if after else "before", after)
                add(
                    s, "S8", "S8 (stock seed step)", s8, "s8_after" if after else "s8_before", after
                )
    (SH / "pareto_rows.json").write_text(json.dumps(rows, indent=1) + "\n")


def pct(x: float) -> float:
    return float(x)


def predictions(st: Dict[str, Any], gate: Dict[str, Any]) -> Dict[str, Any]:
    out: Dict[str, Any] = {}
    ps = gate["per_sample"]
    p1 = {}
    for s in ("mid", "grad", "big"):
        d = ps[s]["prediction1"]
        none_changed = ps[s]["taken_changed"] == 0
        met = bool(
            d["g1"]
            and d["changed_in_range"]
            and d["all_changed_still_095"]
            and d["none_wrong"]
            and d["failed_seed_tests_equal"]
            and (d["mis_in_0_005_to_0_15"] or none_changed)
        )
        p1[s] = dict(met=met, changed=ps[s]["taken_changed"], detail=d)
    out["p1"] = p1
    out["p2"] = {
        s: dict(met=ps[s]["prediction2"]["saving_in_range"], saving=ps[s]["prediction2"]["saving"])
        for s in ("mid", "grad", "big")
    }
    rng = dict(mid=(6, 9), big=(3.5, 5.5), grad=(-0.5, 0.5), held=(2, 6), held2=(2, 8))
    p3 = {}
    for s, (lo, hi) in rng.items():
        r = st[s]["rules"]
        p3[s] = dict(total_cut_pct=r["total_cut_pct"], met=bool(lo <= r["total_cut_pct"] <= hi))
    p3["held2"]["seed_cut_pct"] = st["held2"]["rules"]["seed_cut_pct"]
    p3["held2"]["seed_met"] = bool(6 <= st["held2"]["rules"]["seed_cut_pct"] <= 18)
    p3["held2"]["met"] = bool(p3["held2"]["met"] and p3["held2"]["seed_met"])
    out["p3"] = p3
    p4 = {}
    for s in SAMPLES:
        a = st[s]["arms"]
        pv = st[s]["h_vs_g"]["voxels"]["after"]
        net = pv["broken"] - pv["fixed"]
        dm = a["H"]["after"]["median_err_deg"] - a["G"]["after"]["median_err_deg"]
        p4[s] = dict(net=net, median_shift_deg=dm, met=bool(-3 <= net <= 3 and abs(dm) <= 0.0005))
    out["p4"] = p4
    out["p5"] = dict(
        met=bool(st["mid"]["rules"]["r5"]),
        house_h=st["mid"]["house_h"],
    )
    return out


def house_cell(s: Dict[str, Any]) -> str:
    """Rule (5) numbers of arm H against S8: percent fewer; broken minus fixed; median shift."""
    fewer = SJ.pct_fewer(s["stages"]["total"], s["stages_s8"]["total"])
    net = s["paired_after"]["broken"] - s["paired_after"]["fixed"]
    shift = s["after"]["median_err_deg"] - s["s8_after"]["median_err_deg"]
    return f"{fewer:.1f}; {net}; {shift:.5f}"


def wilson_cell(k: int, total: int) -> str:
    lo, hi = ST.wilson(k, total)
    return f"{n(k)} of {n(total)} ({100 * lo:.1f}% to {100 * hi:.1f}%)"


def tables(st: Dict[str, Any], pred: Dict[str, Any]) -> Dict[str, str]:
    t: Dict[str, str] = {}
    rows = []
    for s in SAMPLES:
        for arm in ("G", "H"):
            sp = st[s]["arms"][arm]["stages_split"]
            rows.append(
                {
                    "sample": s,
                    "arm": arm,
                    "seed step": n(sp["seed_step"]),
                    "rescue": n(sp["rescue"]),
                    "grid": n(sp["grid_fallback"]),
                    "neighbour": n(sp["neighbour"]),
                    "revisit": n(sp["revisit"]),
                    "pass": n(sp["pass_"]),
                    "total": n(sp["total"]),
                }
            )
    t["stallhit_stages"] = DT.markdown_table(
        rows,
        ["sample", "arm", "seed step", "rescue", "grid", "neighbour", "revisit", "pass", "total"],
    )
    rows = []
    for s in SAMPLES:
        d = st[s]
        c = d["cuts"]
        rows.append(
            {
                "sample": s,
                "MC refinements": n(c["mc_refinements"]),
                "stall checks": n(c["stall_checks"]),
                "stopped (cuts)": n(c["stall_checks_stopped"]),
                "went on": n(c["stall_checks_went_on"]),
                "hit ratio of the checks that went on (min / median / max)": (
                    "%.3f / %.3f / %.3f"
                    % (c["went_on_hit_min"], c["went_on_hit_median"], c["went_on_hit_max"])
                    if c["stall_checks_went_on"]
                    else "n/a"
                ),
                "lowest hit ratio at a stop": f(c["stop_hit_min"], 4),
                "evaluations in cut ones": n(c["evals_in_cut"]),
            }
        )
    t["stallhit_cuts"] = DT.markdown_table(rows, list(rows[0].keys()))
    rows = []
    for s in SAMPLES:
        d = st[s]
        ch = d["taken_changed"]
        wg = d["wrong_seeds"]
        rows.append(
            {
                "sample": s,
                "seed steps G / H": f"{n(d['validity']['seed_steps_g'])} / "
                f"{n(d['validity']['seed_steps_h'])}",
                "matched": n(d["validity"]["matched_steps"]),
                "only in H": n(d["unmatched"]["only_h_n"]),
                "only in G": n(d["unmatched"]["only_g_n"]),
                "taken ends changed (matched)": n(d["taken_changed_n"]),
                "median angle to G (deg)": f(d["taken_changed_median_mis_deg"] or 0.0, 4),
                "changed, now at 0.95 or more": n(sum(1 for c in ch if c["hit_h"] >= 0.95)),
                "changed, max error H (deg)": f(max([c["err_h"] for c in ch] or [0.0]), 4),
                "H wrong where G right: matched + only-H": f"{wg['h_wrong_where_g_right']} + "
                f"{wg['unmatched_h_wrong_where_g_map_right']}",
            }
        )
    t["stallhit_changed"] = DT.markdown_table(rows, list(rows[0].keys()))
    rows = []
    for s in SAMPLES:
        d = st[s]
        pg = d["h_vs_g"]["voxels"]["after"]
        for arm in ("G", "H"):
            sd = d["arms"][arm]
            side = sd["after"]
            p8 = sd["paired_after"]
            rows.append(
                {
                    "sample": s,
                    "arm": arm,
                    "evaluations with the pass": n(sd["stages"]["total"]),
                    "seed step": n(sd["stages"]["seed_step"]),
                    "wrong after pass": n(side["wrong"]),
                    "Wilson 95% (%)": "%.2f to %.2f"
                    % (100 * side["wilson_lo"], 100 * side["wilson_hi"]),
                    "median error (deg)": f(side["median_err_deg"], 5),
                    "fixed / broken vs G": (
                        "n/a" if arm == "G" else f"{pg['fixed']} / {pg['broken']}"
                    ),
                    "fixed / broken vs S8": f"{p8['fixed']} / {p8['broken']}",
                }
            )
    t["stallhit_layer"] = DT.markdown_table(rows, list(rows[0].keys()))
    rows = []
    for s in SAMPLES:
        d = st[s]
        for arm in ("G", "H"):
            r = d["rescue"][arm]
            rows.append(
                {
                    "sample": s,
                    "arm": arm,
                    "failed seed tests": n(r["fallback_calls"]),
                    "rescued": n(r["rescued"]),
                    "grid searches": n(r["grid_fallbacks_left"]),
                }
            )
    t["stallhit_fallbacks"] = DT.markdown_table(rows, list(rows[0].keys()))
    rows = []
    for s in SAMPLES:
        d = st[s]
        a = d["arms"]["H"]
        pg = d["h_vs_g"]["voxels"]
        for tag in ("before", "after"):
            rows.append(
                {
                    "sample": s,
                    "pass": tag,
                    "fixed vs G": pg[tag]["fixed"],
                    "broken vs G": pg[tag]["broken"],
                    "both wrong vs G": pg[tag]["both_wrong"],
                    "fixed vs S8": a[f"paired_{tag}"]["fixed"],
                    "broken vs S8": a[f"paired_{tag}"]["broken"],
                    "both wrong vs S8": a[f"paired_{tag}"]["both_wrong"],
                }
            )
    t["stallhit_paired"] = DT.markdown_table(rows, list(rows[0].keys()))
    rows = []
    for s in SAMPLES:
        d = st[s]
        r = d["rules"]
        a = d["arms"]
        pv = d["h_vs_g"]["voxels"]["after"]
        rows.append(
            {
                "sample": s,
                "(1) total fewer vs G (%)": f"{r['total_cut_pct']:.2f}"
                + (f"; seed {r['seed_cut_pct']:.2f}" if s == "held2" else ""),
                "(1)": "pass" if r["r1"] else "fail",
                "(2) broken - fixed vs G": pv["broken"] - pv["fixed"],
                "(2)": "pass" if r["r2"] else "fail",
                "(3) median shift (deg)": f(
                    a["H"]["after"]["median_err_deg"] - a["G"]["after"]["median_err_deg"], 5
                ),
                "(3)": "pass" if r["r3"] else "fail",
                "(4) extra wrong seeds (matched + unmatched); grid H/G": (
                    f"{d['wrong_seeds']['broad']}; {d['rescue']['H']['grid_fallbacks_left']}"
                    f"/{d['rescue']['G']['grid_fallbacks_left']}"
                ),
                "(4)": "pass" if r["r4"] else "fail",
                "(5) vs S8: fewer %; broken - fixed; median shift": house_cell(a["H"]),
                "(5)": "pass" if r["r5"] else "fail",
                "G passes (5)": "yes" if not r["g_fails_5"] else "no",
                "H passes": "pass" if r["passes"] else "FAIL",
            }
        )
    t["stallhit_rules"] = DT.markdown_table(rows, list(rows[0].keys()))
    rows = []
    for s in SAMPLES:
        v = st[s]["validity"]
        rows.append(
            {
                "sample": s,
                "finished": v["h_finished"],
                "voxels processed": n(v["n_processed"]),
                "aligned draws": v["aligned_draws"],
                "fits compared (matched steps)": n(v["fits_compared"]),
                "fits differ": n(v["fits_differ"]),
                "cut MC fits": n(v["fits_differ_cut_mc"]),
                "only one counted check": n(v["fits_differ_checks_only"]),
                "other": n(v["fits_differ_other"]),
                "every stop at 0.95 or more": v["every_stop_at_095_or_more"],
            }
        )
    t["stallhit_validity"] = DT.markdown_table(rows, list(rows[0].keys()))
    rows = []
    tot_ch = 0
    for s in SAMPLES:
        d = st[s]
        rows.append(
            {
                "sample": s,
                "taken seed ends changed, of matched steps (Wilson 95%)": wilson_cell(
                    d["taken_changed_n"], d["validity"]["matched_steps"]
                ),
            }
        )
        tot_ch += d["taken_changed_n"]
    t["stallhit_shares"] = DT.markdown_table(rows, list(rows[0].keys()))
    p1, p2, p3, p4 = pred["p1"], pred["p2"], pred["p3"], pred["p4"]
    rows = [
        {
            "prediction": "1 replay (gate)",
            **{s: f"{p1[s]['met']} ({p1[s]['changed']} changed)" for s in ("mid", "grad", "big")},
            "held": "n/a",
            "held2": "n/a",
        },
        {
            "prediction": "2 replay saving",
            **{s: f"{p2[s]['met']} ({100 * p2[s]['saving']:.2f}%)" for s in ("mid", "grad", "big")},
            "held": "n/a",
            "held2": "n/a",
        },
        {
            "prediction": "3 total fewer (mid 6-9, big 3.5-5.5, grad -0.5..0.5, held 2-6, "
            "held2 2-8; held2 seed 6-18)",
            **{s: f"{p3[s]['met']} ({p3[s]['total_cut_pct']:.2f}%)" for s in SAMPLES},
        },
        {
            "prediction": "4 broken - fixed vs G in -3..3, median shift <= 0.0005",
            **{s: f"{p4[s]['met']} ({p4[s]['net']})" for s in SAMPLES},
        },
        {
            "prediction": "5 mid passes rule 5",
            "mid": f"{pred['p5']['met']}",
            "grad": "n/a",
            "big": "n/a",
            "held": "n/a",
            "held2": "n/a",
        },
    ]
    t["stallhit_predictions"] = DT.markdown_table(
        rows, ["prediction", "held2", "held", "mid", "grad", "big"]
    )
    return t


def main() -> None:
    st = {s: json.loads((SH / f"stats_{s}.json").read_text()) for s in SAMPLES}
    for s in SAMPLES:
        st[s]["unmatched"] = json.loads((SH / f"unmatched_{s}.json").read_text())
        st[s]["unmatched"].pop("only_h", None)
        st[s]["unmatched"].pop("only_g", None)
    gate = json.loads((SH / "gate.json").read_text())
    plot(st)  # the plot comes first
    per = {s: bool(st[s]["rules"]["passes"]) for s in SAMPLES}
    verdict = SJ.step_verdict(per)
    status = (
        STATUS_PENDING
        if verdict
        else "study result of step 12: the rule failed; no key; not a candidate arm"
    )
    pareto(st, status)
    pred = predictions(st, gate)
    DT.write_tables(SH / "tables.md", tables(st, pred))
    failed = [
        f"{s}:({k[1]})"
        for s in SAMPLES
        for k in ("r1", "r2", "r3", "r4", "r5")
        if not st[s]["rules"][k] and not (k == "r5" and st[s]["rules"]["g_fails_5"])
    ]
    out = dict(
        per_sample_passes=per,
        step_passes=bool(verdict),
        verdict="READY FOR KEY" if verdict else "STOP: no key (study merge)",
        failed_rules=failed,
        rules={s: st[s]["rules"] for s in SAMPLES},
        predictions=pred,
        taken_changed={s: st[s]["taken_changed_n"] for s in SAMPLES},
        unmatched={s: st[s]["unmatched"] for s in SAMPLES},
        wrong_seeds_broad={s: st[s]["wrong_seeds"]["broad"] for s in SAMPLES},
        fits_differ_other={s: st[s]["validity"]["fits_differ_other"] for s in SAMPLES},
        status=status,
    )
    (SH / "summary.json").write_text(json.dumps(out, indent=1, default=float) + "\n")
    print(
        json.dumps({k: out[k] for k in ("per_sample_passes", "verdict", "failed_rules")}, indent=1)
    )
    print(json.dumps({k: v for k, v in pred.items() if k not in ("p1",)}, indent=1, default=float))
    print(json.dumps({s: v["met"] for s, v in pred["p1"].items()}))


if __name__ == "__main__":
    main()
