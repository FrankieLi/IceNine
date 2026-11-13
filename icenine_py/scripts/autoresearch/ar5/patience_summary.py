"""AR5 step 11, phase A step 5: validity, plot, tables, pareto rows, predictions, verdict.

  uv run python scripts/autoresearch/ar5/patience_summary.py

Reads ``patience/stats_<sample>.json`` (``patience_stats.py``) and ``patience/gate.json``. Writes
``ar5_patience_evals_wrong.png`` (made first: run evaluations against wrong voxels, S8, G and S),
``tables.md`` (blocks for ``sync_doc_tables.py``), ``summary.json`` (every number of the tables,
the predictions, the rules, the verdict) and ``pareto_rows.json`` (the format of
``ar5/rescue/pareto_rows.json``, arm ``S`` rows). No p-values; paired counts are by voxel.
"""

import json
import sys
from pathlib import Path
from typing import Any, Dict, List

HERE = Path(__file__).resolve().parent
PAT = HERE.parents[2] / "benchmarks" / "autoresearch" / "ar5" / "patience"
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent.parent / "common"))

import doc_tables as DT  # noqa: E402
import patience_judge as PJ  # noqa: E402

SAMPLES = ("held", "mid", "grad", "big")
N = 24570


def f(x: float, d: int = 5) -> str:
    return f"{x:.{d}f}"


def n(x: int) -> str:
    return f"{int(x):,}"


def plot(st: Dict[str, Any]) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.6))
    col = dict(held="tab:red", mid="tab:blue", grad="tab:green", big="tab:orange")
    for ax, tag, title in zip(axes, ("before", "after"), ("before the pass", "after the pass")):
        for s in SAMPLES:
            a = st[s]["arms"]
            key = "before_pass" if tag == "before" else "total"
            pts = [
                ("S8", "^", a["S"]["stages_s8"][key], a["S"][f"s8_{tag}"]),
                ("G", "o", a["G"]["stages"][key], a["G"][tag]),
                ("S", "s", a["S"]["stages"][key], a["S"][tag]),
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
        ax.set_xlabel(
            "run evaluations (millions, log%s)" % (", with the pass" if tag == "after" else "")
        )
        ax.set_ylabel("wrong voxels (more than 1 degree), Wilson 95%")
        ax.set_title(title)
        ax.grid(alpha=0.3)
    fig.legend(loc="lower center", fontsize=7, ncol=6, frameon=False)
    fig.suptitle(
        "Stored S8 (triangle), arm G = step 9 (circle), arm S = G with patience P(3000) (square)",
        fontsize=10,
    )
    fig.tight_layout(rect=(0, 0.12, 1, 1))
    fig.savefig(PAT / "ar5_patience_evals_wrong.png", dpi=110)


def pareto(st: Dict[str, Any]) -> None:
    rows: List[Dict[str, Any]] = []
    for s in SAMPLES:
        d = st[s]["arms"]["S"]
        for after in (False, True):
            side = d["after" if after else "before"]
            total = d["stages"]["total" if after else "before_pass"]
            rows.append(
                dict(
                    sample=s,
                    arm="S",
                    order="interleaved + rescue G1 + patience P(3000)",
                    boundary_pass=after,
                    evals_total=total,
                    evals_per_voxel=round(total / N, 2),
                    wrong=side["wrong"],
                    n=N,
                    interval_pct=[
                        round(100 * side["wilson_lo"], 2),
                        round(100 * side["wilson_hi"], 2),
                    ],
                    source="icenine_py/benchmarks/autoresearch/ar5/patience/stats_%s.json" % s,
                )
            )
    (PAT / "pareto_rows.json").write_text(json.dumps(rows, indent=1) + "\n")


def predictions(st: Dict[str, Any], gate: Dict[str, Any]) -> Dict[str, Any]:
    out: Dict[str, Any] = {}
    # 1 and 2: the replay (gate.json)
    ps = gate["per_sample"]
    out["p1"] = {
        s: dict(
            met=all(v for v in ps[s]["prediction1"].values() if isinstance(v, bool)),
            detail=ps[s]["prediction1"],
        )
        for s in ("mid", "grad", "big")
    }
    out["p2"] = {
        s: dict(met=ps[s]["prediction2"]["saving_in_range"], saving=ps[s]["prediction2"]["saving"])
        for s in ("mid", "grad", "big")
    }
    h, g = st["held"]["arms"], st["held"]["arms"]["G"]
    seed_cut = PJ.pct_fewer(h["S"]["stages"]["seed_step"], g["stages"]["seed_step"])
    tot_cut = PJ.pct_fewer(h["S"]["stages"]["total"], g["stages"]["total"])
    out["p3"] = dict(
        seed_cut_pct=seed_cut,
        total_cut_pct=tot_cut,
        met=bool(10 <= seed_cut <= 30 and 4 <= tot_cut <= 12),
    )
    rng = dict(mid=(8, 12), big=(4, 7), grad=(-1, 1))
    p4 = {}
    for s, (lo, hi) in rng.items():
        a = st[s]["arms"]
        c = PJ.pct_fewer(a["S"]["stages"]["total"], a["G"]["stages"]["total"])
        p4[s] = dict(total_cut_pct=c, met=bool(lo <= c <= hi))
    out["p4"] = p4
    p5 = {}
    for s in SAMPLES:
        a = st[s]["arms"]
        pv = st[s]["s_vs_g"]["voxels"]["after"]
        net = pv["broken"] - pv["fixed"]
        dm = a["S"]["after"]["median_err_deg"] - a["G"]["after"]["median_err_deg"]
        p5[s] = dict(net=net, median_shift_deg=dm, met=bool(-4 <= net <= 4 and abs(dm) <= 0.0005))
    out["p5"] = p5
    return out


def house_cell(s: Dict[str, Any]) -> str:
    """Rule (5) numbers of arm S against S8: percent fewer; broken minus fixed; median shift."""
    fewer = PJ.pct_fewer(s["stages"]["total"], s["stages_s8"]["total"])
    net = s["paired_after"]["broken"] - s["paired_after"]["fixed"]
    shift = s["after"]["median_err_deg"] - s["s8_after"]["median_err_deg"]
    return f"{fewer:.1f}; {net}; {shift:.5f}"


def tables(st: Dict[str, Any], pred: Dict[str, Any]) -> Dict[str, str]:
    t: Dict[str, str] = {}
    rows = []
    for s in SAMPLES:
        a = st[s]["arms"]
        for arm in ("G", "S"):
            sp = a[arm]["stages_split"]
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
    t["patience_stages"] = DT.markdown_table(
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
                "cut": n(c["cut_refinements"]),
                "evaluations in cut ones": n(c["evals_in_cut"]),
                "cut end under 0.95": n(c["cut_end_hit_under_095"]),
                "cut end at least 0.95": n(c["cut_end_hit_at_least_095"]),
                "taken seed ends changed": n(d["taken_changed_n"]),
                "median angle to G (deg)": f(d["taken_changed_median_mis_deg"] or 0.0, 4),
                "S wrong where G right": n(d["wrong_seeds"]["s_wrong_where_g_right"]),
            }
        )
    t["patience_cuts"] = DT.markdown_table(rows, list(rows[0].keys()))
    rows = []
    for s in SAMPLES:
        d = st[s]
        for c in d["taken_changed"]:
            pass
        ch = d["taken_changed"]
        if not ch:
            continue
        import numpy as np

        rows.append(
            {
                "sample": s,
                "changed": len(ch),
                "was at 0.95": sum(1 for c in ch if c["stopped_g"]),
                "now at 0.95": sum(1 for c in ch if c["hit_s"] >= 0.95),
                "median hit G": f(float(np.median([c["hit_g"] for c in ch])), 3),
                "median hit S": f(float(np.median([c["hit_s"] for c in ch])), 3),
                "median cost G": f(float(np.median([c["cost_g"] for c in ch])), 4),
                "median cost S": f(float(np.median([c["cost_s"] for c in ch])), 4),
                "median error G (deg)": f(float(np.median([c["err_g"] for c in ch])), 4),
                "median error S (deg)": f(float(np.median([c["err_s"] for c in ch])), 4),
                "max error S (deg)": f(max(c["err_s"] for c in ch), 4),
            }
        )
    t["patience_changed"] = DT.markdown_table(rows, list(rows[0].keys()))
    rows = []
    for s in SAMPLES:
        d = st[s]
        a = d["arms"]
        r = d["rescue"]
        for arm in ("G", "S"):
            sd = a[arm]
            side = sd["after"]
            rows.append(
                {
                    "sample": s,
                    "arm": arm,
                    "failed seed tests": n(r[arm]["fallback_calls"]),
                    "rescued": n(r[arm]["rescued"]),
                    "grid searches": n(r[arm]["grid_fallbacks_left"]),
                    "wrong after pass": n(side["wrong"]),
                    "Wilson 95% (%)": "%.2f to %.2f"
                    % (100 * side["wilson_lo"], 100 * side["wilson_hi"]),
                    "median error (deg)": f(side["median_err_deg"], 5),
                }
            )
    t["patience_wrong"] = DT.markdown_table(rows, list(rows[0].keys()))
    rows = []
    for s in SAMPLES:
        d = st[s]
        a = d["arms"]["S"]
        pg = d["s_vs_g"]["voxels"]
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
    t["patience_paired"] = DT.markdown_table(rows, list(rows[0].keys()))
    rows = []
    for s in SAMPLES:
        d = st[s]
        r = d["rules"]
        a = d["arms"]
        pv = d["s_vs_g"]["voxels"]["after"]
        rows.append(
            {
                "sample": s,
                "(1) total fewer vs G (%)": f"{r['total_cut_pct']:.2f}"
                + (f"; seed {r['seed_cut_pct']:.2f}" if s == "held" else ""),
                "(1)": "pass" if r["r1"] else "fail",
                "(2) broken - fixed vs G": pv["broken"] - pv["fixed"],
                "(2)": "pass" if r["r2"] else "fail",
                "(3) median shift (deg)": f(
                    a["S"]["after"]["median_err_deg"] - a["G"]["after"]["median_err_deg"], 5
                ),
                "(3)": "pass" if r["r3"] else "fail",
                "(4) extra wrong seeds; grid S/G": f"{d['wrong_seeds']['s_wrong_where_g_right']}; "
                f"{d['rescue']['S']['grid_fallbacks_left']}"
                f"/{d['rescue']['G']['grid_fallbacks_left']}",
                "(4)": "pass" if r["r4"] else "fail",
                "(5) vs S8: fewer %; broken - fixed; median shift": house_cell(a["S"]),
                "(5)": "pass" if r["r5"] else "fail",
                "G passes (5)": "yes" if not r["g_fails_5"] else "no",
                "S passes": "pass" if r["passes"] else "FAIL",
            }
        )
    t["patience_rules"] = DT.markdown_table(rows, list(rows[0].keys()))
    rows = []
    for s in SAMPLES:
        v = st[s]["validity"]
        rows.append(
            {
                "sample": s,
                "finished": v["s_finished"],
                "voxels processed": n(v["n_processed"]),
                "aligned draws": v["aligned_draws"],
                "seed steps G / S": f"{n(v['seed_steps_g'])} / {n(v['seed_steps_s'])}",
                "matched steps": n(v["matched_steps"]),
                "fits compared": n(v["fits_compared"]),
                "fits differ": n(v["fits_differ"]),
                "of them cut MC fits": n(v["fits_differ_cut_mc"]),
                "other": n(v["fits_differ_other"]),
            }
        )
    t["patience_validity"] = DT.markdown_table(rows, list(rows[0].keys()))
    rows = [
        {
            "prediction": "1 replay (gate)",
            "mid": pred["p1"]["mid"]["met"],
            "grad": pred["p1"]["grad"]["met"],
            "big": pred["p1"]["big"]["met"],
            "held": "n/a",
        },
        {
            "prediction": "2 replay saving",
            "mid": pred["p2"]["mid"]["met"],
            "grad": pred["p2"]["grad"]["met"],
            "big": pred["p2"]["big"]["met"],
            "held": "n/a",
        },
        {
            "prediction": "3 held seed 10-30% / total 4-12% fewer",
            "mid": "n/a",
            "grad": "n/a",
            "big": "n/a",
            "held": f"{pred['p3']['met']} ({pred['p3']['seed_cut_pct']:.1f}%; "
            f"{pred['p3']['total_cut_pct']:.1f}%)",
        },
        {
            "prediction": "4 total fewer (mid 8-12, big 4-7, grad -1..1)",
            "mid": f"{pred['p4']['mid']['met']} ({pred['p4']['mid']['total_cut_pct']:.1f}%)",
            "grad": f"{pred['p4']['grad']['met']} ({pred['p4']['grad']['total_cut_pct']:.2f}%)",
            "big": f"{pred['p4']['big']['met']} ({pred['p4']['big']['total_cut_pct']:.1f}%)",
            "held": "n/a",
        },
        {
            "prediction": "5 broken - fixed vs G in -4..4, median shift <= 0.0005",
            "mid": f"{pred['p5']['mid']['met']} ({pred['p5']['mid']['net']})",
            "grad": f"{pred['p5']['grad']['met']} ({pred['p5']['grad']['net']})",
            "big": f"{pred['p5']['big']['met']} ({pred['p5']['big']['net']})",
            "held": f"{pred['p5']['held']['met']} ({pred['p5']['held']['net']})",
        },
    ]
    t["patience_predictions"] = DT.markdown_table(
        rows, ["prediction", "held", "mid", "grad", "big"]
    )
    return t


def main() -> None:
    st = {s: json.loads((PAT / f"stats_{s}.json").read_text()) for s in SAMPLES}
    gate = json.loads((PAT / "gate.json").read_text())
    plot(st)  # the plot comes first
    pareto(st)
    pred = predictions(st, gate)
    DT.write_tables(PAT / "tables.md", tables(st, pred))
    per = {s: bool(st[s]["rules"]["passes"]) for s in SAMPLES}
    verdict = PJ.step_verdict(per)
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
        cut_end_hit_under_095={
            s: dict(
                cut=st[s]["cuts"]["cut_refinements"],
                under=st[s]["cuts"]["cut_end_hit_under_095"],
                at_least=st[s]["cuts"]["cut_end_hit_at_least_095"],
            )
            for s in SAMPLES
        },
    )
    (PAT / "summary.json").write_text(json.dumps(out, indent=1, default=float) + "\n")
    print(
        json.dumps({k: out[k] for k in ("per_sample_passes", "verdict", "failed_rules")}, indent=1)
    )


if __name__ == "__main__":
    main()
