"""AR5 step 8, phase A step 5: validity, tables, plot, predictions, decision rules, verdict.

  uv run python scripts/autoresearch/ar5/order_summary.py

Reads ``benchmarks/autoresearch/ar5/order/stats_<sample>.json`` (made by ``order_stats.py``; no
sample needed here) and the ``order.json`` / ``run.json`` of each run folder for the validity checks
(no ``capped.json``; the commit and the dirty file list the run recorded; the frozen-file check
now). Writes in the same folder: ``summary.json``, ``tables.md`` (blocks for
``scripts/dev/sync_doc_tables.py``), ``pareto_rows.json`` and the plot
``ar5_order_evals_wrong.png``.
Counts and intervals only. The predictions and the decision rules are the ones written before the
runs (PLAN.md, AR5 step 8); nothing is tuned.
"""

import json
import subprocess
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional

HERE = Path(__file__).resolve().parent
ICE = HERE.parents[2]
ORDER = ICE / "benchmarks" / "autoresearch" / "ar5" / "order"
CACHE = ORDER / "cache"
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(ICE / "scripts" / "common"))

import doc_tables as DT  # noqa: E402
import fit_rules as FR  # noqa: E402
import order_freeze as OF  # noqa: E402

SAMPLES = ["held", "mid", "grad", "big"]
N = 24570
EPS = 1e-12
PRED = {
    "held": dict(seed=(20, 45), total=(6, 16)),
    "mid": dict(seed=(30, 45), total=(11, 18)),
    "grad": dict(seed=(10, 25), total=(1, 5)),
    "big": dict(seed=(60, 80), total=(12, 22)),
}


# column headers of tables.md (review fix 5: no unexplained abbreviations; the wall column is
# labelled as contended and not a timing result)
HEADS = dict(
    sample="sample",
    arm="arm",
    steps="seed steps",
    seed_evals="seed-step evaluations",
    mc="MC refinements",
    mc_steps="seed steps with MC",
    nb="neighbour fits",
    rv="revisits",
    grid_e="grid fallback evaluations",
    before="before the pass",
    pas="boundary pass",
    total="total with the pass",
    grid="grid fallbacks (count)",
    wb="wrong before the pass (95% interval)",
    wa="wrong after the pass (95% interval)",
    med="median error after the pass (deg)",
    f8="fixed / broken against S8 (after the pass)",
    fewer="fewer evaluations than S8 (with the pass)",
    r5="rule 5 (house rule against S8)",
    wall="wall, run + pass (minutes; contended, not a timing result)",
    seed="seed-step evaluations, A to Y",
    sc="seed-step cut",
    tc="total cut (with the pass)",
    fb="fixed / broken, Y against A (after the pass)",
    fbb="fixed / broken, Y against A (before the pass)",
    wrong="wrong after the pass, A to Y",
    gridc="grid fallbacks, A to Y",
    medc="median error after the pass (deg), A to Y",
    dif="matched seed steps whose taken end differs (wrong ends in A, in Y)",
    only="seed steps in one arm only (A / Y)",
    r1="rule 1: seed step at least 10% fewer and total fewer",
    r2="rule 2: broken minus fixed at most 5 (net)",
    r3="rule 3: median error at most 0.001 deg higher",
    r4="rule 4: at most 2 more grid fallbacks",
    r5y="rule 5 for Y: house rule against S8",
    a5="rule 5 for A",
    y="Y passes on the sample",
    p12="predictions 1-2: seed-step and total cut in range",
    p3="prediction 3: grid fallbacks within 2, net within -3 to 3",
    p4="prediction 4: median moves at most 0.0005 deg",
    p5="prediction 5: shared fits identical",
)
# the same short key means different columns in two tables; map them per table
_ALIAS = {("grid", "y_vs_a"): "gridc", ("med", "y_vs_a"): "medc", ("r5", "rules"): "r5y"}
TERMS = [
    ("A", "present (interleaved) order: CMA fit, then MC if needed, candidate by candidate"),
    ("Y", "sweep order: CMA fit on every candidate first, MC only if no CMA end reaches 0.95"),
    ("Yl", "sweep order with the library's own random draws (not aligned), big only"),
    ("S8", "the stored house reference run (other seed rule, library draws; not draw-aligned)"),
    ("wrong", "voxel error above 1 degree from the truth (symmetry-reduced angle)"),
    ("fixed / broken", "fixed: wrong in the first run, right in the second; broken: the reverse"),
    ("net", "broken minus fixed"),
    ("matched seed step", "a seed step at the same voxel and occurrence in both arms"),
    ("taken end", "the end of the refinement the seed step kept for the voxel"),
    ("the pass", "the boundary pass after the layer (run D4 of step 2)"),
]


def labelled(
    rows: List[Dict[str, Any]],
    cols: List[str],
    formats: Optional[Dict[str, str]] = None,
    table: str = "",
) -> str:
    """``DT.markdown_table`` with the headers of ``HEADS`` instead of the short row keys
    (``table`` selects the per-table header of a key used in two tables, see ``_ALIAS``)."""
    head = {c: HEADS[_ALIAS.get((c, table), c)] for c in cols}
    new = [{head[c]: r.get(c) for c in cols} for r in rows]
    fm = {head[c]: f for c, f in (formats or {}).items()}
    return DT.markdown_table(new, [head[c] for c in cols], fm)


def pct(x: float) -> float:
    return round(100.0 * x, 1)


def y_passes(rules: Dict[str, bool], a_fails_5: bool) -> bool:
    """The plan's rule for one sample: (1) to (4) against A, and (5), the house rule against S8,
    unless A itself fails (5) there; then Y is judged by (1) to (4) alone."""
    r14 = rules["r1"] and rules["r2"] and rules["r3"] and rules["r4"]
    return bool(r14 and (rules["r5"] or a_fails_5))


def step_passes(per_sample: Dict[str, bool]) -> bool:
    """The step passes if Y passes on held and on all three disclosed samples."""
    return bool(per_sample["held"] and all(per_sample[s] for s in SAMPLES))


def validity(st: Dict[str, Dict[str, Any]]) -> Dict[str, Any]:
    out: Dict[str, Any] = dict(freeze_check_now=OF.problems() or "passed")
    heads: Dict[str, Any] = {}
    capped: List[str] = []
    for s in SAMPLES:
        for a in st[s]["arms"]:
            d = CACHE / f"{a}_{s}"
            if (d / "capped.json").exists() or not (d / "run.json").exists():
                capped.append(f"{a}_{s}")
            o = json.loads((d / "order.json").read_text())
            heads[f"{a}_{s}"] = dict(
                head=o["state"]["head"][:7],
                dirty=[x for x in o["state"]["dirty_library_or_scripts"] if not x.startswith("??")],
            )
    out["capped_or_missing"] = capped
    out["recorded_heads"] = heads
    out["distinct_heads"] = sorted({h["head"] for h in heads.values()})
    # files frozen: unchanged between the earliest and the latest recorded head
    diffs: List[str] = []
    paths = [f"icenine_py/scripts/autoresearch/ar5/{x}" for x in OF.SCRIPTS] + [
        "icenine_py/icenine",
        "icenine_py/pyproject.toml",
        "icenine_py/uv.lock",
    ]
    for h in out["distinct_heads"]:
        r = subprocess.run(
            ["git", "diff", "--stat", h, "HEAD", "--"] + paths,
            cwd=str(ICE),
            capture_output=True,
            text=True,
        )
        diffs.append(f"{h}: {r.stdout.strip() or 'no change in the frozen files or the library'}")
    out["frozen_files_changed_since_head"] = diffs
    out["bad_counts"] = {
        s: {
            a: {
                k: st[s]["arms"][a]["log"][k]
                for k in ("bad_mc_after_cma_pass", "bad_cma_after_cma_pass")
            }
            for a in st[s]["arms"]
        }
        for s in SAMPLES
    }
    out["shared_fits"] = {
        s: (
            st[s]["pairs"]["Y_vs_A"]["seed_steps"]["shared_fits_compared"],
            st[s]["pairs"]["Y_vs_A"]["seed_steps"]["shared_fits_different"],
        )
        for s in SAMPLES
    }
    out["all_ok"] = bool(
        not capped
        and out["freeze_check_now"] == "passed"
        and all(v == 0 for s in out["bad_counts"].values() for a in s.values() for v in a.values())
        and all(x[1] == 0 for x in out["shared_fits"].values())
        and all("no change" in d for d in diffs)
    )
    return out


def main() -> None:
    st = {s: json.loads((ORDER / f"stats_{s}.json").read_text()) for s in SAMPLES}
    val = validity(st)
    rows: Dict[str, Any] = {}
    for s in SAMPLES:
        A, Y = st[s]["arms"]["A"], st[s]["arms"]["Y"]
        pv = st[s]["pairs"]["Y_vs_A"]
        r: Dict[str, Any] = dict(sample=s)
        r["seed_A"], r["seed_Y"] = A["stages"]["seed_step"], Y["stages"]["seed_step"]
        r["total_A"], r["total_Y"] = A["stages"]["total"], Y["stages"]["total"]
        r["before_A"], r["before_Y"] = A["stages"]["before_pass"], Y["stages"]["before_pass"]
        r["seed_cut"] = 1 - r["seed_Y"] / r["seed_A"]
        r["total_cut"] = 1 - r["total_Y"] / r["total_A"]
        r["fixed"], r["broken"] = pv["voxels"]["after"]["fixed"], pv["voxels"]["after"]["broken"]
        r["fixed_b"], r["broken_b"] = (
            pv["voxels"]["before"]["fixed"],
            pv["voxels"]["before"]["broken"],
        )
        r["net"] = r["broken"] - r["fixed"]
        r["wrong_A"], r["wrong_Y"] = A["after"]["wrong"], Y["after"]["wrong"]
        r["wrong_A_b"], r["wrong_Y_b"] = A["before"]["wrong"], Y["before"]["wrong"]
        r["grid_A"], r["grid_Y"] = A["grid_fallbacks"], Y["grid_fallbacks"]
        r["med_A"], r["med_Y"] = A["after"]["median_err_deg"], Y["after"]["median_err_deg"]
        r["differs"] = pv["seed_steps"]["taken_differs"]
        r["differs_wrong_A"] = pv["seed_steps"]["differing_wrong_first"]
        r["differs_wrong_Y"] = pv["seed_steps"]["differing_wrong_second"]
        r["matched"] = pv["seed_steps"]["matched"]
        r["only_A"], r["only_Y"] = pv["seed_steps"]["only_a"], pv["seed_steps"]["only_b"]
        # decision rules, as written before the runs
        rA, rY = A["rule"], Y["rule"]
        r["rule5_A"], r["rule5_Y"] = rA, rY
        d1 = r["seed_cut"] >= 0.10 - EPS and r["total_Y"] < r["total_A"]
        d2 = r["net"] <= 5
        d3 = r["med_Y"] <= r["med_A"] + 0.001 + EPS
        d4 = r["grid_Y"] - r["grid_A"] <= 2
        d5 = bool(rY["passes"])
        r["rules"] = dict(r1=bool(d1), r2=bool(d2), r3=bool(d3), r4=bool(d4), r5=d5)
        r["a_fails_5"] = not bool(rA["passes"])
        r["y_passes"] = y_passes(r["rules"], r["a_fails_5"])
        # predictions
        p = PRED[s]
        r["pred1_2"] = dict(
            seed_cut_pct=pct(r["seed_cut"]),
            total_cut_pct=pct(r["total_cut"]),
            met=bool(
                p["seed"][0] <= 100 * r["seed_cut"] <= p["seed"][1]
                and p["total"][0] <= 100 * r["total_cut"] <= p["total"][1]
            ),
        )
        r["pred3"] = bool(abs(r["grid_Y"] - r["grid_A"]) <= 2 and -3 <= r["net"] <= 3)
        r["pred4"] = bool(abs(r["med_Y"] - r["med_A"]) <= 0.0005 + EPS)
        r["pred5"] = pv["seed_steps"]["shared_fits_different"] == 0
        rows[s] = r
    verdict = step_passes({s: rows[s]["y_passes"] for s in SAMPLES})
    big = st["big"]
    yl = dict(
        seed_Yl=big["arms"]["Yl"]["stages"]["seed_step"],
        seed_Y=big["arms"]["Y"]["stages"]["seed_step"],
        total_Yl=big["arms"]["Yl"]["stages"]["total"],
        total_Y=big["arms"]["Y"]["stages"]["total"],
        wrong_Yl=big["arms"]["Yl"]["after"]["wrong"],
        wrong_Y=big["arms"]["Y"]["after"]["wrong"],
        grid_Yl=big["arms"]["Yl"]["grid_fallbacks"],
        grid_Y=big["arms"]["Y"]["grid_fallbacks"],
        vs_Y=big["pairs"]["Yl_vs_Y"]["voxels"],
        vs_A=big["pairs"]["Yl_vs_A"]["voxels"],
    )
    out = dict(validity=val, rows=rows, step_passes=verdict, yl_vs_y_big=yl)
    (ORDER / "summary.json").write_text(json.dumps(out, indent=1, default=float) + "\n")
    tables(st, rows, val, yl)
    pareto(st)
    plot(st)
    print(
        json.dumps(
            {
                "valid": val["all_ok"],
                "verdict": verdict,
                **{s: {k: rows[s][k] for k in ("seed_cut", "total_cut", "y_passes")} for s in rows},
            },
            indent=1,
        )
    )


def tables(
    st: Dict[str, Any], rows: Dict[str, Any], val: Dict[str, Any], yl: Dict[str, Any]
) -> None:
    T: Dict[str, str] = {}
    main_rows: List[Dict[str, Any]] = []
    for s in SAMPLES:
        for a, d in st[s]["arms"].items():
            lg = d["log"]
            main_rows.append(
                dict(
                    sample=s,
                    arm=a,
                    seed_evals=d["stages"]["seed_step"],
                    mc=lg["mc_fits"],
                    mc_steps=lg["steps_with_mc"],
                    steps=lg["seed_steps"],
                    nb=d["stages"]["neighbour"],
                    rv=d["stages"]["revisit"],
                    grid_e=d["stages"]["grid_fallback"],
                    before=d["stages"]["before_pass"],
                    pas=d["stages"]["pass_"],
                    total=d["stages"]["total"],
                    grid=d["grid_fallbacks"],
                )
            )
    T["order_stages"] = labelled(
        main_rows,
        ["sample", "arm", "steps", "seed_evals", "mc", "mc_steps", "nb", "rv", "grid_e", "before",
         "pas", "total", "grid"],
        {k: "," for k in ("seed_evals", "nb", "rv", "grid_e", "before", "pas", "total")},
    )  # fmt: skip
    wr: List[Dict[str, Any]] = []
    for s in SAMPLES:
        for a, d in st[s]["arms"].items():
            wr.append(
                dict(
                    sample=s,
                    arm=a,
                    wb=f"{d['before']['wrong']} ({100*d['before']['wilson_lo']:.2f} to "
                    f"{100*d['before']['wilson_hi']:.2f}%)",
                    wa=f"{d['after']['wrong']} ({100*d['after']['wilson_lo']:.2f} to "
                    f"{100*d['after']['wilson_hi']:.2f}%)",
                    med=d["after"]["median_err_deg"],
                    f8=f"{d['paired_after']['fixed']} / {d['paired_after']['broken']}",
                    fewer=f"{pct(d['rule']['fewer'])}%",
                    r5="pass" if d["rule"]["passes"] else "fail",
                    wall=f"{(d['wall_s']+d['wall_pass_s'])/60:.0f}",
                )
            )
    T["order_accuracy"] = labelled(
        wr, ["sample", "arm", "wb", "wa", "med", "f8", "fewer", "r5", "wall"], {"med": ".5f"}
    )
    yr = []
    for s in SAMPLES:
        r = rows[s]
        yr.append(
            dict(
                sample=s,
                seed=f"{r['seed_A']:,} to {r['seed_Y']:,}",
                sc=f"{pct(r['seed_cut'])}%",
                tc=f"{pct(r['total_cut'])}%",
                fb=f"{r['fixed']} / {r['broken']}",
                fbb=f"{r['fixed_b']} / {r['broken_b']}",
                wrong=f"{r['wrong_A']} to {r['wrong_Y']}",
                grid=f"{r['grid_A']} to {r['grid_Y']}",
                med=f"{r['med_A']:.5f} to {r['med_Y']:.5f}",
                dif=f"{r['differs']} of {r['matched']} ({r['differs_wrong_A']} wrong in A, "
                f"{r['differs_wrong_Y']} in Y)",
                only=f"{r['only_A']} / {r['only_Y']}",
            )
        )
    T["order_y_vs_a"] = labelled(
        yr,
        ["sample", "seed", "sc", "tc", "fb", "fbb", "wrong", "grid", "med", "dif", "only"],
        table="y_vs_a",
    )
    dr = []
    for s in SAMPLES:
        r = rows[s]
        ru = r["rules"]
        dr.append(
            dict(
                sample=s,
                r1="pass" if ru["r1"] else "fail",
                r2=f"{'pass' if ru['r2'] else 'fail'} ({r['net']})",
                r3="pass" if ru["r3"] else "fail",
                r4="pass" if ru["r4"] else "fail",
                r5=f"{'pass' if ru['r5'] else 'fail'} ({pct(r['rule5_Y']['fewer'])}% fewer, net "
                f"{r['rule5_Y']['net_broken']}, median {r['rule5_Y']['med_diff']:+.5f})",
                a5=f"{'fail' if r['a_fails_5'] else 'pass'} ({pct(r['rule5_A']['fewer'])}% fewer, "
                f"net {r['rule5_A']['net_broken']})",
                y="pass" if r["y_passes"] else "fail",
            )
        )
    T["order_rules"] = labelled(
        dr, ["sample", "r1", "r2", "r3", "r4", "r5", "a5", "y"], table="rules"
    )
    pr = []
    for s in SAMPLES:
        r = rows[s]
        pr.append(
            dict(
                sample=s,
                p12=f"{'met' if r['pred1_2']['met'] else 'not met'} "
                f"(seed {r['pred1_2']['seed_cut_pct']}%"
                f", total {r['pred1_2']['total_cut_pct']}%; range {PRED[s]['seed'][0]}-"
                f"{PRED[s]['seed'][1]} / {PRED[s]['total'][0]}-{PRED[s]['total'][1]})",
                p3=f"{'met' if r['pred3'] else 'not met'} (grid {r['grid_A']} to {r['grid_Y']}, "
                f"net {r['net']})",
                p4=f"{'met' if r['pred4'] else 'not met'} ({r['med_Y']-r['med_A']:+.5f})",
                p5=f"{'met' if r['pred5'] else 'not met'} "
                f"({val['shared_fits'][s][1]} of {val['shared_fits'][s][0]} shared fits differ)",
            )
        )
    T["order_predictions"] = labelled(pr, ["sample", "p12", "p3", "p4", "p5"])
    T["order_terms"] = DT.markdown_table(
        [dict(term=k, meaning=v) for k, v in TERMS], ["term", "meaning"]
    )
    T["order_yl"] = DT.markdown_table(
        [
            dict(
                what="seed-step evaluations",
                Y=f"{yl['seed_Y']:,}",
                Yl=f"{yl['seed_Yl']:,}",
            ),
            dict(what="total with the pass", Y=f"{yl['total_Y']:,}", Yl=f"{yl['total_Yl']:,}"),
            dict(what="wrong after the pass", Y=yl["wrong_Y"], Yl=yl["wrong_Yl"]),
            dict(what="grid fallbacks", Y=yl["grid_Y"], Yl=yl["grid_Yl"]),
            dict(
                what="fixed / broken, Yl against Y (after the pass)",
                Y="",
                Yl=f"{yl['vs_Y']['after']['fixed']} / {yl['vs_Y']['after']['broken']}",
            ),
            dict(
                what="fixed / broken, Yl against A (after the pass)",
                Y="",
                Yl=f"{yl['vs_A']['after']['fixed']} / {yl['vs_A']['after']['broken']}",
            ),
        ],
        ["what", "Y", "Yl"],
    )
    DT.write_tables(ORDER / "tables.md", T)


def pareto(st: Dict[str, Any]) -> None:
    rows: List[Dict[str, Any]] = []
    for s in SAMPLES:
        for a in ("A", "Y"):
            d = st[s]["arms"][a]
            for after in (False, True):
                side = d["after" if after else "before"]
                total = d["stages"]["total" if after else "before_pass"]
                rows.append(
                    dict(
                        sample=s,
                        arm=a,
                        order="interleaved" if a == "A" else "sweep",
                        boundary_pass=after,
                        evals_total=total,
                        evals_per_voxel=round(total / N, 2),
                        wrong=side["wrong"],
                        n=N,
                        interval_pct=[
                            round(100 * side["wilson_lo"], 2),
                            round(100 * side["wilson_hi"], 2),
                        ],
                        source="icenine_py/benchmarks/autoresearch/ar5/order/stats_%s.json" % s,
                    )
                )
    (ORDER / "pareto_rows.json").write_text(json.dumps(rows, indent=1) + "\n")


def plot(st: Dict[str, Any]) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.6))
    col = dict(held="tab:red", mid="tab:blue", grad="tab:green", big="tab:orange")
    # review suggestion: Wilson 95% bars on the wrong counts, and the stored S8 point per sample
    # (S8 is not draw-aligned with A and Y); Yl on big shows one changed draw stream
    for ax, tag, title in zip(axes, ("before", "after"), ("before the pass", "after the pass")):
        for s in SAMPLES:
            arms = st[s]["arms"]
            pts = [("A", "o", arms["A"]["stages"]["seed_step"], arms["A"][tag])]
            pts.append(("Y", "s", arms["Y"]["stages"]["seed_step"], arms["Y"][tag]))
            if "Yl" in arms:
                pts.append(("Yl", "D", arms["Yl"]["stages"]["seed_step"], arms["Yl"][tag]))
            pts.append(("S8", "^", arms["A"]["stages_s8"]["seed_step"], arms["A"][f"s8_{tag}"]))
            for a, mk, x, side in pts:
                k = side["wrong"]
                lo, hi = N * side["wilson_lo"], N * side["wilson_hi"]
                ax.errorbar(
                    x / 1e6,
                    k,
                    yerr=[[k - lo], [hi - k]],
                    fmt=mk,
                    color=col[s],
                    mfc="none" if a in ("S8", "Yl") else col[s],
                    ms=6,
                    capsize=2,
                    lw=0.8,
                    label=f"{s} {a}" if tag == "before" else None,
                )
        ax.set_xscale("log")
        ax.set_xlabel("seed-step evaluations (millions, log)")
        ax.set_ylabel("wrong voxels (more than 1 degree), Wilson 95%")
        ax.set_title(title)
        ax.grid(alpha=0.3)
    fig.legend(loc="lower center", fontsize=7, ncol=7, frameon=False)
    fig.suptitle(
        "Present order A (circle), sweep order Y (square), Yl on big (diamond), stored S8 "
        "(triangle)",
        fontsize=10,
    )
    fig.tight_layout(rect=(0, 0.12, 1, 1))
    fig.savefig(ORDER / "ar5_order_evals_wrong.png", dpi=110)


if __name__ == "__main__":
    main()
