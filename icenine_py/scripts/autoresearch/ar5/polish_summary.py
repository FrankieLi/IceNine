"""AR5 step 14, phase A step 6: the plot, the tables, the pareto rows, the predictions, the verdict.

  uv run python scripts/autoresearch/ar5/polish_summary.py

Reads ``polish/stats_<sample>.json`` and ``polish/unmatched_<sample>.json`` (``polish_stats.py``),
``polish/gate1.json`` and ``polish/gate2.json``, and the stored tables of earlier steps for the
other arms (``ar5/order``, ``ar5/rescue``, ``ar5/stallhit`` ``pareto_rows.json`` and
``stats_<sample>.json``). Writes ``ar5_polish_evals_wrong.png`` (made first: run evaluations
against wrong voxels, S8, A, G, H and J, before and after the pass), ``tables.md`` (blocks for
``sync_doc_tables.py``), ``summary.json`` and ``pareto_rows.json`` (the format of
``ar5/stallhit/pareto_rows.json``, with a status field). ``held2`` after the pass is NOT read or
shown: it waits for step 13's ``held2`` result (plan step 7). No p-values; paired counts are by
voxel. Plain wording: "local optimization", "grid search", "realistic" data.
"""

import json
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional

HERE = Path(__file__).resolve().parent
ICE = HERE.parents[2]
PO = ICE / "benchmarks" / "autoresearch" / "ar5" / "polish"
AR5 = ICE / "benchmarks" / "autoresearch" / "ar5"
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent.parent / "common"))

import doc_tables as DT  # noqa: E402
import polish_judge as PJ  # noqa: E402
import stats as ST  # noqa: E402

SAMPLES = PJ.SAMPLES
FULL = ("held", "mid", "grad", "big")  # samples whose pass is read
N = 24570
STATUS = (
    "study result of step 14 (rule G2b: a second local optimization from a failed rescue's ends "
    "before the grid search); phase A only, held2 after the pass waits for step 13; library key "
    "not written; not a default arm"
)
# H's held2 total with the pass: step 12's tracked stats file (a summary, not a passes_held2 file)
H2_TOTAL = int(
    json.loads((AR5 / "stallhit" / "stats_held2.json").read_text())["arms"]["H"]["stages"]["total"]
)
ORDER_J = "interleaved + rescue G1 + stall check PH(3000) + second local optimization G2b"


def n(x: int) -> str:
    return f"{int(x):,}"


def pct(x: float, d: int = 2) -> str:
    return f"{x:.{d}f}"


def left_cell(d: Dict[str, Any]) -> str:
    """Grid searches left in J, then H's in brackets."""
    j, h = d["rescue"]["J"], d["rescue"]["H"]
    return f"{j['grid_fallbacks_left']} (H {h['grid_fallbacks_left']})"


def wil(k: int, total: int = N) -> str:
    lo, hi = ST.wilson(int(k), total)
    return f"{int(k)} ({100 * lo:.2f}% to {100 * hi:.2f}%)"


def stored_points(st: Dict[str, Any]) -> List[Dict[str, Any]]:
    """Every point of the plot: (sample, arm, after the pass, evaluations, wrong, Wilson lo, hi)."""
    pts: List[Dict[str, Any]] = []

    def add(s: str, arm: str, after: bool, evals: int, k: int, lo: float, hi: float) -> None:
        pts.append(
            dict(sample=s, arm=arm, after=after, evals=int(evals), wrong=int(k), lo=lo, hi=hi)
        )

    for d in ("order", "rescue", "stallhit"):
        for r in json.loads((AR5 / d / "pareto_rows.json").read_text()):
            if r["arm"] in ("A", "G", "H", "S8") and (r["arm"] != "S8" or r["sample"] == "held2"):
                add(
                    r["sample"],
                    r["arm"],
                    bool(r["boundary_pass"]),
                    r["evals_total"],
                    r["wrong"],
                    r["interval_pct"][0] / 100.0,
                    r["interval_pct"][1] / 100.0,
                )
    for s in FULL:  # S8 from the stored part-2 numbers of step 12
        sj = json.loads((AR5 / "stallhit" / f"stats_{s}.json").read_text())["arms"]["H"]
        for after, ek, sk in ((False, "before_pass", "s8_before"), (True, "total", "s8_after")):
            side = sj[sk]
            add(
                s,
                "S8",
                after,
                sj["stages_s8"][ek],
                side["wrong"],
                side["wilson_lo"],
                side["wilson_hi"],
            )
    # J (the held2 after-the-pass point waits)
    for s in SAMPLES:
        d = st[s]
        if s == "held2":
            b = d["before_pass"]["j"]
            add(
                s,
                "J",
                False,
                d["evals_before_pass"]["j"],
                b["wrong"],
                b["wilson_lo"],
                b["wilson_hi"],
            )
            continue
        j = d["arms"]["J"]
        add(
            s,
            "J",
            False,
            j["stages"]["before_pass"],
            j["before"]["wrong"],
            j["before"]["wilson_lo"],
            j["before"]["wilson_hi"],
        )
        add(
            s,
            "J",
            True,
            j["stages"]["total"],
            j["after"]["wrong"],
            j["after"]["wilson_lo"],
            j["after"]["wilson_hi"],
        )
    return pts


def plot(pts: List[Dict[str, Any]]) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.ticker import LogFormatter

    fig, axes = plt.subplots(1, 2, figsize=(11.5, 5.0))
    col = dict(
        held2="tab:purple", held="tab:red", mid="tab:blue", grad="tab:green", big="tab:orange"
    )
    mk = dict(S8="^", A="D", G="o", H="s", J="*")
    for ax, after, title in zip(axes, (False, True), ("before the pass", "after the pass")):
        for p in pts:
            if p["after"] != after:
                continue
            k = p["wrong"]
            lo, hi = N * p["lo"], N * p["hi"]
            ax.errorbar(
                p["evals"] / 1e6,
                k,
                yerr=[[k - lo], [hi - k]],
                fmt=mk[p["arm"]],
                color=col[p["sample"]],
                mfc="none" if p["arm"] == "S8" else col[p["sample"]],
                ms=11 if p["arm"] == "J" else 6,
                capsize=2,
                lw=0.8,
            )
        ax.set_xscale("log")
        ax.xaxis.set_minor_formatter(LogFormatter(labelOnlyBase=False, minor_thresholds=(2, 0.4)))
        ax.tick_params(axis="x", which="minor", labelsize=7)
        ax.set_xlabel("run evaluations (millions, log%s)" % (", with the pass" if after else ""))
        ax.set_ylabel("wrong voxels (more than 1 degree), Wilson 95%")
        ax.set_title(title)
        ax.grid(alpha=0.3)
    axes[1].text(
        0.02,
        0.02,
        "held2 after the pass: waits for step 13",
        transform=axes[1].transAxes,
        fontsize=8,
        color="tab:purple",
    )
    from matplotlib.lines import Line2D

    h = [
        Line2D([], [], marker=m, ls="", color="k", mfc="none" if a == "S8" else "k", label=a)
        for a, m in mk.items()
    ]
    h += [Line2D([], [], marker="o", ls="", color=c, label=s) for s, c in col.items()]
    fig.legend(handles=h, loc="lower center", fontsize=8, ncol=10, frameon=False)
    fig.suptitle(
        "S8 (open triangle), A = step 8 (diamond), G = step 9 (circle), H = step 12 (square), "
        "J = H with a second local optimization (star)",
        fontsize=9,
    )
    fig.tight_layout(rect=(0, 0.07, 1, 0.97))
    fig.savefig(PO / "ar5_polish_evals_wrong.png", dpi=110)


def pareto(st: Dict[str, Any]) -> None:
    rows: List[Dict[str, Any]] = []
    for s in SAMPLES:
        d = st[s]
        for after in (False, True):
            if s == "held2":
                if after:
                    continue  # waits for step 13
                total = d["evals_before_pass"]["j"]
                side = d["before_pass"]["j"]
            else:
                j = d["arms"]["J"]
                total = j["stages"]["total" if after else "before_pass"]
                side = j["after" if after else "before"]
            rows.append(
                dict(
                    sample=s,
                    arm="J",
                    order=ORDER_J,
                    boundary_pass=after,
                    evals_total=int(total),
                    evals_per_voxel=round(total / N, 2),
                    wrong=int(side["wrong"]),
                    n=N,
                    interval_pct=[
                        round(100 * side["wilson_lo"], 2),
                        round(100 * side["wilson_hi"], 2),
                    ],
                    source="icenine_py/benchmarks/autoresearch/ar5/polish/stats_%s.json" % s,
                    status=STATUS,
                )
            )
    (PO / "pareto_rows.json").write_text(json.dumps(rows, indent=1) + "\n")


def inr(x: float, lo: float, hi: float) -> bool:
    return bool(lo <= x <= hi)


def predictions(st: Dict[str, Any], g1: Dict[str, Any], g2: Dict[str, Any]) -> Dict[str, Any]:
    b = g1["forms"][g1["chosen"]]["metrics"]
    a = g1["forms"]["G2a"]["metrics"]
    ang = b["accepted_angle_to_grid_end"]
    p2 = dict(
        accepted=b["accepted"],
        accepted_in_12_to_20=inr(b["accepted"], 12, 20),
        g2b_at_least_g2a=bool(b["accepted"] >= a["accepted"]),
        none_wrong=bool(b["wrong_accepted"] == 0),
        angle_range_deg=[min(ang), max(ang)],
        angle_in_0p002_to_0p2=bool(min(ang) >= 0.002 and max(ang) <= 0.2),
        net_saving_pct=100.0 * b["net_saving"] / 9_595_587,
        net_in_40_to_70=inr(100.0 * b["net_saving"] / 9_595_587, 40, 70),
    )
    m = g2["metrics"]
    net3 = 100.0 * m["net_saving"] / 7_252_187
    p3 = dict(
        accepted=m["accepted"],
        accepted_in_7_to_15=inr(m["accepted"], 7, 15),
        none_wrong=bool(m["wrong_accepted"] == 0),
        net_saving_pct=net3,
        net_in_30_to_70=inr(net3, 30, 70),
    )
    cut = {s: st[s]["rules"]["total_cut_pct"] for s in FULL}
    rng4 = dict(held=(6, 15), mid=(0, 8), grad=(0, 10))
    p4 = {s: dict(cut_pct=cut[s], right=inr(cut[s], *rng4[s])) for s in rng4}
    p4["big"] = dict(cut_pct=cut["big"], right=bool(cut["big"] == 0.0))
    p4["held2"] = dict(
        cut_pct_before_pass=st["held2"]["cut_before_pass_pct"], right="waits (needs the pass)"
    )
    left = dict(
        held2=(st["held2"]["rescue"]["J"]["grid_fallbacks_left"], 6, 14, 21),
        held=(st["held"]["rescue"]["J"]["grid_fallbacks_left"], 5, 12, 19),
    )
    p4["grid_searches_left"] = {
        s: dict(h=v[3], j=v[0], right=inr(v[0], v[1], v[2])) for s, v in left.items()
    }
    p5: Dict[str, Any] = {}
    for s in FULL:
        pa = st[s]["j_vs_h"]["voxels"]["after"]
        net = pa["broken"] - pa["fixed"]
        dm = abs(
            st[s]["arms"]["J"]["after"]["median_err_deg"]
            - st[s]["arms"]["H"]["after"]["median_err_deg"]
        )
        p5[s] = dict(
            broken_minus_fixed=net,
            in_minus2_to_3=inr(net, -2, 3),
            median_move_deg=dm,
            within_0p0005=bool(dm <= 0.0005),
        )
    p5["held2"] = "waits (needs the pass)"
    return dict(p1=True, p2=p2, p3=p3, p4=p4, p5=p5)


def tables(st: Dict[str, Any], pred: Dict[str, Any], g1: Dict[str, Any]) -> Dict[str, str]:
    out: Dict[str, str] = {}
    # main evaluations
    rows = []
    for s in SAMPLES:
        d = st[s]
        if s == "held2":
            rows.append(
                dict(
                    sample=s,
                    j_total="waits (pass)",
                    h_total=n(H2_TOTAL),
                    fewer="waits",
                    j_before=n(d["evals_before_pass"]["j"]),
                    h_before=n(d["evals_before_pass"]["h"]),
                    fewer_before=pct(d["cut_before_pass_pct"]),
                    left=left_cell(d),
                    g2=f"{d['g2_accepted']} of {d['g2_failed_rescues']}",
                )
            )
            continue
        j, h = d["arms"]["J"], d["arms"]["H"]
        rows.append(
            dict(
                sample=s,
                j_total=n(j["stages"]["total"]),
                h_total=n(h["stages"]["total"]),
                fewer=pct(d["rules"]["total_cut_pct"]),
                j_before=n(j["stages"]["before_pass"]),
                h_before=n(h["stages"]["before_pass"]),
                fewer_before=pct(
                    PJ.pct_fewer(j["stages"]["before_pass"], h["stages"]["before_pass"])
                ),
                left=left_cell(d),
                g2=f"{d['g2_accepted']} of {d['g2_failed_rescues']}",
            )
        )
    out["polish_evals"] = DT.markdown_table(
        rows,
        [
            "sample",
            "j_total",
            "h_total",
            "fewer",
            "j_before",
            "h_before",
            "fewer_before",
            "left",
            "g2",
        ],
        labels=None,
    ).replace(
        "| sample | j_total | h_total | fewer | j_before | h_before | fewer_before | left | g2 |",
        "| sample | J total with the pass | H total with the pass | % fewer (J) | "
        "J before the pass | H before the pass | % fewer before the pass | "
        "grid searches left (H) | G2 accepted |",
    )
    # wrong voxels, pairs
    rows = []
    for s in FULL:
        d = st[s]
        j, h = d["arms"]["J"], d["arms"]["H"]
        pa = d["j_vs_h"]["voxels"]["after"]
        rows.append(
            dict(
                sample=s,
                wj=wil(j["after"]["wrong"]),
                wh=wil(h["after"]["wrong"]),
                ws=wil(j["s8_after"]["wrong"]),
                fixed=pa["fixed"],
                broken=pa["broken"],
                both=pa["both_wrong"],
                mj=f"{j['after']['median_err_deg']:.5f}",
                mh=f"{h['after']['median_err_deg']:.5f}",
                seeds=d["seed_steps"]["wrong_seeds"]["broad"],
            )
        )
    rows.append(
        dict(
            sample="held2",
            wj="waits",
            wh="waits",
            ws="waits",
            fixed="waits",
            broken="waits",
            both="waits",
            mj="waits",
            mh="waits",
            seeds=st["held2"]["seed_steps"]["wrong_seeds"]["broad"],
        )
    )
    out["polish_wrong"] = DT.markdown_table(
        rows,
        ["sample", "wj", "wh", "ws", "fixed", "broken", "both", "mj", "mh", "seeds"],
    ).replace(
        "| sample | wj | wh | ws | fixed | broken | both | mj | mh | seeds |",
        "| sample | wrong after the pass, J (Wilson 95%) | wrong, H | wrong, S8 | "
        "fixed by J vs H | broken by J vs H | wrong in both | median error J (degree) | "
        "median error H | extra wrong seeds (matched + unmatched) |",
    )
    # rules
    rows = []
    for s in SAMPLES:
        if s == "held2":
            rows.append(
                dict(
                    sample=s,
                    r1="waits",
                    r2="waits",
                    r3="waits",
                    r4=str(st[s]["seed_steps"]["wrong_seeds"]["broad"] <= 1) + " (seeds only)",
                    r5="waits",
                    passes="waits for step 13",
                )
            )
            continue
        r = st[s]["rules"]
        rows.append(
            dict(
                sample=s,
                r1=f"{r['r1']} ({r['total_cut_pct']:.2f}% fewer)",
                r2=r["r2"],
                r3=r["r3"],
                r4=r["r4"],
                r5=("n/a: H fails it" if r["h_fails_5"] and not r["r5"] else r["r5"]),
                passes=r["passes"],
            )
        )
    out["polish_rules"] = DT.markdown_table(
        rows, ["sample", "r1", "r2", "r3", "r4", "r5", "passes"]
    ).replace(
        "| sample | r1 | r2 | r3 | r4 | r5 | passes |",
        "| sample | (1) total with the pass | (2) broken minus fixed at most 5 | "
        "(3) median error | (4) wrong seeds | (5) house rule against S8 | J passes |",
    )
    # house rule numbers against S8
    rows = []
    for s in FULL:
        j = st[s]["arms"]["J"]
        pa = j["paired_after"]
        rows.append(
            dict(
                sample=s,
                fewer=pct(PJ.pct_fewer(j["stages"]["total"], j["stages_s8"]["total"]), 1),
                fixed=pa["fixed"],
                broken=pa["broken"],
                med=f"{j['after']['median_err_deg'] - j['s8_after']['median_err_deg']:+.5f}",
            )
        )
    out["polish_house"] = DT.markdown_table(
        rows, ["sample", "fewer", "fixed", "broken", "med"]
    ).replace(
        "| sample | fewer | fixed | broken | med |",
        "| sample | % fewer evaluations than S8 (with the pass) | fixed vs S8 | broken vs S8 | "
        "median error minus S8's (degree) |",
    )
    # G2 ends taken
    rows = []
    for s in SAMPLES:
        acc = [r for r in st[s]["g2_cases"] if r["polished"]]
        errs = sorted(r["err_deg"] for r in acc)
        angs = sorted(r["angle_to_h_map_deg"] for r in acc)
        rows.append(
            dict(
                sample=s,
                acc=f"{len(acc)} of {st[s]['g2_failed_rescues']}",
                wrong=sum(1 for e in errs if e > 1.0),
                emax=f"{max(errs):.5f}" if errs else "n/a",
                amed=f"{angs[len(angs) // 2]:.4f}" if angs else "n/a",
                amax=f"{max(angs):.3f}" if angs else "n/a",
                evals=n(sum(r["polish_evals"] for r in st[s]["g2_cases"])),
            )
        )
    out["polish_ends"] = DT.markdown_table(
        rows, ["sample", "acc", "wrong", "emax", "amed", "amax", "evals"]
    ).replace(
        "| sample | acc | wrong | emax | amed | amax | evals |",
        "| sample | G2 accepted (of failed rescues) | accepted ends wrong | "
        "largest error of an accepted end (degree) | median angle to H's map | "
        "largest angle to H's map (degree) | G2 evaluations |",
    )
    # validity
    rows = []
    for s in SAMPLES:
        v = st[s]["validity"]
        rows.append(
            dict(
                sample=s,
                fin=v["j_finished"],
                ali=v["aligned_draws"],
                steps=f"{v['seed_steps_j']} (H {v['seed_steps_h']}); matched {v['matched_steps']}",
                diff=v["fits_differ"],
                grid=f"{v['grid_searches_same_evals']} of {v['grid_searches_both_ran']}",
            )
        )
    out["polish_validity"] = DT.markdown_table(
        rows, ["sample", "fin", "ali", "steps", "diff", "grid"]
    ).replace(
        "| sample | fin | ali | steps | diff | grid |",
        "| sample | J finished | aligned draws | seed steps | "
        "fits that differ in matched steps | "
        "grid searches both arms ran with equal evaluations |",
    )
    # held2: what is allowed now
    d = st["held2"]
    bp = d["before_pass"]
    rows = [
        dict(
            item="run evaluations before the pass, J",
            value=n(d["evals_before_pass"]["j"]),
            status="final",
        ),
        dict(
            item="run evaluations before the pass, H",
            value=n(d["evals_before_pass"]["h"]),
            status="final",
        ),
        dict(item="% fewer before the pass", value=pct(d["cut_before_pass_pct"]), status="final"),
        dict(
            item="grid searches left, J (H)",
            value=left_cell(d),
            status="final",
        ),
        dict(
            item="G2 accepted",
            value=f"{d['g2_accepted']} of {d['g2_failed_rescues']}",
            status="final",
        ),
        dict(
            item="wrong voxels before the pass, J",
            value=wil(bp["j"]["wrong"]),
            status="before the pass only",
        ),
        dict(
            item="wrong voxels before the pass, H",
            value=wil(bp["h"]["wrong"]),
            status="before the pass only",
        ),
        dict(
            item="fixed / broken by J vs H, before the pass",
            value=f"{bp['fixed']} / {bp['broken']}",
            status="before the pass only",
        ),
        dict(
            item="extra wrong seeds (matched + unmatched)",
            value=d["seed_steps"]["wrong_seeds"]["broad"],
            status="final",
        ),
        dict(
            item="total with the pass, wrong after the pass, rules 1, 2, 3, 5, predictions 4 and 5",
            value="not computed",
            status="waits for step 13's held2 result",
        ),
    ]
    out["polish_held2_now"] = DT.markdown_table(rows, ["item", "value", "status"]).replace(
        "| item | value | status |", "| held2 number | value | status |"
    )
    return out


def main() -> None:
    st: Dict[str, Any] = {s: json.loads((PO / f"stats_{s}.json").read_text()) for s in SAMPLES}
    g1 = json.loads((PO / "gate1.json").read_text())
    g2 = json.loads((PO / "gate2.json").read_text())
    pts = stored_points(st)
    plot(pts)  # the plot comes first
    pareto(st)
    pred = predictions(st, g1, g2)
    DT.write_tables(PO / "tables.md", tables(st, pred, g1))
    per = {s: bool(st[s]["rules"]["passes"]) for s in FULL}
    failed = [
        f"{s}:({k[1]})"
        for s in FULL
        for k in ("r1", "r2", "r3", "r4", "r5")
        if not st[s]["rules"][k] and not (k == "r5" and st[s]["rules"]["h_fails_5"])
    ]
    out = dict(
        per_sample_passes_available=per,
        held2=(
            "waits for step 13's held2 result (plan step 7): "
            "after-the-pass numbers and rules not computed"
        ),
        failed_rules=failed,
        verdict=(
            "NOT READY FOR KEY: the four other samples pass; held2's after-the-pass judgement "
            "waits for step 13's held2 result"
            if all(per.values())
            else "STOP: a rule failed on a sample (see failed_rules)"
        ),
        rules={s: st[s]["rules"] for s in FULL},
        predictions=pred,
        wrong_seeds_broad={s: st[s]["seed_steps"]["wrong_seeds"]["broad"] for s in SAMPLES},
        status=STATUS,
    )
    (PO / "summary.json").write_text(json.dumps(out, indent=1, default=float) + "\n")
    print(
        json.dumps(
            {k: out[k] for k in ("per_sample_passes_available", "failed_rules", "verdict")},
            indent=1,
        )
    )
    print(json.dumps(pred, indent=1, default=float))


if __name__ == "__main__":
    main()
