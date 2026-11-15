"""AR5 step 16, phase A step 6: the plot, the tables, the predictions, the pareto rows, the verdict.

  uv run python scripts/autoresearch/ar5/high500_summary.py

Reads ``ar5/high500/stats_<s>.json`` (``high500_stats.py``, five samples) and ``gate.json``. Writes
``ar5_high500_evals_wrong.png`` (run evaluations against wrong voxels after the pass: S8, J, J5,
one panel per sample), ``tables.md`` (blocks for ``sync_doc_tables.py``), ``summary.json`` and
``pareto_rows.json`` (the format of the earlier steps, with a status field). No p-values; wrong
voxels carry a Wilson 95% interval; fixed and broken are counted by voxel.
"""

import json
import sys
from pathlib import Path
from typing import Any, Dict, List

import high500_judge as HJ
import high500_rules as HR

HERE = Path(__file__).resolve().parent
ICE = HERE.parents[2]
HD = ICE / "benchmarks" / "autoresearch" / "ar5" / "high500"
sys.path.insert(0, str(ICE / "scripts" / "common"))
import doc_tables as DT  # noqa: E402
import stats as ST  # noqa: E402

N = 24570
SAMPLES = HJ.SAMPLES


def load() -> Dict[str, Dict[str, Any]]:
    return {s: json.loads((HD / f"stats_{s}.json").read_text()) for s in SAMPLES}


def plot(st: Dict[str, Dict[str, Any]]) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 5, figsize=(17, 3.6))
    col = dict(S8="#888888", J="#1f77b4", J5="#d62728")
    for ax, s in zip(axes, SAMPLES):
        for a in ("S8", "J", "J5"):
            x = st[s]["arms"][a]
            k = x["after"]["wrong"]
            lo, hi = ST.wilson(k, N)
            ev = x["stages"]["total"] / 1e6
            ax.errorbar(
                ev,
                k,
                yerr=[[k - lo * N], [hi * N - k]],
                fmt="o",
                color=col[a],
                label=a,
                capsize=3,
            )
        ax.set_title(s)
        ax.set_xlabel("evaluations with the pass (millions)")
        ax.grid(alpha=0.3)
    axes[0].set_ylabel("wrong voxels after the pass (of 24,570)")
    axes[0].legend()
    fig.tight_layout()
    fig.savefig(HD / "ar5_high500_evals_wrong.png", dpi=110)


def n(x: int) -> str:
    return f"{x:,}"


def wil(k: int) -> str:
    lo, hi = ST.wilson(k, N)
    return f"{k} ({100 * lo:.2f} to {100 * hi:.2f}%)"


def tables(st: Dict[str, Dict[str, Any]]) -> Dict[str, str]:
    t: Dict[str, str] = {}
    rows = []
    for s in SAMPLES:
        for a in ("S8", "J", "J5"):
            x = st[s]["arms"][a]["stages"]
            sp = st[s]["stages_split"].get(a)
            rows.append(
                dict(
                    sample=s,
                    arm=a,
                    seed_step=n(x["seed_step"]),
                    grid=n(sp["grid_fallback"] if sp else x["grid_fallback"]),
                    rescue_g2=n(sp["rescue"] + sp["polish"]) if sp else "0",
                    neighbour=n(x["neighbour"]),
                    revisit=n(x["revisit"]),
                    pass_=n(x["pass_"]),
                    total=n(x["total"]),
                )
            )
    t["high500_stages"] = DT.markdown_table(
        rows,
        [
            "sample",
            "arm",
            "seed_step",
            "grid",
            "rescue_g2",
            "neighbour",
            "revisit",
            "pass_",
            "total",
        ],
    )
    rows = []
    for s in SAMPLES:
        d = st[s]
        j, j5 = d["arms"]["J"], d["arms"]["J5"]
        s8t = d["arms"]["S8"]["stages"]["total"]
        rows.append(
            dict(
                sample=s,
                fewer_vs_j=f"{d['rules']['total_cut_pct']:.2f}%",
                fewer_vs_s8=f"{HJ.pct_fewer(j5['stages']['total'], s8t):.1f}%",
                saved=n(j["stages"]["total"] - j5["stages"]["total"]),
            )
        )
    t["high500_saving"] = DT.markdown_table(rows, ["sample", "fewer_vs_j", "saved", "fewer_vs_s8"])
    rows = []
    for s in SAMPLES:
        m = st[s]["mc"]
        rows.append(
            dict(
                sample=s,
                j_mc=n(m["J"]["mc_evals"]),
                j5_mc=n(m["J5"]["mc_evals"]),
                j_stops=m["J"]["stops"],
                j5_stops=m["J5"]["stops"],
                j_checks=f"{m['J']['checks']} ({m['J']['checks_went_on']} went on)",
                j5_checks=f"{m['J5']['checks']} ({m['J5']['checks_went_on']} went on)",
            )
        )
    t["high500_checks"] = DT.markdown_table(
        rows, ["sample", "j_mc", "j5_mc", "j_stops", "j5_stops", "j_checks", "j5_checks"]
    )
    rows = []
    for s in SAMPLES:
        d = st[s]
        a = d["arms"]
        rows.append(
            dict(
                sample=s,
                s8=wil(a["S8"]["after"]["wrong"]),
                j=wil(a["J"]["after"]["wrong"]),
                j5=wil(a["J5"]["after"]["wrong"]),
                vs_j=f"{d['j5_vs_j']['after']['fixed']} / {d['j5_vs_j']['after']['broken']}",
                vs_s8=f"{d['j5_vs_s8']['after']['fixed']} / {d['j5_vs_s8']['after']['broken']}",
                med_j=f"{a['J']['after']['median_err_deg']:.5f}",
                med_j5=f"{a['J5']['after']['median_err_deg']:.5f}",
            )
        )
    t["high500_wrong"] = DT.markdown_table(
        rows, ["sample", "s8", "j", "j5", "vs_j", "vs_s8", "med_j", "med_j5"]
    )
    rows = []
    for s in SAMPLES:
        d = st[s]
        ws, um = d["seed_steps"]["wrong_seeds"], d["seed_steps"]["unmatched"]
        rows.append(
            dict(
                sample=s,
                matched=d["seed_steps"]["matched_steps"],
                taken_changed=d["seed_steps"]["taken_changed"],
                wrong_j=ws["j"],
                wrong_j5=ws["j5"],
                extra_matched=ws["j5_wrong_where_j_right"],
                extra_unmatched=ws["unmatched_j5_wrong_where_j_map_right"],
                only_j5=um["only_j5_n"],
                only_j=um["only_j_n"],
                grid=f"{d['arms']['J']['grid_searches']} -> {d['arms']['J5']['grid_searches']}",
            )
        )
    t["high500_seeds"] = DT.markdown_table(
        rows,
        [
            "sample",
            "matched",
            "taken_changed",
            "wrong_j",
            "wrong_j5",
            "extra_matched",
            "extra_unmatched",
            "only_j5",
            "only_j",
            "grid",
        ],
    )
    rows = []
    for s in SAMPLES:
        r = st[s]["rules"]
        rows.append(
            dict(
                sample=s,
                r1=r["r1"],
                r2=r["r2"],
                r3=r["r3"],
                r4=r["r4"],
                r5=r["r5"],
                j_fails_5=r["j_fails_5"],
                passes=r["passes"],
            )
        )
    t["high500_rules"] = DT.markdown_table(
        rows, ["sample", "r1", "r2", "r3", "r4", "r5", "j_fails_5", "passes"]
    )
    rows = []
    for s in SAMPLES:
        h = st[s]["house_j5"]
        p = st[s]["j5_vs_s8"]["after"]
        e5, e8 = st[s]["arms"]["J5"]["stages"]["total"], st[s]["arms"]["S8"]["stages"]["total"]
        m5 = st[s]["arms"]["J5"]["after"]["median_err_deg"]
        m8 = st[s]["arms"]["S8"]["after"]["median_err_deg"]
        rows.append(
            dict(
                sample=s,
                fewer=f"{HJ.pct_fewer(e5, e8):.1f}%",
                net=p["broken"] - p["fixed"],
                med_gap=f"{m5 - m8:+.5f}",
                a=h["a"],
                b=h["b"],
                c=h["c"],
            )
        )
    t["high500_house"] = DT.markdown_table(
        rows, ["sample", "fewer", "net", "med_gap", "a", "b", "c"]
    )
    return t


def predictions(st: Dict[str, Dict[str, Any]]) -> Dict[str, Any]:
    rng = dict(held3=(5, 8), held=(5.5, 7), mid=(5.5, 7.5), grad=(1, 2), big=(1.5, 3))
    p3 = {s: dict(value=st[s]["rules"]["total_cut_pct"], range=rng[s]) for s in SAMPLES}
    for s in SAMPLES:
        p3[s]["right"] = bool(rng[s][0] <= p3[s]["value"] <= rng[s][1])
    p4 = {}
    for s in SAMPLES:
        d = st[s]
        net = d["j5_vs_j"]["after"]["broken"] - d["j5_vs_j"]["after"]["fixed"]
        dm = abs(
            d["arms"]["J5"]["after"]["median_err_deg"] - d["arms"]["J"]["after"]["median_err_deg"]
        )
        p4[s] = dict(
            net=net,
            median_move_deg=dm,
            wrong_seed_extra=d["seed_steps"]["wrong_seeds"]["broad"],
            grid_equal=d["arms"]["J5"]["grid_searches"] == d["arms"]["J"]["grid_searches"],
        )
        p4[s]["right"] = bool(
            -2 <= net <= 2
            and dm <= 0.0005
            and p4[s]["wrong_seed_extra"] == 0
            and p4[s]["grid_equal"]
        )
    mid = st["mid"]["j5_vs_s8"]["after"]
    p5 = dict(
        mid_net_vs_s8=mid["broken"] - mid["fixed"],
        mid_house_passes=st["mid"]["house_j5"]["passes"],
        right=bool(st["mid"]["house_j5"]["passes"]),
    )
    return dict(p3=p3, p4=p4, p5=p5)


def pareto(st: Dict[str, Dict[str, Any]]) -> None:
    rows: List[Dict[str, Any]] = []
    order = {
        "S8": "stock seed step",
        "J": "interleaved + rescue G1 + stall check PH(3000) + second local optimization G2b",
        "J5": "interleaved + rescue G1 + stall check PH(500) + second local optimization G2b",
    }
    status = (
        "study result of step 16 (rule PH(500): the high stop alone at a 500-evaluation stall in "
        "the seed step's MC refinement); phase A only, library key not written (step 12's key "
        "SeedIndexMCStallHit already takes k = 500); not a default arm"
    )
    for s in SAMPLES:
        for a in ("S8", "J", "J5"):
            x = st[s]["arms"][a]
            for after in (False, True):
                ev = x["stages"]["total"] if after else x["stages"]["before_pass"]
                k = x["after" if after else "before"]["wrong"]
                lo, hi = ST.wilson(k, N)
                rows.append(
                    dict(
                        sample=s,
                        arm=a,
                        order=order[a],
                        boundary_pass=after,
                        evals_total=int(ev),
                        evals_per_voxel=round(ev / N, 2),
                        wrong=int(k),
                        n=N,
                        interval_pct=[round(100 * lo, 2), round(100 * hi, 2)],
                        source=f"icenine_py/benchmarks/autoresearch/ar5/high500/stats_{s}.json",
                        status=(
                            status if a == "J5" else f"stored or rerun baseline of step 16 ({a})"
                        ),
                    )
                )
    (HD / "pareto_rows.json").write_text(json.dumps(rows, indent=1) + "\n")


def main() -> None:
    st = load()
    plot(st)
    t = tables(st)
    DT.write_tables(HD / "tables.md", t)
    pred = predictions(st)
    verdict = {s: st[s]["rules"]["passes"] for s in SAMPLES}
    summary = dict(
        rule=HR.RULE,
        verdict=verdict,
        step_passes=HJ.step_verdict(verdict),
        predictions=pred,
        validity={s: st[s]["validity"] for s in SAMPLES},
        rules={s: st[s]["rules"] for s in SAMPLES},
    )
    (HD / "summary.json").write_text(json.dumps(summary, indent=1, default=float) + "\n")
    pareto(st)
    print(json.dumps(dict(verdict=verdict, step=summary["step_passes"]), indent=1))
    print(json.dumps(pred, indent=1, default=float)[:2500])


if __name__ == "__main__":
    main()
