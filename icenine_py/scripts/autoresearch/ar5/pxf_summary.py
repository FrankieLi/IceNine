"""AR5 step 28: the plot and the doc tables of phase A (reads ``score_*.json`` of ``pxf_score``).

  uv run python scripts/autoresearch/ar5/pxf_summary.py

Writes ``ar5/pxf/summary.json``, ``ar5/pxf/tables.md`` (``doc_tables.write_tables``) and the plot
``ar5/pxf/pxf_wall.png`` (wall per sample: old, step 27's ``PX`` and ``PXF``, indexing and the rest
stacked; beside it grid searches, and broken and fixed). Whatever score files exist are used; no
truth is read; no p-value. ``mid`` and ``held`` are in-sample for this rule.
"""

import json
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional

HERE = Path(__file__).resolve().parent
ICE = HERE.parents[2]
RD = ICE / "benchmarks" / "autoresearch" / "ar5" / "pxf"
PXD = ICE / "benchmarks" / "autoresearch" / "ar5" / "px"
sys.path.insert(0, str(ICE / "scripts" / "common"))
import doc_tables as DT  # noqa: E402


def load(name: str) -> Optional[Dict[str, Any]]:
    p = RD / f"score_{name}.json"
    return json.loads(p.read_text()) if p.exists() else None


def mtable(rows: List[Dict[str, Any]], spec: List[Any]) -> str:
    heads = [h for _k, h in spec]
    return DT.markdown_table([{h: r.get(k) for k, h in spec} for r in rows], heads)


def pct(x: Optional[float]) -> str:
    return "n/a" if x is None else f"{100 * x:.1f}%"


def tables(live: Any, h1: Any) -> Dict[str, str]:
    out: Dict[str, str] = {}
    if live:
        rows, rows2, rows3 = [], [], []
        for s, r in live["rows"].items():
            t = r["totals"]
            rows.append(
                dict(
                    sample=s,
                    stored=t["stored"]["total"],
                    new=t["new"]["total"],
                    more=f"{100 * (t['new']['total'] / t['stored']['total'] - 1):+.2f}%",
                    wrong_stored=r["wrong"]["stored"]["n"],
                    wrong_new=r["wrong"]["new"]["n"],
                    fixed=r["new_vs_stored"]["fixed"],
                    broken=r["new_vs_stored"]["broken"],
                    med=f"{r['new_vs_stored']['median_higher']:+.5f}",
                    grid_s=r["grid_fallbacks_left"]["stored"],
                    grid_n=r["grid_fallbacks_left"]["new"],
                )
            )
            c = r["steps"]
            rows2.append(
                dict(
                    sample=s,
                    steps=c["seed_steps"],
                    failed=c["failed_first"],
                    fb=c["fallbacks"],
                    fb_ok=c["fallbacks_passed"],
                    rescue=c["to_rescue"],
                    fb_evals=c["fallback_step_evals"],
                    fb_idx=round(c["fallback_index_s"]),
                    eq=f"{r['replay']['equal']} of {r['replay']['seeded_once_in_old']}",
                )
            )
            iw, jo, jn = r["index_work"], r["job"]["old"], r["job"]["new"]
            rows3.append(
                dict(
                    sample=s,
                    pairs_old=iw["old"]["total"]["pairs"],
                    pairs_new=iw["new"]["total"]["pairs"],
                    votes_old=iw["old"]["total"]["votes"],
                    votes_new=iw["new"]["total"]["votes"],
                    idx_old=round(iw["old"]["total"]["seconds"]),
                    idx_new=round(iw["new"]["total"]["seconds"]),
                    wall_old=jo["wall_s"],
                    wall_new=jn["wall_s"],
                    less=pct(r["wall_fewer"]),
                    lab=f"{jo['label']} / {jn['label']}",
                    load=f"{jo['load_mean']:.1f} / {jn['load_mean']:.1f}",
                )
            )
        out["pxf_live"] = mtable(
            rows,
            [
                ("sample", "sample"),
                ("stored", "evaluations, stored"),
                ("new", "evaluations, PXF"),
                ("more", "PXF vs stored"),
                ("wrong_stored", "wrong, stored"),
                ("wrong_new", "wrong, PXF"),
                ("fixed", "fixed"),
                ("broken", "broken"),
                ("med", "median error change (deg)"),
                ("grid_s", "grid searches, stored"),
                ("grid_n", "grid searches, PXF"),
            ],
        )
        out["pxf_steps"] = mtable(
            rows2,
            [
                ("sample", "sample"),
                ("steps", "seed steps"),
                ("failed", "failed the seed test at 50,000"),
                ("fb", "second tries"),
                ("fb_ok", "second try passed"),
                ("rescue", "sent to the rescue"),
                ("fb_evals", "evaluations of second tries"),
                ("fb_idx", "indexing seconds of second tries"),
                ("eq", "second tries equal to the stored old step"),
            ],
        )
        out["pxf_wall"] = mtable(
            rows3,
            [
                ("sample", "sample"),
                ("pairs_old", "pairs, old"),
                ("pairs_new", "pairs, PXF"),
                ("votes_old", "votes, old"),
                ("votes_new", "votes, PXF"),
                ("idx_old", "indexing s, old"),
                ("idx_new", "indexing s, PXF"),
                ("wall_old", "wall s, old"),
                ("wall_new", "wall s, PXF"),
                ("less", "wall less"),
                ("lab", "preflight verdict old / PXF"),
                ("load", "mean load old / PXF"),
            ],
        )
        out["pxf_gates"] = mtable(
            gate_rows(live, h1),
            [("gate", "gate"), ("what", "measured"), ("limit", "limit"), ("verdict", "verdict")],
        )
    if h1:
        r = h1
        rows = []
        for k, name in (("s8", "S8"), ("pxf", "new J5SD (PXF)"), ("old", "plain J5SD")):
            d: Dict[str, Any] = dict(
                arm=name,
                evals=r["totals"][k],
                wrong=r["wrong"][k]["n"],
            )
            if k == "s8":
                d.update(
                    fewer="", fixed="", broken="", med="", grid="", wall=r["job"]["s8"]["wall_s"]
                )
            else:
                pv = r["pxf_vs_s8"] if k == "pxf" else r["old_vs_s8"]
                d.update(
                    fewer=pct(r["fewer"] if k == "pxf" else r["fewer_old"]),
                    fixed=pv["fixed"],
                    broken=pv["broken"],
                    med=f"{pv['median_higher']:+.5f}",
                    grid=r["grid_fallbacks_left"][k],
                    wall=r["job"][k]["wall_s"],
                )
            d["lab"] = r["job"][k]["label"]
            d["load"] = f"{r['job'][k]['load_mean']:.1f}" if r["job"][k]["load_mean"] else ""
            rows.append(d)
        out["pxf_h1"] = mtable(
            rows,
            [
                ("arm", "arm"),
                ("evals", "evaluations"),
                ("fewer", "fewer than S8"),
                ("wrong", "wrong voxels"),
                ("fixed", "fixed vs S8"),
                ("broken", "broken vs S8"),
                ("med", "median error change vs S8 (deg)"),
                ("grid", "grid searches"),
                ("wall", "wall s"),
                ("lab", "preflight verdict"),
                ("load", "mean load"),
            ],
        )
    return out


def gate_rows(live: Any, h1: Any) -> List[Dict[str, Any]]:
    j = live["judge"]
    rows: List[Dict[str, Any]] = []
    c1 = j["c1"]
    rows.append(
        dict(
            gate="c1",
            what=(
                f"pairs {pct(c1['pairs_fewer'])} fewer, votes {pct(c1['votes_fewer'])} fewer "
                "(PXF live runs, second tries included, against the stored seed voxels indexed "
                "again)"
            ),
            limit="at least 50% fewer",
            verdict="met" if c1["ok"] else "not met",
        )
    )
    c2 = j["c2"]
    rows.append(
        dict(
            gate="c2",
            what=f"{c2['new']:,} against {c2['stored']:,} stored ({100 * c2['more']:+.2f}%)",
            limit="at most 2% more",
            verdict="met" if c2["ok"] else "not met",
        )
    )
    for s, c3 in j["c3"].items():
        rows.append(
            dict(
                gate=f"c3 {s}",
                what=(
                    f"broken minus fixed {c3['net']}, median error change "
                    f"{c3['median_higher']:+.5f} deg, grid searches {c3['grid_extra']:+d}"
                ),
                limit="at most 3; at most 0.001 deg; at most 2 more",
                verdict="met" if c3["ok"] else "not met",
            )
        )
    t2 = j["t2"]
    rows.append(
        dict(
            gate="t2",
            what=(
                f"pooled wall {pct(t2['fewer'])} lower ({t2['old_wall_s']} s to "
                f"{t2['new_wall_s']} s; contended, unpaired across sessions; see wall table)"
            ),
            limit="at least 25% lower",
            verdict="met" if t2["ok"] else "not met",
        )
    )
    if h1:
        hh = h1["h1"]
        rows.append(
            dict(
                gate="h1",
                what=(
                    f"held7: {pct(hh['fewer'])} fewer evaluations than S8, broken minus fixed "
                    f"{hh['net']}, median error change {hh['median_higher']:+.5f} deg"
                ),
                limit="at least 20% fewer; at most 5; at most 0.005 deg",
                verdict="met" if hh["ok"] else "not met",
            )
        )
    else:
        rows.append(
            dict(gate="h1", what="held7 not built or not run", limit="", verdict="not evaluated")
        )
    return rows


def plot(live: Any, h1: Any) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    px = json.loads((PXD / "score_live.json").read_text())
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
    if live:
        xs, labels = [], []
        for i, (s, r) in enumerate(live["rows"].items()):
            arms = [
                ("old", r["job"]["old"]["wall_s"], r["index_work"]["old"]["total"]["seconds"]),
                (
                    "PX",
                    px[s]["job"]["new"]["wall_s"],
                    px[s]["index_work"]["new"]["total"]["seconds"],
                ),
                ("PXF", r["job"]["new"]["wall_s"], r["index_work"]["new"]["total"]["seconds"]),
            ]
            for j, (a, wall, idx) in enumerate(arms):
                x = i * 4 + j
                ax[0].bar(x, idx / 60, color="tab:orange", label="indexing" if x == 0 else None)
                ax[0].bar(
                    x,
                    (wall - idx) / 60,
                    bottom=idx / 60,
                    color="tab:blue",
                    label="evaluations and the rest" if x == 0 else None,
                )
                xs.append(x)
                labels.append(f"{s}\n{a}")
        ax[0].set_xticks(xs)
        ax[0].set_xticklabels(labels, fontsize=8)
        ax[0].set_ylabel("minutes of wall (all contended; runs from different sessions, unpaired)")
        ax[0].set_title("Whole J5SD run: old, PX (50,000), PXF (then 600,000)")
        ax[0].legend()
        k = 0.0
        ticks, names = [], []
        for s, r in live["rows"].items():
            v = r["new_vs_stored"]
            ax[1].bar(k, r["grid_fallbacks_left"]["stored"], color="gray", width=0.35,
                      label="grid searches, stored" if k == 0 else None)  # fmt: skip
            ax[1].bar(k + 0.4, r["grid_fallbacks_left"]["new"], color="tab:purple", width=0.35,
                      label="grid searches, PXF" if k == 0 else None)  # fmt: skip
            ax[1].bar(k + 0.8, v["broken"], color="tab:red", width=0.35,
                      label="broken" if k == 0 else None)  # fmt: skip
            ax[1].bar(k + 1.2, v["fixed"], color="tab:green", width=0.35,
                      label="fixed" if k == 0 else None)  # fmt: skip
            ticks.append(k + 0.6)
            names.append(f"{s} vs stored")
            k += 2.2
        if h1:
            v = h1["pxf_vs_s8"]
            ax[1].bar(k + 0.4, h1["grid_fallbacks_left"]["pxf"], color="tab:purple", width=0.35)
            ax[1].bar(k + 0.8, v["broken"], color="tab:red", width=0.35)
            ax[1].bar(k + 1.2, v["fixed"], color="tab:green", width=0.35)
            ticks.append(k + 0.6)
            names.append("held7 vs S8")
        ax[1].set_xticks(ticks)
        ax[1].set_xticklabels(names, fontsize=8)
        ax[1].set_ylabel("voxels or searches")
        ax[1].set_title("Grid searches; paired by voxel (mid, held in-sample)")
        ax[1].legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(RD / "pxf_wall.png", dpi=110)


def derived(live: Any, h1: Any) -> Dict[str, Any]:
    """Numbers quoted in the report, from the score files only (review fixes 4 and 5).

    ``without_grid``: evaluations without the grid-search stage; c2's gain rests on fewer grid
    searches (small counts, route dependent). ``held7_new_vs_plain``: one sample, route effects;
    not evidence of an accuracy gain. t2 and every wall: contended, unpaired across sessions."""
    out: Dict[str, Any] = {}
    if live:
        tot_s = tot_n = g_s = g_n = 0
        per = {}
        for s, r in live["rows"].items():
            st = r["stages"]
            a, b = r["totals"]["stored"]["total"], r["totals"]["new"]["total"]
            ga, gb = st["stored"]["grid"], st["new"]["grid"]
            per[s] = dict(
                stored_no_grid=a - ga,
                new_no_grid=b - gb,
                more_no_grid=(b - gb) / (a - ga) - 1,
                index_more=st["new"]["index"] / st["stored"]["index"] - 1,
                grid_stored=ga,
                grid_new=gb,
            )
            tot_s += a - ga
            tot_n += b - gb
            g_s += ga
            g_n += gb
        out["without_grid"] = dict(
            per_sample=per,
            stored=tot_s,
            new=tot_n,
            more=tot_n / tot_s - 1,
            grid_evals_stored=g_s,
            grid_evals_new=g_n,
        )
    if h1:
        t, st, j = h1["totals"], h1["stages"], h1["job"]
        out["held7_new_vs_plain"] = dict(
            wall_new_s=j["pxf"]["wall_s"],
            wall_plain_s=j["old"]["wall_s"],
            wall_fewer=1 - j["pxf"]["wall_s"] / j["old"]["wall_s"],
            evals_new=t["pxf"],
            evals_plain=t["old"],
            evals_fewer=1 - t["pxf"] / t["old"],
            evals_difference=t["old"] - t["pxf"],
            grid_stage_difference=st["old"]["grid"] - st["pxf"]["grid"],
            index_stage_more=st["pxf"]["index"] / st["old"]["index"] - 1,
            fixed=h1["pxf_vs_old"]["fixed"],
            broken=h1["pxf_vs_old"]["broken"],
            note="one sample, route effects; not evidence that the rule improves accuracy",
        )
    out["timing_note"] = (
        "every wall is contended; t2 compares runs of different sessions (unpaired)"
    )
    return out


def main() -> None:
    live, h1 = load("live"), load("h1")
    tabs = tables(live, h1)
    DT.write_tables(RD / "tables.md", tabs)
    (RD / "summary.json").write_text(
        json.dumps(dict(live=live, h1=h1, derived=derived(live, h1)), indent=1) + "\n"
    )
    plot(live, h1)
    print("wrote", sorted(tabs))


if __name__ == "__main__":
    main()
