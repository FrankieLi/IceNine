"""AR5 step 27: the plot and the doc tables of phase A (reads ``score_*.json`` of ``px_score``).

  uv run python scripts/autoresearch/ar5/px_summary.py

Writes ``ar5/px/summary.json``, ``ar5/px/tables.md`` (``doc_tables.write_tables``) and the plot
``ar5/px/px_wall.png`` (wall per sample, old and new, indexing and the rest stacked; beside it
paired broken and fixed). Whatever score files exist are used; no truth is read; no p-value.
"""

import json
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional

HERE = Path(__file__).resolve().parent
ICE = HERE.parents[2]
RD = ICE / "benchmarks" / "autoresearch" / "ar5" / "px"
sys.path.insert(0, str(ICE / "scripts" / "common"))
import doc_tables as DT  # noqa: E402


def load(name: str) -> Optional[Dict[str, Any]]:
    p = RD / f"score_{name}.json"
    return json.loads(p.read_text()) if p.exists() else None


def mtable(rows: List[Dict[str, Any]], spec: List[Any]) -> str:
    """Rows keyed by ``spec``'s first items; columns headed by its second items."""
    heads = [h for _k, h, _f in spec]
    r2 = [{h: r.get(k) for k, h, _f in spec} for r in rows]
    return DT.markdown_table(r2, heads)


def pct(x: Optional[float]) -> str:
    return "n/a" if x is None else f"{100 * x:.1f}%"


def tables(q1: Any, t1: Any, live: Any, h1: Any) -> Dict[str, str]:
    out: Dict[str, str] = {}
    if q1:
        rows = []
        for s, r in q1["rows"].items():
            rows.append(
                dict(
                    sample=s,
                    n=r["n"],
                    top20_base=r["top20_base"],
                    top20_new=r["top20_new"],
                    lost=r["lost"],
                    gained=r["gained"],
                    top1_base=r["top1_base"],
                    top1_new=r["top1_new"],
                    pairs_less=pct(1 - r["pairs_new"] / r["pairs_base"]),
                )
            )
        j = q1["judge"]
        rows.append(
            dict(
                sample="all",
                n=sum(r["n"] for r in q1["rows"].values()),
                top20_base=sum(r["top20_base"] for r in q1["rows"].values()),
                top20_new=sum(r["top20_new"] for r in q1["rows"].values()),
                lost=j["lost"],
                gained=j["gained"],
                top1_base=sum(r["top1_base"] for r in q1["rows"].values()),
                top1_new=sum(r["top1_new"] for r in q1["rows"].values()),
                pairs_less="",
            )
        )
        out["px_q1"] = mtable(
            rows,
            [
                ("sample", "sample", None),
                ("n", "voxels", None),
                ("top20_base", "in top 20, base", None),
                ("top20_new", "in top 20, new", None),
                ("lost", "lost", None),
                ("gained", "gained", None),
                ("top1_base", "rank 0, base", None),
                ("top1_new", "rank 0, new", None),
                ("pairs_less", "pairs fewer", None),
            ],
        )
    if t1:
        rows = []
        for s, r in t1["rows"].items():
            rows.append(
                dict(
                    sample=s,
                    n=r["n_voxels"],
                    old=round(r["old"]["s"], 1),
                    new=round(r["new"]["s"], 1),
                    less=f"{r['less_pct_seconds']:.1f}%",
                    pairs=f"{r['less_pct_pairs']:.1f}%",
                    votes=f"{r['less_pct_votes']:.1f}%",
                    label=r["timing_label"],
                )
            )
        p = t1["pooled"]
        rows.append(
            dict(
                sample="pooled",
                n=sum(r["n_voxels"] for r in t1["rows"].values()),
                old=round(p["old_s"], 1),
                new=round(p["new_s"], 1),
                less=pct(t1["judge"]["fewer"]),
                pairs=pct(1 - p["new_pairs"] / p["old_pairs"]),
                votes=pct(1 - p["new_votes"] / p["old_votes"]),
                label="",
            )
        )
        out["px_t1"] = mtable(
            rows,
            [
                ("sample", "sample", None),
                ("n", "voxels", None),
                ("old", "old s", None),
                ("new", "new s", None),
                ("less", "seconds less", None),
                ("pairs", "pairs less", None),
                ("votes", "votes less", None),
                ("label", "machine", None),
            ],
        )
    if live:
        rows = []
        for s, r in live.items():
            t = r["totals"]
            rows.append(
                dict(
                    sample=s,
                    stored=t["stored"]["total"],
                    old=t["old"]["total"],
                    new=t["new"]["total"],
                    more=f"{100 * (t['new']['total'] / t['stored']['total'] - 1):+.2f}%",
                    wrong_stored=r["wrong"]["stored"]["n"],
                    wrong_old=r["wrong"]["old"]["n"],
                    wrong_new=r["wrong"]["new"]["n"],
                    fixed=r["new_vs_stored"]["fixed"],
                    broken=r["new_vs_stored"]["broken"],
                    med=f"{r['new_vs_stored']['median_higher']:+.5f}",
                    grid=r["grid_fallbacks_left"]["new"] - r["grid_fallbacks_left"]["stored"],
                )
            )
        out["px_live"] = mtable(
            rows,
            [
                ("sample", "sample", None),
                ("stored", "evaluations, stored", None),
                ("old", "old rerun", None),
                ("new", "new", None),
                ("more", "new vs stored", None),
                ("wrong_stored", "wrong, stored", None),
                ("wrong_old", "wrong, old rerun", None),
                ("wrong_new", "wrong, new", None),
                ("fixed", "fixed", None),
                ("broken", "broken", None),
                ("med", "median error change (deg)", None),
                ("grid", "grid searches vs stored", None),
            ],
        )
        rows = []
        for s, r in live.items():
            iw = r["index_work"]
            jo, jn = r["job"]["old"], r["job"]["new"]
            rows.append(
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
        out["px_wall"] = mtable(
            rows,
            [
                ("sample", "sample", None),
                ("pairs_old", "pairs, old", None),
                ("pairs_new", "pairs, new", None),
                ("votes_old", "votes, old", None),
                ("votes_new", "votes, new", None),
                ("idx_old", "indexing s, old", None),
                ("idx_new", "indexing s, new", None),
                ("wall_old", "wall s, old", None),
                ("wall_new", "wall s, new", None),
                ("less", "wall less", None),
                ("lab", "preflight verdict old / new", None),
                ("load", "mean load old / new", None),
            ],
        )
    if q1 and t1 and live:
        out["px_gates"] = mtable(
            gate_rows(q1, t1, live),
            [
                ("gate", "gate", None),
                ("what", "measured", None),
                ("limit", "limit", None),
                ("verdict", "verdict", None),
            ],
        )
    if h1:
        r = h1
        rows = [
            dict(
                s8=r["totals"]["s8"],
                new=r["totals"]["j5sd_new"],
                fewer=pct(r["fewer"]),
                wrong_s8=r["wrong"]["s8"]["n"],
                wrong_new=r["wrong"]["j5sd_new"]["n"],
                fixed=r["new_vs_s8"]["fixed"],
                broken=r["new_vs_s8"]["broken"],
                med=f"{r['new_vs_s8']['median_higher']:+.5f}",
                ok=str(r["h1"]["ok"]),
            )
        ]
        out["px_h1"] = mtable(
            rows,
            [
                ("s8", "S8 evaluations", None),
                ("new", "new J5SD evaluations", None),
                ("fewer", "fewer", None),
                ("wrong_s8", "wrong, S8", None),
                ("wrong_new", "wrong, new", None),
                ("fixed", "fixed", None),
                ("broken", "broken", None),
                ("med", "median error change (deg)", None),
                ("ok", "house rule met", None),
            ],
        )
    return out


def gate_rows(q1: Any, t1: Any, live: Any) -> List[Dict[str, Any]]:
    """The frozen gates against the measured numbers (verdicts from the frozen judges)."""
    import px_rules as PR

    rows: List[Dict[str, Any]] = []
    j = q1["judge"]
    rows.append(
        dict(
            gate="q1",
            what=(
                f"{j['lost']} of {j['was_in_top20']} voxels lost from the top 20, "
                f"{j['gained']} gained"
            ),
            limit="at most 0.5% lost",
            verdict="met" if j["ok"] else "not met",
        )
    )
    npairs = sum(r["index_work"]["new"]["total"]["pairs"] for r in live.values())
    nvotes = sum(r["index_work"]["new"]["total"]["votes"] for r in live.values())
    p = t1["pooled"]
    c1 = PR.judge_c1(npairs, nvotes, p["old_pairs"], p["old_votes"])
    rows.append(
        dict(
            gate="c1",
            what=(
                f"pairs {pct(c1['pairs_fewer'])} fewer, votes {pct(c1['votes_fewer'])} fewer "
                "(new live runs against the stored seed voxels indexed again)"
            ),
            limit="at least 60% fewer",
            verdict="met" if c1["ok"] else "not met",
        )
    )
    tot = sum(r["totals"]["new"]["total"] for r in live.values())
    c2 = PR.judge_c2(tot)
    rows.append(
        dict(
            gate="c2",
            what=f"{c2['new']:,} against {c2['stored']:,} stored ({100 * c2['more']:+.2f}%)",
            limit="at most 2% more",
            verdict="met" if c2["ok"] else "not met",
        )
    )
    for s, r in live.items():
        c3 = r["c3"]
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
    rows.append(
        dict(
            gate="t1",
            what=f"indexing seconds {pct(t1['judge']['fewer'])} lower (contended)",
            limit="at least 50% lower",
            verdict="met" if t1["judge"]["ok"] else "not met",
        )
    )
    ow = sum(r["job"]["old"]["wall_s"] for r in live.values())
    nw = sum(r["job"]["new"]["wall_s"] for r in live.values())
    t2 = PR.judge_t(PR.fewer(nw, ow), 0.25)
    rows.append(
        dict(
            gate="t2",
            what=f"pooled wall {pct(t2['fewer'])} lower ({ow} s to {nw} s; contended)",
            limit="at least 25% lower",
            verdict="met" if t2["ok"] else "not met",
        )
    )
    rows.append(
        dict(
            gate="h1", what="held7 not built (c2 and c3 not met)", limit="", verdict="not evaluated"
        )
    )
    return rows


def plot(live: Any, h1: Any) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(1, 2, figsize=(10, 4))
    if live:
        xs, labels = [], []
        for i, (s, r) in enumerate(live.items()):
            for j, a in enumerate(("old", "new")):
                x = i * 3 + j
                wall = r["job"][a]["wall_s"] or 0
                idx = r["index_work"][a]["total"]["seconds"]
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
        ax[0].set_xticklabels(labels)
        ax[0].set_ylabel("minutes of wall (contended)")
        ax[0].set_title("Whole J5SD run: old (600,000 pairs) and new (50,000)")
        ax[0].legend()
        k = 0
        for s, r in live.items():
            v = r["new_vs_stored"]
            ax[1].bar(k, v["broken"], color="tab:red", label="broken" if k == 0 else None)
            ax[1].bar(k + 0.4, v["fixed"], color="tab:green", label="fixed" if k == 0 else None)
            k += 1.5
        ticks = [0.2 + 1.5 * i for i in range(len(live))]
        names = [f"{s}\nvs stored" for s in live]
        if h1:
            v = h1["new_vs_s8"]
            ax[1].bar(k, v["broken"], color="tab:red")
            ax[1].bar(k + 0.4, v["fixed"], color="tab:green")
            ticks.append(k + 0.2)
            names.append("held7\nvs S8")
        ax[1].set_xticks(ticks)
        ax[1].set_xticklabels(names)
        ax[1].set_ylabel("voxels")
        ax[1].set_title("Paired by voxel (new against the baseline)")
        ax[1].legend()
    fig.tight_layout()
    fig.savefig(RD / "px_wall.png", dpi=110)


def main() -> None:
    q1, t1, live, h1 = load("q1"), load("t1"), load("live"), load("h1")
    tabs = tables(q1, t1, live, h1)
    DT.write_tables(RD / "tables.md", tabs)
    extra = json.loads((RD / "extra.json").read_text()) if (RD / "extra.json").exists() else None
    split = (
        json.loads((RD / "stage_split.json").read_text())
        if (RD / "stage_split.json").exists()
        else None
    )
    (RD / "summary.json").write_text(
        json.dumps(dict(q1=q1, t1=t1, live=live, h1=h1, extra=extra, stage_split=split), indent=1)
        + "\n"
    )
    plot(live, h1)
    print("wrote", sorted(tabs))


if __name__ == "__main__":
    main()
