"""AR5 step 29: the plot and the doc tables of phase A (reads ``score_*.json`` of ``ru_score``).

  uv run python scripts/autoresearch/ar5/ru_summary.py

Writes ``ar5/ru/summary.json``, ``ar5/ru/tables.md`` (``doc_tables.write_tables``) and the plot
``ar5/ru/ru_wall.png`` (left: wall per sample, step 28's PXF and RU, with indexing, tries and the
rest stacked; right: broken and fixed voxels of RU against PXF). Whatever score files exist are
used; no truth is read; no p-value. ``mid`` and ``held`` are in-sample for this rule.
"""

import json
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional

HERE = Path(__file__).resolve().parent
ICE = HERE.parents[2]
RD = ICE / "benchmarks" / "autoresearch" / "ar5" / "ru"
sys.path.insert(0, str(ICE / "scripts" / "common"))
import doc_tables as DT  # noqa: E402


def load(name: str) -> Optional[Dict[str, Any]]:
    p = RD / f"score_{name}.json"
    return json.loads(p.read_text()) if p.exists() else None


def pct(x: Optional[float]) -> str:
    return "n/a" if x is None else f"{100 * x:.1f}%"


def tables(live: Any, h8: Any) -> Dict[str, str]:
    out: Dict[str, str] = {}
    if live:
        rows, tr, wl = [], [], []
        for s, r in live["rows"].items():
            t, p = r["totals"], r["new_vs_old"]
            rows.append(
                {
                    "sample": s,
                    "PXF evals": f"{t['old']['total']:,}",
                    "RU evals": f"{t['new']['total']:,}",
                    "fewer": pct(r["fewer_evals"]),
                    "wrong PXF / RU": f"{r['wrong']['old']['n']} / {r['wrong']['new']['n']}",
                    "fixed": p["fixed"],
                    "broken": p["broken"],
                    "median higher (deg)": f"{p['median_higher']:+.5f}",
                    "grid PXF / RU": (
                        f"{r['grid_fallbacks_left']['old']} / {r['grid_fallbacks_left']['new']}"
                    ),
                }
            )
            c = r["tries"]
            lo, hi = c["success_wilson"]
            tr.append(
                {
                    "sample": s,
                    "seed steps": f"{c['seed_steps']:,}",
                    "ended by a try": f"{c['success']:,} ({pct(c['success_share'])})",
                    "Wilson 95%": f"{pct(lo)} to {pct(hi)}",
                    "first fit / second fit": (
                        f"{c['success_first_fit']} / {c['success_second_fit']}"
                    ),
                    "failed tries": c["failed_tries"],
                    "try evals": f"{c['try_evals']:,}",
                    "indexing calls": f"{c['index_calls']:,}",
                    "try s": f"{c['try_seconds']:.0f}",
                    "indexing s": f"{c['index_seconds']:.0f}",
                }
            )
            wl.append(
                {
                    "sample": s,
                    "PXF wall s": r["job"]["old"]["wall_s"],
                    "RU wall s": r["job"]["new"]["wall_s"],
                    "RU label": r["job"]["new"]["label"],
                    "RU load mean / max": (
                        f"{r['job']['new']['load_mean']:.1f} / {r['job']['new']['load_max']:.1f}"
                    ),
                    "PXF label": r["job"]["old"]["label"],
                }
            )
        out["ru_live_evals"] = DT.markdown_table(rows, list(rows[0]))
        out["ru_live_tries"] = DT.markdown_table(tr, list(tr[0]))
        out["ru_live_wall"] = DT.markdown_table(wl, list(wl[0]))
        j = live["judge"]
        g = [
            {
                "gate": "r1 evaluations",
                "frozen bar": "at least 15% fewer than 23,402,740",
                "measured": f"{pct(j['r1']['fewer'])} fewer ({j['r1']['new']:,})",
                "verdict": "met" if j["r1"]["ok"] else "not met",
            }
        ]
        for s, v in j["r2"].items():
            g.append(
                {
                    "gate": f"r2 {s}",
                    "frozen bar": "broken minus fixed at most 3; median at most +0.001; grid at "
                    "most +2",
                    "measured": f"net {v['net']:+d}; median {v['median_higher']:+.5f}; "
                    f"grid {v['grid_extra']:+d}",
                    "verdict": "met" if v["ok"] else "not met",
                }
            )
        g.append(
            {
                "gate": "r3 wall (contended)",
                "frozen bar": "at least 15% lower than 4,545 s",
                "measured": f"{pct(j['r3']['fewer'])} lower ({j['r3']['new_wall_s']:,} s)",
                "verdict": "met" if j["r3"]["ok"] else "not met",
            }
        )
        out["ru_gates"] = DT.markdown_table(g, list(g[0]))
    if h8:
        t = h8["totals"]
        rows = [
            {
                "arm": a,
                "evals": f"{t[k]:,}",
                "wrong": h8["wrong"][k]["n"],
                "median error (deg)": f"{h8['median_err_deg'][k]:.5f}",
                "wall s": h8["job"][k]["wall_s"],
                "label": h8["job"][k]["label"],
            }
            for a, k in (("S8", "s8"), ("PXF J5SD", "pxf"), ("RU J5SD", "ru"))
        ]
        out["ru_held8_arms"] = DT.markdown_table(rows, list(rows[0]))
        c = h8["tries"]
        out["ru_held8_gates"] = DT.markdown_table(
            [
                {
                    "gate": "h1 RU vs S8",
                    "frozen bar": "at least 20% fewer; net at most 5; median at most +0.005",
                    "measured": f"{pct(h8['h1']['fewer'])} fewer; fixed {h8['ru_vs_s8']['fixed']}"
                    f" broken {h8['ru_vs_s8']['broken']}; median "
                    f"{h8['ru_vs_s8']['median_higher']:+.5f}",
                    "verdict": "met" if h8["h1"]["ok"] else "not met",
                },
                {
                    "gate": "h2 RU vs PXF",
                    "frozen bar": "at least 10% fewer; net at most 5; median at most +0.001",
                    "measured": f"{pct(h8['h2']['fewer'])} fewer; fixed {h8['ru_vs_pxf']['fixed']}"
                    f" broken {h8['ru_vs_pxf']['broken']}; median "
                    f"{h8['ru_vs_pxf']['median_higher']:+.5f}",
                    "verdict": "met" if h8["h2"]["ok"] else "not met",
                },
            ],
            ["gate", "frozen bar", "measured", "verdict"],
        )
        lo, hi = c["success_wilson"]
        out["ru_held8_tries"] = DT.markdown_table(
            [
                {
                    "seed steps": f"{c['seed_steps']:,}",
                    "ended by a try": f"{c['success']:,} ({pct(c['success_share'])})",
                    "Wilson 95%": f"{pct(lo)} to {pct(hi)}",
                    "failed tries": c["failed_tries"],
                    "try evals": f"{c['try_evals']:,}",
                }
            ],
            ["seed steps", "ended by a try", "Wilson 95%", "failed tries", "try evals"],
        )
    return out


def plot(live: Any, h8: Any) -> Optional[Path]:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    if not live:
        return None
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(11, 4.2), gridspec_kw=dict(width_ratios=[2, 1]))
    labels: List[str] = []
    idx: List[float] = []
    tri: List[float] = []
    rest: List[float] = []
    for s, r in live["rows"].items():
        for arm, wall, ix, ty in (
            ("PXF", r["job"]["old"]["wall_s"], r["index_old"]["seconds"], 0.0),
            (
                "RU",
                r["job"]["new"]["wall_s"],
                r["tries"]["index_seconds"],
                r["tries"]["try_seconds"],
            ),
        ):
            labels.append(f"{s}\n{arm}")
            idx.append(ix)
            tri.append(ty)
            rest.append(max(0.0, wall - ix - ty))
    if h8:
        for arm, k in (("S8", "s8"), ("PXF", "pxf"), ("RU", "ru")):
            labels.append(f"held8\n{arm}")
            w = h8["job"][k]["wall_s"]
            ix = {"s8": 0.0, "pxf": h8["index_pxf"]["seconds"], "ru": h8["tries"]["index_seconds"]}[
                k
            ]
            ty = h8["tries"]["try_seconds"] if k == "ru" else 0.0
            idx.append(ix)
            tri.append(ty)
            rest.append(max(0.0, float(w) - ix - ty))
    x = list(range(len(labels)))
    a1.bar(x, idx, label="indexing", color="tab:red")
    a1.bar(x, tri, bottom=idx, label="tries (RU)", color="tab:orange")
    a1.bar(x, rest, bottom=[i + t for i, t in zip(idx, tri)], label="everything else", color="0.7")
    a1.set_xticks(x)
    a1.set_xticklabels(labels, fontsize=8)
    a1.set_ylabel("wall seconds (contended)")
    a1.set_title("Wall per run, whole J5SD (contended; held8 S8: no split)")
    a1.legend(fontsize=8)
    names = list(live["rows"])
    br = [live["rows"][s]["new_vs_old"] for s in names]
    if h8:
        names.append("held8")
        br.append(h8["ru_vs_pxf"])
    w = 0.38
    a2.bar([i - w / 2 for i in range(len(names))], [b["fixed"] for b in br], w, label="fixed")
    a2.bar([i + w / 2 for i in range(len(names))], [b["broken"] for b in br], w, label="broken")
    a2.set_xticks(range(len(names)))
    a2.set_xticklabels(names)
    a2.set_title("RU against PXF: voxels fixed / broken")
    a2.legend(fontsize=8)
    fig.tight_layout()
    p = RD / "ru_wall.png"
    fig.savefig(p, dpi=130)
    return p


def main() -> None:
    live, h8 = load("live"), load("held8")
    tb = tables(live, h8)
    DT.write_tables(RD / "tables.md", tb)
    png = plot(live, h8)
    (RD / "summary.json").write_text(
        json.dumps(dict(live=live, held8=h8, plot=str(png) if png else None), indent=1) + "\n"
    )
    print("wrote", RD / "tables.md", png)


if __name__ == "__main__":
    main()
