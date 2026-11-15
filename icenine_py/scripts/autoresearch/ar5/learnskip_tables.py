"""AR5 step 18 (stopped at its gate): the plot and the doc tables from stored results only.

  uv run python scripts/autoresearch/ar5/learnskip_tables.py

Reads ``loso.json``, ``g2.json``, ``stats_mid.json``, ``gate_fail.json`` and ``broken_mid.json`` in
``benchmarks/autoresearch/ar5/learnskip/``; writes ``ar5_learnskip_gate.png`` and ``tables.md``
(``doc_tables.write_tables`` blocks). The judge's grid-search comparison is shown under the rule
as it was frozen (one-sided: J5L's count at most J5's plus 1) and as edited (the difference at
most 1 in either direction), to show that the edit changes no verdict.
"""

import json
import sys
from pathlib import Path
from typing import Any, Dict

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

HERE = Path(__file__).resolve().parent
ICE = HERE.parents[2]
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(ICE / "scripts" / "common"))
import doc_tables as DT  # noqa: E402
import learnskip_rules as LR  # noqa: E402

D = ICE / "benchmarks" / "autoresearch" / "ar5" / "learnskip"


def md(rows: Any) -> str:
    head, *body = rows
    out = ["| " + " | ".join(head) + " |", "|" + "---|" * len(head)]
    out += ["| " + " | ".join(str(c) for c in r) + " |" for r in body]
    return "\n".join(out)


def judge_original(j5l: Dict[str, Any], j5: Dict[str, Any], s8: Dict[str, Any]) -> Dict[str, Any]:
    """``LR.judge`` with the grid comparison as frozen (one-sided)."""
    edited = LR.judge(j5l, j5, s8)
    c = LR.JUDGE
    r4 = bool(
        j5l["wrong_seeds"] - j5["wrong_seeds"] <= c["r4_wrong_seeds_more_max"]
        and j5l["grid"] <= j5["grid"] + c["r4_grid_within"]
        and j5l["reseeds"] <= j5["reseeds"] * (1 + c["r4_reseeds_more_max_share"])
    )
    out = dict(edited, r4=r4)
    out["passes"] = bool(
        out["r1"] and out["r2"] and out["r3"] and r4 and (out["r5"] or not out["j5_passes_r5"])
    )
    return out


def main() -> None:
    loso = json.loads((D / "loso.json").read_text())
    g2 = json.loads((D / "g2.json").read_text())
    st = json.loads((D / "stats_mid.json").read_text())
    br = json.loads((D / "broken_mid.json").read_text())
    sm = st["summary"]
    j5l, j5, s8 = sm["j5l"], sm["j5"], sm["s8"]
    mid = loso["samples"]["mid"]
    live_share = 100 * (1 - j5l["total"] / j5["total"])
    # ---- plot: saving, then fixed and broken
    fig, ax = plt.subplots(1, 2, figsize=(10, 4))
    names = ["S(0.3, both)\nreplay", "learned skip\nreplay", "learned skip\nlive"]
    vals = [mid["s_share_pct"], mid["share_pct"], live_share]
    ax[0].bar(names, vals, color=["#aaaaaa", "#4477aa", "#cc6677"])
    for i, v in enumerate(vals):
        ax[0].text(i, v + 0.1, f"{v:.2f}%", ha="center")
    ax[0].axhline(4.5, color="k", ls=":", lw=1)
    ax[0].text(2.45, 4.6, "gate g1: 4.5%", ha="right", fontsize=8)
    ax[0].set_ylabel("evaluations saved, % of J5 with the pass (mid)")
    ax[0].set_title("Saving: the replay and the live run")
    labels = ["J5L vs J5", "J5L vs S8", "J5 vs S8"]
    pairs = [sm["j5l"]["vs_j5"], sm["j5l"]["vs_s8"], sm["j5"]["vs_s8"]]
    x = range(3)
    ax[1].bar(
        [i - 0.2 for i in x], [p["broken"] for p in pairs], 0.4, label="broken", color="#cc6677"
    )
    ax[1].bar(
        [i + 0.2 for i in x], [p["fixed"] for p in pairs], 0.4, label="fixed", color="#4477aa"
    )
    ax[1].set_xticks(list(x))
    ax[1].set_xticklabels(labels)
    ax[1].set_ylabel("voxels (mid, after the pass)")
    ax[1].set_title("Wrong voxels: broken and fixed")
    ax[1].legend()
    fig.tight_layout()
    fig.savefig(D / "ar5_learnskip_gate.png", dpi=120)
    # ---- tables
    t: Dict[str, str] = {}
    rows = [
        [
            "held-out sample",
            "net saving",
            "wrong skips (allowed)",
            "changed final",
            "S(0.3, both)",
            "ceiling",
        ]
    ]
    for s in LR.SAMPLES:
        v = loso["samples"][s]
        rows.append(
            [
                s,
                f"{v['share_pct']:.2f}%",
                f"{v['wrong']} ({v['allowed']})",
                v["changed"],
                f"{v['s_share_pct']:.2f}%, {v['s_wrong']} wrong",
                f"{v['ceiling_pct']:.2f}%",
            ]
        )
    t["loso"] = md(rows)
    t["live_mid"] = md(
        [
            [
                "arm (mid, with the pass)",
                "evaluations",
                "broken / fixed against J5",
                "broken / fixed against S8",
                "median error (degree)",
            ],
            [
                "J5L, live, gate model",
                f"{j5l['total']:,}",
                f"{j5l['vs_j5']['broken']} / {j5l['vs_j5']['fixed']}",
                f"{j5l['vs_s8']['broken']} / {j5l['vs_s8']['fixed']}",
                f"{j5l['median']:.5f}",
            ],
            [
                "J5",
                f"{j5['total']:,}",
                "",
                f"{j5['vs_s8']['broken']} / {j5['vs_s8']['fixed']}",
                f"{j5['median']:.5f}",
            ],
            ["S8", f"{s8['total']:,}", "", "", f"{s8['median']:.5f}"],
        ]
    )
    cmp_rows = [
        ["quantity (mid)", "replay", "live"],
        ["fits", f"{g2['replay_fits']:,}", f"{g2['live_fits']:,}"],
        ["saving as a share of J5 with the pass", f"{mid['share_pct']:.2f}%", f"{live_share:.2f}%"],
        ["skipped fits", f"{mid['skipped']:,}", f"{st['learnskip']['counts']['skipped']:,}"],
        [
            "fits compared up to the replay's first wrong skip",
            f"{g2['compared']:,}",
            f"{g2['compared']:,}",
        ],
        [
            "of these, fits with a feature difference (feature 13 only)",
            "0",
            f"{g2['compared'] - g2['equal']}",
        ],
        ["first fit with a difference", "none", f"{g2['first_difference']['fit']}"],
    ]
    t["replay_vs_live"] = md(cmp_rows)
    o, e = judge_original(j5l, j5, s8), LR.judge(j5l, j5, s8)
    t["judge_edit"] = md(
        [
            ["rule on mid", "frozen grid comparison", "edited grid comparison"],
            ["(1) 3% fewer than J5", o["r1"], e["r1"]],
            ["(2) broken minus fixed against J5 at most 5", o["r2"], e["r2"]],
            ["(3) median at most 0.001 degree higher", o["r3"], e["r3"]],
            [
                f"(4) wrong seeds, grid searches ({j5l['grid']} against {j5['grid']}), reseeds",
                o["r4"],
                e["r4"],
            ],
            ["(5) against S8 (J5 passes it)", o["r5"], e["r5"]],
            ["J5L passes on mid", o["passes"], e["passes"]],
        ]
    )
    bd = br["broken_detail"]
    t["broken"] = md(
        [
            [
                "broken voxel",
                "error in J5 (degree)",
                "error in J5L (degree)",
                "fits skipped there",
                "fits run there",
            ]
        ]
        + [
            [v, f"{d['err_J5_deg']:.3f}", f"{d['err_J5L_deg']:.1f}", d["skipped"], d["ran"]]
            for v, d in bd.items()
        ]
    )
    DT.write_tables(D / "tables.md", t)
    print("wrote plot and", sorted(t))


if __name__ == "__main__":
    main()
