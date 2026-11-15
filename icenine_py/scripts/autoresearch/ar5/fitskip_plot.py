"""AR5 step 17: the plot and the pick table of the findings page (reads ``select.json``).

  uv run python scripts/autoresearch/ar5/fitskip_plot.py

Writes ``ar5/fitskip/ar5_fitskip_gate.png`` (net saving of every admissible rule on each sample, the
pick marked, the gate's line per sample) and adds the table ``pick`` to ``ar5/fitskip/tables.md``.
No run, no recording row is read: only the numbers already in ``select.json``.
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
FS = ICE / "benchmarks" / "autoresearch" / "ar5" / "fitskip"
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(ICE / "scripts" / "common"))
import fitskip_rules as FR  # noqa: E402
import doc_tables  # noqa: E402


def main() -> None:
    sel: Dict[str, Any] = json.loads((FS / "select.json").read_text())
    rules = sel["rules"]
    adm = [x["rule"] for x in sel["ranking"]]
    pick = sel["pick"]
    samples = list(FR.SAMPLES)
    fig, ax = plt.subplots(figsize=(8, 4.5))
    for k, s in enumerate(samples):
        need = 100 * float(FR.GATE["g2_net_saving_min"][s])
        ax.plot([k - 0.35, k + 0.35], [need, need], color="red", ls="--", lw=1.5)
        for name in adm:
            y = 100 * rules[name][s]["net"] / FR.EXPECT[s]["total"]
            ax.plot(k - 0.15, y, "o", color="0.6", ms=5)
        yp = 100 * rules[pick][s]["net"] / FR.EXPECT[s]["total"]
        ax.plot(k + 0.15, yp, "D", color="tab:blue", ms=8)
        ax.annotate(f"{yp:.2f}", (k + 0.15, yp), xytext=(6, -3), textcoords="offset points")
    ax.set_xticks(range(len(samples)))
    ax.set_xticklabels(samples)
    ax.set_ylabel("net saving, % of the run with the pass")
    ax.set_title(f"Admissible rules on the recording (grey), the pick {pick} (blue), gate (red)")
    ax.set_ylim(bottom=0)
    fig.tight_layout()
    fig.savefig(FS / "ar5_fitskip_gate.png", dpi=110)
    lines = [
        "Reading line: the pick on each sample's recording, against the share the gate needs.",
        "",
        "| sample | fits | skipped | wrong skips | changed final answers | net saving % "
        "| gate needs % | met |",
        "|---|---|---|---|---|---|---|---|",
    ]
    for s in samples:
        r = rules[pick][s]
        got = 100 * r["net"] / FR.EXPECT[s]["total"]
        need = 100 * float(FR.GATE["g2_net_saving_min"][s])
        lines.append(
            f"| {s} | {r['n_fits']} | {r['skipped']} | {r['wrong']} | {r['final_changed']} |"
            f" {got:.2f} | {need:.0f} | {'yes' if got >= need else 'NO'} |"
        )
    t = doc_tables.read_tables(FS / "tables.md")
    t["pick"] = "\n".join(lines)
    doc_tables.write_tables(FS / "tables.md", t)
    print("ok")


if __name__ == "__main__":
    main()
