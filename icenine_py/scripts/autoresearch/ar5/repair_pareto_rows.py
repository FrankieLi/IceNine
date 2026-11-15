"""AR5 step 19 (stopped at its gate, D-20261115-N8): pareto rows for arm ``J5R`` on the five
disclosed samples, flagged "failed gate (cost)". Replaces any row with step "AR5 step 19".

  uv run python scripts/autoresearch/ar5/repair_pareto_rows.py

Reads ``benchmarks/autoresearch/ar5/repair/summary.json``. The evaluations are the stored ``J5``
run with the pass plus the repair's counted evaluations.
"""

import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ICE = HERE.parents[2]
sys.path.insert(0, str(ICE / "scripts" / "common"))
import stats  # noqa: E402

ROWS = ICE / "docs" / "autoresearch" / "pareto" / "pareto_rows.json"
SUMMARY = ICE / "benchmarks" / "autoresearch" / "ar5" / "repair" / "summary.json"
STEP = "AR5 step 19"
N = 24570
NAMES = dict(
    held3="held3_24570", held="held_24570", mid="mid_24570", grad="gradient_24570", big="big_24570"
)


def main() -> None:
    d = json.loads(ROWS.read_text())
    summ = json.loads(SUMMARY.read_text())
    new = []
    for k, name in NAMES.items():
        r = summ[k]
        total = int(r["j5_total"]) + int(r["repair_evals"])
        lo, hi = stats.wilson(int(r["wrong_J5R"]), N)
        new.append(
            dict(
                sample=name,
                step=STEP,
                method="BFS indexed, arm J5 + boundary pass + repair of contested voxels (MC local "
                "optimization from the own end and each neighbour candidate; arm J5R)",
                family="BFS with indexed seeds",
                boundary_pass=True,
                evals_total=total,
                evals_per_voxel=round(total / N, 2),
                per_voxel_derived=True,
                wrong=int(r["wrong_J5R"]),
                n=N,
                interval_pct=[round(100 * lo, 2), round(100 * hi, 2)],
                estimated=False,
                what_evals_include="the stored J5 run with its boundary pass plus every counted "
                f"evaluation of the repair ({int(r['repair_evals']):,})",
                source_json="icenine_py/benchmarks/autoresearch/ar5/repair/summary.json :: " + k,
                source_page="icenine_py/MIGRATION_HISTORY.md :: Autoresearch AR5 step 19 (the "
                "findings page text waits in .claude/reports/; no page under docs/ yet)",
                flag="FAILED GATE (cost): the repair's evaluations are "
                f"{r['share_pct']:.2f}% of the J5 run, over the cap of {r['cap_pct']}% "
                "(D-20261115-N8). Accuracy as hoped but in-sample: the flag was chosen on these "
                "samples; held4 was not built. Not a result to build on.",
            )
        )
    d["rows"] = [x for x in d["rows"] if x.get("step") != STEP] + new
    ROWS.write_text(json.dumps(d, indent=1, ensure_ascii=False) + "\n")
    print(
        "rows:", len(d["rows"]), "new:", [(x["sample"], x["evals_total"], x["wrong"]) for x in new]
    )


if __name__ == "__main__":
    main()
