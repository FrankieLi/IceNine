"""AR5 step 20 (phase B): pareto rows for arm ``J5C`` on ``held4`` and the five disclosed samples.
Replaces any row with step "AR5 step 20".

  uv run python scripts/autoresearch/ar5/cmarepair_pareto_rows.py

Reads ``benchmarks/autoresearch/ar5/cmarepair/summary.json``. The evaluations are the stored ``J5``
run with the pass plus the repair's counted evaluations (the pass's 24,570 score evaluations are
left out, as in step 16's totals).
"""

import json
import sys
from pathlib import Path
from typing import Any, Dict

HERE = Path(__file__).resolve().parent
ICE = HERE.parents[2]
sys.path.insert(0, str(ICE / "scripts" / "common"))
import stats  # noqa: E402

ROWS = ICE / "docs" / "autoresearch" / "pareto" / "pareto_rows.json"
SUMMARY = ICE / "benchmarks" / "autoresearch" / "ar5" / "cmarepair" / "summary.json"
STEP = "AR5 step 20"
N = 24570
NAMES = dict(
    held4="held4_24570",
    held3="held3_24570",
    held="held_24570",
    mid="mid_24570",
    grad="gradient_24570",
    big="big_24570",
)


def main() -> None:
    d = json.loads(ROWS.read_text())
    summ = json.loads(SUMMARY.read_text())
    new = []
    for k, name in NAMES.items():
        r: Dict[str, Any] = summ["rows"][k]
        j = summ["judge"][k]
        total = int(j["j5c_total"])
        lo, hi = stats.wilson(int(r["wrong_J5C"]), N)
        clean = k == "held4"
        new.append(
            dict(
                sample=name,
                step=STEP,
                method="BFS indexed, arm J5 + boundary pass + repair of contested voxels (120-"
                "evaluation CMA local optimization from the own end and each neighbour "
                "candidate; arm J5C, key BoundaryPassRepair)",
                family="BFS with indexed seeds",
                boundary_pass=True,
                evals_total=total,
                evals_per_voxel=round(total / N, 2),
                per_voxel_derived=True,
                wrong=int(r["wrong_J5C"]),
                n=N,
                interval_pct=[round(100 * lo, 2), round(100 * hi, 2)],
                estimated=False,
                what_evals_include="the stored J5 run with its boundary pass plus every counted "
                f"evaluation of the repair ({int(r['repair_evals']):,})",
                source_json="icenine_py/benchmarks/autoresearch/ar5/cmarepair/summary.json :: "
                f"rows.{k}",
                source_page="icenine_py/docs/autoresearch/findings_ar5_cmarepair.md",
                flag=(
                    "Clean test sample (built after the freeze): rules 1 to 4 pass. One sample."
                    if clean
                    else "In-sample: the repair's flag and settings were chosen on this sample. "
                    "Rules 1 to 4 pass."
                ),
            )
        )
    d["rows"] = [x for x in d["rows"] if x.get("step") != STEP] + new
    ROWS.write_text(json.dumps(d, indent=1, ensure_ascii=False) + "\n")
    print(
        "rows:", len(d["rows"]), "new:", [(x["sample"], x["evals_total"], x["wrong"]) for x in new]
    )


if __name__ == "__main__":
    main()
