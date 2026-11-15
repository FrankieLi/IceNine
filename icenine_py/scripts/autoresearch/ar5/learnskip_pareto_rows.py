"""AR5 step 18 (stopped at its gate): one pareto row for the live ``J5L`` run on ``mid`` with the
boundary pass, flagged "failed gate". Replaces any row with step "AR5 step 18".

  uv run python scripts/autoresearch/ar5/learnskip_pareto_rows.py

Reads ``benchmarks/autoresearch/ar5/learnskip/stats_mid.json``. The row format has a ``flag`` field
and no status field, so the status goes in ``flag``.
"""

import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
ICE = HERE.parents[2]
ROWS = ICE / "docs" / "autoresearch" / "pareto" / "pareto_rows.json"
STATS = ICE / "benchmarks" / "autoresearch" / "ar5" / "learnskip" / "stats_mid.json"
STEP = "AR5 step 18"


def main() -> None:
    text = ROWS.read_text()
    d = json.loads(text)
    st = json.loads(STATS.read_text())
    a = st["arms"]["J5L"]
    s = a["stages"]
    n = 24570
    row = dict(
        sample="mid_24570",
        step=STEP,
        method="BFS indexed, arm J5 with a learned skip of BFS local optimizations (boosted trees "
        "on 16 features, gate model trained on four other samples) + boundary pass (run J5L)",
        family="BFS with indexed seeds",
        boundary_pass=True,
        evals_total=int(s["total"]),
        evals_per_voxel=round(s["total"] / n, 2),
        per_voxel_derived=True,
        wrong=int(a["after"]["wrong"]),
        n=n,
        interval_pct=[
            round(100 * a["after"]["wilson_lo"], 2),
            round(100 * a["after"]["wilson_hi"], 2),
        ],
        estimated=False,
        what_evals_include="every counted evaluation of the run with the boundary pass, including "
        "one counted check of the start for each of the 52,075 local optimizations; a model call "
        "counts zero",
        source_json="icenine_py/benchmarks/autoresearch/ar5/learnskip/stats_mid.json :: arms.J5L",
        source_page="icenine_py/docs/autoresearch/findings_ar5_learnskip.md",
        flag="FAILED GATE: the step stopped at its gate (D-20261115-N5). Gate model, not the final "
        "model; one run on mid, the sample the gate model was checked on; broken minus fixed "
        "against S8 is 8 (limit 5). Not a result to build on.",
    )
    d["rows"] = [r for r in d["rows"] if r.get("step") != STEP] + [row]
    ROWS.write_text(json.dumps(d, indent=1, ensure_ascii=False) + "\n")
    print("rows:", len(d["rows"]), "row:", row["evals_total"], row["wrong"], row["interval_pct"])


if __name__ == "__main__":
    main()
