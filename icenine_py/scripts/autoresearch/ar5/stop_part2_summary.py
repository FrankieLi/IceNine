"""AR5 step 1, part 2: summary of the ladder on the grid-fallback cases.

  uv run python scripts/autoresearch/ar5/stop_part2_summary.py

Reads ``cache/part2_grad.json`` and ``cache/part2_mid.json`` (``stop_part2_run.py``); writes
``part2_cases.json`` (every case, slim) and ``part2_summary.json`` (counts, totals, the pass rule,
prediction 2 and 3). The cases are few (16 and 15): no setting is picked on them.
"""

import json
from pathlib import Path
from typing import Any, Dict, List

import stop_rules as SR

AR5 = Path(__file__).resolve().parents[3] / "benchmarks" / "autoresearch" / "ar5"
WRONG = 1.0
KEEP = (
    "voxel reseed n_cand rank_near_truth near_truth_deg n_within_2deg in_costed20 in_first100 "
    "in_refined5 replay_evals_match replay_hit_match replay_hit stored_fin_hit R1 R3 grid ladder "
    "grid_evals"
).split()


def stats_for(rows: List[Dict[str, Any]]) -> Dict[str, Any]:
    n = len(rows)
    rung = {k: sum(1 for r in rows if r["ladder"]["rung"] == k) for k in ("R1", "R3", "grid")}
    grid_evals = sum(r["grid_evals"] for r in rows)
    ladder_evals = sum(r["ladder"]["evals"] for r in rows)
    right_grid = sum(1 for r in rows if r["grid"]["err"] <= WRONG)
    right_ladder = sum(1 for r in rows if r["ladder"]["right"])
    lost = [r["voxel"] for r in rows if r["grid"]["err"] <= WRONG and not r["ladder"]["right"]]
    gained = [r["voxel"] for r in rows if r["grid"]["err"] > WRONG and r["ladder"]["right"]]
    return dict(
        n=n,
        n_reseed=sum(1 for r in rows if r["reseed"]),
        n_cand_median=float(sorted(r["n_cand"] for r in rows)[n // 2]),
        near_truth_found=sum(1 for r in rows if r["rank_near_truth"] >= 0),
        rank0=sum(1 for r in rows if r["rank_near_truth"] == 0),
        in_costed20=sum(1 for r in rows if r["in_costed20"]),
        in_refined5=sum(1 for r in rows if r["in_refined5"]),
        in_first100_not_costed20=sum(1 for r in rows if r["in_first100"] and not r["in_costed20"]),
        replay_all_match=all(r["replay_evals_match"] and r["replay_hit_match"] for r in rows),
        r1_pass_seed_test=sum(1 for r in rows if r["R1"]["hit"] >= SR.SEED_TEST - SR.TOL),
        r1_right=sum(1 for r in rows if r["R1"]["err"] <= WRONG),
        r1_pass_and_right=sum(
            1 for r in rows if r["R1"]["hit"] >= SR.SEED_TEST - SR.TOL and r["R1"]["err"] <= WRONG
        ),
        r1_evals=sum(r["R1"]["evals"] for r in rows),
        r3_pass_seed_test=sum(1 for r in rows if r["R3"]["hit"] >= SR.SEED_TEST - SR.TOL),
        r3_right=sum(1 for r in rows if r["R3"]["err"] <= WRONG),
        r3_pass_and_right=sum(
            1 for r in rows if r["R3"]["hit"] >= SR.SEED_TEST - SR.TOL and r["R3"]["err"] <= WRONG
        ),
        r3_evals=sum(r["R3"]["evals"] for r in rows),
        r3_level0_evals=sum(r["R3"]["level0_evals"] for r in rows),
        ladder_rung=rung,
        ladder_evals=ladder_evals,
        grid_evals=grid_evals,
        ladder_share_of_grid=ladder_evals / grid_evals,
        right_grid=right_grid,
        right_ladder=right_ladder,
        grid_right_ladder_wrong=lost,
        grid_wrong_ladder_right=gained,
        verdict=SR.ladder_verdict(ladder_evals, grid_evals, right_ladder, right_grid),
    )


def main() -> None:
    rows: List[Dict[str, Any]] = []
    per: Dict[str, Any] = {}
    cases: Dict[str, Any] = {}
    for s in ("grad", "mid"):
        d = json.loads((AR5 / "cache" / f"part2_{s}.json").read_text())
        rs = d["rows"]
        for r in rs:
            r["sample"] = s
        per[s] = stats_for(rs)
        cases[s] = [{k: r[k] for k in KEEP} for r in rs]
        rows += rs
    pooled = stats_for(rows)
    pred2 = dict(
        cases_with_candidate_in_first100_not_costed20=pooled["in_first100_not_costed20"],
        n=pooled["n"],
        holds=pooled["in_first100_not_costed20"] * 2 >= pooled["n"],
    )
    out = dict(
        per_sample=per,
        pooled=pooled,
        prediction2=pred2,
        prediction3=dict(holds=pooled["verdict"]["passes"], verdict=pooled["verdict"]),
        note="level 0 of the grid search alone is run through the library's own level loop "
        "(max_local_resolution 0, final refinement stubbed); a full grid search after a failed "
        "R3 pays level 0 again",
        timing_label="contended (another agent on the machine); not a timing result",
    )
    (AR5 / "part2_cases.json").write_text(json.dumps(cases, indent=0, default=float) + "\n")
    (AR5 / "part2_summary.json").write_text(json.dumps(out, indent=1, default=float) + "\n")
    print(json.dumps(out, indent=1, default=float))


if __name__ == "__main__":
    main()
