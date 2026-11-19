"""AR5 step 24, phase A step 6 (gate g2) and the gate's verdict: live ``J5SDB`` on ``mid`` and
``held`` against ``J5SD`` (step 22's stored arm).

  uv run python scripts/autoresearch/ar5/bf_gate.py

Reads the runs of ``bf_run_all.sh --stage gate`` (``<cache>/J5DB_<s>``, ``j5sdb/J5S_<s>``), step
22's stored arms (links) and the truth of the disclosed samples (scoring only). Nothing of
``held7`` is read. g2 (frozen, ``bf_rules.GATE``): at least 4% fewer evaluations than ``J5SD``
(``J5DB`` search + the screen pass + the wide repair, against step 22's ``J5SD`` total),
broken minus fixed at most 3 (wrong: over 1 degree from the truth, symmetry-reduced; counted by
voxel), median error over all voxels at most 0.001 degree higher, grid searches at most 2
more than ``J5D`` (= ``J5SD``), seed-test failures at most 10 more; and no stall check and no
pause in the ball pass (stop rule). Writes
``ar5/bf/gate.json`` and ``gate_tables.md``. Exit 1 when a gate fails (the plan's stop rule:
no tuning, no ``held7``). The file ``gate_passed`` is NOT written here.
"""

import json
import sys
from pathlib import Path
from typing import Any, Dict

sys.path.insert(0, str(Path(__file__).resolve().parent))
import bf_rules as BR  # noqa: E402
import bf_score as BS  # noqa: E402
import doc_tables as DT  # noqa: E402

RD = BS.ICE / "benchmarks" / "autoresearch" / "ar5" / "bf"


def g2_row(s: str) -> Dict[str, Any]:
    r = BS.one(s)
    g = BR.GATE
    row: Dict[str, Any] = dict(
        sample=s,
        j5sd_total=r["totals"]["j5sd"],
        j5sdb_total=r["totals"]["j5sdb"],
        fewer=r["fewer"]["j5sd"],
        wrong_j5sd=r["wrong"]["j5sd"]["n"],
        wrong_j5sdb=r["wrong"]["j5sdb"]["n"],
        wrong_j5sdb_wilson=r["wrong"]["j5sdb"]["wilson"],
        fixed=r["pair"]["j5sd"]["fixed"],
        broken=r["pair"]["j5sd"]["broken"],
        broken_minus_fixed=r["pair"]["j5sd"]["net"],
        median_err_higher=r["pair"]["j5sd"]["median_higher"],
        grid_j5d=r["grid"]["j5d"],
        grid_j5db=r["grid"]["j5db"],
        grid_more=r["grid"]["j5db"] - r["grid"]["j5d"],
        seed_test_j5d=r["seed_test"]["j5d"],
        seed_test_j5db=r["seed_test"]["j5db"],
        seed_test_more=r["seed_test"]["j5db"] - r["seed_test"]["j5d"],
        checks=r["checks"],
        pauses=r["pauses"],
        j5sd_total_matches_step22=r["j5sd_total_matches_step22"],
    )
    row["ok_fewer"] = bool(row["fewer"] >= g["g2_fewer_min"])
    row["ok_broken"] = bool(row["broken_minus_fixed"] <= g["g2_broken_minus_fixed_max"])
    row["ok_median"] = bool(row["median_err_higher"] <= g["g2_median_error_higher_max_deg"])
    row["ok_grid"] = bool(row["grid_more"] <= g["g2_grid_more_max"])
    row["ok_seed_test"] = bool(row["seed_test_more"] <= g["g2_seed_test_more_max"])
    row["ok_no_checks"] = bool(row["checks"] == 0 and row["pauses"] == 0)
    row["ok"] = bool(
        row["ok_fewer"]
        and row["ok_broken"]
        and row["ok_median"]
        and row["ok_grid"]
        and row["ok_seed_test"]
        and row["ok_no_checks"]
    )
    row["stages"] = r["stages"]
    row["stages_j5sd"] = r["stages_j5sd"]
    row["ball"] = r["ball"]
    return row


def main() -> int:
    g0 = json.loads((RD / "g0.json").read_text())
    g1 = json.loads((RD / "size.json").read_text())
    out: Dict[str, Any] = dict(rule=BR.GATE, g0=bool(g0["g0_passes"]), g1=bool(g1["g1_passes"]))
    out["g2"] = {s: g2_row(s) for s in BR.GATE_LIVE}
    out["g2_ok"] = all(r["ok"] for r in out["g2"].values())
    out["passed"] = bool(out["g0"] and out["g1"] and out["g2_ok"])
    (RD / "gate.json").write_text(json.dumps(out, indent=1, default=float) + "\n")
    cols = [
        "sample",
        "j5sd_total",
        "j5sdb_total",
        "fewer",
        "wrong_j5sd",
        "wrong_j5sdb",
        "fixed",
        "broken",
        "median_err_higher",
        "grid_more",
        "seed_test_more",
        "ok",
    ]
    fmt = dict(j5sd_total=",d", j5sdb_total=",d", fewer=".4f", median_err_higher=".5f")
    DT.write_tables(
        RD / "gate_tables.md",
        dict(ar5_bf_gate_g2=DT.markdown_table(list(out["g2"].values()), cols, fmt)),
    )
    print(json.dumps({k: v for k, v in out.items() if k != "rule"}, indent=1, default=float))
    return 0 if out["passed"] else 1


if __name__ == "__main__":
    sys.exit(main())
