"""Task A tooling: the C++ ledger parser and the invariants of the stored ledger results."""

import json
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts" / "eval_reduction"))
import cpp_ledger_parse as CLP  # noqa: E402

LEDGER = ROOT / "benchmarks" / "eval_reduction" / "A" / "ledger" / "runs"

LOG = """Normal Reconstruction 0.7 sec, evals=15907
LG L 0 n_fz 10 n_local 9 disc 100 rerecost 8 quick 96 cand 8 kept 2
LG L 1 n_fz 2 n_local 9 disc 30 rerecost 4 quick 48 cand 4 kept 1
LG FIND 20 n_in 30 VAR 10 FINAL 1 converged 1 TOTAL 317
"""


def test_cpp_ledger_parse(tmp_path: Path) -> None:
    p = tmp_path / "cpp_1.log"
    p.write_text(LOG)
    r = CLP.parse(p)
    t = r["totals"]
    # screen = n_fz * n_local, first re-cost = disc - screen
    assert (t["screen"], t["recost"]) == (10 * 9 + 2 * 9, (100 - 90) + (30 - 18))
    assert (t["rerecost"], t["quick"], t["find"], t["variance"], t["final"]) == (12, 144, 20, 10, 1)
    assert t["sum"] == 317 == r["total_counter"]


@pytest.mark.skipif(not LEDGER.exists(), reason="stored ledger results missing")
def test_stored_ledger_stage_sums_and_unit_costs() -> None:
    files = sorted(LEDGER.glob("*_v*.json"))
    assert len(files) >= 24
    for f in files:
        r = json.loads(f.read_text())
        # the stages add up to every evaluation the cost function counted
        assert r["totals"]["sum"] == r["evals_total"] == r["library_eval_counts"][2], f.name
        for lv in r["levels"]:
            assert lv["screen_evals"] == lv["fz_in"] * lv["n_local"], f.name  # 9 per FZ point
            assert lv["recost_evals"] == lv["passed_screen"], f.name  # 1 per passed pair
            assert lv["quick_evals"] == 12 * lv["candidates"], f.name  # 12 per candidate
            assert lv["kept"] + lv["dropped"] == lv["candidates"], f.name
