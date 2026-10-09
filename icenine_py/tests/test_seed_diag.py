"""scripts/phase_d/seed_diag_summary.py: the stage grouping of the per-phase evaluation counts."""

import sys
from pathlib import Path

ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(ROOT / "scripts" / "phase_d"))

import seed_diag_summary as S  # noqa: E402


def test_stage_evals_groups_and_conserves() -> None:
    bp = {
        "discrete_L0|pr3": {"n": 10},
        "discrete_L0|pr0": {"n": 4},
        "quick_L0|pr0": {"n": 44},
        "discrete_L2|pr3": {"n": 3},
        "quick_L3|pr0": {"n": 5},
        "find|pr0": {"n": 7},
        "variance|pr0": {"n": 100},
        "final|pr0": {"n": 1},
    }
    s = S.stage_evals({"by_phase": bp})
    assert s == dict(disc_L0=14, quick_L0=44, disc_L1_3=3, quick_L1_3=5, find=7, variance=100)
    assert sum(s.values()) == sum(v["n"] for v in bp.values()) - 1  # the final evaluation
