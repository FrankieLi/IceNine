"""AR5 step 16, phase B: the library key ``SeedIndexMCStallHit 500`` (rule PH(500)) parses, is the
default 0 when absent, and stops and counts exactly as the script rule of the study
(``stallhit_rules.patch_evaluate`` as ``stallhit_run`` / ``high500_pick`` use it). No new library
code: step 12's key already takes any integer k >= 1; these tests pin k = 500.
"""

import sys
from pathlib import Path
from typing import Dict

import numpy as np
import pytest

from icenine.orientation_search import SearchParameters

# tests/ is on sys.path under pytest
from test_seed_by_indexing import _config_with
from test_seed_index_stall_hit import FakeRec, _hit_by_angle, _script_fit

sys.path.insert(0, str(Path(__file__).parent.parent / "scripts" / "autoresearch" / "ar5"))


def test_k_500_parses_and_absent_or_zero_is_off(tmp_path: Path) -> None:
    on = SearchParameters.from_config(
        _config_with(
            tmp_path, "SeedByIndexing 1\nSeedIndexStopHitRatio 0.95\nSeedIndexMCStallHit 500"
        )
    )
    assert on.seed_index_mc_stall_hit == 500 and isinstance(on.seed_index_mc_stall_hit, int)
    assert SearchParameters().seed_index_mc_stall_hit == 0  # the default stays off
    for extra in ("", "SeedIndexMCStallHit 0"):
        assert (
            SearchParameters.from_config(_config_with(tmp_path, extra)).seed_index_mc_stall_hit == 0
        )


@pytest.mark.parametrize("seed", range(10))
def test_k_500_equals_the_script_rule(seed: int) -> None:
    g = np.random.default_rng(500 + seed)
    n = int(g.integers(700, 1800))
    n_gains = int(g.integers(1, 5))
    gains = sorted({0, *[int(x) for x in g.integers(1, n, n_gains)]})
    costs = [float(0.6 + 0.3 * g.random()) for _ in range(n)]  # never a gain on their own
    for j, i in enumerate(gains):
        costs[i] = 0.5 - 0.04 * j  # strictly decreasing running best
    hits: Dict[int, float] = {}
    for i in gains:
        deg = int(round(np.degrees(np.arctan2(np.sin(np.radians(i)), np.cos(np.radians(i))))))
        hits[deg] = float(g.choice([0.5, 0.94999, 0.95, 0.99]))
    R0, c0, ev0, checks0 = _script_fit(costs, hits, 500)
    lib = FakeRec(costs)
    R1, c1, ev1, log = lib.seed_mc_refine(np.eye(3), 0.6, "V", 0, 500, 0.95, _hit_by_angle(hits))
    assert np.array_equal(np.asarray(R0), R1) and c0 == c1 and ev0 == ev1
    assert checks0 == log["checks"]
    # a check fires only 500 evaluations after the last gain, and never under a low-side rule
    for e, _, _ in log["checks"]:
        assert any(e - gi >= 500 for gi in gains if gi < e)


def test_k_500_stops_at_0_95_and_goes_on_under_it() -> None:
    costs = [0.5] + [0.9] * 700 + [0.4]
    lib = FakeRec(costs)
    end, c, ev, log = lib.seed_mc_refine(
        np.eye(3), 0.6, "V", 0, 500, 0.95, _hit_by_angle({0: 0.97})
    )
    assert log["stopped"] and c == 0.5 and log["checks"] == [[501, 0.97, True]] and ev == 502
    lib = FakeRec(costs)
    end, c, ev, log = lib.seed_mc_refine(
        np.eye(3), 0.6, "V", 0, 500, 0.95, _hit_by_angle({0: 0.10})
    )
    assert not log["stopped"] and log["checks"] == [[501, 0.10, False]]  # no low stop
    assert c == 0.4
