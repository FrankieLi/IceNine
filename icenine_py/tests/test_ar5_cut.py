"""AR5 step 10 (task ar5j): the level-0 cut on a synthetic candidate list. Written before the
scripts ``cut_*.py`` (plan "AR5 step 10", phase A step 2): the first N listed candidates are kept,
level 0 only, the quarter is counted from the cut list, and no draw outside the grid search moves.
"""

import sys
from pathlib import Path
from typing import Any, Dict, List

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "scripts" / "autoresearch" / "ar5"))

import cut_rules as CR  # noqa: E402


def _cands(n: int) -> List[Any]:
    return [dict(i=i, cost=0.001 * i) for i in range(n)]


def _inner(n: int, calls: List[Dict[str, Any]]) -> Any:
    def f(**kw: Any) -> List[Any]:
        calls.append(kw)
        return _cands(n)

    return f


def test_constants_are_the_frozen_values() -> None:
    assert CR.TOP_N == 1000
    assert CR.KEEP_FRACTION == 0.25


def test_cut_list_keeps_the_first_n_in_order() -> None:
    out = CR.cut_list(_cands(10_200), 1000)
    assert len(out) == 1000 and [c["i"] for c in out] == list(range(1000))


def test_cut_list_short_or_off_is_untouched() -> None:
    c = _cands(500)
    assert CR.cut_list(c, 1000) is c
    assert CR.cut_list(c, 0) is c
    c2 = _cands(1000)
    assert CR.cut_list(c2, 1000) is c2


def test_level0_only() -> None:
    full = np.zeros((7, 3, 3))
    calls: List[Dict[str, Any]] = []
    log: List[Dict[str, Any]] = []
    w = CR.make_cut_search(_inner(5000, calls), 1000, lambda: full, log)
    a = w(fz_orientations=full, angular_radius=0.1)
    b = w(fz_orientations=np.zeros((250, 3, 3)), angular_radius=0.05)
    c = w(fz_orientations=full.copy(), angular_radius=0.1)  # equal values, not the full set itself
    assert len(a) == 1000 and len(b) == 5000 and len(c) == 5000
    assert [r["level0"] for r in log] == [True, False, False]
    assert [(r["n_listed"], r["n_run"]) for r in log] == [(5000, 1000), (5000, 5000), (5000, 5000)]


def test_arguments_pass_through_unchanged() -> None:
    full = np.zeros((7, 3, 3))
    calls: List[Dict[str, Any]] = []
    w = CR.make_cut_search(_inner(10, calls), 1000, lambda: full, [])
    marker = object()
    w(fz_orientations=full, local_grid=marker, angular_radius=0.3, phase_index=2)
    assert calls[0]["local_grid"] is marker and calls[0]["phase_index"] == 2
    assert calls[0]["fz_orientations"] is full


def test_off_wrapper_is_the_inner_call() -> None:
    full = np.zeros((7, 3, 3))
    calls: List[Dict[str, Any]] = []
    w = CR.make_cut_search(_inner(5000, calls), 0, lambda: full, [])
    assert len(w(fz_orientations=full)) == 5000


def test_quarter_is_counted_from_the_cut_list() -> None:
    # the library: n_keep = max(1, int(len(candidates) * 0.25)) of the list it is handed
    assert CR.n_keep_after_cut(10_200, 1000) == 250
    assert CR.n_keep_after_cut(10_200, 0) == 2550
    assert CR.n_keep_after_cut(600, 1000) == 150


def test_no_random_draw_in_the_wrapper() -> None:
    # the wrapper takes no generator and the aligned draws of the grid call are those of arm G:
    # the same (voxel, "grid", 0) generator gives the same numbers whatever the cut did
    import order_draws as OD

    full = np.zeros((7, 3, 3))
    w = CR.make_cut_search(_inner(5000, []), 1000, lambda: full, [])
    g1 = OD.make_rng(12345, "grid", 0).random(5)
    w(fz_orientations=full)
    g2 = OD.make_rng(12345, "grid", 0).random(5)
    assert np.array_equal(g1, g2)
