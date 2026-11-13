"""AR5 step 9 (task ar5i): rule G1 on synthetic steps, the wrapper in front of the grid fallback,
and that the aligned draws are unchanged. Written before the scripts ``rescue_*.py`` (plan "AR5
step 9", phase A step 2)."""

import sys
from pathlib import Path
from typing import Any, Dict, List

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "scripts" / "autoresearch" / "ar5"))

import order_draws as OD  # noqa: E402
import rescue_rules as RR  # noqa: E402


def _try(hits: Dict[int, float], evals: int = 1000) -> Any:
    calls: List[int] = []

    def f(k: int) -> Dict[str, Any]:
        calls.append(k)
        return dict(hit=hits.get(k, 0.5), evals=evals)

    f.calls = calls  # type: ignore[attr-defined]
    return f


def test_constants_are_the_frozen_values() -> None:
    assert (RR.TAU, RR.SIGMA0_DEG, RR.EVALS, RR.MAX_TRIES, RR.SEED_TEST) == (
        0.95,
        0.25,
        1000,
        5,
        0.80,
    )


def test_tries_in_order_of_cost_ties_by_position_and_capped() -> None:
    assert RR.rescue_order([0.3, 0.1, 0.2, 0.1]) == [1, 3, 2, 0]
    assert RR.rescue_order([0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3]) == [6, 5, 4, 3, 2]


def test_first_end_at_tau_is_taken_and_later_ends_are_not_tried() -> None:
    f = _try({1: 0.96, 0: 0.99})
    r = RR.rescue_choose([0.3, 0.1, 0.2], f)
    assert r["rescued"] and r["taken"] == 0 and f.calls == [1]
    assert r["evals"] == 1000 and r["tries"][0]["k"] == 1


def test_a_miss_goes_on_to_the_next_cost_and_counts_every_try() -> None:
    f = _try({2: 0.95})  # exactly tau counts
    r = RR.rescue_choose([0.3, 0.1, 0.2], f)
    assert f.calls == [1, 2] and r["taken"] == 1 and r["evals"] == 2000


def test_no_end_reaches_tau_all_tries_counted_and_not_rescued() -> None:
    f = _try({0: 0.949999})
    r = RR.rescue_choose([0.3, 0.1, 0.2], f)
    assert not r["rescued"] and r["taken"] == -1 and len(r["tries"]) == 3
    assert r["evals"] == 3000


def test_no_kept_end_is_not_rescued_and_costs_nothing() -> None:
    r = RR.rescue_choose([], _try({}))
    assert not r["rescued"] and r["evals"] == 0 and r["tries"] == []


class _Rec:
    """Stand-in for the reconstructor's eval counters."""

    def __init__(self) -> None:
        self._last_global_evals = 0
        self._last_local_evals = 0

    @property
    def last_eval_counts(self) -> tuple:
        return (
            self._last_global_evals,
            self._last_local_evals,
            self._last_global_evals + self._last_local_evals,
        )


def _setup(hits: Dict[int, float], kept_costs: List[float], key: int = 7) -> Dict[str, Any]:
    rec = _Rec()
    grid_calls: List[Dict[str, Any]] = []

    def grid(**kw: Any) -> str:
        grid_calls.append(kw)
        rec._last_global_evals, rec._last_local_evals = 300_000, 35_000
        return "GRID"

    step = dict(
        key=key,
        voxel=42,
        used=False,
        kept=[dict(R=np.eye(3), cost=c, hit=0.5) for c in kept_costs],
    )
    log: List[Dict[str, Any]] = []

    def try_end(e: Dict[str, Any], v: Any, ph: int) -> Dict[str, Any]:
        k = [x["cost"] for x in step["kept"]].index(e["cost"])
        return dict(hit=hits.get(k, 0.5), cost=e["cost"], evals=1000, result=f"R{k}")

    w = RR.make_wrapper(grid, lambda: step, lambda v: int(v), try_end, rec, log)
    return dict(w=w, rec=rec, grid=grid_calls, log=log, step=step)


def test_wrapper_rescued_returns_the_rescue_and_no_grid_call_and_counts_only_rescue() -> None:
    s = _setup({1: 0.97}, [0.3, 0.1])
    out = s["w"](voxel_vertices=7, phase_index=0, rng="RNG")
    assert out == "R1" and s["grid"] == []
    assert s["rec"].last_eval_counts[2] == 1000
    row = s["log"][0]
    assert row["attempted"] and row["rescued"] and row["grid_evals"] == 0 and row["voxel"] == 42


def test_wrapper_failed_rescue_runs_grid_unchanged_and_adds_the_rescue_evals() -> None:
    s = _setup({}, [0.3, 0.1])
    out = s["w"](voxel_vertices=7, phase_index=0, rng="RNG")
    assert out == "GRID"
    assert s["grid"] == [dict(voxel_vertices=7, phase_index=0, rng="RNG")]  # same arguments
    assert s["rec"].last_eval_counts[2] == 335_000 + 2000
    row = s["log"][0]
    assert row["attempted"] and not row["rescued"] and row["grid_evals"] == 335_000
    assert row["rescue_evals"] == 2000 and len(row["tries"]) == 2
    assert all("result" not in t for t in row["tries"])


def test_log_rows_are_json_serialisable_with_array_ends() -> None:
    import json

    s = _setup({1: 0.97}, [0.3, 0.1])
    orig = s["w"]
    assert orig is not None

    # try_end that returns an ndarray end, as the real rescue does
    def try_end(e: Dict[str, Any], v: Any, ph: int) -> Dict[str, Any]:
        return dict(hit=0.97, cost=0.0, evals=1000, R=np.eye(3), result="R")

    rec = _Rec()
    log: List[Dict[str, Any]] = []
    w = RR.make_wrapper(lambda **kw: "GRID", lambda: s["step"], lambda v: int(v), try_end, rec, log)
    w(voxel_vertices=7, phase_index=0, rng=None)
    json.dumps(log)
    assert log[0]["tries"][0]["R"] == list(np.eye(3).ravel())


def test_wrapper_without_a_matching_unused_step_is_the_plain_grid_call() -> None:
    for step_edit in ({"key": 8}, {"used": True}, {"kept": []}):
        s = _setup({0: 0.99}, [0.1])
        s["step"].update(step_edit)
        out = s["w"](voxel_vertices=7, phase_index=0, rng=None)
        assert out == "GRID" and len(s["grid"]) == 1
        assert not s["log"][0]["attempted"] and s["log"][0]["rescue_evals"] == 0


def test_a_step_is_rescued_once_the_next_call_is_plain() -> None:
    s = _setup({0: 0.99}, [0.1])
    assert s["w"](voxel_vertices=7, phase_index=0, rng=None) == "R0"
    assert s["w"](voxel_vertices=7, phase_index=0, rng=None) == "GRID"


def test_rescue_draws_nothing_from_the_aligned_book() -> None:
    book = OD.DrawBook()
    before = book.n_issued
    s = _setup({0: 0.99}, [0.1])
    s["w"](voxel_vertices=7, phase_index=0, rng=None)
    assert book.n_issued == before
    # the grid draw of arm G at a voxel is the (voxel, grid, count) draw, as in arm A
    a = OD.make_rng(7, "grid", 0).random(3)
    b = OD.DrawBook().rng(7, "grid").random(3)
    assert np.array_equal(a, b)


def test_aligned_draws_unchanged_by_this_step() -> None:
    g = OD.make_rng(123, "grid", 2)
    assert (
        g.random()
        == np.random.default_rng(np.random.SeedSequence([123, OD.kind_id("grid"), 2])).random()
    )


# ---- the verdict rules (rescue_judge) ---------------------------------------------------------
import rescue_judge as RJ  # noqa: E402

_OK = dict(a=True, b=True, c=True, passes=True)
_BAD = dict(a=False, b=True, c=True, passes=False)


def _rules(sample: str = "held", **kw: Any) -> Dict[str, Any]:
    base: Dict[str, Any] = dict(
        sample=sample,
        total_g=80,
        total_a=100,
        fallbacks_g=1,
        fallbacks_a=5,
        net_vs_a=0,
        med_g=0.02,
        med_a=0.02,
        wrong_rescued=0,
        identical_to_a=True,
        house_g=_OK,
        house_a=_OK,
    )
    base.update(kw)
    return RJ.rules(**base)


def test_pct_fewer_and_house_rule_boundaries() -> None:
    assert RJ.pct_fewer(80, 100) == 20.0 and RJ.pct_fewer(5, 0) == 0.0
    assert RJ.house(80, 100, 5, 0, 0.0, 0.0)["passes"]  # exactly 20% fewer, net 5
    assert not RJ.house(81, 100, 0, 0, 0.0, 0.0)["a"]
    assert not RJ.house(50, 100, 6, 0, 0.0, 0.0)["b"]
    assert not RJ.house(50, 100, 0, 0, 0.0061, 0.0)["c"]


def test_rule_1_by_sample_and_the_10_percent_edge() -> None:
    assert _rules("held")["r1"]
    assert not _rules("held", total_g=91)["r1"]  # 9% fewer
    assert _rules("mid", total_g=90)["r1"]  # exactly 10%
    assert not _rules("held", fallbacks_g=5)["r1"]  # no fewer fallbacks
    assert _rules("grad", total_g=99, fallbacks_g=5)["r1"]
    assert not _rules("big", identical_to_a=False)["r1"]


def test_rules_2_3_4_and_the_a_fails_5_clause() -> None:
    assert not _rules(net_vs_a=6)["passes"]
    assert not _rules(med_g=0.0211)["passes"]
    assert not _rules(wrong_rescued=2)["passes"] and _rules(wrong_rescued=1)["passes"]
    r = _rules(house_g=_BAD)
    assert not r["r5"] and not r["passes"]
    assert _rules(house_g=_BAD, house_a=_BAD)["passes"]  # A fails (5): judged by (1) to (4)


def test_step_verdict_needs_every_sample() -> None:
    assert RJ.step_verdict(dict(held=True, mid=True, grad=True, big=True))
    assert not RJ.step_verdict(dict(held=False, mid=True, grad=True, big=True))
