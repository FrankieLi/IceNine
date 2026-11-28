"""AR5 step 14 (task ar5n): rule G2 on synthetic steps, the wrapper between rule G1 and the grid
search, the aligned draws, the evaluation counts, the read guards, and the gate rules. Written
before the probe scripts ran (plan "AR5 step 14", phase A step 2)."""

import math
import sys
from pathlib import Path
from typing import Any, Dict, List

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "scripts" / "autoresearch" / "ar5"))

import order_draws as OD  # noqa: E402
import polish_gate as PGT  # noqa: E402
import polish_guard as PG  # noqa: E402
import polish_rules as PR  # noqa: E402
import rescue_rules as RR  # noqa: E402


def test_constants_and_bounds_are_the_frozen_values() -> None:
    assert (PR.TAU, PR.RNG_SEED, PR.MAX_TRIES) == (0.95, 1, 5)
    assert PR.FORMS["G2a"] == dict(kind="cma", sigma0_deg=0.5, evals=3000)
    b = PR.FORMS["G2b"]
    assert (b["side_deg"], b["level"], b["n_grid"], b["n_starts"]) == (1.5, 2, 64, 2)
    assert (b["sigma0_deg"], b["evals"]) == (0.1, 500)
    assert PR.case_bound("G2a") == 15_005 and PR.case_bound("G2b") == 5_330
    assert PR.case_bound("G2a") + PR.case_bound("G2b") == 20_335  # the plan's 28 x 20,335


def test_tries_follow_g1s_order_and_the_first_end_at_tau_is_taken() -> None:
    calls: List[int] = []

    def f(j: int) -> Dict[str, Any]:
        calls.append(j)
        return dict(hit={1: 0.95, 3: 0.99}.get(j, 0.5), evals=3001)

    r = PR.polish_choose(5, f)
    assert calls == [0, 1] and r["accepted"] and r["taken"] == 1
    assert r["evals"] == 6002 and [t["j"] for t in r["tries"]] == [0, 1]


def test_all_misses_count_every_try_and_cap_at_five() -> None:
    calls: List[int] = []

    def f(j: int) -> Dict[str, Any]:
        calls.append(j)
        return dict(hit=0.949999, evals=3001)

    r = PR.polish_choose(7, f)
    assert calls == [0, 1, 2, 3, 4] and not r["accepted"] and r["taken"] == -1
    assert r["evals"] == 15_005
    assert PR.polish_choose(0, f)["evals"] == 0


def test_pick_starts_lowest_cost_ties_by_position() -> None:
    assert PR.pick_starts([0.5, 0.1, 0.3, 0.1], 2) == [1, 3]
    assert PR.pick_starts([0.2, 0.2, 0.2], 2) == [0, 1]


def test_g2b_grid_is_64_points_about_a_third_of_a_degree_apart_and_composes_left() -> None:
    g = PR.grid_matrices("G2b")
    assert g.shape == (64, 3, 3)
    ang = [math.degrees(math.acos(max(-1.0, min(1.0, (np.trace(m) - 1) / 2)))) for m in g]
    assert max(ang) < 1.5 * math.sqrt(3)  # inside the side-1.5-degree cube's diagonal
    R = np.array([[0, -1, 0], [1, 0, 0], [0, 0, 1.0]])
    c = PR.candidates(g, R)
    assert np.allclose(c[5], g[5] @ R)


def _fakes(costs: List[float], hits: List[float]) -> Any:
    log: Dict[str, Any] = dict(scores=0, cma=[])

    def score(Rc: np.ndarray) -> float:
        i = log["scores"]
        log["scores"] += 1
        return costs[i]

    def cma(Rs: np.ndarray, sigma: float, n: int, rng: Any) -> Dict[str, Any]:
        k = len(log["cma"])
        log["cma"].append((Rs.copy(), sigma, n))
        return dict(hit=hits[k], cost=0.1 * (k + 1), evals=n + 1, R=Rs, result=f"r{k}")

    return score, cma, log


def test_g2a_is_one_cma_with_its_values_and_counts_its_final_evaluation() -> None:
    score, cma, log = _fakes([], [0.97])
    r = PR.one_try("G2a", np.eye(3), score, cma)
    assert log["scores"] == 0 and len(log["cma"]) == 1
    assert log["cma"][0][1:] == (0.5, 3000)
    assert r["evals"] == 3001 and r["grid_evals"] == 0 and r["hit"] == 0.97


def test_g2b_scores_64_points_once_then_cma_from_the_two_lowest() -> None:
    costs = [1.0] * 64
    costs[10], costs[40] = 0.2, 0.1  # the two lowest: 40 then 10
    score, cma, log = _fakes(costs, [0.5, 0.5])
    r = PR.one_try("G2b", np.eye(3), score, cma)
    assert log["scores"] == 64 and len(log["cma"]) == 2
    g = PR.grid_matrices("G2b")
    assert np.allclose(log["cma"][0][0], g[40]) and np.allclose(log["cma"][1][0], g[10])
    assert [c[1:] for c in log["cma"]] == [(0.1, 500)] * 2
    assert r["evals"] == 64 + 2 * 501 == 1066 and [e["start"] for e in r["ends"]] == [40, 10]


def test_g2b_stops_after_the_first_cma_end_at_tau_and_a_grid_point_is_never_taken() -> None:
    score, cma, log = _fakes([0.1] * 64, [0.96, 0.99])
    r = PR.one_try("G2b", np.eye(3), score, cma)
    assert len(log["cma"]) == 1 and r["hit"] == 0.96 and r["evals"] == 64 + 501
    score, cma, log = _fakes([0.1] * 64, [0.6, 0.7])  # none: the better hit ratio is reported
    r = PR.one_try("G2b", np.eye(3), score, cma)
    assert r["hit"] == 0.7 and r["evals"] == 1066


def test_g2_try_uses_its_own_generator_default_rng_1_for_every_cma() -> None:
    seen: List[Any] = []

    class Res:
        overlap_info = None
        cost = 0.3
        orientation = np.eye(3)

    class Rec:
        last_local_optimization_evals = 3001

        def evaluate_overlap(self, R: Any, v: Any, p: Any) -> Any:
            return type("I", (), {"cost": 0.5})()

    def lo(**kw: Any) -> Any:
        seen.append(kw)
        return Res()

    for form in ("G2a", "G2b"):
        seen.clear()
        PR.g2_try(form, Rec(), lo, lambda info: (0.5, 1.0), np.eye(3), "V", 0)
        assert seen and all(kw["optimizer"] == "cma" for kw in seen)
        for kw in seen:
            assert np.array_equal(kw["rng"].random(4), np.random.default_rng(1).random(4))
            assert kw["voxel_vertices"] == "V" and kw["phase_index"] == 0
    assert seen[0]["cma_sigma0_deg"] == 0.1 and seen[0]["cma_max_evals"] == 500


# ---- the wrapper between G1 and the grid search ----
class _Rec:
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


def _wrap(g1_hits: Dict[int, float], g2_hits: Dict[int, float], aligned: bool = False) -> Any:
    rec = _Rec()
    book = OD.DrawBook()
    seen: List[Any] = []
    gcalls: List[Dict[str, Any]] = []

    def grid(**kw: Any) -> str:
        gcalls.append(kw)
        seen.append(kw.get("rng"))
        rec._last_global_evals, rec._last_local_evals = 300_000, 35_000
        return "GRID"

    g = (
        OD.aligned_reconstruct_voxel(
            book,
            lambda voxel_vertices, phase_index, rng=None: grid(
                voxel_vertices=voxel_vertices, phase_index=phase_index, rng=rng
            ),
        )
        if aligned
        else grid
    )
    step = dict(
        key=7,
        voxel=42,
        used=False,
        kept=[dict(R=np.eye(3) * (i + 1), cost=0.1 * (i + 1), hit=0.5) for i in range(3)],
    )
    log: List[Dict[str, Any]] = []
    g2_calls: List[Any] = []

    def g1_end(e: Dict[str, Any], v: Any, ph: int) -> Dict[str, Any]:
        k = int(round(e["R"][0, 0])) - 1
        return dict(hit=g1_hits.get(k, 0.5), cost=e["cost"], evals=1001, R=e["R"], result=f"G1_{k}")

    def g2_end(R: Any, v: Any, ph: int) -> Dict[str, Any]:
        j = int(round(R[0, 0])) - 1  # the rescue end keeps its kept end's orientation
        g2_calls.append(j)
        return dict(hit=g2_hits.get(j, 0.5), cost=0.0, evals=3001, R=R, result=f"G2_{j}")

    w = PR.make_wrapper(
        "G2a",
        g,
        lambda: step,
        (OD.voxel_key if aligned else (lambda v: int(v))),
        g1_end,
        g2_end,
        rec,
        log,
        RR.rescue_choose,
    )
    return dict(w=w, rec=rec, log=log, gcalls=gcalls, g2=g2_calls, book=book, seen=seen)


def test_g1_rescue_returns_before_g2_and_g2_is_not_called() -> None:
    s = _wrap({1: 0.97}, {})
    assert s["w"](voxel_vertices=7, phase_index=0, rng="RNG") == "G1_1"
    assert s["g2"] == [] and s["gcalls"] == []
    assert s["rec"].last_eval_counts[2] == 2002  # G1 tried ends 0 and 1
    assert not s["log"][0]["polished"]


def test_g2_runs_from_g1s_rescue_ends_in_g1s_order_and_takes_the_first_at_tau() -> None:
    # G1 tries kept ends in cost order 0, 1, 2 (costs 0.1, 0.2, 0.3), all failing
    s = _wrap({}, {1: 0.96, 2: 0.99})
    out = s["w"](voxel_vertices=7, phase_index=0, rng="RNG")
    assert out == "G2_1" and s["g2"] == [0, 1] and s["gcalls"] == []
    assert s["rec"].last_eval_counts[2] == 3 * 1001 + 2 * 3001  # every evaluation counted
    row = s["log"][0]
    assert row["polished"] and row["polish_evals"] == 6002 and row["rescue_evals"] == 3003
    assert row["grid_evals"] == 0


def test_failed_g2_calls_the_grid_with_unchanged_arguments_and_counts_everything() -> None:
    s = _wrap({}, {})
    out = s["w"](voxel_vertices=7, phase_index=0, rng="RNG")
    assert out == "GRID" and s["g2"] == [0, 1, 2]
    assert s["gcalls"] == [dict(voxel_vertices=7, phase_index=0, rng="RNG")]
    assert s["rec"].last_eval_counts[2] == 335_000 + 3 * 1001 + 3 * 3001
    row = s["log"][0]
    assert row["attempted"] and not row["rescued"] and not row["polished"]
    assert row["grid_evals"] == 335_000 and len(row["polish_tries"]) == 3
    assert all("result" not in t for t in row["polish_tries"])


def test_g2_draws_nothing_from_the_aligned_book_and_the_grid_draw_is_that_of_arm_h() -> None:
    s = _wrap({}, {0: 0.99}, aligned=True)  # accepted: no draw at all
    assert s["w"](voxel_vertices=7, phase_index=0, rng="LIBRARY") == "G2_0"
    assert s["book"].n_issued == 0 and s["book"].counts == {} and s["seen"] == []
    s = _wrap({}, {}, aligned=True)  # failed G2: the first grid draw of the voxel, as G1 alone
    assert s["w"](voxel_vertices=7, phase_index=0, rng="LIBRARY") == "GRID"
    assert s["book"].n_issued == 1 and s["book"].counts == {(7, "grid"): 1}
    (g,) = s["seen"]
    assert np.array_equal(g.random(5), OD.make_rng(7, "grid", 0).random(5))


def test_wrapper_without_a_matching_unused_step_is_the_plain_grid_call() -> None:
    s = _wrap({}, {})
    s["w"](voxel_vertices=7, phase_index=0, rng=None)
    assert s["w"](voxel_vertices=7, phase_index=0, rng=None) == "GRID"  # step already used
    assert s["log"][1]["attempted"] is False and s["log"][1]["polish_evals"] == 0


# ---- the read guards ----------------------------------------------------------------------------
def test_guard_refuses_anything_under_passes_held2(tmp_path: Path) -> None:
    for p in (
        tmp_path / "passes_held2" / "D4_H" / "pass.json",
        tmp_path / "x" / "passes_held2",
        Path("a/b/passes_held2/before_H.npy"),
    ):
        with pytest.raises(SystemExit):
            PG.guarded(p, problems=lambda: "")
    # a link into passes_held2 is refused as well
    real = tmp_path / "passes_held2"
    real.mkdir()
    (real / "f").write_text("x")
    link = tmp_path / "innocent"
    link.symlink_to(real)
    with pytest.raises(SystemExit):
        PG.guarded(link / "f", problems=lambda: "")
    # only the explicit step-13 flag opens it
    assert PG.guarded(real / "f", step13_recorded=True, problems=lambda: "") == real / "f"


def test_guard_refuses_held2_run_folders_before_the_freeze_and_opens_them_after(
    tmp_path: Path,
) -> None:
    for name in ("H_held2", "G_held2", "S8_held2", "A_held2"):
        p = tmp_path / name / "rescue_log.json.gz"
        with pytest.raises(SystemExit):
            PG.guarded(p, problems=lambda: "polish frozen.json is not committed")
        assert PG.guarded(p, problems=lambda: "") == p  # freeze agrees: allowed
    # other samples and the held2 sample folder of ar3 are not guarded
    q = tmp_path / "H_held" / "rescue_log.json.gz"
    assert PG.guarded(q, problems=lambda: "any") == q
    r = tmp_path / "ar3" / "held2" / "gradient" / "sample.json"
    assert PG.guarded(r, problems=lambda: "any") == r


def test_default_guard_uses_the_real_freeze_check_and_fails_without_a_commit() -> None:
    # the committed repo has no polish frozen.json at this point or it agrees; either way the call
    # must be a SystemExit or the path, never an exception of another kind
    try:
        PG.guarded(Path("cache/H_held2/rescue_log.json.gz"))
    except SystemExit:
        pass


# ---- the gate rules -----------------------------------------------------------------------------
def _row(
    acc: bool, err: float = 0.01, grid_err: float = 0.01, evals: int = 5000, grid: int = 300000
):
    return dict(
        accepted=acc,
        taken_err=err if acc else None,
        wrong_accepted=bool(acc and err > 1.0),
        grid_end_err=grid_err,
        saved=grid if acc else 0,
        probe_evals=evals,
        grid_evals=grid,
        taken_angle_to_grid_end=0.1 if acc else None,
    )


def test_gate_rules_pass_and_each_clause_can_fail() -> None:
    rows = [_row(True, 0.01, evals=2000, grid=500_000)] * 14 + [_row(False, evals=5330)] * 14
    m = PGT.metrics(rows)
    assert m["accepted"] == 14 and m["net_saving"] == 14 * 500_000 - 14 * 2000 - 14 * 5330
    assert PGT.gate1_form(m, True)["passes"]
    assert not PGT.gate1_form(m, False)["g0"]
    few = PGT.metrics([_row(True, grid=500_000)] * 12 + [_row(False)] * 16)
    assert not PGT.gate1_form(few, True)["g1"]
    wrong = PGT.metrics(
        [_row(True, grid=500_000)] * 13
        + [_row(True, err=1.5, grid_err=0.02, grid=500_000)]
        + [_row(False)] * 14
    )
    assert not PGT.gate1_form(wrong, True)["g2"]
    far = PGT.metrics(
        [_row(True, grid=500_000)] * 13 + [_row(True, err=0.4, grid=500_000)] + [_row(False)] * 14
    )
    assert not PGT.gate1_form(far, True)["g3"]  # one end over 0.3
    poor = PGT.metrics(
        [_row(True, evals=15_000, grid=300_000)] * 13 + [_row(False, evals=15_000)] * 15
    )
    assert not PGT.gate1_form(poor, True)["g4"]


def test_gate_choice_larger_net_saving_tie_goes_to_g2b_none_is_none() -> None:
    def res(a: bool, b: bool, na: int, nb: int) -> Dict[str, Any]:
        return dict(
            G2a=dict(gate=dict(passes=a), metrics=dict(net_saving=na)),
            G2b=dict(gate=dict(passes=b), metrics=dict(net_saving=nb)),
        )

    assert PGT.choose(res(True, True, 5, 4)) == "G2a"
    assert PGT.choose(res(True, True, 4, 5)) == "G2b"
    assert PGT.choose(res(True, True, 5, 5)) == "G2b"
    assert PGT.choose(res(True, False, 1, 9)) == "G2a"
    assert PGT.choose(res(False, False, 1, 9)) is None


def test_gate2_rules() -> None:
    base = PGT.metrics(
        [_row(True, grid=400_000, evals=1000)] * 7 + [_row(False, grid=400_000, evals=5330)] * 14
    )
    assert PGT.gate2(base)["t1"] and PGT.gate2(base)["t3"] is True
    six = PGT.metrics([_row(True, grid=900_000)] * 6 + [_row(False)] * 15)
    assert not PGT.gate2(six)["t1"]
    low = PGT.metrics(
        [_row(True, grid=200_000, evals=1000)] * 7 + [_row(False, grid=300_000, evals=15_005)] * 14
    )
    assert not PGT.gate2(low)["t3"]


# ---- the verdict rules (polish_judge) ----
import polish_judge as PJ  # noqa: E402

_OK = dict(a=True, b=True, c=True, passes=True)
_BAD = dict(a=False, b=True, c=True, passes=False)


def _r(sample: str, tj: int, th: int, **kw: Any) -> Dict[str, Any]:
    base: Dict[str, Any] = dict(
        net_vs_h=0, med_j=0.02, med_h=0.02, extra_wrong_seeds=0, house_j=_OK, house_h=_OK
    )
    base.update(kw)
    return PJ.rules(sample, tj, th, **base)


def test_rule_1_by_sample_and_its_edges() -> None:
    assert _r("held2", 94, 100)["r1"] and not _r("held2", 94_001, 100_000)["r1"]
    assert _r("held", 95, 100)["r1"] and not _r("held", 951, 1000)["r1"]
    assert _r("mid", 1005, 1000)["r1"] and not _r("mid", 1006, 1000)["r1"]
    assert _r("grad", 1005, 1000)["r1"] and not _r("grad", 1006, 1000)["r1"]
    assert _r("big", 100, 100)["r1"] and not _r("big", 99, 100)["r1"]


def test_rules_2_3_4_5_and_the_h_fails_5_clause() -> None:
    assert _r("mid", 100, 100, net_vs_h=5)["r2"] and not _r("mid", 100, 100, net_vs_h=6)["r2"]
    assert _r("mid", 100, 100, med_j=0.021)["r3"] and not _r("mid", 100, 100, med_j=0.0211)["r3"]
    assert _r("mid", 100, 100, extra_wrong_seeds=1)["r4"]
    assert not _r("mid", 100, 100, extra_wrong_seeds=2)["r4"]
    assert not _r("mid", 100, 100, house_j=_BAD)["passes"]  # H passes (5), J fails it
    assert _r("mid", 100, 100, house_j=_BAD, house_h=_BAD)["passes"]  # H fails it too
    assert PJ.house(80, 100, 5, 0, 0.0, 0.0)["passes"]
    assert not PJ.house(81, 100, 0, 0, 0.0, 0.0)["a"]  # less than 20% fewer
    assert PJ.pct_fewer(75, 100) == 25.0 and PJ.pct_fewer(1, 0) == 0.0


def test_step_verdict_needs_every_sample() -> None:
    assert PJ.step_verdict({s: True for s in PJ.SAMPLES})
    assert not PJ.step_verdict({**{s: True for s in PJ.SAMPLES}, "held2": False})


# ---- step 14's freeze check replaces step 12's for held and held2 ----
import polish_replace_freeze as PRF  # noqa: E402


class _SF:
    @staticmethod
    def problems() -> str:
        return "icenine/config_file.py differs from the frozen hash"


def test_held_and_held2_are_refused_when_step_14s_freeze_fails() -> None:
    for s in ("held", "held2"):
        sf = _SF()
        with pytest.raises(SystemExit):
            PRF.install(s, lambda: "polish frozen.json is not committed", sf)
        assert sf.problems() != ""  # step 12's check was not switched off


def test_held_and_held2_are_read_when_step_14s_freeze_passes_and_others_untouched() -> None:
    for s in ("held", "held2"):
        sf = _SF()
        PRF.install(s, lambda: "", sf)
        assert sf.problems() == ""
    sf = _SF()
    PRF.install("mid", lambda: "any", sf)  # not a held sample: nothing happens
    assert sf.problems() != ""
