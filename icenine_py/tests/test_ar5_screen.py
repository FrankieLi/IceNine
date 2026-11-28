"""AR5 step 21 (task ar5u): screen tries and the wide flag, on hand-made maps and a fake
reconstructor. Plan "AR5 step 21", phase A step 2: a screen try costs 2 evaluations and ends at its
start; the wide flag (ratio, 1 degree, the contested part unchanged); the repair is step 20's
(imported); per-fit generators; every evaluation counted; guards on ``held2`` and ``held5``.
"""

import subprocess
import sys
from pathlib import Path
from types import SimpleNamespace
from typing import Any, List

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parent.parent
AR5 = ROOT / "scripts" / "autoresearch" / "ar5"
sys.path.insert(0, str(AR5))

import cmarepair_rules as CR  # noqa: E402
import screen_guard as SG  # noqa: E402
import screen_rules as SR  # noqa: E402

from icenine.orientation_search import CMAOptimizer  # noqa: E402

I3 = np.eye(3)[None]  # identity as the only symmetry operator


def rot_z(deg: float) -> np.ndarray:
    a = np.radians(deg)
    c, s = np.cos(a), np.sin(a)
    return np.array([[c, -s, 0.0], [s, c, 0.0], [0.0, 0.0, 1.0]])


# ---- constants, imports and stored counts ------------------------------------------------------


def test_constants_are_the_plans() -> None:
    assert SR.SCREEN_MAX_EVALS == 1 and SR.EVALS_PER_SCREEN == 2
    assert SR.FULL_MAX_EVALS == 60 and SR.PASS_SIGMA0_DEG == 0.05
    assert SR.RATIO_MAX == 0.95 and SR.WIDE_DEG == 1.0
    assert SR.REPAIR_MAX_EVALS == 120 and SR.EVALS_PER_FIT == 121
    assert SR.fit_kwargs() == dict(optimizer="cma", cma_sigma0_deg=0.05, cma_max_evals=120)
    assert SR.SAMPLES == ("held4", "held3", "held", "mid", "grad", "big")


def test_the_repair_and_the_generators_are_step_20s() -> None:
    for name in ("repair", "fit_rng", "fit_kwargs", "contested", "after_pass_state"):
        assert getattr(SR, name) is getattr(CR, name), name
    src = (AR5 / "screen_rules.py").read_text()
    assert "def repair" not in src and "def fit_rng" not in src and "def contested" not in src


def test_per_fit_generators_do_not_depend_on_other_fits() -> None:
    a = SR.fit_rng(7, None).integers(0, 2**31, 4)
    SR.fit_rng(8, 3).integers(0, 2**31, 100)  # another fit draws: nothing moves
    b = SR.fit_rng(7, None).integers(0, 2**31, 4)
    assert (a == b).all()
    assert not (SR.fit_rng(7, 0).integers(0, 2**31, 4) == a).all()


def test_stored_counts_are_consistent() -> None:
    for s in SR.SAMPLES:
        j5, j5c, s8 = SR.TOTALS[s]
        assert j5c - j5 == SR.G0_J5C_EVALS[s], s  # J5C = J5 plus step 20's repair
        stored = j5 + SR.P1_PASS_EVALS[s] + SR.P1_REPAIR_EVALS[s] - SR.P1_TOTAL[s]
        assert 0.15 < stored / j5 < 0.45, (s, stored / j5)  # the pass is 20% to 41% of J5
        assert SR.P1_PASS_EVALS[s] % 2 == 0  # every try costs 2
    chosen = dict(held4=14.1, held3=14.4, held=15.6, mid=21.3, grad=31.7, big=26.6)
    for s, want in chosen.items():
        got = SR.saving_pct(SR.TOTALS[s][1], SR.P1_TOTAL[s])
        assert abs(got - want) < 0.06, (s, got, want)
    assert SR.P1_FIXED == dict(held4=3, held3=4, held=9, mid=3, grad=2, big=3)
    assert SR.P1_BROKEN == dict(held4=2, held3=2, held=3, mid=3, grad=2, big=1)


def test_total_evals_counts_every_evaluation() -> None:
    assert SR.total_evals(1000, 300, 20, 50) == 1000 - 300 + 20 + 50
    assert SR.pass_evals_of([2, 2, 2]) == 6 and SR.pass_evals_of([]) == 0


# ---- a screen try ------------------------------------------------------------------------------


class FakeInfo:
    def __init__(self, cost: float) -> None:
        self.cost = cost


class FakeRec:
    """Mimics ``Reconstructor.local_optimization``'s CMA path: refuses a budget override under 2,
    the optimizer is the library's ``CMAOptimizer``, the final evaluation is counted."""

    def __init__(self) -> None:
        self.params = SimpleNamespace(local_optimizer="mc")
        self.last_local_optimization_evals = 0
        self.calls: List[Any] = []

    def _make_cma_optimizer(self, cost_fn: Any, vertices: Any, phase: int, rng: Any) -> Any:
        return CMAOptimizer(cost_fn, vertices, phase, rng, max_evals=1000)

    def local_optimization(
        self,
        voxel_vertices: Any,
        phase_index: int,
        initial_orientation: np.ndarray,
        rng: Any = None,
        cma_sigma0_deg: Any = None,
        cma_max_evals: Any = None,
        optimizer: Any = None,
    ) -> Any:
        assert self.params.local_optimizer == "cma"
        n = [0]
        target = rot_z(2.0)

        class CF:
            def evaluate(self, R: np.ndarray, vertices: Any, phase: int) -> Any:
                n[0] += 1
                return FakeInfo(float(np.abs(R - target).sum()))

        cf = CF()
        opt = self._make_cma_optimizer(cf, voxel_vertices, phase_index, rng)
        if cma_sigma0_deg is not None:
            opt.sigma0_deg = float(cma_sigma0_deg)
        if cma_max_evals is not None:
            if cma_max_evals < 2:
                raise ValueError("cma_max_evals must be >= 2")
            opt.max_evals = int(cma_max_evals)
        self.calls.append(opt.max_evals)
        res = opt.optimize(initial_orientation)
        info = cf.evaluate(res.orientation, voxel_vertices, phase_index)  # the final evaluation
        res.overlap_info = info
        res.cost = info.cost
        self.last_local_optimization_evals = n[0]
        return res


def make_env() -> Any:
    rec = FakeRec()
    hooks = SimpleNamespace(
        rec=rec,
        mic=SimpleNamespace(voxels=[SimpleNamespace(phase=0)]),
        vertices=lambda v: np.zeros((3, 3)),
    )
    bfs = SimpleNamespace(_ratios=lambda info: (0.5, 0.25))
    return hooks, bfs, rec


def test_a_screen_try_costs_two_evaluations_and_ends_at_its_start() -> None:
    hooks, bfs, rec = make_env()
    R0 = rot_z(1.3)
    t = SR.pass_try(hooks, bfs, 0, R0, np.random.default_rng(0), SR.SCREEN_MAX_EVALS)
    assert t["evals"] == 2 == SR.EVALS_PER_SCREEN
    assert np.array_equal(t["end"], R0)  # the end is the start, bit for bit
    assert t["cost"] == float(np.abs(R0 - rot_z(2.0)).sum())
    assert rec.calls == [1]
    assert rec.params.local_optimizer == "mc"  # restored
    assert "_make_cma_optimizer" not in vars(rec)  # the patch is gone: the class method is back
    with pytest.raises(ValueError):  # the library's own override still refuses a budget of 1
        rec.params.local_optimizer = "cma"
        rec.local_optimization(np.zeros((3, 3)), 0, R0, np.random.default_rng(0), cma_max_evals=1)


def test_a_full_try_is_the_pass_s_own_budget() -> None:
    hooks, bfs, rec = make_env()
    t = SR.pass_try(hooks, bfs, 0, rot_z(1.3), np.random.default_rng(0), SR.FULL_MAX_EVALS)
    assert t["evals"] == 61  # 60 with the start's own, plus the final evaluation
    assert rec.calls == [60] and t["cost"] <= float(np.abs(rot_z(1.3) - rot_z(2.0)).sum())


def test_a_screen_try_draws_nothing_that_matters() -> None:
    hooks, bfs, _rec = make_env()
    a = SR.pass_try(hooks, bfs, 0, rot_z(1.3), np.random.default_rng(1), 1)
    b = SR.pass_try(hooks, bfs, 0, rot_z(1.3), np.random.default_rng(99), 1)
    assert np.array_equal(a["end"], b["end"]) and a["cost"] == b["cost"]


# ---- the wide flag -----------------------------------------------------------------------------


def hand_state() -> Any:
    """Six voxels. Kept ends: all the identity. Try rows: (voxel, end).
    v0: ratio 0.80, a try 5 degrees away      -> wide (extra) if not contested
    v1: ratio 0.80, a try 0.5 degree away     -> not wide
    v2: ratio 0.99, a try 5 degrees away      -> not wide (ratio not under 0.95)
    v3: ratio 0.80, contested already         -> wide, not extra
    v4: ratio 0.80, no try                    -> not wide
    v5: ratio 0.80, two tries: 0.5 and 3 deg  -> wide (extra)"""
    n = 6
    R = np.tile(np.eye(3), (n, 1, 1))
    ratio = np.array([0.80, 0.80, 0.99, 0.80, 0.80, 0.80])
    flag = np.array([False, False, False, True, False, False])
    t_voxel = np.array([0, 1, 2, 3, 5, 5])
    t_end = np.stack([rot_z(5.0), rot_z(0.5), rot_z(5.0), rot_z(5.0), rot_z(0.5), rot_z(3.0)])
    return R, ratio, flag, t_voxel, t_end


def test_wide_flag_on_hand_made_maps() -> None:
    R, ratio, flag, tv, te = hand_state()
    wide, extra = SR.wide_flag(flag, ratio, tv, te, R, I3)
    assert wide.tolist() == [True, False, False, True, False, True]
    assert extra.tolist() == [True, False, False, False, False, True]


def test_wide_flag_keeps_the_contested_part_unchanged() -> None:
    R, ratio, flag, tv, te = hand_state()
    wide, extra = SR.wide_flag(flag, ratio, tv, te, R, I3)
    assert (wide[flag]).all() and not extra[flag].any()
    assert ((wide & ~extra) == flag).all()  # wide minus the extra voxels is the contested set
    w2, _e = SR.wide_flag(flag, ratio, np.array([], dtype=int), te[:0], R, I3)
    assert (w2 == flag).all()  # no candidates: the flag is the contested one


def test_wide_flag_threshold_is_one_degree_strictly() -> None:
    R = np.tile(np.eye(3), (3, 1, 1))
    ratio = np.array([0.80, 0.80, 0.95])
    flag = np.zeros(3, dtype=bool)
    tv = np.array([0, 1, 2])
    te = np.stack([rot_z(1.0 + 1e-6), rot_z(0.999), rot_z(5.0)])
    wide, _ = SR.wide_flag(flag, ratio, tv, te, R, I3)
    assert wide.tolist() == [True, False, False]  # just over 1 degree; just under; ratio == 0.95


# ---- the repair counts every evaluation --------------------------------------------------------


def test_the_repair_counts_every_evaluation_and_per_fit_rngs() -> None:
    n = 3
    R = np.tile(np.eye(3), (n, 1, 1))
    R[1] = rot_z(10.0)
    cost = np.array([5.0, 5.0, 5.0])
    hit = np.array([0.5, 0.5, np.nan])
    seen: List[Any] = []

    def fit(v: int, R0: np.ndarray, rng: Any) -> Any:
        seen.append((v, int(rng.integers(0, 2**31))))
        return R0, 4.0, 0.6, SR.EVALS_PER_FIT

    def score(v: int, Rv: np.ndarray) -> Any:
        return 5.0, 0.5

    res = SR.repair(R, cost, hit, [0, 2], {0: [1], 2: [0, 1]}, fit, score=score)
    assert res.n_evals_fits == 5 * SR.EVALS_PER_FIT  # 2 + 3 fits
    assert res.n_evals_score == 1  # voxel 2 had no stored hit ratio
    assert res.n_evals == 5 * 121 + 1
    # each fit's draw is that of default_rng([0, v, c + 1]) whatever ran before
    for v, draw in seen:
        assert isinstance(draw, int)
    assert seen[0] == (0, int(np.random.default_rng([0, 0, 0]).integers(0, 2**31)))


# ---- gate and judge ----------------------------------------------------------------------------


def test_gate_g2_rules() -> None:
    ok = SR.gate_g2("held", 80, 100, 9, 3, 0.01, 0.011, 0)
    assert ok["ok"] and ok["ok_saving"] and ok["ok_net"] and ok["ok_median"]
    assert not SR.gate_g2("held", 91, 100, 9, 3, 0.01, 0.011, 0)["ok_saving"]  # 9% < 10%
    assert not SR.gate_g2("held", 80, 100, 2, 3, 0.01, 0.011, 0)["ok_net"]  # +1 > 0
    assert SR.gate_g2("held", 80, 100, 3, 3, 0.01, 0.011, 0)["ok_net"]  # 0 <= 0
    assert not SR.gate_g2("held", 80, 100, 9, 3, 0.012, 0.011, 0)["ok_median"]
    m = SR.gate_g2("mid", 80, 100, 3, 3, 0.01, 0.011, 2)
    assert not m["ok_s8"] and not m["ok"]
    assert SR.gate_g2("mid", 80, 100, 3, 3, 0.01, 0.011, 1)["ok"]


def test_judge_rules() -> None:
    j = SR.judge(900, 1000, 3000, 2, 5, 0.0105, 0.0100, 0.0100)
    assert j["ok"] and j["saving_vs_j5c_pct"] == 10.0
    assert not SR.judge(930, 1000, 3000, 0, 0, 0.01, 0.01, 0.01)["ok1"]  # 7%
    assert not SR.judge(900, 1000, 3000, 3, 0, 0.01, 0.01, 0.01)["ok2"]
    assert not SR.judge(900, 1000, 3000, 0, 0, 0.0111, 0.0100, 0.0100)["ok3"]
    assert not SR.judge(900, 1000, 1100, 0, 0, 0.01, 0.01, 0.01)["ok4"]  # 18% fewer than S8
    assert not SR.judge(900, 1000, 3000, 0, 6, 0.01, 0.01, 0.01)["ok4"]


# ---- guards ------------------------------------------------------------------------------------


def test_guard_refuses_held2_and_unknown_samples() -> None:
    for name in ("held2", "xheld2", "held9"):
        with pytest.raises(SystemExit) as e:
            SG.refuse_sample(name)
        assert e.value.code == 2
    for name in ("", "big", "grad", "mid", "held", "held3", "held4"):
        SG.refuse_sample(name)


def test_guard_refuses_held5_until_the_gate_file_exists(tmp_path: Path) -> None:
    with pytest.raises(SystemExit):
        SG.refuse_sample("held5", gate=tmp_path / "gate.json")
    (tmp_path / "gate.json").write_text("{}")
    SG.refuse_sample("held5", gate=tmp_path / "gate.json")


def test_guard_paths(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    cache = tmp_path / "cache"
    (cache / "logs").mkdir(parents=True)
    monkeypatch.delenv("SCREEN_PRODUCER", raising=False)
    for bad in (cache / "passes_held2" / "x.npy", cache / "J5_held2", tmp_path / "held2_clean.npy"):
        with pytest.raises(SystemExit):
            SG.guarded(bad, cache)
    SG.guarded(cache / "passes_held4" / "after_J5.npy", cache)  # held4 is disclosed now
    SG.guarded(cache / "passes_mid" / "after_J5.npy", cache)
    h5 = cache / "J5_held5" / "run.json"
    with pytest.raises(SystemExit):
        SG.guarded(h5, cache)  # before the marker
    with pytest.raises(SystemExit):
        SG.guarded(cache / "passes_held5" / "after_S8.npy", cache)
    # a producer job opens its own inputs by name, never a log or a result
    monkeypatch.setenv("SCREEN_PRODUCER", "1")
    SG.guarded(cache / "passes_held5" / "before_J5.npy", cache)
    SG.guarded(cache / "passes_held5" / "D4_J5" / "tries.npz", cache)
    with pytest.raises(SystemExit):
        SG.guarded(cache / "J5_held5" / "run.json", cache)
    with pytest.raises(SystemExit):
        SG.guarded(cache / "logs" / "J5_held5.log", cache)
    with pytest.raises(SystemExit):
        SG.guarded(cache / "passes_held5" / "after_S8.npy", cache)
    with pytest.raises(SystemExit):
        SG.guarded(cache / "passes_held2" / "before_J5.npy", cache)  # held2 never
    # after the marker everything of held5 opens
    (cache / "logs" / "screen_runs.done").write_text("exit=0\n")
    monkeypatch.delenv("SCREEN_PRODUCER")
    SG.guarded(h5, cache)
    with pytest.raises(SystemExit):
        SG.guarded(cache / "passes_held2" / "x.npy", cache)


def test_the_run_all_script_parses() -> None:
    p = AR5 / "screen_run_all.sh"
    if not p.exists():
        pytest.skip("screen_run_all.sh not written yet")
    r = subprocess.run(["bash", "-n", str(p)], capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
