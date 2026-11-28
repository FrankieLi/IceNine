"""AR5 step 20 (task ar5t): the repair with a 120-evaluation CMA fit, on hand-made maps and a fake
fit. Plan "AR5 step 20", phase A step 2: a fit's result does not depend on which other fits ran or
how long (per-fit generators); a CMA fit counts 121 evaluations; the keep rule with the current end
first on ties; the flag is step 19's frozen function; every evaluation counted; guards on ``held2``
and ``held4``.
"""

import sys
from pathlib import Path
from typing import Any, Dict, List, Tuple

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "scripts" / "autoresearch" / "ar5"))

import cmarepair_guard as CG  # noqa: E402
import cmarepair_rules as CR  # noqa: E402
import repair_rules as R19  # noqa: E402

from icenine.orientation_search import CMAOptimizer  # noqa: E402


def rot_z(deg: float) -> np.ndarray:
    a = np.radians(deg)
    c, s = np.cos(a), np.sin(a)
    return np.array([[c, -s, 0.0], [s, c, 0.0], [0.0, 0.0, 1.0]])


# ---- constants and the imported flag -----------------------------------------------------------


def test_the_frozen_constants_are_the_plans() -> None:
    assert CR.REPAIR_OPTIMIZER == "cma" and CR.CMA_SIGMA0_DEG == 0.05
    assert CR.CMA_MAX_EVALS == 120 and CR.EVALS_PER_FIT == 121 and CR.REPAIR_SEED == 0
    assert CR.fit_kwargs() == dict(optimizer="cma", cma_sigma0_deg=0.05, cma_max_evals=120)
    assert CR.RATIO_MAX == 0.95 and CR.MARGIN_MAX == 0.15
    assert CR.FAR_DEG == 1.0 and CR.CAND_DEG == 0.5 and CR.WRONG_DEG == 1.0
    assert CR.FLAGGED == dict(held3=194, held=183, mid=142, grad=238, big=56)
    assert CR.FITS == dict(held3=1_214, held=1_026, mid=635, grad=965, big=162)
    assert CR.P1_EVALS == dict(held3=146_894, held=124_146, mid=76_835, grad=116_765, big=19_602)
    assert {s: CR.P1_EVALS[s] for s in CR.FITS} == {s: 121 * n for s, n in CR.FITS.items()}
    assert CR.P1_FIXED == dict(held3=14, held=18, mid=15, grad=31, big=14)
    assert CR.P1_J5L == dict(flagged=145, fits=650, evals=78_650, fixed=18, broken=0)
    assert 121 * CR.J5L["fits"] == CR.P1_J5L["evals"]
    assert CR.max_share("grad") == 2.0 and CR.max_share("mid") == 1.5
    assert CR.J5_TOTAL == dict(
        held3=23_503_385, held=22_939_812, mid=10_589_493, grad=9_041_686, big=3_567_098
    )


def test_the_flag_and_the_keep_rule_are_step_19s_frozen_functions() -> None:
    for name in (
        "after_pass_state",
        "other_end_hit",
        "contested",
        "repair_candidates",
        "wrong_mask",
        "fixed_broken",
    ):
        assert getattr(CR, name) is getattr(R19, name), name
    src = (ROOT / "scripts" / "autoresearch" / "ar5" / "cmarepair_rules.py").read_text()
    assert "def contested" not in src and "def after_pass_state" not in src
    assert "def ratio_of" not in src and "def margin_of" not in src


# ---- per-fit generators ------------------------------------------------------------------------


def test_per_fit_generator_is_keyed_by_voxel_and_start() -> None:
    a = CR.fit_rng(7, None).random(5)
    assert np.array_equal(a, np.random.default_rng([0, 7, 0]).random(5))
    assert np.array_equal(CR.fit_rng(7, 3).random(5), np.random.default_rng([0, 7, 4]).random(5))
    assert not np.array_equal(a, CR.fit_rng(7, 3).random(5))
    assert not np.array_equal(a, CR.fit_rng(8, None).random(5))


def make_fit(draws: Dict[Tuple[int, int], float], extra: int = 0) -> Any:
    """A fake fit: records the first draw of each fit's own generator, then burns ``extra`` more
    (a fit that runs longer). The end is the start, cost and hit 0.5 (nothing is kept)."""

    def fit(v: int, R0: np.ndarray, rng: np.random.Generator) -> Any:
        first = float(rng.random())
        rng.random(extra)
        draws[(v, int(np.rint(np.degrees(np.arctan2(R0[1, 0], R0[0, 0])))))] = first
        return R0, 0.5, 0.5, CR.EVALS_PER_FIT

    return fit


def test_a_fit_does_not_depend_on_which_other_fits_ran_or_how_long() -> None:
    n = 4
    R = np.stack([rot_z(10.0 * i) for i in range(n)])
    cost = np.full(n, 0.9)
    hit = np.full(n, 0.5)
    cands = {1: [0, 2], 2: [3]}
    full: Dict[Tuple[int, int], float] = {}
    CR.repair(R, cost, hit, [1, 2], cands, make_fit(full))
    # longer fits (more draws), another order of the flagged list, an extra flagged voxel
    longer: Dict[Tuple[int, int], float] = {}
    CR.repair(R, cost, hit, [2, 1], cands, make_fit(longer, extra=50))
    assert longer == full
    more: Dict[Tuple[int, int], float] = {}
    CR.repair(R, cost, hit, [0, 1, 2], {**cands, 0: [1]}, make_fit(more))
    for k, x in full.items():
        assert more[k] == x
    # a dropped fit (no candidates for voxel 1) leaves voxel 2's draws unchanged
    dropped: Dict[Tuple[int, int], float] = {}
    CR.repair(R, cost, hit, [1, 2], {2: [3]}, make_fit(dropped))
    assert all(dropped[k] == full[k] for k in dropped if k[0] == 2)


# ---- a CMA fit counts 121 evaluations ----------------------------------------------------------


class _Info:
    def __init__(self, cost: float) -> None:
        self.cost = cost


class _Quadratic:
    """A cost function with no zero (never converges below a cost stop) and an eval counter."""

    def __init__(self) -> None:
        self.eval_count = 0

    def evaluate(self, R: np.ndarray, vertices: Any = None, phase: int = 0) -> _Info:
        self.eval_count += 1
        a = np.degrees(np.arctan2(R[1, 0], R[0, 0]))
        return _Info(1.0 + 0.01 * (a - 0.3) ** 2)


def test_a_cma_fit_counts_121_evaluations() -> None:
    cf = _Quadratic()
    kw = CR.fit_kwargs()
    opt = CMAOptimizer(
        cost_fn=cf,
        voxel_vertices=None,
        phase_index=0,
        rng=CR.fit_rng(5, None),
        sigma0_deg=kw["cma_sigma0_deg"],
        max_evals=kw["cma_max_evals"],
    )
    res = opt.optimize(rot_z(0.0))
    assert cf.eval_count == CR.CMA_MAX_EVALS == 120 and res.n_evals == 120
    cf.evaluate(res.orientation)  # the library's final evaluation (local_optimization)
    assert cf.eval_count == CR.EVALS_PER_FIT == 121


# ---- the keep rule -----------------------------------------------------------------------------


def fit_returning(rows: List[Tuple[float, float]]) -> Any:
    it = iter(rows)

    def fit(v: int, R0: np.ndarray, rng: np.random.Generator) -> Any:
        c, h = next(it)
        return R0 @ rot_z(0.1), c, h, CR.EVALS_PER_FIT

    return fit


def test_keep_rule_current_end_first_on_ties() -> None:
    R = np.stack([rot_z(0.0), rot_z(20.0)])
    cost, hit = np.array([0.5, 0.5]), np.array([0.6, 0.6])
    # a tie on hit and cost: the current end stays
    r = CR.repair(R, cost, hit, [0], {0: [1]}, fit_returning([(0.5, 0.6), (0.5, 0.6)]))
    assert r.changed == [] and np.array_equal(r.R, R)
    # same hit, lower cost: the first such fit (the own end) wins
    r = CR.repair(R, cost, hit, [0], {0: [1]}, fit_returning([(0.4, 0.6), (0.3, 0.6)]))
    assert r.changed == [0] and r.cost[0] == 0.3 and r.hit[0] == 0.6
    # a higher hit wins over a lower cost
    r = CR.repair(R, cost, hit, [0], {0: [1]}, fit_returning([(0.1, 0.6), (0.9, 0.7)]))
    assert r.changed == [0] and r.cost[0] == 0.9 and r.hit[0] == 0.7
    # worse everywhere: unchanged
    r = CR.repair(R, cost, hit, [0], {0: [1]}, fit_returning([(0.6, 0.5), (0.7, 0.4)]))
    assert r.changed == []


def test_candidates_are_read_from_the_snapshot_not_the_updated_map() -> None:
    R = np.stack([rot_z(0.0), rot_z(20.0), rot_z(40.0)])
    cost, hit = np.full(3, 0.5), np.full(3, 0.5)
    starts: List[Tuple[int, float]] = []

    def fit(v: int, R0: np.ndarray, rng: np.random.Generator) -> Any:
        starts.append((v, float(np.degrees(np.arctan2(R0[1, 0], R0[0, 0])))))
        return rot_z(5.0), 0.1, 0.9, CR.EVALS_PER_FIT  # always better: voxel 0 moves to 5 degrees

    CR.repair(R, cost, hit, [0, 1], {0: [], 1: [0]}, fit)
    assert (1, 0.0) in starts  # voxel 1's candidate 0 is read at its snapshot orientation


# ---- every evaluation counted ------------------------------------------------------------------


def test_every_evaluation_is_counted_including_a_missing_hit() -> None:
    R = np.stack([rot_z(0.0), rot_z(20.0), rot_z(40.0)])
    cost = np.array([np.nan, 0.5, 0.5])
    hit = np.array([np.nan, 0.5, 0.5])
    scored: List[int] = []

    def score(v: int, Rv: np.ndarray) -> Tuple[float, float]:
        scored.append(v)
        return 0.5, 0.5

    r = CR.repair(
        R,
        cost,
        hit,
        [0, 1],
        {0: [1, 2], 1: [0]},
        fit_returning([(0.6, 0.4)] * 5),
        score=score,
    )
    assert scored == [0]
    assert len(r.fits) == 3 + 2
    assert r.n_evals_fits == 5 * 121 and r.n_evals_score == 1 and r.n_evals == 5 * 121 + 1
    with pytest.raises(ValueError):
        CR.repair(R, cost, hit, [0], {0: []}, fit_returning([(0.6, 0.4)]))


# ---- guards ------------------------------------------------------------------------------------


def test_guards_on_held2(tmp_path: Path) -> None:
    with pytest.raises(SystemExit) as e:
        CG.refuse_sample("held2")
    assert e.value.code == 2
    with pytest.raises(SystemExit):
        CG.refuse_sample("held5")
    with pytest.raises(SystemExit):
        CG.guarded(tmp_path / "x_held2" / "a.npy", cache=tmp_path)
    CG.refuse_sample("mid")


def test_guards_on_held4(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    gate = tmp_path / "gate.json"
    with pytest.raises(SystemExit):
        CG.refuse_sample("held4", gate=gate)  # no gate file: refused
    gate.write_text("{}")
    CG.refuse_sample("held4", gate=gate)
    cache = tmp_path / "cache"
    (cache / "logs").mkdir(parents=True)
    log = cache / "logs" / "S8_held4.log"
    mp = cache / "passes_held4" / "after_J5.npy"
    truth = cache.parent / "ar3" / "held4" / "gradient" / "sample" / "gradient_truth.mic"
    monkeypatch.delenv("CMAREPAIR_PRODUCER", raising=False)
    for p in (log, mp, truth):
        with pytest.raises(SystemExit):
            CG.guarded(p, cache=cache)
    # the repair job: its own inputs by name, never the truth, a log or a result
    monkeypatch.setenv("CMAREPAIR_PRODUCER", "1")
    CG.guarded(mp, cache=cache)
    CG.guarded(cache / "passes_held4" / "D4_J5" / "tries.npz", cache=cache)
    for p in (
        log,
        truth,
        cache / "J5CR_held4" / "after.npy",
        cache / "passes_held4" / "after_S8.npy",
    ):
        with pytest.raises(SystemExit):
            CG.guarded(p, cache=cache)
    monkeypatch.delenv("CMAREPAIR_PRODUCER")
    # a sample that is not held4 is never blocked by the marker rule
    CG.guarded(cache / "passes_mid" / "after_J5.npy", cache=cache)
    # the marker opens held4
    (cache / "logs" / "cmarepair_runs.done").write_text("exit=0\n")
    for p in (log, mp, truth):
        CG.guarded(p, cache=cache)


def test_freeze_lists_the_imported_flag_and_this_steps_scripts() -> None:
    import cmarepair_freeze as CF

    cur = CF.current()["files"]
    for s in (
        "cmarepair_rules.py",
        "cmarepair_guard.py",
        "cmarepair_data.py",
        "cmarepair_common.py",
        "cmarepair_launch.py",
        "cmarepair_link_cache.sh",
        "repair_rules.py",
    ):
        assert f"scripts/autoresearch/ar5/{s}" in cur, s
    assert CF.problems() != "" or CF.FROZEN.exists()
