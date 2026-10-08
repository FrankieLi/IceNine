"""scripts/finisher_bench/optimizers.py on a smooth synthetic SO(3) quadratic: each local finisher
converges, stops exactly at its evaluation budget, and is deterministic given a seed."""

import math
import sys
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest
from scipy.spatial.transform import Rotation

ROOT = Path(__file__).parent.parent
for sub in ("scripts/finisher_bench", "scripts/common"):
    sys.path.insert(0, str(ROOT / sub))

import optimizers as O  # noqa: E402

TARGET = Rotation.from_rotvec(np.radians([0.12, -0.2, 0.15])).as_matrix()
START = np.eye(3)
START_DEG = float(np.degrees(np.linalg.norm(Rotation.from_matrix(TARGET).as_rotvec())))
BOX = math.radians(0.3292)
STEP0 = math.radians(0.1317)
BOX_A = 1.5 * math.radians(0.25)  # the April MC's box for a told r of 0.25 deg


class Quadratic:
    """cost = (angle to TARGET in degrees)^2: smooth, unique minimum."""

    def evaluate(self, R, vertices, phase):
        ang = np.degrees(np.linalg.norm(Rotation.from_matrix(np.asarray(R) @ TARGET.T).as_rotvec()))
        return SimpleNamespace(cost=float(ang**2))


def err_deg(R: np.ndarray) -> float:
    return float(np.degrees(np.linalg.norm(Rotation.from_matrix(R @ TARGET.T).as_rotvec())))


def methods(seed: int):
    return {
        "mc_local": lambda cc: O.mc_local_restarts(
            cc, START, np.random.default_rng(seed), STEP0, 31
        ),
        "es_box": lambda cc: O.one_plus_one_es(cc, START, np.random.default_rng(seed), STEP0),
        "es_small": lambda cc: O.one_plus_one_es(
            cc, START, np.random.default_rng(seed), math.radians(0.02)
        ),
        "nm": lambda cc: O.nelder_mead_rot(cc, START, 0.13),
        "cma_005": lambda cc: O.cma_local(cc, START, seed, 0.05),
        "cma_02": lambda cc: O.cma_local(cc, START, seed, 0.2),
        "vm_small": lambda cc: O.variance_min_small_box(cc, START, seed, BOX / 4),
        "mc_deployed": lambda cc: O.mc_plain(cc, START, seed, BOX, STEP0, 200, 2, 1e-4),
        "mc_april": lambda cc: O.mc_plain(cc, START, seed, BOX_A, 0.5 * BOX_A, 3500, 2, 0.0),
    }  # fmt: skip


def run(name: str, seed: int, budget: int, ckpts=(50, 137)):
    cc = O.CountingCost(Quadratic(), None, 0, budget, ckpts)
    return O.run_budgeted(methods(seed)[name], cc)


# Tolerances (deg). ES, NM and CMA must reach the optimum. mc_local keeps the deployed rule that
# halves the step at every global improvement, so on a smooth cost it can lock in above the
# optimum (here 0.06 deg from a 0.27 deg start): only an improvement is required of it, and of
# the small-box VM (which moves inside a box of 0.08 deg around its start).
CONVERGE = {"mc_local": 0.1, "es_box": 0.01, "es_small": 0.01, "nm": 0.01, "cma_005": 0.01,
            "cma_02": 0.01, "vm_small": 0.15}  # fmt: skip


@pytest.mark.parametrize("name", sorted(CONVERGE))
def test_converges_on_smooth_quadratic(name):
    cc = run(name, 3, 4000)
    assert START_DEG > 0.2
    assert err_deg(cc.best_R) < CONVERGE[name] < START_DEG
    assert cc.best_cost == pytest.approx(Quadratic().evaluate(cc.best_R, None, 0).cost)


@pytest.mark.parametrize("name", ["mc_deployed", "mc_april"])
def test_plain_mc_improves(name):
    cc = run(name, 3, 4000)
    assert err_deg(cc.best_R) < START_DEG


@pytest.mark.parametrize("name", sorted(methods(0)))
def test_respects_budget_exactly(name):
    # the plain MC runs end by themselves before 4000; give them less than they would use
    budget = 137
    cc = run(name, 5, budget)
    assert cc.n == budget and cc.exhausted
    # one more evaluation is refused and not counted
    with pytest.raises(O.BudgetExhausted):
        cc.evaluate(np.eye(3))
    assert cc.n == budget
    # the checkpoint at the budget equals the final best, the earlier one is a prefix best
    c137, R137, n137 = cc.at(137)
    assert n137 == 137 and c137 == cc.best_cost and np.array_equal(R137, cc.best_R)
    c50, _, n50 = cc.at(50)
    assert n50 == 50 and c50 >= c137


def test_budget_longer_than_natural_run_reports_evals_used():
    cc = run("mc_deployed", 1, 10000, ckpts=(250, 10000))
    assert not cc.exhausted and cc.n < 10000
    c, R, used = cc.at(10000)
    assert used == cc.n and c == cc.best_cost


@pytest.mark.parametrize("name", sorted(methods(0)))
def test_deterministic_given_seed(name):
    a, b = run(name, 11, 300), run(name, 11, 300)
    assert a.n == b.n and np.array_equal(a.best_R, b.best_R) and a.best_cost == b.best_cost
    if name not in ("nm",):  # NM has no random draws
        c = run(name, 12, 300)
        assert not np.array_equal(a.best_R, c.best_R)


@pytest.mark.parametrize("name", ["mc_local", "es_box", "es_small", "nm", "cma_005", "vm_small"])
def test_nested_truncation(name):
    """A run with a larger budget passes through the same states: the best after 100 evaluations
    equals the result of a run with budget 100 (the methods do not look at their budget)."""
    long = run(name, 7, 600, ckpts=(100,))
    short = run(name, 7, 100, ckpts=(100,))
    assert short.n == 100
    assert long.snap[100][0] == short.best_cost
    assert np.array_equal(long.snap[100][1], short.best_R)


def test_es_step_adapts_from_too_small_and_too_large_starts():
    """The 1/5th-type rule grows a step that is far too small and shrinks one that is too large,
    so both starts reach the optimum within 600 evaluations."""
    for s0 in (math.radians(2.0), math.radians(0.0005)):
        cc = O.CountingCost(Quadratic(), None, 0, 600)
        O.run_budgeted(
            lambda c, s0=s0: O.one_plus_one_es(c, START, np.random.default_rng(0), s0), cc
        )
        assert err_deg(cc.best_R) < 0.05


# ---------------------------------------------------------------------------
# summary.py helpers
# ---------------------------------------------------------------------------


def _summary():
    return __import__("finisher_summary")


def test_sign_tests_exclude_ties_and_cluster_by_voxel():
    S = _summary()
    d = np.array([-0.01, -0.01, -0.01, 0.01, 0.0005, -0.0005])  # 3 better, 1 worse, 2 ties
    t = S.sign_test(d)
    assert (t["better"], t["worse"], t["ties"]) == (3, 1, 2)
    assert t["p"] == pytest.approx(0.625)  # two-sided exact binomial, 3 of 4
    # voxel-clustered: 9 cases in 3 voxels; per-voxel medians decide
    d = np.array([-0.1, -0.1, 0.5, -0.1, 0.0, 0.0, 0.2, 0.2, 0.2])
    vox = np.repeat([1, 2, 3], 3)
    v = S.voxel_sign_test(d, vox)
    assert v["n_vox"] == 3 and (v["better"], v["worse"], v["ties"]) == (1, 1, 1)


def test_angles_and_frac_are_consistent():
    S = _summary()
    R = Rotation.from_rotvec(np.radians([0.0, 0.0, 0.5])).as_matrix()
    assert S.angles(R, np.eye(3)) == pytest.approx(0.5)
    f = S.frac(3, 10)
    assert f["k"] == 3 and 0.0 < f["lo"] < 0.3 < f["hi"] < 1.0
