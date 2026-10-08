"""Local finishers for the Phase B3 head-to-head (pure: no pixel data, no physics).

Every optimizer takes a cost object `f` with `f.evaluate(R, vertices, phase)` returning an object
with `.cost` (the CountingCost below wraps VoxelCostFunction; the unit tests use a synthetic SO(3)
quadratic) and a start rotation `R0`. They evaluate the start first and return nothing: the result
is read from the CountingCost (the lowest-cost orientation it evaluated), so a budget stop and a
natural stop are handled the same way.

Conventions:
  * a rotation is a 3x3 matrix; a perturbation is applied on the left (delta * R), as in
    MCOptimizer;
  * ES and MC-local draw proposals exactly as MCOptimizer.optimize does at step s (radius
    tan(s)/sqrt(12); x, y, z ~ U(-radius, radius); delta = get_near_identity_point(x, y, z)),
    which has a median rotation angle of about 0.97 s;
  * NM and CMA work in a rotation vector (degrees) about the centre: R = exp(v) R_centre.
"""

import math
from typing import Any, Callable, Dict, List, Optional, Sequence

import numpy as np
from scipy.spatial.transform import Rotation

from icenine.orientation_search import (
    QuaternionGrid,
    _quat_multiply,
    matrix_to_quaternion,
    quaternion_to_matrix,
)


class BudgetExhausted(Exception):
    """Raised by CountingCost.evaluate when the next evaluation would exceed the budget."""


class CountingCost:
    """Counts evaluations exactly, tracks the best (lowest-cost, first on ties) orientation
    evaluated, snapshots the best at checkpoint counts, and stops the method at `budget`.

    evaluate() has the signature of VoxelCostFunction.evaluate so that MCOptimizer (which calls
    self.cost_fn.evaluate(R, vertices, phase)) can use it; the vertices / phase are those given at
    construction and the call's own are ignored.
    """

    def __init__(
        self,
        fn: Any,
        vertices: Any,
        phase: int,
        budget: int,
        checkpoints: Sequence[int] = (),
    ) -> None:
        self._fn, self._vertices, self._phase = fn, vertices, phase
        self.budget = int(budget)
        self.checkpoints = sorted(int(c) for c in checkpoints if int(c) <= self.budget)
        self._ck = set(self.checkpoints)
        self.n = 0
        self.best_cost = math.inf
        self.best_R: Optional[np.ndarray] = None
        self.snap: Dict[int, Any] = {}  # checkpoint -> (cost, R) of the best at that count
        self.exhausted = False

    def evaluate(self, orientation: np.ndarray, voxel_vertices: Any = None, phase_index: int = 0):
        if self.n >= self.budget:
            self.exhausted = True
            raise BudgetExhausted
        R = np.asarray(orientation, dtype=np.float64)
        info = self._fn.evaluate(R, self._vertices, self._phase)
        self.n += 1
        c = float(info.cost)
        if c < self.best_cost:
            self.best_cost, self.best_R = c, R.copy()
        if self.n in self._ck:
            self.snap[self.n] = (self.best_cost, self.best_R.copy())
        return info

    def __call__(self, R: np.ndarray) -> float:
        return float(self.evaluate(R).cost)

    def at(self, checkpoint: int):
        """(cost, R, evals used) of the method at `checkpoint`: the best after that many
        evaluations, or the final best when the method stopped earlier."""
        if checkpoint in self.snap:
            return (*self.snap[checkpoint], checkpoint)
        assert self.n < checkpoint or checkpoint > self.budget
        assert self.best_R is not None
        return self.best_cost, self.best_R.copy(), self.n


def run_budgeted(method: Callable[[CountingCost], None], cc: CountingCost) -> CountingCost:
    """Run `method(cc)` until it returns or the budget is reached."""
    try:
        method(cc)
    except BudgetExhausted:
        pass
    return cc


# ---------------------------------------------------------------------------
# Shared pieces
# ---------------------------------------------------------------------------


def _propose(
    q: np.ndarray, step_rad: float, rng: np.random.Generator, grid: QuaternionGrid
) -> np.ndarray:
    """One MCOptimizer proposal quaternion at step step_rad around quaternion q."""
    radius = math.tan(step_rad) / math.sqrt(12.0) if step_rad > 0 else 0.01
    x = rng.uniform(-radius, radius)
    y = rng.uniform(-radius, radius)
    z = rng.uniform(-radius, radius)
    return _quat_multiply(grid.get_near_identity_point(x, y, z), q)


def rot_about(R_centre: np.ndarray, v_deg: np.ndarray) -> np.ndarray:
    """exp(v) R_centre for a rotation vector v in degrees."""
    return Rotation.from_rotvec(np.radians(np.asarray(v_deg, dtype=np.float64))).as_matrix() @ (
        np.asarray(R_centre, dtype=np.float64)
    )


# ---------------------------------------------------------------------------
# (ii) MC with local restarts
# ---------------------------------------------------------------------------


def mc_local_restarts(
    f: CountingCost,
    R0: np.ndarray,
    rng: np.random.Generator,
    step0_rad: float,
    stuck_steps: int,
    restart_shrink: float = 0.5,
    min_step_rad: float = math.radians(1e-4),
) -> None:
    """MCOptimizer's zero-temperature greedy walk with two changes to the restart rule:
      * a restart is local: the walker moves to a proposal around the best at the *current* step
        (not a uniform point in the box), and the step is multiplied by restart_shrink (it is not
        reset to step0);
      * "stuck" is a fixed number of steps without a global improvement (stuck_steps), not the
        deployed 2 (box/step)^3, which grows 8-fold at every step halving.
    As in the deployed MC the step still halves at every global improvement. Restarts are
    unlimited (the evaluation budget ends the run)."""
    grid = QuaternionGrid()
    best_q = matrix_to_quaternion(np.asarray(R0, dtype=np.float64))
    cur_q = best_q.copy()
    best = current = f(np.asarray(R0, dtype=np.float64))
    step = step0_rad
    n_since = 0
    while True:
        trial_q = _propose(cur_q, step, rng, grid)
        c = f(quaternion_to_matrix(trial_q))
        if c < current:
            current, cur_q = c, trial_q
            if current < best:
                best, best_q = current, cur_q.copy()
                n_since = 0
                step = max(step * 0.5, min_step_rad)
        else:
            n_since += 1
        if n_since >= stuck_steps:
            step = max(step * restart_shrink, min_step_rad)
            cur_q = _propose(best_q, step, rng, grid)
            current = f(quaternion_to_matrix(cur_q))
            n_since = 0


# ---------------------------------------------------------------------------
# (iii) (1+1)-ES with the 1/5th-type success rule
# ---------------------------------------------------------------------------


def one_plus_one_es(
    f: CountingCost,
    R0: np.ndarray,
    rng: np.random.Generator,
    step0_rad: float,
    p_target: float = 0.25,
    damping: float = 3.0,
    min_step_rad: float = math.radians(2e-4),
    max_step_rad: float = math.radians(2.0),
) -> None:
    """(1+1)-ES on SO(3): propose around the current point at step s (MCOptimizer proposal),
    accept if the cost is strictly lower, and set
        s <- s exp((1[success] - p_target) / (damping (1 - p_target))),
    the success-rate rule whose fixed point is success probability p_target (success x1.40, failure
    x0.894 at the defaults). The step is clipped to [min_step, max_step]."""
    grid = QuaternionGrid()
    q = matrix_to_quaternion(np.asarray(R0, dtype=np.float64))
    cur = f(np.asarray(R0, dtype=np.float64))
    step = step0_rad
    while True:
        trial_q = _propose(q, step, rng, grid)
        c = f(quaternion_to_matrix(trial_q))
        ok = c < cur
        if ok:
            cur, q = c, trial_q
        step *= math.exp((float(ok) - p_target) / (damping * (1.0 - p_target)))
        step = min(max(step, min_step_rad), max_step_rad)


# ---------------------------------------------------------------------------
# (iv) Nelder-Mead on the rotation vector with restarts of a shrinking simplex
# ---------------------------------------------------------------------------


def nelder_mead_rot(
    f: CountingCost,
    R0: np.ndarray,
    size0_deg: float,
    shrink: float = 0.5,
    min_size_deg: float = 0.002,
    xatol_deg: float = 1e-6,
) -> None:
    """scipy Nelder-Mead on v (degrees) with R = exp(v) R_centre, started from a right-angled
    simplex of edge size (deg). NM ends when the simplex is below xatol and all its vertices have
    the same cost; it is then restarted around the best point with an edge `shrink` times the
    previous one (floored at min_size_deg). Deterministic."""
    from scipy.optimize import minimize

    centre = np.asarray(R0, dtype=np.float64)
    f(centre)
    size = size0_deg
    while True:
        c0 = centre

        def g(v: np.ndarray, c0: np.ndarray = c0) -> float:
            return f(rot_about(c0, v))

        simplex = np.vstack([np.zeros(3), size * np.eye(3)])
        minimize(
            g,
            np.zeros(3),
            method="Nelder-Mead",
            options=dict(
                initial_simplex=simplex,
                xatol=xatol_deg,
                fatol=0.0,
                maxiter=10**9,
                maxfev=10**9,
            ),
        )
        assert f.best_R is not None
        centre = f.best_R.copy()
        size = max(size * shrink, min_size_deg)


# ---------------------------------------------------------------------------
# (v) CMA-ES (the `cma` package, already a dependency of the benchmarks extra)
# ---------------------------------------------------------------------------


def cma_local(
    f: CountingCost,
    R0: np.ndarray,
    seed: int,
    sigma0_deg: float,
    popsize: Optional[int] = None,
) -> None:
    """Local CMA-ES on v (degrees), R = exp(v) R0, sigma0 in degrees. The package's flat-fitness
    and function-value stopping tests are switched off (the cost is quantised, so a small
    population is often flat); it stops on its own x-tolerance, otherwise the budget ends it."""
    import cma

    R0 = np.asarray(R0, dtype=np.float64)
    f(R0)
    opts: Dict[str, Any] = dict(
        seed=int(seed) + 1,
        verbose=-9,
        tolfun=0.0,
        tolfunhist=0.0,
        tolflatfitness=10**9,
        tolstagnation=10**9,
        tolx=1e-9,
    )
    if popsize is not None:
        opts["popsize"] = int(popsize)
    es = cma.CMAEvolutionStrategy(np.zeros(3), float(sigma0_deg), opts)
    while not es.stop():
        X: List[np.ndarray] = es.ask()
        costs = [f(rot_about(R0, x)) for x in X]
        es.tell(X, costs)


# ---------------------------------------------------------------------------
# (vi) VarianceMinimizing with a small box
# ---------------------------------------------------------------------------


def variance_min_small_box(f: CountingCost, R0: np.ndarray, seed: int, box_rad: float) -> None:
    """MCOptimizer.variance_minimizing_optimize with search box `box_rad`, as the finisher calls it
    (restarts 2, convergence variance 0.02^2, max convergence cost 0) and an effectively unbounded
    step budget: the CountingCost ends it."""
    from icenine.orientation_search import MCOptimizer

    mc = MCOptimizer(f, None, 0, np.random.default_rng(seed))  # type: ignore[arg-type]
    mc.variance_minimizing_optimize(
        np.asarray(R0, dtype=np.float64),
        box_rad,
        10**9,
        2,
        0.0,
        0.02**2,
    )


# ---------------------------------------------------------------------------
# (i) MC as deployed / as in the April sweep
# ---------------------------------------------------------------------------


def mc_plain(
    f: CountingCost,
    R0: np.ndarray,
    seed: int,
    box_rad: float,
    step_rad: float,
    max_steps: int,
    restarts: int,
    max_convergence_cost: float,
) -> None:
    """The unmodified MCOptimizer.optimize."""
    from icenine.orientation_search import MCOptimizer

    mc = MCOptimizer(f, None, 0, np.random.default_rng(seed))  # type: ignore[arg-type]
    mc.optimize(
        np.asarray(R0, dtype=np.float64),
        box_rad,
        step_rad,
        max_steps,
        restarts,
        max_convergence_cost,
    )
