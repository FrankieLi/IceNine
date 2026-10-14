"""Lockstep quick MC: every candidate of a search level runs its short MC together (S3, option B).

At each level of ``AdaptiveVoxelReconstructor.reconstruct_voxel`` the quick MC runs
``MCOptimizer.optimize(max_mc_steps=10, max_restarts=5)`` on every candidate, one candidate at a
time, each evaluation one ``VoxelCostFunction.evaluate`` call. With the config key
``BatchedQuickMC`` (opt-in, default off) the candidates run in lockstep instead: step k of every
candidate is scored by one ``BatchedCostEvaluator.evaluate_many`` call, so the per-call overhead is
paid once per step and not once per candidate.

**Random draws (option B, a fixed slice per candidate).** Run sequentially, candidate i's 30 trial
uniforms start at stream position 30 i + 3 F_i, F_i the number of failed blocks among candidates
0 .. i-1 (a failed block draws 3 more uniforms for a restart it never uses), so each candidate's
draws depend on the outcomes before it. Reproducing that exactly needs speculative evaluation that
costs 21-67 times the counted evaluations (S3 Step 0), so the lockstep path draws one buffer
``rng.random(33 N)`` per level instead: candidate i owns positions 33 i .. 33 i + 32 (the first 30
are its 10 trials, 3 each in today's x, y, z order; the last 3 are the restart draw of the failure
branch, never used). The generator advances by exactly 33 N per level whatever the outcomes. The
buffer is drawn as ``rng.uniform(-r, r, size=33 N)``, r the trial radius: it consumes the same 33 N
doubles as ``rng.random(33 N)`` and maps each exactly as the scalar ``rng.uniform(-r, r)`` of the
sequential path does (bit for bit; recomputing ``low + (high - low) * u`` with numpy ufuncs is
not the same: numpy's C ``random_uniform`` rounds that expression once (a fused multiply-add) on
this platform, and the two roundings differ by one ulp in about half of the draws). Results
therefore differ from the sequential path's (later levels and FindOptimal also see a different
stream): the key is compared with wrong counts (an A/B), not equality. Given the slices, the result
of each candidate equals what ``MCOptimizer.optimize`` returns when its generator replays that
candidate's slice (tests/test_batched_quick_mc.py).

Only the single-block case is covered (``n_min_ergodic >= max_mc_steps``, the case of every shipped
config, ``MCRadiusScaleFactor`` 0.4: one block of 10 trial steps). Evaluations per candidate: 12,
as sequential (the start as a matrix, the start again after the quaternion round trip, 10 trials).
"""

import math
from typing import Callable, Dict, List, Optional, Sequence

import numpy as np

from .cost_functions import OverlapInfo
from .orientation_search import MCOptimizer, SearchCandidate
from .sampling import (
    QuaternionGrid,
    _quat_multiply,
    matrix_to_quaternion,
    quaternion_to_matrix,
)


def n_min_ergodic(angular_box_side: float, angular_step: float, max_mc_steps: int) -> int:
    """Trial steps per block of ``MCOptimizer.optimize`` (C++ nMinErgodicSteps)."""
    return int(2.0 * (angular_box_side / angular_step) ** 3) if angular_step > 0 else max_mc_steps


def is_single_block(angular_box_side: float, angular_step: float, max_mc_steps: int) -> bool:
    """True when ``optimize`` runs one block of ``max_mc_steps`` trial steps."""
    return n_min_ergodic(angular_box_side, angular_step, max_mc_steps) >= max_mc_steps


def lockstep_applicable(
    optimizer: object, angular_box_side: float, angular_step: float, max_mc_steps: int
) -> bool:
    """Whether the quick MC of this level can run in lockstep: a plain ``MCOptimizer`` (a study
    subclass such as ``LoggedMC`` overrides ``_mc_block`` or ``_on_block_end``) and a single block.
    """
    return type(optimizer) is MCOptimizer and is_single_block(
        angular_box_side, angular_step, max_mc_steps
    )


def uniforms_per_candidate(max_mc_steps: int) -> int:
    """Stream positions owned by one candidate: 3 per trial step plus the unused restart's 3 (33
    for 10 steps)."""
    return 3 * (max_mc_steps + 1)


def trial_radius(angular_step: float) -> float:
    """Half-width r of the trial box, x, y, z ~ U(-r, r) (as ``MCOptimizer._mc_block``)."""
    return math.tan(angular_step) / math.sqrt(12.0) if angular_step > 0 else 0.01


def draw_trial_offsets(
    rng: np.random.Generator, n_candidates: int, n_steps: int, radius: float
) -> np.ndarray:
    """The level's draws (option B): one buffer of ``33 N`` uniforms U(-r, r), shape
    (N, n_steps, 3) for the trial offsets x, y, z of every step of every candidate. Candidate i
    owns positions 33 i .. 33 i + 32; the last 3 (the failure branch's unused restart) are drawn
    and dropped. Equal, bit for bit, to scalar ``rng.uniform(-r, r)`` draws on the same stream."""
    per = uniforms_per_candidate(n_steps)
    buf = rng.uniform(-radius, radius, size=per * n_candidates)
    return buf.reshape(n_candidates, n_steps + 1, 3)[:, :n_steps, :]


def run_batched_quick_mc(
    candidates: Sequence[SearchCandidate],
    evaluate_many: Callable[[np.ndarray], List[OverlapInfo]],
    rng: np.random.Generator,
    angular_box_side: float,
    angular_step: float,
    max_mc_steps: int = 10,
    max_restarts: int = 5,
    max_convergence_cost: float = 0.0,
    grid: Optional[QuaternionGrid] = None,
) -> List[Dict[str, float]]:
    """Run ``MCOptimizer.optimize`` on every candidate in lockstep (draw option B).

    Args:
        candidates: the level's candidates; ``orientation``, ``cost`` and ``overlap_info`` are
            replaced in place by the result, as the reconstructor does with the sequential result.
        evaluate_many: scores a stack of orientations (B, 3, 3) and returns one ``OverlapInfo``
            each (``BatchedCostEvaluator.evaluate_many`` with the voxel and phase bound, or a loop
            of ``VoxelCostFunction.evaluate``); it counts its own evaluations.
        rng: the generator; advances by exactly ``33 N`` draws (N = len(candidates)).
        angular_box_side, angular_step, max_mc_steps, max_restarts, max_convergence_cost: the
            arguments of ``MCOptimizer.optimize``.
        grid: a ``QuaternionGrid`` (default: a new one).

    Returns:
        The ``MCOptimizer.last_run`` dict of each candidate. Raises ValueError for a run with
        more than one block (``is_single_block`` false).
    """
    N = len(candidates)
    if N == 0:
        return []
    n_min = n_min_ergodic(angular_box_side, angular_step, max_mc_steps)
    if n_min < max_mc_steps:
        raise ValueError(
            "the lockstep quick MC covers a single block only (n_min_ergodic >= steps)"
        )
    n_steps = max_mc_steps
    gen = grid if grid is not None else QuaternionGrid()

    # the level's draws: each candidate's trial offsets (x, y, z per step), (N, n_steps, 3)
    offsets = draw_trial_offsets(rng, N, n_steps, trial_radius(angular_step))

    init_R = [np.asarray(c.orientation) for c in candidates]
    init_q = [matrix_to_quaternion(R) for R in init_R]
    round_trip_R = [quaternion_to_matrix(q) for q in init_q]
    # both start evaluations of every candidate in one call: the matrix, then the round trip
    infos = evaluate_many(np.stack(init_R + round_trip_R))
    info_start, info_round = infos[:N], infos[N:]
    cost_start = [i.cost for i in info_start]

    opt_q = [q.copy() for q in init_q]  # best state of each candidate's block
    opt_R = list(round_trip_R)
    opt_info = list(info_round)
    opt_cost = [i.cost for i in info_round]
    for t in range(n_steps):
        trial_q = []
        trial_R = []
        for i in range(N):
            x, y, z = offsets[i, t]
            delta_q = gen.get_near_identity_point(x, y, z)
            q = _quat_multiply(delta_q, opt_q[i])
            trial_q.append(q)
            trial_R.append(quaternion_to_matrix(q))
        trial_infos = evaluate_many(np.stack(trial_R))
        for i in range(N):
            if trial_infos[i].cost < opt_cost[i]:  # strict improvement, as _mc_block
                opt_q[i] = trial_q[i]
                opt_R[i] = trial_R[i]
                opt_cost[i] = trial_infos[i].cost
                opt_info[i] = trial_infos[i]

    runs: List[Dict[str, float]] = []
    for i, cand in enumerate(candidates):
        failed = opt_cost[i] >= cost_start[i]  # the block did not beat the start: a failure
        if failed:
            cand.orientation, cand.cost, cand.overlap_info = (
                round_trip_R[i],
                cost_start[i],
                info_start[i],
            )
            cur_step, n_succ, n_accept, last_accept = angular_step, 1, 0, -1
        else:
            cand.orientation, cand.cost, cand.overlap_info = opt_R[i], opt_cost[i], opt_info[i]
            cur_step, n_succ, n_accept, last_accept = angular_step * 0.5, 0, 1, n_steps
        if cand.cost < max_convergence_cost:
            stop = 2
        elif n_succ > max_restarts:
            stop = 1
        else:
            stop = 0
        runs.append(
            dict(
                stop=stop,
                steps_run=n_steps,
                n_accept=n_accept,
                last_accept=last_accept,
                n_restarts=1 - n_accept,
                n_blocks=1,
                final_step_deg=math.degrees(cur_step),
                min_ergodic=n_min,
                since_improve=n_succ,
                cost_start=cost_start[i],
                cost_end=cand.cost,
            )
        )
    return runs
