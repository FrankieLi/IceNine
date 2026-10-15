"""Lockstep quick MC (S3 Step 2, draw option B): every candidate's result equals
``MCOptimizer.optimize`` replayed on that candidate's slice of the level's uniform buffer; the
generator advances by exactly 33 N; the mapped uniforms are ``Generator.uniform(-r, r)``'s."""

import copy
import math
import warnings
from pathlib import Path
from typing import Any, Callable, Dict, List, Tuple

import numpy as np
import pytest
import torch
from scipy.spatial.transform import Rotation

from icenine import batched_cost as BC
from icenine import batched_quick_mc as BQ
from icenine import cost_functions as CF
from icenine.cost_functions import OverlapInfo, VoxelCostFunction
from icenine.orientation_search import MCOptimizer, SearchCandidate
from icenine.ray_overlap import RayRule

from test_batched_cost import orientations
from test_center_screen import three_voxels  # noqa: F401  (fixture)

needs_c = pytest.mark.skipif(not CF._HAS_C_RASTERIZE, reason="C extension missing")
BOX = math.radians(0.5)
STEP = BOX * 0.4  # MCRadiusScaleFactor 0.4: n_min_ergodic = 31, one block of 10 trial steps


class ReplayRng:
    """Serves the ``uniform`` calls ``MCOptimizer.optimize`` makes from one candidate's slice of
    the level buffer (``rng.uniform(-r, r, size=33 N)``): the 30 trial draws, which must be asked
    for with the trial radius, then the 3 restart draws of a failed block (unused by the run)."""

    def __init__(self, vals: np.ndarray, radius: float) -> None:
        self.vals, self.radius, self.pos = vals, radius, 0

    def uniform(self, low: float, high: float) -> float:
        if self.pos < 30:
            assert (low, high) == (-self.radius, self.radius), self.pos
        v = float(self.vals[self.pos])
        self.pos += 1
        return v


def level_buffer(state: np.random.Generator, n: int, step: float = STEP) -> np.ndarray:
    """The values the lockstep path draws for a level, from a copy of the generator state."""
    return copy.deepcopy(state).uniform(-BQ.trial_radius(step), BQ.trial_radius(step), 33 * n)


def sequential_reference(
    cands: List[SearchCandidate],
    cost_fn: Any,
    verts: torch.Tensor,
    phase: int,
    buf: np.ndarray,
    box: float = BOX,
    step: float = STEP,
    max_conv: float = 0.0,
) -> Tuple[List[SearchCandidate], List[Dict[str, float]], List[int]]:
    """Today's ``optimize`` on each candidate, its generator replaying the candidate's slice.
    Returns (results, last_run dicts, draws consumed per candidate)."""
    per = BQ.uniforms_per_candidate(10)
    out, runs, used = [], [], []
    for i, c in enumerate(cands):
        replay = ReplayRng(buf[i * per : (i + 1) * per], BQ.trial_radius(step))
        mc = MCOptimizer(cost_fn, verts, phase, rng=replay)  # type: ignore[arg-type]
        res = mc.optimize(c.orientation, box, step, 10, 5, max_convergence_cost=max_conv)
        out.append(res)
        runs.append(mc.last_run)
        used.append(replay.pos)
    return out, runs, used


def make_candidates(c: Dict[str, Any], n: int, seed: int = 1) -> List[SearchCandidate]:
    truth = np.asarray(c["voxel"].orientation, dtype=np.float64)
    R = orientations(truth, n_each=n)[:: max(1, (5 * n + n + 1) // n)][:n]
    return [SearchCandidate(orientation=r.copy(), cost=1.0) for r in R]


def assert_same(
    got: List[SearchCandidate],
    ref: List[SearchCandidate],
    got_runs: List[Dict[str, float]],
    ref_runs: List[Dict[str, float]],
) -> None:
    assert len(got) == len(ref)
    for i, (a, b) in enumerate(zip(got, ref)):
        assert np.array_equal(a.orientation, b.orientation), i
        assert a.cost == b.cost, i
        assert a.overlap_info == b.overlap_info, i
        assert got_runs[i] == ref_runs[i], i


def run_lockstep(
    cands: List[SearchCandidate],
    evaluate_many: Callable[[np.ndarray], List[OverlapInfo]],
    rng: np.random.Generator,
    max_conv: float = 0.0,
) -> List[Dict[str, float]]:
    return BQ.run_batched_quick_mc(
        cands, evaluate_many, rng, BOX, STEP, 10, 5, max_convergence_cost=max_conv
    )


# ---------------------------------------------------------------------------------------------
# the draw rule
# ---------------------------------------------------------------------------------------------


@pytest.mark.parametrize("step", [math.radians(0.05), STEP, math.radians(2.0), 0.0])
def test_level_buffer_equals_scalar_generator_uniform(step: float) -> None:
    """The offsets of candidate i, step t are the scalar ``rng.uniform(-r, r)`` draws at stream
    positions 33 i + 3 t + (0, 1, 2), bit for bit; the buffer consumes as ``rng.random(33 N)``."""
    r = BQ.trial_radius(step)
    n = 25
    a, b, c = (np.random.default_rng(123) for _ in range(3))
    offsets = BQ.draw_trial_offsets(a, n, 10, r)
    scalar = np.array([b.uniform(-r, r) for _ in range(33 * n)]).reshape(n, 11, 3)
    assert offsets.shape == (n, 10, 3)
    assert np.array_equal(offsets, scalar[:, :10, :])  # bit for bit
    c.random(33 * n)
    assert a.bit_generator.state == b.bit_generator.state == c.bit_generator.state


def test_trial_radius_and_block_rules() -> None:
    assert BQ.trial_radius(STEP) == math.tan(STEP) / math.sqrt(12.0)
    assert BQ.trial_radius(0.0) == 0.01
    assert BQ.n_min_ergodic(BOX, STEP, 10) == 31
    assert BQ.is_single_block(BOX, STEP, 10)
    assert not BQ.is_single_block(BOX, BOX * 1.0, 10)  # scale factor 1.0: n_min_ergodic 2
    assert BQ.uniforms_per_candidate(10) == 33


def test_applicable_only_to_a_plain_mc_optimizer() -> None:
    class Logged(MCOptimizer):
        pass

    plain = MCOptimizer(None, None, rng=np.random.default_rng(0))  # type: ignore[arg-type]
    sub = Logged(None, None, rng=np.random.default_rng(0))  # type: ignore[arg-type]
    assert BQ.lockstep_applicable(plain, BOX, STEP, 10)
    assert not BQ.lockstep_applicable(sub, BOX, STEP, 10)
    assert not BQ.lockstep_applicable(plain, BOX, BOX, 10)


# ---------------------------------------------------------------------------------------------
# logic, with stub costs that force each outcome
# ---------------------------------------------------------------------------------------------


class StubCost:
    """quality = f(R): a smooth function of the angle to ``target`` (or a constant)."""

    def __init__(self, target: np.ndarray, kind: str) -> None:
        self.target, self.kind, self.eval_count = target, kind, 0

    def evaluate(self, R: np.ndarray, verts: Any = None, phase: int = 0) -> OverlapInfo:
        self.eval_count += 1
        if self.kind == "const":
            q = 0.5
        else:  # decreasing with the misorientation from the target (degrees)
            ang = np.degrees(Rotation.from_matrix(R @ self.target.T).magnitude())
            q = float(np.exp(-ang / 0.3))
        return OverlapInfo(quality=q, n_quality_points=1, peak_overlap=1, peak_on_detector=2)

    def many(self, R: np.ndarray) -> List[OverlapInfo]:
        return [self.evaluate(r) for r in R]


@pytest.mark.parametrize("n", [0, 1, 5, 40])
@pytest.mark.parametrize("kind", ["const", "smooth"])
def test_lockstep_equals_replayed_sequential_with_stub_cost(n: int, kind: str) -> None:
    rng0 = np.random.default_rng(7)
    target = Rotation.random(random_state=3).as_matrix()
    cands = [
        SearchCandidate(
            orientation=(
                target.copy()  # exactly at the optimum: no trial can improve (a failure)
                if i % 4 == 0
                else Rotation.from_rotvec(rng0.normal(size=3) * np.radians(0.4)).as_matrix()
                @ target
            )
        )
        for i in range(n)
    ]
    state = np.random.default_rng(11)
    buf = level_buffer(state, n)
    ref_cost, got_cost = StubCost(target, kind), StubCost(target, kind)
    ref, ref_runs, used = sequential_reference(
        copy.deepcopy(cands), ref_cost, None, 0, buf  # type: ignore[arg-type]
    )
    got_runs = run_lockstep(cands, got_cost.many, state)
    assert_same(cands, ref, got_runs, ref_runs)
    assert got_cost.eval_count == ref_cost.eval_count == 12 * n
    if kind == "const":  # every block fails: n_restarts 1, 33 draws used, the start is kept
        assert all(r["n_restarts"] == 1 and r["last_accept"] == -1 for r in got_runs)
        assert used == [33] * n
    elif n == 40:  # a mix of outcomes
        n_ok = sum(r["n_accept"] for r in got_runs)
        assert 0 < n_ok < n, n_ok
        assert sorted(set(used)) == [30, 33]
        assert [u == 33 for u in used] == [r["n_restarts"] == 1 for r in got_runs]


def test_convergence_stop_code_and_multi_block_refusal() -> None:
    target = np.eye(3)
    cands = [SearchCandidate(orientation=target.copy())]
    cost = StubCost(target, "smooth")  # quality 1 at the start: cost 0
    buf = level_buffer(np.random.default_rng(5), 1)
    ref, ref_runs, _ = sequential_reference(
        copy.deepcopy(cands), StubCost(target, "smooth"), None, 0, buf, max_conv=0.1  # type: ignore
    )
    runs = run_lockstep(cands, cost.many, np.random.default_rng(5), max_conv=0.1)
    assert_same(cands, ref, runs, ref_runs)
    assert runs[0]["stop"] == 2 and runs[0]["n_restarts"] == 1
    with pytest.raises(ValueError, match="single block"):
        BQ.run_batched_quick_mc(
            cands, cost.many, np.random.default_rng(0), BOX, BOX, 10, 5  # n_min_ergodic 2
        )


def test_generator_advances_by_33_n_whatever_the_outcome() -> None:
    for kind in ("const", "smooth"):
        for n in (1, 6, 25):
            target = np.eye(3)
            cands = [
                SearchCandidate(
                    orientation=Rotation.from_rotvec([0.001 * (i + 1), 0.0, 0.0]).as_matrix()
                )
                for i in range(n)
            ]
            rng = np.random.default_rng(99)
            ref_rng = copy.deepcopy(rng)
            run_lockstep(cands, StubCost(target, kind).many, rng)
            ref_rng.random(33 * n)
            assert rng.bit_generator.state == ref_rng.bit_generator.state, (kind, n)


# ---------------------------------------------------------------------------------------------
# the real cost
# ---------------------------------------------------------------------------------------------


def _local(c: Dict[str, Any], **kw: Any) -> VoxelCostFunction:
    return VoxelCostFunction(pixel_radius=0, mode="hard", **kw, **c["kwargs"])


@needs_c
@pytest.mark.parametrize("flags", [{}, {"cpp_overlap_rule": True}], ids=["default", "cpp_rule"])
@pytest.mark.parametrize("use_batched_cost", [False, True], ids=["evaluate", "batched_cost"])
def test_lockstep_equals_replayed_sequential_with_the_real_cost(
    three_voxels: Dict[str, Any], flags: Dict[str, Any], use_batched_cost: bool  # noqa: F811
) -> None:
    c = three_voxels
    phase, verts = c["voxel"].phase, c["vertices"]
    truth = np.asarray(c["voxel"].orientation, dtype=np.float64)
    # the truth and perturbations of 0.1, 0.4, 1 and 3 degrees: blocks that improve and that fail
    cands = [SearchCandidate(orientation=r.copy()) for r in orientations(truth, n_each=5)[:26]]
    n = len(cands)
    state = np.random.default_rng(2)
    buf = level_buffer(state, n)
    ref_cf, got_cf = _local(c, **flags), _local(c, **flags)
    ref, ref_runs, used = sequential_reference(copy.deepcopy(cands), ref_cf, verts, phase, buf)
    if use_batched_cost:
        ev = BC.BatchedCostEvaluator(got_cf, "cpu")
        many = lambda R: ev.evaluate_many(R, verts, phase)  # noqa: E731
    else:
        many = lambda R: [got_cf.evaluate(r, verts, phase) for r in R]  # noqa: E731
    got_runs = run_lockstep(cands, many, state)
    assert_same(cands, ref, got_runs, ref_runs)
    assert got_cf.eval_count == ref_cf.eval_count == 12 * n
    # the test sees both outcomes: blocks that improved and blocks that failed
    n_ok = sum(r["n_accept"] for r in got_runs)
    assert 0 < n_ok < n, f"{n_ok} of {n} blocks improved"
    print(f"\n{n_ok} of {n} blocks improved, used draws {sorted(set(used))}")
    assert [u == 33 for u in used] == [r["n_restarts"] == 1 for r in got_runs]


# ---------------------------------------------------------------------------------------------
# wiring: config keys and the reconstructor (S3 Step 3)
# ---------------------------------------------------------------------------------------------


def test_config_keys_default_off_and_parse(tmp_path: Any) -> None:
    from icenine.config_file import ConfigFile
    from icenine.orientation_search import SearchParameters

    base = Path(__file__).resolve().parents[1] / "scripts" / "phase_d" / "configs"
    text = (base / "ReconstructPhaseD_mc_clean.config").read_text()
    sp0 = SearchParameters.from_config(
        ConfigFile.from_file(str(base / "ReconstructPhaseD_mc_clean.config"))
    )
    assert sp0.batched_quick_mc is False and sp0.search_device == "cpu"
    f = tmp_path / "on.config"
    f.write_text(text + "\nBatchedQuickMC 1\nSearchDevice MPS\n")
    sp = SearchParameters.from_config(ConfigFile.from_file(str(f)))
    assert sp.batched_quick_mc is True and sp.search_device == "mps"
    f.write_text(text + "\nSearchDevice tpu\n")
    with pytest.raises(ValueError, match="SearchDevice"):
        ConfigFile.from_file(str(f))
    with pytest.raises(ValueError, match="search_device"):
        SearchParameters(search_device="tpu")


def _recon(**kw: Any) -> Tuple[Any, Any, np.ndarray]:
    from test_findoptimal_refactor import _build

    rec, voxel, R_true = _build()
    for k, v in kw.items():
        setattr(rec.setup.search_params, k, v)
    return rec, voxel, R_true


@needs_c
def test_key_off_builds_no_evaluator_and_key_on_runs_the_lockstep_quick_mc() -> None:
    from icenine.reconstructor import _get_voxel_vertices

    rec, voxel, R_true = _recon()
    events: List[Any] = []
    rec.recorder = lambda name, data: events.append((name, data))
    rec.reconstruct_voxel(_get_voxel_vertices(voxel), voxel.phase, rng=np.random.default_rng(7))
    assert rec._quick_mc_evaluators == {}  # the default path never touches the batched cost
    rec_on, voxel, R_true = _recon(batched_quick_mc=True)
    events_on: List[Any] = []
    rec_on.recorder = lambda name, data: events_on.append((name, data))
    res = rec_on.reconstruct_voxel(
        _get_voxel_vertices(voxel), voxel.phase, rng=np.random.default_rng(7)
    )
    (ev,) = rec_on._quick_mc_evaluators.values()
    n_cand = len(events_on[0][1]["R"])  # candidates of the one level
    assert n_cand == len(events[0][1]["R"])  # the discrete stage does not depend on the key
    assert ev.counts["orient_compact"] == 12 * n_cand and ev.counts["orient_fixed"] == 0
    assert ev.counts["calls_compact"] == 11  # starts together, then ten trial steps
    ang = np.degrees(Rotation.from_matrix(res.orientation @ R_true.T).magnitude())
    assert ang < 1.0 and np.isfinite(res.cost)  # a sane result (toy problem: one random draw)
    quick = dict(events_on[1][1])
    assert np.all(np.diff(quick["cost"]) >= 0)


@needs_c
def test_unsupported_levels_fall_back_to_the_sequential_quick_mc() -> None:
    from icenine.reconstructor import _get_voxel_vertices

    rec, voxel, _ = _recon(batched_quick_mc=True)
    verts, phase = _get_voxel_vertices(voxel), voxel.phase
    cands = [SearchCandidate(orientation=np.asarray(voxel.orientation, dtype=np.float64))]
    args = (cands, verts, phase, BOX, STEP)
    rec._quick_mc_evaluators = {}

    class Logged(MCOptimizer):
        pass

    cf = rec._make_local_cost_fn()
    with pytest.warns(RuntimeWarning, match="subclass of MCOptimizer"):
        assert not rec._batched_quick_mc(Logged(cf, verts, phase), *args)
    with pytest.warns(RuntimeWarning, match="subclass of MCOptimizer or more than one"):
        assert not rec._batched_quick_mc(
            MCOptimizer(cf, verts, phase), cands, verts, phase, BOX, BOX
        )
    # the ray rule is not covered by the batched cost (CppEdgeRules is, since S1 Step 4)
    ray_cf = rec._make_local_cost_fn(ray_rule=RayRule())
    with pytest.warns(RuntimeWarning, match="ray rule"):
        assert not rec._batched_quick_mc(MCOptimizer(ray_cf, verts, phase), *args)
    assert rec._quick_mc_evaluators == {id(ray_cf): None}  # remembered: one warning per cost fn
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        assert not rec._batched_quick_mc(MCOptimizer(ray_cf, verts, phase), *args)
    rec._quick_mc_evaluators = {}
    rec.setup.search_params.cpp_edge_rules = True  # CppEdgeRules now runs the lockstep quick MC
    edge_cf = rec._make_local_cost_fn()
    rec.setup.search_params.cpp_edge_rules = False
    assert edge_cf.cpp_edge_rules
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        edge_cands = [SearchCandidate(orientation=c.orientation.copy()) for c in cands]
        assert rec._batched_quick_mc(
            MCOptimizer(edge_cf, verts, phase, rng=np.random.default_rng(0)),
            edge_cands,
            *args[1:],
        )
    (ev_edge,) = rec._quick_mc_evaluators.values()
    assert ev_edge is not None and ev_edge.cf is edge_cf
    rec._quick_mc_evaluators = {}
    assert rec._batched_quick_mc(MCOptimizer(cf, verts, phase, rng=np.random.default_rng(0)), *args)
    assert cands[0].cost < 1.0
    rec.setup.search_params.batched_quick_mc = False
    assert not rec._batched_quick_mc(MCOptimizer(cf, verts, phase), *args)


@pytest.mark.skipif(
    not CF._HAS_C_RASTERIZE or not getattr(CF, "_C_HAS_RAY_RULE", False),
    reason="C extension without the ray rule",
)
def test_ray_rule_levels_run_sequentially_and_later_levels_batched() -> None:
    """Mixed mode under the ray rule: levels 0 .. rule.levels-1 score with a per-level ray cost
    function the batched cost does not cover (sequential quick MC, one RuntimeWarning each); the
    later levels use the strict local cost function, which the lockstep quick MC covers."""
    from icenine.ray_overlap import RayRule
    from icenine.reconstructor import _get_voxel_vertices

    rec, voxel, R_true = _recon(
        batched_quick_mc=True, max_local_resolution=2, ray_rule=RayRule(levels=2)
    )
    events: List[Any] = []
    rec.recorder = lambda name, data: events.append((name, data))
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        res = rec.reconstruct_voxel(
            _get_voxel_vertices(voxel), voxel.phase, rng=np.random.default_rng(7)
        )
    msgs = [str(w.message) for w in caught if issubclass(w.category, RuntimeWarning)]
    ray_msgs = [m for m in msgs if "BatchedQuickMC is not supported" in m]
    assert len(ray_msgs) == 2 and all("the ray rule" in m for m in ray_msgs)  # levels 0 and 1
    ev = [e for e in rec._quick_mc_evaluators.values() if e is not None]
    assert len(rec._quick_mc_evaluators) == 3 and len(ev) == 1  # 2 ray cost fns (None) + level 2
    n_cand = [len(d["R"]) for n, d in events if n == "quick_mc"]
    assert len(n_cand) == 3
    assert ev[0].counts["orient_compact"] == 12 * n_cand[2]  # only level 2 went through it
    assert np.isfinite(res.cost) and rec.last_find_optimal["n_candidates"] >= 1
