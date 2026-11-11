"""AR3 step 10: the opt-in keys of the indexed seed step and the boundary pass (all off by default).

Keys: ``BFSNeighborCMAFirst``, ``SeedIndexCandidates n_cost n_refine``, ``SeedIndexStopHitRatio``,
``BoundaryPass``. Two kinds of test, as for ``SeedByIndexing``:
- scripted groups (fake reconstructor, no images): defaults, parsing, bad values, key off equals
  key absent, the call arguments, the early stop, the order of the pass in the BFS;
- the ThreeVoxels toy (real images, ``needs_c``): each new path equals the script code it was
  moved from (``reseed_s8.E1Pick.e1_step``, ``run_wrapped.Wrap.wrapped``, the boundary pass of
  ``boundary_common.Env`` and ``boundary_lib.apply_rule``) written out by hand.
"""

from pathlib import Path
from types import SimpleNamespace
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
import pytest

import icenine.mic_file as mic_module
import icenine.reconstructor as reconstructor_module
from icenine import boundary_pass as BP
from icenine import cost_functions as CF
from icenine import orientation_search as OS
from icenine import ray_indexing as RI
from icenine.mic_file import MicFile, ReconstructionState, Voxel
from icenine.orientation_search import CMAOptimizer, SearchCandidate, SearchParameters
from icenine.reconstructor import BFSReconstruction, BFSVoxelRecord, _get_voxel_vertices

# tests/ is on sys.path under pytest
import bfs_scripted
from test_bfs_refit import (
    GOOD,
    POSITIONS,
    FakeReconstructor,
    RejectingFull,
    _angle,
    _info,
    _run,
    _rz,
    _states,
)
from test_seed_by_indexing import IndexFake, NoIndex, _config_with, _fingerprint, _reset

FITTED = ReconstructionState.FITTED
needs_c = pytest.mark.skipif(not CF._HAS_C_RASTERIZE, reason="C extension missing")

DEFAULTS: Dict[str, Any] = dict(
    bfs_neighbor_cma_first=False,
    seed_index_n_cost=5,
    seed_index_n_refine=2,
    seed_index_stop_hit=0.0,
    boundary_pass=False,
)
ON_LINES = (
    "SeedByIndexing 1\nBFSNeighborCMAFirst 1\nSeedIndexCandidates 20 5\n"
    "SeedIndexStopHitRatio 0.95\nBoundaryPass 1\n"
)


# ------------------------------------------------------------------------------ 1. defaults
def test_defaults_and_parsing(tmp_path: Path) -> None:
    sp = SearchParameters()
    for k, v in DEFAULTS.items():
        assert getattr(sp, k) == v, k
    from_cfg = SearchParameters.from_config(_config_with(tmp_path, ""))
    for k, v in DEFAULTS.items():
        assert getattr(from_cfg, k) == v, k
    written = SearchParameters.from_config(
        _config_with(
            tmp_path,
            "BFSNeighborCMAFirst 0\nSeedIndexCandidates 5 2\nSeedIndexStopHitRatio 0\n"
            "BoundaryPass 0",
        )
    )
    for k, v in DEFAULTS.items():
        assert getattr(written, k) == v, k
    on = SearchParameters.from_config(_config_with(tmp_path, ON_LINES))
    assert on.bfs_neighbor_cma_first is True and on.boundary_pass is True
    assert (on.seed_index_n_cost, on.seed_index_n_refine) == (20, 5)
    assert on.seed_index_stop_hit == 0.95 and on.seed_by_indexing is True


# ------------------------------------------------------------------------------ 2. bad values
@pytest.mark.parametrize(
    "line",
    [
        "SeedIndexCandidates 30 5",
        "SeedIndexCandidates 5 6",
        "SeedIndexCandidates 5 0",
        "SeedIndexCandidates 5",
        "SeedIndexCandidates",
        "SeedIndexCandidates 20 5 7",
        "SeedIndexStopHitRatio 1.5",
        "SeedIndexStopHitRatio -0.1",
        "SeedIndexStopHitRatio",
        "BoundaryPass 2",
        "BoundaryPass",
        "BFSNeighborCMAFirst 2",
    ],
)
def test_bad_values_raise_at_parse_time(tmp_path: Path, line: str) -> None:
    with pytest.raises(ValueError):
        _config_with(tmp_path, line)


@pytest.mark.parametrize(
    "kw",
    [dict(seed_index_n_cost=30), dict(seed_index_n_refine=6), dict(seed_index_stop_hit=1.5)],
)
def test_bad_values_raise_in_search_parameters(kw: Dict[str, Any]) -> None:
    with pytest.raises(ValueError):
        SearchParameters(seed_by_indexing=True, **kw)


@pytest.mark.parametrize(
    "kw",
    [
        dict(seed_index_n_cost=20, seed_index_n_refine=5),  # seed keys without SeedByIndexing
        dict(seed_index_stop_hit=0.9),
        dict(seed_by_indexing=True, paint_grains=True, bfs_neighbor_cma_first=True),
        dict(seed_by_indexing=True, paint_grains=True, seed_index_stop_hit=0.9),
        dict(seed_by_indexing=True, paint_grains=True, seed_index_n_cost=20, seed_index_n_refine=5),
        dict(seed_by_indexing=True, paint_grains=True, boundary_pass=True),
        dict(bfs_neighbor_cma_first=True, local_optimizer="cma"),
        dict(bfs_neighbor_cma_first=True, outward_search=True),
        dict(bfs_neighbor_cma_first=True, use_hybrid_optimizer=True),
    ],
)
def test_invalid_combinations_raise_before_any_search(
    monkeypatch: pytest.MonkeyPatch, kw: Dict[str, Any]
) -> None:
    monkeypatch.setattr(reconstructor_module, "AdaptiveVoxelReconstructor", FakeReconstructor)
    setup = SimpleNamespace(
        sample=SimpleNamespace(get_mic=lambda: None),
        config=SimpleNamespace(min_acceleration_threshold=0.8),
        search_params=SearchParameters(**kw),
    )
    with pytest.raises(ValueError):
        BFSReconstruction(setup)  # type: ignore[arg-type]


def test_keys_written_at_their_defaults_never_raise(monkeypatch: pytest.MonkeyPatch) -> None:
    _reset()
    for extra in ({}, {"seed_by_indexing": True}, {"local_optimizer": "cma"}):
        _run((0.5,), monkeypatch, fake=IndexFake, **DEFAULTS, **extra)


# ------------------------------------------------------------------------------ 3. key off
@pytest.mark.parametrize("seed_idx", [False, True])
@pytest.mark.parametrize("reseed", [False, True])
def test_keys_absent_and_at_defaults_are_identical(
    monkeypatch: pytest.MonkeyPatch, seed_idx: bool, reseed: bool
) -> None:
    """Same orientations, evaluation counts, records, stats (no new key) and ``index_log`` (no new
    field); the ``index_seed`` call is the one of today (two arguments)."""
    extra: Dict[str, Any] = dict(seed_by_indexing=seed_idx, bfs_reseed=reseed)
    runs = []
    for params in ({}, DEFAULTS):
        _reset(bad=(1,) if reseed else ())
        bfs, mic = _run(
            (3.0, 3.5, 3.2),
            monkeypatch,
            fake=IndexFake if seed_idx else NoIndex,
            **extra,
            **params,
        )
        assert bfs.boundary_log == []
        runs.append((_fingerprint(bfs, mic), [tuple(sorted(e)) for e in bfs.index_log]))
    assert runs[0] == runs[1]
    assert not any(
        "cma_first" in k or "boundary" in k or "stopped" in k for k in runs[1][0]["stats"]
    )


# ------------------------------------------------------------------------------ 4. golden
def test_key_off_scripted_grid_still_matches_the_golden() -> None:
    """The scripted 14 x 14 grid (``bfs_grid_golden.npz``) with the new keys written at their
    defaults: the whole fit sequence is unchanged."""
    gold = np.load(Path(__file__).parent / "data" / "bfs_grid_golden.npz")
    out = bfs_scripted.run_grid(
        reconstructor_module, mic_module, SearchCandidate, SearchParameters, dict(DEFAULTS)
    )
    assert np.array_equal(out["done"], gold["mc_done"])
    assert np.array_equal(out["log_idx"], gold["mc_log_idx"])
    assert float(out["rng_after"]) == float(gold["mc_rng_after"])


# ------------------------------------------------------------------------------ 5. call arguments
class RecordingIndex(IndexFake):
    calls: List[Tuple[Any, ...]] = []
    stopped_at = -1

    def index_seed(
        self,
        voxel_vertices: Any,
        phase_index: int,
        n_cost: int = 5,
        n_refine: int = 2,
        stop_hit: float = 0.0,
    ) -> RI.SeedResult:
        RecordingIndex.calls.append((n_cost, n_refine, stop_hit))
        out = super().index_seed(voxel_vertices, phase_index)
        out.stopped_at = RecordingIndex.stopped_at
        return out


@pytest.mark.parametrize(
    "kw,expect",
    [
        (dict(seed_index_n_cost=20, seed_index_n_refine=5), (20, 5, 0.0)),
        (dict(seed_index_stop_hit=0.95), (5, 2, 0.95)),
        (
            dict(seed_index_n_cost=20, seed_index_n_refine=5, seed_index_stop_hit=0.95),
            (20, 5, 0.95),
        ),
    ],
)
def test_index_seed_gets_the_counts_and_the_threshold(
    monkeypatch: pytest.MonkeyPatch, kw: Dict[str, Any], expect: Tuple[Any, ...]
) -> None:
    _reset()
    RecordingIndex.calls = []
    bfs, _mic = _run((0.3,), monkeypatch, fake=RecordingIndex, seed_by_indexing=True, **kw)
    assert RecordingIndex.calls and set(RecordingIndex.calls) == {expect}
    has_stop = "seed_index_stop_hit" in kw
    assert all(("stopped_at" in e) == has_stop for e in bfs.index_log)
    assert ("index_stopped_early" in bfs.stats) == has_stop


def test_stopped_at_is_logged_and_counted_only_with_the_stop_key(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    _reset()
    RecordingIndex.calls, RecordingIndex.stopped_at = [], 1
    try:
        bfs, _ = _run(
            (0.3,), monkeypatch, fake=RecordingIndex, seed_by_indexing=True, seed_index_stop_hit=0.9
        )
    finally:
        RecordingIndex.stopped_at = -1
    assert [e["stopped_at"] for e in bfs.index_log] == [1]
    assert bfs.stats["index_stopped_early"] == 1
    _reset()
    bfs, _ = _run(
        (0.3,), monkeypatch, fake=RecordingIndex, seed_by_indexing=True, seed_index_stop_hit=0.9
    )
    assert bfs.stats["index_stopped_early"] == 0  # the key exists although no seed step stopped


# ------------------------------------------------------------------------------ 6. early stop
K_STRICT = {0: 0.4, 1: 0.1, 2: 0.3, 3: 0.2}  # strict cost of candidate j; order 1, 3, 2, 0
K_FINAL = {1: 0.20, 3: 0.15, 2: 0.05, 0: 0.5}  # final cost of the run on candidate j
K_EVALS = {1: 100, 3: 200, 2: 300, 0: 400}


class ScriptedRec:
    """``refine_from_candidates`` with scripted final costs (candidate j is the rotation j deg)."""

    def __init__(self) -> None:
        self.last_eval_counts = (0, 0, 0)
        self.refined: List[int] = []

    def refine_from_candidates(
        self, cands: List[SearchCandidate], vertices: Any, phase: int, rng: Any = None
    ) -> SearchCandidate:
        j = int(round(_angle(cands[0].orientation, _rz(0.0))))
        self.refined.append(j)
        self.last_eval_counts = (0, K_EVALS[j], 0)
        return SearchCandidate(orientation=cands[0].orientation * 1.0, cost=K_FINAL[j])


def _seed_step(
    monkeypatch: pytest.MonkeyPatch, hits: Dict[int, float], **kw: Any
) -> Tuple[RI.SeedResult, ScriptedRec, List[int]]:
    Rs = np.stack([_rz(float(j)).astype(np.float64) for j in range(4)])
    monkeypatch.setattr(RI, "index_voxel", lambda *a, **k: (SimpleNamespace(R=Rs), np.arange(40)))
    rec = ScriptedRec()
    asked: List[int] = []

    def hit_fn(R: np.ndarray) -> float:
        j = int(round(_angle(R, _rz(0.0))))
        asked.append(j)
        return hits[j]

    def strict(R: np.ndarray) -> float:
        return K_STRICT[int(round(_angle(R, _rz(0.0))))]

    out = RI.seed_by_indexing(
        rec, None, None, None, 0, None, strict, n_cost=4, n_refine=3, hit_fn=hit_fn, **kw
    )
    return out, rec, asked


def test_early_stop_keeps_the_first_end_at_or_above_the_threshold(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    out, rec, asked = _seed_step(monkeypatch, {1: 0.7, 3: 0.95, 2: 0.99, 0: 0.0}, stop_hit=0.95)
    assert rec.refined == [1, 3] and out.refined == [1, 3] and asked == [1, 3]
    assert out.stopped_at == 1 and out.taken == 3 and out.cost == K_FINAL[3]  # 0.05 would be lower
    assert out.evals == 4 + K_EVALS[1] + K_EVALS[3]  # the costing plus the refinements made
    assert out.finisher_costs == [K_FINAL[1], K_FINAL[3]]


def test_without_a_passing_end_all_are_refined_and_the_lowest_cost_is_kept(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    out, rec, asked = _seed_step(monkeypatch, {1: 0.7, 3: 0.9, 2: 0.94, 0: 0.0}, stop_hit=0.95)
    assert rec.refined == [1, 3, 2] and asked == [1, 3, 2]
    assert out.stopped_at == -1 and out.taken == 2 and out.cost == K_FINAL[2]
    assert out.evals == 4 + K_EVALS[1] + K_EVALS[3] + K_EVALS[2]


def test_threshold_zero_is_todays_loop(monkeypatch: pytest.MonkeyPatch) -> None:
    """No hit evaluation, all refined, lowest cost; every field as the old loop made it."""
    out, rec, asked = _seed_step(monkeypatch, {}, stop_hit=0.0)
    assert asked == [] and rec.refined == [1, 3, 2]
    assert out.stopped_at == -1 and out.taken == 2 and out.cost == K_FINAL[2]
    assert out.refined == [1, 3, 2] and out.strict_costs == [0.4, 0.1, 0.3, 0.2]
    assert out.finisher_costs == [K_FINAL[1], K_FINAL[3], K_FINAL[2]]
    assert out.finisher_evals == [K_EVALS[1], K_EVALS[3], K_EVALS[2]]
    assert out.evals == 4 + K_EVALS[1] + K_EVALS[3] + K_EVALS[2] and out.n_cand == 4


def test_a_threshold_without_a_hit_function_is_refused(monkeypatch: pytest.MonkeyPatch) -> None:
    Rs = np.stack([_rz(float(j)).astype(np.float64) for j in range(4)])
    monkeypatch.setattr(RI, "index_voxel", lambda *a, **k: (SimpleNamespace(R=Rs), np.arange(40)))
    with pytest.raises(ValueError):
        RI.seed_by_indexing(ScriptedRec(), None, None, None, 0, None, lambda R: 0.1, stop_hit=0.9)


# ------------------------------------------------------------------------------ toy fixture
@pytest.fixture(scope="module")
def toy() -> Any:
    from test_findoptimal_refactor import _build

    return _build()


@pytest.fixture()
def toy_params(toy: Any) -> Any:
    """The toy's search parameters; the new fields are put back after the test."""
    p = toy[0].setup.search_params
    old = {k: getattr(p, k) for k in (*DEFAULTS, "cma_sigma0_deg", "cma_neighbor_max_evals")}
    yield p
    for k, v in old.items():
        setattr(p, k, v)


# ------------------------------------------------------------------------------ 7. stop on toy
def _e1_by_hand(rec: Any, verts: Any, phase: int, n_cost: int, n_refine: int, tau: float) -> Any:
    """The loop of ``reseed_s8.E1Pick.e1_step`` written out."""
    sc, crystal, params, strict_cf = rec._index_context(phase)
    res, _att = RI.index_voxel(sc, verts, crystal, params, RI.SEED_TOP_K)
    R = res.R
    k = min(n_cost, len(R))
    costs = np.array(
        [float(strict_cf.evaluate(R[j], verts, phase).cost) for j in range(k)], dtype=np.float64
    )
    order = [int(j) for j in np.argsort(costs, kind="stable")[:n_refine]]
    fin, stop = [], None
    for pos, j in enumerate(order):
        r = rec.refine_from_candidates(
            [SearchCandidate(orientation=R[j].copy(), cost=float(costs[j]))],
            verts, phase, rng=np.random.default_rng(0),
        )  # fmt: skip
        fin.append(
            dict(j=j, cost=float(r.cost), R=np.asarray(r.orientation, np.float64),
                 evals=int(rec.last_eval_counts[1]))  # fmt: skip
        )
        info = rec.evaluate_overlap(fin[-1]["R"], verts, phase)
        if BFSReconstruction._ratios(info)[0] >= tau:
            stop = pos
            break
    best = int(np.argmin([f["cost"] for f in fin])) if stop is None else stop
    return dict(
        fin=fin,
        best=best,
        stop=-1 if stop is None else stop,
        evals=k + sum(f["evals"] for f in fin),
    )


@needs_c
@pytest.mark.parametrize("tau", [0.95, 1.0])
def test_index_seed_with_the_stop_equals_the_e1_loop_on_the_toy(toy: Any, tau: float) -> None:
    rec, voxel, _ = toy
    verts = _get_voxel_vertices(voxel)
    lib = rec.index_seed(verts, voxel.phase, 20, 5, stop_hit=tau)
    hand = _e1_by_hand(rec, verts, voxel.phase, 20, 5, tau)
    b = hand["fin"][hand["best"]]
    assert lib.stopped_at == hand["stop"]
    assert lib.refined == [f["j"] for f in hand["fin"]]
    assert lib.finisher_costs == [f["cost"] for f in hand["fin"]]
    assert lib.taken == b["j"] and lib.cost == b["cost"] and np.array_equal(lib.orientation, b["R"])
    assert lib.evals == hand["evals"]


# index_seed(verts, phase) on the toy, recorded from the parent code (commit 4ff3ce6, before the
# counts and the stop were arguments), the parent package on PYTHONPATH.
PARENT_SEED = dict(taken=1, cost=0.180327868852459, evals=549, refined=[1, 0])
PARENT_SEED_FINISHER_COSTS = [0.180327868852459, 0.4479166666666665]
PARENT_SEED_ORIENTATION = np.array(
    [
        [0.9071078420456102, -0.42083749065569104, -0.007153276031125422],
        [0.4185461838800449, 0.9037026031100777, -0.09022581167005998],
        [0.04443483834569397, 0.0788505649357736, 0.9958956439057767],
    ]
)


@needs_c
@pytest.mark.parametrize("explicit", [False, True])
def test_index_seed_without_the_stop_is_the_old_call_on_the_toy(toy: Any, explicit: bool) -> None:
    """The bare call and the call with the defaults written out both equal the parent's call
    (values recorded from commit 4ff3ce6; float fields to 1e-12 relative)."""
    rec, voxel, _ = toy
    verts = _get_voxel_vertices(voxel)
    if explicit:
        a = rec.index_seed(verts, voxel.phase, 5, 2, stop_hit=0.0)
    else:
        a = rec.index_seed(verts, voxel.phase)
    assert a.stopped_at == -1
    assert (a.taken, a.evals, a.refined) == (
        PARENT_SEED["taken"],
        PARENT_SEED["evals"],
        PARENT_SEED["refined"],
    )
    np.testing.assert_allclose(a.cost, PARENT_SEED["cost"], rtol=1e-12, atol=0.0)
    np.testing.assert_allclose(a.finisher_costs, PARENT_SEED_FINISHER_COSTS, rtol=1e-12, atol=0.0)
    np.testing.assert_allclose(a.orientation, PARENT_SEED_ORIENTATION, rtol=1e-12, atol=1e-15)


# ------------------------------------------------------------------------------ 8. CMA first
class CMAFirstFake(FakeReconstructor):
    """``local_optimization`` with the ``optimizer`` argument. The CMA fit lowers the cost unless
    the voxel is in ``stuck``; the MC fit always reaches the truth. Calls are logged."""

    stuck: Tuple[int, ...] = ()
    log: List[Tuple[Any, ...]] = []

    def local_optimization(  # type: ignore[override]
        self,
        voxel_vertices: Any,
        phase_index: int,
        initial_orientation: np.ndarray,
        rng: Any = None,
        cma_sigma0_deg: Optional[float] = None,
        cma_max_evals: Optional[int] = None,
        optimizer: Optional[str] = None,
    ) -> SearchCandidate:
        i = self._idx(voxel_vertices)
        CMAFirstFake.log.append((i, optimizer, cma_sigma0_deg, cma_max_evals))
        R0 = np.asarray(initial_orientation, np.float64)
        if optimizer == "cma":
            self.last_local_optimization_evals = 61
            ok = i not in self.stuck
            R = self.truth[i].astype(np.float64) if ok else R0
            res = SearchCandidate(R, 0.1 if ok else 0.5, _info(0.9 if ok else 0.5))
            res.start_cost = 0.5  # type: ignore[attr-defined]
            return res
        self.last_local_optimization_evals = 150
        return SearchCandidate(self.truth[i].astype(np.float64), 0.05, _info(GOOD))


def test_cma_first_scripted_bookkeeping(monkeypatch: pytest.MonkeyPatch) -> None:
    """Voxel 1 is fitted by the CMA fit alone; voxel 2 is stuck, so the MC fit runs from the same
    start and its lower cost is kept; evaluations add up; the counters exist."""
    CMAFirstFake.stuck, CMAFirstFake.log = (2,), []
    bfs, mic = _run((0.5, 0.4), monkeypatch, fake=CMAFirstFake, bfs_neighbor_cma_first=True)
    assert CMAFirstFake.log[:3] == [
        (1, "cma", 0.2, 250),
        (2, "cma", 0.2, 250),
        (2, "mc", None, None),
    ]
    assert bfs.records[1].n_evals == 61 and bfs.records[2].n_evals == 61 + 150
    assert bfs.stats["cma_first_fallback_runs"] == 1 and bfs.stats["cma_first_fallback_won"] == 1
    assert bfs.stats["n_evals_neighbor"] == 61 + 61 + 150
    assert _states(mic) == [FITTED, FITTED, FITTED]


def test_cma_first_off_calls_local_optimization_as_before(monkeypatch: pytest.MonkeyPatch) -> None:
    CMAFirstFake.stuck, CMAFirstFake.log = (), []
    bfs, _ = _run((0.5,), monkeypatch, fake=CMAFirstFake)
    assert CMAFirstFake.log == [(1, None, None, None)]
    assert "cma_first_fallback_runs" not in bfs.stats


def _hand_wrap(
    bfs: Any, voxel: Any, verts: Any, R0: np.ndarray, seed: int, monkeypatch: Any
) -> Any:
    """``run_wrapped.Wrap.wrapped`` written out: the CMA fit (first step 0.05, 60 evaluations), the
    stock fit from the same start if the end cost is not below the start cost, the lower cost kept
    (the CMA end on a tie). The start cost is the first strict evaluation of the CMA call."""
    rec = bfs.reconstructor
    p = rec.params
    rng = np.random.default_rng(seed)
    first: List[float] = []
    orig_eval = CF.VoxelCostFunction.evaluate
    armed = [True]

    def rec_eval(self: Any, *a: Any, **k: Any) -> Any:
        info = orig_eval(self, *a, **k)
        if armed[0]:
            first.append(float(info.cost))
            armed[0] = False
        return info

    monkeypatch.setattr(CF.VoxelCostFunction, "evaluate", rec_eval)
    old = p.local_optimizer
    p.local_optimizer = "cma"
    try:
        r1 = rec.local_optimization(
            verts, voxel.phase, R0, rng=rng, cma_sigma0_deg=0.05, cma_max_evals=60
        )
    finally:
        p.local_optimizer = old
        monkeypatch.setattr(CF.VoxelCostFunction, "evaluate", orig_eval)
    n1 = int(rec.last_local_optimization_evals)
    if r1.cost < first[0]:
        return r1, n1, False
    r2 = rec.local_optimization(verts, voxel.phase, R0, rng=rng)
    n2 = int(rec.last_local_optimization_evals)
    return (r2 if r2.cost < r1.cost else r1), n1 + n2, True


def _lib_fit(bfs: Any, voxel: Any, verts: Any, seed: int) -> Any:
    bfs.stats = bfs._new_stats()
    rec = BFSVoxelRecord(0)
    res, hit, conf = bfs._local_fit(voxel, verts, np.random.default_rng(seed), rec, "neighbor")
    return res, int(bfs.reconstructor.last_local_optimization_evals), rec, hit


@needs_c
@pytest.mark.parametrize("which", ["near", "truth"])
def test_cma_first_equals_the_wrapper_on_the_toy(
    toy: Any, toy_params: Any, monkeypatch: pytest.MonkeyPatch, which: str
) -> None:
    rec0, voxel, R_true = toy
    p = toy_params
    p.bfs_neighbor_cma_first, p.cma_sigma0_deg, p.cma_neighbor_max_evals = True, 0.05, 60
    try:
        verts = _get_voxel_vertices(voxel)
        R0 = R_true.copy()
        if which == "near":
            R0 = _rot_deg([0.3, -0.2, 0.1]) @ R_true
        old_o = voxel.orientation
        voxel.orientation = R0
        bfs = BFSReconstruction(rec0.setup)
        res, n_lib, brec, _hit = _lib_fit(bfs, voxel, verts, 5)
        calls_cma_lib = bfs.stats.get("cma_first_fallback_runs")
        # the hand-written wrapper on a fresh reconstructor (same setup, same generator)
        p.bfs_neighbor_cma_first, p.cma_sigma0_deg, p.cma_neighbor_max_evals = False, 0.2, 250
        bfs2 = BFSReconstruction(rec0.setup)
        hand, n_hand, fell_back = _hand_wrap(bfs2, voxel, verts, R0, 5, monkeypatch)
        assert fell_back == (calls_cma_lib == 1)
        assert np.array_equal(res.orientation, hand.orientation) and res.cost == hand.cost
        assert n_lib == n_hand and brec.n_evals == n_hand
        if which == "near":
            assert not fell_back  # the CMA fit lowers the cost: its end and its evaluations
            assert n_lib <= 61
        voxel.orientation = old_o
    finally:
        p.bfs_neighbor_cma_first = False


@needs_c
def test_cma_first_falls_back_to_mc_from_the_truth(toy: Any, toy_params: Any) -> None:
    """From the truth the CMA fit does not lower the cost: the MC fit runs as well, the evaluations
    of both are counted, the lower cost is kept."""
    rec0, voxel, R_true = toy
    p = toy_params
    p.bfs_neighbor_cma_first, p.cma_sigma0_deg, p.cma_neighbor_max_evals = True, 0.05, 60
    old_o = voxel.orientation
    voxel.orientation = R_true.copy()
    try:
        bfs = BFSReconstruction(rec0.setup)
        verts = _get_voxel_vertices(voxel)
        res, n, _rec, _hit = _lib_fit(bfs, voxel, verts, 5)
        assert bfs.stats["cma_first_fallback_runs"] == 1
        assert n > 61  # the CMA fit's evaluations plus the MC fit's
    finally:
        voxel.orientation = old_o


@needs_c
def test_refine_from_candidates_makes_no_cma_call_with_cma_first_on(
    toy: Any, toy_params: Any, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Plan test 8, last part: with ``BFSNeighborCMAFirst`` on, ``refine_from_candidates`` (the
    seed step's refinement with ``LocalOptimizer mc``) makes no CMA call. A guard against the key
    being read outside ``_local_fit``; the grid search itself is not run here."""
    rec0, voxel, R_true = toy
    toy_params.bfs_neighbor_cma_first = True
    n = [0]
    orig = CMAOptimizer.optimize

    def counted(self: Any, *a: Any, **k: Any) -> Any:
        n[0] += 1
        return orig(self, *a, **k)

    monkeypatch.setattr(CMAOptimizer, "optimize", counted)
    verts = _get_voxel_vertices(voxel)
    rec0.refine_from_candidates(
        [SearchCandidate(orientation=R_true.copy(), cost=0.5)],
        verts, voxel.phase, rng=np.random.default_rng(0),
    )  # fmt: skip
    assert n[0] == 0


def _rot_deg(v: List[float]) -> np.ndarray:
    from scipy.spatial.transform import Rotation

    return Rotation.from_rotvec(np.radians(v)).as_matrix()


# ------------------------------------------------------------------------------ 9. start_cost
class _Quad:
    def __init__(self) -> None:
        self.costs: List[float] = []

    def evaluate(self, R: np.ndarray, vertices: Any = None, phase: int = 0) -> Any:
        c = float(_angle(R, np.eye(3))) ** 2
        self.costs.append(c)
        return SimpleNamespace(cost=c)


def test_start_cost_is_the_first_evaluation_of_a_cma_run() -> None:
    f = _Quad()
    start = _rot_deg([0.3, -0.2, 0.25])
    res = OS.CMAOptimizer(f, None, sigma0_deg=0.2, max_evals=200).optimize(start, seed=0)
    assert res.start_cost == f.costs[0] and res.start_cost > res.cost
    assert res.n_evals == len(f.costs)
    # the early returns carry it too
    g = _Quad()
    r1 = OS.CMAOptimizer(g, None, max_evals=1).optimize(start, seed=0)
    assert r1.start_cost == g.costs[0] and r1.stop_reason == "max_evals" and r1.n_evals == 1
    h = _Quad()
    r2 = OS.CMAOptimizer(h, None, max_convergence_cost=1e9).optimize(start, seed=0)
    assert r2.start_cost == h.costs[0] and r2.stop_reason == "converged_cost"
    assert OS.CMAResult(np.eye(3), 0.1, None).start_cost is None


# CMAOptimizer(_Quad(), None, sigma0_deg=0.2, max_evals=100).optimize(start, seed=3), recorded from
# the parent code (commit 4ff3ce6, before ``start_cost`` existed), the parent package on PYTHONPATH.
PARENT_CMA_ORIENTATION = np.array(
    [
        [0.9999999852453088, -4.3672235355953386e-05, 0.00016613885314559318],
        [4.374747523552882e-05, 0.9999998964853186, -0.0004528967992064885],
        [-0.0001661190569321757, 0.00045290406067949974, 0.9999998836411789],
    ]
)
PARENT_CMA_COST = 0.0007702389413883503
PARENT_CMA_N_EVALS, PARENT_CMA_STOP = 100, "max_evals"


def test_cma_result_is_otherwise_unchanged() -> None:
    """The branch's CMA run equals the parent's (values recorded from commit 4ff3ce6): the same
    orientation and cost (to 1e-12 relative, for float differences between machines), the same
    evaluations and stop reason. Only ``start_cost`` is new."""
    start = _rot_deg([0.3, -0.2, 0.25])
    a = OS.CMAOptimizer(_Quad(), None, sigma0_deg=0.2, max_evals=100).optimize(start, seed=3)
    np.testing.assert_allclose(a.orientation, PARENT_CMA_ORIENTATION, rtol=1e-12, atol=1e-15)
    np.testing.assert_allclose(a.cost, PARENT_CMA_COST, rtol=1e-12, atol=0.0)
    assert (a.n_evals, a.stop_reason) == (PARENT_CMA_N_EVALS, PARENT_CMA_STOP)
    assert a.start_cost is not None


# ------------------------------------------------------------------------------ 10. candidates
IDENT = np.eye(3)[None]


def _grid_neighbours(n: int = 6) -> List[List[int]]:
    """8-neighbourhood of an n x n grid, index x * n + y, sorted."""
    out = []
    for x in range(n):
        for y in range(n):
            out.append(
                sorted(
                    a * n + b
                    for a in range(max(0, x - 1), min(n, x + 2))
                    for b in range(max(0, y - 1), min(n, y + 2))
                    if (a, b) != (x, y)
                )
            )
    return out


def _two_orientation_grid() -> np.ndarray:
    R = np.stack([_rz(0.0 if i // 6 < 3 else 10.0).astype(np.float64) for i in range(36)])
    R[19] = _rz(10.4)  # (3, 1)
    R[20] = _rz(10.3)  # (3, 2)
    return R


def test_boundary_candidates_on_a_two_orientation_grid() -> None:
    R = _two_orientation_grid()
    isb, cands = BP.boundary_candidates(R, _grid_neighbours(), IDENT)
    expect = {x * 6 + y for x in (2, 3) for y in range(6)}
    assert set(np.where(isb)[0].tolist()) == expect and sorted(cands) == sorted(expect)
    assert cands[12] == [18]  # (2, 0): neighbours 18, 19 are within 0.5 deg of each other
    assert cands[13] == [18]  # (2, 1): 18 (10.0), 19 (10.4), 20 (10.3): 18 kept
    assert cands[14] == [19]  # (2, 2): 19, 20, 21 merge; 19 is first although 21 is closer to B
    # an independent O(n^2) reference
    for v in expect:
        far = [j for j in _grid_neighbours()[v] if _angle(R[v], R[j]) > 0.5]
        keep: List[int] = []
        for j in far:
            if not any(_angle(R[j], R[k]) <= 0.5 for k in keep):
                keep.append(j)
        assert cands[v] == keep, v


def test_boundary_candidates_skip_voxels_that_were_never_fitted() -> None:
    R = _two_orientation_grid()
    active = np.ones(36, bool)
    active[18] = False  # (3, 0): neither a target nor a source
    isb, cands = BP.boundary_candidates(R, _grid_neighbours(), IDENT, active=active)
    assert not isb[18] and 18 not in cands
    assert all(18 not in c for c in cands.values())
    assert cands[12] == [19]


# ------------------------------------------------------------------------------ 11. keep_hit
def test_keep_hit_orders_by_hit_then_cost_then_position() -> None:
    old = (0.30, 0.80)
    assert BP.keep_hit(old, []) == 0
    assert BP.keep_hit(old, [(0.40, 0.85)]) == 1  # higher hit wins although the cost is higher
    assert BP.keep_hit(old, [(0.20, 0.79)]) == 0  # lower cost does not beat a higher hit
    assert BP.keep_hit(old, [(0.25, 0.80)]) == 1  # equal hit: lower cost
    assert BP.keep_hit(old, [(0.35, 0.80)]) == 0
    assert BP.keep_hit(old, [(0.30, 0.80)]) == 0  # a full tie: the old end wins
    assert BP.keep_hit(old, [(0.30, 0.90), (0.30, 0.90), (0.20, 0.90)]) == 3
    assert BP.keep_hit(old, [(0.25, 0.90), (0.25, 0.90)]) == 1  # a tie among tries: the first


# ------------------------------------------------------------------------------ 12. in the BFS
def test_the_pass_runs_once_after_the_reseeds_and_before_the_write(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    events: List[str] = []
    orig_reseed = BFSReconstruction._reseed_clusters
    orig_restart = BFSReconstruction._restart_pass

    def reseed(self: Any, *a: Any, **k: Any) -> Any:
        events.append("reseed")
        return orig_reseed(self, *a, **k)

    def restart(self: Any, *a: Any, **k: Any) -> Any:
        events.append("restart")
        return orig_restart(self, *a, **k)

    monkeypatch.setattr(BFSReconstruction, "_reseed_clusters", reseed)
    monkeypatch.setattr(BFSReconstruction, "_restart_pass", restart)
    monkeypatch.setattr(BP, "run_boundary_pass", lambda bfs, mic: events.append("pass"))
    monkeypatch.setattr(MicFile, "write", lambda self, path, *a, **k: events.append("write"))
    _reset(bad=(1,))
    bfs, _ = _run(
        (3.0, 3.5, 3.2), monkeypatch, fake=IndexFake, seed_by_indexing=True, bfs_reseed=True,
        bfs_restart_pass=True, boundary_pass=True,
    )  # fmt: skip
    assert events == ["reseed", "restart", "pass"]
    events.clear()
    _reset(bad=(1,))
    bfs.reconstruct_sample(
        output_mic=str(tmp_path / "o.mic"), rng=SimpleNamespace(shuffle=lambda x: None)
    )
    assert events == ["reseed", "restart", "pass", "write"]


def test_with_the_key_off_the_pass_is_never_called(monkeypatch: pytest.MonkeyPatch) -> None:
    def boom(bfs: Any, mic: Any) -> None:  # pragma: no cover - must not run
        raise AssertionError("the pass ran with BoundaryPass off")

    monkeypatch.setattr(BP, "run_boundary_pass", boom)
    _run((0.5, -0.5), monkeypatch)
    _run((0.5, -0.5), monkeypatch, boundary_pass=False)


class PassFake:
    """A reconstructor stand-in for ``run_boundary_pass``. Truth of voxel i: ``TRUTH_DEG[i]``
    degrees about z. A strict score is 0.95 - 0.05 per degree of error; a CMA fit returns the
    truth when the start is within 3 degrees of it, else the start."""

    TRUTH_DEG = [2.0, 2.0, 2.0, 2.0]
    fits: List[Tuple[int, float, Optional[str], Optional[float], Optional[int]]] = []

    def __init__(self, setup: Any) -> None:
        self.last_local_optimization_evals = 0

    @staticmethod
    def _idx(v: Any) -> int:
        return FakeReconstructor._idx(v)

    def _hit(self, R: np.ndarray, i: int) -> float:
        return max(0.0, 0.95 - 0.05 * _angle(np.asarray(R, np.float64), _rz(self.TRUTH_DEG[i])))

    def _make_local_cost_fn(self) -> Any:
        me = self

        class CF_:
            def evaluate(self, R: np.ndarray, vertices: Any, phase: int) -> Any:
                return _info(me._hit(R, me._idx(vertices)))

        return CF_()

    def local_optimization(
        self, voxel_vertices: Any, phase_index: int, initial_orientation: np.ndarray,
        rng: Any = None, cma_sigma0_deg: Optional[float] = None,
        cma_max_evals: Optional[int] = None, optimizer: Optional[str] = None,
    ) -> SearchCandidate:  # fmt: skip
        i = self._idx(voxel_vertices)
        R0 = np.asarray(initial_orientation, np.float64)
        PassFake.fits.append((i, _angle(R0, _rz(0.0)), optimizer, cma_sigma0_deg, cma_max_evals))
        rng.integers(0, 10)  # a CMA fit draws its seed from the generator
        truth = _rz(self.TRUTH_DEG[i]).astype(np.float64)
        R = truth if _angle(R0, truth) < 3.0 else R0
        self.last_local_optimization_evals = 61
        return SearchCandidate(R, 1.0 - self._hit(R, i), _info(self._hit(R, i)))


def _pass_bfs(
    monkeypatch: pytest.MonkeyPatch, deg: List[float], fitted: Optional[List[bool]] = None
) -> Tuple[BFSReconstruction, MicFile]:
    monkeypatch.setattr(reconstructor_module, "AdaptiveVoxelReconstructor", PassFake)
    monkeypatch.setattr(BP, "proper_rotations", lambda sym: IDENT)
    n = len(deg)
    voxels = [
        Voxel(
            position=np.array([x, y, 0.0]), orientation=_rz(d).astype(np.float64), side_length=1.0
        )
        for (x, y), d in zip(POSITIONS[:n], deg)
    ]
    mic = MicFile(voxels)
    for i, v in enumerate(mic.voxels):
        v.reconstruction_id = (
            FITTED if (fitted is None or fitted[i]) else (ReconstructionState.NOT_VISITED)
        )
    setup = SimpleNamespace(
        sample=SimpleNamespace(get_mic=lambda: mic),
        config=SimpleNamespace(min_acceleration_threshold=0.8),
        exp_setup=SimpleNamespace(get_sample_symmetry=lambda: None),
        search_params=SearchParameters(boundary_pass=True),
    )
    bfs = BFSReconstruction(setup)  # type: ignore[arg-type]
    bfs.stats = bfs._new_stats()
    PassFake.fits = []
    return bfs, mic


def test_pass_candidates_come_from_the_map_before_the_pass(monkeypatch: pytest.MonkeyPatch) -> None:
    """Voxel 0 sits at 0 degrees, truth 2; voxels 1 to 3 at 2 (right). Voxel 0 changes to 2. Voxels
    1 and 2 still see voxel 0's OLD orientation as a candidate (a live map would drop their tries).
    The tries are CMA fits with the pass's settings."""
    bfs, mic = _pass_bfs(monkeypatch, [0.0, 2.0, 2.0, 2.0])
    res = BP.run_boundary_pass(bfs, mic)
    assert res.boundary_voxels == [0, 1, 2]
    assert [(t["voxel"], t["cand"]) for t in res.tries] == [(0, 1), (1, 0), (2, 0)]
    assert PassFake.fits == [
        (0, pytest.approx(2.0, abs=1e-4), "cma", 0.05, 60),
        (1, pytest.approx(0.0, abs=1e-4), "cma", 0.05, 60),
        (2, pytest.approx(0.0, abs=1e-4), "cma", 0.05, 60),
    ]
    assert res.changed == [0] and _angle(mic.voxels[0].orientation, _rz(2.0)) < 1e-4
    for i in (1, 2, 3):  # a full tie: the old end wins
        assert np.array_equal(mic.voxels[i].orientation, _rz(2.0).astype(np.float64))
    s = bfs.stats
    assert s["boundary_voxels"] == 3 and s["boundary_tries"] == 3 and s["boundary_changed"] == 1
    assert s["n_evals_boundary_pass"] == 3 * 61 and s["n_evals_boundary_score"] == 3
    assert len(bfs.boundary_log) == 1 and bfs.boundary_log[0]["voxel"] == 0
    assert _angle(bfs.boundary_log[0]["old"], _rz(0.0)) < 1e-6
    assert mic.voxels[0].overlap_ratio == pytest.approx(0.95, abs=5e-3)


def test_pass_uses_its_own_generator_and_skips_unfitted_voxels(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """The run's generator is not touched (the pass has its own ``default_rng(0)``). Voxel 1 was
    never fitted: it is neither a target nor a source."""
    bfs, mic = _pass_bfs(monkeypatch, [0.0, 9.0, 2.0, 2.0], fitted=[True, False, True, True])
    res = BP.run_boundary_pass(bfs, mic)
    assert res.boundary_voxels == [0, 2]
    assert [(t["voxel"], t["cand"]) for t in res.tries] == [(0, 2), (2, 0)]
    assert np.array_equal(mic.voxels[1].orientation, _rz(9.0).astype(np.float64))


# ------------------------------------------------------------------------------ 13. toy pass
@needs_c
def test_the_pass_equals_the_script_code_on_the_toy(toy: Any, toy_params: Any) -> None:
    """Voxel 1 turned 5 degrees from voxel 0 (they sit in one place; voxel 2 has no neighbour):
    ``run_boundary_pass`` equals ``Env.boundary``, ``try_fit`` and ``apply_rule(..., "hit")``
    written out: tries, ends, evaluations, the new map."""
    rec0, _voxel, _R = toy
    setup = rec0.setup
    mic = setup.sample.get_mic()
    old_R = [np.array(v.orientation, dtype=np.float64) for v in mic.voxels]
    old_id = [v.reconstruction_id for v in mic.voxels]
    try:
        R = np.stack(old_R)
        R[1] = _rot_deg([0.0, 0.0, 5.0]) @ R[0]
        for v, r in zip(mic.voxels, R):
            v.orientation, v.reconstruction_id = r.copy(), FITTED
        bfs = BFSReconstruction(setup)
        bfs.stats = bfs._new_stats()
        res = BP.run_boundary_pass(bfs, mic)

        # ---- by hand
        sym = RI.proper_rotations(setup.exp_setup.get_sample_symmetry())
        nbr = [
            sorted(int(j) for j in mic.get_neighbors(v, radius=2.0 * mic.voxels[v].side_length))
            for v in range(len(R))
        ]
        cands: Dict[int, List[int]] = {}
        for v in range(len(R)):
            far = [j for j in nbr[v] if RI.misorientation_deg(R[v], R[j], sym) > 0.5]
            keep: List[int] = []
            for j in far:
                if not any(RI.misorientation_deg(R[j], R[k], sym) <= 0.5 for k in keep):
                    keep.append(j)
            if far:
                cands[v] = keep
        assert sorted(cands) == res.boundary_voxels == [0, 1]
        hrec = bfs.reconstructor
        cost_fn = hrec._make_local_cost_fn()
        rng = np.random.default_rng(0)
        want_tries, want_map = [], R.copy()
        for v in sorted(cands):
            verts = _get_voxel_vertices(mic.voxels[v])
            phase = mic.voxels[v].phase
            info = cost_fn.evaluate(R[v], verts, phase)
            oh, _ = BFSReconstruction._ratios(info)
            rows = []
            for c in cands[v]:
                p = hrec.params
                p.local_optimizer = "cma"
                try:
                    r = hrec.local_optimization(
                        voxel_vertices=verts, phase_index=phase, initial_orientation=R[c],
                        rng=rng, cma_sigma0_deg=0.05, cma_max_evals=60,
                    )  # fmt: skip
                finally:
                    p.local_optimizer = "mc"
                h, _cf = BFSReconstruction._ratios(r.overlap_info)
                rows.append(dict(voxel=v, cand=c, end=np.asarray(r.orientation, np.float64),
                                 cost=float(r.cost), hit=h,
                                 evals=int(hrec.last_local_optimization_evals)))  # fmt: skip
            want_tries += rows
            cost = np.array([float(info.cost)] + [t["cost"] for t in rows])
            hit = np.array([oh] + [t["hit"] for t in rows])
            k = int(np.lexsort((np.arange(len(cost)), cost, -hit))[0])
            if k > 0:
                want_map[v] = rows[k - 1]["end"]
        assert len(res.tries) == len(want_tries) >= 2
        for a, b in zip(res.tries, want_tries):
            assert (a["voxel"], a["cand"], a["evals"]) == (b["voxel"], b["cand"], b["evals"])
            assert np.array_equal(a["end"], b["end"])
            assert a["cost"] == b["cost"] and a["hit"] == b["hit"]
        got = np.stack([np.asarray(v.orientation, np.float64) for v in mic.voxels])
        assert np.array_equal(got, want_map)
        assert bfs.stats["n_evals_boundary_pass"] == sum(t["evals"] for t in want_tries)
    finally:
        for v, r, i in zip(mic.voxels, old_R, old_id):
            v.orientation, v.reconstruction_id = r, i
