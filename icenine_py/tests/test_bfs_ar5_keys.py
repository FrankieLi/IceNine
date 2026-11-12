"""AR5 step 5: the two opt-in keys behind run S11 (both off by default).

Keys: ``SeedIndexCMAFirst n`` (the mixed seed refinement) and ``BFSNeighborCMANoFallback``
(the CMA end of a neighbour fit is returned, the MC fit does not follow). Two kinds of test, as in
``test_bfs_ar3_keys.py``:
- scripted groups (fake reconstructor, no images): defaults, parsing, bad values, every refused
  combination, key off equals key absent, the call arguments, the rule on synthetic cases with a
  fixed generator, the counters;
- the ThreeVoxels toy (real images, ``needs_c``): each new path equals the script code it was moved
  from (``mix_common.mixed_refine`` and ``fit_run.FitWrap.wrapped_no_fallback``) written out.
"""

from pathlib import Path
from types import SimpleNamespace
from typing import Any, Dict, List, Tuple

import numpy as np
import pytest

import icenine.reconstructor as reconstructor_module
from icenine import cost_functions as CF
from icenine import ray_indexing as RI
from icenine.orientation_search import SearchCandidate, SearchParameters
from icenine.reconstructor import BFSReconstruction, _get_voxel_vertices

# tests/ is on sys.path under pytest
from test_bfs_ar3_keys import CMAFirstFake, _lib_fit, toy, toy_params  # noqa: F401
from test_bfs_refit import FakeReconstructor, _angle, _info, _run, _rz
from test_seed_by_indexing import IndexFake, NoIndex, _config_with, _fingerprint, _reset

needs_c = pytest.mark.skipif(not CF._HAS_C_RASTERIZE, reason="C extension missing")

DEFAULTS: Dict[str, Any] = dict(seed_index_cma_first=0, bfs_neighbor_cma_no_fallback=False)
ON_LINES = (
    "SeedByIndexing 1\nBFSNeighborCMAFirst 1\nSeedIndexStopHitRatio 0.95\n"
    "SeedIndexCMAFirst 200\nBFSNeighborCMANoFallback 1\n"
)
# the keys a valid S11 parameter set needs
S11 = dict(
    seed_by_indexing=True,
    bfs_neighbor_cma_first=True,
    seed_index_stop_hit=0.95,
    seed_index_cma_first=200,
    bfs_neighbor_cma_no_fallback=True,
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
        _config_with(tmp_path, "SeedIndexCMAFirst 0\nBFSNeighborCMANoFallback 0")
    )
    for k, v in DEFAULTS.items():
        assert getattr(written, k) == v, k
    on = SearchParameters.from_config(_config_with(tmp_path, ON_LINES))
    assert on.seed_index_cma_first == 200 and on.bfs_neighbor_cma_no_fallback is True
    assert isinstance(on.seed_index_cma_first, int)


@pytest.mark.parametrize("n", [2, 3, 60, 200])
def test_cma_first_accepts_two_or_more(tmp_path: Path, n: int) -> None:
    on = SearchParameters.from_config(
        _config_with(
            tmp_path, f"SeedByIndexing 1\nSeedIndexStopHitRatio 0.9\nSeedIndexCMAFirst {n}"
        )
    )
    assert on.seed_index_cma_first == n


# ------------------------------------------------------------------------------ 2. bad values
@pytest.mark.parametrize(
    "line",
    [
        "SeedIndexCMAFirst 1",
        "SeedIndexCMAFirst -1",
        "SeedIndexCMAFirst",
        "SeedIndexCMAFirst x",
        "SeedIndexCMAFirst 2.5",
        "BFSNeighborCMANoFallback 2",
        "BFSNeighborCMANoFallback -1",
        "BFSNeighborCMANoFallback",
    ],
)
def test_bad_values_raise_at_parse_time(tmp_path: Path, line: str) -> None:
    with pytest.raises(ValueError):
        _config_with(tmp_path, line)


@pytest.mark.parametrize("n", [1, -1, -200])
def test_bad_values_raise_in_search_parameters(n: int) -> None:
    with pytest.raises(ValueError):
        SearchParameters(seed_by_indexing=True, seed_index_stop_hit=0.9, seed_index_cma_first=n)


def test_the_valid_set_builds() -> None:
    SearchParameters(**S11)


def _bfs_with(monkeypatch: pytest.MonkeyPatch, kw: Dict[str, Any]) -> None:
    monkeypatch.setattr(reconstructor_module, "AdaptiveVoxelReconstructor", FakeReconstructor)
    setup = SimpleNamespace(
        sample=SimpleNamespace(get_mic=lambda: None),
        config=SimpleNamespace(min_acceleration_threshold=0.8),
        search_params=SearchParameters(**kw),
    )
    BFSReconstruction(setup)  # type: ignore[arg-type]


def _without(*names: str, **change: Any) -> Dict[str, Any]:
    d = {k: v for k, v in S11.items() if k not in names}
    d.update(change)
    return d


@pytest.mark.parametrize(
    "kw",
    [
        _without("seed_by_indexing"),  # the mixed refinement acts on the indexed seed step
        _without("seed_index_stop_hit"),  # needs the stop ratio
        _without(paint_grains=True),
        _without(local_optimizer="cma"),
        _without(use_hybrid_optimizer=True),
        _without("bfs_neighbor_cma_first"),  # no fallback needs BFSNeighborCMAFirst
        _without("bfs_neighbor_cma_first", "seed_index_cma_first", seed_index_stop_hit=0.0),
        _without(outward_search=True),
        dict(bfs_neighbor_cma_no_fallback=True),
    ],
)
def test_invalid_combinations_raise_before_any_search(
    monkeypatch: pytest.MonkeyPatch, kw: Dict[str, Any]
) -> None:
    with pytest.raises(ValueError):
        _bfs_with(monkeypatch, kw)


# Rows that only the AR5 step 5 refusal block refuses: the parent branch's check passes each of
# them (no BFSNeighborCMAFirst, no SeedIndexCandidates, stop ratio 0 where named), so each test
# matches the new message and would fail if that block were deleted.
NEW_ONLY = "SeedIndexCMAFirst and BFSNeighborCMANoFallback"


@pytest.mark.parametrize(
    "kw, message",
    [
        (
            _without(
                "bfs_neighbor_cma_first", "bfs_neighbor_cma_no_fallback", local_optimizer="cma"
            ),
            f"{NEW_ONLY} need LocalOptimizer mc",
        ),
        (
            _without(
                "bfs_neighbor_cma_first", "bfs_neighbor_cma_no_fallback", use_hybrid_optimizer=True
            ),
            f"{NEW_ONLY} cannot be combined with use_hybrid_optimizer",
        ),
        (
            dict(bfs_neighbor_cma_no_fallback=True, paint_grains=True, seed_by_indexing=True),
            f"{NEW_ONLY} cannot be combined with PaintGrains 1",
        ),
        (
            _without("seed_by_indexing", "seed_index_stop_hit", "bfs_neighbor_cma_first"),
            "SeedIndexCMAFirst needs SeedByIndexing 1",
        ),
        (
            _without(
                "seed_index_stop_hit", "bfs_neighbor_cma_first", "bfs_neighbor_cma_no_fallback"
            ),
            r"SeedIndexCMAFirst needs SeedIndexStopHitRatio > 0",
        ),
        (
            dict(bfs_neighbor_cma_no_fallback=True),
            "BFSNeighborCMANoFallback 1 needs BFSNeighborCMAFirst 1",
        ),
    ],
)
def test_refusals_only_the_new_block_makes(
    monkeypatch: pytest.MonkeyPatch, kw: Dict[str, Any], message: str
) -> None:
    with pytest.raises(ValueError) as e:
        _bfs_with(monkeypatch, kw)
    assert str(e.value).startswith(message), str(e.value)


def test_cma_first_refusal_messages_name_the_missing_key(monkeypatch: pytest.MonkeyPatch) -> None:
    with pytest.raises(ValueError, match="SeedByIndexing"):
        _bfs_with(monkeypatch, _without("seed_by_indexing"))
    with pytest.raises(ValueError, match="SeedIndexCMAFirst needs SeedByIndexing 1"):
        _bfs_with(monkeypatch, _without("seed_by_indexing", "seed_index_stop_hit"))
    with pytest.raises(ValueError, match="SeedIndexStopHitRatio"):
        _bfs_with(monkeypatch, _without("seed_index_stop_hit"))
    with pytest.raises(ValueError, match="BFSNeighborCMAFirst"):
        _bfs_with(monkeypatch, _without("bfs_neighbor_cma_first"))


@pytest.mark.parametrize(
    "kw",
    [
        S11,
        _without("bfs_neighbor_cma_no_fallback"),  # the seed key alone
        _without("seed_index_cma_first"),  # the no-fallback key alone (S10)
        _without("bfs_neighbor_cma_first", "bfs_neighbor_cma_no_fallback"),  # seed key, no CMA
    ],
)
def test_valid_combinations_build(monkeypatch: pytest.MonkeyPatch, kw: Dict[str, Any]) -> None:
    _bfs_with(monkeypatch, kw)


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
    field)."""
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
        runs.append((_fingerprint(bfs, mic), [tuple(sorted(e)) for e in bfs.index_log]))
    assert runs[0] == runs[1]
    assert not any("mixed" in k or "not_fell" in k for k in runs[1][0]["stats"])


def test_key_off_cma_first_run_has_no_new_counter(monkeypatch: pytest.MonkeyPatch) -> None:
    CMAFirstFake.stuck, CMAFirstFake.log = (2,), []
    bfs, _ = _run((0.5, 0.4), monkeypatch, fake=CMAFirstFake, bfs_neighbor_cma_first=True)
    assert "cma_not_fell" not in bfs.stats and bfs.stats["cma_first_fallback_runs"] == 1


# ------------------------------------------------------------------------------ 4. call arguments
class RecordingIndex(IndexFake):
    calls: List[Dict[str, Any]] = []
    mixed: List[Dict[str, Any]] = []

    def index_seed(self, voxel_vertices: Any, phase_index: int, **kw: Any) -> RI.SeedResult:
        RecordingIndex.calls.append(dict(kw))
        out = super().index_seed(voxel_vertices, phase_index)
        out.mixed = [dict(m) for m in RecordingIndex.mixed]
        return out


def test_index_seed_gets_the_new_arguments_only_with_the_key(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    _reset()
    RecordingIndex.calls, RecordingIndex.mixed = [], []
    _run((0.3,), monkeypatch, fake=RecordingIndex, seed_by_indexing=True, seed_index_stop_hit=0.9)
    assert RecordingIndex.calls == [dict(stop_hit=0.9)]  # today's call, nothing added
    RecordingIndex.calls = []
    _run(
        (0.3,),
        monkeypatch,
        fake=RecordingIndex,
        seed_by_indexing=True,
        seed_index_stop_hit=0.9,
        seed_index_cma_first=200,
        cma_sigma0_deg=0.05,
    )
    assert RecordingIndex.calls == [dict(stop_hit=0.9, cma_first_evals=200, cma_sigma0_deg=0.05)]


def test_mixed_counters_and_log_exist_only_with_the_key(monkeypatch: pytest.MonkeyPatch) -> None:
    _reset()
    RecordingIndex.mixed = [dict(cma_hit=0.97, mc_ran=False), dict(cma_hit=0.4, mc_ran=True)]
    bfs, _ = _run(
        (0.3, 0.3),
        monkeypatch,
        fake=RecordingIndex,
        seed_by_indexing=True,
        seed_index_stop_hit=0.9,
        seed_index_cma_first=200,
    )
    n = len(bfs.index_log)
    assert n >= 1 and all(len(e["mixed"]) == 2 for e in bfs.index_log)
    assert bfs.stats["index_mixed_cands"] == 2 * n
    assert bfs.stats["index_mixed_cma_reached"] == n and bfs.stats["index_mixed_mc_ran"] == n
    _reset()
    bfs, _ = _run(
        (0.3,), monkeypatch, fake=RecordingIndex, seed_by_indexing=True, seed_index_stop_hit=0.9
    )
    assert all("mixed" not in e for e in bfs.index_log)
    assert not any(k.startswith("index_mixed") for k in bfs.stats)


# ------------------------------------------------------------------------------ 5. the rule
def _fns(
    cma: Tuple[float, float, int], mc: Tuple[float, float, int], log: List[str]
) -> Tuple[Any, Any, Any]:
    """CMA end (hit, cost, evals) and MC end (hit, cost, evals); the end orientation is a rotation
    marked by its angle: 1 deg for the CMA end, 2 deg for the MC end."""
    R_cma, R_mc = _rz(1.0).astype(np.float64), _rz(2.0).astype(np.float64)

    def cma_fn(R: np.ndarray, cost: float) -> Tuple[np.ndarray, float, int]:
        log.append("cma")
        return R_cma, cma[1], cma[2]

    def mc_fn(R: np.ndarray, cost: float) -> Tuple[np.ndarray, float, int]:
        log.append("mc")
        return R_mc, mc[1], mc[2]

    def hit_fn(R: np.ndarray) -> float:
        log.append("hit")
        return cma[0] if np.array_equal(R, R_cma) else mc[0]

    return cma_fn, mc_fn, hit_fn


def _rule(
    cma: Tuple[float, float, int], mc: Tuple[float, float, int], tau: float = 0.95
) -> Tuple[float, float, int, Dict[str, Any], List[str]]:
    log: List[str] = []
    cma_fn, mc_fn, hit_fn = _fns(cma, mc, log)
    R, c, ev, info = RI.mixed_refinement(
        _rz(0.0).astype(np.float64), 0.5, cma_fn, mc_fn, hit_fn, tau
    )
    return round(_angle(R, _rz(0.0))), c, ev, info, log


def test_rule_cma_reaches_the_ratio_so_mc_does_not_run() -> None:
    which, c, ev, info, log = _rule((0.96, 0.2, 201), (0.99, 0.1, 500))
    assert (which, c, ev) == (1, 0.2, 201)
    assert log == ["cma", "hit"]  # MC never called
    assert info["mc_ran"] is False and info["cma_hit"] == 0.96


def test_rule_ratio_exactly_at_the_threshold_counts_as_reached() -> None:
    which, _c, _ev, info, log = _rule((0.95, 0.2, 201), (0.99, 0.1, 500))
    assert which == 1 and info["mc_ran"] is False and "mc" not in log


def test_rule_below_the_ratio_mc_runs_and_the_lower_cost_wins() -> None:
    which, c, ev, info, log = _rule((0.5, 0.2, 201), (0.9, 0.1, 500))
    assert (which, c) == (2, 0.1)
    assert ev == 201 + 500  # both fits are counted
    assert log == ["cma", "hit", "mc"]
    assert info["mc_ran"] is True and info["cma_hit"] == 0.5


def test_rule_mc_wins_only_on_a_strictly_lower_cost() -> None:
    which, c, ev, _i, _l = _rule((0.5, 0.2, 201), (0.9, 0.2, 500))  # tie: the CMA end
    assert (which, c, ev) == (1, 0.2, 701)
    which, c, ev, _i, _l = _rule((0.5, 0.2, 201), (0.9, 0.3, 500))  # MC higher: the CMA end
    assert (which, c, ev) == (1, 0.2, 701)


def test_rule_returns_a_copy_free_pair_not_the_mc_end_on_a_higher_cost() -> None:
    which, c, _ev, info, _l = _rule((0.0, 0.4, 10), (1.0, 0.41, 20), tau=1.0)
    assert which == 1 and c == 0.4 and info["mc_ran"] is True


@pytest.mark.parametrize("hit", [0.0, 0.9499999, 0.95, 1.0])
def test_rule_threshold_boundary(hit: float) -> None:
    *_x, info, log = _rule((hit, 0.2, 201), (0.9, 0.1, 500))
    assert info["mc_ran"] == (hit < 0.95) and (("mc" in log) == (hit < 0.95))


# ------------------------------------------------------------------------------ 6. the seed loop
K_STRICT = {0: 0.4, 1: 0.1, 2: 0.3, 3: 0.2}  # order of refinement: 1, 3, 2, 0
K_HIT = {1: 0.5, 3: 0.97, 2: 0.99, 0: 0.0}
K_EVALS = {1: 701, 3: 201, 2: 201, 0: 201}
K_COST = {1: 0.2, 3: 0.15, 2: 0.05, 0: 0.5}


def _candidate_index(R: np.ndarray) -> int:
    return int(round(_angle(R, _rz(0.0))))


def _loop(
    monkeypatch: pytest.MonkeyPatch, with_refine_fn: bool = True, **kw: Any
) -> Tuple[RI.SeedResult, List[int]]:
    Rs = np.stack([_rz(float(j)).astype(np.float64) for j in range(4)])
    monkeypatch.setattr(RI, "index_voxel", lambda *a, **k: (SimpleNamespace(R=Rs), np.arange(40)))
    refined: List[int] = []

    def refine_fn(R: np.ndarray, cost: float) -> Tuple[np.ndarray, float, int, Dict[str, Any]]:
        j = _candidate_index(R)
        refined.append(j)
        return (
            R * 1.0,
            K_COST[j],
            K_EVALS[j],
            dict(cma_hit=K_HIT[j], mc_ran=K_HIT[j] < 0.95),
        )

    return (
        RI.seed_by_indexing(
            SimpleNamespace(),
            None,
            None,
            None,
            0,
            None,
            lambda R: K_STRICT[_candidate_index(R)],
            n_cost=4,
            n_refine=3,
            stop_hit=0.95,
            hit_fn=lambda R: K_HIT[_candidate_index(R)],
            refine_fn=refine_fn if with_refine_fn else None,
            **kw,
        ),
        refined,
    )


def test_the_loop_uses_refine_fn_and_stops_at_the_first_end_reaching_the_ratio(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    out, refined = _loop(monkeypatch)
    assert refined == [1, 3] and out.refined == [1, 3]
    assert out.stopped_at == 1 and out.taken == 3 and out.cost == K_COST[3]
    assert out.evals == 4 + K_EVALS[1] + K_EVALS[3]  # costing plus the refinements made
    assert out.finisher_evals == [K_EVALS[1], K_EVALS[3]]
    assert [m["mc_ran"] for m in out.mixed] == [True, False]
    assert [m["cma_hit"] for m in out.mixed] == [K_HIT[1], K_HIT[3]]


def test_the_loop_without_refine_fn_has_no_mixed_field(monkeypatch: pytest.MonkeyPatch) -> None:
    class Rec:
        last_eval_counts = (0, 7, 0)

        def refine_from_candidates(self, cands: Any, v: Any, p: Any, rng: Any = None) -> Any:
            return SearchCandidate(orientation=cands[0].orientation * 1.0, cost=0.3)

    Rs = np.stack([_rz(float(j)).astype(np.float64) for j in range(4)])
    monkeypatch.setattr(RI, "index_voxel", lambda *a, **k: (SimpleNamespace(R=Rs), np.arange(40)))
    out = RI.seed_by_indexing(Rec(), None, None, None, 0, None, lambda R: 0.1, n_cost=4, n_refine=2)
    assert out.mixed == [] and out.refined == [0, 1] and out.evals == 4 + 14


def test_a_refine_function_needs_the_stop_ratio(monkeypatch: pytest.MonkeyPatch) -> None:
    Rs = np.stack([_rz(float(j)).astype(np.float64) for j in range(4)])
    monkeypatch.setattr(RI, "index_voxel", lambda *a, **k: (SimpleNamespace(R=Rs), np.arange(40)))
    with pytest.raises(ValueError):
        RI.seed_by_indexing(
            SimpleNamespace(), None, None, None, 0, None, lambda R: 0.1,
            refine_fn=lambda R, c: (R, c, 1, {}),
        )  # fmt: skip


def test_index_seed_refuses_the_cma_budget_without_the_stop_ratio() -> None:
    rec = reconstructor_module.AdaptiveVoxelReconstructor.__new__(
        reconstructor_module.AdaptiveVoxelReconstructor
    )
    with pytest.raises(ValueError):
        rec.index_seed(None, 0, cma_first_evals=200)  # type: ignore[arg-type]
    with pytest.raises(ValueError):
        rec.index_seed(None, 0, stop_hit=0.9, cma_first_evals=1)  # type: ignore[arg-type]


# ------------------------------------------------------------------------------ 7. no MC fallback
def test_no_fallback_returns_the_cma_end_and_never_calls_mc(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Voxel 1 is fitted by the CMA fit; voxel 2 is stuck (its CMA end is not below its start
    cost): with the key on the CMA end is returned as it is and MC does not run."""
    CMAFirstFake.stuck, CMAFirstFake.log = (2,), []
    bfs, _mic = _run(
        (0.5, 0.4),
        monkeypatch,
        fake=CMAFirstFake,
        bfs_neighbor_cma_first=True,
        bfs_neighbor_cma_no_fallback=True,
    )
    assert all(c[1] == "cma" for c in CMAFirstFake.log)  # no MC call at all
    assert bfs.records[1].n_evals == 61 and bfs.records[2].n_evals == 61
    assert bfs.stats["cma_not_fell"] == 1
    assert "cma_first_fallback_runs" not in bfs.stats
    assert bfs.stats["n_evals_neighbor"] == 61 + 61


def test_no_fallback_off_still_falls_back(monkeypatch: pytest.MonkeyPatch) -> None:
    CMAFirstFake.stuck, CMAFirstFake.log = (2,), []
    bfs, _ = _run((0.5, 0.4), monkeypatch, fake=CMAFirstFake, bfs_neighbor_cma_first=True)
    assert [c[1] for c in CMAFirstFake.log] == ["cma", "cma", "mc"]
    assert bfs.records[2].n_evals == 61 + 150


def test_no_fallback_counter_is_zero_when_every_cma_fit_falls(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    CMAFirstFake.stuck, CMAFirstFake.log = (), []
    bfs, _ = _run(
        (0.5, 0.4),
        monkeypatch,
        fake=CMAFirstFake,
        bfs_neighbor_cma_first=True,
        bfs_neighbor_cma_no_fallback=True,
    )
    assert bfs.stats["cma_not_fell"] == 0 and CMAFirstFake.log


def test_no_fallback_returns_the_cma_candidate_object(monkeypatch: pytest.MonkeyPatch) -> None:
    """The result of a stuck fit is the CMA fit's own (cost 0.5, the start orientation)."""
    CMAFirstFake.stuck, CMAFirstFake.log = (1,), []
    bfs, _ = _run(
        (0.5,),
        monkeypatch,
        fake=CMAFirstFake,
        bfs_neighbor_cma_first=True,
        bfs_neighbor_cma_no_fallback=True,
    )
    assert bfs.stats["cma_not_fell"] == 1
    assert bfs.records[1].n_evals == 61


# ------------------------------------------------------------------------------ 8. toy (real)
def _rot_deg(v: List[float]) -> np.ndarray:
    from scipy.spatial.transform import Rotation

    return Rotation.from_rotvec(np.radians(v)).as_matrix()


def _hand_mixed(rec: Any, verts: Any, phase: int, R: np.ndarray, cost: float, tau: float) -> Any:
    """``mix_common.mixed_refine`` written out (arm X200, first step 0.05)."""
    r = rec.local_optimization(
        voxel_vertices=verts,
        phase_index=phase,
        initial_orientation=np.asarray(R, np.float64),
        rng=np.random.default_rng(0),
        cma_sigma0_deg=0.05,
        cma_max_evals=200,
        optimizer="cma",
    )
    rc, cc, ec = (
        np.asarray(r.orientation, np.float64),
        float(r.cost),
        int(rec.last_local_optimization_evals),
    )
    info = rec.evaluate_overlap(rc, verts, phase)
    hc = float(BFSReconstruction._ratios(info)[0])
    if hc >= tau:
        return rc, cc, ec, False
    m = rec.refine_from_candidates(
        [SearchCandidate(orientation=R.copy(), cost=float(cost))],
        verts,
        phase,
        rng=np.random.default_rng(0),
    )
    em = int(rec.last_eval_counts[1])
    if float(m.cost) < cc:
        return np.asarray(m.orientation, np.float64), float(m.cost), ec + em, True
    return rc, cc, ec + em, True


@needs_c
@pytest.mark.parametrize("tau", [0.95, 1.0])
def test_index_seed_with_cma_first_equals_the_script_rule_on_the_toy(toy: Any, tau: float) -> None:
    rec, voxel, _ = toy
    verts = _get_voxel_vertices(voxel)
    lib = rec.index_seed(
        verts, voxel.phase, 20, 5, stop_hit=tau, cma_first_evals=200, cma_sigma0_deg=0.05
    )
    sc, crystal, params, strict_cf = rec._index_context(voxel.phase)
    res, _att = RI.index_voxel(sc, verts, crystal, params, RI.SEED_TOP_K)
    R = res.R
    k = min(20, len(R))
    costs = np.array(
        [float(strict_cf.evaluate(R[j], verts, voxel.phase).cost) for j in range(k)], np.float64
    )
    order = [int(j) for j in np.argsort(costs, kind="stable")[:5]]
    ends, evals, ran, stop = [], 0, [], None
    for pos, j in enumerate(order):
        Re, ce, ev, mc_ran = _hand_mixed(rec, verts, voxel.phase, R[j], float(costs[j]), tau)
        ends.append((j, Re, ce, ev))
        ran.append(mc_ran)
        evals += ev
        info = rec.evaluate_overlap(Re, verts, voxel.phase)
        if BFSReconstruction._ratios(info)[0] >= tau:
            stop = pos
            break
    best = int(np.argmin([e[2] for e in ends])) if stop is None else stop
    assert lib.refined == [e[0] for e in ends] and lib.stopped_at == (-1 if stop is None else stop)
    assert lib.finisher_costs == [e[2] for e in ends]
    assert lib.finisher_evals == [e[3] for e in ends]
    assert lib.taken == ends[best][0] and lib.cost == ends[best][2]
    assert np.array_equal(lib.orientation, ends[best][1])
    assert lib.evals == k + evals
    assert [m["mc_ran"] for m in lib.mixed] == ran
    assert all(set(m) == {"cma_hit", "mc_ran"} for m in lib.mixed)


@needs_c
def test_index_seed_without_cma_first_has_no_mixed_entries_on_the_toy(toy: Any) -> None:
    rec, voxel, _ = toy
    verts = _get_voxel_vertices(voxel)
    a = rec.index_seed(verts, voxel.phase, 20, 5, stop_hit=0.95)
    assert a.mixed == []


def _hand_no_fallback(rec: Any, verts: Any, phase: int, R0: np.ndarray, seed: int) -> Any:
    """``fit_run.FitWrap.wrapped_no_fallback`` written out: the CMA fit (first step 0.05, 60
    evaluations), its end returned whatever its cost."""
    r1 = rec.local_optimization(
        verts, phase, R0, rng=np.random.default_rng(seed),
        cma_sigma0_deg=0.05, cma_max_evals=60, optimizer="cma",
    )  # fmt: skip
    return r1, int(rec.last_local_optimization_evals), bool(r1.cost >= r1.start_cost)


@needs_c
@pytest.mark.parametrize("which", ["near", "truth"])
def test_no_fallback_equals_the_script_wrapper_on_the_toy(
    toy: Any, toy_params: Any, which: str  # noqa: F811
) -> None:
    rec0, voxel, R_true = toy
    p = toy_params
    p.bfs_neighbor_cma_first, p.cma_sigma0_deg, p.cma_neighbor_max_evals = True, 0.05, 60
    p.bfs_neighbor_cma_no_fallback = True
    R0 = R_true.copy() if which == "truth" else _rot_deg([0.3, -0.2, 0.1]) @ R_true
    old_o = voxel.orientation
    voxel.orientation = R0
    try:
        verts = _get_voxel_vertices(voxel)
        bfs = BFSReconstruction(rec0.setup)
        res, n_lib, brec, _hit = _lib_fit(bfs, voxel, verts, 5)
        hand, n_hand, not_fell = _hand_no_fallback(rec0, verts, voxel.phase, R0, 5)
        assert np.array_equal(res.orientation, hand.orientation) and res.cost == hand.cost
        assert n_lib == n_hand == brec.n_evals and n_lib <= 61
        assert bfs.stats["cma_not_fell"] == int(not_fell)
        assert "cma_first_fallback_runs" not in bfs.stats
        if which == "truth":
            assert not_fell  # from the truth the CMA fit cannot lower the cost
    finally:
        voxel.orientation = old_o
        p.bfs_neighbor_cma_first = False
        p.bfs_neighbor_cma_no_fallback = False
        p.cma_sigma0_deg, p.cma_neighbor_max_evals = 0.2, 250


@needs_c
def test_the_stock_paths_are_not_edited(toy: Any) -> None:
    """``local_optimization`` and ``refine_from_candidates`` keep their signatures."""
    import inspect

    from icenine.reconstructor import AdaptiveVoxelReconstructor as A

    assert "cma_first_evals" not in inspect.signature(A.local_optimization).parameters
    assert "cma_first_evals" not in inspect.signature(A.refine_from_candidates).parameters
