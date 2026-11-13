"""AR5 step 9, phase B: the opt-in key ``SeedRescueCMA sigma evals`` (off by default).

After a failed seed test of an indexed seed step, the kept ends of the step (one per refined
candidate) are tried in order of cost: a CMA local optimization from each (first step ``sigma``
degrees, ``evals`` evaluations, a fresh ``default_rng(0)``), at most 5; the first end whose hit
ratio reaches ``SeedIndexStopHitRatio`` is taken; if none does, the grid search runs as today.
Rule ``G1`` of ``scripts/autoresearch/ar5/rescue_rules.py`` (frozen in
``benchmarks/autoresearch/ar5/rescue/frozen.json``), moved into the library.

Scripted tests (fake reconstructor, no images): defaults, parsing, bad values, refused
combinations, the rule on synthetic steps, the call arguments of one rescue, the key off (absent
or written off) equal, the BFS bookkeeping with the key on, and that the rule equals the script's
``rescue_rules.rescue_choose`` on the same synthetic steps.
"""

import sys
from pathlib import Path
from types import SimpleNamespace
from typing import Any, Dict, List, Tuple

import numpy as np
import pytest

import icenine.reconstructor as reconstructor_module
from icenine import ray_indexing as RI
from icenine.orientation_search import SearchCandidate, SearchParameters
from icenine.reconstructor import AdaptiveVoxelReconstructor, BFSReconstruction

# tests/ is on sys.path under pytest
from test_bfs_refit import FakeReconstructor, _info, _run, _rz
from test_seed_by_indexing import IndexFake, _config_with, _fingerprint, _reset

DEFAULTS: Dict[str, Any] = dict(seed_rescue_sigma_deg=0.0, seed_rescue_evals=0)
BASE: Dict[str, Any] = dict(seed_by_indexing=True, seed_index_stop_hit=0.95)
ON: Dict[str, Any] = dict(BASE, seed_rescue_sigma_deg=0.25, seed_rescue_evals=1000)


# ------------------------------------------------------------------------------ 1. defaults
def test_defaults_and_parsing(tmp_path: Path) -> None:
    sp = SearchParameters()
    assert (sp.seed_rescue_sigma_deg, sp.seed_rescue_evals) == (0.0, 0)
    for extra in ("", "SeedRescueCMA 0 0"):
        p = SearchParameters.from_config(_config_with(tmp_path, extra))
        assert (p.seed_rescue_sigma_deg, p.seed_rescue_evals) == (0.0, 0), extra
    on = SearchParameters.from_config(
        _config_with(
            tmp_path, "SeedByIndexing 1\nSeedIndexStopHitRatio 0.95\nSeedRescueCMA 0.25 1000"
        )
    )
    assert (on.seed_rescue_sigma_deg, on.seed_rescue_evals) == (0.25, 1000)
    assert isinstance(on.seed_rescue_sigma_deg, float) and isinstance(on.seed_rescue_evals, int)


@pytest.mark.parametrize(
    "line",
    [
        "SeedRescueCMA",
        "SeedRescueCMA 0.25",
        "SeedRescueCMA 0.25 1000 5",
        "SeedRescueCMA 0.25 1",
        "SeedRescueCMA 0.25 0",
        "SeedRescueCMA 0 1000",
        "SeedRescueCMA -0.25 1000",
        "SeedRescueCMA 0.25 2.5",
        "SeedRescueCMA x 1000",
        "SeedRescueCMA nan 1000",
    ],
)
def test_bad_values_raise_at_parse_time(tmp_path: Path, line: str) -> None:
    with pytest.raises(ValueError):
        _config_with(tmp_path, line)


@pytest.mark.parametrize(
    "sigma, evals", [(0.25, 1), (0.25, 0), (0.0, 1000), (-0.1, 1000), (float("inf"), 1000)]
)
def test_bad_values_raise_in_search_parameters(sigma: float, evals: int) -> None:
    with pytest.raises(ValueError):
        SearchParameters(**dict(BASE, seed_rescue_sigma_deg=sigma, seed_rescue_evals=evals))


# ------------------------------------------------------------------------------ 2. combinations
def _bfs_with(monkeypatch: pytest.MonkeyPatch, kw: Dict[str, Any]) -> None:
    monkeypatch.setattr(reconstructor_module, "AdaptiveVoxelReconstructor", FakeReconstructor)
    setup = SimpleNamespace(
        sample=SimpleNamespace(get_mic=lambda: None),
        config=SimpleNamespace(min_acceleration_threshold=0.8),
        search_params=SearchParameters(**kw),
    )
    BFSReconstruction(setup)  # type: ignore[arg-type]


def _on(*drop: str, **change: Any) -> Dict[str, Any]:
    d = {k: v for k, v in ON.items() if k not in drop}
    d.update(change)
    return d


@pytest.mark.parametrize(
    "kw, message",
    [
        (_on("seed_by_indexing", "seed_index_stop_hit"), "SeedRescueCMA needs SeedByIndexing 1"),
        (_on("seed_index_stop_hit"), "SeedRescueCMA needs SeedIndexStopHitRatio > 0"),
        (_on(local_optimizer="cma"), "SeedRescueCMA needs LocalOptimizer mc"),
        (_on(use_hybrid_optimizer=True), "SeedRescueCMA cannot be combined with"),
    ],
)
def test_refused_combinations_raise_before_any_search(
    monkeypatch: pytest.MonkeyPatch, kw: Dict[str, Any], message: str
) -> None:
    with pytest.raises(ValueError) as e:
        _bfs_with(monkeypatch, kw)
    assert str(e.value).startswith(message), str(e.value)


S11 = dict(
    bfs_neighbor_cma_first=True,
    seed_index_cma_first=200,
    bfs_neighbor_cma_no_fallback=True,
    seed_index_n_cost=20,
    seed_index_n_refine=5,
    boundary_pass=True,
)


@pytest.mark.parametrize("kw", [ON, dict(ON, **S11)])
def test_valid_combinations_build(monkeypatch: pytest.MonkeyPatch, kw: Dict[str, Any]) -> None:
    _bfs_with(monkeypatch, kw)


def test_key_off_never_needs_the_seed_keys(monkeypatch: pytest.MonkeyPatch) -> None:
    _bfs_with(monkeypatch, dict(DEFAULTS))


# ------------------------------------------------------------------------------ 3. the rule
def _tries(hits: Dict[int, float], evals: int = 1001) -> Tuple[Any, Any, List[int]]:
    """cma_fn(R) returns an end marked by the kept end's index (a rotation of 10 + k degrees);
    hit_fn reads the hit of that index from ``hits`` (0.5 when absent) and no overlap info."""
    calls: List[int] = []
    kept = [_rz(float(k)).astype(np.float64) for k in range(7)]

    def cma_fn(R: np.ndarray) -> Tuple[np.ndarray, float, int]:
        k = min(range(7), key=lambda j: float(np.abs(kept[j] - R).max()))
        calls.append(k)
        return _rz(10.0 + k).astype(np.float64), 0.01 * k, evals

    def hit_fn(R: np.ndarray) -> Tuple[float, Any]:
        k = min(range(7), key=lambda j: float(np.abs(_rz(10.0 + j) - R).max()))
        return hits.get(k, 0.5), None

    return cma_fn, hit_fn, calls


def _kept(n: int) -> List[np.ndarray]:
    return [_rz(float(k)).astype(np.float64) for k in range(n)]


def test_rule_tries_in_order_of_cost_ties_by_position() -> None:
    cma_fn, hit_fn, calls = _tries({})
    out = RI.seed_rescue(_kept(4), [0.3, 0.1, 0.2, 0.1], cma_fn, hit_fn, 0.95)
    assert calls == [1, 3, 2, 0] and not out["rescued"] and out["taken"] == -1
    assert out["evals"] == 4 * 1001 and [t["k"] for t in out["tries"]] == [1, 3, 2, 0]


def test_rule_tries_at_most_five() -> None:
    cma_fn, hit_fn, calls = _tries({})
    out = RI.seed_rescue(_kept(7), [0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3], cma_fn, hit_fn, 0.95)
    assert calls == [6, 5, 4, 3, 2] and len(out["tries"]) == RI.SEED_RESCUE_MAX_TRIES == 5


def test_rule_first_end_at_the_ratio_is_taken_and_the_rest_not_tried() -> None:
    cma_fn, hit_fn, calls = _tries({1: 0.96, 0: 0.99})
    out = RI.seed_rescue(_kept(3), [0.3, 0.1, 0.2], cma_fn, hit_fn, 0.95)
    assert out["rescued"] and out["taken"] == 0 and calls == [1] and out["evals"] == 1001
    assert np.array_equal(out["tries"][0]["R"], _rz(11.0).astype(np.float64))
    assert out["tries"][0]["cost"] == 0.01 and out["tries"][0]["hit"] == 0.96


def test_rule_ratio_exactly_at_the_threshold_counts() -> None:
    cma_fn, hit_fn, calls = _tries({2: 0.95})
    out = RI.seed_rescue(_kept(3), [0.3, 0.1, 0.2], cma_fn, hit_fn, 0.95)
    assert calls == [1, 2] and out["taken"] == 1 and out["evals"] == 2002


def test_rule_just_below_the_ratio_is_a_miss() -> None:
    cma_fn, hit_fn, _ = _tries({0: 0.949999, 1: 0.949999, 2: 0.949999})
    out = RI.seed_rescue(_kept(3), [0.3, 0.1, 0.2], cma_fn, hit_fn, 0.95)
    assert not out["rescued"] and len(out["tries"]) == 3 and out["evals"] == 3003


def test_rule_without_a_kept_end_costs_nothing() -> None:
    cma_fn, hit_fn, calls = _tries({})
    out = RI.seed_rescue([], [], cma_fn, hit_fn, 0.95)
    assert out == dict(rescued=False, taken=-1, tries=[], evals=0) and calls == []


@pytest.mark.parametrize(
    "hits, costs",
    [
        ({}, [0.3, 0.1, 0.2]),
        ({2: 0.97}, [0.3, 0.1, 0.2]),
        ({0: 0.95}, [0.1, 0.1, 0.1, 0.1, 0.1, 0.1]),
        ({5: 0.99}, [0.6, 0.5, 0.4, 0.3, 0.2, 0.1]),
        ({6: 0.99}, [0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3]),
    ],
)
def test_rule_equals_the_script_rule(hits: Dict[int, float], costs: List[float]) -> None:
    sys.path.insert(0, str(Path(__file__).parent.parent / "scripts" / "autoresearch" / "ar5"))
    import rescue_rules as RR

    cma_fn, hit_fn, calls = _tries(hits)
    lib = RI.seed_rescue(_kept(len(costs)), costs, cma_fn, hit_fn, 0.95)
    cma2, hit2, calls2 = _tries(hits)

    def one(k: int) -> Dict[str, Any]:
        R, c, ev = cma2(_kept(len(costs))[k])
        return dict(hit=hit2(R)[0], evals=ev)

    scr = RR.rescue_choose(costs, one, tau=0.95, max_tries=RR.MAX_TRIES)
    assert (lib["rescued"], lib["taken"], lib["evals"]) == (
        scr["rescued"],
        scr["taken"],
        scr["evals"],
    )
    assert [t["k"] for t in lib["tries"]] == [t["k"] for t in scr["tries"]] and calls == calls2
    assert (RR.TAU, RR.MAX_TRIES) == (0.95, RI.SEED_RESCUE_MAX_TRIES)


# ------------------------------------------------------------------------------ 4. one rescue
class _LO:
    """Stand-in for the reconstructor in ``AdaptiveVoxelReconstructor.rescue_seed``."""

    def __init__(self, hit: float) -> None:
        self.calls: List[Dict[str, Any]] = []
        self.hit = hit
        self.last_local_optimization_evals = 0
        self.overlaps = 0

    def local_optimization(self, **kw: Any) -> SearchCandidate:
        kw["draw"] = kw["rng"].random()
        self.calls.append(kw)
        self.last_local_optimization_evals = kw["cma_max_evals"] + 1
        return SearchCandidate(np.asarray(kw["initial_orientation"]) @ _rz(1.0), 0.07)

    def evaluate_overlap(self, R: np.ndarray, vertices: Any, phase: int) -> Any:
        self.overlaps += 1
        return _info(self.hit)


def test_one_rescue_calls_cma_with_the_key_values_and_a_fresh_generator() -> None:
    fake = _LO(hit=0.5)
    kept = _kept(3)
    out = AdaptiveVoxelReconstructor.rescue_seed(
        fake, "VERTS", 2, kept, [0.3, 0.1, 0.2], 0.95, 0.25, 1000  # type: ignore[arg-type]
    )
    assert not out["rescued"] and len(fake.calls) == 3 and fake.overlaps == 3
    first = np.random.default_rng(0).random()
    for c, k in zip(fake.calls, [1, 2, 0]):
        assert c["voxel_vertices"] == "VERTS" and c["phase_index"] == 2
        assert np.array_equal(c["initial_orientation"], kept[k])
        assert c["initial_orientation"].dtype == np.float64
        assert (c["cma_sigma0_deg"], c["cma_max_evals"], c["optimizer"]) == (0.25, 1000, "cma")
        assert c["draw"] == first  # each try has its own default_rng(0)
    assert out["evals"] == 3 * 1001 and all(t["evals"] == 1001 for t in out["tries"])
    assert all(t["cost"] == 0.07 for t in out["tries"])


def test_one_rescue_reports_the_taken_end_and_its_overlap() -> None:
    fake = _LO(hit=0.97)
    out = AdaptiveVoxelReconstructor.rescue_seed(
        fake, "V", 0, _kept(2), [0.2, 0.1], 0.95, 0.25, 1000  # type: ignore[arg-type]
    )
    assert out["rescued"] and out["taken"] == 0 and len(fake.calls) == 1
    t = out["tries"][0]
    assert np.allclose(t["R"], _kept(2)[1] @ _rz(1.0)) and t["info"] is not None
    assert BFSReconstruction._ratios(t["info"])[0] == pytest.approx(0.97)


# ------------------------------------------------------------------------------ 5. BFS, key off
class RescueFake(IndexFake):
    """``IndexFake`` with kept ends (two, costs 0.05 and 0.2) and the library's ``rescue_seed``.
    The CMA rescue from a kept end lands on the truth for voxels in ``rescuable`` (hit 0.97 by
    the fake overlap) and on the start otherwise (the seed test's 0.3 for failing voxels)."""

    rescuable: Tuple[int, ...] = ()
    lo_calls: List[Dict[str, Any]] = []
    rescue_seed = AdaptiveVoxelReconstructor.rescue_seed

    def index_seed(self, voxel_vertices: Any, phase_index: int, **kw: Any) -> RI.SeedResult:
        out = super().index_seed(voxel_vertices, phase_index)
        if out.orientation is not None:
            R = np.asarray(out.orientation, np.float64)
            out.finisher_R = [R @ _rz(2.0), R @ _rz(3.0)]
        return out

    def local_optimization(self, voxel_vertices: Any, phase_index: int, **kw: Any) -> Any:
        if kw.get("optimizer") != "cma":
            return super().local_optimization(voxel_vertices, phase_index, **kw)
        i = self._idx(voxel_vertices)
        RescueFake.lo_calls.append(dict(idx=i, **{k: v for k, v in kw.items() if k != "rng"}))
        self.last_local_optimization_evals = int(kw["cma_max_evals"]) + 1
        R0 = np.asarray(kw["initial_orientation"], np.float64)
        R = self.truth[i].astype(np.float64) if i in self.rescuable else R0
        return SearchCandidate(R, 0.04, None)

    def evaluate_overlap(self, R: np.ndarray, vertices: Any, phase: int) -> Any:
        i = self._idx(vertices)
        if i in self.rescuable and np.allclose(R, self.truth[i]) and i in self.failing:
            if not getattr(self, "_grid_last", False) and RescueFake.lo_calls:
                return _info(0.97)
        return super().evaluate_overlap(R, vertices, phase)


def _reset_rescue(failing: Tuple[int, ...], rescuable: Tuple[int, ...] = ()) -> None:
    _reset(failing=failing)
    RescueFake.rescuable, RescueFake.lo_calls = rescuable, []


@pytest.mark.parametrize("reseed", [False, True])
def test_key_absent_and_written_off_are_identical(
    monkeypatch: pytest.MonkeyPatch, reseed: bool
) -> None:
    runs = []
    for params in ({}, DEFAULTS):
        _reset(bad=(1,) if reseed else (), failing=(0,))
        RescueFake.rescuable, RescueFake.lo_calls = (0,), []
        bfs, mic = _run(
            (3.0, 3.5, 3.2), monkeypatch, fake=RescueFake, bfs_reseed=reseed, **BASE, **params
        )
        runs.append((_fingerprint(bfs, mic), [tuple(sorted(e)) for e in bfs.index_log]))
        assert RescueFake.lo_calls == []  # no rescue with the key off
    assert runs[0] == runs[1]
    assert not any("rescue" in k for k in runs[1][0]["stats"])
    assert all("rescue" not in e for e in runs[1][1])


# ------------------------------------------------------------------------------ 6. BFS, key on
def _seed_record(bfs: Any, i: int) -> Any:
    return bfs.records[i]


def test_key_on_a_rescued_seed_skips_the_grid_search(monkeypatch: pytest.MonkeyPatch) -> None:
    _reset_rescue(failing=(0,), rescuable=(0,))
    bfs, mic = _run((0.5,), monkeypatch, fake=RescueFake, **ON)
    s = bfs.stats
    import test_bfs_refit as TBR

    assert TBR.RejectingFull.full_calls == []  # no grid search ran
    assert s["index_fallback_seed_test"] == 1 and s["index_rescue_runs"] == 1
    assert s["index_rescued"] == 1 and s["index_rescue_tries"] == 1
    assert s["index_rescue_evals"] == 1001
    assert s["n_evals_seed"] == 1234 + 1001  # the index step and the rescue are the seed's
    assert _seed_record(bfs, 0).n_evals == 1234 + 1001 and _seed_record(bfs, 0).source == "seed"
    assert "index_fallback_grid_evals" not in s
    e = bfs.index_log[0]
    assert e["fallback"] == "seed_test" and e["rescue"]["rescued"] is True
    assert e["rescue"]["evals"] == 1001 and len(e["rescue"]["tries"]) == 1
    assert set(e["rescue"]["tries"][0]) == {"k", "hit", "cost", "evals"}
    assert e["rescue"]["tries"][0]["k"] == 0  # the cheaper kept end first (cost 0.05)
    (c,) = RescueFake.lo_calls
    assert (c["cma_sigma0_deg"], c["cma_max_evals"], c["optimizer"]) == (0.25, 1000, "cma")
    assert np.allclose(mic.voxels[0].orientation, FakeReconstructor.truth[0])
    assert mic.voxels[0].overlap_ratio == pytest.approx(0.97)


def test_key_on_a_failed_rescue_runs_the_grid_with_the_bfs_generator(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    _reset_rescue(failing=(0,), rescuable=())
    seen: List[Any] = []
    orig = RescueFake.reconstruct_voxel

    def grid(self: Any, voxel_vertices: Any, phase_index: int, rng: Any) -> Any:
        seen.append(rng)
        return orig(self, voxel_vertices, phase_index, rng)

    monkeypatch.setattr(RescueFake, "reconstruct_voxel", grid)
    bfs, _ = _run((0.5,), monkeypatch, fake=RescueFake, **ON)
    s = bfs.stats
    assert len(seen) == 1 and s["index_rescue_runs"] == 1 and "index_rescued" not in s
    assert s["index_rescue_tries"] == 2 and s["index_rescue_evals"] == 2002
    assert s["index_fallback_grid_evals"] == 1000  # the grid alone
    assert s["n_evals_seed"] == 1234 + 2002 + 1000
    assert _seed_record(bfs, 0).n_evals == 1234 + 2002 + 1000
    e = bfs.index_log[0]
    assert e["rescue"]["rescued"] is False and e["grid_evals"] == 1000
    assert [t["k"] for t in e["rescue"]["tries"]] == [0, 1]


def test_key_on_passing_seeds_never_rescue(monkeypatch: pytest.MonkeyPatch) -> None:
    _reset_rescue(failing=(), rescuable=(0, 1))
    bfs, _ = _run((0.5,), monkeypatch, fake=RescueFake, **ON)
    assert RescueFake.lo_calls == [] and not any("rescue" in k for k in bfs.stats)


def test_key_on_no_candidate_goes_to_the_grid_without_a_rescue(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    _reset(no_cand=(0,))
    RescueFake.rescuable, RescueFake.lo_calls = (0,), []
    bfs, _ = _run((0.5,), monkeypatch, fake=RescueFake, **ON)
    assert RescueFake.lo_calls == [] and bfs.stats["index_fallback_no_candidate"] == 1
    assert "index_rescue_runs" not in bfs.stats
