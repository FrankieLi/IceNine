"""AR5 step 28, phase B: the opt-in keys ``RayIndexMaxPairs N`` and ``RayIndexFallbackFull 1``.

``RayIndexMaxPairs N`` (default 600,000, today's setting) is the number of ray pairs the indexer
keeps (a random subsample above it) in a seed step. ``RayIndexFallbackFull 1`` (default 0, off) is
rule ``PXF`` of ``scripts/autoresearch/ar5/pxf_rules.py`` (frozen in
``benchmarks/autoresearch/ar5/pxf/``), moved into the library: a seed step that fails the seed test
(no candidate, or a hit ratio under ``min_acceleration_threshold``) after an indexing that kept a
subsample (``n_pairs >= N``) is run again, the WHOLE step, with 600,000 pairs; the second step is
the result, its kept ends go to the rescue, and both steps' evaluations count as the seed's.

Tests: defaults, parsing and the refused combinations; the keyword wiring of ``index_seed`` and
the pair counts on ``SeedResult``; the BFS with the keys absent, written at their defaults and
on (a scripted reconstructor); that the library path equals the script wrapper's step
(``pxf_rules.run_pxf_step``) on the same scripted steps; a hexagonal fixture with the real indexer.
"""

import dataclasses
import sys
from pathlib import Path
from types import SimpleNamespace
from typing import Any, Dict, List, Tuple

import numpy as np
import pytest

from icenine import ray_indexing as RI
from icenine.orientation_search import (
    SearchParameters,
    check_ray_index_fallback_full,
    check_ray_index_max_pairs,
)
from icenine.reconstructor import AdaptiveVoxelReconstructor

# tests/ is on sys.path under pytest
from test_bfs_refit import _angle, _info, _run, _rz
from test_bfs_seed_rescue import RescueFake, _reset_rescue
from test_seed_by_indexing import _config_with, _fingerprint
from test_seed_index_stall_hit import BASE, RESCUE, _bfs_with

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts" / "autoresearch" / "ar5"))
import pxf_rules as PR  # noqa: E402

FULL = 600_000
SMALL = dict(BASE, ray_index_max_pairs=50_000)
PXF = dict(SMALL, ray_index_fallback_full=1)


# ------------------------------------------------------------------------------ 1. the keys
def test_defaults_are_off_and_the_keys_parse(tmp_path: Path) -> None:
    p = SearchParameters()
    assert (p.ray_index_max_pairs, p.ray_index_fallback_full) == (FULL, 0)
    for extra in ("", f"RayIndexMaxPairs {FULL}", "RayIndexFallbackFull 0"):
        q = SearchParameters.from_config(_config_with(tmp_path, extra))
        assert (q.ray_index_max_pairs, q.ray_index_fallback_full) == (FULL, 0), extra
    on = SearchParameters.from_config(
        _config_with(tmp_path, "SeedByIndexing 1\nRayIndexMaxPairs 50000\nRayIndexFallbackFull 1")
    )
    assert (on.ray_index_max_pairs, on.ray_index_fallback_full) == (50_000, 1)
    assert isinstance(on.ray_index_max_pairs, int) and isinstance(on.ray_index_fallback_full, int)


@pytest.mark.parametrize(
    "line",
    ["RayIndexMaxPairs", "RayIndexMaxPairs 0", "RayIndexMaxPairs -5", "RayIndexMaxPairs 2.5",
     "RayIndexMaxPairs x", "RayIndexMaxPairs 5 6", "RayIndexFallbackFull",
     "RayIndexFallbackFull 2", "RayIndexFallbackFull -1", "RayIndexFallbackFull 1.0",
     "RayIndexFallbackFull x", "RayIndexFallbackFull 1 1"],
)  # fmt: skip
def test_bad_values_raise_at_parse_time(tmp_path: Path, line: str) -> None:
    with pytest.raises(ValueError):
        _config_with(tmp_path, line)


@pytest.mark.parametrize("v", [0, -1, 1.5, True, "5", None])
def test_bad_max_pairs_raise_in_search_parameters(v: Any) -> None:
    with pytest.raises(ValueError):
        check_ray_index_max_pairs(v)
    with pytest.raises(ValueError):
        SearchParameters(ray_index_max_pairs=v)


@pytest.mark.parametrize("v", [2, -1, 1.0, True, "1", None])
def test_bad_fallback_values_raise_in_search_parameters(v: Any) -> None:
    with pytest.raises(ValueError):
        check_ray_index_fallback_full(v)
    with pytest.raises(ValueError):
        SearchParameters(ray_index_fallback_full=v)


def test_refused_combinations_say_why(monkeypatch: pytest.MonkeyPatch) -> None:
    # a changed pair count acts on the indexed seed step only
    with pytest.raises(ValueError, match="RayIndexMaxPairs needs SeedByIndexing 1"):
        _bfs_with(monkeypatch, dict(ray_index_max_pairs=50_000))
    with pytest.raises(ValueError, match="RayIndexFallbackFull needs SeedByIndexing 1"):
        _bfs_with(monkeypatch, dict(ray_index_fallback_full=1))
    # the repeat is the same step when the first one already used 600,000 pairs or more
    with pytest.raises(ValueError, match="RayIndexFallbackFull needs RayIndexMaxPairs below"):
        _bfs_with(monkeypatch, dict(BASE, ray_index_fallback_full=1))
    with pytest.raises(ValueError, match="RayIndexFallbackFull needs RayIndexMaxPairs below"):
        _bfs_with(monkeypatch, dict(BASE, ray_index_fallback_full=1, ray_index_max_pairs=FULL))
    with pytest.raises(ValueError, match="RayIndexFallbackFull needs RayIndexMaxPairs below"):
        _bfs_with(monkeypatch, dict(BASE, ray_index_fallback_full=1, ray_index_max_pairs=FULL + 1))
    # valid: the defaults without SeedByIndexing, a small count alone, both on
    _bfs_with(monkeypatch, dict(ray_index_max_pairs=FULL, ray_index_fallback_full=0))
    _bfs_with(monkeypatch, SMALL)
    _bfs_with(monkeypatch, dict(BASE, ray_index_max_pairs=FULL + 1))  # more pairs: allowed
    _bfs_with(monkeypatch, PXF)


# ------------------------------------------------------------------------------ 2. wiring
class _Ctx:
    index_seed = AdaptiveVoxelReconstructor.index_seed

    def __init__(self) -> None:
        self.params = RI.IndexParams(q_max=8.0)
        self.last_seed_stall_hit: List[Dict[str, Any]] = []

    def _index_context(self, phase: int) -> Any:
        return ("SC", "CRYSTAL", self.params, SimpleNamespace(evaluate=lambda R, v, p: None))


def test_index_seed_replaces_max_pairs_only_when_asked(monkeypatch: pytest.MonkeyPatch) -> None:
    got: List[Any] = []

    def fake(rec: Any, sc: Any, crystal: Any, verts: Any, phase: int, params: Any, *a: Any,
             **kw: Any) -> Any:  # fmt: skip
        got.append(params)
        return "RESULT"

    monkeypatch.setattr(RI, "seed_by_indexing", fake)
    rec = _Ctx()
    assert rec.index_seed("V", 0) == "RESULT"
    assert rec.index_seed("V", 0, max_pairs=50_000) == "RESULT"
    assert rec.index_seed("V", 0, max_pairs=0) == "RESULT"  # 0: the stored setting
    assert got[0] is rec.params and got[2] is rec.params  # not copied, not changed
    assert got[1].max_pairs == 50_000 and got[1].q_max == 8.0
    assert rec.params.max_pairs == FULL  # the shared parameters are never mutated
    assert dataclasses.replace(got[1], max_pairs=FULL) == rec.params


def test_seed_result_carries_the_pair_counts(monkeypatch: pytest.MonkeyPatch) -> None:
    seen: List[int] = []

    def index_voxel(sc: Any, verts: Any, crystal: Any, p: Any, top_k: int = 20, keep: bool = False):
        seen.append(p.max_pairs)
        return SimpleNamespace(R=np.zeros((0, 3, 3)), n_pairs=p.max_pairs - 1, n_votes=11), []

    monkeypatch.setattr(RI, "index_voxel", index_voxel)
    out = RI.seed_by_indexing(None, None, None, "V", 0, RI.IndexParams(max_pairs=77), lambda R: 0.5)
    assert seen == [77] and out.orientation is None
    assert (out.n_pairs, out.n_votes) == (76, 11)
    # a stand-in result without the counts (older fakes) leaves zeros
    monkeypatch.setattr(
        RI, "index_voxel", lambda *a, **k: (SimpleNamespace(R=np.zeros((0, 3, 3))), [])
    )
    out = RI.seed_by_indexing(None, None, None, "V", 0, RI.IndexParams(), lambda R: 0.5)
    assert (out.n_pairs, out.n_votes) == (0, 0)


# ------------------------------------------------------------------------------ 3. the BFS
class PxfFake(RescueFake):
    """A scripted seed step. At ``max_pairs`` m the voxel has ``min(m, available)`` pairs. A step
    with ``m >= good_from`` finds the truth (evaluations 150); a smaller one ends 4 degrees off
    (evaluations 100, hit ratio 0.3 by the fake overlap) or, in ``no_cand_below``, has no
    candidate. Every step leaves its own two kept ends and one stall-check log."""

    calls: List[Dict[str, Any]] = []
    available = 700_000
    good_from = FULL
    no_cand_below = 0

    def index_seed(self, voxel_vertices: Any, phase_index: int, **kw: Any) -> RI.SeedResult:
        PxfFake.calls.append(dict(kw))
        mp = int(kw.get("max_pairs", FULL))
        n_pairs = min(mp, PxfFake.available)
        i = self._idx(voxel_vertices)
        k = len(PxfFake.calls)
        if mp < PxfFake.no_cand_below:
            return RI.SeedResult(n_att=3, evals=7, index_s=0.01, total_s=0.01, n_pairs=n_pairs)
        good = mp >= PxfFake.good_from
        R = self.truth[i].astype(np.float64)
        R = R if good else R @ _rz(4.0).astype(np.float64)
        if kw.get("mc_stall_hit"):
            self.last_seed_stall_hit = [
                dict(n=100 * k, evals=100 * k, checks=[], continued=0, stopped=False)
            ]
        return RI.SeedResult(
            orientation=R, cost=0.05, taken=k, n_cand=3, n_att=40, evals=150 if good else 100,
            index_s=0.02, total_s=0.5, strict_costs=[0.3, 0.1], refined=[1, 0],
            finisher_costs=[0.05, 0.2], finisher_evals=[60, 40], n_pairs=n_pairs,
            n_votes=n_pairs // 10, finisher_R=[R @ _rz(2.0 * k), R @ _rz(2.0 * k + 1)],
        )  # fmt: skip

    def evaluate_overlap(self, R: np.ndarray, vertices: Any, phase: int) -> Any:
        i = self._idx(vertices)
        return _info(0.97 if _angle(np.asarray(R), self.truth[i]) < 1.0 else 0.3)


def _bfs(monkeypatch: pytest.MonkeyPatch, params: Dict[str, Any], rescue: bool = False,
         **fake: Any) -> Tuple[Any, Any]:  # fmt: skip
    _reset_rescue(failing=(), rescuable=())
    PxfFake.calls, PxfFake.available = [], fake.get("available", 700_000)
    PxfFake.good_from, PxfFake.no_cand_below = fake.get("good_from", FULL), fake.get("no_cand", 0)
    extra = RESCUE if rescue else {}
    return _run((0.0, 0.0, 0.0), monkeypatch, fake=PxfFake, **params, **extra)


def _pairs() -> List[int]:
    """The ``max_pairs`` each indexer call was sent (600,000 when none was sent)."""
    return [int(c.get("max_pairs", FULL)) for c in PxfFake.calls]


def _seed0(bfs: Any) -> Dict[str, Any]:
    (e,) = [x for x in bfs.index_log if x["voxel"] == 0]
    return e


def test_keys_absent_and_written_at_their_defaults_are_identical(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    runs = []
    for params in (
        BASE,
        dict(BASE, ray_index_max_pairs=FULL),
        dict(BASE, ray_index_fallback_full=0),
    ):
        bfs, mic = _bfs(monkeypatch, params)
        assert all("max_pairs" not in kw for kw in PxfFake.calls) and len(PxfFake.calls) == 1
        runs.append((_fingerprint(bfs, mic), [tuple(sorted(e)) for e in bfs.index_log]))
    assert runs[0] == runs[1] == runs[2]
    assert not any(
        k.startswith(("index_pairs", "index_votes", "index_fb")) for k in runs[0][0]["stats"]
    )
    assert all("ray_index" not in e for e in bfs.index_log)


def test_a_small_count_alone_sends_the_count_and_never_repeats(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    bfs, _ = _bfs(monkeypatch, SMALL)  # 50,000 pairs miss the truth: the seed test fails
    assert _pairs() == [50_000]
    e = _seed0(bfs)
    assert e["fallback"] == "seed_test" and e["used"] is False and e["evals"] == 100
    assert e["ray_index"] == dict(
        max_pairs=50_000, n_pairs=50_000, n_votes=5_000, fallback_full=False
    )
    s = bfs.stats
    assert s["index_pairs"] == 50_000 and s["index_votes"] == 5_000
    assert "index_fb_full_steps" not in s and s["index_evals"] == 100  # one step, no repeat


def test_the_fallback_repeats_the_whole_step_with_the_full_sample(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    bfs, _ = _bfs(monkeypatch, PXF)
    assert _pairs() == [50_000, FULL]
    e = _seed0(bfs)
    assert e["used"] is True and "fallback" not in e and e["taken"] == 2  # the second step's
    assert e["evals"] == 100 + 150  # both steps count as the seed's
    assert bfs.stats["n_evals_seed"] == 250 and bfs.records[0].n_evals >= 250
    r = e["ray_index"]
    assert r["fallback_full"] is True and (r["n_pairs"], r["n_votes"]) == (FULL, FULL // 10)
    assert r["first"] == dict(evals=100, hit=pytest.approx(0.3), n_pairs=50_000, n_votes=5_000,
                              index_s=pytest.approx(0.02))  # fmt: skip
    s = bfs.stats
    assert s["index_fb_full_steps"] == 1 and s["index_fb_full_passed"] == 1
    assert s["index_fb_full_evals_first"] == 100
    assert s["index_pairs"] == 50_000 + FULL and s["index_votes"] == 5_000 + FULL // 10
    assert s["index_used"] == 1 and s["index_evals_used"] == 250  # one seed step, both costs
    assert s["index_runs"] == 1


def test_no_repeat_when_the_first_step_passes(monkeypatch: pytest.MonkeyPatch) -> None:
    bfs, _ = _bfs(monkeypatch, PXF, good_from=50_000)
    assert _pairs() == [50_000]
    assert _seed0(bfs)["used"] is True and "first" not in _seed0(bfs)["ray_index"]
    assert bfs.stats.get("index_fb_full_steps", 0) == 0


def test_no_repeat_when_the_first_indexing_kept_every_pair(monkeypatch: pytest.MonkeyPatch) -> None:
    bfs, _ = _bfs(monkeypatch, PXF, available=30_000)  # fails, but 30,000 < 50,000: no subsample
    assert _pairs() == [50_000]
    assert _seed0(bfs)["fallback"] == "seed_test" and bfs.stats.get("index_fb_full_steps", 0) == 0


def test_a_voxel_with_exactly_n_pairs_counts_as_a_subsample(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    # the frozen rule: n_pairs == N is a subsample (the real repeat is then the same step, and
    # its evaluations are counted; the scripted fake finds the truth at 600,000 either way)
    bfs, _ = _bfs(monkeypatch, PXF, available=50_000)
    assert _pairs() == [50_000, FULL] and bfs.stats["index_fb_full_steps"] == 1
    assert bfs.stats["index_fb_full_passed"] == 1 and bfs.stats["index_evals"] == 250


def test_no_candidate_after_a_subsample_also_repeats(monkeypatch: pytest.MonkeyPatch) -> None:
    bfs, _ = _bfs(monkeypatch, PXF, no_cand=FULL)  # 50,000 pairs: no candidate; 600,000: truth
    assert _pairs() == [50_000, FULL]
    e = _seed0(bfs)
    assert e["used"] is True and e["ray_index"]["first"]["hit"] is None
    assert e["evals"] == 7 + 150 and bfs.stats["index_fb_full_passed"] == 1
    assert "index_fallback_no_candidate" not in bfs.stats  # the step that counts had a candidate


def test_the_second_step_failing_too_goes_on_to_the_rescue_with_its_ends(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    bfs, _ = _bfs(monkeypatch, PXF, rescue=True, good_from=10 * FULL)  # nothing finds the truth
    assert len(PxfFake.calls) == 2 and bfs.stats["index_fb_full_passed"] == 0
    e = _seed0(bfs)
    assert e["fallback"] == "seed_test" and e["rescue"]["rescued"] is False
    starts = [c["initial_orientation"] for c in PxfFake.lo_calls if c["idx"] == 0]
    truth = bfs_truth(bfs)
    second = [truth @ _rz(4.0) @ _rz(2.0 * 2), truth @ _rz(4.0) @ _rz(2.0 * 2 + 1)]
    assert len(starts) >= 2
    for got, want in zip(starts, second):  # the rescue starts from the SECOND step's kept ends
        assert np.allclose(got, want, atol=1e-5)
    assert e["evals"] == 200  # both steps' evaluations (the rescue's are logged apart)


def bfs_truth(bfs: Any) -> np.ndarray:
    return np.asarray(PxfFake.truth[0], np.float64)


def test_the_stall_logs_of_both_steps_are_kept_in_step_order(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    on = dict(PXF, seed_index_mc_stall_hit=500)
    bfs, _ = _bfs(monkeypatch, on)
    assert [kw.get("mc_stall_hit") for kw in PxfFake.calls] == [500, 500]
    e = _seed0(bfs)
    assert [x["n"] for x in e["stall_hit"]] == [100, 200]  # the first step's, then the second's
    assert bfs.stats["index_stallhit_refinements"] == 2


# ------------------------------------------------------------------------------ 4. = the script
def test_the_library_path_equals_the_script_wrappers_step(monkeypatch: pytest.MonkeyPatch) -> None:
    """The same scripted steps through ``pxf_rules.run_pxf_step`` (the phase A wrapper) and through
    the library's ``_index_seed``: same indexer calls, same step handed on, same evaluations."""
    for kw in (dict(), dict(good_from=50_000), dict(available=30_000), dict(no_cand=FULL),
               dict(good_from=10 * FULL), dict(available=50_000)):  # fmt: skip
        bfs, _ = _bfs(monkeypatch, PXF, **kw)
        lib_calls = [c.get("max_pairs", FULL) for c in PxfFake.calls]
        lib = _seed0(bfs)
        # the script: one step() per call; the controller's max_pairs is the indexer's setting
        PxfFake.calls = []
        ctl = PR.Controller()
        fake = PxfFake(SimpleNamespace())
        verts = [(0.0, 0.0, 0.0)]

        def step() -> Any:
            out = fake.index_seed(verts, 0, max_pairs=ctl.max_pairs)
            ctl.log.append(
                dict(n_rays=0, n_pairs=out.n_pairs, n_votes=out.n_votes, s=0.0,
                     max_pairs=ctl.max_pairs, stage=ctl.stage)
            )  # fmt: skip
            return out

        out, row = PR.run_pxf_step(
            step,
            lambda o: 0.97 if _angle(np.asarray(o.orientation), PxfFake.truth[0]) < 1.0 else 0.3,
            ctl,
            0.8,
        )
        assert [c["max_pairs"] for c in PxfFake.calls] == lib_calls, kw
        assert out.evals == lib["evals"], kw
        assert bool(row["fallback"]) == lib["ray_index"]["fallback_full"], kw
        assert out.n_pairs == lib["ray_index"]["n_pairs"], kw
        if out.orientation is not None:
            assert lib["taken"] == out.taken, kw
            hit = 0.97 if _angle(np.asarray(out.orientation), PxfFake.truth[0]) < 1.0 else 0.3
            assert (hit >= 0.8) == lib["used"], kw


# ------------------------------------------------------------------------------ 5. hexagonal
class _Rec:
    """Just enough reconstructor for ``seed_by_indexing`` with stop_hit off."""

    last_eval_counts = (0, 1, 1)

    def refine_from_candidates(self, cands: Any, verts: Any, phase: int, rng: Any = None) -> Any:
        return SimpleNamespace(cost=0.1, orientation=np.asarray(cands[0].orientation))


def test_hexagonal_runs_through_the_library_step(monkeypatch: pytest.MonkeyPatch) -> None:
    from pymatgen.core.lattice import Lattice
    from pymatgen.core.structure import Structure
    from test_ar5_pxf import clutter, make_rays, noisy, rand_rot

    from icenine.symmetry import CrystalSymmetry

    st = Structure(
        Lattice.hexagonal(2.95, 4.68), ["Ti", "Ti"], [[1 / 3, 2 / 3, 0.25], [2 / 3, 1 / 3, 0.75]]
    )
    rl = st.lattice.reciprocal_lattice
    r = range(-6, 7)
    hkl = np.array([(h, k, ll) for h in r for k in r for ll in r if (h, k, ll) != (0, 0, 0)])
    vec = rl.get_cartesian_coords(hkl.astype(np.float64))
    vec = vec[np.linalg.norm(vec, axis=1) <= 5.0]
    cr = RI.build_crystal(vec, RI.proper_rotations(CrystalSymmetry(st)), q_max=5.0)
    g = np.concatenate([noisy(make_rays(rand_rot(41), cr, 40, 6), 0.3, 7), clutter(cr, 600, 8)])

    def index_voxel(sc: Any, verts: Any, crystal: Any, p: Any, top_k: int = 20, keep: bool = False):
        return RI.index_rays(g, crystal, p, top_k), np.arange(len(g))

    monkeypatch.setattr(RI, "index_voxel", index_voxel)
    outs = {}
    for mp in (50_000, FULL):
        p = dataclasses.replace(RI.IndexParams(q_max=5.0), max_pairs=mp)
        outs[mp] = RI.seed_by_indexing(_Rec(), None, cr, "V", 0, p, lambda R: 0.5, 2, 1)
    assert outs[50_000].n_pairs == 50_000 and outs[FULL].n_pairs > 50_000
    assert outs[FULL].n_votes > outs[50_000].n_votes > 0
    assert outs[50_000].orientation is not None and outs[FULL].orientation is not None
