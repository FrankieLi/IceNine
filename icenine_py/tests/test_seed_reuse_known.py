"""AR5 step 29, phase B: the opt-in key ``SeedReuseKnown 1`` (rule ``RU`` in the library).

Before a seed step (or reseed) at voxel v, the BFS tries the orientations of the 8 nearest voxels
accepted so far (state ``FITTED``): one strict cost each at v, then one CMA fit (first step 0.5
degree, 200 evaluations) from the lowest-cost one, then a second fit from the lowest-cost start
more than 0.5 degree from the first. A fit that ends at hit ratio 0.95 or more is the seed result
and nothing is indexed; otherwise the seed step runs unchanged and the try's evaluations are added
to it. The rule of ``scripts/autoresearch/ar5/ru_rules.py`` (frozen in
``benchmarks/autoresearch/ar5/ru/``), moved into the library. Needs ``SeedByIndexing 1`` and
``RayIndexFallbackFull 1`` (the only setting it was gated on).

Tests: the default and parsing; the refused combinations; key absent and written ``0`` give the
same BFS (map, stats, ``index_log``) with no ``reuse_*`` counter; the pure helpers equal
``ru_rules`` on scripted cases; a success skips the indexing and returns the fit's end; a failure
runs the step unchanged with the try added; the try runs once when the full-sample step repeats;
the fit's overlap equals a fresh ``evaluate_overlap`` (real images); a hexagonal fixture.
"""

import sys
from pathlib import Path
from types import SimpleNamespace
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
import pytest
from scipy.spatial import cKDTree

from icenine import ray_indexing as RI
from icenine.mic_file import ReconstructionState
from icenine.orientation_search import SearchCandidate, SearchParameters

# tests/ is on sys.path under pytest
from test_bfs_refit import _angle, _info, _run, _rz
from test_bfs_seed_rescue import _reset_rescue
from test_ray_index_max_pairs import PXF, PxfFake
from test_seed_by_indexing import _config_with
from test_seed_index_stall_hit import BASE, _bfs_with

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts" / "autoresearch" / "ar5"))
import ru_rules as RR  # noqa: E402

FULL = 600_000
ON = dict(PXF, seed_reuse_known=1)
FITTED = ReconstructionState.FITTED


# ------------------------------------------------------------------------------ 1. the key
def test_default_is_off_and_the_key_parses(tmp_path: Path) -> None:
    assert SearchParameters().seed_reuse_known == 0
    for extra in ("", "SeedReuseKnown 0"):
        assert SearchParameters.from_config(_config_with(tmp_path, extra)).seed_reuse_known == 0
    on = SearchParameters.from_config(
        _config_with(
            tmp_path,
            "SeedByIndexing 1\nRayIndexMaxPairs 50000\nRayIndexFallbackFull 1\nSeedReuseKnown 1",
        )
    )
    assert on.seed_reuse_known == 1 and isinstance(on.seed_reuse_known, int)


@pytest.mark.parametrize(
    "line",
    ["SeedReuseKnown", "SeedReuseKnown 2", "SeedReuseKnown -1", "SeedReuseKnown 1.0",
     "SeedReuseKnown x", "SeedReuseKnown 1 1"],
)  # fmt: skip
def test_bad_values_raise_at_parse_time(tmp_path: Path, line: str) -> None:
    with pytest.raises(ValueError, match="SeedReuseKnown"):
        _config_with(tmp_path, line)


@pytest.mark.parametrize("v", [2, -1, 1.0, "1", None, True, False])
def test_bad_values_raise_in_search_parameters(v: Any) -> None:
    with pytest.raises(ValueError, match="SeedReuseKnown"):
        SearchParameters(seed_reuse_known=v)


def test_refused_combinations_say_why(monkeypatch: pytest.MonkeyPatch) -> None:
    with pytest.raises(ValueError, match="SeedReuseKnown needs SeedByIndexing 1"):
        _bfs_with(monkeypatch, dict(seed_reuse_known=1))
    with pytest.raises(ValueError, match="SeedReuseKnown needs SeedByIndexing 1"):
        _bfs_with(
            monkeypatch,
            dict(ray_index_max_pairs=50_000, ray_index_fallback_full=1, seed_reuse_known=1),
        )
    # the fallback key (so RayIndexMaxPairs under 600,000) is the setting the rule was gated on
    with pytest.raises(ValueError, match="SeedReuseKnown needs RayIndexFallbackFull 1"):
        _bfs_with(monkeypatch, dict(BASE, seed_reuse_known=1))
    with pytest.raises(ValueError, match="SeedReuseKnown needs RayIndexFallbackFull 1"):
        _bfs_with(monkeypatch, dict(BASE, seed_reuse_known=1, ray_index_max_pairs=50_000))
    # the paint path calls index_seed directly and would never reach the block
    paint = dict(
        seed_by_indexing=True,
        paint_grains=True,
        ray_index_max_pairs=50_000,
        ray_index_fallback_full=1,
    )  # no other seed key, so no older refusal can stand in for this one
    with pytest.raises(ValueError, match="SeedReuseKnown cannot be combined with PaintGrains"):
        _bfs_with(monkeypatch, dict(paint, seed_reuse_known=1))
    with pytest.raises(ValueError, match="SeedReuseKnown cannot be combined with PaintGrains"):
        _bfs_with(monkeypatch, dict(seed_by_indexing=True, paint_grains=True, seed_reuse_known=1))
    # the gated arms ran the MC local optimizer
    with pytest.raises(ValueError, match="SeedReuseKnown cannot be combined with use_hybrid"):
        _bfs_with(monkeypatch, dict(ON, use_hybrid_optimizer=True))
    with pytest.raises(ValueError, match="SeedReuseKnown needs LocalOptimizer mc"):
        _bfs_with(monkeypatch, dict(ON, local_optimizer="cma"))
    _bfs_with(monkeypatch, ON)  # valid
    _bfs_with(monkeypatch, dict(ON, seed_reuse_known=0))
    _bfs_with(monkeypatch, dict(PXF, seed_reuse_known=0, local_optimizer="mc"))


# ------------------------------------------------------------------------------ 2. the helpers
def _grid(n: int = 7) -> np.ndarray:
    """A square grid: many equal distances, so the tie rule matters."""
    return np.array([(float(i), float(j)) for j in range(n) for i in range(n)])


def test_constants_are_the_frozen_ones() -> None:
    assert (RI.REUSE_K_NEAREST, RI.REUSE_SIGMA0_DEG, RI.REUSE_MAX_EVALS) == (
        RR.K_NEAREST,
        RR.SIGMA0_DEG,
        RR.MAX_EVALS,
    )
    assert (RI.REUSE_SEP_DEG, RI.REUSE_ACCEPT) == (RR.SEP_DEG, RR.ACCEPT)
    assert (RR.K_NEAREST, RR.SIGMA0_DEG, RR.MAX_EVALS, RR.SEP_DEG, RR.ACCEPT) == (
        8, 0.5, 200, 0.5, 0.95,
    )  # fmt: skip


@pytest.mark.parametrize("seed", [0, 1, 2])
def test_nearest_accepted_equals_the_script_helper(seed: int) -> None:
    pos = _grid()
    tree = cKDTree(pos)
    rng = np.random.default_rng(seed)
    accepted = rng.random(len(pos)) < 0.4
    for v in (0, 3, 24, 25, 48):
        lib = RI.nearest_accepted(tree, pos, v, lambda j: bool(accepted[j]))
        ref = RR.nearest_accepted(tree, pos, v, lambda j: bool(accepted[j]))
        assert lib == ref and len(lib) <= 8 and v not in lib and all(accepted[j] for j in lib)


def test_nearest_accepted_ties_by_index_rejected_left_out_and_fewer() -> None:
    pos = _grid()
    tree = cKDTree(pos)
    v = 24  # the middle of the 7 x 7 grid
    # all accepted: the 4 at distance 1 (indices 17, 23, 25, 31), then the 4 at sqrt 2
    # (16, 18, 30, 32): ties are broken by voxel index
    out = RI.nearest_accepted(tree, pos, v, lambda j: True)
    assert out == [17, 23, 25, 31, 16, 18, 30, 32]
    # a rejected voxel is left out (and the next one in order takes its place)
    out = RI.nearest_accepted(tree, pos, v, lambda j: j != 23)
    assert 23 not in out and out == [17, 25, 31, 16, 18, 30, 32, 10]
    # the voxel itself is never a start; fewer than 8 when fewer are accepted; none
    assert sorted(RI.nearest_accepted(tree, pos, v, lambda j: j in (v, 3, 40))) == [3, 40]
    assert RI.nearest_accepted(tree, pos, v, lambda j: False) == []
    assert RI.nearest_accepted(tree, pos, v, lambda j: j == v) == []


def _sym_cubic() -> np.ndarray:
    from icenine.symmetry import CrystalSymmetry
    from pymatgen.core.lattice import Lattice
    from pymatgen.core.structure import Structure

    fcc = [[0.0, 0.0, 0.0], [0.0, 0.5, 0.5], [0.5, 0.0, 0.5], [0.5, 0.5, 0.0]]
    st = Structure(Lattice.cubic(3.6), ["Cu"] * 4, fcc)
    return RI.proper_rotations(CrystalSymmetry(st))


def test_second_start_is_the_cheapest_over_half_a_degree_from_the_first() -> None:
    sym = _sym_cubic()

    def mis(a: np.ndarray, b: np.ndarray) -> float:
        return float(RI.misorientation_deg(a[None], b[None], sym)[0])

    # costs: starts 0 (lowest), 1 (0.2 deg from 0), 2 (0.4), 3 (0.6 deg), 4 (a symmetry copy of 0)
    R = [
        _rz(0.0).astype(np.float64), _rz(0.2).astype(np.float64), _rz(0.4).astype(np.float64),
        _rz(0.6).astype(np.float64), _rz(90.0).astype(np.float64),
    ]  # fmt: skip
    order = [0, 1, 2, 3, 4]
    assert RI.second_start(order, R, mis) == 3 == RR.second_start(order, R, mis)
    # the exact half degree is not "more than"; a start 90 degrees about z is a symmetry copy
    assert RI.second_start([0, 1, 2], R[:3], mis) is None
    assert RI.second_start([0, 1], [R[0], R[4]], mis) is None
    assert RR.second_start([0, 1], [R[0], R[4]], mis) is None
    # the order is the cost order, not the position
    assert RI.second_start([3, 2, 1, 0], R, mis) == RR.second_start([3, 2, 1, 0], R, mis) == 0
    assert RI.second_start([0], R[:1], mis) is None


# ------------------------------------------------------------------------------ 3. the BFS
class ReuseFake(PxfFake):
    """``PxfFake`` with the pieces of the try. A CMA fit from a start within ``capture`` degrees of
    the voxel's truth ends on the truth (hit ratio 0.97 by the fake overlap), else on its start
    (0.3). A fit costs ``cma_max_evals + 1`` evaluations (the final overlap evaluation included)."""

    capture = 5.0
    sym = np.eye(3)[None]
    fit_calls: List[Dict[str, Any]] = []
    cost_calls: List[Tuple[int, float]] = []
    scripted_costs: Optional[List[float]] = None  # the strict costs in call order, if scripted

    def _index_context(self, phase_index: int) -> Any:
        def evaluate(R: np.ndarray, verts: Any, phase: int) -> Any:
            i = self._idx(verts)
            c = 0.01 * _angle(np.asarray(R), self.truth[i])
            if ReuseFake.scripted_costs is not None:
                c = ReuseFake.scripted_costs[len(ReuseFake.cost_calls)]
            ReuseFake.cost_calls.append((i, c))
            return SimpleNamespace(cost=c)

        return None, SimpleNamespace(sym=ReuseFake.sym), None, SimpleNamespace(evaluate=evaluate)

    def local_optimization(self, voxel_vertices: Any, phase_index: int, **kw: Any) -> Any:
        if kw.get("optimizer") != "cma":
            return super().local_optimization(voxel_vertices, phase_index, **kw)
        i = self._idx(voxel_vertices)
        R0 = np.asarray(kw["initial_orientation"], np.float64)
        ReuseFake.fit_calls.append(
            dict(idx=i, start=R0, sigma0=kw["cma_sigma0_deg"], max_evals=kw["cma_max_evals"],
                 draw=float(kw["rng"].random()))
        )  # fmt: skip
        close = _angle(R0, self.truth[i]) < ReuseFake.capture
        self.last_local_optimization_evals = int(kw["cma_max_evals"]) + 1
        R = self.truth[i].astype(np.float64) if close else R0
        return SearchCandidate(R, 0.04 if close else 0.5, _info(0.97 if close else 0.3))


def _reuse(monkeypatch: pytest.MonkeyPatch, params: Dict[str, Any], capture: float = 5.0,
           **fake: Any) -> Tuple[Any, Any]:  # fmt: skip
    ReuseFake.capture, ReuseFake.fit_calls, ReuseFake.cost_calls = capture, [], []
    ReuseFake.sym = fake.pop("sym", np.eye(3)[None])
    _reset_rescue(failing=(), rescuable=())
    PxfFake.calls, PxfFake.available = [], fake.get("available", 700_000)
    PxfFake.good_from, PxfFake.no_cand_below = fake.get("good_from", FULL), fake.get("no_cand", 0)
    return _run((3.0, 3.5, 3.2), monkeypatch, fake=ReuseFake, **params)


def _entry(bfs: Any, voxel: int) -> Dict[str, Any]:
    (e,) = [x for x in bfs.index_log if x["voxel"] == voxel]
    return e


def _fp(bfs: Any, mic: Any) -> Dict[str, Any]:
    from dataclasses import asdict

    stats = {k: v for k, v in bfs.stats.items() if "wall" not in k and k != "reseeds"}
    records = {}
    for k, r in bfs.records.items():
        d = asdict(r)
        d.pop("wall_s")
        records[k] = d
    log = [{k: v for k, v in e.items() if k not in ("index_s", "path_s")} for e in bfs.index_log]
    return dict(
        stats=stats,
        records=records,
        log=log,
        orientations=[np.asarray(v.orientation).tolist() for v in mic.voxels],
        states=[v.reconstruction_id for v in mic.voxels],
    )


def test_absent_and_written_zero_are_identical_and_leave_no_trace(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    runs = []
    for params in (PXF, dict(PXF, seed_reuse_known=0)):
        bfs, mic = _reuse(monkeypatch, params)
        assert ReuseFake.fit_calls == [] and ReuseFake.cost_calls == []  # no try, no cost
        assert not hasattr(bfs, "_reuse_geo") or bfs._reuse_geo is None  # no tree built
        assert not any(k.startswith("reuse_") for k in bfs.stats)
        assert all("reuse" not in e for e in bfs.index_log)
        runs.append(_fp(bfs, mic))
    assert runs[0] == runs[1]
    # this compares key absent with key written 0 on the new code; identity with the build before
    # the change rests on code inspection plus the full suite and its goldens, not on this test
    assert [e["voxel"] for e in runs[0]["log"]] == [0, 3]


def test_a_success_skips_the_indexing_and_returns_the_fit_end(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    bfs, mic = _reuse(monkeypatch, ON)
    # voxel 0: no accepted voxel, so no try, and the PXF step runs (50,000 pairs fail, then the
    # full sample); voxel 3: the try from voxel 0's orientation reaches the truth, no indexer call
    assert [int(c.get("max_pairs", FULL)) for c in PxfFake.calls] == [50_000, FULL]
    (call,) = ReuseFake.fit_calls
    assert call["idx"] == 3 and call["sigma0"] == 0.5 and call["max_evals"] == 200
    assert call["draw"] == float(np.random.default_rng(0).random())  # a fresh default_rng(0)
    assert ReuseFake.cost_calls == [(3, pytest.approx(0.032))]  # one strict cost at the voxel
    e0, e3 = _entry(bfs, 0), _entry(bfs, 3)
    assert e0["reuse"]["tried"] is False and e0["reuse"]["starts"] == []
    assert e0["used"] is True and e0["evals"] == 100 + 150
    r = e3["reuse"]
    assert r["tried"] and r["success"] and r["second"] is False and r["starts"] == [0]
    assert [f["evals"] for f in r["fits"]] == [201] and r["evals"] == 1 + 201
    assert r["fits"][0]["hit"] == pytest.approx(0.97) and r["costs"] == [pytest.approx(0.032)]
    assert e3["used"] is True and e3["taken"] == -1 and e3["n_cand"] == 0 and e3["n_att"] == 0
    assert e3["evals"] == 202 and e3["strict_costs"] == [] and e3["refined"] == []
    assert e3["hit_ratio"] == pytest.approx(0.97) and "fallback" not in e3
    assert e3["ray_index"]["n_pairs"] == 0 and e3["ray_index"]["fallback_full"] is False
    assert np.allclose(mic.voxels[3].orientation, PxfFake.truth[3])
    assert mic.voxels[3].reconstruction_id == FITTED
    s = bfs.stats
    assert s["n_evals_seed"] == 250 + 202 and bfs.records[3].n_evals >= 202
    assert (s["reuse_steps"], s["reuse_tried"], s["reuse_success"]) == (2, 1, 1)
    assert s["reuse_second_fits"] == 0 and s["reuse_evals"] == 202 and s["reuse_evals_wasted"] == 0
    # phase A's accounting: a step ended by a try is an index run and an index use
    assert s["index_runs"] == 2 and s["index_used"] == 2 and s["index_evals"] == 250 + 202
    assert s["index_evals_used"] == 250 + 202
    assert s.get("index_evals_wasted", 0) == 0


def test_a_failure_runs_the_step_unchanged_with_the_try_added(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    ref, _ = _reuse(monkeypatch, PXF)  # key off: voxel 3's step alone
    ref3 = _entry(ref, 3)
    bfs, mic = _reuse(monkeypatch, ON, capture=1.0)  # the fit stays at voxel 0's orientation
    e3 = _entry(bfs, 3)
    r = e3["reuse"]
    assert r["tried"] and not r["success"] and r["starts"] == [0] and len(r["fits"]) == 1
    assert r["evals"] == 202
    # the indexing ran as without the key: same pairs, same steps, evaluations add up
    for key in ("n_cand", "n_att", "taken", "strict_costs", "refined", "finisher_costs",
                "finisher_evals", "used", "hit_ratio"):  # fmt: skip
        assert e3[key] == ref3[key], key
    assert e3["evals"] == ref3["evals"] + 202
    assert e3["ray_index"]["first"] == ref3["ray_index"]["first"]
    assert e3["ray_index"]["n_pairs"] == ref3["ray_index"]["n_pairs"]
    s, s0 = bfs.stats, ref.stats
    assert s["n_evals_seed"] == s0["n_evals_seed"] + 202
    assert s["reuse_success"] == 0 and s["reuse_evals"] == 202 and s["reuse_evals_wasted"] == 202
    assert s["index_evals"] == s0["index_evals"] + 202 and s["index_used"] == s0["index_used"]
    assert np.allclose(mic.voxels[3].orientation, PxfFake.truth[3])


def test_the_try_runs_once_when_the_full_sample_step_repeats(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    bfs, _ = _reuse(monkeypatch, ON, capture=1.0)
    # voxel 0: 50,000 then 600,000 pairs; voxel 3: the same two steps after ONE try
    assert [int(c.get("max_pairs", FULL)) for c in PxfFake.calls] == [50_000, FULL] * 2
    assert len(ReuseFake.fit_calls) == 1 and len(ReuseFake.cost_calls) == 1
    e3 = _entry(bfs, 3)
    assert e3["ray_index"]["fallback_full"] is True
    assert e3["evals"] == 202 + 100 + 150  # try, then both PXF steps
    assert bfs.stats["reuse_tried"] == 1 and bfs.stats["index_fb_full_steps"] == 2


def test_a_second_fit_starts_from_a_far_enough_voxel(monkeypatch: pytest.MonkeyPatch) -> None:
    """Two accepted voxels, 2 degrees apart. The scripted strict costs rank voxel 0's orientation
    first, but it is 2.1 degrees from voxel 2's truth (capture 0.5), so the first fit fails; the
    second start is voxel 1's orientation (more than 0.5 degree from the first), 0.1 degree from
    the truth, and its fit succeeds."""
    import icenine.reconstructor as RM
    from icenine.mic_file import MicFile, Voxel
    from test_bfs_refit import POSITIONS

    ReuseFake.capture, ReuseFake.fit_calls, ReuseFake.cost_calls = 0.5, [], []
    ReuseFake.sym = np.eye(3)[None]
    ReuseFake.scripted_costs = [0.01, 0.02]
    _reset_rescue(failing=(), rescuable=())
    PxfFake.calls, PxfFake.available, PxfFake.good_from, PxfFake.no_cand_below = (
        [],
        700_000,
        FULL,
        0,
    )
    ReuseFake.truth = [_rz(0.0), _rz(2.0), _rz(2.1), _rz(0.0)]
    voxels = [
        Voxel(position=np.array([x, y, 0.0]), orientation=np.eye(3), side_length=1.0)
        for x, y in POSITIONS[:4]
    ]
    mic = MicFile(voxels)
    monkeypatch.setattr(RM, "AdaptiveVoxelReconstructor", ReuseFake)
    setup = SimpleNamespace(
        sample=SimpleNamespace(get_mic=lambda: mic),
        config=SimpleNamespace(min_acceleration_threshold=0.8),
        search_params=SearchParameters(**ON),
    )
    bfs = RM.BFSReconstruction(setup)  # type: ignore[arg-type]
    bfs.stats = bfs._new_stats()
    for j in (0, 1):  # accepted with their truths; voxels 2 and 3 are not
        mic.voxels[j].orientation = np.asarray(ReuseFake.truth[j], np.float64)
        mic.voxels[j].reconstruction_id = FITTED
    rec = SimpleNamespace(voxel_idx=2, n_evals=0, wall_s=0.0)
    try:
        out = bfs._index_seed(mic.voxels[2], [(0.75, 1.3, 0.0)], rec, False)  # type: ignore
    finally:
        ReuseFake.scripted_costs = None
    assert out is not None and PxfFake.calls == []  # no indexer call
    row = bfs.index_log[-1]["reuse"]
    assert row["starts"] == [0, 1] and row["costs"] == [0.01, 0.02]
    assert [(f["start"], f["second"]) for f in row["fits"]] == [(0, False), (1, True)]
    assert row["second"] is True and row["success"] is True
    assert row["evals"] == 2 + 201 + 201 and bfs.stats["reuse_second_fits"] == 1
    assert np.allclose(out[0].orientation, ReuseFake.truth[2])


# ------------------------------------------------------------------------------ 4. real images
def test_the_fit_overlap_equals_a_fresh_evaluate_overlap() -> None:
    from test_findoptimal_refactor import _build
    from test_seed_by_indexing import needs_c

    if needs_c.args[0]:
        pytest.skip("C extension missing")
    rec, voxel, R_true = _build()
    from scipy.spatial.transform import Rotation

    from icenine.reconstructor import _get_voxel_vertices

    verts = _get_voxel_vertices(voxel)
    R0 = Rotation.from_rotvec(np.radians([0.3, -0.2, 0.1])).as_matrix() @ R_true
    r = rec.local_optimization(
        voxel_vertices=verts,
        phase_index=voxel.phase,
        initial_orientation=R0,
        rng=np.random.default_rng(0),
        cma_sigma0_deg=RI.REUSE_SIGMA0_DEG,
        cma_max_evals=RI.REUSE_MAX_EVALS,
        optimizer="cma",
    )
    fresh = rec.evaluate_overlap(np.asarray(r.orientation, np.float64), verts, voxel.phase)
    from icenine.reconstructor import BFSReconstruction

    assert BFSReconstruction._ratios(r.overlap_info) == BFSReconstruction._ratios(fresh)
    assert r.overlap_info.pixel_overlap == fresh.pixel_overlap and r.cost == fresh.cost


def test_the_block_on_the_real_toy_ends_on_a_nearby_orientation(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Real images (ThreeVoxels): voxels 1 and 2 are marked accepted with orientations 0.3 and
    0.5 degree from voxel 0's truth; the try at voxel 0 costs both (2 evaluations), fits from the
    cheaper one and, with the real overlap, ends at the accepted hit ratio: no indexing. The
    returned overlap is the same evaluation a fresh ``evaluate_overlap`` gives."""
    import icenine.reconstructor as RM
    from scipy.spatial.transform import Rotation
    from test_findoptimal_refactor import _build
    from test_seed_by_indexing import needs_c

    if needs_c.args[0]:
        pytest.skip("C extension missing")
    rec, voxel, R_true = _build()
    mic = rec.setup.sample.get_mic()
    for j, deg in ((1, 0.3), (2, 0.5)):
        R = Rotation.from_rotvec(np.radians([deg, 0.0, 0.0])).as_matrix() @ R_true
        mic.voxels[j].orientation = R
        mic.voxels[j].reconstruction_id = FITTED
    assert mic.voxels[0].reconstruction_id != FITTED
    monkeypatch.setattr(RM, "AdaptiveVoxelReconstructor", lambda setup: rec)
    bfs = RM.BFSReconstruction(rec.setup)
    verts = RM._get_voxel_vertices(voxel)
    end, row = bfs._reuse_try(0, verts, voxel.phase)
    assert row["tried"] and row["starts"] == [1, 2] and len(row["costs"]) == 2
    assert row["evals"] == 2 + sum(f["evals"] for f in row["fits"])
    assert row["success"] and end is not None and row["fits"][-1]["hit"] >= RI.REUSE_ACCEPT
    fresh = rec.evaluate_overlap(end[0], verts, voxel.phase)
    assert RM.BFSReconstruction._ratios(end[2]) == RM.BFSReconstruction._ratios(fresh)
    assert row["fits"][-1]["hit"] == RM.BFSReconstruction._ratios(fresh)[0]
    sym = RI.proper_rotations(rec.setup.exp_setup.get_sample_symmetry())
    assert RI.misorientation_deg(end[0][None], R_true[None], sym)[0] < 0.5
    # the crystal symmetry the block uses is the one the scripts used
    assert np.array_equal(rec._index_context(voxel.phase)[1].sym, sym)


# ------------------------------------------------------------------------------ 5. hexagonal
def test_hexagonal_symmetry_runs_through_the_library_block(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    from pymatgen.core.lattice import Lattice
    from pymatgen.core.structure import Structure

    from icenine.symmetry import CrystalSymmetry

    st = Structure(
        Lattice.hexagonal(2.95, 4.68), ["Ti", "Ti"], [[1 / 3, 2 / 3, 0.25], [2 / 3, 1 / 3, 0.75]]
    )
    sym = RI.proper_rotations(CrystalSymmetry(st))
    assert len(sym) == 12
    bfs, _ = _reuse(monkeypatch, ON, sym=sym)
    assert _entry(bfs, 3)["reuse"]["success"] is True and bfs.stats["reuse_success"] == 1
    # a rotation of 60 degrees about the hexagonal axis is a symmetry copy for the second start
    ez = np.array([0.0, 0.0, 1.0])
    six = [s for s in sym if abs(np.trace(s) - 2.0) < 1e-9 and np.allclose(s @ ez, ez)]
    assert six  # the 60 degree rotations about z
    a, b = np.eye(3), six[0]
    assert RI.misorientation_deg(a[None], b[None], sym)[0] < 1e-6
    assert RI.second_start([0, 1], [a, b], lambda x, y: float(
        RI.misorientation_deg(x[None], y[None], sym)[0])) is None  # fmt: skip
