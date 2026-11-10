"""SeedByIndexing (autoresearch AR2 step 2, opt-in, default off): BFS seeds and reseeds are
indexed from the observed rays first; the grid search is the fallback.

Two kinds of test, as for OutwardSearch:
- The scripted voxel groups of ``test_bfs_refit`` (fake reconstructor, no images) with a fake
  ``index_seed``: the key's bookkeeping (what is taken, what falls back, the counts), and that the
  key off (absent or written off) never calls it and leaves orientations, evaluation counts,
  records and stats unchanged.
- ``AdaptiveVoxelReconstructor.index_seed`` on the ThreeVoxels toy problem (real images): it equals
  the script path of ``scripts/autoresearch/ar2/holdout_run.py`` written out by hand.
"""

from dataclasses import asdict
from pathlib import Path
from typing import Any, Dict, List, Tuple

import numpy as np
import pytest

import icenine.mic_file as mic_module
import icenine.reconstructor as reconstructor_module
from icenine import cost_functions as CF
from icenine import ray_indexing as RI
from icenine.config_file import ConfigFile
from icenine.mic_file import ReconstructionState
from icenine.orientation_search import SearchCandidate, SearchParameters
from icenine.reconstructor import _get_voxel_vertices

# tests/ is on sys.path under pytest
import bfs_scripted
from test_bfs_refit import FakeReconstructor, RejectingFull, _info, _run, _rz, _states

FITTED, REFIT = ReconstructionState.FITTED, ReconstructionState.REFIT
needs_c = pytest.mark.skipif(not CF._HAS_C_RASTERIZE, reason="C extension missing")


class IndexFake(RejectingFull):
    """A full search is logged (and rejected for ``bad``) as in ``RejectingFull``. ``index_seed``
    returns the voxel's truth with 1,234 evaluations; for the voxels in ``no_cand`` it returns no
    candidate (no evaluation spent); the seed test (a fake overlap) fails for the voxels in
    ``failing`` until a grid search has run."""

    index_calls: List[int] = []
    no_cand: Tuple[int, ...] = ()
    failing: Tuple[int, ...] = ()

    def index_seed(self, voxel_vertices: Any, phase_index: int) -> RI.SeedResult:
        i = self._idx(voxel_vertices)
        IndexFake.index_calls.append(i)
        if i in self.no_cand:
            return RI.SeedResult(n_att=3, evals=0, index_s=0.01, total_s=0.01)
        R = self.truth[i].astype(np.float64)
        return RI.SeedResult(
            orientation=R, cost=0.05, taken=1, n_cand=3, n_att=40, evals=1234, index_s=0.02,
            total_s=0.5, strict_costs=[0.3, 0.1], refined=[1, 0], finisher_costs=[0.05, 0.2],
            finisher_evals=[600, 600],
        )  # fmt: skip

    def evaluate_overlap(self, R: np.ndarray, vertices: Any, phase: int) -> Any:
        if self._idx(vertices) in self.failing and not getattr(self, "_grid_last", False):
            return _info(0.3)
        return super().evaluate_overlap(R, vertices, phase)

    def reconstruct_voxel(self, voxel_vertices: Any, phase_index: int, rng: Any) -> Any:
        self._grid_last = True
        return super().reconstruct_voxel(voxel_vertices, phase_index, rng)


def _reset(bad: Tuple[int, ...] = (), no_cand: Tuple[int, ...] = (), failing: Tuple[int, ...] = ()):
    IndexFake.index_calls = []
    IndexFake.no_cand, IndexFake.failing = no_cand, failing
    RejectingFull.bad, RejectingFull.full_calls = bad, []


class NoIndex(RejectingFull):
    def index_seed(self, *a: Any, **k: Any) -> Any:  # pragma: no cover - must not run
        raise AssertionError("index_seed ran with SeedByIndexing off")


# ------------------------------------------------------------------------------ config key
_CFG = Path(__file__).parent.parent.parent / "Examples" / "Example2.ThreeVoxels"
_CFG = _CFG / "ConfigFiles" / "ReconstructQ8.config"


def _config_with(tmp_path: Path, extra: str) -> ConfigFile:
    if not _CFG.exists():
        pytest.skip("ThreeVoxels config not available")
    p = tmp_path / "x.config"
    p.write_text(_CFG.read_text() + "\n" + extra + "\n")
    return ConfigFile.from_file(str(p))


def test_seed_by_indexing_key_default_off_and_parsed(tmp_path: Path) -> None:
    assert SearchParameters().seed_by_indexing is False
    assert SearchParameters.from_config(_config_with(tmp_path, "")).seed_by_indexing is False
    on = SearchParameters.from_config(_config_with(tmp_path, "SeedByIndexing 1"))
    off = SearchParameters.from_config(_config_with(tmp_path, "SeedByIndexing 0"))
    assert on.seed_by_indexing is True and off.seed_by_indexing is False


@pytest.mark.parametrize("line", ["SeedByIndexing 2", "SeedByIndexing"])
def test_seed_by_indexing_key_bad_values(tmp_path: Path, line: str) -> None:
    with pytest.raises(ValueError):
        _config_with(tmp_path, line)


# ------------------------------------------------------------------------------ key off
def _fingerprint(bfs: Any, mic: Any) -> Dict[str, Any]:
    stats = {k: v for k, v in bfs.stats.items() if "wall" not in k and k != "reseeds"}
    records = {}
    for k, r in bfs.records.items():
        d = asdict(r)
        d.pop("wall_s")
        records[k] = d
    return dict(
        stats=stats,
        records=records,
        orientations=[np.asarray(v.orientation).tolist() for v in mic.voxels],
        states=_states(mic),
        reseeds=[
            {k: v for k, v in e.items() if "seconds" not in k} for e in bfs.stats.get("reseeds", [])
        ],
    )


@pytest.mark.parametrize("extra", [{}, {"bfs_reseed": True}])
def test_key_absent_and_key_off_are_identical_and_never_index(
    monkeypatch: pytest.MonkeyPatch, extra: Dict[str, Any]
) -> None:
    """Same orientations, evaluation counts, records and stats with the key absent and written off;
    ``index_seed`` is never called and no ``index_*`` counter exists."""
    runs = []
    for params in ({}, {"seed_by_indexing": False}):
        _reset(bad=(1,) if extra else ())
        bfs, mic = _run((3.0, 3.5, 3.2), monkeypatch, fake=NoIndex, **extra, **params)
        assert bfs.index_log == []
        assert not any(k.startswith("index_") for k in bfs.stats)
        runs.append(_fingerprint(bfs, mic))
    assert runs[0] == runs[1]
    assert runs[0]["stats"]["n_evals_seed"] > 0


def test_key_off_scripted_grid_still_matches_the_golden() -> None:
    """The scripted 14 x 14 grid (golden recorded before the Phase D options) with the key
    written explicitly off: the whole fit sequence is unchanged."""
    gold = np.load(Path(__file__).parent / "data" / "bfs_grid_golden.npz")
    out = bfs_scripted.run_grid(
        reconstructor_module, mic_module, SearchCandidate, SearchParameters,
        {"seed_by_indexing": False},
    )  # fmt: skip
    assert np.array_equal(out["done"], gold["mc_done"])
    assert np.array_equal(out["log_idx"], gold["mc_log_idx"])
    assert float(out["rng_after"]) == float(gold["mc_rng_after"])


# ------------------------------------------------------------------------------ key on
def test_on_takes_the_indexing_result_and_skips_the_grid_search(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Voxel 0 is a seed: it is indexed, the result (the truth) passes the seed test, the grid
    search never runs, the seed's evaluations are the indexing path's, and the BFS goes on."""
    _reset()
    bfs, mic = _run((0.2, 0.3, 0.2), monkeypatch, fake=IndexFake, seed_by_indexing=True)
    assert _states(mic) == [FITTED] * 4
    assert IndexFake.index_calls == [0] and RejectingFull.full_calls == []
    s = bfs.stats
    assert s["n_seeds"] == 1 and s["n_evals_seed"] == 1234
    assert s["index_runs"] == 1 and s["index_used"] == 1 and s["index_evals"] == 1234
    assert s["index_evals_used"] == 1234 and "index_evals_wasted" not in s
    assert "index_fallback_no_candidate" not in s and "index_fallback_seed_test" not in s
    assert bfs.records[0].n_evals == 1234 and bfs.records[0].source == "seed"
    (e,) = bfs.index_log
    assert e["voxel"] == 0 and e["phase"] == "seed" and e["used"] and e["taken"] == 1
    assert e["hit_ratio"] is not None and e["refined"] == [1, 0]


def test_on_without_a_candidate_falls_back_to_the_grid_search(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """No candidate: the grid search runs as today, its evaluations are counted as the fallback's
    and the indexing path's own evaluations as wasted."""
    _reset(no_cand=(0,))
    bfs, mic = _run((0.2, 0.3, 0.2), monkeypatch, fake=IndexFake, seed_by_indexing=True)
    assert _states(mic) == [FITTED] * 4
    assert IndexFake.index_calls == [0] and RejectingFull.full_calls == [0]
    s = bfs.stats
    assert s["index_fallback_no_candidate"] == 1 and "index_used" not in s
    assert s["index_fallback_grid_evals"] == 1000 and s["n_evals_seed"] == 1000
    (e,) = bfs.index_log
    assert e["fallback"] == "no_candidate" and not e["used"] and e["grid_evals"] == 1000


def test_on_failing_the_seed_test_falls_back_and_counts_the_wasted_evaluations(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """The indexing result fails the seed test (hit ratio 0.3 < 0.8): the grid search runs, its
    result is taken, and the 1,234 evaluations of the indexing path are in the seed total."""
    _reset(failing=(0,))
    bfs, mic = _run((0.2, 0.3, 0.2), monkeypatch, fake=IndexFake, seed_by_indexing=True)
    assert _states(mic) == [FITTED] * 4
    assert RejectingFull.full_calls == [0]
    s = bfs.stats
    assert s["index_fallback_seed_test"] == 1 and "index_used" not in s
    assert s["index_evals_wasted"] == 1234 and s["index_fallback_grid_evals"] == 1000
    assert s["n_evals_seed"] == 1234 + 1000 and bfs.records[0].n_evals == 1234 + 1000
    (e,) = bfs.index_log
    assert e["fallback"] == "seed_test" and e["hit_ratio"] == pytest.approx(0.3)


def test_on_indexes_reseeds_and_keeps_their_cost_apart(monkeypatch: pytest.MonkeyPatch) -> None:
    """Voxels 1, 2 are rejected from seed 0 (3 deg), voxel 3 is a seed. The reseed of the cluster
    {1, 2} is indexed (phase "reseed"); its evaluations go to the reseed bucket, not the seeds'."""
    _reset()
    bfs, mic = _run(
        (3.0, 3.5, 3.2), monkeypatch, fake=IndexFake, seed_by_indexing=True, bfs_reseed=True
    )
    assert _states(mic) == [FITTED] * 4
    assert RejectingFull.full_calls == []
    assert sorted(IndexFake.index_calls) == [0, 1, 3]
    s = bfs.stats
    assert s["n_reseeds"] == 1 and s["n_reseed_accepted"] == 1
    assert s["n_evals_seed"] == 2 * 1234 and s["n_evals_reseed"] == 1234
    assert [e["phase"] for e in bfs.index_log] == ["seed", "seed", "reseed"]
    assert s["index_runs"] == 3 and s["index_used"] == 3


# ------------------------------------------------------------------------------ real images
@pytest.fixture(scope="module")
def toy() -> Any:
    from test_findoptimal_refactor import _build

    return _build()


@needs_c
def test_library_path_equals_the_script_path_on_the_toy(toy: Any) -> None:
    """``index_seed`` equals the path of ``holdout_run.py`` written out by hand: index, strict cost
    on the first 5 candidates, ``refine_from_candidates`` (rng seed 0) on the 2 lowest, the lowest
    final cost; and its evaluations are 5 + the finisher runs'."""
    rec, voxel, R_true = toy
    verts = _get_voxel_vertices(voxel)
    lib = rec.index_seed(verts, voxel.phase)
    assert rec.index_pairing_s > 0 and lib.orientation is not None

    sc, crystal, p, _ = rec._index_context(voxel.phase)
    res, _att = RI.index_voxel(sc, verts, crystal, p, 20)
    lf = rec._make_local_cost_fn()
    k5 = min(5, len(res.R))
    costs = np.array([float(lf.evaluate(res.R[j], verts, voxel.phase).cost) for j in range(k5)])
    order = [int(j) for j in np.argsort(costs, kind="stable")[:2]]
    fins = []
    for j in order:
        r = rec.refine_from_candidates(
            [SearchCandidate(orientation=res.R[j].copy(), cost=float(costs[j]))],
            verts, voxel.phase, rng=np.random.default_rng(0),
        )  # fmt: skip
        fins.append(dict(cand=j, cost=float(r.cost), R=np.asarray(r.orientation),
                         evals=int(rec.last_eval_counts[1])))  # fmt: skip
    best = min(fins, key=lambda f: f["cost"])
    assert np.array_equal(lib.candidates, res.R)
    assert lib.strict_costs == [float(c) for c in costs] and lib.refined == order
    assert lib.taken == best["cand"] and lib.cost == best["cost"]
    assert np.array_equal(lib.orientation, best["R"])
    assert lib.evals == k5 + sum(f["evals"] for f in fins)
    # the key's result is within the toy's resolution of the truth (symmetry-reduced)
    sym = RI.proper_rotations(rec.setup.exp_setup.get_sample_symmetry())
    assert RI.misorientation_deg(lib.orientation[None], R_true[None], sym)[0] < 1.0


@needs_c
def test_no_candidate_when_the_voxel_has_no_rays(toy: Any, monkeypatch: pytest.MonkeyPatch) -> None:
    """Fewer than two attached rays: no candidate, no cost evaluation."""
    rec, voxel, _ = toy
    verts = _get_voxel_vertices(voxel)
    sc, *_ = rec._index_context(voxel.phase)
    monkeypatch.setattr(sc, "attached_rays", lambda v: np.zeros(1, dtype=int))
    out = rec.index_seed(verts, voxel.phase)
    assert out.orientation is None and out.evals == 0 and out.n_cand == 0
