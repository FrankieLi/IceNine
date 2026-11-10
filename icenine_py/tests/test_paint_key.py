"""PaintGrains (autoresearch AR3 step 2, opt-in, default off): ``grain_paint.paint_grains`` runs in
place of the BFS.

* the key parses and is off by default; the two refused settings raise;
* with the key off the paint loop is never called and a small BFS run is unchanged;
* with the key on, a small real case (ThreeVoxels, real images) gives the same map and the same
  evaluations by stage as a direct call of ``paint_grains`` with the same order and generator, and
  the stages sum to the evaluations the cost functions counted.
"""

from pathlib import Path
from types import SimpleNamespace
from typing import Any, Dict, List

import numpy as np
import pytest

import icenine.grain_paint as GP
from icenine import cost_functions as CF
from icenine.config_file import ConfigFile
from icenine.mic_file import MicFile
from icenine.orientation_search import SearchParameters
from icenine.reconstructor import BFSReconstruction, check_paint_grains_keys

# tests/ is on sys.path under pytest
from test_bfs_refit import _run

needs_c = pytest.mark.skipif(not CF._HAS_C_RASTERIZE, reason="C extension missing")
_CFG = Path(__file__).parent.parent.parent / "Examples" / "Example2.ThreeVoxels"
_CFG = _CFG / "ConfigFiles" / "ReconstructQ8.config"


def _config_with(tmp_path: Path, extra: str) -> ConfigFile:
    if not _CFG.exists():
        pytest.skip("ThreeVoxels config not available")
    p = tmp_path / "x.config"
    p.write_text(_CFG.read_text() + "\n" + extra + "\n")
    return ConfigFile.from_file(str(p))


# ------------------------------------------------------------------------------ the key
def test_paint_grains_key_default_off_and_parsed(tmp_path: Path) -> None:
    assert SearchParameters().paint_grains is False
    assert SearchParameters.from_config(_config_with(tmp_path, "")).paint_grains is False
    on = SearchParameters.from_config(_config_with(tmp_path, "SeedByIndexing 1\nPaintGrains 1"))
    off = SearchParameters.from_config(_config_with(tmp_path, "PaintGrains 0"))
    assert on.paint_grains is True and off.paint_grains is False


@pytest.mark.parametrize("line", ["PaintGrains 2", "PaintGrains"])
def test_paint_grains_key_bad_values(tmp_path: Path, line: str) -> None:
    with pytest.raises(ValueError):
        _config_with(tmp_path, line)


def test_refused_settings_name_the_keys() -> None:
    with pytest.raises(ValueError, match="PaintGrains.*SeedByIndexing"):
        check_paint_grains_keys(SearchParameters(paint_grains=True))
    with pytest.raises(ValueError, match="PaintGrains.*OutwardSearch"):
        check_paint_grains_keys(
            SearchParameters(paint_grains=True, seed_by_indexing=True, outward_search=True)
        )
    check_paint_grains_keys(SearchParameters(paint_grains=True, seed_by_indexing=True))
    check_paint_grains_keys(SearchParameters())
    with pytest.raises(ValueError, match="SeedByIndexing"):
        fake: Any = SimpleNamespace(search_params=SearchParameters(paint_grains=True))
        BFSReconstruction(fake)


# ------------------------------------------------------------------------------ key off
def _snapshot(bfs: Any, mic: Any) -> Dict[str, Any]:
    stats = {k: v for k, v in bfs.stats.items() if "wall" not in k}
    return dict(
        stats=stats,
        sources=[bfs.records[i].source for i in sorted(bfs.records)],
        n_evals=[bfs.records[i].n_evals for i in sorted(bfs.records)],
        orientations=[np.asarray(v.orientation).tolist() for v in mic.voxels],
        states=[v.reconstruction_id for v in mic.voxels],
    )


def test_key_off_never_calls_the_paint_loop_and_changes_nothing(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    runs = []
    for params in ({}, {"paint_grains": False}):
        bfs, mic = _run((3.0, 3.5, 3.2), monkeypatch, **params)
        runs.append(_snapshot(bfs, mic))
        assert bfs.paint_result is None and "paint_evals" not in bfs.stats
    assert runs[0] == runs[1]

    def boom(*a: Any, **k: Any) -> Any:  # pragma: no cover - must not run
        raise AssertionError("the paint loop ran with PaintGrains off")

    monkeypatch.setattr(GP, "paint_grains", boom)
    bfs, mic = _run((3.0, 3.5, 3.2), monkeypatch, paint_grains=False)
    assert _snapshot(bfs, mic) == runs[0]


# ------------------------------------------------------------------------------ key on
@pytest.fixture(scope="module")
def toy() -> Any:
    from test_findoptimal_refactor import _build

    return _build()


@needs_c
def test_key_on_equals_a_direct_call(toy: Any) -> None:
    """The key runs the loop of ``paint_grains`` through ``ReconstructorHooks`` with the BFS's
    order and generator: same map and same evaluations by stage as a direct call, and the stages
    sum to the evaluations counted by the cost function and the seed paths."""
    rec, _voxel, _R = toy
    setup = rec.setup
    mic = setup.sample.get_mic()
    n = len(mic.voxels)
    p = setup.search_params
    old = (p.seed_by_indexing, p.paint_grains)
    p.seed_by_indexing, p.paint_grains = True, True
    try:
        # direct call: a fresh reconstructor (the key's own), the same order and generator
        bfs_d = BFSReconstruction(setup)
        rng = np.random.default_rng(3)
        order: List[int] = list(range(n))
        rng.shuffle(order)
        mic_d = MicFile.read(str(Path(_CFG).parents[1] / "SimInput" / "three_voxels.mic"))
        hooks = GP.ReconstructorHooks(bfs_d, mic_d, rng)
        direct = GP.paint_grains(order, n, hooks)

        bfs = BFSReconstruction(setup)
        done = bfs.reconstruct_sample(rng=np.random.default_rng(3))
        assert sorted(done) == list(range(n))
        res = bfs.paint_result
        assert res is not None and bfs.stats["paint_evals"] == res.evals
        assert res.evals == direct.evals and res.counts == direct.counts
        for v in range(n):
            Rk, src = res.final_orientation(v)
            Rd, srcd = direct.final_orientation(v)
            assert src == srcd and np.array_equal(Rk, Rd)
            assert np.array_equal(mic.voxels[v].orientation, Rk)
        # the stages sum to the cost functions' and seed paths' own counts
        seed_paths = sum(r["evals"] for r in bfs.paint_hooks.seed_records)
        assert res.evals["score"] == bfs.paint_hooks.cost_fn.eval_count
        assert res.evals["seed_index"] + res.evals["index_rejected"] == seed_paths
        assert bfs.stats["paint_evals_total"] == sum(res.evals.values())
        assert all(bfs.records[v].source in ("seed", "paint", "fit", "failed") for v in range(n))
    finally:
        p.seed_by_indexing, p.paint_grains = old
