"""AR5 step 20, phase B: the opt-in key ``BoundaryPassRepair`` (off by default).

After the boundary pass, each contested voxel gets a 120-evaluation CMA local optimization from its
own end and from each candidate, with a generator of its own, and keeps the best end. Scripted
groups with the fake reconstructor of ``test_bfs_ar3_keys`` (no images): defaults and parsing,
invalid combinations, the flag functions equal the script rule's (step 19's frozen flag), key off
equals key absent, key on repairs a contested voxel and counts every evaluation, per-fit
generators, and the call from the BFS.
"""

import sys
from pathlib import Path
from types import SimpleNamespace
from typing import Any, Dict, List

import numpy as np
import pytest

import icenine.reconstructor as reconstructor_module
from icenine import boundary_pass as BP
from icenine import boundary_repair as BR
from icenine.orientation_search import SearchCandidate, SearchParameters
from icenine.reconstructor import BFSReconstruction

# tests/ is on sys.path under pytest
from test_bfs_ar3_keys import IDENT, PassFake, _info, _pass_bfs
from test_bfs_refit import FakeReconstructor, _angle, _run, _rz
from test_seed_by_indexing import _config_with

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "scripts" / "autoresearch" / "ar5"))
import cmarepair_rules as CR  # noqa: E402

KEY = "BoundaryPassRepair"


# ------------------------------------------------------------------------------ 1. parsing
def test_default_off_and_parsed(tmp_path: Path) -> None:
    assert SearchParameters().boundary_pass_repair is False
    assert SearchParameters.from_config(_config_with(tmp_path, "")).boundary_pass_repair is False
    on = SearchParameters.from_config(_config_with(tmp_path, f"BoundaryPass 1\n{KEY} 1"))
    assert on.boundary_pass_repair is True and on.boundary_pass is True
    off = SearchParameters.from_config(_config_with(tmp_path, f"BoundaryPass 1\n{KEY} 0"))
    assert off.boundary_pass_repair is False


def test_bad_value_raises_at_parse_time(tmp_path: Path) -> None:
    with pytest.raises(ValueError):
        _config_with(tmp_path, f"{KEY} maybe")


@pytest.mark.parametrize(
    "kw",
    [
        dict(boundary_pass_repair=True),  # needs BoundaryPass 1
        dict(
            seed_by_indexing=True, paint_grains=True, boundary_pass=True, boundary_pass_repair=True
        ),
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


def test_module_constants_are_the_plans_and_the_scripts() -> None:
    assert (BR.RATIO_MAX, BR.MARGIN_MAX, BR.FAR_DEG) == (0.95, 0.15, 1.0)
    assert (BR.SIGMA0_DEG, BR.MAX_EVALS, BR.SEED) == (0.05, 120, 0)
    assert (BR.RATIO_MAX, BR.MARGIN_MAX, BR.FAR_DEG) == (CR.RATIO_MAX, CR.MARGIN_MAX, CR.FAR_DEG)
    assert (BR.SIGMA0_DEG, BR.MAX_EVALS, BR.SEED) == (
        CR.CMA_SIGMA0_DEG,
        CR.CMA_MAX_EVALS,
        CR.REPAIR_SEED,
    )


# ------------------------------------------------------------------------------ 2. flag ports
def test_the_flag_functions_equal_the_script_rules() -> None:
    rng = np.random.default_rng(5)
    n = 60
    hit = rng.uniform(0.3, 1.0, n)
    other = np.where(rng.random(n) < 0.3, -1.0, rng.uniform(0.2, 1.0, n))
    nbrs = [
        sorted(int(j) for j in rng.choice(n, size=rng.integers(0, 6), replace=False))
        for _ in range(n)
    ]
    assert np.array_equal(BR.margin_of(hit, other), CR.R19.margin_of(hit, other))
    assert np.array_equal(BR.ratio_of(hit, nbrs), CR.R19.ratio_of(hit, nbrs))
    for a, b in zip(BR.contested(hit, other, nbrs), CR.contested(hit, other, nbrs)):
        assert np.array_equal(a, b)
    # other_end_hit on random maps and tries
    R = np.stack([_rz(d).astype(np.float64) for d in rng.uniform(0, 5, n)])
    Rb = np.stack([_rz(d).astype(np.float64) for d in rng.uniform(0, 5, n)])
    tv = rng.integers(0, n, 40)
    te = np.stack([_rz(d).astype(np.float64) for d in rng.uniform(0, 5, 40)])
    th = rng.uniform(0.2, 1.0, 40)
    y = np.eye(3)[None]
    assert np.array_equal(
        BR.other_end_hit(R, Rb, hit, tv, te, th, y), CR.other_end_hit(R, Rb, hit, tv, te, th, y)
    )


def test_per_fit_generator_equals_the_scripts() -> None:
    for v, c in ((7, None), (7, 3), (0, 0), (12, 11)):
        assert np.array_equal(BR.fit_rng(v, c).random(4), CR.fit_rng(v, c).random(4))


# ------------------------------------------------------------------------------ 3. the pass
class RepairFake(PassFake):
    """``PassFake`` whose repair fits (budget 120) return the truth, as a good local optimization
    would, with 121 evaluations; it records each fit's generator state. Truth of voxel 0: 10
    degrees (its candidates, at 2 degrees, cannot reach it in the pass: a try within 3 degrees
    only)."""

    TRUTH_DEG = [10.0, 2.0, 2.0, 2.0]
    states: List[Any] = []

    def local_optimization(
        self, voxel_vertices: Any, phase_index: int, initial_orientation: np.ndarray,
        rng: Any = None, cma_sigma0_deg: Any = None, cma_max_evals: Any = None,
        optimizer: Any = None,
    ) -> SearchCandidate:  # fmt: skip
        if cma_max_evals != BR.MAX_EVALS:
            return super().local_optimization(
                voxel_vertices, phase_index, initial_orientation, rng, cma_sigma0_deg,
                cma_max_evals, optimizer,
            )  # fmt: skip
        i = self._idx(voxel_vertices)
        RepairFake.states.append((i, rng.bit_generator.state))
        R = _rz(self.TRUTH_DEG[i]).astype(np.float64)
        self.last_local_optimization_evals = BR.MAX_EVALS + 1
        return SearchCandidate(R, 1.0 - self._hit(R, i), _info(self._hit(R, i)))


def _repair_bfs(monkeypatch: pytest.MonkeyPatch, deg: List[float]) -> Any:
    bfs, mic = _pass_bfs(monkeypatch, deg)
    monkeypatch.setattr(reconstructor_module, "AdaptiveVoxelReconstructor", RepairFake)
    bfs.reconstructor = RepairFake(None)
    RepairFake.fits = []
    RepairFake.states = []
    return bfs, mic


def _summary(res: BP.BoundaryPassResult, mic: Any) -> Dict[str, Any]:
    return dict(
        boundary_voxels=res.boundary_voxels,
        tries=[(t["voxel"], t["cand"], t["evals"], t["cost"], t["hit"]) for t in res.tries],
        old_score=res.old_score,
        changed=res.changed,
        n_evals_tries=res.n_evals_tries,
        n_evals_score=res.n_evals_score,
    )


def test_key_off_equals_key_absent_and_adds_nothing(monkeypatch: pytest.MonkeyPatch) -> None:
    outs = []
    for kw in ({}, {"repair": False}):
        bfs, mic = _repair_bfs(monkeypatch, [0.0, 2.0, 2.0, 2.0])
        res = BP.run_boundary_pass(bfs, mic, **kw)
        outs.append((_summary(res, mic), [np.asarray(v.orientation).tobytes() for v in mic.voxels]))
        assert res.repair is None and bfs.repair_log == []
        assert not any("repair" in k for k in bfs.stats)
    assert outs[0] == outs[1]


def test_key_on_leaves_the_pass_itself_unchanged(monkeypatch: pytest.MonkeyPatch) -> None:
    bfs0, mic0 = _repair_bfs(monkeypatch, [0.0, 2.0, 2.0, 2.0])
    r0 = BP.run_boundary_pass(bfs0, mic0)
    bfs1, mic1 = _repair_bfs(monkeypatch, [0.0, 2.0, 2.0, 2.0])
    r1 = BP.run_boundary_pass(bfs1, mic1, repair=True)
    assert _summary(r0, mic0) == _summary(r1, mic1)
    for k in ("boundary_voxels", "boundary_tries", "boundary_changed", "n_evals_boundary_pass"):
        assert bfs0.stats[k] == bfs1.stats[k]


def test_a_contested_voxel_is_repaired_and_every_evaluation_counted(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Voxel 0 starts at 0 degrees (truth 10). The pass moves it to 2 degrees (a neighbour's end;
    hit 0.55 against the neighbours' 0.95; its old end, 2 degrees away, scored 0.45: margin 0.10).
    It is contested; its own-end fit reaches the truth with 121 evaluations."""
    bfs, mic = _repair_bfs(monkeypatch, [0.0, 2.0, 2.0, 2.0])
    res = BP.run_boundary_pass(bfs, mic, repair=True)
    assert res.changed == [0]
    rep = res.repair
    assert rep.flagged == [0]
    assert [(f["voxel"], f["start"]) for f in rep.fits] == [(0, None)]  # no candidate: all at 2
    assert rep.changed == [0] and rep.n_evals_fits == 121 and rep.n_evals_score == 0
    assert _angle(mic.voxels[0].orientation, _rz(10.0)) < 1e-4
    assert mic.voxels[0].overlap_ratio == pytest.approx(0.95)
    s = bfs.stats
    assert (s["repair_flagged"], s["repair_fits"], s["repair_changed"]) == (1, 1, 1)
    assert s["n_evals_boundary_repair"] == 121
    assert len(bfs.repair_log) == 1 and bfs.repair_log[0]["voxel"] == 0
    assert _angle(bfs.repair_log[0]["old"], _rz(2.0)) < 1e-4  # the state after the pass
    # the state before the repair is kept, and the pass's own counts did not move
    assert _angle(rep.R_pre[0], _rz(2.0)) < 1e-4
    assert s["n_evals_boundary_pass"] == 3 * 61 and s["n_evals_boundary_score"] == 3


def test_each_fit_has_the_generator_of_its_voxel_and_start(monkeypatch: pytest.MonkeyPatch) -> None:
    bfs, mic = _repair_bfs(monkeypatch, [0.0, 2.0, 2.0, 2.0])
    BP.run_boundary_pass(bfs, mic, repair=True)
    assert RepairFake.states == [(0, BR.fit_rng(0, None).bit_generator.state)]


def test_nothing_contested_means_no_fit(monkeypatch: pytest.MonkeyPatch) -> None:
    bfs, mic = _repair_bfs(monkeypatch, [2.0, 2.0, 2.0, 2.0])  # a clean map: no boundary
    res = BP.run_boundary_pass(bfs, mic, repair=True)
    assert res.repair.flagged == [] and res.repair.fits == [] and res.repair.n_evals == 0
    assert bfs.stats["n_evals_boundary_repair"] == 0


def test_a_missing_cost_costs_one_counted_evaluation(monkeypatch: pytest.MonkeyPatch) -> None:
    """A contested voxel whose stored cost is NaN (its hit ratio finite) is scored once, counted.
    A NaN hit ratio never reaches the scoring: its ratio is NaN and ``NaN < 0.95`` is False, so
    such a voxel is not contested (the script rule behaves the same)."""
    bfs, mic = _repair_bfs(monkeypatch, [0.0, 2.0, 2.0, 2.0])
    res = BP.run_boundary_pass(bfs, mic)
    nbrs = [
        sorted(int(j) for j in mic.get_neighbors(v, radius=2.0) if int(j) != v) for v in range(4)
    ]
    calls: List[int] = []
    mic.voxels[0].cost = float("nan")  # voxel 0 (changed by the pass) keeps its finite hit ratio
    cf = bfs.reconstructor._make_local_cost_fn()
    orig_eval = cf.evaluate

    def counting(*a: Any, **k: Any) -> Any:
        calls.append(1)
        return orig_eval(*a, **k)

    cf.evaluate = counting  # type: ignore[method-assign]
    R_before = np.stack([_rz(d) for d in [0.0, 2.0, 2.0, 2.0]])
    rep = BR.run_repair(bfs, mic, res, R_before, nbrs, IDENT, None, cf)
    assert rep.flagged == [0]
    assert rep.n_evals_score == 1 == len(calls)
    assert rep.n_evals == rep.n_evals_fits + 1
    # a NaN hit ratio is never contested, so it costs nothing
    bfs2, mic2 = _repair_bfs(monkeypatch, [0.0, 2.0, 2.0, 2.0])
    res2 = BP.run_boundary_pass(bfs2, mic2)
    mic2.voxels[0].overlap_ratio = float("nan")
    rep2 = BR.run_repair(bfs2, mic2, res2, R_before, nbrs, IDENT, None, cf)
    assert rep2.flagged == [] and rep2.n_evals_score == 0


# ------------------------------------------------------------------------------ 4. the BFS call
def test_the_bfs_passes_the_key_only_when_it_is_on(monkeypatch: pytest.MonkeyPatch) -> None:
    calls: List[Dict[str, Any]] = []

    def spy(bfs: Any, mic: Any, **kw: Any) -> None:
        calls.append(kw)

    monkeypatch.setattr(BP, "run_boundary_pass", spy)
    _run((0.5, -0.5), monkeypatch, boundary_pass=True)
    _run((0.5, -0.5), monkeypatch, boundary_pass=True, boundary_pass_repair=True)
    _run(
        (0.5, -0.5), monkeypatch, boundary_pass=True, boundary_pass_repair=True,
        boundary_pass_skip_hit=0.9,
    )  # fmt: skip
    _run((0.5, -0.5), monkeypatch, boundary_pass=False)
    assert calls == [
        {},
        {"skip_hit": 0.0, "repair": True},
        {"skip_hit": 0.9, "repair": True},
    ]
