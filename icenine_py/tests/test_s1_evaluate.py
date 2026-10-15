"""S1 (docs/todo_s1_evaluate_loops.md): the per-call loops are out of ``VoxelCostFunction.evaluate``
and its results are bit-identical.

* ``frame_of_omegas`` against the Python loop it replaced (kept here as the reference), on synthetic
  ranges with gaps, edge omegas and both edge rules;
* ``FrameTable``: lazy fill, nothing touched at construction, size check, and rebuild after the
  ``exp_data`` of a live cost function is swapped;
* ``scattering_omegas_fixed_beam`` against ``get_scattering_omegas_torch`` on observable rows
  (random g, chi 0 and not 0, two energies), and the constants' rebuild when the beam changes;
* ``evaluate`` against ``tests/data/s1_evaluate_ref.npz`` (made on the base commit by
  ``scripts/speed_parity/s1_parity.py dump --fixture``): every OverlapInfo field of every cost mode,
  exactly. Single-thread torch, as the fixture; a fixture made on another CPU or torch build could
  differ in the last bit of a float32 geometry value.
"""

import os
import sys
from pathlib import Path
from types import SimpleNamespace
from typing import Any, Dict, Iterator, List

import numpy as np
import pytest
import torch

from icenine import cost_functions as CF
from icenine.cost_functions import FrameTable, VoxelCostFunction
from icenine.diffraction_core import (
    BraggConstants,
    get_scattering_omegas_torch,
    scattering_omegas_fixed_beam,
)
from icenine.simulation_range import OmegaRange, SimulationRange

from test_center_screen import three_voxels  # noqa: F401  (fixture)

ROOT = Path(__file__).resolve().parents[1]
FIXTURE = ROOT / "tests" / "data" / "s1_evaluate_ref.npz"
sys.path.insert(0, str(ROOT / "scripts" / "speed_parity"))
import s1_parity as S1  # noqa: E402

needs_c = pytest.mark.skipif(not CF._HAS_C_RASTERIZE, reason="C extension missing")


# ----------------------------------------------------------------------------------------------
# stage A: the frame lookup
# ----------------------------------------------------------------------------------------------
def _loop_reference(omegas: np.ndarray, rng: Any, edge: bool) -> np.ndarray:
    """The per-peak Python loop of calculate_diffraction_overlap_batched before S1."""
    bins = CF._omega_bins(omegas, rng.low, rng.width, edge)
    out = np.full(len(omegas), -1, dtype=np.int64)
    n_bins = len(rng.index_list)
    for i in range(len(omegas)):
        n = bins[i]
        if 0 <= n < n_bins and rng.index_list[n] is not None:
            out[i] = rng.index_list[n]
    return out


def _wedge_range(low: float, high: float, width: float, skip: Any) -> SimulationRange:
    """One wedge per bin (the center bin is the one that is marked), the bins in ``skip`` left
    out, so ``index_list`` has None entries and frame numbers differ from bin numbers."""
    n = round(abs(high - low) / width)
    sgn = 1.0 if high > low else -1.0
    wedges = [
        OmegaRange(low + sgn * k * width, low + sgn * (k + 1) * width)
        for k in range(n)
        if k not in skip
    ]
    return SimulationRange(low, high, width, wedges)


def _ranges() -> List[SimulationRange]:
    return [
        _wedge_range(-1.0, 1.0, 0.1, {3, 4, 11, 12, 13, 19}),  # gaps (None entries)
        _wedge_range(-1.0, 1.0, 0.1, set()),  # no gap
        _wedge_range(0.3, -0.7, 0.05, {2, 9, 10}),  # backwards range: negative width
        _wedge_range(-np.pi / 2, np.pi / 2, np.pi / 180, set(range(0, 180, 7))),
    ]


def _test_omegas(rng: SimulationRange) -> np.ndarray:
    lo, hi, w = rng.low, rng.high, abs(rng.width)
    r = np.random.default_rng(3)
    edges = np.array(
        [lo, hi, lo - w, hi + w, lo - 0.5 * w, lo - 1e-6, lo + 1e-6, hi - 1e-6, hi + 1e-6]
        + [lo + k * rng.width for k in range(len(rng.index_list) + 3)]
        + [lo - k * rng.width for k in range(1, 4)]
    )
    span = abs(hi - lo)
    rand = r.uniform(min(lo, hi) - 3 * w, max(lo, hi) + 3 * w, 3000)
    grid = np.linspace(min(lo, hi) - 2 * span, max(lo, hi) + 2 * span, 997)
    return np.concatenate([edges, rand, grid]).astype(np.float32)


@pytest.mark.parametrize("edge", [False, True], ids=["truncate", "cpp_edge_rules"])
@pytest.mark.parametrize("which", range(4))
def test_frame_of_omegas_equals_the_loop(which: int, edge: bool) -> None:
    rng = _ranges()[which]
    index = CF.frame_index_array(rng)
    om = _test_omegas(rng)
    assert index.dtype == np.int64 and len(index) == len(rng.index_list) + 1 and index[-1] == -1
    got = CF.frame_of_omegas(om, rng.low, rng.width, index, edge)
    want = _loop_reference(om, rng, edge)
    assert got.dtype == np.int64 and np.array_equal(got, want)
    # not vacuous: frames, gaps, below and above the range all occur
    assert (want >= 0).sum() > 100 and (want < 0).sum() > 100
    if which == 0:
        assert any(i is None for i in rng.index_list)
        bins = CF._omega_bins(om, rng.low, rng.width, edge)
        in_gap = [(0 <= b < len(index) - 1) and rng.index_list[b] is None for b in bins]
        assert any(in_gap)


def test_frame_of_omegas_edge_bins_follow_the_rule() -> None:
    """An omega in [low - width, low): bin 0 by truncation (historical), no frame with the floor."""
    rng = _ranges()[1]
    index = CF.frame_index_array(rng)
    om = np.array([rng.low - 0.5 * rng.width, rng.low + 0.5 * rng.width], dtype=np.float32)
    assert CF.frame_of_omegas(om, rng.low, rng.width, index, False).tolist() == [0, 0]
    assert CF.frame_of_omegas(om, rng.low, rng.width, index, True).tolist() == [-1, 0]


# ----------------------------------------------------------------------------------------------
# stage D: the frame table
# ----------------------------------------------------------------------------------------------
class _Exp:
    """A stand-in with the ExperimentalData interface the table uses; counts get_image calls."""

    def __init__(self, n_frames: int = 6, n_det: int = 2, rows: int = 8, cols: int = 10) -> None:
        self.calls: List[Any] = []
        self.shape = (rows, cols)
        self.frames = {
            (w, d): np.full((rows, cols), 1 + w * n_det + d, dtype=np.uint8)
            for w in range(n_frames)
            for d in range(n_det)
        }

    def get_image(self, w: int, d: int) -> Any:
        self.calls.append((w, d))
        arr = self.frames[(w, d)]
        return SimpleNamespace(
            num_rows=arr.shape[0], num_cols=arr.shape[1], get_binary_numpy=lambda: arr
        )


def test_frame_table_fills_lazily_and_once() -> None:
    exp = _Exp()
    t = FrameTable(exp, 6, 2)
    assert exp.calls == [] and len(t.images) == 12 and t.rows == 0  # nothing touched yet
    assert all(im is CF._PLACEHOLDER_FRAME for im in t.images)
    t.ensure(np.array([3, 1, 3], dtype=np.int32))
    assert sorted(exp.calls) == [(1, 0), (1, 1), (3, 0), (3, 1)]
    assert (t.rows, t.cols) == (8, 10)
    assert t.images[3 * 2 + 1] is exp.frames[(3, 1)] and t.images[1 * 2 + 0] is exp.frames[(1, 0)]
    assert t.images[0] is CF._PLACEHOLDER_FRAME  # frame 0 never used
    n = len(exp.calls)
    t.ensure(np.array([1, 3], dtype=np.int32))  # all filled: no get_image call
    assert len(exp.calls) == n
    t.ensure(np.array([1, 5], dtype=np.int32))  # only frame 5 is new
    assert sorted(exp.calls[n:]) == [(5, 0), (5, 1)]


def test_frame_table_is_built_without_touching_the_data() -> None:
    """A cost function may be built with exp_data=None and get its data attached later."""
    t = FrameTable(None, 4, 2)  # type: ignore[arg-type]
    assert t.source is None and len(t.images) == 8


def test_frame_table_checks_the_detector_0_frame_size() -> None:
    exp = _Exp()
    exp.frames[(4, 0)] = np.zeros((8, 11), dtype=np.uint8)  # a frame of another width
    t = FrameTable(exp, 6, 2)
    t.ensure(np.array([1], dtype=np.int32))
    with pytest.raises(ValueError, match="size"):
        t.ensure(np.array([4], dtype=np.int32))
    assert not t._filled[4]  # nothing marked as filled for the failed frame
    t.ensure(np.array([2], dtype=np.int32))  # the table is still usable


@needs_c
def test_table_passed_to_the_overlap_function_is_equivalent_to_none(
    three_voxels: Dict[str, Any], monkeypatch: pytest.MonkeyPatch  # noqa: F811
) -> None:
    """A caller that passes no index and no table gets the same OverlapInfo; a table of another
    exp_data is ignored."""
    c = three_voxels
    cf = VoxelCostFunction(pixel_radius=0, mode="hard", cpp_overlap_rule=True, **c["kwargs"])
    cap: List[Dict[str, Any]] = []
    orig = CF.calculate_diffraction_overlap_batched

    def spy(**kw: Any) -> Any:
        cap.append(kw)
        return orig(**kw)

    monkeypatch.setattr(CF, "calculate_diffraction_overlap_batched", spy)
    truth = np.asarray(c["voxel"].orientation, dtype=np.float64)
    full = cf.evaluate(truth, c["vertices"], c["voxel"].phase)
    kw = cap[0]
    assert kw["frame_table"] is cf.frame_table() and kw["frame_index"] is cf._frame_index
    bare = {k: v for k, v in kw.items() if k not in ("frame_index", "frame_table")}
    assert orig(**bare) == full
    foreign = FrameTable(_Exp(), 180, 2)  # another exp_data: not used, results unchanged
    assert orig(**bare, frame_index=kw["frame_index"], frame_table=foreign) == full
    assert not foreign._filled.any()
    assert full.quality > 0.5


@needs_c
def test_swapping_exp_data_of_a_live_cost_function_follows_the_new_data(
    three_voxels: Dict[str, Any],  # noqa: F811
) -> None:
    """optimizer_sweep / findoptimal_sweep attach new images to a live cost function; the table is
    rebuilt (identity check) and the results are those of a cost function built on the new data."""
    c = three_voxels
    kw = dict(c["kwargs"])
    cf = VoxelCostFunction(pixel_radius=0, mode="hard", **kw)
    truth = np.asarray(c["voxel"].orientation, dtype=np.float64)
    ph, vt = c["voxel"].phase, c["vertices"]
    before = cf.evaluate(truth, vt, ph)
    t_before = cf.frame_table()
    assert before.quality > 0.5 and t_before.source is kw["exp_data"] and t_before._filled.any()

    class Dark:  # every frame dark, one shared zero array (as optimizer_sweep's PixelSetData)
        n_omega_intervals, n_detectors = 180, 2
        zero = np.zeros((2048, 2048), dtype=np.uint8)

        def get_image(self, w: int, d: int) -> Any:
            return SimpleNamespace(num_rows=2048, num_cols=2048, get_binary_numpy=lambda: self.zero)

    dark = Dark()
    cf.exp_data = dark  # type: ignore[assignment]
    swapped = cf.evaluate(truth, vt, ph)
    fresh = VoxelCostFunction(pixel_radius=0, mode="hard", **{**kw, "exp_data": dark})
    assert swapped == fresh.evaluate(truth, vt, ph)
    assert swapped.quality == 0.0 and swapped.pixel_overlap == 0 and swapped.peak_on_detector > 0
    assert cf.frame_table() is not t_before and cf.frame_table().source is dark
    cf.exp_data = kw["exp_data"]  # and back
    assert cf.evaluate(truth, vt, ph) == before
    # attach later: exp_data=None at construction, data assigned before the first evaluate
    late = VoxelCostFunction(pixel_radius=0, mode="hard", **{**kw, "exp_data": None})
    late.exp_data = kw["exp_data"]
    assert late.evaluate(truth, vt, ph) == before


# ----------------------------------------------------------------------------------------------
# Bragg: hoisted constants
# ----------------------------------------------------------------------------------------------
@pytest.mark.parametrize("chi", [0.0, 0.0173, -0.31], ids=["chi0", "chi_small", "chi_neg"])
@pytest.mark.parametrize("energy", [64.351, 50.02], ids=["64keV", "50keV"])
def test_scattering_omegas_fixed_beam_equals_the_reference_on_observable_rows(
    chi: float, energy: float
) -> None:
    r = np.random.default_rng(11)
    K = 112
    mag = torch.tensor(r.uniform(0.5, 8.0, K), dtype=torch.float32)
    n_obs = 0
    consts = BraggConstants(mag, energy, chi)
    for _ in range(60):
        u = r.normal(size=(K, 3))
        u[::3] = [0.0, 0.0, 1.0] + 0.04 * u[::3]  # near the rotation axis: not observable
        u /= np.linalg.norm(u, axis=1, keepdims=True)
        g = torch.from_numpy(u).float() * mag[:, None]
        g = g.T.contiguous().T  # a (K, 3) view of a (3, K) array, as evaluate's (R @ g_hkl.T).T
        ref = get_scattering_omegas_torch(g, mag, energy, chi)
        got = scattering_omegas_fixed_beam(g, consts)
        assert torch.equal(got.observable, ref.observable)
        obs = ref.observable
        n_obs += int(obs.sum())
        assert torch.equal(got.omega1[obs], ref.omega1[obs])
        assert torch.equal(got.omega2[obs], ref.omega2[obs])
    assert 1000 < n_obs < 60 * K - 1000  # both observable and unobservable rows occur


def test_bragg_constants_know_what_they_were_built_from() -> None:
    mag = torch.ones(5)
    c = BraggConstants(mag, 64.351, 0.0)
    assert c.is_for(mag, 64.351, 0.0)
    assert not c.is_for(mag.clone(), 64.351, 0.0)  # another |g| tensor
    assert not c.is_for(mag, 64.4, 0.0) and not c.is_for(mag, 64.351, 0.001)


@needs_c
def test_evaluate_rebuilds_the_bragg_constants_when_the_beam_changes(
    three_voxels: Dict[str, Any],  # noqa: F811
) -> None:
    c = three_voxels
    kw = c["kwargs"]
    sim = kw["simulator"]
    cf = VoxelCostFunction(pixel_radius=0, mode="hard", **kw)
    truth = np.asarray(c["voxel"].orientation, dtype=np.float64)
    ph, vt = c["voxel"].phase, c["vertices"]
    e0, chi0 = sim.beam_energy, sim.beam_deflection_chi
    try:
        base = cf.evaluate(truth, vt, ph)
        for energy, chi in ((e0 * 1.01, chi0), (e0, chi0 + 0.01)):
            sim.beam_energy, sim.beam_deflection_chi = energy, chi
            fresh = VoxelCostFunction(pixel_radius=0, mode="hard", **kw)
            assert cf.evaluate(truth, vt, ph) == fresh.evaluate(truth, vt, ph)
            assert cf.evaluate(truth, vt, ph) != base  # the change is visible in the cost
    finally:
        sim.beam_energy, sim.beam_deflection_chi = e0, chi0
    assert cf.evaluate(truth, vt, ph) == base


# ----------------------------------------------------------------------------------------------
# evaluate against the fixture of the base commit
# ----------------------------------------------------------------------------------------------
@pytest.fixture(scope="module")
def single_thread() -> Iterator[None]:
    prev = torch.get_num_threads()
    torch.set_num_threads(1)
    yield
    torch.set_num_threads(prev)


@pytest.fixture(scope="module")
def fixture_dataset(single_thread: None) -> Any:
    if not FIXTURE.exists():
        pytest.skip("tests/data/s1_evaluate_ref.npz missing")
    data_dir = S1.THREE / "ScatteringData_Python"
    if not data_dir.exists() or len(list(data_dir.glob("*.d*"))) < 360:
        pytest.skip(f"Forward sim output not found or incomplete: {data_dir}")
    old = os.getcwd()
    try:
        return S1.load_dataset("three0", "sparse")  # chdir to the example directory
    finally:
        os.chdir(old)


def test_fixture_covers_the_modes_and_is_small() -> None:
    ref = np.load(FIXTURE, allow_pickle=False)
    modes = sorted({k.split("/")[1] for k in ref.files if k.count("/") == 2})
    assert modes == sorted(S1.MODES)
    assert FIXTURE.stat().st_size < 200 * 1024
    assert len(ref["three0/R"]) == 64
    # the edge rules change counts on this orientation set, so the edge modes are not vacuous
    ndiff = (ref["three0/local/peak_on_detector"] != ref["three0/edge/peak_on_detector"]).sum()
    assert ndiff >= 10


@needs_c
@pytest.mark.parametrize("mode", sorted(S1.MODES))
def test_evaluate_equals_the_base_commit_fixture_exactly(fixture_dataset: Any, mode: str) -> None:
    ref = np.load(FIXTURE, allow_pickle=False)
    R = ref["three0/R"][ref[f"three0/{mode}/index"]]
    got = S1.run_mode(fixture_dataset, mode, R)
    for f, v in got.items():
        want = ref[f"three0/{mode}/{f}"]
        assert v.dtype == want.dtype and np.array_equal(v, want), (mode, f)
