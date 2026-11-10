"""Task B tooling: event classification and the attribution of a change of the cost (synthetic
traces), the per-peak trace of a real ``VoxelCostFunction`` on the ThreeVoxels images
(recombination, restored module functions, attribution of a real change of the cost), and the
invariants of a stored ray walk when its cache is present."""

import os
import sys
from pathlib import Path
from types import SimpleNamespace
from typing import Any, Dict

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts" / "eval_reduction"))
import cost_edges as CE  # noqa: E402

DIMS = (2048, 2048)
EXAMPLE_DIR = ROOT.parent / "Examples" / "Example2.ThreeVoxels"
_CACHE: Dict[str, Any] = {}


def trace(P: int = 4) -> Dict[str, Any]:
    """P peaks, all evaluated on 2 detectors, one vertex-free (screen) signature each, q = 1."""
    return dict(
        status=np.full(P, CE.PRESENT, np.int8),
        frame=np.arange(P, dtype=np.int32),
        hit=np.ones((P, 2), np.uint8),
        sig=np.full((P, 2, 1, 2), 100, np.int32),
        cnt=np.tile(np.array([[5, 5, 1, 1, 2]], np.int32), (P, 1)),
        q=np.ones(P),
        rays=np.zeros((P, 8), np.int32),
        q_sum=float(P),
        n_q=P,
    )


def test_no_change_is_no_event() -> None:
    a, b = trace(), trace()
    ev = CE.classify(a, b, DIMS, 3)
    assert (ev == -1).all()
    d, parts, flags = CE.attribute(a, b, ev)
    assert d == 0 and parts == {} and flags["n_den_only"] == 0


def test_event_types_and_priority() -> None:
    a, b = trace(5), trace(5)
    b["frame"][0] += 1  # frame
    b["sig"][1, 0, 0, 0] += 1  # pixel
    b["sig"][2, 1, 0, 1] = 5000  # leaves the detector area: detector (wins over the pixel change)
    b["status"][3] = 1  # eta cut: peak not evaluated any more
    b["cnt"][4, 0] += 1  # counts only: subpixel
    ev = CE.classify(a, b, DIMS, 3)
    names = [CE.EVENT_TYPES[i] if i >= 0 else None for i in ev]
    assert names == ["frame", "pixel", "detector", "eta_cut", "subpixel"]


def test_attribution_sums_to_the_change_of_the_cost() -> None:
    a, b = trace(4), trace(4)
    # peak 0: frame crossing drops its quality to 0; peak 1: pixel event lowers it to 0.5
    b["frame"][0] += 1
    b["q"][0] = 0.0
    b["sig"][1, 0, 0, 0] += 1
    b["q"][1] = 0.5
    b["q_sum"], b["n_q"] = 2.5, 4
    ev = CE.classify(a, b, DIMS, 3)
    d, parts, _ = CE.attribute(a, b, ev)
    assert d == pytest.approx(1 - 2.5 / 4 - 0.0)  # cost 0 -> 0.375
    assert sum(parts.values()) == pytest.approx(d)
    assert parts["frame"] == pytest.approx(0.25) and parts["pixel"] == pytest.approx(0.125)


def test_denominator_change_with_no_numerator_change() -> None:
    a, b = trace(5), trace(5)
    a["status"][4] = 1  # the fifth peak is eta-cut in a ...
    a["q"][4] = 0.0
    a["q_sum"], a["n_q"] = 4.0, 4
    b["q_sum"], b["n_q"] = (
        5.0,
        5,
    )  # ... and evaluated in b with quality 1: the mean stays 1, N changes
    ev = CE.classify(a, b, DIMS, 3)
    d, parts, flags = CE.attribute(a, b, ev)
    assert d == pytest.approx(0.0)
    assert flags["dN"] == 1
    assert sum(parts.values()) == pytest.approx(0.0, abs=1e-12)
    assert set(parts) == {"eta_cut", "denominator"}


def test_ray_radii_windows() -> None:
    w = CE.ray_radii()
    assert set(w) == {"fin", "mc", "grid"}
    assert w["fin"].shape == (1, 101) and w["grid"].shape == (4, 101)
    assert w["grid"][-1, -1] == pytest.approx(5.0)


@pytest.fixture(scope="module")
def three_voxels() -> Dict[str, Any]:
    """ThreeVoxels forward-simulated images (sparse) and the base namespace ``CE.make_cf`` needs
    (same setup as the fixture of tests/test_center_screen.py, plus eta limit and registration)."""
    if _CACHE:
        return _CACHE
    from icenine.config_file import ConfigFile
    from icenine.experiment_setup import XDMExperimentSetup
    from icenine.experimental_data import ExperimentalData
    from icenine.mic_file import MicFile
    from icenine.reconstructor import _get_voxel_vertices
    from icenine.sample import Sample
    from icenine.simulation import Simulation

    config_path = EXAMPLE_DIR / "ConfigFiles" / "Example2.Simulation.config"
    data_dir = EXAMPLE_DIR / "ScatteringData_Python"
    mic_path = EXAMPLE_DIR / "SimInput" / "three_voxels.mic"
    if not config_path.exists() or not mic_path.exists():
        pytest.skip("ThreeVoxels config or .mic not found")
    if not data_dir.exists() or len(list(data_dir.glob("*.d*"))) < 360:
        pytest.skip(f"Forward sim output not found or incomplete: {data_dir}")
    old_cwd = os.getcwd()
    os.chdir(EXAMPLE_DIR)
    try:
        config = ConfigFile.from_file(str(config_path))
        exp_data = ExperimentalData.from_image_directory(
            directory=str(data_dir),
            basename="3Grains.sim",
            ext="d",
            serial_length=5,
            n_omega=180,
            n_detectors=2,
            num_rows=2048,
            num_cols=2048,
            mode="sparse",
        )
        setup = XDMExperimentSetup(config)
        setup.initialize_experiment()
        detector_list = setup.get_detector_list()
        sample = Sample()
        setup.initialize_sample(sample, detector_list[0])
        mic = MicFile.read(str(mic_path))
        voxel = mic.voxels[0]
        _CACHE.update(
            base=SimpleNamespace(
                simulator=Simulation(setup),
                detector_list=detector_list,
                range_map=setup.get_range_to_index_map(),
                exp_data=exp_data,
                sample=sample,
                structure_list=sample.get_structure_list(),
                eta_limit=setup.get_eta_limit(),
                reg_mm=float(config.detector_dist_deviation),
            ),
            R=np.asarray(voxel.orientation, dtype=np.float64),
            phase=int(voxel.phase),
            verts=_get_voxel_vertices(voxel),
        )
    finally:
        os.chdir(old_cwd)
    return _CACHE


def _module_functions() -> Any:
    from icenine import cost_functions as CF

    return (
        CF.get_scattering_omegas_torch,
        CF.scattering_omegas_fixed_beam,
        CF.calculate_diffraction_overlap_batched,
        CF._c_stage_d_overlap,
    )


@pytest.mark.parametrize("setting", ["full", "screen"])
@pytest.mark.parametrize("rule", ["old", "cpp", "r2"])
def test_trace_recombines_the_library_cost(
    three_voxels: Dict[str, Any], setting: str, rule: str
) -> None:
    """The per-peak trace recombines to the library's own quality and counts (float precision),
    and the wrapped module functions are the original objects afterwards."""
    c = three_voxels
    before = _module_functions()
    cf = CE.make_cf(c["base"], setting, rule)
    with CE.quiet():
        tr = CE.trace_eval(cf, c["R"], c["verts"], c["phase"])
    assert tr["agg"]["nq"] > 0 and tr["agg"]["pkov"] > 0
    assert CE.check_trace(tr, tol=1e-12) < 1e-12
    after = _module_functions()
    assert all(x is y for x, y in zip(before, after))


def test_trace_restores_module_functions_after_an_exception(three_voxels: Dict[str, Any]) -> None:
    c = three_voxels
    before = _module_functions()
    cf = CE.make_cf(c["base"], "full", "cpp")
    original = cf.evaluate

    def evaluate_then_raise(*a: Any, **k: Any) -> Any:
        original(*a, **k)  # runs through the wrapped module functions
        raise RuntimeError("raised inside evaluate")

    cf.evaluate = evaluate_then_raise
    with pytest.raises(RuntimeError, match="inside evaluate"), CE.quiet():
        CE.trace_eval(cf, c["R"], c["verts"], c["phase"])
    after = _module_functions()
    assert all(x is y for x, y in zip(before, after))


@pytest.mark.parametrize("setting", ["full", "screen"])
def test_parts_sum_to_a_real_change_of_the_cost(three_voxels: Dict[str, Any], setting: str) -> None:
    """classify + attribute on two nearby orientations: the parts sum to the change of the cost."""
    c = three_voxels
    cf = CE.make_cf(c["base"], setting, "cpp")
    d = CE.unit_directions(1, seed=7)[0]
    Rs = CE.rotate_about(c["R"], np.array([0.0, 0.3]), d)
    with CE.quiet():
        a, b = (CE.trace_eval(cf, R, c["verts"], c["phase"]) for R in Rs)
    ev = CE.classify(a, b, DIMS, int(cf.pixel_radius))
    dcost, parts, _ = CE.attribute(a, b, ev)
    expected = (1.0 - b["agg"]["quality"]) - (1.0 - a["agg"]["quality"])
    assert dcost == pytest.approx(expected, abs=1e-12)
    assert (ev >= 0).any() and dcost != 0.0
    assert abs(sum(parts.values()) - dcost) < 1e-12


STORED = ROOT / "benchmarks" / "eval_reduction" / "B" / "cache" / "full" / "v15901_d0.npz"


@pytest.mark.skipif(not STORED.exists(), reason="stored ray walk (gitignored cache) missing")
def test_stored_walk_parts_sum_to_cost_change() -> None:
    z = np.load(STORED)
    for combo in ("full|cpp", "screen|old"):
        c = z[f"{combo}|cost|dense"].astype(np.float64)
        parts = z[f"{combo}|parts|dense"].astype(np.float64)
        assert np.abs(parts.sum(1) - np.diff(c)).max() < 1e-4  # float32 storage
    assert float(z["maxdiff"]) < 1e-12
