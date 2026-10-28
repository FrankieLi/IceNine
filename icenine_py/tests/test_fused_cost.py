"""The fused evaluate of the opt-in key FusedCost (S2 Steps 1-2, docs/todo_s2_fused_c_eval.md).

1. Classification: every argument of ``VoxelCostFunction.__init__`` is classified in
   ``FUSED_ARGUMENTS``, and each "fallback" argument has a value that ``fused_unsupported_reason``
   rejects (a new cost mode cannot reach the fused path untested).
2. Key off: the C entry is never called and the results are those of a cost function built without
   the argument (the base-commit comparison is tests/test_s1_evaluate.py and s1_parity.py).
3. Key on: on ThreeVoxels voxel 0 every covered mode agrees with the torch evaluate except on
   orientations whose differences are all pixel-edge flips (parity rule D-20261027-I2: a traced
   ``pixel_integer`` row with the torch or the fused coordinate within 0.01 px of the integer).
4. Fallbacks (soft mode, the ray rule, a build without HAS_FUSED, a mode changed on the live
   object): the torch result, one RuntimeWarning per reason.
5. The lazy FrameTable fill and the context rebuild when the geometry, the data or a mode changes.
6. Wiring: the config key, ``SearchParameters`` and the reconstructor's cost flags.
"""

import contextlib
import inspect
import io
import os
import sys
import warnings
from pathlib import Path
from typing import Any, Dict, Iterator, List

import numpy as np
import pytest
import torch

from icenine import cost_functions as CF
from icenine.cost_functions import FUSED_ARGUMENTS, VoxelCostFunction, fused_unsupported_reason

S2 = Path(__file__).resolve().parents[1] / "scripts" / "speed_parity" / "s2"
SPEED = S2.parent
EXAMPLE = Path(__file__).resolve().parents[2] / "Examples" / "Example2.ThreeVoxels"
N = 24
BOUND_PX = 0.01

needs_fused = pytest.mark.skipif(
    not (CF._HAS_C_RASTERIZE and CF._C_HAS_FUSED),
    reason="needs the _rasterize extension built with the fused evaluate (HAS_FUSED)",
)

# covered modes (the names of scripts/speed_parity/s2/s2_parity.py)
MODES: Dict[str, Dict[str, Any]] = {
    "local": dict(pixel_radius=0, max_q=8.0),
    "cpp_rule": dict(pixel_radius=0, max_q=8.0, cpp_overlap_rule=True),
    "edge": dict(pixel_radius=0, max_q=8.0, cpp_edge_rules=True),
    "edge_cpp_rule": dict(pixel_radius=0, max_q=8.0, cpp_edge_rules=True, cpp_overlap_rule=True),
    "min_sin_eta": dict(pixel_radius=0, max_q=8.0, min_sin_eta=0.2),
    "screen_default": dict(pixel_radius=3, max_q=5.0, center_screen=True),
    "screen_center": dict(pixel_radius=3, max_q=5.0, cpp_overlap_rule=True, center_screen=True),
    "screen_vertex0": dict(pixel_radius=3, max_q=5.0, cpp_overlap_rule=True, center_screen=False),
    "screen_edge": dict(pixel_radius=3, max_q=5.0, cpp_overlap_rule=True, cpp_edge_rules=True),
}


@pytest.fixture(scope="module")
def s2() -> Iterator[Dict[str, Any]]:
    data = EXAMPLE / "ScatteringData_Python"
    if not data.exists() or len(list(data.glob("*.d*"))) < 360:
        pytest.skip("ThreeVoxels forward-simulated images not found")
    if not (CF._HAS_C_RASTERIZE and CF._C_HAS_FUSED):
        pytest.skip("needs the _rasterize extension built with the fused evaluate (HAS_FUSED)")
    names = ("s1_parity", "s2fused", "torch_mirror", "s2_parity")
    saved = {n: sys.modules.pop(n) for n in names if n in sys.modules}
    sys.path.insert(0, str(S2))
    sys.path.insert(0, str(SPEED))
    old = os.getcwd()
    threads = torch.get_num_threads()
    torch.set_num_threads(1)
    try:
        import s1_parity as S1
        import s2_parity as SP
        import torch_mirror as TM

        with contextlib.redirect_stdout(io.StringIO()):
            ds = S1.load_dataset("three0")
        R = S1.orientation_set(ds, N // 2, N // 2)
        yield dict(S1=S1, SP=SP, TM=TM, ds=ds, R=R)
    finally:
        torch.set_num_threads(threads)
        os.chdir(old)
        for p in (str(S2), str(SPEED)):
            sys.path.remove(p)
        for n in names:
            sys.modules.pop(n, None)
        sys.modules.update(saved)


def _tuple(oi: Any) -> tuple:
    return (
        oi.pixel_overlap,
        oi.pixel_on_detector,
        oi.peak_overlap,
        oi.peak_on_detector,
        oi.quality,
        oi.n_quality_points,
        oi.detectors_overlap,
    )


@pytest.fixture()
def fresh_warnings() -> Iterator[None]:
    saved = set(CF._WARNED_FUSED)
    CF._WARNED_FUSED.clear()
    try:
        yield
    finally:
        CF._WARNED_FUSED.clear()
        CF._WARNED_FUSED.update(saved)


# ---------------------------------------------------------------------------------------------
# 1. classification of the VoxelCostFunction arguments
# ---------------------------------------------------------------------------------------------


def test_every_init_argument_is_classified() -> None:
    params = set(inspect.signature(VoxelCostFunction.__init__).parameters) - {"self"}
    assert params == set(FUSED_ARGUMENTS), (
        "VoxelCostFunction.__init__ and FUSED_ARGUMENTS differ: classify the new argument and "
        f"make fused_unsupported_reason reject its values that C does not implement "
        f"(unclassified {sorted(params - set(FUSED_ARGUMENTS))}, "
        f"stale {sorted(set(FUSED_ARGUMENTS) - params)})"
    )
    assert set(FUSED_ARGUMENTS.values()) <= {"context", "covered", "fallback", "key"}


@needs_fused
def test_every_fallback_argument_has_a_rejected_value(s2: Dict[str, Any]) -> None:
    from icenine.ray_overlap import RayRule

    ds = s2["ds"]
    # a rejected value per "fallback" argument, set on the live object (Dataset.cost_fn fixes some
    # arguments, e.g. the setup's eta limit)
    rejected: Dict[str, Any] = {
        "mode": "soft",
        "ray_rule": RayRule(arm="R", radius_px=2, registration_mm=0.0),
        "eta_limit": None,
        "detector_list": list(ds.setup.detector_list) * 5,  # 10 detectors
    }
    for name, kind in FUSED_ARGUMENTS.items():
        if kind != "fallback":
            continue
        assert name in rejected, f"no rejected value tested for the fallback argument {name}"
        cf = ds.cost_fn(dict(MODES["local"], fused_cost=True))
        setattr(cf, name, rejected[name])
        assert fused_unsupported_reason(cf) is not None, name
    assert fused_unsupported_reason(ds.cost_fn(dict(MODES["local"], fused_cost=True))) is None


# ---------------------------------------------------------------------------------------------
# 2. key off
# ---------------------------------------------------------------------------------------------


@needs_fused
def test_key_off_never_calls_the_fused_entry(s2: Dict[str, Any], monkeypatch: Any) -> None:
    ds, R = s2["ds"], s2["R"]

    def boom(*a: Any) -> Any:
        raise AssertionError("fused entry called with the key off")

    monkeypatch.setattr(CF, "_c_fused_evaluate", boom)
    monkeypatch.setattr(CF, "_c_fused_context", boom)
    for mode in ("local", "screen_center"):
        a = ds.cost_fn(MODES[mode])
        b = ds.cost_fn(dict(MODES[mode], fused_cost=False))
        assert a.fused_cost is False and b.fused_cost is False
        for i in range(len(R)):
            assert _tuple(a.evaluate(R[i], ds.verts, ds.phase)) == _tuple(
                b.evaluate(R[i], ds.verts, ds.phase)
            )
        assert b._fused_ctx == {}


# ---------------------------------------------------------------------------------------------
# 3. key on: the parity rule
# ---------------------------------------------------------------------------------------------


@pytest.mark.parametrize("mode", list(MODES))
def test_key_on_differences_are_pixel_edge_flips(s2: Dict[str, Any], mode: str) -> None:
    ds, R, SP, TM = s2["ds"], s2["R"], s2["SP"], s2["TM"]
    kw = MODES[mode]
    lib, mir = ds.cost_fn(kw), ds.cost_fn(kw)
    fused = ds.cost_fn(dict(kw, fused_cost=True))
    n_det = len(ds.setup.detector_list)
    n_diff = 0
    with warnings.catch_warnings():
        warnings.simplefilter("error", RuntimeWarning)  # no fallback in a covered mode
        for i in range(len(R)):
            a = _tuple(lib.evaluate(R[i], ds.verts, ds.phase))
            b = _tuple(fused.evaluate(R[i], ds.verts, ds.phase))
            d = fused.fused_debug(R[i], ds.verts, ds.phase)
            assert tuple(d["result"]) == b  # the debug call computes what evaluate returns
            if a == b:
                continue
            n_diff += 1
            m = TM.mirror(mir, R[i], ds.verts, ds.phase)
            ref = TM.reference64(mir, R[i], ds.verts, ds.phase)
            rows = SP.trace(m, d, ref, mir, n_det)
            assert rows, f"orientation {i}: results differ but no stage-D input differs"
            for r in rows:
                assert r["cause"] == "pixel_integer", (i, r)
                assert min(r["torch_px_from_edge"], r["fused_px_from_edge"]) <= BOUND_PX, (i, r)
    assert fused.eval_count == len(R)  # fused_debug does not count
    assert n_diff <= max(2, len(R) // 10)


# ---------------------------------------------------------------------------------------------
# 4. fallbacks
# ---------------------------------------------------------------------------------------------


def _fallback_case(s2: Dict[str, Any], kw: Dict[str, Any], match: str, n: int = 4) -> None:
    ds, R = s2["ds"], s2["R"]
    torch_cf = ds.cost_fn(kw)
    cf = ds.cost_fn(dict(kw, fused_cost=True))
    with pytest.warns(RuntimeWarning, match=match) as rec:
        for i in range(n):
            assert _tuple(cf.evaluate(R[i], ds.verts, ds.phase)) == _tuple(
                torch_cf.evaluate(R[i], ds.verts, ds.phase)
            )
    assert sum("FusedCost" in str(w.message) for w in rec) == 1  # one warning per reason
    assert cf.eval_count == n


@needs_fused
def test_fallback_soft_mode(s2: Dict[str, Any], fresh_warnings: None) -> None:
    _fallback_case(s2, dict(MODES["local"], mode="soft"), "soft mode", n=2)


@needs_fused
def test_fallback_ray_rule(s2: Dict[str, Any], fresh_warnings: None) -> None:
    from icenine.ray_overlap import RayRule

    rule = RayRule(arm="R", radius_px=2, registration_mm=0.0)
    _fallback_case(s2, dict(MODES["screen_center"], ray_rule=rule), "ray rule")


@needs_fused
def test_fallback_build_without_fused(
    s2: Dict[str, Any], fresh_warnings: None, monkeypatch: Any
) -> None:
    monkeypatch.setattr(CF, "_C_HAS_FUSED", False)
    _fallback_case(s2, MODES["local"], "HAS_FUSED")


@needs_fused
def test_mode_changed_on_live_object(s2: Dict[str, Any], fresh_warnings: None) -> None:
    """A cost mode changed after the context was built rebuilds it; a change to a mode C does not
    cover falls back."""
    from icenine.ray_overlap import RayRule

    ds, R = s2["ds"], s2["R"]
    cf = ds.cost_fn(dict(MODES["local"], fused_cost=True))
    cf.evaluate(R[0], ds.verts, ds.phase)
    ctx0 = cf._fused_ctx[ds.phase]
    cf.cpp_edge_rules = True
    cf.cpp_overlap_rule = True
    ref = ds.cost_fn(dict(MODES["edge_cpp_rule"], fused_cost=True))
    for i in range(len(R)):
        assert _tuple(cf.evaluate(R[i], ds.verts, ds.phase)) == _tuple(
            ref.evaluate(R[i], ds.verts, ds.phase)
        )
    assert cf._fused_ctx[ds.phase] is not ctx0
    cf.ray_rule = RayRule(arm="R", radius_px=2, registration_mm=0.0)
    with pytest.warns(RuntimeWarning, match="ray rule"):
        cf.evaluate(R[0], ds.verts, ds.phase)


# ---------------------------------------------------------------------------------------------
# 5. the lazy frame fill and the context rebuild
# ---------------------------------------------------------------------------------------------


@needs_fused
def test_lazy_frame_fill(s2: Dict[str, Any]) -> None:
    ds, R = s2["ds"], s2["R"]
    kw = MODES["screen_center"]
    cf = ds.cost_fn(dict(kw, fused_cost=True))
    table = cf.frame_table()
    assert not table._filled.any()
    ctx = cf._fused_context(ds.phase)
    assert ctx is not None
    verts = cf._fused_vertices(ds.verts)
    first = CF._c_fused_evaluate(ctx.capsule, R[0], verts)
    assert isinstance(first, list) and first  # the frames to fill
    assert not table._filled.any()  # C only names them
    b = _tuple(cf.evaluate(R[0], ds.verts, ds.phase))
    assert table._filled[np.asarray(first)].all()
    assert b == _tuple(ds.cost_fn(dict(kw, fused_cost=True)).evaluate(R[0], ds.verts, ds.phase))
    assert cf.eval_count == 1


def _results(cf: Any, ds: Any, R: np.ndarray) -> List[tuple]:
    return [_tuple(cf.evaluate(R[i], ds.verts, ds.phase)) for i in range(len(R))]


@needs_fused
def test_context_rebuilt_when_the_geometry_changes(s2: Dict[str, Any]) -> None:
    """Moving a detector, moving the sample (in place) or replacing the data rebuilds the context:
    the results equal those of a cost function built after the change."""
    ds, R = s2["ds"], s2["R"]
    kw = dict(MODES["local"], fused_cost=True)
    det = ds.setup.detector_list[0]
    sample = ds.setup.sample
    pos0 = det._position.clone()
    s2l0 = sample.sample_to_lab_matrix.clone()
    try:
        cf = ds.cost_fn(kw)
        before = _results(cf, ds, R)
        ctx0 = cf._fused_ctx[ds.phase]
        det.set_position(pos0 + torch.tensor([0.0, 0.05, 0.0]))
        moved = _results(cf, ds, R)
        assert cf._fused_ctx[ds.phase] is not ctx0
        assert moved == _results(ds.cost_fn(kw), ds, R)
        assert moved != before  # the move is large enough to change some result
        ctx1 = cf._fused_ctx[ds.phase]
        sample.sample_to_lab_matrix[1, 3] += 0.01  # in place: the tensor version changes
        shifted = _results(cf, ds, R)
        assert cf._fused_ctx[ds.phase] is not ctx1
        assert shifted == _results(ds.cost_fn(kw), ds, R)
        ctx2 = cf._fused_ctx[ds.phase]
        assert _results(cf, ds, R[:2]) and cf._fused_ctx[ds.phase] is ctx2  # unchanged: kept
    finally:
        det.set_position(pos0)
        sample.sample_to_lab_matrix.copy_(s2l0)


@needs_fused
def test_context_rebuilt_when_the_data_is_replaced(s2: Dict[str, Any]) -> None:
    import copy

    ds, R = s2["ds"], s2["R"]
    cf = ds.cost_fn(dict(MODES["local"], fused_cost=True))
    cf.evaluate(R[0], ds.verts, ds.phase)
    ctx0, table0 = cf._fused_ctx[ds.phase], cf.frame_table()
    cf.exp_data = copy.copy(cf.exp_data)  # a new object with the same images
    r = _tuple(cf.evaluate(R[0], ds.verts, ds.phase))
    assert cf.frame_table() is not table0 and cf._fused_ctx[ds.phase] is not ctx0
    assert r == _tuple(
        ds.cost_fn(dict(MODES["local"], fused_cost=True)).evaluate(R[0], ds.verts, ds.phase)
    )


@needs_fused
def test_vertices_modified_in_place_are_reread(s2: Dict[str, Any]) -> None:
    ds, R = s2["ds"], s2["R"]
    kw = dict(MODES["local"], fused_cost=True)
    cf = ds.cost_fn(kw)
    verts = ds.verts.clone()
    a = [_tuple(cf.evaluate(R[i], verts, ds.phase)) for i in range(len(R))]
    with torch.no_grad():
        verts += 0.01
    b = [_tuple(cf.evaluate(R[i], verts, ds.phase)) for i in range(len(R))]
    c = [_tuple(ds.cost_fn(kw).evaluate(R[i], verts, ds.phase)) for i in range(len(R))]
    assert b == c and a != b


# ---------------------------------------------------------------------------------------------
# 6. wiring
# ---------------------------------------------------------------------------------------------


def test_config_key_default_off_and_parse(tmp_path: Any) -> None:
    from icenine.config_file import ConfigFile
    from icenine.orientation_search import SearchParameters

    base = Path(__file__).resolve().parents[1] / "scripts" / "phase_d" / "configs"
    text = (base / "ReconstructPhaseD_mc_clean.config").read_text()
    sp0 = SearchParameters.from_config(
        ConfigFile.from_file(str(base / "ReconstructPhaseD_mc_clean.config"))
    )
    assert sp0.fused_cost is False
    f = tmp_path / "on.config"
    f.write_text(text + "\nFusedCost 1\n")
    assert SearchParameters.from_config(ConfigFile.from_file(str(f))).fused_cost is True
    f.write_text(text + "\nFusedCost 2\n")
    with pytest.raises(ValueError, match="FusedCost"):
        ConfigFile.from_file(str(f))


@needs_fused
def test_reconstructor_passes_the_key(monkeypatch: Any) -> None:
    """Key on: every cost function the reconstructor builds through ``_cost_flags`` has it and the
    search runs through the fused entry; key off: the entry is never called."""
    from scipy.spatial.transform import Rotation
    from test_findoptimal_refactor import _build

    from icenine.reconstructor import _get_voxel_vertices

    calls = {"n": 0}
    real = CF._c_fused_evaluate

    def counting(*a: Any) -> Any:
        calls["n"] += 1
        return real(*a)

    monkeypatch.setattr(CF, "_c_fused_evaluate", counting)
    rec, voxel, _ = _build()
    assert "fused_cost" not in rec._cost_flags()
    rec.reconstruct_voxel(_get_voxel_vertices(voxel), voxel.phase, rng=np.random.default_rng(7))
    assert calls["n"] == 0
    rec, voxel, R_true = _build()
    rec.setup.search_params.fused_cost = True
    assert rec._cost_flags()["fused_cost"] is True
    assert rec._make_local_cost_fn().fused_cost is True
    res = rec.reconstruct_voxel(
        _get_voxel_vertices(voxel), voxel.phase, rng=np.random.default_rng(7)
    )
    assert calls["n"] > 0
    ang = np.degrees(Rotation.from_matrix(res.orientation @ R_true.T).magnitude())
    assert ang < 1.0 and np.isfinite(res.cost)  # a sane result (toy problem)
