"""Smoke test of the S2 Step 0 fused-evaluate prototype (scripts/speed_parity/s2/, ThreeVoxels).

Pins the facts the Step 0 measurements rest on (docs/todo_s2_fused_c_eval.md, Step 0):
1. the torch mirror (``torch_mirror.mirror``) returns the library ``evaluate`` result, so its
   intermediates are the library's;
2. the fused prototype agrees with ``evaluate`` except on orientations whose differing stage-D
   input is a peak the trace finds (every differing orientation has at least one such peak);
3. the lazy frame fill: a fresh cost function's first fused call fills the frames it needs;
4. numpy's promotion of ``(omega - low) / width`` is classified as the C code assumes;
5. on macOS, the Accelerate build gives torch's g_lab bit for bit.
"""

import contextlib
import io
import os
import shutil
import sys
from pathlib import Path
from typing import Any, Dict, Iterator

import numpy as np
import pytest
import torch

S2 = Path(__file__).resolve().parents[1] / "scripts" / "speed_parity" / "s2"
SPEED = S2.parent
EXAMPLE = Path(__file__).resolve().parents[2] / "Examples" / "Example2.ThreeVoxels"
N = 24


@pytest.fixture(scope="module")
def s2() -> Iterator[Dict[str, Any]]:
    data = EXAMPLE / "ScatteringData_Python"
    if not data.exists() or len(list(data.glob("*.d*"))) < 360:
        pytest.skip("ThreeVoxels forward-simulated images not found")
    if shutil.which("cc") is None:
        pytest.skip("needs a C compiler (cc) to build the scratch extension")
    from icenine import cost_functions as CF

    if not CF._HAS_C_RASTERIZE:
        pytest.skip("needs the _rasterize C extension")
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
        import s2fused as SF
        import torch_mirror as TM

        with contextlib.redirect_stdout(io.StringIO()):
            ds = S1.load_dataset("three0")
        R = S1.orientation_set(ds, N // 2, N // 2)
        yield dict(S1=S1, SF=SF, TM=TM, SP=SP, ds=ds, R=R)
    finally:
        torch.set_num_threads(threads)
        os.chdir(old)
        for p in (str(S2), str(SPEED)):
            sys.path.remove(p)
        for n in names:
            sys.modules.pop(n, None)
        sys.modules.update(saved)


MODES = ("local", "cpp_rule", "edge_cpp_rule", "min_sin_eta", "screen_center", "screen_vertex0")


@pytest.mark.parametrize("mode", MODES)
def test_mirror_equals_library(s2: Dict[str, Any], mode: str) -> None:
    ds, R, SF, TM, SP = s2["ds"], s2["R"], s2["SF"], s2["TM"], s2["SP"]
    kw = SP.MODES[mode]
    lib, mir = ds.cost_fn(kw), ds.cost_fn(kw)
    for i in range(len(R)):
        a = SF.result_tuple(lib.evaluate(R[i], ds.verts, ds.phase))
        assert tuple(TM.mirror(mir, R[i], ds.verts, ds.phase)["result"]) == a, i


@pytest.mark.parametrize("mode", MODES)
def test_fused_differences_are_traced_peaks(s2: Dict[str, Any], mode: str) -> None:
    ds, R, SF, TM, SP = s2["ds"], s2["R"], s2["SF"], s2["TM"], s2["SP"]
    kw = SP.MODES[mode]
    lib, mir = ds.cost_fn(kw), ds.cost_fn(kw)
    fe = SF.FusedEval(ds.cost_fn(kw))
    n_det = len(ds.setup.detector_list)
    n_diff = 0
    for i in range(len(R)):
        a = SF.result_tuple(lib.evaluate(R[i], ds.verts, ds.phase))
        b = SF.result_tuple(fe.evaluate(R[i], ds.verts, ds.phase))
        d = fe.debug(R[i], ds.verts, ds.phase)
        assert tuple(d["result"]) == b  # the debug call computes what the plain call returns
        if a != b:
            n_diff += 1
            m = TM.mirror(mir, R[i], ds.verts, ds.phase)
            ref = TM.reference64(mir, R[i], ds.verts, ds.phase)
            rows = SP.trace(m, d, ref, mir, n_det)
            assert rows, f"orientation {i}: results differ but no stage-D input differs"
    assert n_diff <= max(2, len(R) // 10)


def test_reference64_edge_rule(s2: Dict[str, Any]) -> None:
    """With CppEdgeRules the float64 reference applies the -1 rule as the library does, so its
    pixel error against torch stays at float32 size (without the rule it is about 1 px or more).
    ThreeVoxels has no negative detector coordinate, so detector 0's origin is moved by the median
    column coordinate of the first orientation for the test (restored afterwards)."""
    ds, R, TM, SP = s2["ds"], s2["R"], s2["TM"], s2["SP"]
    mir = ds.cost_fn(SP.MODES["edge_cpp_rule"])
    det = mir.detector_list[0]
    m0 = TM.mirror(mir, R[0], ds.verts, ds.phase)
    shift = float(np.nanmedian(m0["jk"][m0["frame"] >= 0, 0, :, 0]))
    origin = det._lab_frame_coord_origin
    det._lab_frame_coord_origin = origin + shift * det._lab_frame_basis_j
    try:
        total, n, neg = 0.0, 0, 0
        for i in range(len(R)):
            m = TM.mirror(mir, R[i], ds.verts, ds.phase)
            ref = TM.reference64(mir, R[i], ds.verts, ds.phase)
            rows = m["frame"] >= 0
            e = np.abs(m["pix"][rows].astype(np.float64) - ref["pix"][rows])
            total, n = total + float(e.sum()), n + e.size
            neg += int((ref["pix"][rows] == -1.0).sum())
    finally:
        det._lab_frame_coord_origin = origin
    assert neg > 0  # the rule is exercised
    assert total / n < 0.01


def test_lazy_frame_fill(s2: Dict[str, Any]) -> None:
    ds, R, SF, SP = s2["ds"], s2["R"], s2["SF"], s2["SP"]
    kw = SP.MODES["screen_center"]
    cf = ds.cost_fn(kw)
    fe = SF.FusedEval(cf)
    table = cf.frame_table()
    assert not table._filled.any()
    ctx = fe.context(ds.phase)
    first = fe._evaluate(ctx.capsule, R[0], fe._verts(ds.verts))
    assert isinstance(first, list) and first  # the frames to fill
    b = SF.result_tuple(fe.evaluate(R[0], ds.verts, ds.phase))
    assert table._filled[np.asarray(first)].all()
    a = SF.result_tuple(ds.cost_fn(kw).evaluate(R[0], ds.verts, ds.phase))
    assert a == b
    assert cf.eval_count == 1


def test_frame_mode(s2: Dict[str, Any]) -> None:
    # The expected modes follow NumPy 2 (NEP 50) promotion of a Python float with a numpy
    # scalar, as under the lock (numpy 2.0.2); NumPy 1 value-based casting would differ.
    SF = s2["SF"]
    assert SF.frame_mode(np.float64(-1.5), 0.01) == 0
    assert SF.frame_mode(-1.5, np.float64(0.01)) == 1
    assert SF.frame_mode(-1.5, 0.01) == 2


@pytest.mark.skipif(sys.platform != "darwin", reason="Accelerate is macOS only")
def test_accelerate_glab_bit_identical(s2: Dict[str, Any]) -> None:
    ds, R, SF, TM, SP = s2["ds"], s2["R"], s2["SF"], s2["TM"], s2["SP"]
    kw = SP.MODES["local"]
    mir = ds.cost_fn(kw)
    fe = SF.FusedEval(ds.cost_fn(kw), "accelerate")
    for i in range(len(R)):
        m = TM.mirror(mir, R[i], ds.verts, ds.phase)
        d = fe.debug(R[i], ds.verts, ds.phase)
        assert np.array_equal(m["g_lab"].view(np.int32), d["g_lab"].view(np.int32)), i
