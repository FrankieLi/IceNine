"""Smoke test of the S2 Step 0 fused-evaluate prototype (scripts/speed_parity/s2/, ThreeVoxels).

Pins the facts the Step 0 measurements rest on (docs/todo_s2_fused_c_eval.md, Step 0):
1. the torch mirror (``torch_mirror.mirror``) returns the library ``evaluate`` result, so its
   intermediates are the library's;
2. the fused prototype agrees with ``evaluate`` except on orientations whose differing stage-D
   input is a peak the trace finds (every differing orientation has at least one such peak);
3. the lazy frame fill: a fresh cost function's first fused call fills the frames it needs;
4. numpy's promotion of ``(omega - low) / width`` is classified as the C code assumes;
5. on macOS, the Accelerate build gives torch's g_lab bit for bit.

It also tests the pure helpers of the Task A check (``taska_check.py``): the C++ float32 product
order and the parity-rule classification.
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


# ------------------------------------------------------------------------------------------------
# S2 Task A check helpers (scripts/speed_parity/s2/taska_check.py): the C++ float32 product order
# and the parity-rule classification. Pure functions; no dataset needed.
# ------------------------------------------------------------------------------------------------
@pytest.fixture(scope="module")
def taska() -> Iterator[Any]:
    names = ("taska_check", "parity_common", "s1_parity", "s2fused", "torch_mirror", "s2_parity")
    saved = {n: sys.modules.pop(n) for n in names if n in sys.modules}
    paths = [str(S2), str(SPEED), str(SPEED.parent / "eval_reduction")]
    for p in paths:
        sys.path.insert(0, p)
    try:
        import taska_check as TC

        yield TC
    finally:
        for p in paths:
            sys.path.remove(p)
        for n in names:
            sys.modules.pop(n, None)
        sys.modules.update(saved)


# a = 1 + 2^-12: a*a = 1 + 2^-11 + 2^-24 rounds (tie to even) to 1 + 2^-11 in float32, so
# a*a - 1 is 2^-11 + 2^-24 when the product is fused into the add and 2^-11 when it is rounded first.
_A = np.float32(1 + 2**-12)
_FUSED = np.float32(2**-11 + 2**-24)
_ROUNDED = np.float32(2**-11)


def _dot(TC: Any, slots: Any, order: str) -> np.float32:
    """Element (0, 0) of ``product`` with row 0 of local and column 0 of centre from ``slots``."""
    L = np.zeros((3, 3), np.float32)
    C = np.zeros((3, 3), np.float32)
    for k, (x, y) in enumerate(slots):
        L[0, k], C[k, 0] = x, y
    return np.float32(TC.product(L, C, order)[0, 0])


def test_taska_fmaf_is_fused(taska: Any) -> None:
    assert taska.fmaf(_A, _A, np.float32(-1)) == _FUSED
    assert np.float32(np.float32(_A * _A) - np.float32(1)) == _ROUNDED


def test_taska_product_clang_order(taska: Any) -> None:
    # clang -ffp-contract=on on m[i][0]*R[0][j] + m[i][1]*R[1][j] + m[i][2]*R[2][j]:
    # fmaf(a2*b2, fmaf(a0*b0, a1*b1)); the a*a term is fused when it is in slot 0 or 2, not in 1.
    neg1, zero = (np.float32(-1), np.float32(1)), (np.float32(0), np.float32(0))
    cases = {
        "slot0": [(_A, _A), neg1, zero],
        "slot2": [neg1, zero, (_A, _A)],
        "slot1": [zero, (_A, _A), neg1],
    }
    expect = {
        "clang_contract": (_FUSED, _FUSED, _ROUNDED),
        "plain": (_ROUNDED, _ROUNDED, _ROUNDED),
        "fma_chain": (_ROUNDED, _FUSED, _ROUNDED),
        "fma_rev": (_FUSED, _ROUNDED, _FUSED),
    }
    assert taska.PRODUCT_ORDER == "clang_contract"
    for order, want in expect.items():
        got = tuple(_dot(taska, cases[c], order) for c in ("slot0", "slot2", "slot1"))
        assert got == want, order


def test_taska_product_identity_local(taska: Any) -> None:
    rng = np.random.default_rng(0)
    C = rng.standard_normal((3, 3)).astype(np.float32)
    for order in taska.ORDERS:
        assert np.array_equal(taska.product(np.eye(3), C, order), C), order
    m = dict(local=np.eye(3).tolist(), centre=C.astype(np.float64).tolist())
    assert np.array_equal(taska.orientation_of(m), C.astype(np.float64))


def test_taska_classify(taska: Any) -> None:
    pw = 0.0015  # detector pixel width (any positive value)

    def pint(t: float, f: float) -> Dict[str, Any]:
        return dict(cause="pixel_integer", torch_px_from_edge=t, fused_px_from_edge=f)

    v = taska.classify([pint(0.002, 0.004), pint(0.01, 0.03)], pw)
    assert v == dict(causes=["pixel_integer"], worst_near_px=0.01, within_rule=True)
    v = taska.classify([pint(0.002, 0.004), pint(0.02, 0.03)], pw)
    assert v["within_rule"] is False and v["worst_near_px"] == 0.02
    # an empty trace does not explain a difference
    assert taska.classify([], pw) == dict(causes=[], worst_near_px=0.0, within_rule=False)
    # the -1 rule: the distance is the detector coordinate (torch) over the pixel width
    neg = dict(cause="pixel_negative", torch=-5.25e-05 * pw)
    v = taska.classify([neg], pw)
    assert v["causes"] == ["pixel_negative"] and v["within_rule"] is True
    assert v["worst_near_px"] == pytest.approx(5.25e-05)
    v = taska.classify([dict(cause="pixel_negative", torch=-0.02 * pw)], pw)
    assert v["within_rule"] is False
    # any other cause (frame, eta, observability, t > 0) is outside the rule
    v = taska.classify([pint(0.001, 0.001), dict(cause="frame")], pw)
    assert v == dict(causes=["frame", "pixel_integer"], worst_near_px=0.001, within_rule=False)
