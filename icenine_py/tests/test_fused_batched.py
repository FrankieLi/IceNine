"""The batched fused entry (S2 Step 3, docs/todo_s2_fused_c_eval.md): ``BatchedCostEvaluator`` with
the key FusedCost on scores a batch in one C call (``_rasterize.fused_evaluate_many``).

1. Key on: every orientation of a batch equals the single-orientation fused ``evaluate`` (every
   ``OverlapInfo`` field, ``==``), for every mode the batched cost supports, at batch sizes 1, 7
   and 64, including the first calls that fill the FrameTable.
2. Key off: the fused entry is never called and the results are the torch ``evaluate``'s.
3. Dispatch: ``variant`` forces a variant; on a GPU device at or above ``GPU_MIN_BATCH`` the fixed
   variant keeps priority over the fused one.
4. Fallbacks (a build without ``fused_evaluate_many`` or without the fused evaluate): the compact
   variant, equal to the torch ``evaluate``, with one RuntimeWarning per reason.
"""

import warnings
from typing import Any, Dict, Iterator, List

import numpy as np
import pytest

from icenine import batched_cost as BC
from icenine import cost_functions as CF
from icenine.cost_functions import OverlapInfo

from test_batched_cost import FIELDS, _local, differing, orientations, run_chunks
from test_center_screen import three_voxels  # noqa: F401  (fixture)

needs_fused_many = pytest.mark.skipif(
    not (CF._HAS_C_RASTERIZE and CF._C_HAS_FUSED and BC._C_HAS_FUSED_MANY),
    reason="needs the _rasterize extension built with the batched fused entry",
)

FLAGS = {
    "default": {},
    "cpp_rule": {"cpp_overlap_rule": True},
    "min_sin_eta": {"min_sin_eta": 0.4},
    "cpp_edge_rules": {"cpp_edge_rules": True},
    "cpp_edge_rules+cpp_rule": {"cpp_edge_rules": True, "cpp_overlap_rule": True},
}


@pytest.fixture()
def fresh_warnings() -> Iterator[None]:
    saved = set(CF._WARNED_FUSED)
    CF._WARNED_FUSED.clear()
    try:
        yield
    finally:
        CF._WARNED_FUSED.clear()
        CF._WARNED_FUSED.update(saved)


def _single(c: Dict[str, Any], R: np.ndarray, **kw: Any) -> List[OverlapInfo]:
    """The single-orientation fused evaluate of every orientation (a fresh cost function)."""
    cf = _local(c, fused_cost=True, **kw)
    out = [cf.evaluate(r, c["vertices"], c["voxel"].phase) for r in R]
    assert cf._fused_ctx and all(x.reason is None for x in cf._fused_ctx.values())
    return out


@needs_fused_many
@pytest.mark.parametrize("flags", list(FLAGS.values()), ids=list(FLAGS))
@pytest.mark.parametrize("max_q", [0.0, 5.0], ids=["all_refl", "qmax5"])
@pytest.mark.parametrize("B", [1, 7, 64])
def test_batched_fused_equals_single_fused(
    three_voxels: Dict[str, Any], flags: Dict[str, Any], max_q: float, B: int  # noqa: F811
) -> None:
    c = three_voxels
    phase = c["voxel"].phase
    R = orientations(np.asarray(c["voxel"].orientation, dtype=np.float64))
    ref = _single(c, R, max_q=max_q, **flags)
    cf = _local(c, fused_cost=True, max_q=max_q, **flags)
    ev = BC.BatchedCostEvaluator(cf, "cpu")
    got = run_chunks(ev, R, c["vertices"], phase, B)  # variant None: fused with the key on
    assert differing(got, ref) == {f: 0 for f in FIELDS}
    assert got == ref
    assert cf.eval_count == len(R)
    assert ev.counts["orient_fused"] == len(R) and ev.counts["calls_fused"] == -(-len(R) // B)
    assert ev.counts["orient_compact"] == 0 and ev.counts["orient_fixed"] == 0
    assert ev.timing["fused"] > 0.0 and ev.timing["torch"] == 0.0
    # not vacuous: the truth is bright, the far orientations are not, and peaks exist throughout
    assert ref[0].quality > 0.5 and min(i.quality for i in ref) < 0.1
    assert sum(i.peak_on_detector > 0 for i in ref) > len(R) // 2


@needs_fused_many
def test_batched_fused_fills_the_frame_table_once(
    three_voxels: Dict[str, Any], monkeypatch: pytest.MonkeyPatch  # noqa: F811
) -> None:
    """A fresh cost function: the first call names the frames of the whole batch (each once), the
    FrameTable fills them, the call is repeated; later calls need no fill."""
    c = three_voxels
    phase = c["voxel"].phase
    R = orientations(np.asarray(c["voxel"].orientation, dtype=np.float64))[:40]
    cf = _local(c, fused_cost=True)
    ev = BC.BatchedCostEvaluator(cf, "cpu")
    calls: List[Any] = []
    real = BC._c_fused_evaluate_many

    def spy(*a: Any) -> Any:
        r = real(*a)
        calls.append(r[1])
        return r

    monkeypatch.setattr(BC, "_c_fused_evaluate_many", spy)
    first = ev.evaluate_many(R, c["vertices"], phase)
    assert len(calls) == 2 and calls[1] is None
    missing = calls[0]
    assert missing and len(missing) == len(set(missing))  # each frame once
    assert bool(cf.frame_table()._filled[missing].all())
    calls.clear()
    assert ev.evaluate_many(R, c["vertices"], phase) == first
    assert calls == [None]
    assert first == _single(c, R)


@needs_fused_many
def test_key_off_never_calls_the_fused_entry(
    three_voxels: Dict[str, Any], monkeypatch: pytest.MonkeyPatch  # noqa: F811
) -> None:
    c = three_voxels
    phase = c["voxel"].phase

    def boom(*a: Any) -> Any:
        raise AssertionError("the fused entry was called with the key off")

    monkeypatch.setattr(BC, "_c_fused_evaluate_many", boom)
    monkeypatch.setattr(CF, "_c_fused_evaluate", boom)
    R = orientations(np.asarray(c["voxel"].orientation, dtype=np.float64))
    cf, ref_cf = _local(c), _local(c)
    ev = BC.BatchedCostEvaluator(cf, "cpu")
    got = run_chunks(ev, R, c["vertices"], phase, 64)
    assert got == [ref_cf.evaluate(r, c["vertices"], phase) for r in R]
    assert ev.counts["orient_compact"] == len(R) and ev.counts["orient_fused"] == 0
    assert not cf._fused_ctx


@needs_fused_many
def test_variant_forces_a_variant(three_voxels: Dict[str, Any]) -> None:  # noqa: F811
    c = three_voxels
    phase = c["voxel"].phase
    R = orientations(np.asarray(c["voxel"].orientation, dtype=np.float64))[:30]
    single = _single(c, R)
    torch_ref = [_local(c).evaluate(r, c["vertices"], phase) for r in R]
    on = _local(c, fused_cost=True)
    ev = BC.BatchedCostEvaluator(on, "cpu")
    assert ev.evaluate_many(R, c["vertices"], phase, "compact") == torch_ref
    assert ev.evaluate_many(R, c["vertices"], phase, "fixed") == torch_ref
    assert ev.evaluate_many(R, c["vertices"], phase) == single
    off = _local(c)
    ev_off = BC.BatchedCostEvaluator(off, "cpu")
    assert ev_off.evaluate_many(R, c["vertices"], phase, "fused") == single
    assert ev_off.counts["orient_fused"] == len(R)
    with pytest.raises(ValueError, match="'fused'"):
        ev.evaluate_many(R, c["vertices"], phase, "gpu")


class _FixedCalled(Exception):
    pass


@needs_fused_many
def test_gpu_fixed_variant_keeps_priority(
    three_voxels: Dict[str, Any], monkeypatch: pytest.MonkeyPatch  # noqa: F811
) -> None:
    """On ``mps`` or ``cuda`` a batch of at least ``GPU_MIN_BATCH`` goes to the fixed variant even
    with the key on; a smaller batch to the fused one. (The device is set after construction and
    the fixed variant's first step, the GPU constants, is stubbed to raise, so the dispatch is
    checked without a GPU.)"""
    c = three_voxels
    phase = c["voxel"].phase
    R = orientations(np.asarray(c["voxel"].orientation, dtype=np.float64), n_each=20)
    assert len(R) >= BC.GPU_MIN_BATCH
    cf = _local(c, fused_cost=True)
    ev = BC.BatchedCostEvaluator(cf, "cpu")
    ev.device = "cuda"
    real_constants = ev._constants

    def constants(dev_name: str) -> Any:  # the fixed variant's first step on the GPU device
        if dev_name != "cpu":
            raise _FixedCalled
        return real_constants(dev_name)

    monkeypatch.setattr(ev, "_constants", constants)
    with pytest.raises(_FixedCalled):
        ev.evaluate_many(R[: BC.GPU_MIN_BATCH], c["vertices"], phase)
    small = R[: BC.GPU_MIN_BATCH - 1]
    assert ev.evaluate_many(small, c["vertices"], phase) == _single(c, small)
    assert ev.counts["orient_fused"] == len(small)


@needs_fused_many
@pytest.mark.skipif(not BC.device_available("mps"), reason="MPS not available")
def test_mps_fixed_variant_keeps_priority(three_voxels: Dict[str, Any]) -> None:  # noqa: F811
    c = three_voxels
    phase = c["voxel"].phase
    R = orientations(np.asarray(c["voxel"].orientation, dtype=np.float64), n_each=20)
    cf = _local(c, fused_cost=True)
    ev = BC.BatchedCostEvaluator(cf, "mps")
    ev.evaluate_many(R[: BC.GPU_MIN_BATCH], c["vertices"], phase)
    assert ev.counts["calls_fixed"] == 1 and ev.counts["calls_fused"] == 0
    small = R[: BC.GPU_MIN_BATCH - 1]
    assert ev.evaluate_many(small, c["vertices"], phase) == _single(c, small)
    assert ev.counts["calls_fused"] == 1


@needs_fused_many
def test_build_without_the_batched_entry_falls_back(
    three_voxels: Dict[str, Any],
    monkeypatch: pytest.MonkeyPatch,
    fresh_warnings: None,  # noqa: F811
) -> None:
    c = three_voxels
    phase = c["voxel"].phase
    R = orientations(np.asarray(c["voxel"].orientation, dtype=np.float64))[:30]
    torch_ref = [_local(c).evaluate(r, c["vertices"], phase) for r in R]
    monkeypatch.setattr(BC, "_C_HAS_FUSED_MANY", False)
    cf = _local(c, fused_cost=True)
    ev = BC.BatchedCostEvaluator(cf, "cpu")
    with pytest.warns(RuntimeWarning, match="fused_evaluate_many"):
        assert ev.evaluate_many(R, c["vertices"], phase) == torch_ref
    with warnings.catch_warnings():  # one warning per reason
        warnings.simplefilter("error")
        assert ev.evaluate_many(R, c["vertices"], phase) == torch_ref
    assert ev.counts["orient_compact"] == 2 * len(R) and ev.counts["orient_fused"] == 0
    assert cf.eval_count == 2 * len(R)
    # a build without the fused evaluate at all: the reason the single path gives
    monkeypatch.setattr(CF, "_C_HAS_FUSED", False)
    with pytest.warns(RuntimeWarning, match="HAS_FUSED"):
        got = BC.BatchedCostEvaluator(_local(c, fused_cost=True)).evaluate_many(
            R, c["vertices"], phase
        )
    assert got == torch_ref
    monkeypatch.setattr(BC, "_C_HAS_FUSED_MANY", True)  # the batched entry, but HAS_FUSED off
    CF._WARNED_FUSED.clear()
    with pytest.warns(RuntimeWarning, match="HAS_FUSED"):
        got = BC.BatchedCostEvaluator(_local(c, fused_cost=True)).evaluate_many(
            R, c["vertices"], phase
        )
    assert got == torch_ref


@needs_fused_many
def test_bad_orientation_shape_raises(three_voxels: Dict[str, Any]) -> None:  # noqa: F811
    c = three_voxels
    cf = _local(c, fused_cost=True)
    ctx = cf._fused_context(c["voxel"].phase)
    verts = cf._fused_vertices(c["vertices"])
    with pytest.raises(ValueError, match="B, 3, 3"):
        BC._c_fused_evaluate_many(ctx.capsule, np.eye(3), verts)
    with pytest.raises(TypeError, match="verts"):
        BC._c_fused_evaluate_many(ctx.capsule, np.eye(3)[None], verts.astype(np.float64))
    res, missing = BC._c_fused_evaluate_many(ctx.capsule, np.zeros((0, 3, 3)), verts)
    assert res == [] and missing is None


@needs_fused_many
def test_float32_orientations_match_single(three_voxels: Dict[str, Any]) -> None:  # noqa: F811
    """Orientations given as float32 go through float64 to float32 unchanged, as in the single
    path."""
    c = three_voxels
    phase = c["voxel"].phase
    R32 = orientations(np.asarray(c["voxel"].orientation, dtype=np.float64))[:20].astype(np.float32)
    cf = _local(c, fused_cost=True)
    ev = BC.BatchedCostEvaluator(cf, "cpu")
    assert ev.evaluate_many(R32, c["vertices"], phase) == _single(c, R32.astype(np.float64))
