"""Batched local cost (S3 Step 1, BatchedCostEvaluator): the CPU compact variant and the fixed-shape
variant run on the CPU equal ``VoxelCostFunction.evaluate`` in every OverlapInfo field; an MPS run
is only compared (counts reported, no equality)."""

import json
import os
import warnings
from pathlib import Path
from types import SimpleNamespace
from typing import Any, Dict, List

import numpy as np
import pytest
import torch
from scipy.spatial.transform import Rotation

from icenine import batched_cost as BC
from icenine import cost_functions as CF
from icenine.cost_functions import OverlapInfo, VoxelCostFunction

from test_center_screen import three_voxels  # noqa: F401  (fixture)

needs_c = pytest.mark.skipif(not CF._HAS_C_RASTERIZE, reason="C extension missing")
FIELDS = (
    "pixel_overlap",
    "pixel_on_detector",
    "peak_overlap",
    "peak_on_detector",
    "detectors_overlap",
    "quality",
    "n_quality_points",
    "ray_cases",
)


def _local(c: Dict[str, Any], **kw: Any) -> VoxelCostFunction:
    return VoxelCostFunction(pixel_radius=0, mode="hard", **kw, **c["kwargs"])


def orientations(truth: np.ndarray, n_each: int = 24, seed: int = 0) -> np.ndarray:
    """The truth, perturbations of 0.1 to 8 degrees about random axes, and random rotations."""
    rng = np.random.default_rng(seed)
    out = [truth]
    for deg in (0.1, 0.4, 1.0, 3.0, 8.0):
        for _ in range(n_each):
            ax = rng.normal(size=3)
            ax /= np.linalg.norm(ax)
            out.append(Rotation.from_rotvec(ax * np.radians(deg)).as_matrix() @ truth)
    out += [
        Rotation.random(random_state=int(s)).as_matrix() for s in rng.integers(1 << 30, size=n_each)
    ]
    return np.array(out)


def differing(a: List[OverlapInfo], b: List[OverlapInfo]) -> Dict[str, int]:
    """Number of orientations in which each OverlapInfo field differs (exact comparison)."""
    return {f: sum(getattr(x, f) != getattr(y, f) for x, y in zip(a, b)) for f in FIELDS}


def run_chunks(
    ev: BC.BatchedCostEvaluator, R: np.ndarray, verts: torch.Tensor, phase: int, B: int, **kw: Any
) -> List[OverlapInfo]:
    out: List[OverlapInfo] = []
    for s in range(0, len(R), B):
        out += ev.evaluate_many(R[s : s + B], verts, phase, **kw)
    return out


@needs_c
@pytest.mark.parametrize(
    "flags",
    [
        {},
        {"cpp_overlap_rule": True},
        {"min_sin_eta": 0.4},
        {"cpp_edge_rules": True},
        {"cpp_edge_rules": True, "cpp_overlap_rule": True},
    ],
    ids=["default", "cpp_rule", "min_sin_eta", "cpp_edge_rules", "cpp_edge_rules+cpp_rule"],
)
@pytest.mark.parametrize("max_q", [0.0, 5.0], ids=["all_refl", "qmax5"])
@pytest.mark.parametrize("variant", ["compact", "fixed"])
@pytest.mark.parametrize("B", [1, 7, 64])
def test_cpu_variants_equal_evaluate(
    three_voxels: Dict[str, Any], flags: Dict[str, Any], max_q: float, variant: str, B: int
) -> None:  # noqa: F811
    c = three_voxels
    phase = c["voxel"].phase
    cf = _local(c, max_q=max_q, **flags)
    ref_cf = _local(c, max_q=max_q, **flags)
    truth = np.asarray(c["voxel"].orientation, dtype=np.float64)
    R = orientations(truth)
    ref = [ref_cf.evaluate(r, c["vertices"], phase) for r in R]
    ev = BC.BatchedCostEvaluator(cf, "cpu")
    got = run_chunks(ev, R, c["vertices"], phase, B, variant=variant)
    assert len(got) == len(ref)
    assert differing(got, ref) == {f: 0 for f in FIELDS}
    assert got == ref  # the dataclass comparison, every field at once
    assert cf.eval_count == ref_cf.eval_count == len(R)
    # not vacuous: the truth is bright, the far orientations are not, and peaks exist throughout
    assert ref[0].quality > 0.5 and min(i.quality for i in ref) < 0.1
    assert sum(i.peak_on_detector > 0 for i in ref) > len(R) // 2
    assert ev.counts["orient_" + variant] == len(R)


@needs_c
def test_default_variant_is_compact_on_cpu_and_counts_evaluations(
    three_voxels: Dict[str, Any],
) -> None:  # noqa: F811
    c = three_voxels
    cf = _local(c)
    ev = BC.BatchedCostEvaluator(cf, "cpu")
    R = orientations(np.asarray(c["voxel"].orientation, dtype=np.float64))[:100]
    out = ev.evaluate_many(R, c["vertices"], c["voxel"].phase)
    assert len(out) == 100 and cf.eval_count == 100
    assert ev.counts["calls_compact"] == 1 and ev.counts["calls_fixed"] == 0
    assert ev.evaluate_many(R[:0], c["vertices"], c["voxel"].phase) == []
    assert cf.eval_count == 100


@needs_c
def test_unknown_phase_falls_back_to_evaluate(three_voxels: Dict[str, Any]) -> None:  # noqa: F811
    c = three_voxels
    cf = _local(c)
    ev = BC.BatchedCostEvaluator(cf, "cpu")
    R = orientations(np.asarray(c["voxel"].orientation, dtype=np.float64))[:3]
    out = ev.evaluate_many(R, c["vertices"], 7)
    assert out == [OverlapInfo()] * 3 and cf.eval_count == 3


@needs_c
def test_sample_matrix_change_is_seen(three_voxels: Dict[str, Any]) -> None:  # noqa: F811
    """The constants follow an in-place change of the sample's matrix (as evaluate reads it)."""
    c = three_voxels
    cf = _local(c)
    ev = BC.BatchedCostEvaluator(cf, "cpu")
    R = orientations(np.asarray(c["voxel"].orientation, dtype=np.float64))[:5]
    ev.evaluate_many(R, c["vertices"], c["voxel"].phase)
    m = cf.sample.sample_to_lab_matrix
    old = m.clone()
    try:
        m[:3, 3] += torch.tensor([0.0, 0.0, 0.2])
        want = [cf.evaluate(r, c["vertices"], c["voxel"].phase) for r in R]
        got = ev.evaluate_many(R, c["vertices"], c["voxel"].phase)
        assert got == want
    finally:
        m.copy_(old)


class _DarkExp:
    """A stand-in ``exp_data`` whose every frame is dark (one shared zero array)."""

    n_omega_intervals, n_detectors = 180, 2
    zero = np.zeros((2048, 2048), dtype=np.uint8)

    def get_image(self, w: int, d: int) -> Any:
        return SimpleNamespace(num_rows=2048, num_cols=2048, get_binary_numpy=lambda: self.zero)


@needs_c
def test_exp_data_swap_of_a_live_cost_function_is_followed(
    three_voxels: Dict[str, Any],  # noqa: F811
) -> None:
    """S1: the evaluator takes its frames from ``cf.frame_table()``, which follows a swap of
    ``cf.exp_data`` (optimizer_sweep, findoptimal_sweep). Frames filled from the old data before
    the swap are not kept: after the swap ``evaluate_many`` equals a fresh cost function's
    ``evaluate`` on the new data, and after swapping back it equals the first results. (The
    S3 evaluator before S1 read ``cf.exp_data`` at fill time and kept the frames it had already
    filled, so it mixed old and new frames.)"""
    c = three_voxels
    phase, verts = c["voxel"].phase, c["vertices"]
    cf = _local(c)
    ev = BC.BatchedCostEvaluator(cf, "cpu")
    R = orientations(np.asarray(c["voxel"].orientation, dtype=np.float64))[:40]
    first = ev.evaluate_many(R, verts, phase)  # fills frames from the real images
    assert first[0].quality > 0.5 and sum(i.pixel_overlap > 0 for i in first) > 20
    original = c["kwargs"]["exp_data"]
    dark = _DarkExp()
    cf.exp_data = dark  # type: ignore[assignment]
    try:
        swapped = ev.evaluate_many(R, verts, phase)
        fresh = VoxelCostFunction(pixel_radius=0, mode="hard", **{**c["kwargs"], "exp_data": dark})
        assert swapped == [fresh.evaluate(r, verts, phase) for r in R]
        assert all(i.pixel_overlap == 0 and i.quality == 0.0 for i in swapped)
        assert all(i.peak_on_detector > 0 for i in swapped[:5])  # peaks exist, images dark
    finally:
        cf.exp_data = original
    assert ev.evaluate_many(R, verts, phase) == first


@needs_c
def test_unsupported_cost_functions(
    three_voxels: Dict[str, Any], monkeypatch: pytest.MonkeyPatch  # noqa: F811
) -> None:
    from icenine.ray_overlap import RayRule

    c = three_voxels
    assert BC.unsupported_reason(_local(c)) is None
    assert BC.unsupported_reason(_local(c, cpp_overlap_rule=True)) is None
    assert BC.unsupported_reason(_local(c, cpp_edge_rules=True)) is None  # S1 Step 4
    assert BC.unsupported_reason(_local(c, cpp_edge_rules=True, cpp_overlap_rule=True)) is None
    bad = {
        "screen": VoxelCostFunction(pixel_radius=3, **c["kwargs"]),
        "soft": VoxelCostFunction(pixel_radius=0, mode="soft", **c["kwargs"]),
        "ray": _local(c, ray_rule=RayRule()),
    }
    reasons = {"screen": "pixel_radius", "soft": "soft", "ray": "ray rule"}
    for name, cf in bad.items():
        assert reasons[name] in str(BC.unsupported_reason(cf)), name
        with pytest.warns(RuntimeWarning, match="sequential quick MC"):
            assert BC.make_batched_cost_evaluator(cf) is None
        with warnings.catch_warnings():
            warnings.simplefilter("error")
            assert BC.make_batched_cost_evaluator(cf, warn=False) is None
        with pytest.raises(NotImplementedError):
            BC.BatchedCostEvaluator(cf)
    assert isinstance(BC.make_batched_cost_evaluator(_local(c)), BC.BatchedCostEvaluator)
    # a build without the _rasterize extension, or one older than CppOverlapRule: a supported cost
    # function becomes unsupported (the reconstructor then runs the sequential quick MC)
    ok, ok_rule = _local(c), _local(c, cpp_overlap_rule=True)
    ok_edge = _local(c, cpp_edge_rules=True)
    monkeypatch.setattr(CF, "_C_HAS_VERTEX_LIT", False)  # an extension before the exact edge rules
    assert BC.unsupported_reason(ok) is None and BC.unsupported_reason(ok_rule) is None
    assert "predates the exact CppEdgeRules" in str(BC.unsupported_reason(ok_edge))
    with pytest.warns(RuntimeWarning, match="sequential quick MC"):
        assert BC.make_batched_cost_evaluator(ok_edge) is None
    monkeypatch.setattr(CF, "_C_HAS_VERTEX_LIT", True)
    monkeypatch.setattr(CF, "_C_HAS_CPP_RULE", False)
    assert BC.unsupported_reason(ok) is None
    assert "predates CppOverlapRule" in str(BC.unsupported_reason(ok_rule))
    monkeypatch.setattr(CF, "_HAS_C_RASTERIZE", False)
    for cf in (ok, ok_rule):
        assert "_rasterize C extension is missing" in str(BC.unsupported_reason(cf))
        with pytest.warns(RuntimeWarning, match="sequential quick MC"):
            assert BC.make_batched_cost_evaluator(cf) is None
        with pytest.raises(NotImplementedError):
            BC.BatchedCostEvaluator(cf)


def test_device_validation(three_voxels: Dict[str, Any]) -> None:  # noqa: F811
    c = three_voxels
    with pytest.raises(ValueError):
        BC.BatchedCostEvaluator(_local(c), "tpu")
    assert BC.device_available("cpu") and not BC.device_available("tpu")
    if not BC.device_available("cuda"):
        with pytest.raises(RuntimeError):
            BC.BatchedCostEvaluator(_local(c), "cuda")


@needs_c
@pytest.mark.skipif(not BC.device_available("mps"), reason="MPS not available")
@pytest.mark.parametrize("max_q", [0.0, 5.0], ids=["all_refl", "qmax5"])
def test_mps_is_compared_not_equal(
    three_voxels: Dict[str, Any], max_q: float, capsys: Any
) -> None:  # noqa: F811
    """MPS differs from evaluate in the last bits of the transcendental functions: report how many
    orientations differ in each field; assert only that the costs stay close."""
    c = three_voxels
    phase = c["voxel"].phase
    cf = _local(c, max_q=max_q)
    ref_cf = _local(c, max_q=max_q)
    R = orientations(np.asarray(c["voxel"].orientation, dtype=np.float64), n_each=60)
    ref = [ref_cf.evaluate(r, c["vertices"], phase) for r in R]
    ev = BC.BatchedCostEvaluator(cf, "mps")
    got = ev.evaluate_many(R, c["vertices"], phase)  # >= GPU_MIN_BATCH: the fixed-shape variant
    assert ev.counts["calls_fixed"] == 1
    n_diff = differing(got, ref)
    dcost = np.array([abs(a.cost - b.cost) for a, b in zip(got, ref)])
    with capsys.disabled():
        print(
            f"\nMPS vs evaluate ({len(R)} orientations, max_q {max_q}): fields differing {n_diff}; "
            f"cost differs in {(dcost > 0).sum()}, max |dcost| {dcost.max():.4f}"
        )
    save = os.environ.get("ICENINE_S3_MPS_JSON")  # optional: keep the printed counts (a JSON)
    if save:
        path = Path(save)
        rec = json.loads(path.read_text()) if path.exists() else {}
        rec[f"max_q_{max_q:g}"] = dict(
            n_orientations=len(R),
            max_q=max_q,
            fields_differing=n_diff,
            cost_differs=int((dcost > 0).sum()),
            max_abs_dcost=float(dcost.max()),
        )
        rec["torch"] = torch.__version__
        path.write_text(json.dumps(rec, indent=1, sort_keys=True) + "\n")
    assert dcost.max() < 0.1
    small = ev.evaluate_many(R[: BC.GPU_MIN_BATCH - 1], c["vertices"], phase)
    assert ev.counts["calls_compact"] == 1  # below GPU_MIN_BATCH: the CPU compact variant
    assert small == ref[: BC.GPU_MIN_BATCH - 1]  # ... which is exact
