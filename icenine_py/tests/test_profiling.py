"""scripts/profiling: the stage instrumentation leaves reconstruct_voxel bit-identical and is
removed again on exit; the interleaving rotation; the verdict rules on synthetic numbers."""

import os
import sys
from pathlib import Path
from typing import Any, Dict

import numpy as np
import pytest

ROOT = Path(__file__).parent.parent
PROF = ROOT / "scripts" / "profiling"
for p in (PROF, ROOT / "scripts"):
    sys.path.insert(0, str(p))

from test_findoptimal_refactor import _build  # noqa: E402  (tests/ is on sys.path under pytest)


@pytest.fixture(scope="module")
def pc():
    """scripts/profiling/prof_common without leaking its OMP / MKL / torch thread settings."""
    import torch

    keep, threads = dict(os.environ), torch.get_num_threads()
    try:
        import prof_common

        yield prof_common
    finally:
        for k in set(os.environ) - set(keep):
            del os.environ[k]
        torch.set_num_threads(threads)


def test_rotated_order_is_a_rotation_and_balanced(pc):
    pipes = ["A", "B", "C", "D", "E"]
    firsts = []
    for i in range(10):
        o = pc.rotated(pipes, i)
        assert sorted(o) == sorted(pipes)
        k = pipes.index(o[0])
        assert o == pipes[k:] + pipes[:k]  # cyclic order kept (ABAB-style interleaving)
        firsts.append(o[0])
    assert all(firsts.count(p) == 2 for p in pipes)


def test_summarize_times_and_wilson(pc):
    s = pc.summarize_times([1.0, 2.0, 3.0, 4.0, 5.0])
    assert s["median"] == 3.0 and s["n"] == 5 and s["mean"] == 3.0
    assert s["p10"] == pytest.approx(1.4) and s["p90"] == pytest.approx(4.6)
    p, lo, hi = pc.wilson(0, 200)
    assert p == 0 and lo == 0 and hi == pytest.approx(0.0188, abs=5e-4)  # as in the studies


def test_instrumented_reconstruct_voxel_is_bit_identical_and_uninstalled(pc):
    """One reconstruct_voxel run with and without the stage wrappers: same orientation, same cost,
    same evaluation counts; the evaluate calls counted by the wrappers equal the reconstructor's
    own counts; every patched attribute is restored."""
    import icenine.reconstructor as RC
    from icenine.cost_functions import VoxelCostFunction
    from icenine.orientation_search import MCOptimizer
    from icenine.reconstructor import _get_voxel_vertices

    rec, voxel, _ = _build()
    verts = _get_voxel_vertices(voxel)
    res0 = rec.reconstruct_voxel(verts, voxel.phase, rng=np.random.default_rng(7))
    counts0 = tuple(rec.last_eval_counts)
    orig = (
        RC.run_discrete_search_spaced,
        MCOptimizer.__dict__["optimize"],
        MCOptimizer.__dict__["variance_minimizing_optimize"],
        VoxelCostFunction.__dict__["evaluate"],
        RC.AdaptiveVoxelReconstructor.__dict__["reconstruct_voxel"],
    )
    T = pc.ST.StageTimer()
    inst = pc.Instrument(T, rec.params.max_mc_steps, net=False)
    with inst.installed():
        before = T.snapshot()
        res1 = rec.reconstruct_voxel(verts, voxel.phase, rng=np.random.default_rng(7))
        ev = inst.eval_counts(before)
    assert np.array_equal(np.asarray(res0.orientation), np.asarray(res1.orientation))
    assert res0.cost == res1.cost
    assert tuple(rec.last_eval_counts) == counts0
    assert (ev["global"], ev["local"]) == (counts0[0], counts0[1])
    assert T.calls["reconstruct"] == 1 and T.calls["discrete_L0"] == 1
    assert T.calls["find_optimal"] >= 1 and T.calls["variance"] >= 1
    assert T.inclusive["reconstruct"] >= T.inclusive["discrete_L0"] > 0
    after = (
        RC.run_discrete_search_spaced,
        MCOptimizer.__dict__["optimize"],
        MCOptimizer.__dict__["variance_minimizing_optimize"],
        VoxelCostFunction.__dict__["evaluate"],
        RC.AdaptiveVoxelReconstructor.__dict__["reconstruct_voxel"],
    )
    assert all(a is b for a, b in zip(orig, after))


def test_network_patches_are_installed_where_looked_up_and_removed(pc):
    import icenine.orientation_eval as OE
    import perturbation_sweep as ps

    T = pc.ST.StageTimer()
    orig = (ps.prepare_nominal, OE.decode_windows, OE.render_windows, OE.make_realistic_dataset)
    inst = pc.Instrument(T, 3500, net=True)
    with inst.installed():
        assert ps.prepare_nominal is not orig[0] and OE.decode_windows is not orig[1]
        assert OE.render_windows is not orig[2] and OE.make_realistic_dataset is not orig[3]
    assert (
        ps.prepare_nominal,
        OE.decode_windows,
        OE.render_windows,
        OE.make_realistic_dataset,
    ) == orig


def _row(rate: float, lo: float, hi: float, med: float) -> Dict[str, Any]:
    return dict(rate=rate, lo=lo, hi=hi, median_err_right=med)


def test_u0_verdict_rule_on_synthetic_numbers(pc):
    import prof_summary as S

    def mk(times: Dict[str, float]) -> Dict[str, Any]:
        return {v: {p: dict(time=dict(median=t)) for p, t in times.items()} for v in S.VARIANTS}

    names = dict(S.U0_ROW)
    rows = {
        "baseline": _row(0.30, 0.25, 0.36, 0.030),
        "F1": _row(0.03, 0.01, 0.06, 0.027),
        "F1b": _row(0.0, 0.0, 0.019, 0.026),
        "e2": _row(0.04, 0.02, 0.077, 0.0285),
        "proxy": _row(0.07, 0.04, 0.11, 0.028),
        "proxy+F1": _row(0.0, 0.0, 0.02, 0.028),
    }
    t2 = dict(
        end_to_end={v: [dict(name=names[p], **r) for p, r in rows.items()] for v in S.VARIANTS},
        multi_seed={
            f"p_i/{v}": dict(rows=[dict(name=S.U0_ROW3[p], **rows[p]) for p in S.U0_ROW3])
            for v in S.VARIANTS
        },
    )
    # proxy+F1 21 s vs F1b 28 s (25% saving) at matched accuracy: helps; vs F1 24 s only 12.5%
    # faster and the intervals overlap (0.0-0.02 vs 0.01-0.06): neither (a) nor (b) holds
    res = mk(
        {"baseline": 20.0, "F1": 24.0, "F1b": 28.0, "e2": 21.0, "proxy": 20.5, "proxy+F1": 21.0}
    )
    out, _ = S.u0_verdicts(res, t2)
    assert out["proxy+F1 vs F1b"]["helps"] is True  # 25% faster at matched accuracy
    assert out["proxy+F1 vs F1"]["helps"] is False
    assert out["proxy vs baseline"]["helps"] is True  # lower wrong, non-overlapping, +2.5% time
    slow = mk(
        {"baseline": 20.0, "F1": 24.0, "F1b": 28.0, "e2": 23.0, "proxy": 20.5, "proxy+F1": 29.0}
    )
    out2, _ = S.u0_verdicts(slow, t2)
    assert out2["e2 vs baseline"]["helps"] is False  # +15% time exceeds the 10% allowance
    assert out2["proxy+F1 vs F1b"]["helps"] is False  # slower than F1b


def test_h1_case_is_bit_identical_with_and_without_the_patches(pc):
    """One U1/U2 case (voxel 0, r = 0.5 deg, direction 0, both variants) of the H1 pipeline:
    network stage + FindOptimal give the same orientation and evaluation count with the stage
    wrappers installed and after they are uninstalled."""
    sweep_raw = ROOT / "benchmarks" / "toy_orientation_sweep" / "perturbation_sweep_raw.npz"
    model = ROOT / "scripts" / "toy_orientation_sweep_model_realistic_s0.pt"
    if not sweep_raw.exists() or not model.exists():
        pytest.skip("sweep raw / network not present")
    import prof_seeded as SD

    wargs, sweep = SD.NH.worker_args([SD.MODEL])
    SD.init_worker(wargs)  # installs the wrappers
    vidx = int(sweep["voxel_indices"][0])
    item = (vidx, 0, 3, 0, [0], 0, sweep["n_roi"][0, 3], sweep["fail_pass1"][0, 3])
    keep = SD.pipes_for
    SD.pipes_for = lambda ri: ["H1"]
    try:
        with_patches = SD.run_task(item)
        SD._S["I"].uninstall()
        without = SD.run_task(item)
    finally:
        SD.pipes_for = keep
        SD._S["T"].uninstall()
    assert len(with_patches) == len(without) == 2
    for a, b in zip(with_patches, without):
        assert a["R_final"] == b["R_final"] and a["evals"] == b["evals"]
        assert a["err"] == b["err"]
    assert sum(a["stages"]["inclusive"].get("finisher", 0) for a in with_patches) > 0
    assert all("find_optimal" not in b["stages"]["inclusive"] for b in without)  # patches gone
