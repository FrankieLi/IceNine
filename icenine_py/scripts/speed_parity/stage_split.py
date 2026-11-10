"""Split one VoxelCostFunction.evaluate() into stages, without profiler overhead.

Times the real library calls on captured inputs (Example2.ThreeVoxels voxel 0, one orientation near
the truth, single thread) plus the two Python loops that S1 removed from
calculate_diffraction_overlap_batched (kept here as reference loops) and what replaced them. Runs on
the base commit and on S1 (the S1 lines are skipped when the library has no FrameTable).

Usage (from icenine_py/): uv run python scripts/speed_parity/stage_split.py
"""

import sys
import time
from pathlib import Path

import numpy as np
import torch

sys.path.insert(0, str(Path(__file__).resolve().parent))
import prof_eval as P  # noqa: E402  (builds the ThreeVoxels setup)
from icenine import cost_functions as cf  # noqa: E402
from icenine.diffraction_core import get_scattering_omegas_torch  # noqa: E402


def t_us(f, n=400):
    for _ in range(20):
        f()
    t0 = time.perf_counter()
    for _ in range(n):
        f()
    return (time.perf_counter() - t0) / n * 1e6


def capture_overlap_inputs(fn, o):
    """Keyword arguments evaluate() passes to calculate_diffraction_overlap_batched."""
    cap = {}
    orig = cf.calculate_diffraction_overlap_batched

    def spy(**kw):
        cap.update(kw)
        return orig(**kw)

    cf.calculate_diffraction_overlap_batched = spy
    try:
        fn.evaluate(o, P.verts, P.vox.phase)
    finally:
        cf.calculate_diffraction_overlap_batched = orig
    return cap


def capture_stage_d_args(fn, o):
    """Positional arguments of the stage-D C call made by one evaluate()."""
    cap = {}
    orig = cf._c_stage_d_overlap

    def spy(*a):
        cap["a"] = a
        return orig(*a)

    cf._c_stage_d_overlap = spy
    try:
        fn.evaluate(o, P.verts, P.vox.phase)
    finally:
        cf._c_stage_d_overlap = orig
    return cap["a"]


def main():
    o = P.orients[0]
    fn = P.make(0, 8.0, True)  # final cost: triangle, Q-max 8, C++ rule
    g, m = fn._phase_recip_vecs[P.vox.phase]
    cap = capture_overlap_inputs(fn, o)

    total = t_us(lambda: fn.evaluate(o, P.verts, P.vox.phase))
    bragg = t_us(
        lambda: get_scattering_omegas_torch(
            (torch.from_numpy(o).float() @ g.T).T, m, P.sim.beam_energy, P.sim.beam_deflection_chi
        )
    )
    overlap = t_us(lambda: cf.calculate_diffraction_overlap_batched(**cap))

    # Stage A: per-peak Python loop vs a lookup table
    om = cap["peak_omegas"].numpy()
    rm = cap["range_map"]
    bins = ((om - rm.low) / rm.width).astype(int)
    il = rm.index_list
    nb = len(il)

    def stage_a_loop():
        w = np.full(len(om), -1, dtype=np.int64)
        for i in range(len(om)):
            n = bins[i]
            if 0 <= n < nb and il[n] is not None:
                w[i] = il[n]
        return w

    lut = np.array([-1 if x is None else x for x in il] + [-1], dtype=np.int64)

    def stage_a_lut():
        b = ((om - rm.low) / rm.width).astype(int)
        return lut[np.where((b >= 0) & (b < nb), b, nb)]

    w = stage_a_loop()
    assert (stage_a_lut() == w).all()
    uw = np.unique(w[w >= 0])
    ed = cap["exp_data"]
    n_det = len(cap["detector_list"])

    def image_list_loop():
        return [ed.get_image(int(ww), d).get_binary_numpy() for ww in uw for d in range(n_det)]

    stage_d = capture_stage_d_args(fn, o)
    fs = P.make(3, 5.0, True)  # coarse screen: pixel radius 3, Q-max 5, C++ rule
    stage_d_screen = capture_stage_d_args(fs, o)
    s1 = hasattr(cf, "FrameTable")

    print(f"peaks: N={len(om)} valid M={int((w >= 0).sum())} unique wedges={len(uw)}")
    print(f"evaluate() total                   {total:7.1f} us")
    print(f"  Bragg solve (torch)              {bragg:7.1f} us")
    print(f"  rest of evaluate body (eta ...)  {total - overlap - bragg:7.1f} us")
    print(f"  calculate_overlap_batched        {overlap:7.1f} us")
    loop_us, lut_us = t_us(stage_a_loop), t_us(stage_a_lut)
    images_us = t_us(image_list_loop)
    stage_d_us = t_us(lambda: cf._c_stage_d_overlap(*stage_d))
    print(f"    stage A Python loop            {loop_us:7.1f} us  (lookup table: {lut_us:.1f} us)")
    n_calls = n_det * len(uw)
    print(f"    image-list Python loop         {images_us:7.1f} us  ({n_calls} get_image calls)")
    print(f"    stage D C call                 {stage_d_us:7.1f} us")
    print(
        f"screen arm: evaluate {t_us(lambda: fs.evaluate(o, P.verts, P.vox.phase)):.1f} us,"
        f" stage D C call {t_us(lambda: cf._c_stage_d_overlap(*stage_d_screen)):.1f} us"
    )
    if s1:
        from icenine.diffraction_core import BraggConstants, scattering_omegas_fixed_beam

        gl = (torch.from_numpy(o).float() @ g.T).T
        consts = BraggConstants(m, P.sim.beam_energy, P.sim.beam_deflection_chi)
        fixed = t_us(lambda: scattering_omegas_fixed_beam(gl, consts))
        # the same reflections through the original function, on the same g (no matmul)
        ref = t_us(
            lambda: get_scattering_omegas_torch(gl, m, P.sim.beam_energy, P.sim.beam_deflection_chi)
        )
        table = fn.frame_table()
        table.ensure(np.asarray(w[w >= 0], dtype=np.int32))  # filled, as in a running search
        vw = np.asarray(w[w >= 0], dtype=np.int32)
        index = fn._frame_index
        lookup = t_us(lambda: cf.frame_of_omegas(om, rm.low, rm.width, index, fn.cpp_edge_rules))
        check = t_us(lambda: (fn.frame_table(), table.ensure(vw)))
        print("S1 replacements (same inputs, same machine state):")
        print(f"    Bragg solve, original function  {ref:7.1f} us (on g only, no matmul)")
        print(f"    Bragg solve, hoisted constants  {fixed:7.1f} us")
        print(f"    stage A, frame_of_omegas        {lookup:7.1f} us")
        print(f"    frame table: identity + filled  {check:7.1f} us")
        # pointer prefetch of the C call: the whole table (n_frames * n_det entries) against a
        # compact list of the frames in use, for the same peaks (stage_d[1] are global indices)
        images, wedge_g = stage_d[0], stage_d[1]
        n_det = len(P.dets)
        uw2 = np.unique(wedge_g)
        compact = [images[int(u) * n_det + d] for u in uw2 for d in range(n_det)]
        compact_args = (compact, np.searchsorted(uw2, wedge_g).astype(np.int32)) + tuple(
            stage_d[2:]
        )
        t_table = t_us(lambda: cf._c_stage_d_overlap(*stage_d))
        t_compact = t_us(lambda: cf._c_stage_d_overlap(*compact_args))
        print(
            f"    stage D C call, table of {len(images)} entries {t_table:.1f} us, compact list of"
            f" {len(compact)} entries {t_compact:.1f} us"
        )


if __name__ == "__main__":
    main()
