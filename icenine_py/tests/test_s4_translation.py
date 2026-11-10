"""Smoke test of the S4 translation-screen prototype (scripts/speed_parity/s4/, ThreeVoxels).

Pins three facts the S4 study rests on (docs/todo_translation_screen_prototype.md):
1. the float64 affine law: for a fixed orientation each peak's projected vertices move by A_p dx
   (u(x + dx) = u(x) + A_p dx), and the float32 path of ``evaluate`` agrees to float32 rounding;
2. the dilated bitmap plus the box-bounds test equals the C pixel-radius test
   (``pixel_radius_overlap_uint8`` through ``stage_d_overlap``), border band included;
3. the table-driven screen equals ``BatchedScreen`` on a few FZ points.
"""

import os
import sys
from pathlib import Path
from typing import Any, Dict, Iterator

import numpy as np
import pytest
import torch

S4 = Path(__file__).resolve().parents[1] / "scripts" / "speed_parity" / "s4"
EXAMPLE = Path(__file__).resolve().parents[2] / "Examples" / "Example2.ThreeVoxels"


@pytest.fixture(scope="module")
def s4() -> Iterator[Dict[str, Any]]:
    data = EXAMPLE / "ScatteringData_Python"
    if not data.exists() or len(list(data.glob("*.d*"))) < 360:
        pytest.skip("ThreeVoxels forward-simulated images not found")
    from icenine import cost_functions as CF

    if not CF._HAS_C_RASTERIZE:
        pytest.skip("needs the _rasterize C extension")
    # other script directories also have a top-level ``common`` module: import the S4 modules with
    # those names hidden, then restore sys.modules (the S4 modules keep their own references)
    names = ("common", "translation_correctness")
    saved = {n: sys.modules.pop(n) for n in names if n in sys.modules}
    sys.path.insert(0, str(S4))
    old = os.getcwd()
    try:
        import common as C
        import translation_correctness as TC

        ds = C.load_set("three", threads=1)
        cap = C.capture_level0(ds, ds.voxels[0])
        yield dict(C=C, TC=TC, ds=ds, cap=cap)
    finally:
        os.chdir(old)
        sys.path.remove(str(S4))
        for n in names:
            sys.modules.pop(n, None)
        sys.modules.update(saved)


def test_affine_law_float64(s4: Dict[str, Any]) -> None:
    C, TC, ds, cap = s4["C"], s4["TC"], s4["ds"], s4["cap"]

    lcf = cap["local_cost_fn"]
    geoms = C.det_geoms(lcf)
    t_s = C.sample_translation(lcf)
    v = ds.voxels[0]
    verts = C.voxel_vertices(v)
    X = verts.double().numpy()
    rng = np.random.default_rng(0)
    R_all = C.all_orientations(cap)
    with TC.Spy() as spy:
        for R in [ds.truth[0], R_all[5], R_all[30001]]:
            c0 = spy.run(lcf, R, verts, v.phase)
            keep = TC.frame_valid_mask(lcf, c0["omegas"]) >= 0
            full, refl = TC.f64_geometry(lcf, c0["omegas"][keep], c0["normals"][keep])
            A, c, _, _, _ = C.affine_maps(full, refl, t_s, geoms)
            U0, _ = C.project_direct64(full, refl, t_s, geoms, X)
            for _ in range(2):
                dx = np.array([*rng.uniform(-0.3, 0.3, 2), 0.0])
                U1, _ = C.project_direct64(full, refl, t_s, geoms, X + dx)
                pred = U0 + np.einsum("pdaj,j->pda", A, dx)[:, :, None, :]
                assert np.abs(U1 - pred).max() < 1e-9
                c1 = spy.run(lcf, R, verts + torch.from_numpy(dx).float(), v.phase)
                assert np.array_equal(c1["wedge"], c0["wedge"])  # same peaks, same frames
                assert np.abs(c1["verts"] - pred).max() < 5e-3  # float32 path


def test_dilated_bitmap_equals_pixel_radius(s4: Dict[str, Any]) -> None:
    C, cap = s4["C"], s4["cap"]
    from icenine import cost_functions as CF

    gcf = cap["global_cost_fn"]
    r = int(gcf.pixel_radius)
    img = None
    for f in range(0, 180, 7):
        cand = gcf.exp_data.get_image(f, 0).get_binary_numpy()
        if cand.any():
            img = cand
            break
    assert img is not None
    rows, cols = img.shape
    dil = C.dilate(img, r)
    assert np.array_equal(dil, C.dilate_sparse(img, r))
    lit_rows = np.unique(np.nonzero(img)[0])
    ys = np.unique(
        np.concatenate([np.arange(-r - 2, 6), np.arange(rows - 6, rows + r + 2), lit_rows[:20] + 1])
    )
    cxs = np.arange(-r - 2, cols + r + 2)
    n_checked = 0
    for cy in ys:
        inb = (-r <= cy < rows + r) & (cxs >= -r) & (cxs < cols + r)
        br = np.zeros(len(cxs), bool)
        if -r <= cy < rows + r:
            br[inb] = dil[cy + r, cxs[inb] + r] != 0
        code = inb.astype(np.int8) * 2 + br
        cuts = np.nonzero(np.diff(code))[0] + 1
        for s, e in zip(np.r_[0, cuts], np.r_[cuts, len(cxs)]):
            M = int(e - s)
            centers = np.zeros((M, 1, 2), np.int32)
            centers[:, 0, 0] = cxs[s:e]
            centers[:, 0, 1] = cy
            out = CF._c_stage_d_overlap(
                [img], np.zeros(M, np.int32), np.ones((M, 1), np.uint8), centers, None,
                1, r, M, rows, cols,
            )  # fmt: skip
            assert out[3] == M * int(inb[s]) and out[2] == M * int(br[s])
            n_checked += M
    assert n_checked > 50_000


def test_table_screen_equals_batched_screen(s4: Dict[str, Any]) -> None:
    C, ds, cap = s4["C"], s4["ds"], s4["cap"]
    gcf = cap["global_cost_fn"]
    v = ds.voxels[0]
    verts = C.voxel_vertices(v)
    fz = cap["fz_orientations"][:64]  # two chunks of 32 FZ points
    blocks = list(C.chunk_orientations(fz, cap["local_grid"], 32))
    po_ref = np.concatenate([cap["screen"].overlaps_exact(b, verts, v.phase)[0] for b in blocks])
    tab = C.build_table(gcf, blocks, v.phase, extent=C.sample_extent(ds.voxels))
    bm = C.build_bitmaps(
        gcf.exp_data, int(tab.frame.max()) + 1, 2, int(gcf.pixel_radius), np.unique(tab.frame),
        sparse=True,
    )  # fmt: skip
    out = C.kernel_numpy(tab, C.centroid64(verts)[None], bm, bool(gcf.cpp_overlap_rule))
    assert np.array_equal(out["po"][:, 0], po_ref)
    assert (po_ref > 0).any()


def test_qualified_rules_match_library_truth_table(s4: Dict[str, Any]) -> None:
    """The kernels' default and C++ qualified-peak rules (``common.qualified``) against the
    library's ``count_qualified_peaks`` and ``count_qualified_peaks_cpp`` on every (lit, bright)
    pattern of 2 and 3 detectors (bright implies lit); with 2 detectors the default overlap is
    'bright on detector 0' and the C++ one 'bright on 0 and (bright on 1 or not lit on 1)'."""
    import itertools

    from icenine import cost_functions as CF

    C = s4["C"]
    for nd in (2, 3):
        for lit in itertools.product((False, True), repeat=nd):
            for bright in itertools.product((False, True), repeat=nd):
                if any(b and not li for b, li in zip(bright, lit)):
                    continue
                L, B = np.array([lit]), np.array([bright])
                on0, ov0 = C.qualified(L, B, False)
                on1, ov1 = C.qualified(L, B, True)
                ref0 = CF.count_qualified_peaks(list(lit), list(bright), nd)
                ref1 = CF.count_qualified_peaks_cpp(list(lit), list(bright), nd)
                assert (int(on0[0]), int(ov0[0])) == (ref0[0], ref0[1]), (lit, bright)
                assert (int(on1[0]), int(ov1[0])) == (ref1[0], ref1[1]), (lit, bright)
                if nd == 2:
                    assert bool(ov0[0]) == bright[0]
                    assert bool(ov1[0]) == (bright[0] and (bright[1] or not lit[1]))
