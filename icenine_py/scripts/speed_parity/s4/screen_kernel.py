"""S4 Step 2: the table-driven level-0 screen and its equality with BatchedScreen.

The peak table (``common.build_table``) holds every peak of the 43,974 level-0 orientations: the
float32 torch operations of ``ScreenEvaluator._counts`` up to the frame lookup, run in the same
blocks as the reference (32 FZ points times the local grid), then ``u = A x + c`` in float64 from
those directions. The kernel takes ``trunc(A x + c)`` (the float -> int32 cast of the screen), the
plane-hit flag, the lit and bright flags from the dilated bitmap, and the config's qualified-peak
rule; ``peak_overlap[o, v]`` is the number of qualified peaks.

Equality with the Step 0 reference (BatchedScreen, float32) per orientation, for the float64 table
and for a float32 table (A, c and x in float32). Every differing orientation is listed with the
peak, detector and coordinate whose pixel index flipped and that coordinate's float64 distance to
the nearest pixel boundary. Peak-level flips (centre pixel index differs from the reference's,
whether or not the count changed) are counted too.

Usage (from icenine_py/):
  uv run python scripts/speed_parity/s4/screen_kernel.py table  --set phase_d
  uv run python scripts/speed_parity/s4/screen_kernel.py equal  --set phase_d
  uv run python scripts/speed_parity/s4/screen_kernel.py equal  --set three
Output: benchmarks/s4_translation/step2_equality_<set>.json; tables and bitmaps in cache/.
"""

import argparse
import json
import sys
import time
from pathlib import Path
from typing import Any, Dict, List, Tuple

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import common as C  # noqa: E402

N_DET_MAX = 2


def table_path(name: str) -> Path:
    return C.CACHE / f"table_{name}.npz"


def bitmap_path(name: str, r: int) -> Path:
    return C.CACHE / f"bitmaps_{name}_r{r}.npy"


def save_table(tab: C.PeakTable, path: Path, meta: Dict[str, Any]) -> None:
    np.savez(
        path,
        seg=tab.seg.astype(np.int32),
        frame=tab.frame.astype(np.int16),
        A=tab.A,
        c=tab.c,
        hit=tab.hit,
        n_orient=tab.n_orient,
        meta=json.dumps(meta),
    )


def load_table(path: Path) -> C.PeakTable:
    z = np.load(path)
    t = C.PeakTable(
        seg=z["seg"].astype(np.int64),
        frame=z["frame"].astype(np.int64),
        A=z["A"],
        c=z["c"],
        hit=z["hit"],
        la=np.zeros(0),
        lb=np.zeros(0),
        full_rot=np.zeros(0),
        refl=np.zeros(0),
        n_orient=int(z["n_orient"]),
    )
    t.meta = json.loads(str(z["meta"]))
    return t


def get_capture(ds: C.DataSet) -> Dict[str, Any]:
    vid = 0 if ds.name == "three" else C.reference_voxels(ds)[0]
    return C.capture_level0(ds, ds.voxels[vid])


def make_table(ds: C.DataSet, cap: Dict[str, Any]) -> Tuple[C.PeakTable, float]:
    gcf = cap["global_cost_fn"]
    blocks = list(
        C.chunk_orientations(cap["fz_orientations"], cap["local_grid"], int(cap["screen_chunk"]))
    )
    phase = int(ds.voxels[0].phase)
    t0 = time.perf_counter()
    tab = C.build_table(gcf, blocks, phase, extent=C.sample_extent(ds.voxels))
    return tab, time.perf_counter() - t0


def ensure_table(ds: C.DataSet, cap: Dict[str, Any]) -> C.PeakTable:
    p = table_path(ds.name)
    if p.exists():
        return load_table(p)
    tab, dt = make_table(ds, cap)
    gcf = cap["global_cost_fn"]
    meta = dict(
        n_orient=tab.n_orient,
        n_peaks=tab.n_peaks,
        build_s_untimed=dt,
        pixel_radius=int(gcf.pixel_radius),
        cpp_overlap_rule=bool(gcf.cpp_overlap_rule),
        cpp_edge_rules=bool(gcf.cpp_edge_rules),
        center_screen=bool(gcf.center_screen),
        screen_chunk=int(cap["screen_chunk"]),
        n_fz=int(len(cap["fz_orientations"])),
        n_local=int(len(cap["local_grid"])),
    )
    save_table(tab, p, meta)
    tab.meta = meta
    return tab


def ensure_bitmaps(ds: C.DataSet, cap: Dict[str, Any], tab: C.PeakTable) -> C.Bitmaps:
    gcf = cap["global_cost_fn"]
    r = int(gcf.pixel_radius)
    p = bitmap_path(ds.name, r)
    first = gcf.exp_data.get_image(0, 0).get_binary_numpy()
    rows, cols = first.shape
    if p.exists():
        return C.Bitmaps(np.load(p, mmap_mode="r"), r, rows, cols)
    n_frames = int(tab.frame.max()) + 1
    bm = C.build_bitmaps(gcf.exp_data, n_frames, len(gcf.detector_list), r, np.unique(tab.frame))
    np.save(p, bm.packed)
    return bm


def cmd_table(a: argparse.Namespace) -> None:
    ds = C.load_set(a.set)
    cap = get_capture(ds)
    tab = ensure_table(ds, cap)
    bm = ensure_bitmaps(ds, cap, tab)
    per_o = np.bincount(tab.seg, minlength=tab.n_orient)
    info = dict(
        set=a.set,
        n_orient=tab.n_orient,
        n_peaks=tab.n_peaks,
        peaks_per_orientation_mean=float(per_o.mean()),
        peaks_per_orientation_median=float(np.median(per_o)),
        peaks_per_orientation_min=int(per_o.min()),
        peaks_per_orientation_max=int(per_o.max()),
        hit_fraction=tab.hit.mean(axis=0).tolist(),
        frames_occurring=int(len(np.unique(tab.frame))),
        table_bytes_float64=int(tab.A.nbytes + tab.c.nbytes + tab.hit.nbytes + 4 * tab.n_peaks * 2),
        bitmap_packed_bytes=int(bm.packed.nbytes),
        meta=tab.meta,
    )
    (C.OUT / f"step2_table_{a.set}.json").write_text(json.dumps(info, indent=1) + "\n")
    print(json.dumps(info, indent=1))


def f32_kernel(tab: C.PeakTable, X: np.ndarray, bm: C.Bitmaps, cpp_rule: bool) -> Dict[str, Any]:
    """The same kernel with a float32 table: A, c, x and the product in float32."""
    t32 = C.PeakTable(
        seg=tab.seg,
        frame=tab.frame,
        A=tab.A.astype(np.float32),
        c=tab.c.astype(np.float32),
        hit=tab.hit,
        la=tab.la,
        lb=tab.lb,
        full_rot=tab.full_rot,
        refl=tab.refl,
        n_orient=tab.n_orient,
    )
    X32 = X.astype(np.float32)
    U = (np.einsum("pdaj,vj->pdva", t32.A, X32) + t32.c[:, :, None, :]).astype(np.float32)
    return _kernel_from_U(t32, U.astype(np.float64), bm, cpp_rule)


def _kernel_from_U(
    tab: C.PeakTable, U: np.ndarray, bm: C.Bitmaps, cpp_rule: bool
) -> Dict[str, np.ndarray]:
    P, nd, V, _ = U.shape
    lit = np.zeros((P, V, nd), bool)
    bright = np.zeros((P, V, nd), bool)
    for d in range(nd):
        cx = C.trunc_index(U[:, d, :, 0])
        cy = C.trunc_index(U[:, d, :, 1])
        inb, br = bm.lookup(tab.frame[:, None], d, cx, cy)
        lit[:, :, d] = inb & tab.hit[:, d][:, None]
        bright[:, :, d] = br & tab.hit[:, d][:, None]
    on_det, ovl = C.qualified(lit, bright, cpp_rule)
    po = np.zeros((tab.n_orient, V), np.int64)
    pon = np.zeros((tab.n_orient, V), np.int64)
    np.add.at(po, tab.seg, ovl.astype(np.int64))
    np.add.at(pon, tab.seg, on_det.astype(np.int64))
    return dict(po=po, pon=pon, U=U)


def reference_centres(
    cap: Dict[str, Any], verts: Any, phase: int, tab: C.PeakTable
) -> Tuple[np.ndarray, np.ndarray]:
    """The reference's float32 centre pixels and coordinates of every table peak (same blocks)."""
    gcf = cap["global_cost_fn"]
    cen, cont = [], []
    for blk in C.chunk_orientations(
        cap["fz_orientations"], cap["local_grid"], int(cap["screen_chunk"])
    ):
        pk = C.block_peaks_f32(gcf, blk, phase)
        c_, _, u_ = C.centre_pixels_f32(gcf, pk, verts)
        cen.append(c_)
        cont.append(u_)
    return np.concatenate(cen), np.concatenate(cont)


def cmd_equal(a: argparse.Namespace) -> None:
    ds = C.load_set(a.set)
    cap = get_capture(ds)
    gcf = cap["global_cost_fn"]
    assert not gcf.cpp_edge_rules, "CppEdgeRules form not implemented (the config leaves it off)"
    tab = ensure_table(ds, cap)
    bm = ensure_bitmaps(ds, cap, tab)
    cpp = bool(gcf.cpp_overlap_rule)
    vids = list(range(len(ds.voxels))) if ds.name == "three" else C.reference_voxels(ds)
    rows: List[Dict[str, Any]] = []
    listed: List[Dict[str, Any]] = []
    t0 = time.time()
    for vid in vids:
        z = np.load(C.CACHE / "ref" / f"{ds.name}_{vid}.npz")
        po_ref, pon_ref = z["po"].astype(np.int64), z["pon"].astype(np.int64)
        v = ds.voxels[vid]
        verts = C.voxel_vertices(v)
        X = C.centroid64(verts)[None]
        out = {}
        for prec in ("f64", "f32"):
            if prec == "f64":
                U = np.einsum("pdaj,vj->pdva", tab.A, X) + tab.c[:, :, None, :]
                out[prec] = _kernel_from_U(tab, U, bm, cpp)
            else:
                out[prec] = f32_kernel(tab, X, bm, cpp)
        cen_ref, cont_ref = reference_centres(cap, verts, v.phase, tab)
        U64 = out["f64"]["U"][:, :, 0, :]  # (P, nd, 2)
        row: Dict[str, Any] = dict(vid=int(vid), points_up=bool(v.points_up))
        for prec in ("f64", "f32"):
            po = out[prec]["po"][:, 0]
            pon = out[prec]["pon"][:, 0]
            Uc = out[prec]["U"][:, :, 0, :]
            idx = C.trunc_index(Uc)
            peak_flip = idx != cen_ref  # (P, nd, 2)
            row[f"{prec}_peak_overlap_diff"] = int((po != po_ref).sum())
            row[f"{prec}_peak_on_det_diff"] = int((pon != pon_ref).sum())
            row[f"{prec}_pass_diff"] = int(((po > 0) != (po_ref > 0)).sum())
            row[f"{prec}_peak_coord_flips"] = int(peak_flip.sum())
            fl = U64[peak_flip]
            row[f"{prec}_flip_dist_max"] = (
                float(np.abs(fl - np.round(fl)).max()) if fl.size else None
            )
            for o in np.nonzero(po != po_ref)[0][:200]:
                ps = np.nonzero(tab.seg == o)[0]
                fp = np.argwhere(peak_flip[ps])
                items = []
                for pi, d, ax in fp:
                    u = float(U64[ps[pi], d, ax])
                    items.append(
                        dict(
                            peak=int(ps[pi]),
                            det=int(d),
                            axis="col" if ax == 0 else "row",
                            u64=u,
                            u_ref32=float(cont_ref[ps[pi], d, ax]),
                            dist_to_boundary=abs(u - round(u)),
                        )
                    )
                listed.append(
                    dict(
                        table=prec,
                        vid=int(vid),
                        orientation=int(o),
                        po_ref=int(po_ref[o]),
                        po_new=int(po[o]),
                        pass_flip=bool((po[o] > 0) != (po_ref[o] > 0)),
                        flipped=items,
                    )
                )
        row["passers_ref"] = int((po_ref > 0).sum())
        rows.append(row)
        print(json.dumps(row), flush=True)
    res = dict(
        set=a.set,
        n_voxels=len(rows),
        n_orient=tab.n_orient,
        n_peaks=tab.n_peaks,
        rows=rows,
        listed=listed,
        totals={
            k: int(sum(r[k] for r in rows))
            for k in rows[0]
            if k.endswith("_diff") or k.endswith("_flips")
        },
        wall_s=time.time() - t0,
    )
    (C.OUT / f"step2_equality_{a.set}.json").write_text(json.dumps(res, indent=1) + "\n")
    print(json.dumps(res["totals"], indent=1))


def main() -> None:
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    for name in ("table", "equal"):
        p = sub.add_parser(name)
        p.add_argument("--set", choices=("three", "phase_d"), required=True)
    a = ap.parse_args()
    C.CACHE.mkdir(parents=True, exist_ok=True)
    dict(table=cmd_table, equal=cmd_equal)[a.cmd](a)


if __name__ == "__main__":
    main()
