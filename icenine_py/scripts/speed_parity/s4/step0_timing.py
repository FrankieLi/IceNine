"""S4 Step 0 timings and data facts (run alone; preflight first, label quiet or contended).

* ``BatchedScreen`` per voxel as the reconstructor runs it (chunks of 32 FZ points), 10 seeded
  Phase D voxels: wall time split into the torch geometry (stages A-C) and the per-orientation
  stage-D loop, from a timed copy of ``ScreenEvaluator._counts`` whose output is checked equal to
  the library's on every chunk. Peaks per orientation (after the eta filter and the frame lookup)
  and screen passers per voxel.
* Images: size, lit pixels per frame raw and dilated by the pixel radius, frames that occur.
* Re-cost: ``evaluate`` of every passer with the level-0 local cost function, 3 voxels.
* A no-start search (``reconstruct_voxel``) on the same 3 voxels: screen, re-cost and the rest
  (quick MC per level from the reconstructor's own level prints; FindOptimal, VarianceMinimizing
  and the final evaluation as the remainder). One BFS neighbour step (``local_optimization`` with
  the BFS's CMA budget from a start 0.3 deg off the truth) on 10 voxels.

Usage (from icenine_py/): uv run python scripts/speed_parity/s4/step0_timing.py [--n-search 3]
Output: benchmarks/s4_translation/step0_timing.json, preflight_step0.json.
"""

import argparse
import contextlib
import io
import json
import re
import sys
import time
from pathlib import Path
from typing import Any, Dict, List, Tuple

import numpy as np
import torch

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import common as C  # noqa: E402
import screen_kernel as SK  # noqa: E402

from icenine import cost_functions as CF  # noqa: E402
from icenine import reconstructor as REC  # noqa: E402


def timed_counts(ev: Any, omegas: Any, normals: Any, seg: Any, n_seg: int, verts: Any) -> Any:
    """ScreenEvaluator._counts with timers: (po, pon, t_geom, t_loop). Same operations."""
    t0 = time.perf_counter()
    cf = ev.cf
    po = np.zeros(n_seg, np.int64)
    pon = np.zeros(n_seg, np.int64)
    sim, sample = cf.simulator, cf.sample
    beam_dir = sim.beam_direction
    base_rot = sample.sample_to_lab_matrix[:3, :3]
    cos_w, sin_w = torch.cos(omegas), torch.sin(omegas)
    N = omegas.shape[0]
    Rz = torch.zeros(N, 3, 3)
    Rz[:, 0, 0], Rz[:, 0, 1], Rz[:, 1, 0], Rz[:, 1, 1], Rz[:, 2, 2] = (
        cos_w,
        -sin_w,
        sin_w,
        cos_w,
        1.0,
    )
    full_rot = Rz @ base_rot
    lab_normals = torch.bmm(full_rot, normals.unsqueeze(-1)).squeeze(-1)
    dot = (beam_dir.unsqueeze(0) * lab_normals).sum(dim=1, keepdim=True)
    reflected = beam_dir.unsqueeze(0) - 2.0 * dot * lab_normals
    rd_norms = torch.norm(reflected, dim=1)
    safe = torch.where(rd_norms > 0, rd_norms, torch.ones_like(rd_norms))
    eta = torch.atan2(torch.abs(reflected[:, 1]) / safe, torch.abs(reflected[:, 2]) / safe)
    valid = (eta < cf.eta_limit) & (rd_norms > 0)
    if cf.min_sin_eta > 0.0:
        valid = valid & (torch.sin(eta) >= cf.min_sin_eta)
    bins = CF._omega_bins(omegas.numpy(), cf.range_map.low, cf.range_map.width, cf.cpp_edge_rules)
    ok_bin = (bins >= 0) & (bins < ev._n_bins)
    wedge = np.where(ok_bin, ev._index[np.clip(bins, 0, ev._n_bins - 1)], -1)
    keep = np.nonzero(valid.numpy() & (wedge >= 0))[0]
    order = keep[np.argsort(seg.numpy()[keep], kind="stable")]
    o_idx = torch.from_numpy(order)
    v_wedge = wedge[order]
    v_seg = seg.numpy()[order]
    M = len(order)
    pk = dict(full_rot=full_rot[o_idx], reflected=reflected[o_idx])
    centers_a, det_hit_a, _ = C.centre_pixels_f32(cf, pk, verts)
    n_det = len(cf.detector_list)
    unique = np.unique(v_wedge)
    local_of = np.full(int(v_wedge.max()) + 1, -1, dtype=np.int32)
    local_of[unique] = np.arange(len(unique), dtype=np.int32)
    local_wedge = local_of[v_wedge].astype(np.int32)
    image_list = [
        cf.exp_data.get_image(int(w), d).get_binary_numpy() for w in unique for d in range(n_det)
    ]
    first = cf.exp_data.get_image(int(unique[0]), 0)
    rows, cols = first.num_rows, first.num_cols
    bounds = np.nonzero(np.diff(v_seg))[0] + 1
    starts = np.concatenate([[0], bounds])
    ends = np.concatenate([bounds, [M]])
    t1 = time.perf_counter()
    for s, e in zip(starts, ends):
        res = CF._c_stage_d_overlap(
            image_list,
            local_wedge[s:e],
            det_hit_a[s:e],
            centers_a[s:e],
            None,
            n_det,
            cf.pixel_radius,
            int(e - s),
            rows,
            cols,
            *((ev._rule_bits,) if ev._rule_bits else ()),
        )
        po[v_seg[s]] = res[2]
        pon[v_seg[s]] = res[3]
    t2 = time.perf_counter()
    return po, pon, t1 - t0, t2 - t1


def bragg_part(ev: Any, orientations: np.ndarray, phase: int) -> Tuple[Any, ...]:
    """overlaps_exact up to the call of _counts (timed with the geometry)."""
    from icenine.diffraction_core import get_scattering_omegas_torch

    cf = ev.cf
    K = len(orientations)
    g_hkl, g_mag = cf._phase_recip_vecs[phase]
    n = g_hkl.shape[0]
    R = torch.from_numpy(np.asarray(orientations, dtype=np.float64)).float()
    g = torch.matmul(R, g_hkl.T).transpose(1, 2).reshape(K * n, 3)
    br = get_scattering_omegas_torch(
        g, g_mag.repeat(K), cf.simulator.beam_energy, cf.simulator.beam_deflection_chi
    )
    obs = br.observable
    og = g[obs]
    normals = og / torch.norm(og, dim=1, keepdim=True)
    seg_one = torch.arange(K).repeat_interleave(n)[obs]
    omegas = torch.cat([br.omega1[obs], br.omega2[obs]])
    normals = torch.cat([normals, normals])
    seg = torch.cat([seg_one, seg_one])
    return omegas, normals, seg, K


def screen_split(cap: Dict[str, Any], verts: Any, phase: int) -> Dict[str, Any]:
    ev = cap["screen"]
    t_bragg = t_geom = t_loop = 0.0
    po_all = []
    t_lib = 0.0
    for blk in C.chunk_orientations(cap["fz_orientations"], cap["local_grid"], 32):
        t0 = time.perf_counter()
        om, nr, sg, K = bragg_part(ev, blk, phase)
        t_bragg += time.perf_counter() - t0
        po, pon, tg, tl = timed_counts(ev, om, nr, sg, K, verts)
        t_geom += tg
        t_loop += tl
        t1 = time.perf_counter()
        po_l, pon_l = ev.overlaps_exact(blk, verts, phase)
        t_lib += time.perf_counter() - t1
        assert np.array_equal(po, po_l) and np.array_equal(pon, pon_l)
        po_all.append(po)
    po = np.concatenate(po_all)
    return dict(
        t_library_s=t_lib,
        t_copy_total_s=t_bragg + t_geom + t_loop,
        t_bragg_s=t_bragg,
        t_geometry_s=t_geom,
        t_stage_d_loop_s=t_loop,
        passers=int((po > 0).sum()),
    )


def search_split(ds: C.DataSet, vid: int) -> Dict[str, Any]:
    """One no-start reconstruct_voxel with the screen and re-cost timed inside the discrete
    stage of every level (proxies around the screen object and the local evaluate)."""
    v = ds.voxels[vid]
    acc = dict(screen=0.0, recost=0.0, recost_calls=0, discrete=0.0)
    orig = REC.run_discrete_search_spaced

    class ScreenProxy:
        def __init__(self, ev: Any) -> None:
            self.ev = ev

        def overlaps_exact(self, *a: Any) -> Any:
            t0 = time.perf_counter()
            r = self.ev.overlaps_exact(*a)
            acc["screen"] += time.perf_counter() - t0
            return r

    class LocalProxy:
        def __init__(self, cf: Any) -> None:
            self.cf = cf

        def evaluate(self, **kw: Any) -> Any:
            t0 = time.perf_counter()
            r = self.cf.evaluate(**kw)
            acc["recost"] += time.perf_counter() - t0
            acc["recost_calls"] += 1
            return r

    def wrapped(**kw: Any) -> Any:
        if "screen" in kw:
            kw["screen"] = ScreenProxy(kw["screen"])
        kw["local_cost_fn"] = LocalProxy(kw["local_cost_fn"])
        t0 = time.perf_counter()
        r = orig(**kw)
        acc["discrete"] += time.perf_counter() - t0
        return r

    REC.run_discrete_search_spaced = wrapped  # type: ignore[assignment]
    buf = io.StringIO()
    try:
        t0 = time.perf_counter()
        with contextlib.redirect_stdout(buf):
            res = ds.rec.reconstruct_voxel(
                voxel_vertices=C.voxel_vertices(v),
                phase_index=v.phase,
                rng=np.random.default_rng(0),
            )
        wall = time.perf_counter() - t0
    finally:
        REC.run_discrete_search_spaced = orig  # type: ignore[assignment]
    txt = buf.getvalue()
    quick = [float(x) for x in re.findall(r"quick MC \(([\d.]+)s\)", txt)]
    import stats

    err = float(stats.misorientation_deg_cubic(res.orientation[None], ds.truth[vid][None])[0])
    return dict(
        vid=vid,
        wall_s=wall,
        screen_s=acc["screen"],
        recost_s=acc["recost"],
        recost_calls=acc["recost_calls"],
        discrete_other_s=acc["discrete"] - acc["screen"] - acc["recost"],
        quick_mc_s=sum(quick),
        quick_mc_levels_s=quick,
        rest_s=wall - acc["discrete"] - sum(quick),
        err_deg=err,
        cost=float(res.cost),
    )


def neighbour_steps(ds: C.DataSet, vids: List[int]) -> List[float]:
    from scipy.spatial.transform import Rotation

    rng = np.random.default_rng(3)
    out = []
    for vid in vids:
        v = ds.voxels[vid]
        ax = rng.normal(size=3)
        ax /= np.linalg.norm(ax)
        start = Rotation.from_rotvec(np.radians(0.3) * ax).as_matrix() @ ds.truth[vid]
        t0 = time.perf_counter()
        with contextlib.redirect_stdout(io.StringIO()):
            ds.rec.local_optimization(
                C.voxel_vertices(v),
                v.phase,
                start.astype(np.float32),
                rng,
                cma_max_evals=int(ds.rec.params.cma_neighbor_max_evals),
            )
        out.append(time.perf_counter() - t0)
    return out


def image_facts(ds: C.DataSet, cap: Dict[str, Any]) -> Dict[str, Any]:
    gcf = cap["global_cost_fn"]
    r = int(gcf.pixel_radius)
    packed = np.load(SK.bitmap_path(ds.name, r), mmap_mode="r")
    raw, dil = [], []
    for f in range(packed.shape[0]):
        for d in range(packed.shape[1]):
            raw.append(int((gcf.exp_data.get_image(f, d).get_binary_numpy() != 0).sum()))
            dil.append(int(np.unpackbits(packed[f, d], axis=-1).sum()))
    img = gcf.exp_data.get_image(0, 0).get_binary_numpy()
    n_px = img.size
    return dict(
        rows=int(img.shape[0]),
        cols=int(img.shape[1]),
        n_frames=int(packed.shape[0]),
        n_det=int(packed.shape[1]),
        lit_raw_per_frame_median=float(np.median(raw)),
        lit_raw_fraction_median=float(np.median(raw) / n_px),
        lit_dilated_per_frame_median=float(np.median(dil)),
        lit_dilated_fraction_median=float(np.median(dil) / n_px),
        lit_dilated_fraction_max=float(max(dil) / n_px),
        pixel_radius=r,
    )


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-screen", type=int, default=10)
    ap.add_argument("--n-search", type=int, default=3)
    ap.add_argument("--n-neighbour", type=int, default=10)
    a = ap.parse_args()
    pf = C.preflight_json(C.OUT / "preflight_step0.json")
    label = C.label_of(pf)
    ds = C.load_set("phase_d", threads=1)
    vids = C.timed_voxels(ds, a.n_screen)
    res: Dict[str, Any] = dict(label=label, preflight=pf.get("loadavg"), threads=1, voxels=vids)
    cap0 = C.capture_level0(ds, ds.voxels[vids[0]])
    res["images"] = image_facts(ds, cap0)
    tab = SK.load_table(SK.table_path("phase_d"))
    per_o = np.bincount(tab.seg, minlength=tab.n_orient)
    res["peaks_per_orientation"] = dict(
        mean=float(per_o.mean()),
        median=float(np.median(per_o)),
        min=int(per_o.min()),
        max=int(per_o.max()),
    )
    rows = []
    for vid in vids:
        v = ds.voxels[vid]
        cap = C.capture_level0(ds, v)
        r = screen_split(cap, C.voxel_vertices(v), v.phase)
        r["vid"] = vid
        rows.append(r)
        print(json.dumps(r), flush=True)
    res["screen"] = rows
    # re-cost of every passer, 3 voxels
    rc = []
    for vid in vids[: a.n_search]:
        v = ds.voxels[vid]
        cap = C.capture_level0(ds, v)
        z = (
            np.load(C.CACHE / "ref" / f"phase_d_{vid}.npz")
            if (C.CACHE / "ref" / f"phase_d_{vid}.npz").exists()
            else None
        )
        verts = C.voxel_vertices(v)
        if z is not None:
            passers = z["passers"]
        else:
            po, _ = C.screen_reference(cap, verts, v.phase)
            passers = np.nonzero(po > 0)[0]
        R = C.all_orientations(cap)
        lcf = cap["local_cost_fn"]
        t0 = time.perf_counter()
        for k in passers:
            lcf.evaluate(orientation=R[k], voxel_vertices=verts, phase_index=v.phase)
        dt = time.perf_counter() - t0
        rc.append(
            dict(vid=vid, n=int(len(passers)), s=dt, us_per_eval=1e6 * dt / max(len(passers), 1))
        )
        print(json.dumps(rc[-1]), flush=True)
    res["recost"] = rc
    ss = []
    for vid in vids[: a.n_search]:
        s = search_split(ds, vid)
        ss.append(s)
        print(json.dumps(s), flush=True)
    res["search"] = ss
    nb = neighbour_steps(ds, vids[: a.n_neighbour])
    res["neighbour_s"] = nb
    res["neighbour_median_s"] = float(np.median(nb))
    print("neighbour", nb, flush=True)
    res["preflight_end"] = C.preflight_json(C.OUT / "preflight_step0_end.json").get("loadavg")
    (C.OUT / "step0_timing.json").write_text(json.dumps(res, indent=1) + "\n")


if __name__ == "__main__":
    main()
