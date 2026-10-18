"""Task D step 1: per-spot provenance of the clean full-sample render (eval-reduction).

A sibling of render_full.py: it re-renders the same clean images (same geometry, beam, omega
range, Q-max and sample as `render_full.py --variants clean`) but through a recording image, so
that every lit pixel is attributed to the (voxel, reflection) pairs that drew it. It writes NO
image files and never touches clean/ or realistic/ (render_full.py is unchanged, so its output is
byte-identical). Script-level instrumentation: the only library change is the opt-in
`ForwardSimulation.peak_hook` (default None), which names the voxel and reflection of each
rasterized peak; the physics and the rasterizer are the stock ones.

Truth written, per (frame, detector):
  * triangles: one record per rasterized (voxel, reflection, omega) triangle on that detector:
    its pixel count, the mean (column, row) index of its pixels, the exact omega, and the id of the
    connected component of lit pixels (4-connected, like noise.py) it lies in;
  * components: pixel count, binary centroid, intensity-weighted centroid, bounding box, number of
    distinct (voxel, reflection) sources and of distinct voxels.

Outputs (under <out-root>/<tag>/provenance/, gitignored with the images): provenance.npz and
provenance_meta.json; the run record and preflight go to scripts/phase_d/results/.

Usage (from icenine_py/):
  uv run python scripts/phase_d/render_provenance.py --workers 10 [--n-voxels 500] [--tag pilot]
"""

import argparse
import json
import multiprocessing as mp
import os
import resource
import sys
import time
from pathlib import Path
from typing import Any, Dict, List

import numpy as np
from scipy import ndimage

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent / "common"))

ROOT = HERE.parents[2]
EX = ROOT / "Examples" / "Example2.ManyGrains"
BASENAME = "500Grains.sim"


def _recording_image_class() -> Any:
    """The recording ImageData subclass (built lazily so that importing this file is cheap)."""
    import torch
    from icenine.image_data import ImageData

    class RecordingImage(ImageData):
        """Dense image that also logs, per triangle, which pixels it filled."""

        def __init__(self, rows: int, cols: int, owner: Dict[str, Any], key: Any) -> None:
            super().__init__(rows, cols)
            self._owner = owner
            self._key = key

        def add_triangle_scanline(self, v0: Any, v1: Any, v2: Any, intensity: Any = 1.0) -> None:
            # steps 0-2 of ImageData.add_triangle_scanline, verbatim (float32 vertices,
            # truncation, Sutherland-Hodgman clip, rounding), then the stock scanline fill
            def trunc(val: float) -> float:
                return -1.0 if val < 0 else float(int(val))

            t = [self._to_tensor(v, dtype=torch.float32) for v in (v0, v1, v2)]
            polygon = [(trunc(x[0].item()), trunc(x[1].item())) for x in t]
            clipped = self._sutherland_hodgman_clip(
                polygon, 0.0, float(self.num_cols - 1), 0.0, float(self.num_rows - 1)
            )
            if len(clipped) < 3:
                return
            pixels = self._scanline_fill([(round(x), round(y)) for x, y in clipped])
            if not pixels:
                return
            val = float(intensity) if isinstance(intensity, (int, float)) else intensity.item()
            arr = np.asarray(pixels, dtype=np.int64)
            ok = (
                (arr[:, 0] >= 0)
                & (arr[:, 0] < self.num_cols)
                & (arr[:, 1] >= 0)
                & (arr[:, 1] < self.num_rows)
            )
            arr = arr[ok]
            if len(arr) == 0:
                return
            # np.add.at-free: a triangle's pixels are unique, so fancy-index add is exact
            self._pixels_dense[arr[:, 1], arr[:, 0]] += val
            cur = self._owner["current"]
            self._owner["tri"].append(
                (
                    self._key[0],
                    self._key[1],
                    cur[0],
                    cur[1],
                    cur[2],
                    len(arr),
                    float(arr[:, 0].mean()),
                    float(arr[:, 1].mean()),
                    int(arr[0, 0]),
                    int(arr[0, 1]),
                )
            )

    return RecordingImage


class _NullImage:
    def add_triangle_scanline(self, *a: Any, **k: Any) -> None:
        return None


def _worker(args: Dict[str, Any]) -> Dict[str, Any]:
    import torch

    torch.set_num_threads(1)
    os.chdir(EX)
    from icenine.config_file import ConfigFile
    from icenine.forward_simulation import ForwardSimulation
    from icenine.sample import Sample
    from icenine.simulation import Simulation

    t0 = time.time()
    w, W = args["worker"], args["workers"]
    cfg = ConfigFile.from_file(str(EX / "ConfigFiles" / "Example2.Simulation.config"))
    cfg.sample_filename = args["mic"]
    if args["max_q"]:
        cfg.max_q = args["max_q"]
    sim = ForwardSimulation(cfg)
    sim.exp_setup.initialize_experiment()
    omega_ranges = sim.exp_setup.get_omega_range_list()
    file_ranges = sim.exp_setup.get_file_range_list()
    dets = sim.exp_setup.get_detector_list()
    range_map = sim.exp_setup.get_range_to_index_map()
    sim.simulator = Simulation(sim.exp_setup)
    sample = Sample()
    sim.exp_setup.initialize_sample(sample, dets[0])
    mic = sample.get_mic()
    if args["voxel_idx"] is not None:
        mic.voxels = [mic.voxels[i] for i in args["voxel_idx"]]
    n_vox = len(mic.voxels)
    voxel_ids = (
        np.asarray(args["voxel_idx"], dtype=np.int64)
        if args["voxel_idx"] is not None
        else np.arange(n_vox, dtype=np.int64)
    )

    n_om = len(omega_ranges)
    owned = [i for i in range(n_om) if i % W == w]
    owned_set = set(owned)
    owner: Dict[str, Any] = {"current": (-1, -1, 0.0), "tri": []}
    rec_cls = _recording_image_class()
    images = [
        [
            rec_cls(d.num_rows, d.num_cols, owner, (i, di)) if i in owned_set else _NullImage()
            for di, d in enumerate(dets)
        ]
        for i in range(n_om)
    ]

    def hook(voxel: int, refl: int, oi: int, omega: float) -> None:
        owner["current"] = (int(voxel_ids[voxel]), refl, omega)

    sim.peak_hook = hook
    t_setup = time.time() - t0
    t1 = time.time()
    sim._simulate_peaks_batched(images, dets, sample, range_map, batch_size=args["batch"])
    t_sim = time.time() - t1

    t2 = time.time()
    tri = np.array(
        owner["tri"],
        dtype=[
            ("frame", "i4"),
            ("det", "i1"),
            ("voxel", "i4"),
            ("refl", "i4"),
            ("omega", "f8"),
            ("npix", "i4"),
            ("c_col", "f8"),
            ("c_row", "f8"),
            ("p_col", "i4"),
            ("p_row", "i4"),
        ],
    )
    comp_rows: List[tuple] = []
    comp_of_tri = np.full(len(tri), -1, dtype=np.int64)
    n_checked = 0
    check_bad = 0
    for di in range(len(dets)):
        for i in owned:
            img = images[i][di]._pixels_dense.numpy()
            lab, n = ndimage.label(img > 0)  # 4-connected, as in noise.py
            if args["check_dir"]:
                from noise import read_ascii_image

                name = f"{BASENAME}{str(file_ranges[di].low + i).zfill(5)}.d{di}"
                ref = read_ascii_image(str(Path(args["check_dir"]) / name), *img.shape)
                n_checked += 1
                check_bad += int(
                    not (
                        np.array_equal(ref > 0, img > 0)
                        and np.allclose(ref, img, rtol=1e-4, atol=1e-5)
                    )
                )
            rr, cc = np.nonzero(lab)
            lb = lab[rr, cc]
            vv = img[rr, cc].astype(np.float64)
            cnt = np.bincount(lb, minlength=n + 1)
            sc = np.bincount(lb, weights=cc, minlength=n + 1)
            sr = np.bincount(lb, weights=rr, minlength=n + 1)
            si = np.bincount(lb, weights=vv, minlength=n + 1)
            sic = np.bincount(lb, weights=vv * cc, minlength=n + 1)
            sir = np.bincount(lb, weights=vv * rr, minlength=n + 1)
            cmin = np.full(n + 1, 1 << 30)
            cmax = np.full(n + 1, -1)
            rmin = np.full(n + 1, 1 << 30)
            rmax = np.full(n + 1, -1)
            np.minimum.at(cmin, lb, cc)
            np.maximum.at(cmax, lb, cc)
            np.minimum.at(rmin, lb, rr)
            np.maximum.at(rmax, lb, rr)
            base = len(comp_rows)
            for L in range(1, n + 1):
                comp_rows.append(
                    (
                        i,
                        di,
                        cnt[L],
                        sc[L] / cnt[L],
                        sr[L] / cnt[L],
                        sic[L] / si[L],
                        sir[L] / si[L],
                        si[L],
                        cmin[L],
                        cmax[L],
                        rmin[L],
                        rmax[L],
                    )
                )
            sel = np.nonzero((tri["frame"] == i) & (tri["det"] == di))[0]
            comp_of_tri[sel] = base + lab[tri["p_row"][sel], tri["p_col"][sel]] - 1
            images[i][di] = None  # type: ignore[assignment]
    comp = np.array(
        comp_rows,
        dtype=[
            ("frame", "i4"),
            ("det", "i1"),
            ("npix", "i4"),
            ("c_col", "f8"),
            ("c_row", "f8"),
            ("i_col", "f8"),
            ("i_row", "f8"),
            ("sum_int", "f8"),
            ("cmin", "i4"),
            ("cmax", "i4"),
            ("rmin", "i4"),
            ("rmax", "i4"),
        ],
    )
    t_lab = time.time() - t2
    rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e9
    return {
        "worker": w,
        "n_voxels": n_vox,
        "tri": tri,
        "comp": comp,
        "comp_of_tri": comp_of_tri,
        "t_setup_s": t_setup,
        "t_sim_s": t_sim,
        "t_label_s": t_lab,
        "peak_rss_gb": rss,
        "n_checked": n_checked,
        "check_bad": check_bad,
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=10)
    ap.add_argument("--seed", type=int, default=0, help="new-orientation sample seed")
    ap.add_argument("--n-voxels", type=int, default=0, help="pilot: voxels nearest the centre")
    ap.add_argument("--tag", default="full")
    ap.add_argument("--batch", type=int, default=2000)
    ap.add_argument("--max-q", type=float, default=0.0, help="override Q-max (0 = config)")
    ap.add_argument("--out-root", default=str(EX / "ScatteringData_PhaseD"))
    ap.add_argument(
        "--check",
        action="store_true",
        help="compare every recorded image with the existing <out-root>/<tag>/clean/ file",
    )
    args = ap.parse_args()
    assert args.workers <= 10

    import preflight

    res_dir = HERE / "results"
    res_dir.mkdir(exist_ok=True)
    info = preflight.preflight()
    (res_dir / f"preflight_provenance_{args.tag}.json").write_text(
        json.dumps(info, indent=2) + "\n"
    )
    try:
        preflight.require_quiet(info=info)
        quiet = True
    except preflight.MachineBusyError as e:
        print("WARNING:", e)
        quiet = False

    mic_rel = f"SimInput/rand_500grains_1mm_neworient_s{args.seed}.mic"
    voxel_idx = None
    if args.n_voxels:
        raw = np.loadtxt(EX / mic_rel, skiprows=1)
        d = np.linalg.norm(raw[:, :2], axis=1)
        voxel_idx = np.argsort(d)[: args.n_voxels].tolist()
    out = Path(args.out_root) / args.tag
    jobs = [
        {
            "worker": w,
            "workers": args.workers,
            "mic": mic_rel,
            "voxel_idx": voxel_idx,
            "batch": args.batch,
            "max_q": args.max_q,
            "check_dir": str(out / "clean") if args.check else None,
        }
        for w in range(args.workers)
    ]
    t0 = time.time()
    ctx = mp.get_context("spawn")
    with ctx.Pool(args.workers) as pool:
        results = pool.map(_worker, jobs)
    wall = time.time() - t0

    # merge: component ids become global
    tri_parts, comp_parts, off = [], [], 0
    for r in results:
        t = r["tri"]
        cid = r["comp_of_tri"] + off
        tri_parts.append((t, cid))
        comp_parts.append(r["comp"])
        off += len(r["comp"])
    tri = np.concatenate([t for t, _ in tri_parts])
    comp_of_tri = np.concatenate([c for _, c in tri_parts])
    comp = np.concatenate(comp_parts)
    order = np.lexsort((comp["det"], comp["frame"]))  # stable component order by (frame, det)
    inv = np.empty_like(order)
    inv[order] = np.arange(len(order))
    comp = comp[order]
    comp_of_tri = inv[comp_of_tri]
    # distinct (voxel, refl) sources and distinct voxels per component
    key = (
        comp_of_tri.astype(np.int64) * (1 << 40)
        + tri["voxel"].astype(np.int64) * (1 << 12)
        + tri["refl"]
    )
    uk = np.unique(key)
    n_src = np.bincount(uk >> 40, minlength=len(comp))
    vk = np.unique(comp_of_tri.astype(np.int64) * (1 << 40) + tri["voxel"].astype(np.int64))
    n_vox_c = np.bincount(vk >> 40, minlength=len(comp))
    pdir = out / "provenance"
    pdir.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        pdir / "provenance.npz",
        **{f"tri_{k}": tri[k] for k in tri.dtype.names},
        tri_comp=comp_of_tri,
        **{f"comp_{k}": comp[k] for k in comp.dtype.names},
        comp_n_src=n_src,
        comp_n_vox=n_vox_c,
    )
    meta = {
        "tag": args.tag,
        "sample_seed": args.seed,
        "n_voxels": results[0]["n_voxels"],
        "workers": args.workers,
        "timing_label": f"contended ({args.workers} workers; machine quiet at start: {quiet})",
        "wall_s": wall,
        "worker_t_sim_s_max": max(r["t_sim_s"] for r in results),
        "worker_t_label_s_max": max(r["t_label_s"] for r in results),
        "worker_peak_rss_gb_max": max(r["peak_rss_gb"] for r in results),
        "n_triangles": int(len(tri)),
        "n_components": int(len(comp)),
        "n_frames_checked_against_clean": int(sum(r["n_checked"] for r in results)),
        "n_frames_check_mismatch": int(sum(r["check_bad"] for r in results)),
        "connectivity": 4,
        "output_dir": str((pdir).relative_to(ROOT)),
    }
    (pdir / "provenance_meta.json").write_text(json.dumps(meta, indent=2) + "\n")
    (res_dir / f"render_provenance_{args.tag}.json").write_text(json.dumps(meta, indent=2) + "\n")
    print(json.dumps(meta, indent=2))


if __name__ == "__main__":
    main()
