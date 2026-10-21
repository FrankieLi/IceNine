"""S4 Step 3: loop orderings of the level-0 screen of every voxel (owner question 2).

All orderings compute ``peak_overlap[o, v]`` (or the pass bit) with the Step 2 kernel semantics
(``s4kernel.c``; ``common.kernel_numpy`` is the reference):

  O1  voxel-major: each voxel walks the whole table (today's order, table-driven)
  O2  orientation-major: each orientation's peaks over every voxel (stored table, or ``stream``:
      the peaks of a block of b FZ points computed with the torch ops, used and discarded)
  O3  frame-major: the table sorted by frame, each frame's peaks over a voxel tile
  O5  orientation block x voxel tile (tiled O2)
  O4  pixel voting: cost model only (``model``), built only if the model predicts >= 2x
  MPS O2 on the GPU (torch, float32 with u = n_ref + (phi_ref + A (x - x_ref)), n_ref integer)

Subcommands (from icenine_py/):
  check   every ordering, bitmap layout and rule against the numpy reference on the 20 Phase D
          reference voxels and the 3 ThreeVoxels voxels; the streamed O2 and MPS flips
  time    timings on a contiguous 2,000-voxel patch (and 500, 1,000 for linearity), threads
          1..10, block and tile sweeps, uint8 and bit-packed bitmaps, counts and pass-only
  full    the best ordering on the whole layer (24,570 voxels), tiled over voxels
  model   the cost model of every ordering (O4 included) from measured inputs
  rules   pass rates under the default and the C++ qualified-peak rule (reference voxels, patch,
          layer; counts only, not timed)
  time_o3 O3 tile 2,000 and 256 against O2 at 1 and 10 threads (interleaved; review follow-up)
Every timing is preceded by the preflight and labelled quiet or contended.
"""

import argparse
import json
import resource
import sys
import time
from pathlib import Path
from typing import Any, Callable, Dict, List, Optional, Tuple

import numpy as np
import torch

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import common as C  # noqa: E402
import s4kernel as K  # noqa: E402
import screen_kernel as SK  # noqa: E402

PATCH_SIZES = (500, 1000, 2000)
THREADS = (1, 2, 4, 6, 8, 10)


# ---------------------------------------------------------------------------------------------
# set-up
# ---------------------------------------------------------------------------------------------


class Ctx:
    def __init__(self, name: str = "phase_d", threads: int = 1) -> None:
        self.ds = C.load_set(name, threads=threads)
        self.cap = SK.get_capture(self.ds)
        self.gcf = self.cap["global_cost_fn"]
        self.r = int(self.gcf.pixel_radius)
        self.rule = 1 if self.gcf.cpp_overlap_rule else 0
        self.tab = SK.ensure_table(self.ds, self.cap)
        self.packed_all = np.load(SK.bitmap_path(name, self.r), mmap_mode="r")
        img = self.gcf.exp_data.get_image(0, 0).get_binary_numpy()
        self.rows, self.cols = img.shape
        self.bm = C.Bitmaps(self.packed_all, self.r, self.rows, self.cols)
        self.X = np.stack([C.centroid64(C.voxel_vertices(v)) for v in self.ds.voxels])
        z = np.unique(self.X[:, 2])
        assert len(z) == 1, "one layer expected"
        self.z = float(z[0])
        self._ct: Dict[Tuple[str, bool, int], K.CTable] = {}
        self._bms: Dict[str, List[np.ndarray]] = {}

    def bitmaps(self, layout: str) -> List[np.ndarray]:
        if layout not in self._bms:
            if layout == "packed":
                self._bms[layout] = [np.ascontiguousarray(self.packed_all[:, d]) for d in range(2)]
            else:  # uint8; detector 1 only when the rule needs it
                b0 = K.unpack_bitmaps(self.packed_all, 0)
                b1 = K.unpack_bitmaps(self.packed_all, 1) if self.rule == 1 else b0
                self._bms[layout] = [b0, b1]
        return self._bms[layout]

    def ctab(
        self, layout: str = "u8", by_frame: bool = False, rule: Optional[int] = None
    ) -> K.CTable:
        rule = self.rule if rule is None else rule
        key = (layout, by_frame, rule)
        if key not in self._ct:
            if rule == 1 and layout == "u8" and self.rule == 0:
                layout_bms = [
                    K.unpack_bitmaps(self.packed_all, 0),
                    K.unpack_bitmaps(self.packed_all, 1),
                ]
            else:
                layout_bms = self.bitmaps(layout)
            self._ct[key] = K.CTable(
                self.tab,
                self.z,
                layout_bms,
                layout == "packed",
                self.r,
                self.rows,
                self.cols,
                rule,
                order=K.frame_order(self.tab) if by_frame else None,
            )
        return self._ct[key]

    def patch(self, n: int) -> np.ndarray:
        """The n voxels nearest the sample centre (nested patches)."""
        d = np.linalg.norm(self.X[:, :2] - self.X[:, :2].mean(axis=0), axis=1)
        return np.sort(np.argsort(d, kind="stable")[:n])


def run_kind(
    ctx: Ctx,
    kind: str,
    vids: np.ndarray,
    threads: int,
    layout: str,
    block: int = 288,
    tile: int = 256,
) -> np.ndarray:
    ct_ = ctx.ctab(layout, by_frame=(kind == "o3"))
    return K.run(kind, ct_, ctx.X[vids, 0], ctx.X[vids, 1], threads=threads, block=block, tile=tile)


# ---------------------------------------------------------------------------------------------
# O2 streamed (no stored table) and MPS
# ---------------------------------------------------------------------------------------------


def stream_o2(
    ctx: Ctx, vids: np.ndarray, b_fz: int, threads: int = 1
) -> Tuple[np.ndarray, Dict[str, Any]]:
    """O2 without a stored table: per block of b_fz FZ points the float32 torch peaks, the float64
    affine maps and the C kernel; returns the counts and the frame/peak-set changes against the
    chunk-32 table (float32 batch-shape effect)."""
    cap, gcf = ctx.cap, ctx.gcf
    geoms = C.det_geoms(gcf)
    t_s = C.sample_translation(gcf)
    phase = int(ctx.ds.voxels[0].phase)
    xs, ys = ctx.X[vids, 0], ctx.X[vids, 1]
    out = np.zeros((ctx.tab.n_orient, len(vids)), np.uint8)
    t_peaks = t_kernel = 0.0
    n_peaks = 0
    o = 0
    blocks = list(C.chunk_orientations(cap["fz_orientations"], cap["local_grid"], b_fz))
    bms = ctx.bitmaps("u8")
    peak_frames: List[np.ndarray] = []
    for blk in blocks:
        t0 = time.perf_counter()
        pk = C.block_peaks_f32(gcf, blk, phase)
        fr = pk["full_rot"].double().numpy()
        rf = pk["reflected"].double().numpy()
        A, c, par, _, _ = C.affine_maps(fr, rf, t_s, geoms)
        tab = C.PeakTable(
            seg=pk["seg"],
            frame=pk["frame"],
            A=A,
            c=c,
            hit=~par,
            la=np.zeros(0),
            lb=np.zeros(0),
            full_rot=np.zeros(0),
            refl=np.zeros(0),
            n_orient=len(blk),
        )
        ct_ = K.CTable(tab, ctx.z, bms, False, ctx.r, ctx.rows, ctx.cols, ctx.rule)
        t1 = time.perf_counter()
        res = K.run("o2", ct_, xs, ys, threads=threads, block=max(1, len(blk) // max(threads, 1)))
        t2 = time.perf_counter()
        out[o : o + len(blk)] = res
        o += len(blk)
        t_peaks += t1 - t0
        t_kernel += t2 - t1
        n_peaks += tab.n_peaks
        peak_frames.append(pk["frame"])
    fr_all = np.concatenate(peak_frames)
    same_set = len(fr_all) == ctx.tab.n_peaks and bool(np.array_equal(fr_all, ctx.tab.frame))
    frame_changes = int((fr_all != ctx.tab.frame).sum()) if len(fr_all) == ctx.tab.n_peaks else None
    return out, dict(
        b_fz=b_fz,
        n_blocks=len(blocks),
        t_peaks_s=t_peaks,
        t_kernel_s=t_kernel,
        n_peaks=n_peaks,
        peak_set_and_frames_equal_chunk32=same_set,
        frames_changed_vs_chunk32=frame_changes,
    )


class MPSO2:
    """O2 on the GPU (torch MPS), counts, rule 0 (detector 0)."""

    def __init__(self, ctx: Ctx, vids: np.ndarray, block: int = 576) -> None:
        assert ctx.rule == 0, "MPS arm implements the default rule only"
        self.dev = torch.device("mps")
        tab = ctx.tab
        self.ctx = ctx
        self.block = block
        X = ctx.X[vids]
        self.x_ref = X[:, :2].mean(axis=0)
        A = tab.A[:, 0, :, :2]  # (P, 2, 2) det 0, in-layer
        u_ref = np.einsum("paj,j->pa", tab.A[:, 0], np.array([*self.x_ref, ctx.z])) + tab.c[:, 0]
        n_ref = np.floor(u_ref)
        self.A = torch.from_numpy(A.astype(np.float32)).to(self.dev)
        self.n_ref = torch.from_numpy(n_ref.astype(np.int32)).to(self.dev)
        self.phi = torch.from_numpy((u_ref - n_ref).astype(np.float32)).to(self.dev)
        self.hit = torch.from_numpy(tab.hit[:, 0]).to(self.dev)
        b0 = ctx.bitmaps("u8")[0]
        self.rs = int(b0.shape[2])
        self.fstride = int(b0.shape[1] * b0.shape[2])
        self.bm = torch.from_numpy(b0.reshape(-1)).to(self.dev)
        self.base = torch.from_numpy(tab.frame.astype(np.int64) * self.fstride).to(self.dev)
        self.dx = torch.from_numpy((X[:, :2] - self.x_ref).astype(np.float32)).to(self.dev)
        self.pstart = np.searchsorted(tab.seg, np.arange(tab.n_orient + 1))
        self.seg = torch.from_numpy(tab.seg.astype(np.int64)).to(self.dev)
        self.V = len(vids)
        self.n_orient = tab.n_orient

    def run(self) -> np.ndarray:
        r, rows, cols = self.ctx.r, self.ctx.rows, self.ctx.cols
        out = torch.zeros((self.n_orient, self.V), dtype=torch.int32, device=self.dev)
        for o0 in range(0, self.n_orient, self.block):
            o1 = min(o0 + self.block, self.n_orient)
            p0, p1 = int(self.pstart[o0]), int(self.pstart[o1])
            A = self.A[p0:p1]
            w = self.phi[p0:p1, None, :] + torch.einsum("pab,vb->pva", A, self.dx)
            fl = torch.floor(w)
            i = self.n_ref[p0:p1, None, :] + fl.to(torch.int32)
            i = i + ((i < 0) & (w != fl)).to(torch.int32)  # truncation toward zero
            cx, cy = i[..., 0], i[..., 1]
            inb = (cx >= -r) & (cx < cols + r) & (cy >= -r) & (cy < rows + r)
            inb = inb & self.hit[p0:p1, None]
            flat = self.base[p0:p1, None] + ((cy + r) * self.rs + (cx + r)).to(torch.int64)
            flat = torch.where(inb, flat, torch.zeros_like(flat))
            br = (self.bm[flat.reshape(-1)].reshape(flat.shape) != 0) & inb
            out.index_add_(0, self.seg[p0:p1], br.to(torch.int32))
        torch.mps.synchronize()
        return out.to(torch.uint8).cpu().numpy()


# ---------------------------------------------------------------------------------------------
# check
# ---------------------------------------------------------------------------------------------


def cmd_check(a: argparse.Namespace) -> None:
    res: Dict[str, Any] = {}
    for name in ("phase_d", "three"):
        ctx = Ctx(name)
        vids = np.array(C.reference_voxels(ctx.ds) if name == "phase_d" else [0, 1, 2])
        X = ctx.X[vids]
        want = {rule: C.kernel_numpy(ctx.tab, X, ctx.bm, bool(rule))["po"] for rule in (0, 1)}
        rows = {}
        for layout in ("u8", "packed"):
            for rule in (0, 1):
                for kind in ("o1", "o1p", "o2", "o2p", "o3", "o5"):
                    ct_ = ctx.ctab(layout, by_frame=(kind == "o3"), rule=rule)
                    out = K.run(kind, ct_, X[:, 0], X[:, 1], threads=4, block=288, tile=7)
                    w = want[rule] if not kind.endswith("p") else (want[rule] > 0).astype(np.uint8)
                    rows[f"{kind}_{layout}_rule{rule}"] = bool(np.array_equal(out, w))
        if name == "phase_d":  # O4 needs the voxel lattice (the 3 ThreeVoxels voxels are not one)
            up = np.array([bool(v.points_up) for v in ctx.ds.voxels])
            lat = K.Lattice(X, up[vids], float(ctx.ds.voxels[0].side_length))
            o4, _ = K.run_o4(ctx.ctab("packed", rule=0), lat, X[:, 0], X[:, 1], threads=4)
            rows["o4_packed_rule0"] = bool(np.array_equal(o4, want[0]))
        flips = {}
        for b in (1, 32, 128):
            out, info = stream_o2(ctx, vids, b)
            info["po_diff_vs_table"] = int((out != want[ctx.rule]).sum())
            flips[f"stream_b{b}"] = info
        if torch.backends.mps.is_available():
            m = MPSO2(ctx, vids).run()
            w = want[ctx.rule]
            flips["mps"] = dict(
                po_diff=int((m != w).sum()),
                pass_diff=int(((m > 0) != (w > 0)).sum()),
                pairs=int(w.size),
            )
        res[name] = dict(voxels=[int(v) for v in vids], equal=rows, variants=flips)
        print(name, json.dumps(res[name], indent=1), flush=True)
        del ctx
    res["all_equal"] = all(all(r["equal"].values()) for r in res.values() if isinstance(r, dict))
    (C.OUT / "step3_check.json").write_text(json.dumps(res, indent=1) + "\n")


# ---------------------------------------------------------------------------------------------
# time
# ---------------------------------------------------------------------------------------------


def timed(fn: Callable[[], Any]) -> Tuple[float, Any]:
    t0 = time.perf_counter()
    out = fn()
    return time.perf_counter() - t0, out


def cmd_time(a: argparse.Namespace) -> None:
    pf = C.preflight_json(C.OUT / "preflight_step3.json")
    label = C.label_of(pf)
    torch.set_num_threads(1)
    ctx = Ctx("phase_d")
    patches = {n: ctx.patch(n) for n in PATCH_SIZES}
    v2 = patches[2000]
    ctx.bitmaps("u8")
    ctx.bitmaps("packed")
    # table build (untimed in ensure_table) and frame sort, timed here
    t_build, _ = timed(lambda: SK.make_table(ctx.ds, ctx.cap))
    t_sort, _ = timed(
        lambda: K.CTable(
            ctx.tab,
            ctx.z,
            ctx.bitmaps("u8"),
            False,
            ctx.r,
            ctx.rows,
            ctx.cols,
            ctx.rule,
            order=K.frame_order(ctx.tab),
        )
    )
    ref_counts = run_kind(ctx, "o2", v2, 1, "u8")
    ref_pass = (ref_counts > 0).astype(np.uint8)
    n_tests_2000 = int(ctx.tab.hit[:, 0].sum()) * len(v2)

    arms: List[Dict[str, Any]] = []

    def arm(
        name: str,
        kind: str,
        threads: int,
        layout: str = "u8",
        block: int = 288,
        tile: int = 256,
        n: int = 2000,
        reps: int = a.reps,
    ) -> None:
        arms.append(
            dict(
                name=name,
                kind=kind,
                threads=threads,
                layout=layout,
                block=block,
                tile=tile,
                n=n,
                reps=reps,
            )
        )

    # single thread, every ordering
    for layout in ("u8", "packed"):
        arm(f"O1 counts {layout}", "o1", 1, layout)
        arm(f"O1 pass {layout}", "o1p", 1, layout)
        arm(f"O2 counts {layout} b=32FZ", "o2", 1, layout)
        arm(f"O2 pass {layout} b=32FZ", "o2p", 1, layout)
        for t in (256, 2000):
            arm(f"O3 counts {layout} tile={t}", "o3", 1, layout, tile=t)
        for t in (64, 256, 1024):
            arm(f"O5 counts {layout} tile={t}", "o5", 1, layout, tile=t)
    for b in (9, 2304):
        arm(f"O2 counts u8 b={b // 9}FZ", "o2", 1, "u8", block=b)
    # threads
    for th in THREADS[1:]:
        arm(f"O2 counts u8 b=32FZ th={th}", "o2", th)
        arm(f"O2 pass u8 b=32FZ th={th}", "o2p", th)
        arm(f"O1 pass packed th={th}", "o1p", th, "packed")
        arm(f"O3 counts u8 tile=256 th={th}", "o3", th, tile=256)
        arm(f"O5 counts u8 tile=256 th={th}", "o5", th, tile=256)
    for th in (10,):
        arm(f"O2 counts packed b=32FZ th={th}", "o2", th, "packed")
        arm(f"O2 counts u8 b=1FZ th={th}", "o2", th, block=9)
        arm(f"O2 counts u8 b=4FZ th={th}", "o2", th, block=36)
    # linearity
    for n in (500, 1000):
        for name, kind, th, layout in (
            ("O1 pass packed", "o1p", 1, "packed"),
            ("O2 counts u8 b=32FZ", "o2", 1, "u8"),
            ("O3 counts u8 tile=256", "o3", 1, "u8"),
            ("O5 counts u8 tile=256", "o5", 1, "u8"),
            ("O2 counts u8 b=32FZ th=10", "o2", 10, "u8"),
        ):
            arm(f"{name} n={n}", kind, th, layout, n=n)

    results: Dict[str, Dict[str, Any]] = {x["name"]: dict(x, times=[], equal=True) for x in arms}
    t_all = time.time()
    for rep in range(a.reps):
        for x in arms:
            vids = patches[x["n"]]
            dt, out = timed(
                lambda: run_kind(
                    ctx, x["kind"], vids, x["threads"], x["layout"], x["block"], x["tile"]
                )
            )
            if x["n"] == 2000:
                w = ref_pass if x["kind"].endswith("p") else ref_counts
                results[x["name"]]["equal"] &= bool(np.array_equal(out, w))
            results[x["name"]]["times"].append(dt)
            del out
        print(f"rep {rep} done {time.time() - t_all:.0f}s", flush=True)
    # streamed O2 (torch peaks per block) and MPS
    stream = {}
    for b in (1, 32, 256):
        ts = []
        for _ in range(a.reps):
            dt, (out, info) = timed(lambda: stream_o2(ctx, v2, b, threads=1))
            ts.append(dt)
            info["equal_table_kernel"] = bool(np.array_equal(out, ref_counts))
        info["times"] = ts
        stream[f"b={b}FZ"] = info
        print("stream", b, info, flush=True)
    mps: Dict[str, Any] = {}
    if torch.backends.mps.is_available():
        for blk in (288, 576, 1152):
            m = MPSO2(ctx, v2, block=blk)
            m.run()  # warm-up
            ts = []
            for _ in range(a.reps):
                dt, out = timed(m.run)
                ts.append(dt)
            mps[f"block={blk}"] = dict(
                times=ts,
                po_diff=int((out != ref_counts).sum()),
                pass_diff=int(((out > 0) != (ref_counts > 0)).sum()),
            )
            print("mps", blk, mps[f"block={blk}"], flush=True)
            del m
    for x in results.values():
        x["median_s"] = float(np.median(x["times"]))
        n_tests = int(ctx.tab.hit[:, 0].sum()) * x["n"]
        x["ns_per_test"] = 1e9 * x["median_s"] / n_tests
        x["pairs_per_s"] = ctx.tab.n_orient * x["n"] / x["median_s"]
        x["projected_full_layer_s"] = x["median_s"] * len(ctx.X) / x["n"]
    res = dict(
        label=label,
        preflight_load=pf.get("loadavg"),
        preflight_end=C.preflight_json(C.OUT / "preflight_step3_end.json").get("loadavg"),
        reps=a.reps,
        n_tests_2000=n_tests_2000,
        n_orient=ctx.tab.n_orient,
        n_voxels_layer=int(len(ctx.X)),
        table_build_s=t_build,
        frame_sort_s=t_sort,
        table_bytes=int(
            ctx.tab.A.nbytes + ctx.tab.c.nbytes + ctx.tab.hit.nbytes + 12 * ctx.tab.n_peaks
        ),
        bitmap_bytes=dict(
            u8_det0=int(ctx.bitmaps("u8")[0].nbytes), packed_both=int(ctx.packed_all.nbytes)
        ),
        out_bytes_2000=int(ctx.tab.n_orient * 2000),
        peak_rss_gb=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e9,
        arms=results,
        stream=stream,
        mps=mps,
        passers_patch=int(ref_pass.sum()),
    )
    (C.OUT / "step3_timing.json").write_text(json.dumps(res, indent=1) + "\n")
    print(
        json.dumps(
            {
                k: (v["median_s"], round(v["ns_per_test"], 3), v["equal"])
                for k, v in results.items()
            },
            indent=1,
        )
    )


# ---------------------------------------------------------------------------------------------
# full layer
# ---------------------------------------------------------------------------------------------


def cmd_full(a: argparse.Namespace) -> None:
    pf = C.preflight_json(C.OUT / "preflight_step3_full.json")
    label = C.label_of(pf)
    ctx = Ctx("phase_d")
    ctx.bitmaps(a.layout)
    ct_ = ctx.ctab(a.layout, by_frame=(a.kind == "o3"))
    V = len(ctx.X)
    tile = a.voxel_tile
    n_pass = 0
    per_voxel = np.zeros(V, np.int64)
    t_kernel = t_reduce = 0.0
    t0 = time.perf_counter()
    for s in range(0, V, tile):
        vids = np.arange(s, min(s + tile, V))
        t1 = time.perf_counter()
        out = K.run(
            a.kind,
            ct_,
            ctx.X[vids, 0],
            ctx.X[vids, 1],
            threads=a.threads,
            block=a.block,
            tile=a.tile,
        )
        t2 = time.perf_counter()
        pm = out > 0
        per_voxel[vids] = pm.sum(axis=0)
        n_pass += int(pm.sum())
        t3 = time.perf_counter()
        t_kernel += t2 - t1
        t_reduce += t3 - t2
        del out, pm
    wall = time.perf_counter() - t0
    # check a few voxels against the Step 0 reference
    checks = {}
    for vid in C.reference_voxels(ctx.ds)[:5]:
        z = np.load(C.CACHE / "ref" / f"phase_d_{vid}.npz")
        checks[int(vid)] = dict(ref_passers=int(len(z["passers"])), passers=int(per_voxel[vid]))
    res = dict(
        label=label,
        preflight_load=pf.get("loadavg"),
        kind=a.kind,
        layout=a.layout,
        threads=a.threads,
        block=a.block,
        tile=a.tile,
        voxel_tile=tile,
        n_voxels=V,
        n_orient=ctx.tab.n_orient,
        wall_s=wall,
        kernel_s=t_kernel,
        reduce_s=t_reduce,
        passers_total=n_pass,
        passers_per_voxel_median=float(np.median(per_voxel)),
        passers_per_voxel_min=int(per_voxel.min()),
        passers_per_voxel_max=int(per_voxel.max()),
        sparse_list_bytes_int32_pairs=int(8 * n_pass),
        pass_bitset_bytes=int(ctx.tab.n_orient * V / 8),
        counts_matrix_bytes_uint8=int(ctx.tab.n_orient * V),
        reference_check=checks,
        peak_rss_gb=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e9,
        preflight_end=C.preflight_json(C.OUT / "preflight_step3_full_end.json").get("loadavg"),
    )
    np.save(C.CACHE / "full_layer_passers_per_voxel.npy", per_voxel)
    (C.OUT / "step3_full_layer.json").write_text(json.dumps(res, indent=1) + "\n")
    print(json.dumps(res, indent=1))


# ---------------------------------------------------------------------------------------------
# cost model
# ---------------------------------------------------------------------------------------------


def cmd_model(a: argparse.Namespace) -> None:
    """Inputs measured from the data (2,000-voxel patch): per peak the distinct pixels the patch
    voxels land on, the bright ones among them and their bounding box; per test the bright rate.
    Per-operation costs are taken from the Step 3 timing (O2 ns per test) and from a random-gather
    micro-measure over the bitmap working set (in step3_timing.json when present)."""
    ctx = Ctx("phase_d")
    v2 = ctx.patch(2000)
    X = ctx.X[v2]
    rng = np.random.default_rng(5)
    hit0 = np.nonzero(ctx.tab.hit[:, 0])[0]
    sample = rng.choice(hit0, a.n_peaks, replace=False)
    U = np.einsum("paj,vj->pva", ctx.tab.A[sample, 0], X) + ctx.tab.c[sample, 0][:, None, :]
    cx, cy = C.trunc_index(U[..., 0]), C.trunc_index(U[..., 1])
    distinct, bright_d, bbox, tests_bright, inb_tests = [], [], [], [], []
    for i, p in enumerate(sample):
        inb, br = ctx.bm.lookup(np.full(len(v2), ctx.tab.frame[p]), 0, cx[i], cy[i])
        key = cx[i][inb] * 100000 + cy[i][inb]
        u, first = np.unique(key, return_index=True)
        distinct.append(len(u))
        bright_d.append(int(br[inb][first].sum()))
        tests_bright.append(int(br.sum()))
        inb_tests.append(int(inb.sum()))
        if inb.any():
            bbox.append(int((np.ptp(cx[i][inb]) + 1) * (np.ptp(cy[i][inb]) + 1)))
        else:
            bbox.append(0)
    distinct_a, bright_a, bbox_a = np.array(distinct), np.array(bright_d), np.array(bbox)
    V = len(v2)
    tim = (
        json.loads((C.OUT / "step3_timing.json").read_text())
        if (C.OUT / "step3_timing.json").exists()
        else None
    )
    t_c = None
    if tim:
        t_c = tim["arms"]["O2 counts u8 b=32FZ"]["ns_per_test"]
    res = dict(
        n_peaks_sampled=int(a.n_peaks),
        patch=V,
        tests_per_peak=V,
        bright_rate_per_test=float(np.sum(tests_bright) / (V * a.n_peaks)),
        in_bounds_rate_per_test=float(np.sum(inb_tests) / (V * a.n_peaks)),
        distinct_pixels_per_peak_median=float(np.median(distinct_a)),
        voxels_per_distinct_pixel=float(np.sum(inb_tests) / max(np.sum(distinct_a), 1)),
        bright_distinct_pixels_per_peak_mean=float(bright_a.mean()),
        bbox_pixels_per_peak_median=float(np.median(bbox_a)),
        o2_ns_per_test=t_c,
    )
    # O4 model: per peak, scan the footprint bounding box of the packed dilated bitmap (64-bit
    # words) and, per bright pixel, back-project and enumerate its voxels
    for t_bp in (5.0, 10.0, 20.0):
        t_word = 0.5
        o4_ns = np.mean(bbox_a / 64.0 * t_word + bright_a * t_bp)
        res[f"o4_ns_per_peak_tbp{int(t_bp)}"] = float(o4_ns)
        if t_c:
            res[f"o2_over_o4_tbp{int(t_bp)}"] = float(V * t_c / o4_ns)
    (C.OUT / "step3_model.json").write_text(json.dumps(res, indent=1) + "\n")
    print(json.dumps(res, indent=1))


def cmd_time_o4(a: argparse.Namespace) -> None:
    """O4 (built because the model of ``model`` predicted >= 2x at its central cost assumption)
    timed against O2 on the same patch, same bitmaps (packed), interleaved repeats."""
    pf = C.preflight_json(C.OUT / "preflight_step3_o4.json")
    label = C.label_of(pf)
    ctx = Ctx("phase_d")
    up = np.array([bool(v.points_up) for v in ctx.ds.voxels])
    v2 = ctx.patch(2000)
    X = ctx.X[v2]
    lat = K.Lattice(X, up[v2], float(ctx.ds.voxels[0].side_length))
    ct_ = ctx.ctab("packed")
    arms = {
        "O2 counts packed th=1": lambda: (K.run("o2", ct_, X[:, 0], X[:, 1], threads=1), None),
        "O4 counts packed th=1": lambda: K.run_o4(ct_, lat, X[:, 0], X[:, 1], threads=1),
        "O2 counts packed th=10": lambda: (K.run("o2", ct_, X[:, 0], X[:, 1], threads=10), None),
        "O4 counts packed th=10": lambda: K.run_o4(ct_, lat, X[:, 0], X[:, 1], threads=10),
    }
    ref = K.run("o2", ct_, X[:, 0], X[:, 1], threads=4)
    res: Dict[str, Dict[str, Any]] = {k: dict(times=[], equal=True) for k in arms}
    for _ in range(a.reps):
        for k, fn in arms.items():
            dt, (out, st) = timed(fn)
            res[k]["times"].append(dt)
            res[k]["equal"] &= bool(np.array_equal(out, ref))
            if st:
                res[k]["stats"] = st
    for k, v in res.items():
        v["median_s"] = float(np.median(v["times"]))
    st = res["O4 counts packed th=1"]["stats"]
    n_hit = int(ctx.tab.hit[:, 0].sum())
    out = dict(
        label=label,
        preflight_load=pf.get("loadavg"),
        reps=a.reps,
        arms=res,
        bright_pixels_scanned=st["bright_pixels"],
        bright_pixels_per_peak=st["bright_pixels"] / n_hit,
        candidates=st["candidates"],
        candidates_per_bright_pixel=st["candidates"] / max(st["bright_pixels"], 1),
        o2_tests=n_hit * len(v2),
        ns_per_bright_pixel_1thread=1e9
        * res["O4 counts packed th=1"]["median_s"]
        / st["bright_pixels"],
        o4_over_o2_1thread=res["O4 counts packed th=1"]["median_s"]
        / res["O2 counts packed th=1"]["median_s"],
        preflight_end=C.preflight_json(C.OUT / "preflight_step3_o4_end.json").get("loadavg"),
    )
    (C.OUT / "step3_o4.json").write_text(json.dumps(out, indent=1) + "\n")
    print(json.dumps(out, indent=1))


# ---------------------------------------------------------------------------------------------
# rules: pass rates under the default and the C++ qualified-peak rule (no timing)
# ---------------------------------------------------------------------------------------------


def pass_split(p0: np.ndarray, p1: np.ndarray) -> Dict[str, Any]:
    """Paired pass counts of two rules on the same (orientation, voxel) pairs; p0, p1 bool."""
    n = int(p0.size)
    per0, per1 = p0.sum(axis=0), p1.sum(axis=0)
    return dict(
        pairs=n,
        voxels=int(p0.shape[1]),
        pass_rule0=int(p0.sum()),
        pass_rule1=int(p1.sum()),
        rate_rule0=float(p0.sum() / n),
        rate_rule1=float(p1.sum() / n),
        rule0_only=int((p0 & ~p1).sum()),
        rule1_only=int((p1 & ~p0).sum()),
        both=int((p0 & p1).sum()),
        per_voxel_rate_rule0_min=float(per0.min() / p0.shape[0]),
        per_voxel_rate_rule0_max=float(per0.max() / p0.shape[0]),
        per_voxel_rate_rule1_min=float(per1.min() / p1.shape[0]),
        per_voxel_rate_rule1_max=float(per1.max() / p1.shape[0]),
        per_voxel_passers_rule1_median=float(np.median(per1)),
    )


def cmd_rules(a: argparse.Namespace) -> None:
    """Screen pass rates under rule 0 (the config's default qualified-peak rule) and rule 1
    (``CppOverlapRule``: bright on detector 0 and, on detector 1, bright or not lit) on the 20
    reference voxels (also checked against the numpy kernel), the 2,000-voxel patch and the whole
    layer. Counts only: nothing here is timed."""
    ctx = Ctx("phase_d")
    assert ctx.rule == 0, "the config is expected to use the default rule"
    cts = {r: ctx.ctab("u8", rule=r) for r in (0, 1)}
    res: Dict[str, Any] = dict(n_orient=int(ctx.tab.n_orient))

    def passes(vids: np.ndarray, rule: int) -> np.ndarray:
        X = ctx.X[vids]
        return K.run("o2p", cts[rule], X[:, 0], X[:, 1], threads=a.threads).astype(bool)

    ref = np.array(C.reference_voxels(ctx.ds))
    p = {r: passes(ref, r) for r in (0, 1)}
    want = {r: C.kernel_numpy(ctx.tab, ctx.X[ref], ctx.bm, bool(r))["po"] > 0 for r in (0, 1)}
    res["reference"] = pass_split(p[0], p[1])
    res["reference"]["equal_numpy"] = bool(all(np.array_equal(p[r], want[r]) for r in (0, 1)))
    res["reference"]["voxel_ids"] = [int(v) for v in ref]
    res["reference"]["per_voxel"] = [
        dict(vid=int(v), rule0=int(p[0][:, j].sum()), rule1=int(p[1][:, j].sum()))
        for j, v in enumerate(ref)
    ]
    print("reference", json.dumps({k: v for k, v in res["reference"].items() if k != "per_voxel"}))
    del want
    v2 = ctx.patch(2000)
    res["patch_2000"] = pass_split(passes(v2, 0), passes(v2, 1))
    print("patch", json.dumps(res["patch_2000"]), flush=True)
    V = len(ctx.X)
    per = {r: np.zeros(V, np.int64) for r in (0, 1)}
    tot = dict(pairs=0, pass_rule0=0, pass_rule1=0, rule0_only=0, rule1_only=0, both=0)
    for s in range(0, V, a.voxel_tile):
        vids = np.arange(s, min(s + a.voxel_tile, V))
        q0, q1 = passes(vids, 0), passes(vids, 1)
        per[0][vids], per[1][vids] = q0.sum(axis=0), q1.sum(axis=0)
        tot["pairs"] += int(q0.size)
        tot["pass_rule0"] += int(q0.sum())
        tot["pass_rule1"] += int(q1.sum())
        tot["rule0_only"] += int((q0 & ~q1).sum())
        tot["rule1_only"] += int((q1 & ~q0).sum())
        tot["both"] += int((q0 & q1).sum())
        del q0, q1
    n_o = ctx.tab.n_orient
    res["layer"] = dict(
        tot,
        voxels=V,
        rate_rule0=tot["pass_rule0"] / tot["pairs"],
        rate_rule1=tot["pass_rule1"] / tot["pairs"],
        per_voxel_passers_rule0_median=float(np.median(per[0])),
        per_voxel_passers_rule1_median=float(np.median(per[1])),
        per_voxel_passers_rule1_min=int(per[1].min()),
        per_voxel_passers_rule1_max=int(per[1].max()),
        per_voxel_rate_rule1_min=float(per[1].min() / n_o),
        per_voxel_rate_rule1_max=float(per[1].max() / n_o),
    )
    print("layer", json.dumps(res["layer"]), flush=True)
    (C.OUT / "step3_rule_rates.json").write_text(json.dumps(res, indent=1) + "\n")


# ---------------------------------------------------------------------------------------------
# time_o3: O3 with tile 2,000 on threads, interleaved with O2 (review follow-up)
# ---------------------------------------------------------------------------------------------


def cmd_time_o3(a: argparse.Namespace) -> None:
    """O3 tile 2,000 (the fastest one-thread arm of ``time``) and tile 256 against O2 at 1 and
    10 threads on the 2,000-voxel patch, interleaved repeats in one process."""
    pf = C.preflight_json(C.OUT / "preflight_step3_o3.json")
    label = C.label_of(pf)
    torch.set_num_threads(1)
    ctx = Ctx("phase_d")
    v2 = ctx.patch(2000)
    ctx.bitmaps("u8")
    ref = run_kind(ctx, "o2", v2, 4, "u8")
    arms: Dict[str, Tuple[str, int, int]] = {}
    for th in (1, 10):
        sfx = "" if th == 1 else f" th={th}"
        arms[f"O2 counts u8 b=32FZ{sfx}"] = ("o2", th, 256)
        arms[f"O3 counts u8 tile=2000{sfx}"] = ("o3", th, 2000)
        arms[f"O3 counts u8 tile=256{sfx}"] = ("o3", th, 256)
    res: Dict[str, Dict[str, Any]] = {k: dict(times=[], equal=True) for k in arms}
    for rep in range(a.reps):
        for k, (kind, th, tile) in arms.items():
            dt, out = timed(lambda: run_kind(ctx, kind, v2, th, "u8", tile=tile))
            res[k]["times"].append(dt)
            res[k]["equal"] &= bool(np.array_equal(out, ref))
            del out
        print(f"rep {rep} done", flush=True)
    for v in res.values():
        v["median_s"] = float(np.median(v["times"]))
    out = dict(
        label=label,
        preflight_load=pf.get("loadavg"),
        reps=a.reps,
        n_voxels=int(len(v2)),
        arms=res,
        preflight_end=C.preflight_json(C.OUT / "preflight_step3_o3_end.json").get("loadavg"),
    )
    (C.OUT / "step3_o3_threads.json").write_text(json.dumps(out, indent=1) + "\n")
    print(json.dumps({k: (v["median_s"], v["equal"]) for k, v in res.items()}, indent=1))


def main() -> None:
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("check")
    p = sub.add_parser("rules")
    p.add_argument("--threads", type=int, default=8)
    p.add_argument("--voxel-tile", type=int, default=4096)
    p = sub.add_parser("time_o3")
    p.add_argument("--reps", type=int, default=5)
    p = sub.add_parser("time_o4")
    p.add_argument("--reps", type=int, default=3)
    p = sub.add_parser("time")
    p.add_argument("--reps", type=int, default=3)
    p = sub.add_parser("full")
    p.add_argument("--kind", default="o2")
    p.add_argument("--layout", default="u8")
    p.add_argument("--threads", type=int, default=10)
    p.add_argument("--block", type=int, default=288)
    p.add_argument("--tile", type=int, default=256)
    p.add_argument("--voxel-tile", type=int, default=4096)
    p = sub.add_parser("model")
    p.add_argument("--n-peaks", type=int, default=3000)
    a = ap.parse_args()
    dict(
        check=cmd_check,
        rules=cmd_rules,
        time=cmd_time,
        time_o3=cmd_time_o3,
        full=cmd_full,
        model=cmd_model,
        time_o4=cmd_time_o4,
    )[a.cmd](a)


if __name__ == "__main__":
    main()
