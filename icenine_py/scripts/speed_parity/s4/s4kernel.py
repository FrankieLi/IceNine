"""ctypes wrapper of the S4 scratch C kernels (s4kernel.c), compiled with the system cc.

``CTable`` holds the per-peak arrays of one layer (z folded into c), the bitmaps and the ctypes
``Table`` struct; ``run`` executes one ordering over orientation and voxel ranges with Python
threads (the ctypes call releases the GIL, so the threads run the C code in parallel).
"""

import ctypes as ct
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import common as C  # noqa: E402

SRC = HERE / "s4kernel.c"
LIB = C.CACHE / "s4kernel.so"
CFLAGS = ["-O3", "-ffp-contract=off", "-shared", "-fPIC"]


def compile_lib(force: bool = False) -> Path:
    LIB.parent.mkdir(parents=True, exist_ok=True)
    if force or not LIB.exists() or LIB.stat().st_mtime < SRC.stat().st_mtime:
        subprocess.run(["cc", *CFLAGS, "-o", str(LIB), str(SRC)], check=True)
    return LIB


P_D = ct.POINTER(ct.c_double)
P_U8 = ct.POINTER(ct.c_uint8)
P_I32 = ct.POINTER(ct.c_int32)
P_I64 = ct.POINTER(ct.c_int64)


class _Det(ct.Structure):
    _fields_ = [
        ("a0", P_D),
        ("a1", P_D),
        ("c", P_D),
        ("b0", P_D),
        ("b1", P_D),
        ("d", P_D),
        ("hit", P_U8),
        ("bm", P_U8),
    ]


class _Table(ct.Structure):
    _fields_ = [
        ("det", _Det * 2),
        ("frame", P_I32),
        ("frame_stride", ct.c_int64),
        ("row_stride", ct.c_int32),
        ("rows", ct.c_int32),
        ("cols", ct.c_int32),
        ("r", ct.c_int32),
        ("packed", ct.c_int32),
        ("rule", ct.c_int32),
    ]


def _ptr(a: np.ndarray, typ: Any) -> Any:
    assert a.flags["C_CONTIGUOUS"]
    return a.ctypes.data_as(typ)


_lib: Optional[ct.CDLL] = None


def lib() -> ct.CDLL:
    global _lib
    if _lib is None:
        L = ct.CDLL(str(compile_lib()))
        args_counts = [ct.c_void_p, P_I64, ct.c_int64, ct.c_int64, P_D, P_D, ct.c_int64,
                       ct.c_int64, P_U8, ct.c_int64]  # fmt: skip
        for name in ("o2_counts", "o1_counts", "o1_pass"):
            getattr(L, name).argtypes = args_counts
            getattr(L, name).restype = None
        L.o2_pass.argtypes = args_counts + [P_I32]
        L.o2_pass.restype = None
        L.o5_counts.argtypes = args_counts[:8] + [ct.c_int64, P_U8, ct.c_int64]
        L.o5_counts.restype = None
        L.o3_counts.argtypes = [ct.c_void_p, P_I64, ct.c_int64, ct.c_int64, P_I32, ct.c_int64,
                                P_D, P_D, ct.c_int64, ct.c_int64, ct.c_int64, P_U8,
                                ct.c_int64]  # fmt: skip
        L.o3_counts.restype = None
        L.table_size.restype = ct.c_int64
        assert L.table_size() == ct.sizeof(_Table)
        _lib = L
    return _lib


class CTable:
    """Arrays and struct of the kernel for one peak order (by orientation, or by frame)."""

    def __init__(
        self,
        tab: C.PeakTable,
        z: float,
        bitmaps: List[np.ndarray],
        packed: bool,
        r: int,
        rows: int,
        cols: int,
        rule: int,
        order: Optional[np.ndarray] = None,
    ) -> None:
        idx = np.arange(tab.n_peaks) if order is None else order
        self.keep: Dict[str, np.ndarray] = {}
        dets = []
        for d in range(2):
            A, c = tab.A[idx, d], tab.c[idx, d]
            # u = A0 x + A1 y + (A2 z + c); for z = 0 the folded constant is c exactly
            cz = c + A[:, :, 2] * z if z != 0.0 else c.copy()
            arrs = dict(
                a0=np.ascontiguousarray(A[:, 0, 0]),
                a1=np.ascontiguousarray(A[:, 0, 1]),
                c=np.ascontiguousarray(cz[:, 0]),
                b0=np.ascontiguousarray(A[:, 1, 0]),
                b1=np.ascontiguousarray(A[:, 1, 1]),
                d=np.ascontiguousarray(cz[:, 1]),
                hit=np.ascontiguousarray(tab.hit[idx, d].astype(np.uint8)),
            )
            for k, v in arrs.items():
                self.keep[f"{k}{d}"] = v
            bm = bitmaps[d] if d < len(bitmaps) else bitmaps[0]
            self.keep[f"bm{d}"] = bm
            dets.append(
                _Det(
                    _ptr(arrs["a0"], P_D),
                    _ptr(arrs["a1"], P_D),
                    _ptr(arrs["c"], P_D),
                    _ptr(arrs["b0"], P_D),
                    _ptr(arrs["b1"], P_D),
                    _ptr(arrs["d"], P_D),
                    _ptr(arrs["hit"], P_U8),
                    _ptr(bm, P_U8),
                )
            )
        self.frame = np.ascontiguousarray(tab.frame[idx].astype(np.int32))
        self.seg = np.ascontiguousarray(tab.seg[idx].astype(np.int32))
        bm0 = bitmaps[0]
        self.struct = _Table(
            (_Det * 2)(*dets),
            _ptr(self.frame, P_I32),
            int(bm0.shape[1] * bm0.shape[2]),
            int(bm0.shape[2]),
            rows,
            cols,
            r,
            1 if packed else 0,
            rule,
        )
        self.n_orient = tab.n_orient
        self.pstart = np.ascontiguousarray(
            np.searchsorted(tab.seg[idx], np.arange(tab.n_orient + 1)).astype(np.int64)
        )
        n_frames = int(bm0.shape[0])
        self.fstart = np.ascontiguousarray(
            np.searchsorted(self.frame, np.arange(n_frames + 1)).astype(np.int64)
        )
        self.n_frames = n_frames
        self.n_peaks = len(idx)

    @property
    def ref(self) -> Any:
        return ct.byref(self.struct)


def frame_order(tab: C.PeakTable) -> np.ndarray:
    return np.argsort(tab.frame, kind="stable")


def unpack_bitmaps(packed: np.ndarray, det: int) -> np.ndarray:
    """(n_frames, rows + 2r, 8 * stride) uint8 one byte per pixel from the packed bitmaps."""
    return np.ascontiguousarray(np.unpackbits(packed[:, det], axis=-1))


def run(
    kind: str,
    ctab: CTable,
    xs: np.ndarray,
    ys: np.ndarray,
    threads: int = 1,
    block: int = 288,
    tile: int = 256,
    o_range: Optional[Tuple[int, int]] = None,
) -> np.ndarray:
    """Run ordering ``kind`` (o1, o1p, o2, o2p, o3, o5) over all orientations (or o_range) and
    the voxels xs, ys: (n_o, V) uint8 counts (pass flags for the *p kinds). Work units: blocks of
    ``block`` orientations (o1: voxel chunks; o3: voxel tiles) shared by ``threads`` threads."""
    L = lib()
    V = len(xs)
    xs = np.ascontiguousarray(xs, np.float64)
    ys = np.ascontiguousarray(ys, np.float64)
    o0, o1 = o_range or (0, ctab.n_orient)
    n_o = o1 - o0
    out = np.zeros((n_o, V), np.uint8)
    xp, yp = _ptr(xs, P_D), _ptr(ys, P_D)
    ps = _ptr(ctab.pstart, P_I64)
    T = ctab.ref

    def out_at(o: int) -> Any:
        return ct.cast(out.ctypes.data + (o - o0) * V, P_U8)

    jobs: List[Any] = []
    if kind in ("o2", "o2p", "o5"):
        for s in range(o0, o1, block):
            e = min(s + block, o1)
            jobs.append((s, e))

        def work(j: Tuple[int, int]) -> None:
            s, e = j
            if kind == "o2":
                L.o2_counts(T, ps, s, e, xp, yp, 0, V, out_at(s), V)
            elif kind == "o5":
                L.o5_counts(T, ps, s, e, xp, yp, 0, V, tile, out_at(s), V)
            else:
                alive = np.empty(V, np.int32)
                L.o2_pass(T, ps, s, e, xp, yp, 0, V, out_at(s), V, _ptr(alive, P_I32))

    elif kind in ("o1", "o1p"):
        step = max(1, min(tile, (V + threads - 1) // threads))
        for s in range(0, V, step):
            jobs.append((s, min(s + step, V)))

        def work(j: Tuple[int, int]) -> None:
            s, e = j
            fn = L.o1_counts if kind == "o1" else L.o1_pass
            fn(T, ps, o0, o1, xp, yp, s, e, ct.cast(out.ctypes.data + s, P_U8), V)

    elif kind == "o3":
        assert o_range is None
        fs = _ptr(ctab.fstart, P_I64)
        sg = _ptr(ctab.seg, P_I32)
        for s in range(0, V, tile):
            jobs.append((s, min(s + tile, V)))

        def work(j: Tuple[int, int]) -> None:
            s, e = j
            L.o3_counts(T, fs, 0, ctab.n_frames, sg, 0, xp, yp, s, e, e - s,
                        ct.cast(out.ctypes.data + s, P_U8), V)  # fmt: skip

    else:
        raise ValueError(kind)
    if threads == 1:
        for j in jobs:
            work(j)
    else:
        with ThreadPoolExecutor(threads) as ex:
            list(ex.map(work, jobs))
    return out


class _Lattice(ct.Structure):
    _fields_ = [
        ("grid", P_I32 * 2),
        ("i0", ct.c_int32 * 2),
        ("j0", ct.c_int32 * 2),
        ("ni", ct.c_int32 * 2),
        ("nj", ct.c_int32 * 2),
        ("org", (ct.c_double * 2) * 2),
        ("e1", ct.c_double * 2),
        ("e2", ct.c_double * 2),
        ("bx0", ct.c_double),
        ("bx1", ct.c_double),
        ("by0", ct.c_double),
        ("by1", ct.c_double),
        ("slack", ct.c_double),
    ]


class Lattice:
    """The two centroid lattices (one per triangle type) of the patch, for O4."""

    def __init__(self, X: np.ndarray, up: np.ndarray, side: float, slack: float = 0.01) -> None:
        e1 = np.array([side, 0.0])
        e2 = np.array([0.5 * side, 0.5 * side * np.sqrt(3.0)])
        G = np.stack([e1, e2], axis=1)
        self.keep: List[np.ndarray] = []
        grids, i0, j0, ni, nj, org = [], [], [], [], [], []
        for t in (True, False):
            sel = np.nonzero(up == t)[0]
            o = X[sel[0], :2]
            ij = np.linalg.solve(G, (X[sel, :2] - o).T).T
            assert np.abs(ij - np.round(ij)).max() < 1e-3, "centroids are not on a lattice"
            ij = np.round(ij).astype(np.int64)
            lo, hi = ij.min(axis=0), ij.max(axis=0)
            g = np.full((hi[0] - lo[0] + 1, hi[1] - lo[1] + 1), -1, np.int32)
            g[ij[:, 0] - lo[0], ij[:, 1] - lo[1]] = sel.astype(np.int32)
            g = np.ascontiguousarray(g)
            self.keep.append(g)
            grids.append(_ptr(g, P_I32))
            i0.append(int(lo[0]))
            j0.append(int(lo[1]))
            ni.append(int(g.shape[0]))
            nj.append(int(g.shape[1]))
            org.append((ct.c_double * 2)(*o))
        bx0, by0 = X[:, :2].min(axis=0)
        bx1, by1 = X[:, :2].max(axis=0)
        self.struct = _Lattice(
            (P_I32 * 2)(*grids),
            (ct.c_int32 * 2)(*i0),
            (ct.c_int32 * 2)(*j0),
            (ct.c_int32 * 2)(*ni),
            (ct.c_int32 * 2)(*nj),
            ((ct.c_double * 2) * 2)(*org),
            (ct.c_double * 2)(*e1),
            (ct.c_double * 2)(*e2),
            float(bx0),
            float(bx1),
            float(by0),
            float(by1),
            slack,
        )


def run_o4(
    ctab: CTable, lat: Lattice, xs: np.ndarray, ys: np.ndarray, threads: int = 1, block: int = 288
) -> Tuple[np.ndarray, Dict[str, int]]:
    """O4 counts (n_orient, V) uint8 over all orientations (packed bitmaps, rule 0)."""
    L = lib()
    if not hasattr(L.o4_counts, "_set"):
        L.o4_counts.argtypes = [ct.c_void_p, ct.c_void_p, P_I64, ct.c_int64, ct.c_int64, P_D,
                                P_D, P_U8, ct.c_int64, P_I64]  # fmt: skip
        L.o4_counts.restype = None
        L.o4_counts._set = True
    assert ctab.struct.packed == 1 and ctab.struct.rule == 0
    V = len(xs)
    xs = np.ascontiguousarray(xs, np.float64)
    ys = np.ascontiguousarray(ys, np.float64)
    out = np.zeros((ctab.n_orient, V), np.uint8)
    jobs = [(s, min(s + block, ctab.n_orient)) for s in range(0, ctab.n_orient, block)]
    stats = np.zeros((len(jobs), 2), np.int64)

    def work(k: int) -> None:
        s, e = jobs[k]
        L.o4_counts(ctab.ref, ct.byref(lat.struct), _ptr(ctab.pstart, P_I64), s, e,
                    _ptr(xs, P_D), _ptr(ys, P_D), ct.cast(out.ctypes.data + s * V, P_U8), V,
                    ct.cast(stats.ctypes.data + 16 * k, P_I64))  # fmt: skip

    if threads == 1:
        for k in range(len(jobs)):
            work(k)
    else:
        with ThreadPoolExecutor(threads) as ex:
            list(ex.map(work, range(len(jobs))))
    return out, dict(bright_pixels=int(stats[:, 0].sum()), candidates=int(stats[:, 1].sum()))
