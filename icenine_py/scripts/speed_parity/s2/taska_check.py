"""S2 Task A check: the Step 0 fused prototype against C++ on Task A's 2,230 orientations.

Plan ``docs/todo_s2_fused_c_eval.md``, section "Task A check with the prototype" (D-20261027-I4),
Parity item 2. Scripts only; no library change. Three subcommands:

``inputs``  checks that ``scripts/eval_reduction/parity_inputs.py`` regenerates the stored inputs.
    The C++ driver output and the input files of Task A were not kept; what is stored is
    ``benchmarks/eval_reduction/A/step0/results.json.gz`` (per orientation: labels, the FZ / grid
    indices of the pairs, whether the float64 product cast to float32 equals the C++ float32 product
    and their largest difference, and the C++ and Python aggregates of every arm). The check
    compares the regenerated labels and indices with the stored ones and finds the float32
    operation order that reproduces the C++ product ``oLocal * oCenter`` (``SMatrix3x3::operator*``,
    ``XDM++/libXDM/3dMath.cpp:410-430``, Apple clang -O3, FP contraction on) from the stored
    ``product_bit_equal`` and ``product_max_abs_diff`` of every pair record.

``run``     evaluates every record of the given voxels on the C++ product orientation: the library
    ``evaluate`` (torch) on Task A's four Python arms (screen and complete cost, ``CppEdgeRules``
    off and on: a regeneration check against the stored Python aggregates), and on the two
    ``CppEdgeRules`` arms also ``torch_mirror.mirror`` (must equal the library), the float64
    reference and the fused prototype (``s2fused.FusedEval``, variants ``portable`` and
    ``accelerate``). Every orientation whose fused ``OverlapInfo`` differs from torch is traced
    with ``s2_parity.trace``. Writes one JSON per call. S2 Step 4 (D-20261027-I7): ``--variants
    library`` runs the library kernel instead (``VoxelCostFunction`` with the key FusedCost, the
    shipped g_lab order, intermediates from ``fused_debug``), and ``--every-coord`` traces every
    flipped pixel coordinate, not only the worst per row and detector.

``outlier`` traces the 1.50 px float64-reference outlier of Step 0 (Phase D voxel 18558, mode
    ``edge``, S1's 1,000 orientations): every stage-D pixel coordinate whose torch or fused value
    is more than ``--min-px`` from the float64 reference, with the detector coordinate in float32
    (torch) and float64.

  uv run python scripts/speed_parity/s2/taska_check.py inputs --dir DIR --out OUT.json
  uv run python scripts/speed_parity/s2/taska_check.py run --dir DIR --voxels 15901 --out OUT.json
  uv run python scripts/speed_parity/s2/taska_check.py outlier --out OUT.json

``DIR`` is the output directory of ``parity_inputs.py --out DIR``. Paths are made absolute first
(the Task A context changes the working directory to ``Examples/Example2.ManyGrains``).
"""

import argparse
import contextlib
import ctypes
import ctypes.util
import gzip
import io
import json
import sys
import time
from collections import Counter
from pathlib import Path
from typing import Any, Callable, Dict, List, Tuple

import numpy as np
import torch

HERE = Path(__file__).resolve().parent
ICE = HERE.parents[2]
sys.path.insert(0, str(ICE))
sys.path.insert(0, str(HERE.parent))
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(ICE / "scripts" / "eval_reduction"))

import parity_common as PC  # noqa: E402
import s1_parity as S1  # noqa: E402
import s2_parity as S2P  # noqa: E402
import s2fused as SF  # noqa: E402
import torch_mirror as TM  # noqa: E402

STORED = ICE / "benchmarks" / "eval_reduction" / "A" / "step0" / "results.json.gz"
META_KEYS = ("deg", "fz", "j", "why", "cpp", "py")
#: Task A arm names -> (Python cost function arm, CppEdgeRules)
ARMS: Dict[str, Tuple[str, bool]] = {
    "EXS": ("screen", False),
    "EXF": ("full", False),
    "EXS_edge": ("screen", True),
    "EXF_edge": ("full", True),
}
FUSED_ARMS = ("EXS_edge", "EXF_edge")
AGG = ("pixov", "pixon", "pkov", "pkon", "quality", "conf")
PIXEL_EDGE_PX = 0.01  # parity rule D-20261027-I2


def load_stored() -> Dict[Tuple[int, int], Dict[str, Any]]:
    d = json.loads(gzip.open(STORED).read())
    return {(int(r["voxel"]), int(r["id"])): r for r in d["rows"]}


def load_inputs(d: Path, voxel: int) -> List[Dict[str, Any]]:
    return list(json.loads((d / f"in_{voxel}.json").read_text()))


def git_meta() -> Dict[str, Any]:
    return dict(
        head=S1.git_head(),
        dirty=S1.git_dirty(),
        dirty_files=S2P.dirty_files(),
        torch=torch.__version__,
        numpy=np.__version__,
    )


# ----------------------------------------------------------------------------------------------
# the C++ float32 product
# ----------------------------------------------------------------------------------------------
_libm = ctypes.CDLL(ctypes.util.find_library("m") or "libSystem.B.dylib")
_libm.fmaf.restype = ctypes.c_float
_libm.fmaf.argtypes = [ctypes.c_float, ctypes.c_float, ctypes.c_float]


def fmaf(a: np.float32, b: np.float32, c: np.float32) -> np.float32:
    return np.float32(_libm.fmaf(float(a), float(b), float(c)))


def _plain(a: np.ndarray, b: np.ndarray) -> np.float32:
    return np.float32(np.float32(a[0] * b[0]) + np.float32(a[1] * b[1])) + np.float32(a[2] * b[2])


def _clang(a: np.ndarray, b: np.ndarray) -> np.float32:
    # clang -ffp-contract=on: a0*b0 + a1*b1 fuses the left product, then the outer + fuses a2*b2
    return fmaf(a[2], b[2], fmaf(a[0], b[0], np.float32(a[1] * b[1])))


def _chain(a: np.ndarray, b: np.ndarray) -> np.float32:
    return fmaf(a[2], b[2], fmaf(a[1], b[1], np.float32(a[0] * b[0])))


def _rev(a: np.ndarray, b: np.ndarray) -> np.float32:
    return fmaf(a[0], b[0], fmaf(a[1], b[1], np.float32(a[2] * b[2])))


ORDERS: Dict[str, Callable[[np.ndarray, np.ndarray], np.float32]] = {
    "plain": _plain,
    "clang_contract": _clang,
    "fma_chain": _chain,
    "fma_rev": _rev,
}
PRODUCT_ORDER = "clang_contract"  # the order the ``inputs`` check selects (see inputs.json)


def product(local: np.ndarray, centre: np.ndarray, order: str = PRODUCT_ORDER) -> np.ndarray:
    """float32 ``local @ centre`` in the given operation order (C++ ``oLocal * oCenter``)."""
    L = np.asarray(local, np.float32)
    C = np.asarray(centre, np.float32)
    f = ORDERS[order]
    out = np.zeros((3, 3), np.float32)
    for i in range(3):
        for j in range(3):
            out[i, j] = f(L[i, :], C[:, j])
    return out


def orientation_of(m: Dict[str, Any]) -> np.ndarray:
    """The orientation the C++ driver evaluated (its float32 product), as float64."""
    return product(np.asarray(m["local"]), np.asarray(m["centre"])).astype(np.float64)


def cmd_inputs(a: argparse.Namespace) -> int:
    d = Path(a.dir).resolve()
    stored = load_stored()
    res: Dict[str, Any] = dict(meta=git_meta(), voxels={}, orders={})
    meta_bad: List[Any] = []
    order_ok = {o: 0 for o in ORDERS}
    order_bitflag_ok = {o: 0 for o in ORDERS}
    order_maxdiff_ok = {o: 0 for o in ORDERS}
    # the limit of the evidence: how often another order gives a different product yet matches both
    # statistics, and how often the statistics leave more than the largest element unpinned
    order_differs = {o: 0 for o in ORDERS}
    order_differs_both_match = {o: 0 for o in ORDERS}
    pp_elements_differ: Counter = Counter()  # float64-cast vs selected order: elements differing
    kinds: Counter = Counter()
    single_identity_local = 0
    n_pairs = 0
    n_seen = 0
    for v in PC.PARITY_VOXELS:
        recs = load_inputs(d, v)
        res["voxels"][str(v)] = len(recs)
        for m in recs:
            key = (v, int(m["id"]))
            if key not in stored:
                meta_bad.append(dict(voxel=v, id=m["id"], why="not stored"))
                continue
            n_seen += 1
            s = stored[key]
            same = s["label"] == m["label"] and int(s["kind"]) == int(m["kind"])
            for k in META_KEYS:
                same = same and s.get("meta_" + k) == m.get(k)
            if not same:
                meta_bad.append(dict(voxel=v, id=m["id"], label=m["label"]))
            kinds[int(m["kind"])] += 1
            if int(m["kind"]) != 1:
                single_identity_local += int(
                    np.array_equal(np.asarray(m["local"], np.float64), np.eye(3))
                )
                continue
            n_pairs += 1
            L = np.asarray(m["local"], np.float64)
            C = np.asarray(m["centre"], np.float64)
            pp = (L @ C).astype(np.float32)
            P_sel = product(L, C, PRODUCT_ORDER)
            pp_elements_differ[int((pp != P_sel).sum())] += 1
            for o in ORDERS:
                P = product(L, C, o)
                be = bool(np.array_equal(pp, P))
                md = float(np.abs(pp.astype(np.float64) - P.astype(np.float64)).max())
                ok_b = be == bool(s["product_bit_equal"])
                ok_m = md == float(s["product_max_abs_diff"])
                order_bitflag_ok[o] += int(ok_b)
                order_maxdiff_ok[o] += int(ok_m)
                order_ok[o] += int(ok_b and ok_m)
                if not np.array_equal(P, P_sel):
                    order_differs[o] += 1
                    order_differs_both_match[o] += int(ok_b and ok_m)
    for o in ORDERS:
        res["orders"][o] = dict(
            n_pairs=n_pairs,
            bit_equal_flag_matches=order_bitflag_ok[o],
            max_abs_diff_matches=order_maxdiff_ok[o],
            both_match=order_ok[o],
            product_differs_from_selected=order_differs[o],
            product_differs_from_selected_yet_both_match=order_differs_both_match[o],
        )
    res["n_records"] = n_seen
    res["records_by_kind"] = {str(k): kinds[k] for k in sorted(kinds)}
    res["single_records_identity_local"] = single_identity_local
    # per pair record, the elements (of 9) where the float64 product cast to float32 differs from
    # the selected order's product; the stored largest difference pins only one of them
    res["pair_records_by_elements_differing"] = {
        str(k): pp_elements_differ[k] for k in sorted(pp_elements_differ)
    }
    res["pair_records_two_or_more_elements_differing"] = sum(
        c for k, c in pp_elements_differ.items() if k >= 2
    )
    res["n_stored"] = len(stored)
    res["labels_and_indices_differ"] = meta_bad
    res["selected_order"] = PRODUCT_ORDER
    res["selected_order_matches_all"] = order_ok[PRODUCT_ORDER] == n_pairs
    Path(a.out).resolve().write_text(json.dumps(res, indent=1) + "\n")
    print(json.dumps({k: res[k] for k in ("voxels", "orders", "n_records", "n_stored")}, indent=1))
    print(f"labels/indices differing: {len(meta_bad)}")
    return 0 if (not meta_bad and res["selected_order_matches_all"]) else 1


# ----------------------------------------------------------------------------------------------
# run: torch and the fused prototype on every record
# ----------------------------------------------------------------------------------------------
def agg(oi: Any) -> Dict[str, Any]:
    return dict(
        pixov=int(oi.pixel_overlap),
        pixon=int(oi.pixel_on_detector),
        pkov=int(oi.peak_overlap),
        pkon=int(oi.peak_on_detector),
        quality=float(oi.quality),
        conf=float(oi.confidence),
    )


def agg_of_tuple(t: Tuple[Any, ...]) -> Dict[str, Any]:
    from icenine.cost_functions import OverlapInfo

    oi = OverlapInfo()
    (
        oi.pixel_overlap,
        oi.pixel_on_detector,
        oi.peak_overlap,
        oi.peak_on_detector,
        oi.quality,
        oi.n_quality_points,
        oi.detectors_overlap,
    ) = t
    return agg(oi)


def classify(rows: List[Dict[str, Any]], pixel_width: float) -> Dict[str, Any]:
    """Whether every traced row is a pixel-edge flip within the parity rule's 0.01 px."""
    causes = sorted({r["cause"] for r in rows})
    worst_px = 0.0
    ok = bool(rows)
    for r in rows:
        if r["cause"] == "pixel_integer":
            near = min(r["torch_px_from_edge"], r["fused_px_from_edge"])
        elif r["cause"] == "pixel_negative":
            near = abs(r["torch"]) / pixel_width  # detector coordinate 0 is the edge of pixel -1
        else:
            ok = False
            continue
        worst_px = max(worst_px, near)
        ok = ok and near <= PIXEL_EDGE_PX
    return dict(causes=causes, worst_near_px=worst_px, within_rule=ok)


class LibraryFused:
    """Variant ``library`` (S2 Step 4): the library kernel (``VoxelCostFunction`` with the key
    FusedCost on, the shipped g_lab order) in the prototype's ``debug`` interface. The
    intermediates come from ``fused_debug``, the instrumented fused path (``parity_common``'s
    probes patch torch functions the fused path never calls); ``debug`` also checks that the
    plain ``evaluate`` gives the same result and counts every call."""

    def __init__(self, cf: Any) -> None:
        cf.fused_cost = True
        self.cf = cf

    def debug(self, R: np.ndarray, verts: Any, phase: int) -> Dict[str, Any]:
        n0 = self.cf.eval_count
        oi = self.cf.evaluate(R, verts, phase)
        if self.cf.eval_count != n0 + 1:
            raise SystemExit("library fused evaluate did not count one evaluation")
        d = dict(self.cf.fused_debug(R, verts, phase))
        r = SF.result_tuple(oi)
        if tuple(d["result"]) != r:
            raise SystemExit(f"fused_debug {tuple(d['result'])} differs from evaluate {r}")
        d["K"] = int(d["g_lab"].shape[0])
        return d


def make_fused(cf: Any, variant: str) -> Any:
    return LibraryFused(cf) if variant == "library" else SF.FusedEval(cf, variant)


def cmd_run(a: argparse.Namespace) -> int:
    torch.set_num_threads(1)
    d = Path(a.dir).resolve()
    out = Path(a.out).resolve()
    with contextlib.redirect_stdout(io.StringIO()):
        ctx = PC.PyContext()
    cfs = {arm: ctx.cost_fn(*ARMS[arm]) for arm in ARMS}
    mirror_cf = {arm: ctx.cost_fn(*ARMS[arm]) for arm in FUSED_ARMS}
    fused = {
        (arm, v): make_fused(ctx.cost_fn(*ARMS[arm]), v) for arm in FUSED_ARMS for v in a.variants
    }
    n_det = len(ctx.setup.detector_list)
    pw = float(ctx.setup.detector_list[0].pixel_width)
    rows: List[Dict[str, Any]] = []
    t0 = time.perf_counter()
    for v in a.voxels:
        recs = load_inputs(d, v)
        if a.limit:
            recs = recs[: a.limit]
        verts = ctx.verts(v)
        phase = ctx.mic.voxels[v].phase
        for m in recs:
            R = orientation_of(m)
            row: Dict[str, Any] = dict(voxel=v, id=int(m["id"]), label=m["label"], kind=m["kind"])
            with contextlib.redirect_stdout(io.StringIO()):
                for arm, cf in cfs.items():
                    row[arm] = dict(torch=agg(cf.evaluate(R, verts, phase)))
                for arm in FUSED_ARMS:
                    mcf = mirror_cf[arm]
                    mm = TM.mirror(mcf, R, verts, phase)
                    lib_t = (
                        row[arm]["torch"]["pixov"],
                        row[arm]["torch"]["pixon"],
                        row[arm]["torch"]["pkov"],
                        row[arm]["torch"]["pkon"],
                        row[arm]["torch"]["quality"],
                    )
                    mt = tuple(mm["result"])
                    row[arm]["mirror_equal"] = bool(
                        tuple(mt[:4]) == lib_t[:4] and float(mt[4]) == lib_t[4]
                    )
                    ref = TM.reference64(mcf, R, verts, phase)
                    for var in a.variants:
                        dd = fused[(arm, var)].debug(R, verts, phase)
                        fa = agg_of_tuple(tuple(dd["result"]))
                        tr = S2P.trace(mm, dd, ref, mcf, n_det, every=a.every_coord)
                        ent: Dict[str, Any] = dict(agg=fa)
                        differs = tuple(dd["result"]) != tuple(mt)
                        ent["differs_from_torch"] = differs
                        if differs:
                            ent["trace"] = tr
                            ent["verdict"] = classify(tr, pw)
                        elif tr:  # a stage-D input differs, the OverlapInfo does not
                            ent["silent_rows"] = len(tr)
                            ent["silent_trace"] = tr
                            ent["silent_verdict"] = classify(tr, pw)
                        row[arm][var] = ent
            rows.append(row)
            if len(rows) % 200 == 0:
                print(f"  {len(rows)} records, {time.perf_counter() - t0:.0f} s", flush=True)
    res = dict(
        meta=dict(
            git_meta(),
            voxels=a.voxels,
            variants=a.variants,
            every_coord=a.every_coord,
            product_order=PRODUCT_ORDER,
            threads=torch.get_num_threads(),
            seconds=time.perf_counter() - t0,
        ),
        rows=rows,
    )
    out.parent.mkdir(parents=True, exist_ok=True)
    with gzip.GzipFile(filename="", mode="wb", fileobj=open(out, "wb"), mtime=0) as gz:
        gz.write(json.dumps(res, default=float).encode())
    print(f"wrote {len(rows)} rows to {out.name} in {time.perf_counter() - t0:.0f} s")
    return 0


# ----------------------------------------------------------------------------------------------
# outlier: the 1.50 px float64-reference outlier of Step 0
# ----------------------------------------------------------------------------------------------
def cmd_outlier(a: argparse.Namespace) -> int:
    torch.set_num_threads(1)
    out = Path(a.out).resolve()  # before loading: the dataset may change the working directory
    with contextlib.redirect_stdout(io.StringIO()):
        ds = S1.load_dataset(a.dataset)
    R = S1.orientation_set(ds, a.n_near, a.n_fz)
    res: Dict[str, Any] = dict(
        meta=dict(git_meta(), dataset=a.dataset, n=len(R), min_px=a.min_px), modes={}
    )
    for mode in a.modes:
        kw = S2P.MODES[mode]
        cf = ds.cost_fn(kw)
        fe = SF.FusedEval(ds.cost_fn(kw), "portable")
        cf_raw = ds.cost_fn(kw)
        found: List[Dict[str, Any]] = []
        n_peaks = 0
        for i in range(len(R)):
            with contextlib.redirect_stdout(io.StringIO()):
                m = TM.mirror(cf, R[i], ds.verts, ds.phase)
                d = fe.debug(R[i], ds.verts, ds.phase)
                ref = TM.reference64(cf, R[i], ds.verts, ds.phase)
                cf_raw.cpp_edge_rules = False  # the float64 coordinate without the -1 rule
                raw = TM.reference64(cf_raw, R[i], ds.verts, ds.phase)
            ft = m["frame"] >= 0
            n_peaks += int(ft.sum())
            et = np.abs(m["pix"].astype(np.float64) - ref["pix"])
            ef = np.abs(d["pix"].astype(np.float64) - ref["pix"])
            bad = np.argwhere(((et > a.min_px) | (ef > a.min_px)) & ft[:, None, None, None])
            for r, det, vx, ax in bad:
                dev = ds.setup.detector_list[int(det)]
                half = float(np.float32(dev.pixel_half_width if ax == 0 else dev.pixel_half_height))
                width = float(np.float32(dev.pixel_width if ax == 0 else dev.pixel_height))
                j32 = float(m["jk"][r, det, vx, ax])
                j64 = float(raw["pix"][r, det, vx, ax]) * width - half
                found.append(
                    dict(
                        orientation=int(i),
                        row=int(r),
                        det=int(det),
                        vertex=int(vx),
                        axis="col" if ax == 0 else "row",
                        hit_torch=int(m["hit"][r, det]),
                        hit_fused=int(d["hit"][r, det]),
                        torch_pix=float(m["pix"][r, det, vx, ax]),
                        fused_pix=float(d["pix"][r, det, vx, ax]),
                        ref_pix=float(ref["pix"][r, det, vx, ax]),
                        torch_err=float(et[r, det, vx, ax]),
                        fused_err=float(ef[r, det, vx, ax]),
                        coord_torch_f32=j32,
                        coord_f64=j64,
                        coord_torch_ulp_from_0=int(
                            S2P.ulp_diff(np.array([j32], np.float32), np.zeros(1, np.float32))[0]
                        ),
                        pixel_size=width,
                        coord_torch_over_pixel=j32 / width,
                        coord_f64_over_pixel=j64 / width,
                        torch_stage_d_int=int(
                            -1
                            if m["pix"][r, det, vx, ax] < 0
                            else np.floor(m["pix"][r, det, vx, ax])
                        ),
                        fused_stage_d_int=int(
                            -1
                            if d["pix"][r, det, vx, ax] < 0
                            else np.floor(d["pix"][r, det, vx, ax])
                        ),
                        ref_stage_d_int=int(
                            -1
                            if ref["pix"][r, det, vx, ax] < 0
                            else np.floor(ref["pix"][r, det, vx, ax])
                        ),
                        omega_torch=float(m["omega"][r]),
                        omega_fused=float(d["omega"][r]),
                        omega_ref=float(ref["omega"][r]),
                        asin_arg_torch=float(m["bragg_num"][r % m["K"]])
                        / float(m["bragg_den"][r % m["K"]]),
                        torch_result=list(m["result"]),
                        fused_result=list(d["result"]),
                    )
                )
        res["modes"][mode] = dict(n_stage_d_rows=n_peaks, entries=found)
        print(f"{mode}: {len(found)} coordinates beyond {a.min_px} px", flush=True)
        for e in found:
            print(
                "  ",
                {
                    k: e[k]
                    for k in (
                        "orientation",
                        "row",
                        "det",
                        "vertex",
                        "axis",
                        "torch_pix",
                        "fused_pix",
                        "ref_pix",
                        "coord_torch_f32",
                        "coord_f64",
                    )
                },
            )
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(res, indent=1, default=float) + "\n")
    return 0


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("inputs")
    p.add_argument("--dir", required=True)
    p.add_argument("--out", required=True)
    p = sub.add_parser("run")
    p.add_argument("--dir", required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--voxels", type=int, nargs="+", default=PC.PARITY_VOXELS)
    p.add_argument(
        "--variants",
        nargs="+",
        default=["portable", "accelerate"],
        choices=["portable", "accelerate", "library"],
        help="library: the S2 library kernel (key FusedCost; Step 4)",
    )
    p.add_argument("--every-coord", action="store_true", help="trace every flipped coordinate")
    p.add_argument("--limit", type=int, default=0, help="first N records per voxel (pilot)")
    p = sub.add_parser("outlier")
    p.add_argument("--out", required=True)
    p.add_argument("--dataset", default="full18558")
    p.add_argument("--modes", nargs="+", default=["edge", "screen_edge"])
    p.add_argument("--n-near", type=int, default=500)
    p.add_argument("--n-fz", type=int, default=500)
    p.add_argument("--min-px", type=float, default=0.1)
    a = ap.parse_args()
    return {"inputs": cmd_inputs, "run": cmd_run, "outlier": cmd_outlier}[a.cmd](a)


if __name__ == "__main__":
    sys.exit(main())
