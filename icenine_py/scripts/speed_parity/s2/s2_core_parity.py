"""S2 Steps 1-2: the library fused evaluate (key FusedCost) against the torch ``evaluate``, with
every differing orientation traced to its peaks under the parity rule (D-20261027-I2), for each
g_lab operation order (D-20261027-I5), and against the Step 0 prototype bit for bit.

For each dataset and covered mode (``s2_parity.MODES``) on S1's orientation set:

* the torch ``evaluate`` (the reference result) and ``torch_mirror.mirror`` (its intermediates;
  must equal the result, else the run stops) and ``torch_mirror.reference64`` (float64);
* per g_lab order (``fma``: the fmaf chain of the prototype; ``plain``: ``(a*b + c*d) + e*f``,
  the library default since Step 1) a ``VoxelCostFunction(fused_cost=True)``: its ``evaluate``
  result and ``fused_debug`` intermediates (the debug result must equal the evaluate result);
* order ``fma`` only: the Step 0 prototype (``s2fused.FusedEval``, portable build), every
  intermediate and the result compared bit for bit with the library's (the port must be exact).

A differing orientation is traced with ``s2_parity.trace`` (the rows whose stage-D input differs
and their cause). The parity rule: every traced row is a pixel-edge flip, cause ``pixel_integer``
with the torch or the fused pixel coordinate within ``--bound-px`` (0.01) of the integer; a frame,
eta, observability, ``min_sin_eta``, t > 0 or ``pixel_negative`` row, or a differing orientation
without a traced row, fails.

  uv run python scripts/speed_parity/s2/s2_core_parity.py --out benchmarks/s2_core/glab_order.json
  uv run python scripts/speed_parity/s2/s2_core_parity.py --out /tmp/p.json --datasets three0 \
      --n-near 20 --n-fz 20                                   # pilot
"""

import argparse
import contextlib
import io
import json
import sys
import time
from pathlib import Path
from typing import Any, Dict, List

import numpy as np
import torch

HERE = Path(__file__).resolve().parent
ICE = HERE.parents[2]
sys.path.insert(0, str(ICE))
sys.path.insert(0, str(HERE.parent))
sys.path.insert(0, str(HERE))

import s1_parity as S1  # noqa: E402
import s2_parity as SP  # noqa: E402
import s2fused as SF  # noqa: E402
import torch_mirror as TM  # noqa: E402
from icenine import cost_functions as CF  # noqa: E402

ORDERS = {"fma": 0, "plain": 1}
DEBUG_KEYS = (
    "g_lab",
    "obs",
    "omega1",
    "omega2",
    "omega",
    "cos",
    "sin",
    "eta",
    "valid",
    "frame",
    "lab_normal",
    "reflected",
    "pix",
    "hit",
    "t",
)


def fused_cost_fn(ds: Any, kw: Dict[str, Any], order: str) -> Any:
    """A FusedCost cost function whose C context is built with ``order`` (the module switch is read
    when the context is built, at the first evaluation)."""
    cf = ds.cost_fn(dict(kw, fused_cost=True))
    saved = CF._FUSED_GLAB_ORDER
    CF._FUSED_GLAB_ORDER = ORDERS[order]
    try:
        cf.fused_debug(ds.truth, ds.verts, ds.phase)  # builds the context
    finally:
        CF._FUSED_GLAB_ORDER = saved
    return cf


def bit_equal(a: np.ndarray, b: np.ndarray) -> bool:
    a, b = np.asarray(a), np.asarray(b)
    return a.shape == b.shape and a.dtype == b.dtype and a.tobytes() == b.tobytes()


def classify(rows: List[Dict[str, Any]], bound_px: float) -> List[str]:
    """Why each traced row fails the parity rule ([] = every row is an accepted pixel-edge flip)."""
    bad = []
    for r in rows:
        if r["cause"] != "pixel_integer":
            bad.append(r["cause"])
        elif min(r["torch_px_from_edge"], r["fused_px_from_edge"]) > bound_px:
            bad.append("pixel_integer beyond the bound")
    return bad


def run_dataset(
    ds: Any, R: np.ndarray, modes: List[str], orders: List[str], bound_px: float
) -> Dict[str, Any]:
    res: Dict[str, Any] = {}
    n_det = len(ds.setup.detector_list)
    for mode in modes:
        kw = SP.MODES[mode]
        cf_lib, cf_mir = ds.cost_fn(kw), ds.cost_fn(kw)
        fused = {o: fused_cost_fn(ds, kw, o) for o in orders}
        proto = SF.FusedEval(ds.cost_fn(kw), "portable") if "fma" in orders else None
        st: Dict[str, Dict[str, Any]] = {
            o: dict(
                differing=[],
                field=dict.fromkeys(SF.FIELDS, 0),
                rows=0,
                silent=0,
                ref_side=dict(torch=0, fused=0),
                max_px_from_edge=0.0,
                max_ulp=0,
                g_lab=SP.Stat(),
                pixel=SP.Stat(),
            )
            for o in orders
        }
        proto_diff = 0
        t0 = time.perf_counter()
        for i in range(len(R)):
            with contextlib.redirect_stdout(io.StringIO()):
                lib = SF.result_tuple(cf_lib.evaluate(R[i], ds.verts, ds.phase))
                m = TM.mirror(cf_mir, R[i], ds.verts, ds.phase)
            if tuple(m["result"]) != lib:
                raise SystemExit(f"mirror differs from the library: {mode} {i}")
            ref = TM.reference64(cf_mir, R[i], ds.verts, ds.phase)
            for o in orders:
                cf = fused[o]
                res_f = SF.result_tuple(cf.evaluate(R[i], ds.verts, ds.phase))
                d = cf.fused_debug(R[i], ds.verts, ds.phase)
                if tuple(d["result"]) != res_f:
                    raise SystemExit(f"fused_debug result differs from evaluate: {mode} {i} {o}")
                if o == "fma" and proto is not None:
                    p = proto.debug(R[i], ds.verts, ds.phase)
                    same = tuple(p["result"]) == res_f and all(
                        bit_equal(p[k], d[k]) for k in DEBUG_KEYS
                    )
                    proto_diff += not same
                s = st[o]
                s["g_lab"].add(m["g_lab"], d["g_lab"])
                fb = (m["frame"] >= 0) & (d["frame"] >= 0)
                s["pixel"].add(m["pix"][fb], d["pix"][fb])
                rows = SP.trace(m, d, ref, cf_mir, n_det)
                if res_f != lib:
                    for j, f in enumerate(SF.FIELDS):
                        s["field"][f] += int(res_f[j] != lib[j])
                    bad = classify(rows, bound_px) if rows else ["no traced row"]
                    s["differing"].append(dict(index=i, rows=rows, failures=bad))
                    s["rows"] += len(rows)
                    for r in rows:
                        side = r.get("ref_side")
                        if side in s["ref_side"]:
                            s["ref_side"][side] += 1
                        if r["cause"] == "pixel_integer":
                            s["max_px_from_edge"] = max(
                                s["max_px_from_edge"],
                                min(r["torch_px_from_edge"], r["fused_px_from_edge"]),
                            )
                        s["max_ulp"] = max(s["max_ulp"], int(r["margin_ulp"]))
                elif rows:
                    s["silent"] += 1
        dt = time.perf_counter() - t0
        res[mode] = dict(
            n=len(R),
            prototype_differs=proto_diff if proto is not None else None,
            orders={
                o: dict(
                    n_orientations_differing=len(s["differing"]),
                    n_failing=sum(bool(x["failures"]) for x in s["differing"]),
                    failures=[x for x in s["differing"] if x["failures"]],
                    differing_index=[x["index"] for x in s["differing"]],
                    flipped_rows=s["rows"],
                    field_orientations_differing=s["field"],
                    silent_input_differences=s["silent"],
                    ref_side=s["ref_side"],
                    max_px_from_edge=s["max_px_from_edge"],
                    max_margin_ulp=s["max_ulp"],
                    g_lab=s["g_lab"].as_dict(),
                    pixel=s["pixel"].as_dict(),
                )
                for o, s in st.items()
            },
        )
        line = "  ".join(
            f"{o}: {len(s['differing'])} differ "
            f"({sum(bool(x['failures']) for x in s['differing'])} fail, {s['rows']} rows)"
            for o, s in st.items()
        )
        pd = f"  prototype differs {proto_diff}" if proto is not None else ""
        print(f"  {mode:15s} n={len(R)} {dt:5.1f}s  {line}{pd}", flush=True)
    return res


def main() -> int:
    torch.set_num_threads(1)
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--out", required=True)
    ap.add_argument(
        "--datasets", nargs="+", default=["three0"] + [f"full{v}" for v in S1.FULL_VOXELS]
    )
    ap.add_argument("--modes", nargs="+", default=list(SP.MODES), choices=list(SP.MODES))
    ap.add_argument("--orders", nargs="+", default=list(ORDERS), choices=list(ORDERS))
    ap.add_argument("--n-near", type=int, default=500)
    ap.add_argument("--n-fz", type=int, default=500)
    ap.add_argument("--bound-px", type=float, default=0.01)
    a = ap.parse_args()
    out_path = Path(a.out).resolve()
    result: Dict[str, Any] = dict(
        meta=dict(
            head=S1.git_head(),
            dirty=S1.git_dirty(),
            dirty_files=SP.dirty_files(),
            torch=torch.__version__,
            numpy=np.__version__,
            threads=torch.get_num_threads(),
            seed=S1.SEED,
            n_near=a.n_near,
            n_fz=a.n_fz,
            bound_px=a.bound_px,
            orders=a.orders,
            modes={k: SP.MODES[k] for k in a.modes},
        ),
        datasets={},
    )
    for name in a.datasets:
        t0 = time.perf_counter()
        with contextlib.redirect_stdout(io.StringIO()):
            ds = S1.load_dataset(name)
        R = S1.orientation_set(ds, a.n_near, a.n_fz)
        print(f"[{name}] loaded in {time.perf_counter() - t0:.1f} s, {len(R)} orientations")
        result["datasets"][name] = run_dataset(ds, R, a.modes, a.orders, a.bound_px)
    totals = {
        o: dict(
            differing=sum(
                r["orders"][o]["n_orientations_differing"]
                for d in result["datasets"].values()
                for r in d.values()
            ),
            failing=sum(
                r["orders"][o]["n_failing"] for d in result["datasets"].values() for r in d.values()
            ),
            flipped_rows=sum(
                r["orders"][o]["flipped_rows"]
                for d in result["datasets"].values()
                for r in d.values()
            ),
        )
        for o in a.orders
    }
    proto = [
        r["prototype_differs"]
        for d in result["datasets"].values()
        for r in d.values()
        if r["prototype_differs"] is not None
    ]
    result["totals"] = dict(orders=totals, prototype_differs=sum(proto) if proto else None)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(result, indent=1) + "\n")
    print(f"wrote {out_path.name}; totals {json.dumps(result['totals'])}")
    bad = sum(t["failing"] for t in totals.values()) + (sum(proto) if proto else 0)
    return 0 if bad == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
