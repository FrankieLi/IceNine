"""S2 Step 4: the library fused evaluate (key FusedCost, the shipped g_lab order) against the torch
``evaluate``, every flipped pixel coordinate traced under the parity rule (D-20261027-I2).

Plan ``docs/todo_s2_fused_c_eval.md``, "Steps 4-6", item 4 and Parity item 1. For each dataset
and covered mode (``s2_parity.MODES``) on S1's orientation set (500 near the truth + 500 FZ points):

* torch: the library ``evaluate`` with the key off, and ``torch_mirror.mirror`` (its intermediates;
  must equal the result, else the run stops) and ``torch_mirror.reference64`` (float64);
* fused: a ``VoxelCostFunction(fused_cost=True)``: its ``evaluate`` and its ``fused_debug``
  intermediates (the instrumented fused path; the debug result must equal ``evaluate``). The
  module switch ``cost_functions._FUSED_GLAB_ORDER`` is not touched (the shipped order).

Every orientation with a differing stage-D input (``s2_parity.trace(..., every=True)``: one entry
per flipped pixel coordinate, not only the worst per row and detector) is traced, whether its
``OverlapInfo`` differs ("differing") or not ("silent"). The rule: every entry is a pixel-edge
flip (cause ``pixel_integer``, or ``pixel_negative``: the ``CppEdgeRules`` -1 rule at detector
coordinate 0, which is the edge of pixel 0) with the torch or the fused coordinate within
``--bound-px`` (0.01) of the edge; a frame, eta, observability, ``min_sin_eta`` or t > 0 entry,
or a differing orientation without an entry, fails. Per flipped coordinate: which value is nearer
the float64 reference ("nearer": torch, fused or tie) and on whose side of the edge the float64
value lies ("side"). Per intermediate the ulp statistics as in Step 0d.

  uv run python scripts/speed_parity/s2/s2_lib_parity.py --out benchmarks/s2_results/parity.json
  uv run python scripts/speed_parity/s2/s2_lib_parity.py --out /tmp/p.json --datasets three0 \
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

ACCEPTED = ("pixel_integer", "pixel_negative")


def near_px(e: Dict[str, Any], ds: Any) -> float:
    """Distance (px) of the nearer of the torch and fused coordinates from the pixel edge."""
    if e["cause"] == "pixel_integer":
        return float(min(e["torch_px_from_edge"], e["fused_px_from_edge"]))
    # pixel_negative: e["torch"] is torch's detector coordinate (mm); the edge is coordinate 0
    dev = ds.setup.detector_list[int(e["det"])]
    width = float(np.float32(dev.pixel_width if e["coord"] % 2 == 0 else dev.pixel_height))
    return abs(float(e["torch"])) / width


def check(entries: List[Dict[str, Any]], ds: Any, bound_px: float) -> List[str]:
    """Why the traced entries fail the parity rule ([] = all accepted pixel-edge flips)."""
    bad = []
    for e in entries:
        if e["cause"] not in ACCEPTED:
            bad.append(e["cause"])
        elif near_px(e, ds) > bound_px:
            bad.append(f"{e['cause']} beyond the bound")
    return bad


class Acc:
    """Counts of one (dataset, mode)."""

    def __init__(self) -> None:
        self.differing: List[int] = []
        self.silent: List[int] = []
        self.failures: List[Dict[str, Any]] = []
        self.field = dict.fromkeys(SF.FIELDS, 0)
        self.coords = dict(differing=0, silent=0)
        self.rows_dets = dict(differing=0, silent=0)  # distinct (row, detector) pairs
        self.multi_coord_rows_dets = (
            0  # (row, detector) pairs with more than one flipped coordinate
        )
        self.causes: Dict[str, int] = {}
        self.nearer = dict(torch=0, fused=0, tie=0)
        self.side = dict(torch=0, fused=0)
        self.max_near_px = 0.0
        self.max_ulp = 0
        self.max_fused_ulp = 0
        self.stats = {k: SP.Stat() for k in SP.INTERMEDIATES}
        self.frame_rows = 0
        self.frame_differ = 0


def run_dataset(ds: Any, R: np.ndarray, modes: List[str], bound_px: float) -> Dict[str, Any]:
    out: Dict[str, Any] = {}
    n_det = len(ds.setup.detector_list)
    for mode in modes:
        kw = SP.MODES[mode]
        cf_lib, cf_mir = ds.cost_fn(kw), ds.cost_fn(kw)
        cf_f = ds.cost_fn(dict(kw, fused_cost=True))
        a = Acc()
        t0 = time.perf_counter()
        for i in range(len(R)):
            with contextlib.redirect_stdout(io.StringIO()):
                lib = SF.result_tuple(cf_lib.evaluate(R[i], ds.verts, ds.phase))
                m = TM.mirror(cf_mir, R[i], ds.verts, ds.phase)
            if tuple(m["result"]) != lib:
                raise SystemExit(f"mirror differs from the library: {ds.name} {mode} {i}")
            res_f = SF.result_tuple(cf_f.evaluate(R[i], ds.verts, ds.phase))
            d = cf_f.fused_debug(R[i], ds.verts, ds.phase)
            if tuple(d["result"]) != res_f:
                raise SystemExit(f"fused_debug differs from evaluate: {ds.name} {mode} {i}")
            # intermediates, with Step 0d's masks (s2_parity.run_dataset): rows both sides computed
            a.stats["g_lab"].add(m["g_lab"], d["g_lab"])
            ob = (m["obs"] == 1) & (d["obs"] == 1)
            a.stats["omega_bragg"].add(m["omega1"][ob], d["omega1"][ob])
            a.stats["omega_bragg"].add(m["omega2"][ob], d["omega2"][ob])
            vb = (m["valid"] == 1) & (d["valid"] == 1)
            a.stats["omega"].add(m["omega"][vb], d["omega"][vb])
            a.stats["eta"].add(m["eta"][vb], d["eta"][vb])
            fb = (m["frame"] >= 0) & (d["frame"] >= 0)
            a.frame_rows += int(fb.sum())
            a.frame_differ += int((m["frame"][fb] != d["frame"][fb]).sum())
            a.stats["pixel"].add(m["pix"][fb], d["pix"][fb])
            a.stats["t"].add(m["t"][fb], d["t"][fb])
            ref = TM.reference64(cf_mir, R[i], ds.verts, ds.phase)
            entries = SP.trace(m, d, ref, cf_mir, n_det, every=True)
            differs = res_f != lib
            if differs:
                for j, f in enumerate(SF.FIELDS):
                    a.field[f] += int(res_f[j] != lib[j])
            if not differs and not entries:
                continue
            kind = "differing" if differs else "silent"
            (a.differing if differs else a.silent).append(i)
            bad = check(entries, ds, bound_px) if entries else ["no traced entry"]
            if bad:
                a.failures.append(dict(index=i, kind=kind, failures=bad, entries=entries))
            a.coords[kind] += len(entries)
            rd: Dict[Any, int] = {}
            for e in entries:
                a.causes[e["cause"]] = a.causes.get(e["cause"], 0) + 1
                if "det" in e:
                    rd[(e["row"], e["det"])] = rd.get((e["row"], e["det"]), 0) + 1
                if e["cause"] in ACCEPTED:
                    a.max_near_px = max(a.max_near_px, near_px(e, ds))
                a.max_ulp = max(a.max_ulp, int(e.get("margin_ulp", 0)))
                if e["cause"] == "pixel_integer":
                    a.max_fused_ulp = max(a.max_fused_ulp, int(e["fused_margin_ulp"]))
                    a.side[e["ref_side"]] += 1
                    et = abs(e["torch"] - e["ref"])
                    ef = abs(e["fused"] - e["ref"])
                    a.nearer["torch" if et < ef else "fused" if ef < et else "tie"] += 1
            a.rows_dets[kind] += len(rd)
            a.multi_coord_rows_dets += sum(c > 1 for c in rd.values())
        dt = time.perf_counter() - t0
        out[mode] = dict(
            n=len(R),
            seconds=dt,
            n_differing=len(a.differing),
            n_silent=len(a.silent),
            n_failing=len(a.failures),
            failures=a.failures,
            differing_index=a.differing,
            silent_index=a.silent,
            field_orientations_differing=a.field,
            flipped_coords=a.coords,
            flipped_rows_dets=a.rows_dets,
            rows_dets_with_several_coords=a.multi_coord_rows_dets,
            causes=a.causes,
            nearer_f64=a.nearer,
            f64_side=a.side,
            max_near_px=a.max_near_px,
            max_margin_ulp=a.max_ulp,
            max_fused_margin_ulp=a.max_fused_ulp,
            frame_rows_both=a.frame_rows,
            frame_differ=a.frame_differ,
            intermediates={k: s.as_dict() for k, s in a.stats.items()},
        )
        print(
            f"  {mode:15s} n={len(R)} {dt:5.1f}s  differ {len(a.differing)}  silent "
            f"{len(a.silent)}  fail {len(a.failures)}  coords {a.coords}  causes {a.causes}"
            f"  max {a.max_near_px:.4f} px",
            flush=True,
        )
    return out


def main() -> int:
    torch.set_num_threads(1)
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--out", required=True)
    ap.add_argument(
        "--datasets", nargs="+", default=["three0"] + [f"full{v}" for v in S1.FULL_VOXELS]
    )
    ap.add_argument("--modes", nargs="+", default=list(SP.MODES), choices=list(SP.MODES))
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
            glab_order=CF._FUSED_GLAB_ORDER,
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
        result["datasets"][name] = run_dataset(ds, R, a.modes, a.bound_px)
    recs = [r for d in result["datasets"].values() for r in d.values()]
    result["totals"] = dict(
        differing=sum(r["n_differing"] for r in recs),
        silent=sum(r["n_silent"] for r in recs),
        failing=sum(r["n_failing"] for r in recs),
        flipped_coords=sum(sum(r["flipped_coords"].values()) for r in recs),
        frame_differ=sum(r["frame_differ"] for r in recs),
    )
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(result, indent=1, default=float) + "\n")
    print(f"wrote {out_path.name}; totals {json.dumps(result['totals'])}")
    return 0 if result["totals"]["failing"] == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
