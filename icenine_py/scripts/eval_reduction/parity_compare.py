"""Task A step 0: compare the C++ driver output with the Python cost, orientation by orientation.

  uv run python scripts/eval_reduction/parity_compare.py --dir DIR \
      --out benchmarks/eval_reduction/A/step0

Reads DIR/in_<v>.json (labels) and DIR/cpp_<v>.txt (C++ driver output), evaluates the same float32
orientations (the C++ product local @ centre) with the Python screen (pixel radius 3, centre point,
Q-max 5, C++ overlap rule) and the Python complete cost (radius 0, all reflections, C++ rule), and
writes results.json.gz (every orientation) and summary.json (counts and attributions).

Arms: EXS = exact screen, EXF = exact complete cost (both compared C++ vs Python);
INT = the C++ interpolated screen (compared with C++ EXS to isolate the interpolation, and with
Python EXS for the pairs where the stored traces differ).
"""

import argparse
import json
import sys
import time
from collections import Counter
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import parity_common as PC  # noqa: E402


def PC_pix(x: float) -> int:
    """C++ ToColPixel for a coordinate already in pixels (the half pixel is inside x)."""
    return -1 if x < 0 else int(np.floor(x))


def ulp_diff(a: float, b: float) -> int:
    ia = int(np.float32(a).view(np.int32))
    ib = int(np.float32(b).view(np.int32))
    if ia < 0:
        ia = -(1 << 31) - ia
    if ib < 0:
        ib = -(1 << 31) - ib
    return abs(ia - ib)


def trunc_pixel(v: float) -> int:
    return -1 if v < 0 else int(np.floor(v))


def match_peaks(
    cpp: List[Dict[str, Any]], py: List[Dict[str, Any]], tol: float = 2e-3
) -> Tuple[List[Tuple[Dict[str, Any], Dict[str, Any]]], List[Dict[str, Any]], List[Dict[str, Any]]]:
    """Pair peaks of the two sides by branch and unit direction (the reflection order differs:
    C++ sorts by |G| with an unstable sort, Python keeps the structure order)."""
    pdir = np.array([p["dir"] for p in py]) if py else np.zeros((0, 3))
    used = set()
    pairs = []
    cpp_only = []
    pbr = np.array([p["branch"] for p in py]) if py else np.zeros(0, dtype=int)
    pmag = np.array([p["mag"] for p in py]) if py else np.zeros(0)
    for c in cpp:
        # several reflections of one family share a direction (hkl, 2h2k2l): add |G| to the key
        d = np.linalg.norm(pdir - np.array(c["dir"]), axis=1) if len(pdir) else np.array([])
        d = d + 10.0 * (np.abs(pmag - c["mag"]) > 1e-3) if len(d) else d
        d = np.where(pbr == c["branch"], d, np.inf) if len(d) else d
        hit = None
        if len(d):
            i = int(np.argmin(d))
            if d[i] < tol and i not in used:
                hit = i
        if hit is None:
            cpp_only.append(c)
        else:
            used.add(hit)
            pairs.append((c, py[hit]))
    py_only = [p for i, p in enumerate(py) if i not in used]
    return pairs, cpp_only, py_only


def attach_dirs(cpp_arm: Dict[str, Any]) -> None:
    """NOMATCH peaks carry no direction in the C++ output; nothing to do (kept for clarity)."""
    return


def _counts(x: Dict[str, Any], radius: int) -> Tuple[int, ...]:
    """The counts one peak adds to the cost: peak counts for the screen, plus pixel counts for the
    complete cost (the screen's pixel counts are not used by anything)."""
    g = x.get
    if radius > 0:
        return (g("pkov", 0), g("pkon", 0), g("ndet", 0))
    return (g("pixov", 0), g("pixon", 0), g("pkov", 0), g("pkon", 0), g("ndet", 0))


def _has_flat_lit_detector(c: Dict[str, Any]) -> bool:
    """True if a detector the C++ code calls lit without overlap carries a zero-area triangle
    (integer pixel vertices collinear or equal), i.e. it lit only through the in-bounds vertex rule.
    """
    for d in c["dets"]:
        if d["hit"] == 3 and d["lit"] and not d["ov"]:
            (x0, y0), (x1, y1), (x2, y2) = [tuple(map(float, v)) for v in d["xy"]]
            if (x1 - x0) * (y2 - y0) - (x2 - x0) * (y1 - y0) == 0.0:
                return True
    return False


def compare_peak(
    c: Dict[str, Any], p: Dict[str, Any], radius: int, n_vert: int, low: float
) -> Dict[str, Any]:
    """Per-peak comparison; ``stage`` is the earliest differing stage, '' if all equal.
    ``affects`` is True when the peak's counts differ (the peak changes the cost)."""
    out: Dict[str, Any] = dict(refl=c["refl"], refl_py=p["refl"], branch=c["branch"], stage="")
    out["omega_ulp"] = ulp_diff(c["omega"], p["omega"])
    out["omega_abs"] = abs(c["omega"] - p["omega"])
    if c["nomatch"] or p["nomatch"]:
        if c["nomatch"] != p["nomatch"]:
            # C++ drops an omega below the range (AngleToIndex < 0); the Python batched lookup
            # truncates toward zero and keeps it as frame 0
            out["stage"] = (
                "frame_below_range" if (c["nomatch"] and c["omega"] < low) else "frame_boundary"
            )
            out["affects"] = _counts(c, radius) != _counts(p, radius)
        elif out["omega_ulp"]:
            out["stage"] = "omega_benign"  # both dropped, omega differs
        return out
    if c["frame"] != p["frame"]:
        out["stage"] = "frame_boundary"
        out["affects"] = _counts(c, radius) != _counts(p, radius)
        return out
    # detector hit flags: C++ h == n_vert means all vertices hit the plane
    c_hit = [d["hit"] == n_vert for d in c["dets"]]
    if c_hit != list(p.get("hit", c_hit)):
        out["stage"] = "hit"
        out["affects"] = _counts(c, radius) != _counts(p, radius)
        return out
    # pixel integers
    if radius > 0:
        pc = [(int(d["xy"][0][0]), int(d["xy"][0][1])) if d["hit"] else None for d in c["dets"]]
        pp = [tuple(p["centers"][k]) if c_hit[k] else None for k in range(len(c["dets"]))]
        pp = [None if x is None else (int(x[0]), int(x[1])) for x in pp]
    else:
        pc = [
            tuple((int(v[0]), int(v[1])) for v in d["xy"]) if d["hit"] == n_vert else None
            for d in c["dets"]
        ]
        pp = [
            (
                tuple((trunc_pixel(v[0]), trunc_pixel(v[1])) for v in p["verts"][k])
                if c_hit[k]
                else None
            )
            for k in range(len(c["dets"]))
        ]
    if pc != pp:
        out["stage"] = "pixel"
        diffs = []
        for k in range(len(pc)):
            if pc[k] == pp[k]:
                continue
            a_, b_ = pc[k], pp[k]
            if radius > 0:
                a_, b_ = (a_,), (b_,)  # one point
            if a_ is None or b_ is None:
                diffs.append((k, -1, -1, None, None))
                continue
            for vi, (u, w) in enumerate(zip(a_, b_)):
                for ax in range(2):
                    if u[ax] != w[ax]:
                        diffs.append((k, vi, ax, int(u[ax]), int(w[ax])))
        out["pix_diffs"] = diffs
        out["affects"] = _counts(c, radius) != _counts(p, radius)
        return out
    cc, pc2 = _counts(c, radius), _counts(p, radius)
    if cc != pc2:
        out["stage"] = "counts"
        out["affects"] = True
        out["lit_no_ov"] = any(d["hit"] == n_vert and d["lit"] and not d["ov"] for d in c["dets"])
        out["cpp_pixon_zero"] = c["pixon"] == 0
        out["cnt_cpp"], out["cnt_py"] = cc, pc2
        return out
    if out["omega_ulp"]:
        out["stage"] = "omega_benign"
    return out


def compare_arm(
    cpp_arm: Dict[str, Any], py: Dict[str, Any], radius: int, n_vert: int, low: float
) -> Dict[str, Any]:
    pairs, cpp_only, py_only = match_peaks(cpp_arm["peaks"], py["peaks"])
    # C++ peaks with no hit on any detector are eta-filtered (no pixel) and absent from Python
    cpp_only_live = [c for c in cpp_only if not c["nomatch"] and any(d["hit"] for d in c["dets"])]
    cpp_only_dead = len(cpp_only) - len(cpp_only_live)
    per = [compare_peak(c, p, radius, n_vert, low) for c, p in pairs]
    stages = Counter(x["stage"] for x in per)
    keys = ("pkov", "pkon") if radius > 0 else ("pixov", "pixon", "pkov", "pkon")
    agg_equal = all(cpp_arm[k] == py[k] for k in keys)
    q_close = abs(cpp_arm["quality"] - py["quality"]) < 1e-6
    return dict(
        n_cpp=len(cpp_arm["peaks"]),
        n_py=len(py["peaks"]),
        n_matched=len(pairs),
        cpp_only_live=len(cpp_only_live),
        cpp_only_dead=cpp_only_dead,
        py_only=len(py_only),
        py_only_live=sum(1 for p in py_only if not p["nomatch"]),
        stages=dict(stages),
        counts_equal=agg_equal,
        conf_equal=(cpp_arm["pkov"], cpp_arm["pkon"]) == (py["pkov"], py["pkon"]),
        quality_close=q_close,
        quality_bits_equal=bool(np.float32(cpp_arm["quality"]) == np.float32(py["quality"])),
        cpp=dict(
            pixov=cpp_arm["pixov"],
            pixon=cpp_arm["pixon"],
            pkov=cpp_arm["pkov"],
            pkon=cpp_arm["pkon"],
            quality=cpp_arm["quality"],
            conf=cpp_arm["conf"],
        ),
        py=dict(
            pixov=py["pixov"],
            pixon=py["pixon"],
            pkov=py["pkov"],
            pkon=py["pkon"],
            quality=py["quality"],
            conf=py["conf"],
        ),
        affecting=dict(Counter(x["stage"] for x in per if x.get("affects"))),
        degenerate_lit=sum(
            1
            for (c, p), x in zip(pairs, per)
            if x["stage"] == "counts" and radius == 0 and _has_flat_lit_detector(c)
        ),
        counts_lit_no_ov=sum(1 for x in per if x["stage"] == "counts" and x.get("lit_no_ov")),
        counts_cpp_pixon_zero=sum(
            1 for x in per if x["stage"] == "counts" and x.get("cpp_pixon_zero")
        ),
        counts_pairs=dict(
            Counter(f"{x['cnt_cpp']}|{x['cnt_py']}" for x in per if x["stage"] == "counts")
        ),
        pixel_omega_ulp=[x["omega_ulp"] for x in per if x["stage"] == "pixel"],
        pixel_peaks=[
            dict(refl=x["refl_py"], branch=x["branch"], diffs=x["pix_diffs"], affects=x["affects"])
            for x in per
            if x["stage"] == "pixel"
        ],
        bad_peaks=[x for x in per if x.get("affects")][:6],
        max_omega_ulp=max([x["omega_ulp"] for x in per], default=0),
        omega_abs_max=max([x["omega_abs"] for x in per], default=0.0),
        omega_abs_list=[x["omega_abs"] for x in per] if radius == 0 else [],
        n_omega_bits_equal=sum(1 for x in per if x["omega_ulp"] == 0),
        n_omega_differs=sum(1 for x in per if x["omega_ulp"]),
    )


def interp_vs_exact(int_arm: Dict[str, Any], exs_arm: Dict[str, Any]) -> Dict[str, Any]:
    """C++ interpolated vs C++ exact screen of the same orientation: omega error and peak sets."""
    pairs, a_only, b_only = match_peaks(
        int_arm["peaks"], [dict(x, dir=x.get("dir", [0, 0, 0])) for x in exs_arm["peaks"]]
    )
    err = []
    for i, e in pairs:
        if not i["nomatch"] or not e["nomatch"]:
            d = (i["omega"] - e["omega"] + np.pi) % (2 * np.pi) - np.pi
            err.append(abs(float(np.degrees(d))))
    live_missing = [e for e in b_only if not e["nomatch"] and any(d["hit"] for d in e["dets"])]
    live_extra = [i for i in a_only if not i["nomatch"] and any(d["hit"] for d in i["dets"])]
    return dict(
        omega_err_deg_max=max(err, default=0.0),
        omega_err_deg_med=float(np.median(err)) if err else 0.0,
        n_peaks_int=len(int_arm["peaks"]),
        n_peaks_exs=len(exs_arm["peaks"]),
        exact_peaks_missing_in_int=len(b_only),
        exact_live_missing_in_int=len(live_missing),
        int_peaks_not_in_exact=len(a_only),
        int_live_not_in_exact=len(live_extra),
        same_counts=all(int_arm[k] == exs_arm[k] for k in ("pkov", "pkon")),
        same_pass=(int_arm["pkov"] > 0) == (exs_arm["pkov"] > 0),
        pkov=(int_arm["pkov"], exs_arm["pkov"]),
        nobs=int_arm.get("nobs"),
    )


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--dir", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--voxels", type=int, nargs="*", default=PC.PARITY_VOXELS)
    a = ap.parse_args()
    d, out = Path(a.dir), Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    ctx = PC.PyContext()
    screen, full = ctx.cost_fn("screen"), ctx.cost_fn("full")
    screen_e, full_e = ctx.cost_fn("screen", True), ctx.cost_fn("full", True)
    low = float(ctx.setup.range_map.low)
    results: List[Dict[str, Any]] = []
    t0 = time.time()
    extra: Dict[str, Any] = {}
    for v in a.voxels:
        meta = {m["id"]: m for m in json.load(open(d / f"in_{v}.json"))}
        cpp = PC.parse_cpp(d / f"cpp_{v}.txt")
        verts = ctx.verts(v)
        phase = ctx.mic.voxels[v].phase
        # voxel geometry: C++ vertices/centre against the Python ones
        vv = np.asarray(verts.numpy(), dtype=np.float32)
        cv = np.array([cpp["voxel"][f"vertex{k}"] for k in range(3)], dtype=np.float32)
        cc = np.array(cpp["voxel"]["center"], dtype=np.float32)
        extra[str(v)] = dict(
            vertices_bit_equal=bool(np.array_equal(vv, cv)),
            vertices_max_abs_diff=float(np.abs(vv - cv).max()),
            centre_py_mean=vv.mean(0).tolist(),
            centre_cpp=cc.tolist(),
            centre_max_abs_diff=float(np.abs(vv.mean(0) - cc).max()),
            nrecip_cpp=[cpp["head"]["nrecip_all"], cpp["head"]["nrecip_screen"]],
        )
        for rid, rec in sorted(cpp["recs"].items()):
            m = meta[rid]
            P = rec["prod"].astype(np.float64)
            row: Dict[str, Any] = dict(voxel=v, id=rid, label=m["label"], kind=m["kind"])
            for k in ("deg", "fz", "j", "why", "cpp", "py"):
                if k in m:
                    row["meta_" + k] = m[k]
            # product: C++ float32 product against the Python float64 product cast to float32
            if m["kind"] == 1:
                pp = (
                    np.asarray(m["local"], dtype=np.float64)
                    @ np.asarray(m["centre"], dtype=np.float64)
                ).astype(np.float32)
                row["product_bit_equal"] = bool(np.array_equal(pp, rec["prod"]))
                row["product_max_abs_diff"] = float(np.abs(pp.astype(np.float64) - P).max())
            pys = PC.probe_eval(ctx, screen, P, verts, phase)
            pyf = PC.probe_eval(ctx, full, P, verts, phase)
            row["EXS"] = compare_arm(rec["arms"]["EXS"], pys, PC.SCREEN_RADIUS, 1, low)
            row["EXF"] = compare_arm(rec["arms"]["EXF"], pyf, 0, 3, low)
            pys2 = PC.probe_eval(ctx, screen_e, P, verts, phase)
            pyf2 = PC.probe_eval(ctx, full_e, P, verts, phase)
            row["EXS_edge"] = compare_arm(rec["arms"]["EXS"], pys2, PC.SCREEN_RADIUS, 1, low)
            row["EXF_edge"] = compare_arm(rec["arms"]["EXF"], pyf2, 0, 3, low)
            # which side of the pixel boundary does the double-precision position lie on?
            for arm_key, cf, vv_ in (("EXS_edge", screen_e, "centre"), ("EXF_edge", full_e, "tri")):
                vts = verts.mean(dim=0, keepdim=True).repeat(3, 1) if vv_ == "centre" else verts
                for pk in row[arm_key]["pixel_peaks"]:
                    t64 = PC.pix64(ctx, cf, P, vts, pk["refl"], pk["branch"], phase)
                    det_info = []
                    for k, vi, ax, ci, pi_ in pk["diffs"]:
                        if ci is None:
                            continue
                        x = float(t64[k, vi, ax])
                        det_info.append(
                            dict(
                                x64=x,
                                dist=float(abs(x - round(x))),
                                cpp=ci,
                                py=pi_,
                                cpp_right=bool(PC_pix(x) == ci),
                                py_right=bool(PC_pix(x) == pi_),
                            )
                        )
                    pk["pix64"] = det_info
            row["INT_vs_EXS_cpp"] = interp_vs_exact(rec["arms"]["INT"], rec["arms"]["EXS"])
            ia = rec["arms"]["INT"]
            row["INT"] = dict(
                pkov=ia["pkov"], pkon=ia["pkon"], conf=ia["conf"], nobs=ia.get("nobs")
            )
            row["INT_vs_PyEXS"] = dict(
                same_pass=(ia["pkov"] > 0) == (pys["pkov"] > 0), pkov=(ia["pkov"], pys["pkov"])
            )
            row["INT_vs_PyEXS_edge"] = dict(
                same_pass=(ia["pkov"] > 0) == (pys2["pkov"] > 0), pkov=(ia["pkov"], pys2["pkov"])
            )
            row["cost_cpp_exf"] = 1.0 - rec["arms"]["EXF"]["conf"]
            row["cost_py_exf"] = 1.0 - pyf["conf"]
            results.append(row)
        print(f"voxel {v}: done {len(cpp['recs'])} records, {time.time() - t0:.0f}s", flush=True)
    import gzip

    with gzip.GzipFile(
        filename="", mode="wb", fileobj=open(out / "results.json.gz", "wb"), mtime=0
    ) as gz:
        gz.write(json.dumps(dict(rows=results, voxel_check=extra), default=float).encode())
    print(f"wrote {len(results)} rows")


if __name__ == "__main__":
    main()
