"""Task D step 2: geometric spot pairing on the clean full-sample provenance.

For every frame: the detector-0 and detector-1 components, candidate pairs by the line/hull test
(optionally with the ring-angle and pixel-count gates), greedy assignment. Writes the predicted
pairs of each preset to benchmarks/eval_reduction/D/cache/pairs_<preset>.npz. The truth, the
metrics and the figures are in pairing_summary.py.

Usage (from icenine_py/):
  uv run python scripts/eval_reduction/pairing_run.py --presets all
"""

import argparse
import json
import time
from typing import Any, Dict

import numpy as np

import pairing_common as pc
from icenine import spot_pairing as sp

# preset name -> overrides of the base parameters. tight: pixel tolerance only; plan: the arm A
# terms (registration, voxel size x (L1/L0 - 1), arc motion). hull / ring / shape: the gates.
PRESETS: Dict[str, Dict[str, Any]] = {
    "hull_tight": dict(tol="tight"),
    "hull_plan": dict(tol="plan"),
    "ring_tight": dict(tol="tight", ring=True),
    "ring_plan": dict(tol="plan", ring=True),
    "shape_tight": dict(tol="tight", ring=True, shape=1.0),
    "shape_plan": dict(tol="plan", ring=True, shape=1.0),
    # sensitivity to the centroid tolerance (pixels), tight terms otherwise
    "ring_tight_c0.5": dict(tol="tight", ring=True, ctol=0.5),
    "ring_tight_c2": dict(tol="tight", ring=True, ctol=2.0),
    "shape_tight_c0.5": dict(tol="tight", ring=True, shape=1.0, ctol=0.5),
    "shape_tight_c2": dict(tol="tight", ring=True, shape=1.0, ctol=2.0),
}


def make_params(geo: Dict[str, Any], spec: Dict[str, Any]) -> sp.PairingParams:
    om_step = float(np.mean(geo["omega_hi"] - geo["omega_lo"]))
    kw: Dict[str, Any] = dict(
        voxel_size_mm=float(geo["voxel_side"]),
        omega_step=om_step,
        centroid_tol_px=float(spec.get("ctol", 1.0)),
    )
    if spec["tol"] == "plan":
        kw.update(registration_mm=float(geo["registration_mm"]), voxel_term=True, arc_term=True)
    else:
        kw.update(registration_mm=0.0, voxel_term=False, arc_term=True)
    if spec.get("ring"):
        kw["ring_two_theta"] = np.asarray(geo["rings_two_theta"])
    if spec.get("shape"):
        kw["shape_weight"] = float(spec["shape"])
    # a spot of more than 2.5 voxel footprints may be several sources
    side_px = float(geo["voxel_side"]) / float(geo["det0_pw"])
    kw.update(multi_npix=int(2.5 * 0.433 * side_px**2), max_partners=1, multi_partners=3)
    return sp.PairingParams(**kw)


def run_preset(name: str, geo: Dict[str, Any], ct: Dict[str, np.ndarray], centroid: str) -> Dict:
    params = make_params(geo, PRESETS[name])
    gs = pc.detector_geoms(geo)
    hull = sp.SampleHull(geo["voxel_tri"].reshape(-1, 2))
    om_c = 0.5 * (geo["omega_lo"] + geo["omega_hi"])
    ck, rk = ("i_col", "i_row") if centroid == "intensity" else ("c_col", "c_row")
    frame, det = ct["frame"], ct["det"]
    order = np.lexsort((np.arange(len(frame)), det, frame))
    out: Dict[str, list] = {
        k: []
        for k in (
            "c0",
            "c1",
            "cost",
            "chosen",
            "sigma_mm",
            "cond",
            "dist_out",
            "frame",
        )
    }
    ol, os_, dn = [], [], []
    n_cand = np.zeros(len(frame), dtype=np.int64)
    for f in range(len(om_c)):
        sel = {d: np.nonzero((frame == f) & (det == d))[0] for d in (0, 1)}
        comps = []
        for d in (0, 1):
            ix = sel[d]
            comps.append(
                {
                    "c_col": ct[ck][ix],
                    "c_row": ct[rk][ix],
                    "npix": ct["npix"][ix].astype(np.float64),
                    "edge": ct["edge"][ix],
                    "w": (ct["cmax"][ix] - ct["cmin"][ix]).astype(np.float64),
                    "h": (ct["rmax"][ix] - ct["rmin"][ix]).astype(np.float64),
                }
            )
        r = sp.pair_frame(gs, hull, float(om_c[f]), comps[0], comps[1], params)
        out["c0"].append(sel[0][r["i0"].astype(int)])
        out["c1"].append(sel[1][r["i1"].astype(int)])
        for k in ("cost", "chosen", "sigma_mm", "cond", "dist_out"):
            out[k].append(r[k])
        out["frame"].append(np.full(len(r["i0"]), f))
        ol.append(r["origin_lab"])
        os_.append(r["origin_sample"])
        dn.append(r["dirn"])
        n_cand[sel[0]] = r["n_cand_i0"]
    res = {k: np.concatenate(v) for k, v in out.items()}
    res["origin_lab"] = np.concatenate(ol)
    res["origin_sample"] = np.concatenate(os_)
    res["dirn"] = np.concatenate(dn)
    res["n_cand_c0"] = n_cand
    return res


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--presets", default="all")
    ap.add_argument("--tag", default="full")
    ap.add_argument("--centroid", default="intensity", choices=["intensity", "binary"])
    args = ap.parse_args()
    geo = pc.load_geometry()
    prov = pc.load_provenance(args.tag)
    ct = pc.component_table(prov, geo)
    names = list(PRESETS) if args.presets == "all" else args.presets.split(",")
    timing = {}
    for n in names:
        t0 = time.time()
        res = run_preset(n, geo, ct, args.centroid)
        timing[n] = time.time() - t0
        sfx = "" if args.centroid == "intensity" else "_binary"
        np.savez_compressed(pc.CACHE / f"pairs_{args.tag}_{n}{sfx}.npz", **res)
        print(
            f"{n}{sfx}: {len(res['c0'])} candidates, {int(res['chosen'].sum())} chosen, "
            f"{timing[n]:.1f}s",
            flush=True,
        )
    (pc.CACHE / f"pairing_run_timing_{args.tag}.json").write_text(
        json.dumps(
            {"seconds_per_preset_single_process": timing, "centroid": args.centroid}, indent=2
        )
    )


if __name__ == "__main__":
    main()
