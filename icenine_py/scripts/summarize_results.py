#!/usr/bin/env python3
"""Tabulate train_toy_orientation_nn.py --results-json files, averaging seeds.

Usage: summarize_results.py LABEL=glob [LABEL=glob ...] [--variant NAME]
Each glob matches the per-seed json files of one run; values are means over the matching files.
Prints, per group (in-dist / held-out) and metric, the four |delta| bins (0.1/0.25/0.5/1.0),
then corr(voxel error, r_perp) and the median per-voxel error. The 'gn' extra row is shown once.

Reads both layouts of the results json: the current one, {variant: rows}, and the older one
written for a single eval variant, where the rows were stored directly. --variant picks the
test-set corruption (default clean).
"""

import glob
import json
import sys
from typing import Any, Dict, List

import numpy as np

VARIANTS = {"clean", "neighbours", "noise", "all"}


def load(pat: str) -> List[Dict[str, Any]]:
    out = []
    for f in sorted(glob.glob(pat)):
        with open(f) as fh:
            out.append(json.load(fh))
    return out


def rows_of(r: Dict[str, Any], variant: str) -> Dict[str, Any]:
    """Rows of one eval variant from either results-json layout."""
    if set(r) <= VARIANTS:  # {variant: rows}
        return r[variant]
    return r  # older single-variant layout: rows stored directly


def mag_keys(rr: Dict[str, Any], group: str) -> List[str]:
    """The |delta| bin keys of a group as written in the json (e.g. '0.10000000149011612')."""
    keys = [k.split("/", 1)[1] for k in rr if k.startswith(group + "/")]
    return sorted(keys, key=float)


def get(rs: List[Dict[str, Any]], variant: str, group: str, field: str, which: str) -> np.ndarray:
    out = []
    for r in rs:
        r = rows_of(r, variant)
        row = []
        for m in mag_keys(r, group):
            e = r[f"{group}/{m}"]
            if which == "net":
                name = [k for k in e if k.startswith("net")][0]
                row.append(e[name][field] if field != "maha" else e["mean_mahalanobis_sq"])
            else:
                row.append(e[which][field])
        out.append(row)
    return np.mean(out, axis=0)


def fmt(v: np.ndarray) -> str:
    return "/".join(f"{x:.3f}" for x in v)


def main() -> None:
    args = [a for a in sys.argv[1:] if "=" in a]
    variant = "clean"
    if "--variant" in sys.argv:
        variant = sys.argv[sys.argv.index("--variant") + 1]

    print(f"{'run':<26} {'group':<9} {'median angle':<28} {'z RMS':<28} {'perp RMS':<28} Mahal^2")
    gn_done = set()
    for a in args:
        label, pat = a.split("=", 1)
        rs = load(pat)
        if not rs:
            print(label, "no files")
            continue
        for group in ("in-dist", "held-out"):
            ma = get(rs, variant, group, "median_angle", "net")
            z = get(rs, variant, group, "rms_z", "net")
            p = get(rs, variant, group, "rms_perp", "net")
            mh = get(rs, variant, group, "maha", "net")
            name = f"{label} (n={len(rs)})"
            print(f"{name:<26} {group:<9} {fmt(ma):<28} {fmt(z):<28} {fmt(p):<28} {fmt(mh)}")
            rr = rows_of(rs[0], variant)
            if group not in gn_done and "gn" in rr[f"{group}/{mag_keys(rr, group)[1]}"]:
                gn_done.add(group)
                g = [
                    get(rs, variant, group, f, "gn") for f in ("median_angle", "rms_z", "rms_perp")
                ]
                print(
                    f"{'Gauss-Newton':<26} {group:<9} "
                    f"{fmt(g[0]):<28} {fmt(g[1]):<28} {fmt(g[2]):<28}"
                )
        sm = [rows_of(r, variant)["summary"] for r in rs]
        c = np.mean([s["net"]["corr_err_rperp"] for s in sm])
        print(
            f"{'':<26} corr(err, r_perp) {c:+.2f}  per-voxel median err: in-dist "
            f"{np.mean([s['net']['median_voxel_err_in_dist'] for s in sm]):.4f} held-out "
            f"{np.mean([s['net']['median_voxel_err_held_out'] for s in sm]):.4f}"
            + (f"  [GN corr {sm[0]['gn']['corr_err_rperp']:+.2f}]" if "gn" in sm[0] else "")
        )


if __name__ == "__main__":
    main()
