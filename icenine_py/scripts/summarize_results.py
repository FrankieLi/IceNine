#!/usr/bin/env python3
"""Tabulate train_toy_orientation_nn.py --results-json files, averaging seeds.

Usage: summarize_results.py LABEL=glob [LABEL=glob ...] [--variant PREFIX]
Each glob matches the per-seed json files of one run; values are means over the matching files.
Prints, per group (in-dist / held-out) and metric, the four |delta| bins (0.1/0.25/0.5/1.0),
then corr(voxel error, r_perp) and the median per-voxel error. The 'gn' extra row is shown once.
"""
import glob
import json
import sys

import numpy as np

args = [a for a in sys.argv[1:] if "=" in a]
prefix = ""
if "--variant" in sys.argv:
    prefix = sys.argv[sys.argv.index("--variant") + 1] + "/"
mags = ["0.10000000149011612", "0.25", "0.5", "1.0"]


def load(pat):
    return [json.load(open(f)) for f in sorted(glob.glob(pat))]


def get(rs, group, key, field, which):
    out = []
    for r in rs:
        r = r.get(prefix.rstrip("/"), r) if prefix else r
        row = []
        for m in mags:
            e = r[f"{group}/{m}"]
            if which == "net":
                name = [k for k in e if k.startswith("net")][0]
                row.append(e[name][field] if field != "maha" else e["mean_mahalanobis_sq"])
            else:
                row.append(e[which][field])
        out.append(row)
    return np.mean(out, axis=0)


def fmt(v):
    return "/".join(f"{x:.3f}" for x in v)


print(f"{'run':<26} {'group':<9} {'median angle':<28} {'z RMS':<28} {'perp RMS':<28} Mahal^2")
gn_done = set()
for a in args:
    label, pat = a.split("=", 1)
    rs = load(pat)
    if not rs:
        print(label, "no files")
        continue
    for group in ("in-dist", "held-out"):
        ma = get(rs, group, None, "median_angle", "net")
        z = get(rs, group, None, "rms_z", "net")
        p = get(rs, group, None, "rms_perp", "net")
        mh = get(rs, group, None, "maha", "net")
        print(
            f"{label + f' (n={len(rs)})':<26} {group:<9} {fmt(ma):<28} {fmt(z):<28} {fmt(p):<28} {fmt(mh)}"
        )
        gkey = (group, pat.split("/")[-1][:0])
        rr = rs[0].get(prefix.rstrip("/"), rs[0]) if prefix else rs[0]
        if group not in gn_done and "gn" in rr[f"{group}/0.25"]:
            gn_done.add(group)
            g = [get(rs, group, None, f, "gn") for f in ("median_angle", "rms_z", "rms_perp")]
            print(
                f"{'Gauss-Newton':<26} {group:<9} {fmt(g[0]):<28} {fmt(g[1]):<28} {fmt(g[2]):<28}"
            )
    sm = [(r.get(prefix.rstrip("/"), r) if prefix else r)["summary"] for r in rs]
    c = np.mean([s["net"]["corr_err_rperp"] for s in sm])
    print(
        f"{'':<26} corr(err, r_perp) {c:+.2f}  per-voxel median err: in-dist "
        f"{np.mean([s['net']['median_voxel_err_in_dist'] for s in sm]):.4f} held-out "
        f"{np.mean([s['net']['median_voxel_err_held_out'] for s in sm]):.4f}"
        + (f"  [GN corr {sm[0]['gn']['corr_err_rperp']:+.2f}]" if "gn" in sm[0] else "")
    )
