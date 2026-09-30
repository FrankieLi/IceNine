#!/usr/bin/env python3
"""Summarise benchmarks/toy_orientation_arch/diag_multi.json: per r_perp bin and correlations."""
import json
import sys

import numpy as np

rows = json.load(open(sys.argv[1]))
r = np.array([x["r_perp_um"] for x in rows])
bins = [(0, 130), (130, 260), (260, 390), (390, 510)]
for key in ("both", "det0", "det1", "linear_both"):
    print(f"GN {key}: median over voxels of (median angle / z rms / perp rms), by r_perp bin")
    for lo, hi in bins:
        m = (r >= lo) & (r < hi)
        g = [x["gn"][key] for x, k in zip(rows, m) if k]
        print(
            f"  {lo:3d}-{hi:3d} um (n={m.sum():2d}) {np.median([a['median_angle'] for a in g]):.4f} "
            f"{np.median([a['rms_z'] for a in g]):.4f} {np.median([a['rms_perp'] for a in g]):.4f}"
        )
    for f in ("median_angle", "rms_z", "rms_perp"):
        v = np.array([x["gn"][key][f] for x in rows])
        print(f"  corr({f}, r_perp) = {np.corrcoef(r, v)[0, 1]:+.2f}")
print(
    "J^T W J at nominal (config both): median over bins of sigma_z, sigma_perp (deg), cond, z-share"
)
for cfg in ("both", "det0", "det1"):
    print(f" {cfg}")
    for lo, hi in bins:
        m = (r >= lo) & (r < hi)
        i = [x["info"][cfg] for x, k in zip(rows, m) if k]
        sz = np.median([a["sigma_xyz_deg"][2] for a in i])
        sp = np.median([np.hypot(a["sigma_xyz_deg"][0], a["sigma_xyz_deg"][1]) / 2**0.5 for a in i])
        print(
            f"  {lo:3d}-{hi:3d} sigma_z {sz:.4f} sigma_perp {sp:.4f} cond {np.median([a['cond'] for a in i]):.1f}"
            f" zshare {np.median([a['weakest_eigvec_z_share'] for a in i]):.2f}"
        )
    sz = np.array([x["info"][cfg]["sigma_xyz_deg"][2] for x in rows])
    print(f"  corr(sigma_z, r_perp) = {np.corrcoef(r, sz)[0, 1]:+.2f}")
n = np.array([x["n_entries"] for x in rows])
p = np.array([x["n_pairs"] for x in rows])
b = np.array([x["mean_both_present_entries"] for x in rows])
d0 = np.array([x["n_det"][0] for x in rows])
d1 = np.array([x["n_det"][1] for x in rows])
print(
    f"entries {n.min()}-{n.max()} (mean {n.mean():.1f}); det0 {d0.mean():.1f} det1 {d1.mean():.1f}"
)
print(
    f"pairs {p.min()}-{p.max()} (mean {p.mean():.1f}); paired entries = {2*p.mean()/n.mean():.0%}"
)
print(f"mean entries whose partner is also recorded per sample: {b.mean():.1f} of {n.mean():.1f}")
