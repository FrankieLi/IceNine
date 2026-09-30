"""Summarise the exact Bayes npz of the Stage 3 far-voxel test set with error_summary (writes far_res_bayes.json)."""

import numpy as np, json
from icenine.orientation_eval import error_summary

bz = np.load("benchmarks/toy_orientation_stage3/far_test_bayes.npz")
import torch

te = torch.load("scripts/toy_orientation_stage3_far_test.pt")
mags = te["magnitudes_deg"].numpy()
tr = bz["offsets_deg"]
ok = np.isfinite(bz["mean"]).all(1)
print("bayes finite", ok.sum(), len(ok))
out = {}
for m in sorted(set(mags.tolist())):
    k = (mags == m) & ok
    s = error_summary(bz["mean"][k], tr[k])
    fl = float(np.nanmean(np.sqrt(np.trace(bz["cov"][mags == m], axis1=1, axis2=2))))
    out[str(m)] = dict(s, floor=fl)
    print(m, s["median_angle"], s["rms_z"], s["rms_perp"], fl)
json.dump(out, open("benchmarks/toy_orientation_stage3/far_res_bayes.json", "w"), indent=1)
