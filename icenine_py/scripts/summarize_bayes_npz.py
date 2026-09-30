"""Summarise an exact-Bayes npz on a test set with error_summary (same output as the GN table).

Usage:
  cd icenine_py
  uv run python scripts/summarize_bayes_npz.py \
      --bayes benchmarks/toy_orientation_stage3/far_test_bayes.npz \
      --test scripts/toy_orientation_stage3_far_test.pt \
      --out benchmarks/toy_orientation_stage3/far_res_bayes.json
"""

import argparse
import json

import numpy as np
import torch

from icenine.orientation_eval import error_summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bayes", required=True, help="Bayes npz (mean, cov, offsets_deg)")
    parser.add_argument("--test", required=True, help="test-set .pt the Bayes run was made on")
    parser.add_argument("--out", required=True, help="output JSON")
    args = parser.parse_args()

    bz = np.load(args.bayes)
    te = torch.load(args.test)
    mags = te["magnitudes_deg"].numpy()
    tr = bz["offsets_deg"]
    assert np.allclose(
        tr, te["offsets_deg"].double().numpy(), atol=1e-6
    ), "Bayes offsets_deg do not match the test set"
    ok = np.isfinite(bz["mean"]).all(1)
    print("bayes finite", ok.sum(), len(ok))
    out = {}
    for m in sorted(set(mags.tolist())):
        k = (mags == m) & ok
        s = error_summary(bz["mean"][k], tr[k])
        fl = float(np.nanmean(np.sqrt(np.trace(bz["cov"][mags == m], axis1=1, axis2=2))))
        out[str(m)] = dict(s, floor=fl)
        print(m, s["median_angle"], s["rms_z"], s["rms_perp"], fl)
    with open(args.out, "w") as f:
        json.dump(out, f, indent=1)


if __name__ == "__main__":
    main()
