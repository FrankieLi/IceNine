#!/usr/bin/env python3
"""
Gauss-Newton centroid baseline (Stage 2) on a frame-coded (observer-rendered) test set.

Fits the orientation offset to each test sample's frame indices and lit-pixel
centroids (icenine/orientation_baselines.py) and saves predictions in the format
train_toy_orientation_nn.py's --extra option reads.

Usage:
  cd icenine_py
  uv run python scripts/gauss_newton_baseline.py --test scripts/toy_orientation_stage1_test.pt \
      --out benchmarks/toy_orientation_stage2/pred_gauss_newton.npz
"""

import argparse
import sys
import time
from pathlib import Path

import numpy as np
import torch

sys.path.insert(0, str(Path(__file__).parent))
from generate_toy_orientation_dataset import DEFAULT_EXAMPLE, build_problem  # noqa: E402


def main():
    from icenine.orientation_baselines import CentroidGaussNewton, extract_measurements
    from icenine.orientation_eval import BatchedObserver, WindowSpec, error_summary

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--test", required=True)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    test_path, out_path = Path(args.test).resolve(), Path(args.out).resolve()  # before chdir
    data = torch.load(test_path)
    if data.get("renderer") != "observer":
        raise SystemExit("needs an observer-rendered (frame-coded) dataset")
    max_q = data.get("max_q", float("nan"))
    problem = build_problem(
        DEFAULT_EXAMPLE,
        data["voxel_index"],
        max_q=None if max_q != max_q else max_q,
        detectors=data.get("detectors", "first"),
        min_sin_eta=float(data.get("min_sin_eta", 0.0)),
    )
    assert len(problem["roi_list"]) == data["n_peaks"], "ROI set differs from the dataset's"
    obs = BatchedObserver(
        problem["R_nom"],
        problem["vertices"],
        problem["sample"],
        problem["detector_list"],
        problem["range_map"],
        problem["exp_setup"],
        problem["roi_list"],
    )
    spec = WindowSpec.from_nominal(obs, data["window_size"], data["frame_half_width"])
    gn = CentroidGaussNewton(obs)

    truth = data["offsets_deg"].double().numpy()
    mags = data["magnitudes_deg"].numpy()
    preds, covs, info = [], [], []
    t0 = time.time()
    for n in range(len(truth)):
        meas = extract_measurements(data["windows"][n], spec, obs)
        r = gn.solve(meas)
        preds.append(r["delta"])
        covs.append(r["cov"])
        info.append((r["n_used"], r["n_iter"], r["converged"], r["chi2"]))
    dt = time.time() - t0
    preds, covs = np.array(preds), np.array(covs)
    info = np.array(info, dtype=float)
    print(
        f"{len(truth)} cases in {dt:.1f}s ({dt / len(truth) * 1e3:.0f} ms each); "
        f"converged {int(info[:, 2].sum())}/{len(truth)}; median spots used {np.median(info[:, 0]):.0f}, "
        f"median iterations {np.median(info[:, 1]):.0f}, median chi2/dof {np.nanmedian(info[:, 3]):.2f}"
    )
    print(
        f"{'|delta|':>8} {'rms_z':>10} {'rms_perp':>10} {'median_ang':>11} {'<0.1deg':>8}   predicted sigma (x,y,z)"
    )
    for mag in sorted(set(mags.tolist())):
        m = mags == mag
        s = error_summary(preds[m], truth[m])
        sd = np.sqrt(np.nanmean(np.diagonal(covs[m], axis1=1, axis2=2), axis=0))
        print(
            f"{mag:8.2f} {s['rms_z']:10.5f} {s['rms_perp']:10.5f} {s['median_angle']:11.5f} "
            f"{s['success_0p1']:8.0%}   {np.round(sd, 5)}"
        )
    out_path.parent.mkdir(parents=True, exist_ok=True)
    np.savez(out_path, pred_deg=preds, cov=covs, truth_deg=truth, magnitudes_deg=mags, info=info)
    print(f"saved {out_path}")


if __name__ == "__main__":
    main()
