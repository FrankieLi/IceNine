#!/usr/bin/env python3
"""
Exact Bayes baseline for a toy-NN test set.

For every test offset, computes the posterior over the orientation offset given
the noise-free thresholded data (the prior restricted to the set of offsets that
reproduce exactly the same observed peaks, frames and lit pixels), by importance
sampling with the batched observer. Its mean is the best possible estimate under
squared error and sqrt(trace(cov)) is the error floor for that test case.

Usage:
  cd icenine_py
  uv run python scripts/exact_bayes_baseline.py --test scripts/toy_orientation_stage0_test.pt
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
    from icenine.orientation_eval import BatchedObserver, ExactBayes

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--test", required=True, help="Test .pt from generate_toy_orientation_dataset.py"
    )
    parser.add_argument("--out", default=None, help="Output .npz (default: <test>_bayes.npz)")
    parser.add_argument(
        "--frames-only", action="store_true", help="Ignore pixel data (frames and presence only)"
    )
    parser.add_argument(
        "--limit", type=int, default=None, help="Only the first N test cases (debugging)"
    )
    parser.add_argument("--seed", type=int, default=7)
    args = parser.parse_args()

    # build_problem() changes directory, so resolve every path first
    test_path = Path(args.test).resolve()
    out_arg = Path(args.out).resolve() if args.out else None
    data = torch.load(test_path)
    offsets = data["offsets_deg"].double().numpy()
    if args.limit:
        offsets = offsets[: args.limit]
    problem = build_problem(DEFAULT_EXAMPLE, data["voxel_index"])
    assert len(problem["roi_list"]) == data["n_peaks"], "ROI set differs from the dataset's"
    assert np.allclose(
        problem["R_nom"], data["R_nom"].numpy()
    ), "nominal orientation differs from the dataset's"

    observer = BatchedObserver(
        problem["R_nom"],
        problem["vertices"],
        problem["sample"],
        problem["detector_list"],
        problem["range_map"],
        problem["exp_setup"],
        problem["roi_list"],
    )
    bayes = ExactBayes(
        observer, prior_radius_deg=float(data["prior_radius_deg"]), use_pixels=not args.frames_only
    )
    rng = np.random.default_rng(args.seed)

    means, covs, ess, n_present = [], [], [], []
    t0 = time.time()
    for i, d in enumerate(offsets):
        r = bayes.posterior(d, rng)
        means.append(r["mean"])
        covs.append(r["cov"])
        ess.append(r["ess"])
        n_present.append(r["n_present"])
        if (i + 1) % 10 == 0 or i + 1 == len(offsets):
            print(
                f"  {i + 1}/{len(offsets)}  ({time.time() - t0:.0f}s)  last: ess {r['ess']:.0f}, "
                f"sqrt(tr cov) {np.sqrt(np.trace(r['cov'])):.2e} deg"
            )

    out = (
        out_arg
        if out_arg
        else Path(
            str(test_path.with_suffix(""))
            + ("_bayes_frames.npz" if args.frames_only else "_bayes.npz")
        )
    )
    np.savez(
        out,
        mean=np.array(means),
        cov=np.array(covs),
        ess=np.array(ess),
        n_present=np.array(n_present),
        offsets_deg=offsets,
        use_pixels=not args.frames_only,
    )
    print(f"Saved {out}")
    n_failed = int(np.sum(~np.isfinite(np.array(means)).all(axis=1)))
    if n_failed:
        print(
            f"WARNING: the sampler found no consistent offsets for {n_failed} case(s); their rows are NaN"
        )
    if min(ess) < 100:
        print(
            f"WARNING: {int(np.sum(np.array(ess) < 100))} cases have ESS < 100; treat those as unreliable"
        )


if __name__ == "__main__":
    main()
