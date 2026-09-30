#!/usr/bin/env python3
"""
MC and Riemannian Adam baselines (Stage 2) at the toy test set's perturbation sizes.

Protocol (bench_hp_sweep.py): the experimental images are the Python-simulated
Example2.ThreeVoxels data, generated at the voxel's ground-truth orientation R_nom.
Each optimizer starts from exp([delta]x) R_nom for a test offset delta (same
magnitudes and directions as the toy test set) and searches for the ground truth;
its error is the rotation-vector offset of the result from R_nom. This is the same
local problem as the networks' (data displaced from a known nominal), viewed with the
data fixed and the start displaced, and equivalent to first order.

  mc:   MCOptimizer on the hard (binary pixel overlap) cost, reconstructor defaults
        (3500 steps, 2 restarts), search box 1.5x the perturbation as in the sweep.
  adam: RiemannianAdam (geoopt) on the differentiable cost, scale 2, omega window 1,
        lr 1e-4, 100 steps (the sweep's best setting at 1 deg).

Both use only reflections with |g| <= --max-q (default 8), like the networks.

Usage:
  cd icenine_py
  uv run python scripts/optimizer_baselines.py --test scripts/toy_orientation_stage1_test.pt \
      --out-dir benchmarks/toy_orientation_stage2
"""

import argparse
import math
import sys
import time
from pathlib import Path

import numpy as np
import torch
from scipy.spatial.transform import Rotation

ICENINE_PY = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ICENINE_PY / "benchmarks"))
sys.path.insert(0, str(Path(__file__).parent))

DEG = math.pi / 180.0


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--test", required=True)
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--max-q", type=float, default=8.0)
    parser.add_argument("--which", nargs="+", default=["mc", "adam"], choices=["mc", "adam"])
    parser.add_argument("--mc-steps", type=int, default=3500)
    parser.add_argument("--mc-restarts", type=int, default=2)
    parser.add_argument("--mc-step-frac", type=float, default=0.5)
    parser.add_argument("--adam-lr", type=float, default=1e-4)
    parser.add_argument("--adam-steps", type=int, default=100)
    parser.add_argument("--limit", type=int, default=None)
    parser.add_argument("--seed", type=int, default=0)
    args = parser.parse_args()

    test_path, out_dir = Path(args.test).resolve(), Path(args.out_dir).resolve()  # before chdir
    data = torch.load(test_path)
    offsets = data["offsets_deg"].double().numpy()
    mags = data["magnitudes_deg"].numpy()
    if args.limit:
        offsets, mags = offsets[: args.limit], mags[: args.limit]

    import bench_hp_sweep as hp  # benchmarks/bench_hp_sweep.py
    from icenine.cost_functions import VoxelCostFunction
    from icenine.differentiable_cost import DifferentiableCostFunction

    example = ICENINE_PY.parent / "Examples" / "Example2.ThreeVoxels"
    mic, hard_fn, diff_fn, get_vertices = hp.setup_example(example, "3Grains.sim")
    # Restrict both cost functions to the same reflections as the networks.
    hard_fn = VoxelCostFunction(
        simulator=hard_fn.simulator,
        detector_list=hard_fn.detector_list,
        range_map=hard_fn.range_map,
        exp_data=hard_fn.exp_data,
        sample=hard_fn.sample,
        structure_list=hard_fn.structure_list,
        mode="hard",
        max_q=args.max_q,
    )
    diff_fn = DifferentiableCostFunction(
        simulator=diff_fn.simulator,
        detector_list=diff_fn.detector_list,
        range_map=diff_fn.range_map,
        image_stack=diff_fn.image_stack,
        sample=diff_fn.sample,
        structure_list=diff_fn.structure_list,
        max_q=args.max_q,
    )
    voxel = mic.voxels[int(data["voxel_index"])]
    vertices = get_vertices(voxel)
    R_nom = voxel.orientation.astype(np.float64)
    assert np.allclose(R_nom, data["R_nom"].numpy()), "voxel differs from the dataset's"

    def offset_of(R_final):
        return (
            Rotation.from_matrix(np.asarray(R_final, dtype=np.float64) @ R_nom.T).as_rotvec() / DEG
        )

    out_dir.mkdir(parents=True, exist_ok=True)
    for name in args.which:
        preds, times, quality = [], [], []
        for n, (delta, mag) in enumerate(zip(offsets, mags)):
            R_start = (Rotation.from_rotvec(delta * DEG).as_matrix() @ R_nom).astype(np.float32)
            t0 = time.time()
            if name == "mc":
                np.random.seed(args.seed + n)
                res = hp.run_one_mc(
                    hard_fn,
                    voxel,
                    vertices,
                    R_start,
                    R_nom,
                    args.mc_steps,
                    args.mc_restarts,
                    args.mc_step_frac,
                    mag * DEG,
                    record_traj=False,
                )
            else:
                res = hp.run_one_riemannian_adam_geoopt(
                    diff_fn,
                    voxel,
                    vertices,
                    R_start,
                    hp.FIXED_SCALE,
                    args.adam_steps,
                    args.adam_lr,
                    traj_subsample=0,
                )
            times.append(time.time() - t0)
            preds.append(offset_of(res["R_final"]))
            quality.append(res["final_quality"])
            if (n + 1) % 20 == 0 or n + 1 == len(offsets):
                print(f"  [{name}] {n + 1}/{len(offsets)}  ({np.sum(times):.0f}s)", flush=True)
        errs = np.array(preds)  # final offset from the ground truth (which is at offset 0 here)
        path = out_dir / f"pred_{'mc' if name == 'mc' else 'riemannian_adam'}.npz"
        # pred_deg = the equivalent estimate of the test offset (delta_true + error), so it can be
        # tabulated next to the networks; error_deg is the raw final offset from the ground truth.
        np.savez(
            path,
            pred_deg=offsets + errs,
            error_deg=errs,
            truth_deg=offsets,
            magnitudes_deg=mags,
            seconds=np.array(times),
            quality=np.array(quality),
        )
        err = np.linalg.norm(errs, axis=1)
        print(
            f"{name}: median final error {np.median(err):.4f} deg, success(<0.5deg) {np.mean(err < 0.5):.0%}, "
            f"mean time {np.mean(times):.2f}s  -> {path}"
        )


if __name__ == "__main__":
    main()
