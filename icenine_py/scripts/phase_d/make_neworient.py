"""Phase D step 1: new orientations for the 500-grain sample.

Keeps the voxel geometry and the voxel -> grain map of rand_500grains_1mm_inFZ.mic and draws one
new orientation per grain (uniform on SO(3), fixed seed, reduced to the fundamental zone of the
sample symmetry from the simulation config). Writes, next to the original:
  rand_500grains_1mm_neworient_s<seed>.mic, ..._grainmap.npy, ..._stats.json

Usage (from icenine_py/): uv run python scripts/phase_d/make_neworient.py [--seed 0]
"""

import argparse
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
import grains as G  # noqa: E402

from icenine.config_file import ConfigFile  # noqa: E402
from icenine.experiment_setup import XDMExperimentSetup  # noqa: E402
from icenine.mic_file import MicFile  # noqa: E402
from icenine.orientation_search import get_symmetry_quaternions  # noqa: E402

ROOT = Path(__file__).resolve().parents[3]
EX = ROOT / "Examples" / "Example2.ManyGrains"


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--min-sep-deg", type=float, default=1.0)
    ap.add_argument("--config", default=str(EX / "ConfigFiles" / "Example2.Simulation.config"))
    args = ap.parse_args()

    src = EX / "SimInput" / "rand_500grains_1mm_inFZ.mic"
    stem = f"rand_500grains_1mm_neworient_s{args.seed}"
    dst = EX / "SimInput" / f"{stem}.mic"

    raw = np.loadtxt(src, skiprows=1)
    positions, euler = raw[:, :3], raw[:, 6:9]
    grain, uniq_euler = G.group_grains(euler)
    n_grains = len(uniq_euler)

    sym = get_symmetry_quaternions(
        XDMExperimentSetup(ConfigFile.from_file(args.config)).get_sample_symmetry()
    )
    old_q = G.euler_to_quat(uniq_euler)

    # near-duplicate check among the OLD grains (different Euler rows, same orientation up to
    # symmetry or within 0.1 deg)
    d_old_old = G.misorientation_matrix_deg(old_q, old_q, sym)
    np.fill_diagonal(d_old_old, np.inf)
    side = float(open(src).readline().split()[0]) / 2 ** int(raw[0, 4])
    # centroid of each triangle (left vertex + direction flag)
    s3 = np.sqrt(3.0)
    cy = raw[:, 1] + np.where(raw[:, 3] == 1, 1.0, -1.0) * side * s3 / 6.0
    cent = np.stack([raw[:, 0] + side / 2.0, cy], axis=1)
    adjacency = G.grain_adjacency(cent, grain, 1.01 * side)

    new_q, stats = G.draw_separated(
        n_grains, old_q, adjacency, sym, seed=args.seed, min_sep_deg=args.min_sep_deg
    )
    new_euler = np.round(G.quat_to_euler(new_q), G.EULER_DECIMALS)
    G.write_mic_with_euler(str(src), str(dst), grain, new_euler)
    np.save(EX / "SimInput" / f"{stem}_grainmap.npy", grain)

    # round trip through the Python reader
    mic = MicFile.read(str(dst))
    from icenine.geometry import matrix_to_euler

    back = np.array([matrix_to_euler(v.orientation) for v in mic.voxels])
    q_back = G.euler_to_quat(back)
    q_exp = new_q[grain]
    rt = np.degrees(2 * np.arccos(np.clip(np.abs((q_back * q_exp).sum(1)), 0, 1))).max()

    sizes = np.bincount(grain)
    stats.update(
        {
            "seed": args.seed,
            "n_voxels": int(len(grain)),
            "n_grains": int(n_grains),
            "grain_size_min_median_max": [
                int(sizes.min()),
                float(np.median(sizes)),
                int(sizes.max()),
            ],
            "min_old_vs_old_deg": float(d_old_old.min()),
            "n_old_grain_pairs_below_1deg": int(
                (d_old_old[np.triu_indices(n_grains, 1)] < 1.0).sum()
            ),
            "n_old_in_fz": int(
                sum(
                    G.reduce_quats(old_q[i : i + 1], sym)[0] @ old_q[i] > 0.9999
                    for i in range(n_grains)
                )
            ),
            "side_length_m": side,
            "adjacency_radius_over_side": 1.01,
            "roundtrip_max_misorientation_deg": float(rt),
            "mic_symmetry_ops": int(len(sym)),
            "new_euler_range_deg": [float(new_euler.min()), float(new_euler.max())],
        }
    )
    with open(EX / "SimInput" / f"{stem}_stats.json", "w") as f:
        json.dump(stats, f, indent=2)
    print(json.dumps(stats, indent=2))


if __name__ == "__main__":
    main()
