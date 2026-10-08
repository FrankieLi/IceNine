"""Phase D step 3: full-sample forward simulation of the 500-grain sample (new orientations).

Renders both detectors over the whole omega range with the Python forward model (batched path,
exactly the Example2.Simulation.config geometry/beam/omega/Q-max), writes the clean images and a
realistic copy (noise.py, detector noise only) as ASCII 'j, k, intensity' files that both the
Python ExperimentalData loader and the C++ ASCII reader take.

Work is split by omega frame over worker processes: every worker runs the (cheap) per-voxel
physics for all voxels and rasterizes only the frames it owns (the others are null sinks), so no
worker holds more than ~1/W of the 360 dense images and nothing needs merging.

Usage (from icenine_py/):
  uv run python scripts/phase_d/render_full.py --workers 10 [--n-voxels 500] [--tag pilot]
"""

import argparse
import json
import multiprocessing as mp
import os
import resource
import sys
import time
from pathlib import Path
from typing import Any, Dict, List

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent / "common"))

ROOT = HERE.parents[2]
EX = ROOT / "Examples" / "Example2.ManyGrains"
NOISE_SEED = 12345
BASENAME = "500Grains.sim"


class _NullImage:
    """Sink for frames owned by another worker."""

    def add_triangle_scanline(self, *a: Any, **k: Any) -> None:
        return None


def _worker(args: Dict[str, Any]) -> Dict[str, Any]:
    import torch

    torch.set_num_threads(1)
    os.chdir(EX)
    from icenine.config_file import ConfigFile
    from icenine.forward_simulation import ForwardSimulation
    from icenine.image_data import ImageData
    from icenine.sample import Sample
    from noise import NoiseParams, add_detector_noise, write_ascii_image

    t0 = time.time()
    w, W = args["worker"], args["workers"]
    cfg = ConfigFile.from_file(str(EX / "ConfigFiles" / "Example2.Simulation.config"))
    cfg.sample_filename = args["mic"]
    sim = ForwardSimulation(cfg)
    sim.exp_setup.initialize_experiment()
    omega_ranges = sim.exp_setup.get_omega_range_list()
    file_ranges = sim.exp_setup.get_file_range_list()
    dets = sim.exp_setup.get_detector_list()
    range_map = sim.exp_setup.get_range_to_index_map()
    from icenine.simulation import Simulation

    sim.simulator = Simulation(sim.exp_setup)
    sample = Sample()
    sim.exp_setup.initialize_sample(sample, dets[0])
    mic = sample.get_mic()
    if args["voxel_idx"] is not None:
        mic.voxels = [mic.voxels[i] for i in args["voxel_idx"]]
    n_vox = len(mic.voxels)

    n_om = len(omega_ranges)
    owned = [i for i in range(n_om) if i % W == w]
    owned_set = set(owned)
    images = [
        [ImageData(d.num_rows, d.num_cols) if i in owned_set else _NullImage() for d in dets]
        for i in range(n_om)
    ]
    t_setup = time.time() - t0
    t1 = time.time()
    sim._simulate_peaks_batched(images, dets, sample, range_map, batch_size=args["batch"])
    t_sim = time.time() - t1

    t2 = time.time()
    out_clean, out_real = Path(args["out"]) / "clean", Path(args["out"]) / "realistic"
    out_clean.mkdir(parents=True, exist_ok=True)
    out_real.mkdir(parents=True, exist_ok=True)
    params = NoiseParams()
    frames: List[Dict[str, Any]] = []
    for di in range(len(dets)):
        for i in owned:
            img = images[i][di]._pixels_dense.numpy()
            name = f"{BASENAME}{str(file_ranges[di].low + i).zfill(5)}" f".d{di}"
            n_clean = write_ascii_image(str(out_clean / name), img)
            rng = np.random.default_rng([args["noise_seed"], di, i])
            noisy, cnt = add_detector_noise(img, rng, params)
            n_real = write_ascii_image(str(out_real / name), noisy)
            frames.append(
                {"det": di, "omega_idx": i, "lit_clean": n_clean, "lit_real": n_real, **cnt}
            )
            images[i][di] = None  # type: ignore[assignment]
    t_write = time.time() - t2
    rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e9  # GB (macOS reports bytes)
    return {
        "worker": w,
        "n_voxels": n_vox,
        "n_frames": len(owned),
        "t_setup_s": t_setup,
        "t_sim_s": t_sim,
        "t_write_noise_s": t_write,
        "peak_rss_gb": rss,
        "frames": frames,
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=10)
    ap.add_argument("--seed", type=int, default=0, help="new-orientation sample seed")
    ap.add_argument("--n-voxels", type=int, default=0, help="pilot: voxels nearest the centre")
    ap.add_argument("--tag", default="full")
    ap.add_argument("--batch", type=int, default=2000)
    ap.add_argument("--noise-seed", type=int, default=NOISE_SEED)
    ap.add_argument("--out-root", default=str(EX / "ScatteringData_PhaseD"))
    args = ap.parse_args()
    assert args.workers <= 10

    sys.path.insert(0, str(ROOT / "icenine_py"))
    import preflight

    res_dir = HERE / "results"
    res_dir.mkdir(exist_ok=True)
    info = preflight.preflight()
    (res_dir / f"preflight_render_{args.tag}.json").write_text(json.dumps(info, indent=2) + "\n")
    try:
        preflight.require_quiet(info=info)
        quiet = True
    except preflight.MachineBusyError as e:
        print("WARNING:", e)
        quiet = False

    mic_rel = f"SimInput/rand_500grains_1mm_neworient_s{args.seed}.mic"
    voxel_idx = None
    if args.n_voxels:
        raw = np.loadtxt(EX / mic_rel, skiprows=1)
        d = np.linalg.norm(raw[:, :2] - np.array([0.0, 0.0]), axis=1)
        voxel_idx = np.argsort(d)[: args.n_voxels].tolist()
    out = Path(args.out_root) / args.tag
    jobs = [
        {
            "worker": w,
            "workers": args.workers,
            "mic": mic_rel,
            "out": str(out),
            "voxel_idx": voxel_idx,
            "batch": args.batch,
            "noise_seed": args.noise_seed,
        }
        for w in range(args.workers)
    ]
    t0 = time.time()
    ctx = mp.get_context("spawn")
    with ctx.Pool(args.workers) as pool:
        results = pool.map(_worker, jobs)
    wall = time.time() - t0
    frames = [f for r in results for f in r["frames"]]
    summary = {
        "tag": args.tag,
        "sample_seed": args.seed,
        "noise_seed": args.noise_seed,
        "n_voxels": results[0]["n_voxels"],
        "workers": args.workers,
        "timing_label": f"contended ({args.workers} workers; machine quiet at start: {quiet})",
        "wall_s": wall,
        "worker_peak_rss_gb_max": max(r["peak_rss_gb"] for r in results),
        "worker_peak_rss_gb_sum": sum(r["peak_rss_gb"] for r in results),
        "worker_t_setup_s_max": max(r["t_setup_s"] for r in results),
        "worker_t_sim_s_max": max(r["t_sim_s"] for r in results),
        "worker_t_write_noise_s_max": max(r["t_write_noise_s"] for r in results),
        "n_frames": len(frames),
        "lit_clean_per_frame_median": float(np.median([f["lit_clean"] for f in frames])),
        "lit_clean_per_frame_min_max": [
            min(f["lit_clean"] for f in frames),
            max(f["lit_clean"] for f in frames),
        ],
        "lit_real_total": int(sum(f["lit_real"] for f in frames)),
        "lit_clean_total": int(sum(f["lit_clean"] for f in frames)),
        "spots_total": int(sum(f["n_spots"] for f in frames)),
        "missed_total": int(sum(f["n_missed"] for f in frames)),
        "hot_total": int(sum(f["n_hot"] for f in frames)),
        "blob_total": int(sum(f["n_blob"] for f in frames)),
        "output_dir": str(out.relative_to(ROOT)),
    }
    (res_dir / f"render_{args.tag}.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
