"""One Phase D pilot BFS run: one arm (mc / cma / cma_noretry) on one image variant, restricted
to the pilot region (benchmarks/phase_d_pilot/region_grid.mic), single process, memmap images,
BFS rng seed 0. Writes <out>/<arm>_<variant>/{recon.mic, records.json, run.json}.

The configs (configs/ReconstructPhaseD_<arm>_<variant>.config, with BFSRevisitRefit 1) are used
unchanged except for SampleFilename, which is pointed at the region grid in a copy under <out>.

Usage (from icenine_py/):
  uv run python scripts/phase_d/pilot_bfs.py --arm mc --variant clean --out DIR [--max-voxels 50]
"""

import os

# One thread per process (set before numpy / torch are imported). NOTE: the five pilot runs
# (benchmarks/phase_d_pilot) were launched BEFORE this header existed, so they ran with torch's
# default 8 threads (and the launcher without the env settings); run.json of later runs records
# `torch_threads`. A contended A/B on 2026-10-08 (the two arms ran at the same time as each other
# and as the 5 pilot runs; ab_threads.json) found no detectable difference, 178.8 s per seed
# voxel and 0.379 s per local fit for both; their agreement to 0.001 s is unexplained. The setting
# is kept so that the full runs cannot oversubscribe the machine.
for _k in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ[_k] = "1"

import argparse  # noqa: E402
import json  # noqa: E402
import re  # noqa: E402
import time  # noqa: E402
from dataclasses import asdict  # noqa: E402
from pathlib import Path  # noqa: E402

import numpy as np  # noqa: E402
import torch  # noqa: E402

torch.set_num_threads(1)

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
EX = ROOT.parent / "Examples" / "Example2.ManyGrains"
REGION = ROOT / "benchmarks" / "phase_d_pilot"


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--arm", required=True, choices=["mc", "cma", "cma_noretry"])
    ap.add_argument("--variant", required=True, choices=["clean", "realistic", "realistic_q16"])
    ap.add_argument("--out", required=True)
    ap.add_argument("--max-voxels", type=int, default=None)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument(
        "--timing-label",
        default="contended: other jobs on the machine (state the actual number of concurrent runs)",
    )
    a = ap.parse_args()
    out = Path(a.out).resolve() / f"{a.arm}_{a.variant}"
    out.mkdir(parents=True, exist_ok=True)
    text = (HERE / "configs" / f"ReconstructPhaseD_{a.arm}_{a.variant}.config").read_text()
    text, n = re.subn(
        r"(?m)^SampleFilename\s+\S+", f"SampleFilename {(REGION / 'region_grid.mic')}", text
    )
    assert n == 1 and "BFSRevisitRefit 1" in text
    cfg_path = out / "run.config"
    cfg_path.write_text(text)

    os.chdir(EX)
    from icenine.config_file import ConfigFile
    from icenine.experimental_data import ExperimentalData
    from icenine.reconstructor import BFSReconstruction, setup_reconstruction

    cfg = ConfigFile.from_file(str(cfg_path))
    stack = EX / "ScatteringData_PhaseD" / "stacks" / f"{a.variant}.npy"
    t0 = time.time()
    setup = setup_reconstruction(cfg, exp_data=ExperimentalData.from_binary_memmap(stack, 180, 2))
    p = setup.search_params
    t_setup = time.time() - t0
    bfs = BFSReconstruction(setup)
    rng = np.random.default_rng(a.seed)
    t1 = time.time()
    order = bfs.reconstruct_sample(
        output_mic=str(out / "recon.mic"), max_voxels=a.max_voxels, rng=rng
    )
    wall = time.time() - t1
    (out / "records.json").write_text(
        json.dumps({str(k): asdict(v) for k, v in bfs.records.items()}) + "\n"
    )
    info = {
        "arm": a.arm,
        "variant": a.variant,
        "bfs_rng_seed": a.seed,
        "max_voxels": a.max_voxels,
        "n_processed": len(order),
        "local_optimizer": p.local_optimizer,
        "cma_neighbor_max_evals": p.cma_neighbor_max_evals,
        "cma_retry_sigma0_deg": p.cma_retry_sigma0_deg,
        "bfs_revisit_refit": bool(p.bfs_revisit_refit),
        "bfs_revisit_max": p.bfs_revisit_max,
        "wall_total_s": wall,
        "setup_s": t_setup,
        "timing_label": a.timing_label,
        "torch_threads": torch.get_num_threads(),
        "stats": bfs.stats,
    }
    (out / "run.json").write_text(json.dumps(info, indent=2, default=float) + "\n")
    print(json.dumps(info, indent=2, default=float), flush=True)


if __name__ == "__main__":
    main()
