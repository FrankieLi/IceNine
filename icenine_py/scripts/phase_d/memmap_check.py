"""Dense vs memmap-uint8 loader: bit-identical hard-cost results and measured RSS per process.

For one image variant: build the uint8 stack (once), then in two fresh subprocesses (dense
`setup_reconstruction(cfg)` vs `setup_reconstruction(cfg, exp_data=from_binary_memmap(...))`)
evaluate the hard cost at the true and a 0.5 deg-off orientation of N voxels and report every
OverlapInfo field, plus the process peak RSS (ru_maxrss) and load time. The parent compares.

Usage (from icenine_py/): uv run python scripts/phase_d/memmap_check.py --variant clean
"""

import argparse
import json
import os
import resource
import subprocess
import sys
import time
from pathlib import Path
from typing import Any, Dict

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
EX = ROOT / "Examples" / "Example2.ManyGrains"
IMG = {"clean": "full/clean", "realistic": "full/realistic", "realistic_q16": "full_q16/realistic"}
STACKS = EX / "ScatteringData_PhaseD" / "stacks"


def child(variant: str, mode: str, n_vox: int, out: str) -> None:
    from scipy.spatial.transform import Rotation

    from icenine.config_file import ConfigFile
    from icenine.experimental_data import ExperimentalData
    from icenine.mic_file import MicFile
    from icenine.reconstructor import (
        AdaptiveVoxelReconstructor,
        _get_voxel_vertices,
        setup_reconstruction,
    )

    os.chdir(EX)
    cfg = ConfigFile.from_file(str(HERE / "configs" / f"ReconstructPhaseD_mc_{variant}.config"))
    t0 = time.time()
    exp = None
    if mode == "memmap":
        exp = ExperimentalData.from_binary_memmap(STACKS / f"{variant}.npy", 180, 2)
    setup = setup_reconstruction(cfg, exp_data=exp)
    t_load = time.time() - t0
    rec = AdaptiveVoxelReconstructor(setup)
    mic = setup.sample.get_mic()
    truth = MicFile.read("SimInput/rand_500grains_1mm_neworient_s0.mic")
    rng = np.random.default_rng(1)
    rows = []
    t1 = time.time()
    for i in rng.choice(len(mic.voxels), n_vox, replace=False):
        v = mic.voxels[int(i)]
        R = np.asarray(truth.voxels[int(i)].orientation, dtype=np.float64)
        ax = rng.normal(size=3)
        ax /= np.linalg.norm(ax)
        R2 = Rotation.from_rotvec(np.radians(0.5) * ax).as_matrix() @ R
        for name, Rx in (("true", R), ("off0.5", R2)):
            oi = rec.evaluate_overlap(Rx.astype(np.float32), _get_voxel_vertices(v), v.phase)
            rows.append({"voxel": int(i), "which": name, **vars(oi)})
    res = {
        "mode": mode,
        "load_s": t_load,
        "eval_s": time.time() - t1,
        "peak_rss_gb": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e9,
        "rows": rows,
    }
    Path(out).write_text(json.dumps(res))


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--variant", default="clean")
    ap.add_argument("--n-vox", type=int, default=20)
    ap.add_argument("--child", default=None)
    ap.add_argument("--out", default=None)
    a = ap.parse_args()
    if a.child:
        child(a.variant, a.child, a.n_vox, a.out)
        return
    sys.path.insert(0, str(ROOT / "icenine_py"))
    from icenine.experimental_data import write_binary_stack

    stack = STACKS / f"{a.variant}.npy"
    t_write = None
    if not stack.exists():
        t0 = time.time()
        write_binary_stack(
            EX / "ScatteringData_PhaseD" / IMG[a.variant],
            "500Grains.sim",
            "d",
            5,
            180,
            2,
            2048,
            2048,
            stack,
        )
        t_write = time.time() - t0
    res: Dict[str, Any] = {}
    for mode in ("dense", "memmap"):
        out = HERE / "cache" / f"memmap_{a.variant}_{mode}.json"
        subprocess.run(
            [
                sys.executable,
                __file__,
                "--variant",
                a.variant,
                "--n-vox",
                str(a.n_vox),
                "--child",
                mode,
                "--out",
                str(out),
            ],
            check=True,
        )
        res[mode] = json.loads(out.read_text())
    identical = res["dense"]["rows"] == res["memmap"]["rows"]
    summ = {
        "variant": a.variant,
        "n_voxels": a.n_vox,
        "n_evaluations": len(res["dense"]["rows"]),
        "all_overlap_fields_identical": identical,
        "stack_write_s": t_write,
        "stack_size_gb": stack.stat().st_size / 1e9,
        **{
            f"{m}_{k}": res[m][k]
            for m in ("dense", "memmap")
            for k in ("load_s", "eval_s", "peak_rss_gb")
        },
    }
    (HERE / "results" / f"memmap_check_{a.variant}.json").write_text(
        json.dumps(summ, indent=2) + "\n"
    )
    print(json.dumps(summ, indent=2))
    assert identical


if __name__ == "__main__":
    main()
