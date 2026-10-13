"""Frame-split invariance: rendering ~50 voxels with 1 and with 3 frame-owner workers must give
exactly equal images (same pixels, same float32 intensities). Runs the workers in-process, one
after the other. Writes results/frame_split_check.json.

Usage (from icenine_py/): uv run python scripts/phase_d/frame_split_check.py [--n-voxels 50]
"""

import argparse
import json
import sys
import tempfile
from pathlib import Path
from typing import Any, Dict, Tuple

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import render_full as RF  # noqa: E402


def render(workers: int, voxel_idx: list, out: str) -> Dict[Tuple[int, int], Tuple[Any, ...]]:
    kept: Dict[Tuple[int, int], Tuple[Any, ...]] = {}
    for w in range(workers):
        job = {
            "worker": w,
            "workers": workers,
            "mic": "SimInput/rand_500grains_1mm_neworient_s0.mic",
            "out": out,
            "voxel_idx": voxel_idx,
            "batch": 2000,
            "noise_seed": RF.NOISE_SEED,
            "max_q": 0.0,
            "variants": ["clean", "realistic"],
            "keep": True,
        }
        kept.update(RF._worker(job)["kept"])
    return kept


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-voxels", type=int, default=50)
    a = ap.parse_args()
    raw = np.loadtxt(RF.EX / "SimInput" / "rand_500grains_1mm_neworient_s0.mic", skiprows=1)
    idx = np.argsort(np.linalg.norm(raw[:, :2], axis=1))[: a.n_voxels].tolist()
    with tempfile.TemporaryDirectory() as t1, tempfile.TemporaryDirectory() as t3:
        k1, k3 = render(1, idx, t1), render(3, idx, t3)
        names = sorted(p.name for p in Path(t1, "realistic").iterdir())
        same_files = all(
            (Path(t1, v, n).read_bytes() == Path(t3, v, n).read_bytes())
            for v in ("clean", "realistic")
            for n in names
        )
    same_frames = k1.keys() == k3.keys() and all(
        all(np.array_equal(x, y) for x, y in zip(k1[k], k3[k])) for k in k1
    )
    res = {
        "n_voxels": a.n_voxels,
        "n_frames": len(k1),
        "lit_pixels": int(sum(len(v[0]) for v in k1.values())),
        "arrays_identical_W1_vs_W3": bool(same_frames),
        "written_files_identical_clean_and_realistic": bool(same_files),
    }
    (HERE / "results" / "frame_split_check.json").write_text(json.dumps(res, indent=2) + "\n")
    print(json.dumps(res, indent=2))
    assert same_frames and same_files


if __name__ == "__main__":
    main()
