"""Task G3: trace the BFS flips of the SmoothCostRank-only arm to their full searches.

A BFS voxel that one run has wrong and the other right carries, in the run that has it wrong, the
orientation of a neighbouring grain (``g3_summary.py`` part ``trace``). Where that orientation came
from a full search (the voxel's own, or a wrong full search the BFS carried to it), the smooth key
could be the cause. This script re-runs those full searches through Task E's probe
(``smooth_rank_probe.task``: the Task C probe, clean full-sample images, the same Q-max 8 mc
configuration as the BFS, seed 1, ``min_sin_eta`` 0), with the strict ranking (``--spec off``,
re-ranking skipped) and with ``SmoothCostRank E3+E4,tau=0.25``, lenient rule. Task E's rule
(``smooth_rank_summary.trace`` / ``fixed_trace``, D-20261014-N13) is then applied by
``g3_summary.py`` part ``probe``. The probe's random draws are not the BFS's, so a probe run is
the same search on the same images, not a replay of the BFS one.

  uv run python scripts/eval_reduction/g3_trace_probe.py list           # voxels from summary.json
  uv run python scripts/eval_reduction/g3_trace_probe.py run --spec off --tag base --workers 3
  uv run python scripts/eval_reduction/g3_trace_probe.py run --spec "E3+E4,tau=0.25" \
      --tag e34t25 --workers 3

One npz per full-sample voxel under benchmarks/eval_reduction/G/cache/g3/probe/full_old_<tag>/
(gitignored), a ``<tag>.done`` marker under G/cache/logs. Wall times are contended.
"""

import argparse
import json
import multiprocessing as mp
import os
import sys
import time
from pathlib import Path
from typing import Any, List, Tuple

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")

HERE = Path(__file__).resolve().parent
ICE = HERE.parents[1]
sys.path.insert(0, str(HERE))
import smooth_rank_probe as SRP  # noqa: E402

G = ICE / "benchmarks" / "eval_reduction" / "G"
CACHE = G / "cache" / "g3" / "probe"
LOGS = G / "cache" / "logs"
VOXELS = G / "g3" / "trace_probe_voxels.json"
FLIP = "keyoff_rerun->smooth_only"


def voxel_list() -> List[int]:
    """Full-sample voxels of the full searches behind the flips of the SmoothCostRank-only arm
    against the key-off rerun: for a voxel whose own record is a full search, that voxel; for one
    carried from a wrong full search, the wrong full-search voxels (region rows -> full index)."""
    import numpy as np

    s = json.loads((G / "g3" / "summary.json").read_text())["trace"][FLIP]
    idx = np.load(ICE / "benchmarks" / "phase_d_pilot" / "region_index.npy")
    out = []
    for kind in ("broken", "fixed"):
        for x in s[kind]:
            if x["origin"] == "own full search":
                out.append(int(x["full_sample_voxel"]))
            elif x["origin"] == "carried from a wrong full search":
                out += [int(idx[c]) for c in x["carriers_wrong"]]
    return sorted(set(out))


def run_dir(tag: str) -> Path:
    return CACHE / f"full_old_{tag}"


def _init() -> None:
    SRP._init("full")
    SRP._W["no_rerank"] = True


def task(item: Tuple[Any, ...]) -> str:
    t0 = time.time()
    msg = SRP.task(item)
    path = Path(item[-1])
    side = path.with_suffix(".json")
    if path.exists() and not side.exists():
        side.write_text(json.dumps(dict(wall_s=time.time() - t0, label="contended")) + "\n")
    return msg


def run(spec: str, tag: str, workers: int) -> None:
    vox = json.loads(VOXELS.read_text())["voxels"]
    root = run_dir(tag)
    root.mkdir(parents=True, exist_ok=True)
    LOGS.mkdir(parents=True, exist_ok=True)
    items = [
        ("full", v, p, "old", spec, str(root / f"v{v}.npz"))
        for p, v in enumerate(vox)
        if not (root / f"v{v}.npz").exists()
    ]
    print(f"{root.name} ({spec}): {len(items)} tasks, {workers} workers", flush=True)
    t0 = time.time()
    if workers <= 1:
        _init()
        for i, it in enumerate(items):
            print(f"[{i + 1}/{len(items)} {time.time() - t0:.0f}s] {task(it)}", flush=True)
    else:
        with mp.get_context("spawn").Pool(workers, initializer=_init) as pool:
            for i, msg in enumerate(pool.imap_unordered(task, items)):
                print(f"[{i + 1}/{len(items)} {time.time() - t0:.0f}s] {msg}", flush=True)
    (LOGS / f"g3_probe_{tag}.done").write_text(f"{time.time() - t0:.0f}s\n")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd", choices=["list", "run"])
    ap.add_argument("--spec", default="off")
    ap.add_argument("--tag", default="base")
    ap.add_argument("--workers", type=int, default=3)
    a = ap.parse_args()
    if a.cmd == "list":
        vox = voxel_list()
        VOXELS.write_text(json.dumps(dict(flip=FLIP, voxels=vox), indent=1) + "\n")
        print(f"{len(vox)} voxels: {vox}")
        return
    run(a.spec, a.tag, a.workers)


if __name__ == "__main__":
    main()
