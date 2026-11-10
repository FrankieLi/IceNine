"""Task E measurement 1: the cost along Task B's rays, for the strict cost and each smooth arm.

The rays, voxels and radii are Task B's (``cost_edges_run.py``: six interior full-sample voxels,
0 to 5 deg at 0.0005 deg; the first 30 E0 voxels, 0 to 0.5 deg at 0.0005 deg plus four 0.05-deg
windows and a log profile to 5 deg; 20 directions per voxel, seed 1000 + voxel). The strict cost
along each ray is read from Task B's cache (``full|cpp|cost``: complete cost, C++ rule;
``full|old|cost``: lenient rule); the pilot re-evaluates a few points to check the alignment.

Arms (``icenine.smooth_cost``; fixed widths, since a walk has no search level):
  S            sampled binary overlap (control: the strict cost's pixel test at 4 sample points)
  E1s3, E1s12  distance-transform overlap, kernel width 3 and 12 pixels
  E2           omega blending, kernel reach 1 frame
  E3           denominator stabilisation (soft peak weights, reachable detectors)
  E4t05, E4t25 observed-ray residual, tau 0.05 and 0.25 deg
  E12          E1 (3 px) + E2;  E123  E1 (3 px) + E2 + E3;  E23  E2 + E3;  E34  E3 + E4 (0.25 deg)

  uv run python scripts/eval_reduction/smooth_landscape.py pilot --kind full
  nohup uv run python scripts/eval_reduction/smooth_landscape.py run --kind full --workers 6 &
  nohup uv run python scripts/eval_reduction/smooth_landscape.py run --kind e0 --workers 6 &

Outputs: benchmarks/eval_reduction/E/cache/landscape/<kind>/v<voxel>_d<dir>.npz (gitignored).
"""

import argparse
import multiprocessing as mp
import os
import sys
import time
from pathlib import Path
from typing import Any, Dict, List, Tuple

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")

import numpy as np

HERE = Path(__file__).resolve().parent
ICE = HERE.parents[1]
sys.path.insert(0, str(HERE))
import cost_edges as CE  # noqa: E402
import cost_edges_run as RUN  # noqa: E402

OUT = ICE / "benchmarks" / "eval_reduction" / "E"
CACHE = OUT / "cache" / "landscape"
CACHE_L0 = OUT / "cache" / "landscape_L0"
ARMS: Dict[str, str] = {
    "S": "S",
    "E1s3": "E1,sigma=3",
    "E1s12": "E1,sigma=12",
    "E2": "E2,omega=1",
    "E3": "E3",
    "E4t05": "E4,tau=0.05",
    "E4t25": "E4,tau=0.25",
    "E12": "E1+E2,sigma=3,omega=1",
    "E123": "E1+E2+E3,sigma=3,omega=1",
    "E23": "E2+E3,omega=1",
    "E34": "E3+E4,tau=0.25",
}
# second pass: the "auto" widths of the ranking key at level 0 (diameter 5 deg: kernel 16.5 px,
# omega reach 1.67 frames, tau 0.42 deg), the arms as the search uses them at level 0
ARMS_L0: Dict[str, str] = {
    "E1a": "E1",
    "E2a": "E2",
    "E4a": "E4",
    "E12a": "E1+E2",
    "E123a": "E1+E2+E3",
    "E23a": "E2+E3",
    "E34a": "E3+E4",
    "E1234a": "E1+E2+E3+E4",
}
L0_DIAMETER_DEG = 5.0
_W: Dict[str, Any] = {}


def smooth_cost_for(base: Any, voxel: int) -> Any:
    from icenine.smooth_cost import SmoothCost

    key = ("full",) if _W.get("kind") == "full" else ("e0", voxel)
    if _W.get("sc_key") != key:
        _W["sc"] = SmoothCost(
            base.simulator,
            base.detector_list,
            base.range_map,
            base.exp_data,
            base.sample,
            base.structure_list,
            eta_limit=base.eta_limit,
        )
        _W["sc_key"] = key
    return _W["sc"]


def walk_task(item: Tuple[str, int, List[int]]) -> str:
    from icenine.smooth_cost import SmoothSpec

    kind, voxel, dirs = item
    t0 = time.time()
    base, Rt, verts, phase = RUN.case_data(voxel)
    sc = smooth_cost_for(base, voxel)
    arms = ARMS_L0 if _W.get("set") == "L0" else ARMS
    diam = np.radians(L0_DIAMETER_DEG) if _W.get("set") == "L0" else None
    specs = [SmoothSpec.parse(s) for s in arms.values()]
    cache = CACHE if _W.get("set") != "L0" else CACHE_L0
    all_dirs = CE.unit_directions(RUN.N_DIR, 1000 + voxel)
    done = 0
    for d in dirs:
        path = cache / kind / f"v{voxel}_d{d}.npz"
        if path.exists():
            continue
        out: Dict[str, Any] = {}
        for name, (typ, radii) in RUN.segments(kind).items():
            Rs = CE.rotate_about(Rt, radii, all_dirs[d])
            res = sc.evaluate_specs(Rs, verts, phase, specs, diameter_rad=diam, chunk=512)
            out[f"radii|{name}"] = radii.astype(np.float32)
            for arm, c in zip(arms, res):
                out[f"{arm}|{name}"] = c.astype(np.float32)
        out["n_attached"] = np.int64(len(sc.attached_rays(verts)))
        out["n_rays_obs"] = np.int64(len(sc.observed_rays().frame))
        path.parent.mkdir(parents=True, exist_ok=True)
        np.savez_compressed(path, **out)
        done += 1
    return f"{kind} v{voxel} dirs {dirs[0]}..{dirs[-1]}: {done} written, {time.time() - t0:.0f}s"


def _init(kind: str, arm_set: str = "fixed") -> None:
    import torch

    torch.set_num_threads(1)
    _W["kind"] = kind
    _W["set"] = arm_set
    (RUN.init_full if kind == "full" else RUN.init_e0)()


def items_for(
    kind: str, limit_cases: int = 0, limit_dirs: int = 0
) -> List[Tuple[str, int, List[int]]]:
    vox = RUN.cases(kind)[: limit_cases or None]
    nd = limit_dirs or RUN.N_DIR
    if kind == "full":  # one task per ray (one smooth-cost build per worker)
        return [(kind, v, [d]) for d in range(nd) for v in vox]
    return [(kind, v, list(range(nd))) for v in vox]  # E0: one task per voxel (its own images)


def pilot(kind: str) -> None:
    """One ray of the first voxel, timed, plus a strict-cost alignment check against Task B."""
    _init(kind)
    v = RUN.cases(kind)[0]
    t0 = time.time()
    print(walk_task((kind, v, [0])), f"{time.time() - t0:.1f}s")
    base, Rt, verts, phase = RUN.case_data(v)
    cf = CE.make_cf(base, "full", "cpp")
    zb = np.load(RUN.CACHE / kind / f"v{v}_d0.npz")
    r = zb["radii|dense"].astype(np.float64)
    idx = np.array([0, 10, 100, 500, 900])
    Rs = CE.rotate_about(Rt, r[idx], CE.unit_directions(RUN.N_DIR, 1000 + v)[0])
    mine = np.array([cf.evaluate(R.astype(np.float32), verts, phase).cost for R in Rs])
    theirs = zb["full|cpp|cost|dense"][idx].astype(np.float64)
    print("strict cost now:", np.round(mine, 5), "Task B cache:", np.round(theirs, 5))
    print("max abs difference", float(np.abs(mine - theirs).max()))


def run(kind: str, workers: int, limit_cases: int, limit_dirs: int, arm_set: str) -> None:
    log = OUT / "cache" / "logs"
    log.mkdir(parents=True, exist_ok=True)
    items = items_for(kind, limit_cases, limit_dirs)
    print(f"{kind}: {len(items)} tasks, {workers} workers", flush=True)
    t0 = time.time()
    with mp.get_context("spawn").Pool(workers, initializer=_init, initargs=(kind, arm_set)) as pool:
        for i, msg in enumerate(pool.imap_unordered(walk_task, items)):
            print(f"[{i + 1}/{len(items)} {time.time() - t0:.0f}s] {msg}", flush=True)
    sfx = "" if arm_set == "fixed" else "_L0"
    (log / f"landscape_{kind}{sfx}.done").write_text(f"{time.time() - t0:.0f}s\n")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd", choices=["pilot", "run"])
    ap.add_argument("--kind", default="full", choices=["full", "e0"])
    ap.add_argument("--workers", type=int, default=6)
    ap.add_argument("--cases", type=int, default=0)
    ap.add_argument("--dirs", type=int, default=0)
    ap.add_argument("--set", default="fixed", choices=["fixed", "L0"])
    a = ap.parse_args()
    if a.cmd == "pilot":
        pilot(a.kind)
        return
    run(a.kind, a.workers, a.cases, a.dirs, a.set)


if __name__ == "__main__":
    main()
