"""S4 Step 0 reference: the level-0 screen of every orientation and the re-cost of its passers.

For each reference voxel (the 3 ThreeVoxels voxels and 20 seeded random Phase D voxels, both
triangle types): ``peak_overlap`` and ``peak_on_detector`` of all 43,974 level-0 orientations
from ``BatchedScreen`` exactly as ``run_discrete_search_spaced`` runs it, and for every passer
(``peak_overlap > 0``) the re-cost with the level's local cost function (every OverlapInfo
field). An equality reference, not a timing: it runs in worker processes in the background.

Usage (from icenine_py/):
  uv run python scripts/speed_parity/s4/step0_reference.py --set three --workers 1
  uv run python scripts/speed_parity/s4/step0_reference.py --set phase_d --workers 4
Output: benchmarks/s4_translation/cache/ref/<set>_<voxel>.npz (+ ref_index.json).
"""

import argparse
import json
import multiprocessing as mp
import sys
import time
from pathlib import Path
from typing import Any, Dict, List

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import common as C  # noqa: E402

REF_DIR = C.CACHE / "ref"
INFO_FIELDS = (
    "pixel_overlap",
    "pixel_on_detector",
    "peak_overlap",
    "peak_on_detector",
    "quality",
    "n_quality_points",
    "detectors_overlap",
)


def voxel_ids(ds: C.DataSet) -> List[int]:
    return list(range(len(ds.voxels))) if ds.name == "three" else C.reference_voxels(ds)


def run_voxel(ds: C.DataSet, vid: int) -> Dict[str, Any]:
    out = REF_DIR / f"{ds.name}_{vid}.npz"
    if out.exists():
        return dict(vid=vid, skipped=True)
    v = ds.voxels[vid]
    cap = C.capture_level0(ds, v)
    verts = C.voxel_vertices(v)
    t0 = time.time()
    po, pon = C.screen_reference(cap, verts, v.phase)
    t_screen = time.time() - t0
    R = C.all_orientations(cap)
    lcf = cap["local_cost_fn"]
    passers = np.nonzero(po > 0)[0]
    info = np.zeros((len(passers), len(INFO_FIELDS)), np.float64)
    t1 = time.time()
    for i, k in enumerate(passers):
        oi = lcf.evaluate(orientation=R[k], voxel_vertices=verts, phase_index=v.phase)
        info[i] = [float(getattr(oi, f)) for f in INFO_FIELDS]
    t_recost = time.time() - t1
    np.savez_compressed(
        out,
        po=po.astype(np.int16),
        pon=pon.astype(np.int16),
        passers=passers.astype(np.int32),
        recost=info,
        fields=np.array(INFO_FIELDS),
        verts=verts.numpy(),
        points_up=bool(v.points_up),
        position=np.asarray(v.position, float),
        truth=ds.truth[vid] if ds.truth is not None else np.zeros((3, 3)),
    )
    return dict(
        vid=vid,
        n_passers=int(len(passers)),
        t_screen_s=t_screen,
        t_recost_s=t_recost,
        points_up=bool(v.points_up),
    )


def worker(args: Dict[str, Any]) -> List[Dict[str, Any]]:
    ds = C.load_set(args["set"])
    res = []
    for vid in args["vids"]:
        r = run_voxel(ds, vid)
        print(json.dumps(r), flush=True)
        res.append(r)
    return res


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--set", choices=("three", "phase_d"), required=True)
    ap.add_argument("--workers", type=int, default=1)
    ap.add_argument("--limit", type=int, default=0, help="pilot: first N voxels only")
    a = ap.parse_args()
    assert a.workers <= 6
    REF_DIR.mkdir(parents=True, exist_ok=True)
    ds = C.load_set(a.set)
    vids = voxel_ids(ds)
    if a.limit:
        vids = vids[: a.limit]
    del ds
    groups = [vids[i :: a.workers] for i in range(a.workers)]
    t0 = time.time()
    if a.workers == 1:
        res = worker(dict(set=a.set, vids=groups[0]))
    else:
        with mp.get_context("spawn").Pool(a.workers) as pool:
            res = [
                r for g in pool.map(worker, [dict(set=a.set, vids=g) for g in groups]) for r in g
            ]
    idx_path = REF_DIR / f"ref_index_{a.set}.json"
    old = json.loads(idx_path.read_text()) if idx_path.exists() else {}
    for r in res:
        if not r.get("skipped"):
            old[str(r["vid"])] = r
    idx_path.write_text(json.dumps(old, indent=1) + "\n")
    print(f"done {a.set}: {len(res)} voxels in {time.time() - t0:.0f} s (not a timing)")


if __name__ == "__main__":
    main()
