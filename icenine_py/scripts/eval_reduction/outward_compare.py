"""Task F comparison (D-20261021-N1): every carry option on voxels Task F did not use.

Runs, from the same starts, the methods behind the eight options of the comparison:

  f4     F4, 1-frame matching window (library default)
  f4w    F4, 5-frame window
  f5     F5 expanding shells
  f54    F5, then F4 (1-frame) from the point F5 accepted
  f54w   F5, then F4 with the 5-frame window (the fall-back of the 5-frame policy)
  f2_1, f2_2   CMA-ES with sigma0 1 and 2 deg (stops as in Task F: x-tolerance 0.005 deg or 4000)

The policies "F4, then F5+F4 where F4's own acceptance test rejects it" are not run separately: they
are formed per run in ``outward_compare_summary.py`` (F4's record, and F5+F4's where F4 was
rejected; cost = F4's plus F5+F4's), as Task F's ``policy`` table does. f5, f54 and f54w
share one F5 run per start.

Voxels (``sets``; written to benchmarks/eval_reduction/Fcompare/voxel_sets.json):
  e0    24 E0 voxels with clean per-voxel images that Task F did not use (its 24 benchmark and 4
        pilot voxels, F/voxel_sets.json). Only 10 voxels of the 200 have r_perp below 120 um and
        Task F used all 10, so the "near-axis" group is the 8 unused voxels of smallest r_perp
        (122-150 um); the far group (r_perp >= 250 um) is 16 spread evenly over the unused ones.
  full  the 24 voxels of Task E's 30-voxel list (phase_d_seed_diag/ray_overlap/voxels.json) that are
        not the six interior ones; clean full-sample images.
Offsets 1, 2, 3, 5 deg (indices 2..5 of ``outward_bench.OFFSETS``), 20 directions (seed 7000 + voxel
as in Task F), CMA-ES seed from (voxel, offset index, direction) as in Task F.

  uv run python scripts/eval_reduction/outward_compare.py sets
  nohup uv run python scripts/eval_reduction/outward_compare.py run --kind e0 --workers 3 &
  uv run python scripts/eval_reduction/outward_compare.py timing --kind e0   # single worker, gated

Outputs: benchmarks/eval_reduction/Fcompare/cache/<tag>/v<voxel>_o<offset index>.npz (gitignored),
Fcompare/timing/ (single-worker pass and its preflight JSON).
"""

import argparse
import json
import multiprocessing as mp
import sys
import time
from pathlib import Path
from typing import Any, Dict, List, Optional

import numpy as np

HERE = Path(__file__).resolve().parent
ICE = HERE.parents[1]
sys.path.insert(0, str(HERE))
import outward_bench as B  # noqa: E402  (sets OMP/MKL threads to 1 and the sys.path)

OUT = ICE / "benchmarks" / "eval_reduction" / "Fcompare"
F_OUT = ICE / "benchmarks" / "eval_reduction" / "F"
OFFSET_IDX = [2, 3, 4, 5]  # 1, 2, 3, 5 deg
METHODS = ["f2_1", "f2_2", "f4", "f4w", "f5", "f54", "f54w"]
TIMING_N_VOXELS = 2  # per set, the first two of each list
TIMING_N_DIRS = 5  # the first five of the 20 directions


def fresh_sets() -> Dict[str, Any]:
    """The comparison's voxel lists (see the module docstring)."""
    import common as C

    info = dict(np.load(C.OUT_DIR / "voxels.npz"))
    vox = [int(v) for v in info["voxel_indices"]][: C.N_VOXELS]
    rp = np.asarray(info["voxel_r_perp_um"])[: C.N_VOXELS]
    ok = [(C.CACHE_DIR / "images" / f"v{v}_clean.npz").exists() for v in vox]
    used_json = json.loads((F_OUT / "voxel_sets.json").read_text())
    used = set(used_json["near"]) | set(used_json["far"]) | set(used_json["pilot"])
    rest = [(float(r), v) for v, r, o in zip(vox, rp, ok) if o and v not in used]
    n_near_total = sum(1 for v, r, o in zip(vox, rp, ok) if o and r < B.NEAR_UM)
    n_near_unused = sum(1 for r, v in rest if r < B.NEAR_UM)
    near = [v for r, v in sorted(rest)[:8]]
    far_pool = [v for r, v in rest if r >= 250.0]  # in the order of the E0 voxel list
    far_pool = [v for v in vox if v in set(far_pool)]
    idx = np.linspace(0, len(far_pool) - 1, 16)
    far = [far_pool[int(i)] for i in idx]
    r_of = {v: float(r) for v, r in zip(vox, rp)}
    from smooth_rank_probe import SIX, VOXELS_30

    full = [v for v in VOXELS_30 if v not in SIX]
    assert len(full) == 24 and not (set(near + far) & used)
    return dict(
        e0=dict(near=near, far=far),
        full=full,
        notes=dict(
            e0_near_axis_voxels_total_with_images=n_near_total,
            e0_near_axis_voxels_unused_by_F=n_near_unused,
            near_group="8 unused voxels of smallest r_perp (none left below %g um)" % B.NEAR_UM,
            near_r_perp_um={str(v): round(r_of[v], 1) for v in near},
            far_pool_size=len(far_pool),
            far_r_perp_um={str(v): round(r_of[v], 1) for v in far},
        ),
    )


def run_group_c(
    env: Dict[str, Any], oi: int, nd: int, take: int, methods: List[str]
) -> Dict[str, np.ndarray]:
    """One (voxel, offset): ``take`` of the ``nd`` directions, every method from every start."""
    from icenine.outward_search import angle_deg

    offset = B.OFFSETS[oi]
    starts = B.starts_for(env, offset, nd)[:take]
    lf, R_true = env["lf"], env["R_true"]
    rows: Dict[str, List[Dict[str, Any]]] = {m: [] for m in methods}
    for d in range(take):
        R0 = starts[d]
        seed = int(np.random.default_rng([env["voxel"], oi, d]).integers(2**31))
        for m, sg in B.SIGMA0.items():
            if m in methods:
                rows[m].append(B.run_f2(env, R0, sg, seed))
        if "f4" in methods:
            rows["f4"].append(B.run_f4(env, R0))
        if "f4w" in methods:
            rows["f4w"].append(B.run_f4(env, R0, B.F4W_FRAMES))
        if any(m in methods for m in ("f5", "f54", "f54w")):
            r5, res5 = B.run_f5(env, R0)
            if "f5" in methods:
                rows["f5"].append(r5)
            if "f54" in methods:
                rows["f54"].append(B.run_f54(env, R0, r5, res5))
            if "f54w" in methods:
                rows["f54w"].append(B.run_f54(env, R0, r5, res5, window=B.F4W_FRAMES))
    out: Dict[str, np.ndarray] = dict(
        err_start=np.array([angle_deg(s, R_true) for s in starts]),
        cost_start=np.array(
            [float(lf.evaluate(s, env["verts"], env["phase"]).cost) for s in starts]
        ),
        cost_true=np.array(float(lf.evaluate(R_true, env["verts"], env["phase"]).cost)),
        offset=np.array(offset),
    )
    for m, rs in rows.items():
        for k in rs[0]:
            out[f"{m}|{k}"] = np.array([r[k] for r in rs])
    return out


def task_c(item: Any) -> str:
    kind, voxel, tag, nd, take, methods = item
    t0 = time.time()
    done = 0
    for oi in OFFSET_IDX:
        path = OUT / "cache" / tag / f"v{voxel}_o{oi}.npz"
        if path.exists():
            continue
        env = B.voxel_env(voxel)
        res = run_group_c(env, oi, nd, take, methods)
        path.parent.mkdir(parents=True, exist_ok=True)
        np.savez_compressed(path, **res)
        done += 1
    return f"{kind} v{voxel}: {done} offsets written, {time.time() - t0:.0f}s"


def _init(kind: str) -> None:
    B.init_worker(kind)


def run_pool(kind: str, voxels: List[int], tag: str, nd: int, take: int, workers: int,
             methods: List[str]) -> None:  # fmt: skip
    (OUT / "cache" / "logs").mkdir(parents=True, exist_ok=True)
    items = [(kind, v, tag, nd, take, methods) for v in voxels]
    print(f"{tag}: {len(items)} voxels, {workers} workers, {take} of {nd} dirs", flush=True)
    t0 = time.time()
    if workers <= 1:
        _init(kind)
        for i, it in enumerate(items):
            print(f"[{i + 1}/{len(items)} {time.time() - t0:.0f}s] {task_c(it)}", flush=True)
    else:
        with mp.get_context("spawn").Pool(workers, initializer=_init, initargs=(kind,)) as p:
            for i, msg in enumerate(p.imap_unordered(task_c, items)):
                print(f"[{i + 1}/{len(items)} {time.time() - t0:.0f}s] {msg}", flush=True)
    (OUT / "cache" / "logs" / f"{tag}.done").write_text(f"{time.time() - t0:.0f}s\n")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd", choices=["sets", "run", "timing"])
    ap.add_argument("--kind", choices=["e0", "full"], default="e0")
    ap.add_argument("--workers", type=int, default=3)
    ap.add_argument("--dirs", type=int, default=20)
    ap.add_argument("--take", type=int, default=0, help="run only the first N directions")
    ap.add_argument("--tag", default="")
    ap.add_argument("--voxels", default="", help="comma separated; default: the set's list")
    ap.add_argument("--n-voxels", type=int, default=0, help="use only the first N of the list")
    ap.add_argument("--methods", default=",".join(METHODS))
    a = ap.parse_args()
    sets_path = OUT / "voxel_sets.json"
    if a.cmd == "sets":
        s = fresh_sets()
        sets_path.write_text(json.dumps(s, indent=1) + "\n")
        print(json.dumps(s["notes"]))
        print("e0:", s["e0"]["near"], s["e0"]["far"], "full:", s["full"])
        return
    s = json.loads(sets_path.read_text())
    vox = s["full"] if a.kind == "full" else s["e0"]["near"] + s["e0"]["far"]
    if a.voxels:
        vox = [int(x) for x in a.voxels.split(",")]
    if a.n_voxels:
        vox = vox[: a.n_voxels]
    methods = a.methods.split(",")
    take = a.take or a.dirs
    if a.cmd == "timing":
        import preflight

        pf = preflight.preflight()
        (OUT / "timing").mkdir(parents=True, exist_ok=True)
        (OUT / "timing" / f"preflight_{a.kind}.json").write_text(json.dumps(pf, indent=1) + "\n")
        preflight.require_quiet(info=pf)
        if a.kind == "e0":  # one near-group and one far-group voxel
            vox = [s["e0"]["near"][0], s["e0"]["far"][0]]
        else:
            vox = s["full"][:TIMING_N_VOXELS]
        take = a.take or TIMING_N_DIRS
        run_pool(a.kind, vox, a.tag or f"timing_{a.kind}", a.dirs, take, 1, methods)
        return
    run_pool(a.kind, vox, a.tag or a.kind, a.dirs, take, a.workers, methods)


if __name__ == "__main__":
    main()
