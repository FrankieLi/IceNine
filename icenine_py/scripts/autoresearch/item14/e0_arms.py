"""Item 14: the Task E E0 harness (smooth_rank_probe.py) with ``FusedCost`` on or off.

E0 = the 200 hard test voxels, one global search each (clean images, seed 0, lenient rule, ``mc``
arm, min_sin_eta 0). The harness is imported, not edited. ``fused_cost`` is set on the worker's
search parameters before each voxel and reset after it, in a pool initializer of its own (spawned
workers do not see a patch made in the parent). The re-ranking arms of the baseline run are off.

  uv run python scripts/autoresearch/item14/e0_arms.py run --spec off --fused --workers 4
  uv run python scripts/autoresearch/item14/e0_arms.py run --spec E3+E4,tau=0.25 --fused
  uv run python scripts/autoresearch/item14/e0_arms.py check --n 20

``run`` writes one npz per voxel under benchmarks/autoresearch/item14/cache/e0/e0_<tag>/
(gitignored) and prints the wrong count (error over 1 degree). ``check`` reruns the first N
voxels of the stored
arms ``off`` and ``E3+E4,tau=0.25`` without the fused key and compares every voxel with the stored
file (evaluations, smooth evaluations, final orientation and cost, exactly); exit 1 on a difference.
Wall times are contended multi-worker and are not timing claims.
"""

import argparse
import json
import multiprocessing as mp
import os
import sys
import time
from pathlib import Path
from typing import Any, Dict, List, Optional

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")

import numpy as np

HERE = Path(__file__).resolve().parent
ICE = HERE.parents[2]
sys.path.insert(0, str(ICE / "scripts" / "eval_reduction"))
sys.path.insert(0, str(ICE / "scripts" / "findoptimal_robustness"))
sys.path.insert(0, str(ICE / "scripts" / "common"))
sys.path.insert(0, str(ICE / "scripts"))

import smooth_rank_probe as P  # noqa: E402

SPEC_S = "E3+E4,tau=0.25"
STORED = P.CACHE  # benchmarks/eval_reduction/E/cache/rank (read only here)
OUT = ICE / "benchmarks" / "autoresearch" / "item14" / "cache" / "e0"
WRONG_DEG = 1.0
_F: Dict[str, bool] = {"fused": False}


def stored_tag(spec: str) -> str:
    """Tag of the stored E0 arm (lenient rule) with this spec."""
    return {"off": "old_base", SPEC_S: "old_e34t25"}[spec]


def arm_tag(spec: str, fused: bool) -> str:
    """Folder tag of an arm of this step: the stored tag, or ``old_off_fused`` /
    ``old_e34t25_fused`` with the fused key."""
    if not fused:
        return stored_tag(spec)
    return {"off": "old_off_fused", SPEC_S: "old_e34t25_fused"}[spec]


def n_wrong(err_deg: np.ndarray) -> int:
    return int((np.asarray(err_deg) > WRONG_DEG).sum())


def diff_voxel(new: Dict[str, Any], stored: Dict[str, Any]) -> Optional[str]:
    """First field in which two voxel results differ (exact equality), or None."""
    for k in ("voxel", "evals", "smooth_evals", "err_deg", "cost_final"):
        if float(new[k]) != float(stored[k]):
            return f"{k}: {float(new[k])!r} vs {float(stored[k])!r}"
    if not np.array_equal(np.asarray(new["R_final"]), np.asarray(stored["R_final"])):
        return "R_final differs"
    return None


def _init(fused: bool) -> None:
    P._init("e0", "")
    P._W["no_rerank"] = True  # the re-ranking arms of the baseline run are off
    _F["fused"] = fused


def task(item: Any) -> str:
    """One voxel through the harness with ``fused_cost`` set for it and reset after."""
    sp = P._W["C"].get_worker().rec.setup.search_params
    sp.fused_cost = _F["fused"]
    try:
        return P.task(item)
    finally:
        sp.fused_cost = False


def run_arm(
    spec: str, fused: bool, workers: int, root: Path, n: int = 0, skip: int = 0
) -> List[Dict[str, Any]]:
    """Run one arm into ``root`` (existing voxel files are kept); return the voxel results."""
    root.mkdir(parents=True, exist_ok=True)
    (root / "spec.txt").write_text(f"{spec}\nfused={fused}\n")
    vox = P.voxel_list("e0", "30")[skip:][: n or None]
    items = [("e0", v, p, "old", spec, str(root / f"v{v}.npz")) for v, p in vox]
    todo = [it for it in items if not Path(it[5]).exists()]
    print(
        f"{root.name} (spec {spec}, fused {fused}): {len(todo)} tasks, {workers} workers",
        flush=True,
    )
    t0 = time.time()
    if workers <= 1:
        _init(fused)
        for i, it in enumerate(todo):
            print(f"[{i + 1}/{len(todo)} {time.time() - t0:.0f}s] {task(it)}", flush=True)
    else:
        with mp.get_context("spawn").Pool(workers, initializer=_init, initargs=(fused,)) as pool:
            for i, msg in enumerate(pool.imap_unordered(task, todo)):
                print(f"[{i + 1}/{len(todo)} {time.time() - t0:.0f}s] {msg}", flush=True)
    (root / "run_meta.json").write_text(
        json.dumps(
            dict(workers=workers, wall_s=time.time() - t0, tasks=len(todo), spec=spec, fused=fused,
                 note="contended multi-worker wall, not a timing claim")
        )  # fmt: skip
        + "\n"
    )
    return [dict(np.load(it[5])) for it in items]


def cmd_run(a: argparse.Namespace) -> int:
    root = OUT / f"e0_{arm_tag(a.spec, a.fused)}{a.suffix}"
    res = run_arm(a.spec, a.fused, a.workers, root, a.limit)
    err = np.array([float(r["err_deg"]) for r in res])
    ev = np.array([int(r["evals"]) for r in res])
    sm = np.array([int(r["smooth_evals"]) for r in res])
    print(
        f"RESULT {root.name}: {len(res)} voxels, wrong {n_wrong(err)}, "
        f"evals/voxel {ev.mean():.0f}, smooth/voxel {sm.mean():.0f}",
        flush=True,
    )
    return 0


def cmd_check(a: argparse.Namespace) -> int:
    bad = 0
    for spec in ("off", SPEC_S):
        root = OUT.parent / "check" / f"e0_{stored_tag(spec)}_n{a.n}"
        res = run_arm(spec, False, a.workers, root, a.n)
        for r in res:
            sp = STORED / f"e0_{stored_tag(spec)}" / f"v{int(r['voxel'])}.npz"
            d = diff_voxel(r, dict(np.load(sp)))
            if d:
                bad += 1
                print(f"CHECK MISMATCH {spec} v{int(r['voxel'])}: {d}", flush=True)
        print(f"CHECK {spec}: {len(res)} voxels, {bad} mismatches so far", flush=True)
    print("CHECK PASSED" if bad == 0 else "CHECK FAILED", flush=True)
    return 1 if bad else 0


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    r = sub.add_parser("run")
    r.add_argument("--spec", required=True, choices=["off", SPEC_S])
    r.add_argument("--fused", action="store_true")
    r.add_argument("--workers", type=int, default=4)
    r.add_argument("--limit", type=int, default=0, help="first N voxels only (0 = all 200)")
    r.add_argument("--suffix", default="", help="folder name suffix (short checks)")
    c = sub.add_parser("check")
    c.add_argument("--n", type=int, default=20)
    c.add_argument("--workers", type=int, default=4)
    a = ap.parse_args()
    return cmd_run(a) if a.cmd == "run" else cmd_check(a)


if __name__ == "__main__":
    sys.exit(main())
