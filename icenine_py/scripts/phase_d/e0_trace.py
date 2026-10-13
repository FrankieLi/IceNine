#!/usr/bin/env python3
"""Per-level Python trace of E0 cases with the C++ overlap rule off and on (existing recorders).

For each (voxel, variant) of the list (default: cases the switch made wrong, i.e. wrong with the
switch on and right in the stored faithful-MC run) and each arm, runs reconstruct_voxel as E0 does
(seed 0, per-voxel images) with `rec.recorder` set and stores, per level, the discrete candidates
(R, score), the post-quick-MC sorted list (R, cost, perm, n_keep), the FindOptimal candidates
(R_in, R_out, cost), the variance result and the final result.

  uv run python scripts/phase_d/e0_trace.py list           # write the case list
  uv run python scripts/phase_d/e0_trace.py run --workers 10
"""

import argparse
import json
import sys
import time
from pathlib import Path
from typing import Any, Dict, List, Tuple

import numpy as np

HERE = Path(__file__).resolve().parent
ICE = HERE.parents[1]
sys.path.insert(0, str(ICE / "scripts" / "findoptimal_robustness"))
sys.path.insert(0, str(ICE / "scripts" / "common"))
import common as C  # noqa: E402

CACHE = HERE / "cache"
BASE = ICE / "scripts" / "cma_finisher" / "cache" / "e0_run_mcfaithful"
LIST = ICE / "benchmarks" / "phase_d_seed_diag" / "level0" / "e0" / "newly_wrong_cases.txt"
TR = CACHE / "e0_trace"


def newly_wrong() -> List[Tuple[int, str]]:
    out = []
    for f in sorted((CACHE / "on").glob("v*.npz")):
        z, b = np.load(f), np.load(BASE / f.name)
        Rt = z["R_true"][None]
        if C.err_deg(z["R"][None], Rt)[0] > C.WRONG_DEG and C.err_deg(b["mc_R"][None], Rt)[0] <= 1:
            v, var = f.stem[1:].split("_", 1)
            out.append((int(v), var))
    return out


def task(item: Tuple[int, int, str, str]) -> str:
    from optimizer_sweep import voxel_context

    vidx, vpos, variant, arm = item
    path = TR / f"v{vidx}_{variant}_{arm}.npz"
    if path.exists():
        return "skip"
    W = C.get_worker()
    t0 = time.time()
    keys = np.load(C.CACHE_DIR / "images" / f"v{vidx}_{variant}.npz")["keys"]
    C.attach(keys)
    vctx = voxel_context(W.ctx, vidx)
    rec = W.rec
    ev: List[Tuple[str, Dict[str, Any]]] = []
    rec.recorder = lambda n, d: ev.append((n, {k: np.array(x) for k, x in d.items()}))
    rec.setup.search_params.cpp_overlap_rule = arm == "on"
    try:
        with C.quiet():
            res = rec.reconstruct_voxel(vctx.vertices, vctx.voxel.phase, rng=C.run_seed(vpos, 0))
    finally:
        rec.recorder = None
        rec.setup.search_params.cpp_overlap_rule = False
    out: Dict[str, Any] = dict(R_true=np.asarray(vctx.R_true, dtype=np.float64))
    out["R_final"] = np.asarray(res.orientation, dtype=np.float64)
    out["cost_final"] = float(res.cost)
    fi = 0
    for name, d in ev:
        if name in ("discrete", "quick_mc"):
            lv = int(d["level"])
            for k, x in d.items():
                if k != "level":
                    out[f"{name}_{lv}_{k}"] = x
        elif name == "find_candidate":
            for k in ("R_in", "R_out", "cost"):
                out[f"find_{fi}_{k}"] = d[k]
            fi += 1
        elif name in ("variance", "final"):
            for k, x in d.items():
                out[f"{name}_{k}"] = x
    out["n_find"] = fi
    np.savez_compressed(path, **out)
    return f"{vidx} {variant} {arm} {time.time() - t0:.0f}s"


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    ap.add_argument("cmd", choices=["list", "run"])
    ap.add_argument("--workers", type=int, default=10)
    a = ap.parse_args()
    if a.cmd == "list":
        cs = newly_wrong()
        LIST.write_text("\n".join(f"{v}:{var}" for v, var in cs) + "\n")
        print(len(cs), "cases ->", LIST)
        return
    info = dict(np.load(C.OUT_DIR / "voxels.npz"))
    vox = [int(v) for v in info["voxel_indices"]][: C.N_VOXELS]
    vpos = {v: i for i, v in enumerate(vox)}
    cs = [tuple(x.split(":")) for x in LIST.read_text().split()]
    TR.mkdir(parents=True, exist_ok=True)
    items = [(int(v), vpos[int(v)], var, arm) for v, var in cs for arm in ("off", "on")]
    t0 = time.time()
    C.run_pool(task, items, a.workers, "e0trace")
    (LIST.parent / "run_meta_trace.json").write_text(
        json.dumps(dict(workers=a.workers, wall_s=time.time() - t0, tasks=len(items),
                        note="contended multi-worker wall, not a timing claim"))
    )  # fmt: skip


if __name__ == "__main__":
    main()
