#!/usr/bin/env python3
"""Low-Q feature tables (F-lowQ) for every E2 candidate.

One task = one (voxel, variant): the candidates of cache/e2 (harvested search candidates,
FindOptimal results, synthetic perturbed truths / CSL relatives; R matrices as stored there) get
the features of `FeatureExtractor(q_max=Q)` for Q = 4 and 5 (|q| <= Q Angstrom^-1; Q = 3 holds
no reflection, see the report): one forward pass over the reflections with |q| <= Q plus the
two cost columns at max_q = Q (pixel radius 0 and 3).
Output cache/lowq/v{voxel}_{variant}.npz with X4, X5.

  uv run python scripts/coarse_proxy/lowq_dataset.py run --workers 10
"""

import argparse
import sys
import time
from pathlib import Path
from typing import Any, Dict, Tuple

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import base as B  # noqa: E402

C = B.C
import features as F  # noqa: E402

_FE: Dict[Tuple[float, int], Any] = {}


def task(item: Tuple[Any, ...]) -> str:
    vidx, variant, path = item
    from optimizer_sweep import voxel_context

    W = C.get_worker()
    t0 = time.time()
    keys = np.load(C.CACHE_DIR / "images" / f"v{vidx}_{variant}.npz")["keys"]
    C.attach(keys)
    vctx = voxel_context(W.ctx, vidx)
    vertices, phase = vctx.vertices, vctx.voxel.phase
    R = np.load(C.CACHE_DIR / "e2" / f"v{vidx}_{variant}.npz")["R"]
    out = {}
    for q in B.Q_LEVELS:
        if (q, phase) not in _FE:
            _FE[(q, phase)] = F.FeatureExtractor(W.local_fn, W.ctx.geo, phase, q_max=q)
        fe = _FE[(q, phase)]
        fe.set_image(keys)
        out[f"X{int(q)}"] = np.stack([fe.features(r, vertices, phase) for r in R])
    np.savez_compressed(path, **out)
    return f"voxel {vidx} {variant}: {len(R)} candidates in {time.time() - t0:.0f}s"


def write_qlevels() -> None:
    """Record the |q| families of the structure and which fall under each cut-off Q = 3, 4, 5."""
    from optimizer_sweep import voxel_context

    W = C.get_worker()
    vctx = voxel_context(W.ctx, B.voxels()[0])
    fe = F.FeatureExtractor(W.local_fn, W.ctx.geo, vctx.voxel.phase)
    fam_n = [int((fe.fam == i).sum()) for i in range(fe.n_fam)]
    out = dict(
        q_levels=[float(x) for x in fe.q_levels], reflections_per_family=fam_n,
        n_reflections=int(len(fe.fam)),
        per_cutoff={
            str(q): dict(
                families=[i for i, x in enumerate(fe.q_levels) if x <= q],
                n_reflections=int(sum(n for x, n in zip(fe.q_levels, fam_n) if x <= q)),
            )
            for q in (3, 4, 5)
        },
    )  # fmt: skip
    B.OUT.mkdir(parents=True, exist_ok=True)
    C.save_json(B.OUT / "qlevels.json", out)
    print(out)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    ap.add_argument("cmd", choices=["run", "qlevels"])
    ap.add_argument("--workers", type=int, default=10)
    ap.add_argument("--limit", type=int, default=0)
    a = ap.parse_args()
    if a.cmd == "qlevels":
        C.init_worker(C.worker_args())
        return write_qlevels()
    cache = B.CACHE / "lowq"
    cache.mkdir(parents=True, exist_ok=True)
    its = [
        (v, var, str(cache / f"v{v}_{var}.npz"))
        for v, _, var in B.tasks()
        if not (cache / f"v{v}_{var}.npz").exists()
    ]
    if a.limit:
        its = its[: a.limit]
    print(len(its), "tasks", flush=True)
    C.run_pool(task, its, a.workers, "lowq")


if __name__ == "__main__":
    main()
