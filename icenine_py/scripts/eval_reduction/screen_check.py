"""Task A step 3: level-0 screen of one voxel, four ways, against the stored C++ trace.

Arms (all with the C++ overlap rule and the centroid screen, optionally the edge rules):
  unbatched   VoxelCostFunction.evaluate per orientation (today)
  batched     ScreenEvaluator.overlaps_exact, chunks of 32 FZ points
  interp1     ScreenEvaluator.overlaps_interp, one FZ point per call (CppOmegaInterp alone)
  interp32    ScreenEvaluator.overlaps_interp, 32 FZ points per call (both switches)
Reports the per-FZ-point agreement of the pass sets (peak_overlap > 0) with the stored C++ trace
(benchmarks/phase_d_seed_diag/level0/cpp_v15901.json), the agreement of the peak_overlap values,
and, between arms, the number of orientations whose peak_overlap differs.

  uv run python scripts/eval_reduction/screen_check.py --voxel 15901 --out DIR [--edge] [--limit N]
"""

import argparse
import json
import sys
import time
from pathlib import Path
from typing import Any, Dict

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import parity_common as PC  # noqa: E402

TRACE = PC.ICE / "benchmarks" / "phase_d_seed_diag" / "level0" / "cpp_v15901.json"


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--voxel", type=int, default=15901)
    ap.add_argument("--out", required=True)
    ap.add_argument("--edge", action="store_true")
    ap.add_argument("--limit", type=int, default=0, help="first N FZ points only")
    ap.add_argument("--arms", nargs="*", default=["unbatched", "batched", "interp1", "interp32"])
    a = ap.parse_args()
    out = Path(a.out).resolve()
    ctx = PC.PyContext()
    from icenine.batched_screen import ScreenEvaluator

    cf = ctx.cost_fn("screen", cpp_edge_rules=a.edge)
    fz = np.asarray(ctx.setup.fz_orientations, dtype=np.float64)
    if a.limit:
        fz = fz[: a.limit]
    lg = ctx.local_grid()
    verts, phase = ctx.verts(a.voxel), ctx.mic.voxels[a.voxel].phase
    ev = ScreenEvaluator(cf)
    n_fz, n_lg = len(fz), len(lg)
    po: Dict[str, np.ndarray] = {}
    wall: Dict[str, float] = {}
    for arm in a.arms:
        t0 = time.perf_counter()
        if arm == "unbatched":
            m = np.zeros((n_fz, n_lg), np.int64)
            for i in range(n_fz):
                for j in range(n_lg):
                    m[i, j] = cf.evaluate(lg[j] @ fz[i], verts, phase).peak_overlap
        elif arm == "batched":
            m = np.zeros((n_fz, n_lg), np.int64)
            for s in range(0, n_fz, 32):
                ch = fz[s : s + 32]
                prod = np.matmul(lg[None], ch[:, None])
                m[s : s + len(ch)] = ev.overlaps_exact(prod.reshape(-1, 3, 3), verts, phase)[
                    0
                ].reshape(len(ch), n_lg)
        else:
            step = 1 if arm == "interp1" else 32
            m = np.zeros((n_fz, n_lg), np.int64)
            for s in range(0, n_fz, step):
                ch = fz[s : s + step]
                m[s : s + len(ch)] = ev.overlaps_interp(ch, lg, verts, phase)[0]
        wall[arm] = time.perf_counter() - t0
        po[arm] = m
        print(f"{arm}: {wall[arm]:.1f}s", flush=True)

    res: Dict[str, Any] = dict(voxel=a.voxel, edge=a.edge, n_fz=n_fz, n_local=n_lg, wall_s=wall)
    if a.voxel == 15901:
        cpp = json.load(open(TRACE))["fz"][:n_fz]
        cm = np.array([c["peak_overlap"] for c in cpp])
        res["cpp_pass_pairs"] = int((cm > 0).sum())
        for arm, m in po.items():
            same_fz = int(((m > 0) == (cm > 0)).all(axis=1).sum())
            res[arm] = dict(
                pass_pairs=int((m > 0).sum()),
                fz_identical_pass_set=same_fz,
                pairs_pass_differs=int(((m > 0) != (cm > 0)).sum()),
                pairs_pass_cpp_only=int(((cm > 0) & (m == 0)).sum()),
                pairs_pass_py_only=int(((m > 0) & (cm == 0)).sum()),
                overlap_equal_pairs=int((m == cm).sum()),
                overlap_equal_pct=100.0 * float((m == cm).mean()),
            )
    arms = list(po)
    res["between_arms_overlap_differs"] = {
        f"{x}|{y}": int((po[x] != po[y]).sum()) for i, x in enumerate(arms) for y in arms[i + 1 :]
    }
    res["between_arms_pass_differs"] = {
        f"{x}|{y}": int(((po[x] > 0) != (po[y] > 0)).sum())
        for i, x in enumerate(arms)
        for y in arms[i + 1 :]
    }
    out.mkdir(parents=True, exist_ok=True)
    tag = f"v{a.voxel}{'_edge' if a.edge else ''}{'_lim' + str(a.limit) if a.limit else ''}"
    np.savez_compressed(out / f"screen_{tag}.npz", **po)
    (out / f"screen_{tag}.json").write_text(json.dumps(res, indent=1) + "\n")
    print(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()
