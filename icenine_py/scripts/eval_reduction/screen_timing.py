"""Task A step 3: what the level-0 screen costs per orientation, and how much of it the switches
remove. Single worker, preflight-gated, one voxel, the first N FZ points x the 9 local-grid points.

Arms (C++ overlap rule + centroid screen, CppEdgeRules on in all):
  unbatched   VoxelCostFunction.evaluate per orientation (today)
  batched     ScreenEvaluator.overlaps_exact, 32 FZ points per call
  interp1     ScreenEvaluator.overlaps_interp, one FZ point per call (CppOmegaInterp alone)
  interp32    ScreenEvaluator.overlaps_interp, 32 FZ points per call (both switches)

Per arm: microseconds per orientation in total, and the split into the exact omega solve
(get_scattering_omegas_torch), the C overlap stage (stage_d_overlap) and the rest (stages A to C of
the geometry in torch, the frame lookup, Python glue). The timers wrap those functions and add about
a microsecond per call; each arm is run twice, the first run warms caches and is not reported.

  uv run python scripts/eval_reduction/screen_timing.py \
      --out benchmarks/eval_reduction/A/step3 [--n-fz 600]
"""

import argparse
import json
import sys
import time
from pathlib import Path
from typing import Any, Callable, Dict

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import parity_common as PC  # noqa: E402

sys.path.insert(0, str(PC.ICE / "scripts" / "common"))
import preflight  # noqa: E402


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("--voxel", type=int, default=15901)
    ap.add_argument("--n-fz", type=int, default=600)
    ap.add_argument("--allow-contended", action="store_true")
    a = ap.parse_args()
    out = Path(a.out).resolve()
    out.mkdir(parents=True, exist_ok=True)
    pf = preflight.preflight()
    label = "single-worker"
    try:
        preflight.require_quiet(info=pf)
    except preflight.MachineBusyError as e:
        if not a.allow_contended:
            raise
        label = "contended"
        pf = dict(pf, not_quiet=str(e)[:300])
    (out / "preflight_screen_timing.json").write_text(json.dumps(pf, indent=1))

    ctx = PC.PyContext()
    import torch  # noqa: F401
    from icenine import batched_screen as BS
    from icenine import cost_functions as CF

    cf = ctx.cost_fn("screen", cpp_edge_rules=True)
    fz = np.asarray(ctx.setup.fz_orientations, dtype=np.float64)[: a.n_fz]
    lg = ctx.local_grid()
    verts, phase = ctx.verts(a.voxel), ctx.mic.voxels[a.voxel].phase
    ev = BS.ScreenEvaluator(cf)
    n_orient = len(fz) * len(lg)

    acc = {"solve": 0.0, "stage_d": 0.0}

    def timed(fn: Callable[..., Any], key: str) -> Callable[..., Any]:
        def w(*args: Any, **kw: Any) -> Any:
            t0 = time.perf_counter()
            r = fn(*args, **kw)
            acc[key] += time.perf_counter() - t0
            return r

        return w

    orig = dict(
        solve_cf=CF.get_scattering_omegas_torch,
        solve_fixed=CF.scattering_omegas_fixed_beam,  # what evaluate() solves with since S1
        solve_bs=BS.get_scattering_omegas_torch,
        stage_d=CF._c_stage_d_overlap,
    )

    def arm_unbatched() -> None:
        for i in range(len(fz)):
            for j in range(len(lg)):
                cf.evaluate(lg[j] @ fz[i], verts, phase)

    def arm_batched() -> None:
        for s in range(0, len(fz), 32):
            ch = fz[s : s + 32]
            ev.overlaps_exact(np.matmul(lg[None], ch[:, None]).reshape(-1, 3, 3), verts, phase)

    def arm_interp(step: int) -> Callable[[], None]:
        def f() -> None:
            for s in range(0, len(fz), step):
                ev.overlaps_interp(fz[s : s + step], lg, verts, phase)

        return f

    arms = {
        "unbatched": arm_unbatched,
        "batched": arm_batched,
        "interp1": arm_interp(1),
        "interp32": arm_interp(32),
    }
    res: Dict[str, Any] = dict(
        voxel=a.voxel,
        n_fz=len(fz),
        n_local=len(lg),
        n_orientations=n_orient,
        label=label,
        c_ext=bool(CF._HAS_C_RASTERIZE),
        torch_threads=int(torch.get_num_threads()),
    )
    for name, fn in arms.items():
        fn()  # warm-up run (image binarisation caches, torch)
        CF.get_scattering_omegas_torch = timed(orig["solve_cf"], "solve")
        CF.scattering_omegas_fixed_beam = timed(orig["solve_fixed"], "solve")
        BS.get_scattering_omegas_torch = timed(orig["solve_bs"], "solve")
        CF._c_stage_d_overlap = timed(orig["stage_d"], "stage_d")
        acc["solve"] = acc["stage_d"] = 0.0
        t0 = time.perf_counter()
        fn()
        total = time.perf_counter() - t0
        CF.get_scattering_omegas_torch = orig["solve_cf"]
        CF.scattering_omegas_fixed_beam = orig["solve_fixed"]
        BS.get_scattering_omegas_torch = orig["solve_bs"]
        CF._c_stage_d_overlap = orig["stage_d"]
        # untimed repeat for the headline number (no wrappers)
        t0 = time.perf_counter()
        fn()
        clean = time.perf_counter() - t0
        res[name] = dict(
            total_s=clean,
            us_per_orientation=1e6 * clean / n_orient,
            us_solve=1e6 * acc["solve"] / n_orient,
            us_stage_d=1e6 * acc["stage_d"] / n_orient,
            us_other=1e6 * (total - acc["solve"] - acc["stage_d"]) / n_orient,
            us_per_orientation_with_timers=1e6 * total / n_orient,
        )
        print(name, {k: round(v, 2) for k, v in res[name].items()}, flush=True)
    (out / "screen_timing.json").write_text(json.dumps(res, indent=1) + "\n")


if __name__ == "__main__":
    main()
