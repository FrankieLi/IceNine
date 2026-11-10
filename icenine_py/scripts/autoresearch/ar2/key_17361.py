"""AR2 step 2, named test case: voxel 17361 (full-sample index), outside the pilot region.

With the library path (``AdaptiveVoxelReconstructor.index_seed``), record what the finisher does
from the truth candidate (about 0.045 degrees from the truth, strict cost 0.162; the strict cost at
the exact truth is 0.0949; the finisher did not move in 924 evaluations in the holdout): which
finisher stage ran, how many steps were accepted, and whether a restart of the finisher from its own
end point, or the ray Gauss-Newton step (``RayGaussNewton``) from that point, reaches a cost near
0.095. A diagnosis, not a new arm: nothing in the library changes because of it.

  uv run python scripts/autoresearch/ar2/key_17361.py
Writes benchmarks/autoresearch/ar2/key/voxel_17361.json.
"""

import json
import math
import sys
import time
from pathlib import Path
from typing import Any, Dict

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[0] / "ar1"))
import ar1_common as AC  # noqa: E402

from icenine import ray_indexing as RI  # noqa: E402
from icenine.orientation_search import SearchCandidate  # noqa: E402

OUT = AC.ICE / "benchmarks" / "autoresearch" / "ar2" / "key"
VOXEL = 17361


def ang(a: np.ndarray, b: np.ndarray, sym: np.ndarray) -> float:
    return float(RI.misorientation_deg(np.asarray(a)[None], np.asarray(b)[None], sym)[0])


def main() -> None:
    AC.init("full")
    e = AC.env("full", VOXEL, need_sc=False)
    rec, verts, phase, Rt = e["rec"], e["verts"], e["phase"], e["R_true"]
    sym = RI.proper_rotations(AC.symmetry_of(e))
    t0 = time.perf_counter()
    lib = rec.index_seed(verts, phase)  # the library path (builds the SmoothCost, pairs the rays)
    out: Dict[str, Any] = dict(
        voxel=VOXEL, label="contended", library_seconds=time.perf_counter() - t0
    )
    d = RI.misorientation_deg(lib.candidates, Rt[None], sym)
    jt = int(np.nonzero(d < 2.0)[0][0])
    lf = rec._make_local_cost_fn()

    def strict(R: np.ndarray) -> float:
        return float(lf.evaluate(np.asarray(R, np.float64), verts, phase).cost)

    out.update(
        n_cand=int(lib.n_cand),
        truth_rank=jt + 1,
        truth_start_deg=float(d[jt]),
        strict_costs_first5=lib.strict_costs,
        refined=lib.refined,
        finisher_costs=lib.finisher_costs,
        finisher_evals=lib.finisher_evals,
        taken=lib.taken,
        taken_err_deg=ang(lib.orientation, Rt, sym),
        cost_at_exact_truth=strict(Rt),
        truth_candidate_cost=strict(lib.candidates[jt]),
    )
    # the finisher from the truth candidate, stage by stage (as refine_from_candidates builds it)
    from icenine.orientation_search import MCOptimizer  # noqa: F401

    params = rec.params
    lcf = rec._make_local_cost_fn()
    mc, find = rec._make_optimizers(lcf, verts, phase, np.random.default_rng(0))
    out["hybrid_optimizer"] = find is not None
    out["local_optimizer"] = params.local_optimizer
    diameter = params.local_grid_radius / 1.5 ** (params.max_local_resolution + 1)
    final_radius = max(diameter / 3.0, math.radians(0.2))
    box = final_radius / (2**params.min_local_resolution)
    step = box * params.mc_radius_scale_factor
    R0 = lib.candidates[jt].copy()
    n0 = lcf.eval_count
    r1 = mc.optimize(R0, box, step, params.max_mc_steps, params.successive_restarts,
                     params.max_convergence_cost)  # fmt: skip
    run1 = dict(mc.last_run)
    n1 = lcf.eval_count
    r2 = mc.variance_minimizing_optimize(
        r1.orientation, box, params.max_mc_steps, params.successive_restarts, 0.0, 0.02**2
    )
    n2 = lcf.eval_count
    out["stages"] = dict(
        box_deg=math.degrees(box),
        step_deg=math.degrees(step),
        findoptimal=dict(
            ran=True,
            evals=int(n1 - n0),
            cost_before=float(run1["cost_start"]),
            cost_after=float(r1.cost),
            moved_deg=ang(r1.orientation, R0, sym),
            **{
                k: (float(v) if isinstance(v, (int, float, np.floating)) else v)
                for k, v in run1.items()
            },
        ),  # fmt: skip
        variance_minimizing=dict(
            ran=True,
            evals=int(n2 - n1),
            cost_before=float(r1.cost),
            cost_after=float(r2.cost),
            moved_deg=ang(r2.orientation, r1.orientation, sym),
        ),
    )
    best = r2 if r2.cost < r1.cost else r1
    out["replica_end_cost"] = float(best.cost)
    out["replica_matches_library"] = bool(
        abs(best.cost - lib.finisher_costs[lib.refined.index(jt)]) < 1e-12
        if jt in lib.refined
        else False
    )
    out["replica_end_moved_from_start_deg"] = ang(best.orientation, R0, sym)
    # restart of the finisher from its own end point (a fresh random stream, as the library's)
    Rend = np.asarray(best.orientation, dtype=np.float64)
    rr = rec.refine_from_candidates(
        [SearchCandidate(orientation=Rend.copy(), cost=float(best.cost))],
        verts, phase, rng=np.random.default_rng(0),
    )  # fmt: skip
    out["restart_from_end"] = dict(
        cost=float(rr.cost), evals=int(rec.last_eval_counts[1]),
        moved_deg=ang(rr.orientation, Rend, sym), err_deg=ang(rr.orientation, Rt, sym),
    )  # fmt: skip
    # the ray Gauss-Newton step from the end point (Task F's F4 settings, as OutwardSearch)
    from icenine.outward_search import RayGaussNewton
    from icenine.reconstructor import (
        OUTWARD_ACCEPT_COST,
        OUTWARD_FRAME_WINDOW,
        OUTWARD_TAU_DEG,
    )

    _, _, sc = rec._outward_costs()
    gn = RayGaussNewton(
        sc, verts, strict, phase, tau_schedule_deg=OUTWARD_TAU_DEG,
        accept_cost=OUTWARD_ACCEPT_COST, frame_window=OUTWARD_FRAME_WINDOW,
    )  # fmt: skip
    t1 = time.perf_counter()
    g = gn.solve(Rend)
    out["ray_gauss_newton_from_end"] = dict(
        cost=float(g.strict_cost), accepted=bool(g.accepted), n_strict=int(g.n_strict),
        n_smooth=int(g.n_smooth), n_jac=int(g.n_jac), moved_deg=ang(g.R, Rend, sym),
        err_deg=ang(g.R, Rt, sym), seconds=time.perf_counter() - t1,
    )  # fmt: skip
    # and from the exact truth, to show what the cost does around it
    from scipy.spatial.transform import Rotation

    out["strict_cost_at_truth_perturbed"] = {
        f"{a}deg": strict(
            Rotation.from_rotvec(np.radians(a) * np.array([1.0, 0, 0])).as_matrix() @ Rt
        )
        for a in (0.02, 0.045, 0.1)
    }
    (OUT / "voxel_17361.json").write_text(json.dumps(out, indent=1, default=float) + "\n")
    print(json.dumps(out, indent=1, default=float))


if __name__ == "__main__":
    main()
