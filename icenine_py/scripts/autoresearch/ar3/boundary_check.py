"""AR3 step 5, practical check: one try reproduces a fit of step 4's run C.

  uv run python scripts/autoresearch/ar3/boundary_check.py

The first 200 painted voxels of step 4's run C (arm C: CMA fit, first step 0.05 degrees, 60
evaluations, from the painted start; voxel order and generator state as in ``run_fine.py``) are
fitted again with ``boundary_common.try_fit``. Same end orientation (rounded to 9 digits as stored)
and same evaluation count are required for all 200. Also: a strict re-score of each end with the
same cost function that ``Env.score3`` uses gives the stored end cost (so the keep rules compare
like with like). Writes ``boundary/check_try.json``.
"""

import copy
import json

import boundary_common as B
import numpy as np
import run_fine as RF


def main() -> None:
    out = B.CACHE / "check"
    out.mkdir(parents=True, exist_ok=True)
    hooks, _T, voxels, starts, _sc, _sh, state0 = RF.setup(out)
    bfs = hooks.bfs
    rows = B.G.load_gz(B.FINE_RUNS / "C" / "refine.json.gz")
    assert [r["voxel"] for r in rows] == voxels
    hooks.rng.bit_generator.state = copy.deepcopy(state0)
    same = 0
    rescore_equal = 0
    max_cost_diff = 0.0
    k = 200
    for r in rows[:k]:
        v = r["voxel"]
        t = B.try_fit(hooks, bfs, v, starts[v], hooks.rng)
        if np.round(t["end"], 9).tolist() == r["end"] and t["evals"] == r["evals"]:
            same += 1
        info = hooks.cost_fn.evaluate(t["end"], hooks.vertices(v), hooks.mic.voxels[v].phase)
        d = abs(float(info.cost) - t["cost"])
        max_cost_diff = max(max_cost_diff, d)
        rescore_equal += int(d == 0.0)
    res = dict(
        n_fits=k,
        same_end_and_evals=same,
        rescore_cost_equal=rescore_equal,
        rescore_max_abs_diff=max_cost_diff,
        reference="fine/runs/C/refine.json.gz",
    )
    (B.BND / "check_try.json").write_text(json.dumps(res) + "\n")
    print(res, flush=True)
    if same != k:
        raise SystemExit("one try does not reproduce a fit of run C: stop")


if __name__ == "__main__":
    main()
