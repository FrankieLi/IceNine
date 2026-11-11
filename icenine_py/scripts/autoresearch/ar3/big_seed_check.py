"""AR3 step 7: the truth passes the seed test on 200 voxels of ``big`` (drawn with seed 0).

  AR3_SAMPLE=big uv run python scripts/autoresearch/ar3/big_seed_check.py

The seed test: the hit ratio of an orientation is at least ``MinAccelerationThreshold`` (0.80).
Here the orientation is the voxel's true one (analysis only, no decision reads it). The same 200
voxel indices (``default_rng(0).choice(24570, 200, replace=False)``) are read from step 3's stored
Reading 0 (``gradient/reading0.json.gz``, hit ratio at the truth of every gradient-sample voxel),
so nothing is rerun on the gradient sample. Writes ``big/seed_check.json``.
"""

import json

import ar3_common as C
import gradient_common as G
import numpy as np

assert G.BIG, "run with AR3_SAMPLE=big"
GRAD_READING0 = C.ICE / "benchmarks" / "autoresearch" / "ar3" / "gradient" / "reading0.json.gz"


def main() -> None:
    out = G.GRAD / "runs" / "seed_check"
    out.mkdir(parents=True, exist_ok=True)
    setup, bfs, _t = C.build(
        C.CONFIG, out, sample=G.GRID_MIC, extra="PaintGrains 0", stack_name=G.STACK
    )
    from icenine.reconstructor import _get_voxel_vertices

    T = G.Truth()
    mic = setup.sample.get_mic()
    n = len(mic.voxels)
    assert n == T.n == 24570
    thr = float(setup.config.min_acceleration_threshold)
    idx = np.sort(np.random.default_rng(0).choice(n, 200, replace=False))
    cf = bfs.reconstructor._make_local_cost_fn()
    hit_big = []
    for v in idx:
        vox = mic.voxels[int(v)]
        info = cf.evaluate(T.R[int(v)], _get_voxel_vertices(vox), vox.phase)
        hit_big.append(float(bfs._ratios(info)[0]))
    hit_big_a = np.array(hit_big)
    r0 = G.load_gz(GRAD_READING0)
    col = r0["cols"].index("hit_true")
    hit_grad = np.array([r0["rows"][int(v)][col] for v in idx])
    assert all(r0["rows"][int(v)][0] == int(v) for v in idx)
    res = dict(
        threshold=thr,
        n=200,
        voxels_seed=0,
        big_pass=int((hit_big_a >= thr).sum()),
        big_hit_median=float(np.median(hit_big_a)),
        big_hit_min=float(hit_big_a.min()),
        gradient_pass=int((hit_grad >= thr).sum()),
        gradient_hit_median=float(np.median(hit_grad)),
        gradient_hit_min=float(hit_grad.min()),
        gradient_pass_all_voxels=int(
            sum(1 for r in r0["rows"] if r[r0["cols"].index("hit_true")] >= thr)
        ),
        n_gradient_all=len(r0["rows"]),
    )
    (G.GRAD.parent / "seed_check.json").write_text(json.dumps(res, indent=1) + "\n")
    print(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()
