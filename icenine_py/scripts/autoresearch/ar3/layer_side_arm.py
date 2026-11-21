"""AR3 step 2, the rule "no fit on a contested voxel": the side arm of step 1 on the whole layer.

  layer_side_arm.py [RUN]      # RUN = C (default)

Reads ``<OUT>/layer/runs/RUN/paint.json``. On every contested voxel (two or more grains pass)
whose passing grains are more than 1 degree apart, and which no fit settled: one local fit from
each of the two lowest-cost grains; the lower final cost is kept. Writes
``<OUT>/layer/runs/RUN/side/side.json`` (per voxel: grains, costs, margin, fit costs and hits,
evaluations, kept fit) and ``side_fits.npz`` (the kept orientation per voxel). Evaluations are
counted apart from the paint run's. The generator is ``default_rng(0)``, as in step 1.
"""

import json
import sys
import time
from pathlib import Path
from typing import Any, Dict, List

import ar3_common as C
import numpy as np
from run_layer import LAYER_OUT, LAYER_SAMPLE


def main() -> None:
    name = sys.argv[1] if len(sys.argv) > 1 else "C"
    run = LAYER_OUT / "runs" / name
    paint = json.loads((run / "paint.json").read_text())
    ents = paint["entries"]
    Re = np.array([e["R"] for e in ents], dtype=np.float64)
    passing: Dict[int, List[Any]] = {}
    for g, rows in enumerate(paint["scores"]):
        for v, c, h, ok in rows:
            if ok:
                passing.setdefault(int(v), []).append((g, c, h))
    fit_settled = {int(v) for v, s in paint["stage"].items() if s == "fit"}

    out = run / "side"
    out.mkdir(exist_ok=True)
    setup, bfs, _ = C.build(Path(C.CONFIG), out, sample=LAYER_SAMPLE, extra="PaintGrains 1")
    mic = setup.sample.get_mic()
    from icenine import grain_paint as GP

    rng = np.random.default_rng(0)
    hooks = GP.ReconstructorHooks(bfs, mic, rng)
    rows: List[Dict[str, Any]] = []
    kept_v: List[int] = []
    kept_R: List[np.ndarray] = []
    evals = 0
    t0 = time.perf_counter()
    for v in sorted(passing):
        r = sorted(passing[v], key=lambda t: (t[1], t[0]))
        if len(r) < 2 or v in fit_settled:
            continue
        dmax = max(
            hooks.misorientation_deg(Re[a[0]], Re[b[0]]) for a in r for b in r if a[0] < b[0]
        )
        if dmax <= 1.0:
            continue  # passing entries of one grain: nothing to decide
        fits = []
        for g, cost, hit in r[:2]:
            R, c, h, ne = hooks.local_fit(v, Re[g])
            fits.append((g, R, c, h, ne))
            evals += ne
        k = int(np.argmin([f[2] for f in fits]))
        kept_v.append(v)
        kept_R.append(fits[k][1])
        rows.append(
            dict(
                voxel=int(v),
                grains=[int(f[0]) for f in fits],
                paint_costs=[float(r[0][1]), float(r[1][1])],
                margin=float(r[1][1] - r[0][1]),
                n_passing=len(r),
                fit_costs=[float(f[2]) for f in fits],
                fit_hits=[float(f[3]) for f in fits],
                fit_evals=[int(f[4]) for f in fits],
                kept=k,
            )
        )
        if len(rows) % 200 == 0:
            print(
                f"  {len(rows)} voxels, {evals} evals, {time.perf_counter() - t0:.0f}s", flush=True
            )
    (out / "side.json").write_text(
        json.dumps(dict(rows=rows, evals=int(evals), wall_s=time.perf_counter() - t0)) + "\n"
    )
    np.savez(out / "side_fits.npz", voxel=np.array(kept_v), R=np.array(kept_R))
    print(json.dumps(dict(n_voxels=len(rows), evals=int(evals))), flush=True)


if __name__ == "__main__":
    main()
