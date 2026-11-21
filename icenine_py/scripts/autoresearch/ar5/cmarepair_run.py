"""AR5 step 20, phase A: the CMA repair (arm ``J5C``) on one stored map after the pass.

  RESEED_SAMPLE=<sample> CMAREPAIR_CACHE=<dir> uv run python \
      scripts/autoresearch/ar5/cmarepair_run.py [--arm J5|J5L] [--limit N]

Reads the stored arm's map before the pass and its pass tries from the cache (``cmarepair_data``,
without opening the truth map: the neighbours come from the sample's grid), builds the sample's
reconstructor as the pass does (``boundary_common.Env``), and applies the frozen
``cmarepair_rules.repair`` with one local CMA run of the library per fit
(``Reconstructor.local_optimization`` with ``cmarepair_rules.fit_kwargs()``) and the fit's own
generator. BFS, seeds and the pass are not run: arm ``J5C`` is the stored ``J5`` run plus this
repair. Writes ``<cache>/<arm>C_<sample>/``: ``repair.json`` (counts, per-fit evaluations, changed
voxels, wall time, label "contended" unless the machine was quiet: this is not a timing result),
``after.npy`` (the map after the repair) and ``fits.npz`` (one row per fit). ``--limit N`` takes
the first N flagged voxels (a smoke test; nothing is scored). The freeze check runs first for every
sample (the stop rule: no repair before ``frozen.json`` is committed).
"""

import argparse
import json
import time
from typing import Any

import cmarepair_common as RC  # noqa: E402  (first: guards, freeze for the held samples, cache)

RC.CF.check_or_exit()  # the repair runs only from the committed freeze, for every sample
import cmarepair_data as RD  # noqa: E402
import cmarepair_rules as CR  # noqa: E402
import numpy as np  # noqa: E402


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--arm", default="J5", choices=("J5", "J5L"))
    ap.add_argument("--limit", type=int, default=0)
    a = ap.parse_args()
    s = RC.SAMPLE
    if not s:
        raise SystemExit("cmarepair_run.py: set RESEED_SAMPLE")
    out = RC.RCACHE / (f"{a.arm}C_{s}" + ("_smoke" if a.limit else ""))
    if (out / "repair.json").exists():
        raise SystemExit(f"refusing: {out / 'repair.json'} exists")
    out.mkdir(parents=True, exist_ok=True)
    d = RD.load(s, a.arm, truth=False)
    flagged = d["flagged"][: a.limit] if a.limit else d["flagged"]
    import boundary_common as B  # noqa: E402  (after cmarepair_common: the sample is set)

    env = B.Env("grad")
    assert env.n == d["n"], (env.n, d["n"])
    for v in flagged:  # the reconstructor's neighbour relation is the one the flag used
        assert sorted(int(j) for j in env.nbr[v]) == d["nbrs"][v], v
    rec = env.hooks.rec
    kw = CR.fit_kwargs()

    def fit(v: int, R0: np.ndarray, rng: np.random.Generator) -> Any:
        r = rec.local_optimization(
            voxel_vertices=env.hooks.vertices(v),
            phase_index=env.mic.voxels[v].phase,
            initial_orientation=np.asarray(R0, dtype=np.float64),
            rng=rng,
            **kw,
        )
        hit, _conf = env.bfs._ratios(r.overlap_info)
        return (
            np.asarray(r.orientation, dtype=np.float64),
            float(r.cost),
            float(hit),
            int(rec.last_local_optimization_evals),
        )

    def score(v: int, R: np.ndarray) -> Any:
        c, h, _f = env.score3(R, v)
        return c, h

    t0 = time.time()
    res = CR.repair(d["R"], d["cost"], d["hit"], flagged, d["cands"], fit, score=score)
    wall = time.time() - t0
    np.save(out / "after.npy", res.R)
    np.savez_compressed(
        out / "fits.npz",
        voxel=np.array([r["voxel"] for r in res.fits], np.int32),
        start=np.array([-1 if r["start"] is None else r["start"] for r in res.fits], np.int32),
        end=np.array([r["end"] for r in res.fits], np.float64).reshape(-1, 3, 3),
        cost=np.array([r["cost"] for r in res.fits]),
        hit=np.array([r["hit"] for r in res.fits]),
        evals=np.array([r["evals"] for r in res.fits], np.int32),
    )
    ev = np.array([r["evals"] for r in res.fits])
    info = dict(
        sample=s,
        arm=a.arm,
        limit=a.limit,
        n_flagged=len(res.flagged),
        n_fits=len(res.fits),
        n_changed=len(res.changed),
        changed=res.changed,
        n_evals=int(res.n_evals),
        n_evals_fits=int(res.n_evals_fits),
        n_evals_score=int(res.n_evals_score),
        median_evals_per_fit=float(np.median(ev)) if len(ev) else 0.0,
        wall_s=wall,
        timing_label="contended (not a timing result)",
        optimizer=CR.REPAIR_OPTIMIZER,
        cma_sigma0_deg=CR.CMA_SIGMA0_DEG,
        cma_max_evals=CR.CMA_MAX_EVALS,
    )
    (out / "repair.json").write_text(json.dumps(info, indent=1) + "\n")
    print("done", json.dumps({k: v for k, v in info.items() if k != "changed"}), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
