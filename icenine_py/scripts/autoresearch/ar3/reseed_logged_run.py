"""AR3 step 8, parts 1.1 to 1.3: the one logged rerun of S5 on ``big`` (the same voxel order).

  uv run python scripts/autoresearch/ar3/reseed_logged_run.py     # cap 90 minutes

This calls step 6's ``run_seed.run("S5", ...)`` unchanged (the same code path as the stored run) and
only adds recorders around instance methods of the run's objects; none of them changes a value the
run uses or consumes the run's random generator:

  ``bfs._fit_from_seed``   the flood counter (one flood per first seed or reseed)
  ``bfs._insert_seed``     who handed its orientation to whom (the expander of each queued voxel)
  ``bfs._local_fit``       one row per local optimization: voxel, flood, expander, phase, start and
                           end orientation, hit ratio, the flood's best hit ratio before the fit,
                           evaluations
  ``bfs._refit_clusters``  the clusters of rejected voxels at each reseed
  ``rec.index_seed``       (below the seed-step wrapper of S5) the seed step's own detail and, for
                           each refined candidate, the hit ratio and the error to the truth
                           (evaluated after the seed step; no evaluation is added to it)

Writes ``reseed/cache/S5_big_logged/`` (``fits.npz``, ``seeds.json.gz``, ``clusters.json.gz``,
``floods.json.gz``, plus the run's own files). Truth is used for analysis only.
"""

import os
import sys
import threading
import time
from pathlib import Path
from typing import Any, Dict, List

os.environ["RESEED_SAMPLE"] = "big"

import numpy as np  # noqa: E402
import reseed_common as RC  # noqa: E402

import ar3_common as C  # noqa: E402
import gradient_common as G  # noqa: E402
import run_seed as RS  # noqa: E402

CAP_MIN = 90
OUTD = RC.CACHE / "S5_big_logged"
state: Dict[str, Any] = {}


class Recorder:
    def __init__(self, bfs: Any, T: Any) -> None:
        self.bfs = bfs
        self.T = T
        self.flood = -1
        self.floods: List[Dict[str, Any]] = []
        self.inh: Dict[int, int] = {}
        self.fit_meta: List[List[float]] = []  # v, flood, expander, phase, hit, best, evals, conf
        self.fit_R: List[np.ndarray] = []  # start (9) and end (9)
        self.clusters: List[Dict[str, Any]] = []
        self.seeds: List[Dict[str, Any]] = []
        self.n_refit_calls = 0
        mic = bfs.setup.sample.get_mic()
        self.mic = mic

    def install(self) -> None:
        b = self.bfs
        o_fs, o_ins, o_lf, o_rc = b._fit_from_seed, b._insert_seed, b._local_fit, b._refit_clusters
        rec = b.reconstructor
        o_is = rec.index_seed
        me = self

        def fit_from_seed(mic: Any, seed_idx: int, rng: Any, force_revisit: bool = False, **k: Any):
            me.flood += 1
            me.floods.append(dict(flood=me.flood, voxel=int(seed_idx), reseed=bool(force_revisit)))
            return o_fs(mic, seed_idx, rng, force_revisit=force_revisit, **k)

        def insert_seed(mic: Any, voxel_idx: int, queue: Any, *a: Any, **k: Any):
            n0 = len(queue)
            out = o_ins(mic, voxel_idx, queue, *a, **k)
            for n_idx in list(queue)[n0:]:
                me.inh[int(n_idx)] = int(voxel_idx)
            return out

        def local_fit(voxel: Any, vertices: Any, rng: Any, rec_: Any, phase: str, *a: Any, **k):
            fr = sys._getframe(1)
            best_before = float(fr.f_locals.get("best_conf", -1.0))
            start = np.asarray(voxel.orientation, np.float64).copy()
            n0 = rec_.n_evals
            res = o_lf(voxel, vertices, rng, rec_, phase, *a, **k)
            r, hit, conf = res
            v = int(rec_.voxel_idx)
            me.fit_meta.append(
                [
                    v,
                    me.flood,
                    me.inh.get(v, -1),
                    {"neighbor": 0, "revisit": 1, "restart": 2}[phase],
                    hit,
                    best_before,
                    rec_.n_evals - n0,
                    conf,
                ]
            )
            me.fit_R.append(
                np.concatenate([start.ravel(), np.asarray(r.orientation, np.float64).ravel()])
            )
            return res

        def refit_clusters(*a: Any, **k: Any):
            todo = o_rc(*a, **k)
            if todo:
                me.clusters.append(
                    dict(call=len(me.clusters), members=[int(i) for i in todo[0]], n_all=len(todo))
                )
            return todo

        def index_seed(vertices: Any, phase_index: int, *a: Any, **k: Any):
            res = o_is(vertices, phase_index, *a, **k)
            v = me.floods[-1]["voxel"]
            hits, errs = [], []
            for R in res.finisher_R:
                info = rec.evaluate_overlap(np.asarray(R, np.float64), vertices, phase_index)
                hits.append(float(b._ratios(info)[0]))
                errs.append(float(me.T.err(np.asarray(R, np.float64)[None], [v])[0]))
            me.seeds.append(
                dict(
                    voxel=v,
                    flood=me.flood,
                    reseed=me.floods[-1]["reseed"],
                    evals=int(res.evals),
                    costed=len(res.strict_costs),
                    strict_costs=[float(c) for c in res.strict_costs],
                    refined=[int(j) for j in res.refined],
                    fin_evals=[int(e) for e in res.finisher_evals],
                    fin_costs=[float(c) for c in res.finisher_costs],
                    fin_hit=hits,
                    fin_err=errs,
                    taken=int(res.taken),
                    n_cand=int(res.n_cand),
                    n_att=int(res.n_att),
                )
            )
            return res

        b._fit_from_seed = fit_from_seed
        b._insert_seed = insert_seed
        b._local_fit = local_fit
        b._refit_clusters = refit_clusters
        rec.index_seed = index_seed

    def dump(self) -> None:
        OUTD.mkdir(parents=True, exist_ok=True)
        np.savez_compressed(
            OUTD / "fits.npz",
            meta=np.asarray(self.fit_meta, np.float64).reshape(-1, 8),
            R=np.asarray(self.fit_R, np.float64).reshape(-1, 18).astype(np.float32),
        )
        RC.write_gz(self.seeds, OUTD / "seeds.json.gz")
        RC.write_gz(self.clusters, OUTD / "clusters.json.gz")
        RC.write_gz(self.floods, OUTD / "floods.json.gz")


def patch_build() -> None:
    orig_build = C.build

    def build(*a: Any, **k: Any) -> Any:
        setup, bfs, t = orig_build(*a, **k)
        state["bfs"] = bfs
        state["rec"] = Recorder(bfs, G.Truth())
        state["rec"].install()
        return setup, bfs, t

    C.build = build
    assert RS.C is C


def watchdog(cap_s: float, t0: float) -> None:
    def loop() -> None:
        time.sleep(cap_s)
        r = state.get("rec")
        if r is not None:
            r.dump()
        (OUTD / "capped.json").write_text(
            f'{{"cap_s": {cap_s}, "elapsed_s": {time.time() - t0}}}\n'
        )
        print("CAP reached; stopped", flush=True)
        os._exit(3)

    threading.Thread(target=loop, daemon=True).start()

    def periodic() -> None:
        while True:
            time.sleep(300.0)
            r = state.get("rec")
            try:
                if r is not None:
                    r.dump()
            except Exception:  # a progress dump must never stop the run
                pass

    threading.Thread(target=periodic, daemon=True).start()


def main() -> None:
    OUTD.mkdir(parents=True, exist_ok=True)
    patch_build()
    t0 = time.time()
    watchdog(CAP_MIN * 60.0, t0)
    RS.run("S5", OUTD, None)
    state["rec"].dump()
    (OUTD / "logged.done").write_text(f"{time.time() - t0:.1f}\n")
    assert Path(OUTD).exists()


if __name__ == "__main__":
    main()
