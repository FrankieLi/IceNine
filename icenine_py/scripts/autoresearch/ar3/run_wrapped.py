"""AR3 step 4, part 3: whole-layer runs with the picked fit ("fine first, stock as fallback").

  run_wrapped.py --kind D4    # BFS with indexed seeds, gradient sample
  run_wrapped.py --kind P4    # paint through the key, gradient sample, then R4
  run_wrapped.py --kind E4    # BFS with indexed seeds, step 2's sample (one orientation per grain)

No library change. ``Wrap`` replaces the instance method ``local_optimization`` of the run's
reconstructor object, and nothing else. A call of the wrapper runs the picked fit first: the CMA
local fit (``local_optimizer`` "cma" switched on only for that call, first step 0.05 degrees, at
most 60 evaluations, the start's own evaluation included; Part 1's arm C). When the end cost is not
below the start cost, the stock fit (variance-minimizing MC, box 5 degrees, the unchanged code
path) runs from the same start and the end with the lower cost is returned (the picked fit's on a
tie). The evaluations of both fits are counted (``last_local_optimization_evals``).

The caller's acceptance test (BFS: ``hit_ratio / best_conf > 0.9`` in ``_expand``, where
``best_conf`` already includes the fit's own hit ratio; paint: ``hit / max(entry best hit, hit) >
0.9`` in ``grain_paint``) was not evaluated in the wrapper: it would need the caller's local
variables (``best_conf`` of ``_expand``; the entry ``g`` of the paint loop). So the trigger of the
stock fit is "the cost did not fall".

Counts that show the wrapper is the only change: calls of the wrapper by caller; calls of
``CMAOptimizer.optimize`` (all of them are the wrapper's, the run's own local_optimizer being
"classic"); calls of ``MCOptimizer.variance_minimizing_optimize`` inside and outside the wrapper
(outside: the seed path, the grid search and the finisher, unchanged).

Writes ``<OUT>/fine/runs/<kind>/``: ``recon.mic.gz``, ``records.json.gz``, ``run.json``,
``progress.json`` (written every 5 minutes, so a killed run leaves a record), ``preflight.json``,
and for P4 also ``paint.json.gz``, ``refine_R4.json.gz``.
"""

import argparse
import gzip
import json
import shutil
import sys
import threading
import time
from collections import Counter
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import ar3_common as C
import gradient_common as G
import numpy as np
import run_fine as RF
import run_gradient as RG
import run_layer as RL

FRUNS = RF.FRUNS


class Wrap:
    def __init__(self, bfs: Any) -> None:
        from icenine.orientation_search import CMAOptimizer, MCOptimizer

        self.rec = bfs.reconstructor
        self.orig = self.rec.local_optimization
        self.cr = RF.CostRecorder()
        self.calls: Counter = Counter()
        self.n_picked_fell = 0
        self.n_fallback = 0
        self.fallback_won = 0
        self.evals_picked = 0
        self.evals_fallback = 0
        self.n_cma_total = 0
        self.n_vmo_total = 0
        self.n_vmo_in_wrapper = 0
        self.in_wrapper = False
        me = self
        cma_orig = CMAOptimizer.optimize
        vmo_orig = MCOptimizer.variance_minimizing_optimize

        def cma_counted(self_: Any, *a: Any, **k: Any) -> Any:
            me.n_cma_total += 1
            return cma_orig(self_, *a, **k)

        def vmo_counted(self_: Any, *a: Any, **k: Any) -> Any:
            me.n_vmo_total += 1
            if me.in_wrapper:
                me.n_vmo_in_wrapper += 1
            return vmo_orig(self_, *a, **k)

        CMAOptimizer.optimize = cma_counted
        MCOptimizer.variance_minimizing_optimize = vmo_counted
        self.rec.local_optimization = self.wrapped

    def cma_fit(
        self, voxel_vertices: Any, phase_index: int, R0: np.ndarray, rng: Any
    ) -> Tuple[Any, int, float]:
        """The picked fit. Returns (result, evaluations, start cost)."""
        params = self.rec.params
        old = params.local_optimizer
        params.local_optimizer = "cma"
        self.cr.calls = []
        self.cr.on = True
        try:
            r = self.orig(
                voxel_vertices=voxel_vertices,
                phase_index=phase_index,
                initial_orientation=R0,
                rng=rng,
                cma_sigma0_deg=RF.CMA_SIGMA0_DEG,
                cma_max_evals=RF.CMA_MAX_EVALS,
            )
        finally:
            params.local_optimizer = old
            self.cr.on = False
        n = int(self.rec.last_local_optimization_evals)
        return r, n, float(self.cr.calls[0][0])

    def wrapped(
        self,
        voxel_vertices: Any,
        phase_index: int,
        initial_orientation: np.ndarray,
        rng: Optional[Any] = None,
        **kw: Any,
    ) -> Any:
        assert not kw, "the wrapper expects the classic-mode call (no cma keywords)"
        self.calls[sys._getframe(1).f_code.co_name] += 1
        R0 = np.asarray(initial_orientation, np.float64)
        r1, n1, c0 = self.cma_fit(voxel_vertices, phase_index, R0, rng)
        self.evals_picked += n1
        if r1.cost < c0:
            self.n_picked_fell += 1
            self.rec.last_local_optimization_evals = n1
            return r1
        self.n_fallback += 1
        self.in_wrapper = True
        try:
            r2 = self.orig(
                voxel_vertices=voxel_vertices,
                phase_index=phase_index,
                initial_orientation=R0,
                rng=rng,
            )
        finally:
            self.in_wrapper = False
        n2 = int(self.rec.last_local_optimization_evals)
        self.evals_fallback += n2
        self.rec.last_local_optimization_evals = n1 + n2
        if r2.cost < r1.cost:
            self.fallback_won += 1
            return r2
        return r1

    def counts(self) -> Dict[str, Any]:
        return dict(
            wrapper_calls=int(sum(self.calls.values())),
            wrapper_calls_by_caller=dict(self.calls),
            picked_fell=self.n_picked_fell,
            fallback_runs=self.n_fallback,
            fallback_won=self.fallback_won,
            evals_picked=self.evals_picked,
            evals_fallback=self.evals_fallback,
            cma_calls_total=self.n_cma_total,
            vmo_calls_total=self.n_vmo_total,
            vmo_calls_in_wrapper=self.n_vmo_in_wrapper,
            vmo_calls_outside_wrapper=self.n_vmo_total - self.n_vmo_in_wrapper,
        )


def progress_thread(out: Path, w: Wrap, bfs: Any, t0: float, every: float = 300.0) -> None:
    def loop() -> None:
        while True:
            time.sleep(every)
            try:
                st = {k: v for k, v in bfs.stats.items() if isinstance(v, (int, float))}
                (out / "progress.json").write_text(
                    json.dumps(
                        dict(elapsed_s=time.time() - t0, wrapper=w.counts(), stats=st),
                        default=float,
                    )
                    + "\n"
                )
            except Exception:  # a progress write must never stop the run
                pass

    threading.Thread(target=loop, daemon=True).start()


def gz_file(path: Path) -> None:
    with open(path, "rb") as fi, open(str(path) + ".gz", "wb") as raw:
        with gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as fo:
            shutil.copyfileobj(fi, fo)
    path.unlink()


def run_bfs(kind: str, out: Path, max_voxels: Optional[int] = None) -> None:
    label = C.preflight_label(out)
    if kind == "D4":
        setup, bfs, t_setup = C.build(C.CONFIG, out, sample=G.GRID_MIC, stack_name=G.STACK)
    else:
        setup, bfs, t_setup = C.build(C.CONFIG, out, sample=RL.LAYER_SAMPLE)
    mic = setup.sample.get_mic()
    p = setup.search_params
    assert not p.paint_grains and p.seed_by_indexing and p.fused_cost
    assert p.local_optimizer != "cma"
    w = Wrap(bfs)
    rng = np.random.default_rng(0)
    t0 = time.time()
    progress_thread(out, w, bfs, t0)
    order = bfs.reconstruct_sample(
        output_mic=str(out / "recon.mic"), max_voxels=max_voxels, rng=rng
    )
    wall = time.time() - t0
    (out / "records.json").write_text(
        json.dumps({str(k): vars(v) for k, v in bfs.records.items()}) + "\n"
    )
    gz_file(out / "records.json")
    gz_file(out / "recon.mic")
    st = bfs.stats
    info = dict(
        name=kind,
        n_voxels=len(mic.voxels),
        n_processed=len(order),
        wall_total_s=wall,
        setup_s=t_setup,
        timing_label=label,
        stats=st,
        wrapper=w.counts(),
        evals_total=int(
            sum(st[k] for k in ("n_evals_seed", "n_evals_neighbor", "n_evals_revisit"))
            + st.get("n_evals_reseed", 0)
            + st.get("n_evals_restart", 0)
        ),
        index_pairing_s=float(st.get("index_pairing_s", 0.0)),
    )
    (out / "run.json").write_text(json.dumps(info, indent=2, default=float) + "\n")
    print(kind, "done", wall, "s", flush=True)


def run_p4(out: Path) -> None:
    from icenine import grain_paint as GP

    label = C.preflight_label(out)
    setup, bfs, t_setup = C.build(
        C.CONFIG, out, sample=G.GRID_MIC, extra="PaintGrains 1", stack_name=G.STACK
    )
    mic = setup.sample.get_mic()
    n = len(mic.voxels)
    p = setup.search_params
    assert bool(p.paint_grains) and p.seed_by_indexing and p.fused_cost
    assert p.local_optimizer != "cma"
    T = G.Truth()
    assert n == T.n
    w = Wrap(bfs)
    rng = np.random.default_rng(0)
    t0 = time.time()
    progress_thread(out, w, bfs, t0)
    done = bfs.reconstruct_sample(output_mic=str(out / "recon.mic"), rng=rng)
    wall_p = time.time() - t0
    res, hooks = bfs.paint_result, bfs.paint_hooks
    counts_after_paint = w.counts()
    G.dump_gz(RG.paint_dict(res), out / "paint.json.gz")
    (out / "records.json").write_text(
        json.dumps({str(k): vars(v) for k, v in bfs.records.items()}) + "\n"
    )
    gz_file(out / "records.json")
    # ---- R4: one picked fit (plain CMA, no fallback) from each painted start
    order = list(range(n))
    np.random.default_rng(0).shuffle(order)
    starts: Dict[int, np.ndarray] = {}
    start_cost: Dict[int, float] = {}
    for v in range(n):
        _Rv, src = res.final_orientation(v)
        if src == "grain" and res.stage.get(v) != "seed":
            g = res.painted_grain(v)
            assert g is not None
            starts[v] = np.asarray(res.entries[g].R, np.float64)
            start_cost[v] = res.scores[g][v][0]
    voxels = [v for v in order if v in starts]
    final_map = np.zeros((n, 3, 3))
    for v in range(n):
        Rv, _src = res.final_orientation(v)
        final_map[v] = Rv if Rv is not None else np.eye(3)
    rows: List[Dict[str, Any]] = []
    t1 = time.time()
    for v in voxels:
        vox = mic.voxels[v]
        r, nev, c0 = w.cma_fit(hooks.vertices(v), vox.phase, starts[v], hooks.rng)
        hit, _ = bfs._ratios(r.overlap_info)
        rows.append(
            dict(
                voxel=int(v),
                start_cost=start_cost[v],
                end=np.round(r.orientation, 9).tolist(),
                end_cost=float(r.cost),
                end_hit=float(hit),
                evals=int(nev),
            )
        )
        final_map[v] = r.orientation
    wall_r = time.time() - t1
    G.dump_gz(rows, out / "refine_R4.json.gz")
    np.savez_compressed(out / "maps.npz", P4=final_map.astype(np.float32))
    info = dict(
        name="P4",
        n_voxels=n,
        n_processed=len(done),
        wall_p_s=wall_p,
        wall_r4_s=wall_r,
        setup_s=t_setup,
        timing_label=label,
        n_entries=len(res.entries),
        evals=res.evals,
        evals_total=int(sum(res.evals.values())),
        evals_r4=int(sum(r["evals"] for r in rows)),
        n_refined_r4=len(rows),
        counts=res.counts,
        wrapper_after_paint=counts_after_paint,
        wrapper_final=w.counts(),
        index_pairing_s=float(bfs.stats.get("index_pairing_s", 0.0)),
        truth_sha256=T.sha,
    )
    (out / "run.json").write_text(json.dumps(info, indent=2, default=float) + "\n")
    print("P4 done", wall_p, wall_r, flush=True)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--kind", required=True, choices=["D4", "P4", "E4"])
    ap.add_argument("--smoke", type=int, default=0, help="pilot: first N voxels, own folder")
    a = ap.parse_args()
    out = FRUNS / (f"smoke_{a.kind}" if a.smoke else a.kind)
    out.mkdir(parents=True, exist_ok=True)
    if a.kind == "P4":
        run_p4(out)
    else:
        run_bfs(a.kind, out, a.smoke or None)


if __name__ == "__main__":
    main()
