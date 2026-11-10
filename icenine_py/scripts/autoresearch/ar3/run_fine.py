"""AR3 step 4, part 1: local fits with a smaller box, from step 3's stored painted starts.

  run_fine.py --check                 # stock box on the first 200 painted voxels vs stored run R
  run_fine.py --arm B06 [--limit N]   # one arm: B06, B03, B015, L (ladder) or C (CMA)

Every arm starts every painted voxel of the gradient sample from the same stored start and from the
same generator state (``runs/P/rng_state.json``), one run per arm, in the voxel order of
``default_rng(0)``. Arm B5 (the stock fit, box 5 degrees) is run R of step 3, read from its stored
record, and is not rerun here.

The box. ``local_optimization`` passes ``search_box_side`` (radians; the optimizer takes
``tan(search_box_side)``) to ``MCOptimizer.variance_minimizing_optimize``. This script replaces
that argument by ``math.radians(box_deg)``, and only while ``FitBox.active`` is set, which the
fit function sets around its own ``local_optimization`` call. No config key is used.

Writes ``<OUT>/fine/runs/<arm>/refine.json.gz`` (one row per voxel: end orientation, cost, hit
ratio, evaluations, the trace of the fit) and ``run.json``.
"""

import argparse
import copy
import json
import math
import time
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import ar3_common as C
import gradient_common as G
import numpy as np
import run_gradient as RG

FINE = C.OUT / "fine"
FRUNS = FINE / "runs"
STOCK_DEG = 5.0
# box (degrees) of the single-box arms
BOX = {"B06": 5.0 / 8, "B03": 5.0 / 16, "B015": 5.0 / 32}
LADDER = (5.0 / 32, 5.0 / 8, STOCK_DEG)  # the ladder's boxes, finest first
CMA_SIGMA0_DEG = 0.05
CMA_MAX_EVALS = 60


class FitBox:
    """Substitutes ``search_box_side`` of ``variance_minimizing_optimize`` while ``active``."""

    def __init__(self) -> None:
        from icenine.orientation_search import MCOptimizer

        self.active = False
        self.box_rad: Optional[float] = None
        self.n_in = 0  # variance_minimizing_optimize calls with the substitution on
        self.n_out = 0  # calls with it off
        self.orig = MCOptimizer.variance_minimizing_optimize
        me = self

        def patched(self_: Any, *args: Any, **kw: Any) -> Any:
            if me.active and me.box_rad is not None:
                me.n_in += 1
                kw["search_box_side"] = me.box_rad
                if len(args) > 1:  # search_box_side given by position
                    args = (args[0], me.box_rad) + tuple(args[2:])
            else:
                me.n_out += 1
            return me.orig(self_, *args, **kw)

        MCOptimizer.variance_minimizing_optimize = patched

    def set_deg(self, deg: float) -> None:
        self.box_rad = math.radians(deg)


class CostRecorder:
    """Records (cost, orientation) of every ``VoxelCostFunction.evaluate`` call while ``on``. It
    draws no random number and changes no value."""

    def __init__(self) -> None:
        from icenine.cost_functions import VoxelCostFunction

        self.on = False
        self.calls: List[Tuple[float, np.ndarray]] = []
        orig = VoxelCostFunction.evaluate
        me = self

        def patched(self_: Any, orientation: Any, *args: Any, **kw: Any) -> Any:
            out = orig(self_, orientation, *args, **kw)
            if me.on:
                me.calls.append((float(out.cost), np.array(orientation, np.float64)))
            return out

        VoxelCostFunction.evaluate = patched


CMA_POP = 7  # the cma package's default population in 3 dimensions (4 + int(3 ln 3))


class Fitter:
    """One local fit of a voxel from a start, with a box, as the ladder or as CMA."""

    def __init__(self, hooks: Any, tracer: Any, T: Any, fb: FitBox, cr: Any = None) -> None:
        self.h, self.tr, self.T, self.fb, self.cr = hooks, tracer, T, fb, cr

    def one(
        self,
        v: int,
        R0: np.ndarray,
        c0: float,
        box_deg: Optional[float],
        cma: bool = False,
        trace: bool = True,
    ) -> Dict[str, Any]:
        h = self.h
        rec = h.rec
        vox = h.mic.voxels[v]
        self.tr.reset()
        self.fb.active = box_deg is not None
        if box_deg is not None:
            self.fb.set_deg(box_deg)
        kw: Dict[str, Any] = {}
        old = rec.params.local_optimizer
        if cma:
            rec.params.local_optimizer = "cma"
            kw = dict(cma_sigma0_deg=CMA_SIGMA0_DEG, cma_max_evals=CMA_MAX_EVALS)
        use_cr = bool(cma and trace and self.cr is not None)
        if use_cr:
            self.cr.calls = []
            self.cr.on = True
        try:
            r = rec.local_optimization(
                voxel_vertices=h.vertices(v),
                phase_index=vox.phase,
                initial_orientation=np.asarray(R0, np.float64),
                rng=h.rng,
                **kw,
            )
        finally:
            rec.params.local_optimizer = old
            self.fb.active = False
            if use_cr:
                self.cr.on = False
        hit, _ = h.bfs._ratios(r.overlap_info)
        Re = np.asarray(r.orientation, np.float64)
        out: Dict[str, Any] = dict(
            end=Re,
            end_cost=float(r.cost),
            end_hit=float(hit),
            evals=int(rec.last_local_optimization_evals),
        )
        if trace and not cma:
            out["trace"] = self.tr.trace(R0, c0, Re, self.T, v)
        elif use_cr:
            # blocks of CMA_POP evaluations after the start's own (the call after the optimizer,
            # the final overlap evaluation, is dropped): as the MC fit's blocks, (evaluations so
            # far, best cost in the block, its orientation)
            calls = self.cr.calls[:-1]
            assert len(calls) == out["evals"] - 1
            self.tr.reset()
            for j in range(1, len(calls), CMA_POP):
                blk = calls[j : j + CMA_POP]
                k = min(range(len(blk)), key=lambda i: blk[i][0])
                self.tr.rows.append((j + len(blk), blk[k][0], blk[k][1]))
            out["trace"] = self.tr.trace(R0, c0, Re, self.T, v)
        return out

    def row(
        self, v: int, R0: np.ndarray, c0: float, h0: float, arm: str, trace: bool = True
    ) -> Dict[str, Any]:
        """The fit of arm ``arm`` on voxel ``v``. The ladder runs rungs from the same start; a
        rung failed when its end cost is not below the start cost."""
        if arm == "L":
            rungs: List[Dict[str, Any]] = []
            for box in LADDER:
                f = self.one(v, R0, c0, box, trace=trace)
                f["box_deg"] = box
                rungs.append(f)
                if f["end_cost"] < c0:
                    break
            last = rungs[-1]
            best = min(rungs, key=lambda f: f["end_cost"])  # first on ties
            fit = dict(best)
            fit["evals"] = int(sum(f["evals"] for f in rungs))
            fit["rungs"] = [
                dict(
                    box_deg=f["box_deg"],
                    evals=f["evals"],
                    end_cost=f["end_cost"],
                    trace=f.get("trace"),
                )
                for f in rungs
            ]
            fit["last_rung"] = last["box_deg"]
            fit["trace"] = best.get("trace")
        elif arm == "C":
            fit = self.one(v, R0, c0, None, cma=True, trace=trace)
        elif arm == "B5":
            fit = self.one(v, R0, c0, None, trace=trace)
        else:
            fit = self.one(v, R0, c0, BOX[arm], trace=trace)
        row: Dict[str, Any] = dict(
            voxel=int(v),
            start_cost=float(c0),
            start_hit=float(h0),
            end=np.round(fit["end"], 9).tolist(),
            end_cost=fit["end_cost"],
            end_hit=fit["end_hit"],
            evals=fit["evals"],
        )
        for k in ("rungs", "last_rung", "trace"):
            if fit.get(k) is not None:
                row[k] = fit[k]
        return row


def setup(
    out: Path,
) -> Tuple[Any, Any, List[int], Dict[int, Any], Dict[int, float], Dict[int, float], Any]:
    """(hooks, truth, voxels in R's order, starts, start cost, start hit, rng state0)."""
    from icenine import grain_paint as GP

    P = RG.RUNS / "P"
    b, bfs, _t = C.build(
        C.CONFIG, out, sample=G.GRID_MIC, extra="PaintGrains 1", stack_name=G.STACK
    )
    mic = b.sample.get_mic()
    n = len(mic.voxels)
    T = G.Truth()
    assert n == T.n
    order = list(range(n))
    np.random.default_rng(0).shuffle(order)
    rng = np.random.default_rng(0)
    hooks = GP.ReconstructorHooks(bfs, mic, rng)
    res = RG.load_paint(P)
    state0 = json.loads((P / "rng_state.json").read_text())
    hooks.rng.bit_generator.state = copy.deepcopy(state0)
    starts: Dict[int, np.ndarray] = {}
    sc: Dict[int, float] = {}
    sh: Dict[int, float] = {}
    for v in range(n):
        _Rv, src = res.final_orientation(v)
        if src == "grain" and res.stage.get(v) != "seed":
            g = res.painted_grain(v)
            assert g is not None
            starts[v] = np.asarray(res.entries[g].R, np.float64)
            sc[v], sh[v] = res.scores[g][v][0], res.scores[g][v][1]
    voxels = [v for v in order if v in starts]
    return hooks, T, voxels, starts, sc, sh, state0


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--arm", choices=sorted(list(BOX) + ["L", "C"]))
    ap.add_argument("--check", action="store_true")
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--name", default="")
    ap.add_argument("--trace-cma", action="store_true", help="arm C with a per-generation trace")
    ap.add_argument("--second", action="store_true", help="then a second fit, every fourth voxel")
    a = ap.parse_args()
    assert a.check or a.arm
    FINE.mkdir(parents=True, exist_ok=True)
    out = FRUNS / (a.name or (a.arm or "check"))
    out.mkdir(parents=True, exist_ok=True)
    label = "not a timing run" if a.check else C.preflight_label(out)
    hooks, T, voxels, starts, sc, sh, state0 = setup(out)
    print(len(voxels), "painted voxels", flush=True)
    tracer = RG.Tracer()
    tracer.on = True
    fb = FitBox()
    cr = CostRecorder() if a.trace_cma else None
    fit = Fitter(hooks, tracer, T, fb, cr)
    ebh = None
    rR = G.load_gz(RG.RUNS / "P" / "refine_R.json.gz")
    ebh = {r["voxel"]: r["entry_best_hit"] for r in rR}
    assert [r["voxel"] for r in rR] == voxels, "voxel order differs from stored run R"

    if a.check:
        k = 200
        hooks.rng.bit_generator.state = copy.deepcopy(state0)
        same = 0
        for r in rR[:k]:
            v = r["voxel"]
            row = fit.row(v, starts[v], sc[v], sh[v], "B5")
            if row["end"] == r["end"] and row["evals"] == r["evals"]:
                same += 1
        res = dict(
            n_fits=k,
            same_end_and_evals=same,
            box_deg=STOCK_DEG,
            n_vmo_in=fb.n_in,
            n_vmo_out=fb.n_out,
        )
        (out / "check_b5.json").write_text(json.dumps(res) + "\n")
        print("check B5:", res, flush=True)
        if same != k:
            raise SystemExit("baseline does not reproduce: stop")
        return

    todo = voxels[: a.limit] if a.limit else voxels
    # trace check: first 200 fits of this arm with the trace off and on from the same state
    if a.arm != "C":
        k = min(200, len(todo))
        offs, ons = [], []
        for on, store in ((False, offs), (True, ons)):
            hooks.rng.bit_generator.state = copy.deepcopy(state0)
            tracer.on = on
            for v in todo[:k]:
                store.append(fit.row(v, starts[v], sc[v], sh[v], a.arm, trace=on))
        same = sum(
            1 for x, y in zip(offs, ons) if x["end"] == y["end"] and x["evals"] == y["evals"]
        )
        tracer.on = True
        (out / "trace_check.json").write_text(
            json.dumps(dict(arm=a.arm, n_fits=k, same_end_and_evals=same)) + "\n"
        )
        print("trace check:", same, "of", k, flush=True)
    hooks.rng.bit_generator.state = copy.deepcopy(state0)
    fb.n_in = fb.n_out = 0
    rows: List[Dict[str, Any]] = []
    t0 = time.time()
    for i, v in enumerate(todo):
        rows.append(fit.row(v, starts[v], sc[v], sh[v], a.arm))
        rows[-1]["entry_best_hit"] = ebh[v]
        if (i + 1) % 1000 == 0:
            print(f"  fit {i + 1}/{len(todo)} {time.time() - t0:.0f}s", flush=True)
    wall = time.time() - t0
    G.dump_gz(rows, out / "refine.json.gz")
    if a.trace_cma:
        old = {r["voxel"]: r for r in G.load_gz(FRUNS / "C" / "refine.json.gz")}
        same = sum(
            1
            for r in rows
            if r["end"] == old[r["voxel"]]["end"] and r["evals"] == old[r["voxel"]]["evals"]
        )
        (out / "trace_check.json").write_text(
            json.dumps(dict(arm="C", n_fits=len(rows), same_end_and_evals_as_run_C=same)) + "\n"
        )
        print("trace vs untraced C run:", same, "of", len(rows), flush=True)
    rows2: List[Dict[str, Any]] = []
    if a.second:
        order_all = list(range(T.n))
        np.random.default_rng(0).shuffle(order_all)
        done = {r["voxel"]: r for r in rows}
        v2 = [v for i, v in enumerate(order_all) if i % 4 == 0 and v in done]
        t2 = time.time()
        for v in v2:
            d = done[v]
            r2 = fit.row(v, np.array(d["end"]), d["end_cost"], d["end_hit"], a.arm)
            r2["entry_best_hit"] = ebh[v]
            rows2.append(r2)
        G.dump_gz(rows2, out / "refine2.json.gz")
        print(f"second fits {len(rows2)} {time.time() - t2:.0f}s", flush=True)
    info = dict(
        arm=a.arm,
        n_fits2=len(rows2),
        evals2=int(sum(r["evals"] for r in rows2)),
        n_fits=len(rows),
        evals=int(sum(r["evals"] for r in rows)),
        wall_s=wall,
        timing_label=label,
        box_deg=BOX.get(a.arm),
        ladder_deg=list(LADDER) if a.arm == "L" else None,
        cma=dict(sigma0_deg=CMA_SIGMA0_DEG, max_evals=CMA_MAX_EVALS) if a.arm == "C" else None,
        n_vmo_with_box=fb.n_in,
        n_vmo_without_box=fb.n_out,
        limit=a.limit,
    )
    (out / "run.json").write_text(json.dumps(info, indent=2) + "\n")
    print("done", a.arm, f"{wall:.0f}s", info["evals"], "evaluations", flush=True)


if __name__ == "__main__":
    main()
