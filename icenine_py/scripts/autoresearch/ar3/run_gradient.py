"""AR3 step 3: runs P (paint through the key), R and R2 (refinement of P's map), D (BFS with
indexed seeds) on the gradient sample. Clean images, Q-max 8, FusedCost 1, one worker, the voxel
order of ``default_rng(0)``.

  run_gradient.py --name P [--limit-r N]    # paint, then R, then R2 (each with a time cap)
  run_gradient.py --name D                  # BFS, indexed seeds, with step 2's fit log

Writes ``<OUT>/gradient/runs/<name>/``. P: ``recon.mic`` (the painted map), ``records.json``,
``run.json``, ``preflight.json``, ``paint.json.gz`` (entries, stages, evaluations, scores of each
grain, the loop's log), ``seed_records.json.gz``, ``refine_R.json.gz`` and ``refine_R2.json.gz``
(one row per refined voxel: start cost and hit ratio, the entry's best hit ratio, end orientation,
cost and hit ratio, evaluations, and the trace of the fit), ``maps.npz`` (P, R, R2 orientation
maps), ``trace_check.json``. D: as run C/D of step 2 (``fit_log.json.gz``).

The trace of a fit. ``variance_minimizing_optimize`` (the BFS's local fit) has no ``trajectory``
argument, so the trace wraps ``MCOptimizer._zero_temp_with_variance`` (one call per block of the
Monte Carlo fit): after each block, the evaluations so far (the block's own evaluations and the
earlier ones, not the final overlap evaluation), the best cost, whether the block failed (its cost
not below the best so far; the optimizer's own test), and the misorientation (degrees) of the best
orientation at that block to the fit's end. The wrapper draws no random number and changes no
value; ``trace_check.json`` shows the first 200 fits with it on and off from the same generator
state.
"""

import argparse
import copy
import json
import time
from pathlib import Path
from typing import Any, Dict, List, Tuple

import ar3_common as C
import gradient_common as G
import numpy as np
import run_layer as RL

RUNS = G.GRAD / "runs"
R_CAP_S = 60 * 60
R2_CAP_S = 30 * 60
P_CAP_S = 90 * 60
D_CAP_S = 150 * 60


class Tracer:
    """Wraps ``MCOptimizer._zero_temp_with_variance`` to record one row per block."""

    def __init__(self) -> None:
        from icenine.orientation_search import MCOptimizer

        self.cls = MCOptimizer
        self.orig = MCOptimizer._zero_temp_with_variance
        self.rows: List[Tuple[int, float, np.ndarray]] = []
        self.on = False
        tracer = self

        def traced(self_: Any, init: np.ndarray, radius: float, n_steps: int) -> Any:
            out = tracer.orig(self_, init, radius, n_steps)
            if tracer.on:
                tracer.rows.append(
                    (int(self_.cost_fn.eval_count), float(out[1]), np.asarray(out[0], np.float64))
                )
            return out

        MCOptimizer._zero_temp_with_variance = traced

    def reset(self) -> None:
        self.rows = []

    def trace(
        self, R0: np.ndarray, c0: float, R_end: np.ndarray, T: G.Truth, voxel: int
    ) -> Dict[str, Any]:
        """Per block: evals, best cost, failed, degrees of the best orientation to the end (de)
        and to the voxel's truth (te); te0: the start's error."""
        gmin, best = c0, np.asarray(R0, np.float64)
        ev, bc, fl, bests = [], [], [], []
        for e, c, R in self.rows:
            if c >= gmin:
                fl.append(1)
            else:
                fl.append(0)
                gmin, best = c, R
            ev.append(e)
            bc.append(gmin)
            bests.append(best)
        de = (
            T.mis(np.stack(bests), np.broadcast_to(R_end, (len(bests), 3, 3))).tolist()
            if bests
            else []
        )
        te = T.err(np.stack(bests), voxel) if bests else np.zeros(0)
        return dict(
            te=[round(float(x), 5) for x in te],
            te0=round(float(T.err(np.asarray(R0, np.float64)[None], voxel)[0]), 5),
            ev=ev,
            bc=[round(x, 6) for x in bc],
            fl=fl,
            de=[round(x, 5) for x in de],
            end_ok=bool(bests and T.mis(best[None], R_end[None])[0] < 1e-4),
        )


def paint_dict(res: Any) -> Dict[str, Any]:
    return dict(
        entries=[RL.entry_row(e) for e in res.entries],
        stage={str(k): v for k, v in res.stage.items()},
        fits={
            str(v): dict(R=f.R.tolist(), cost=f.cost, hit=f.hit, from_grain=f.from_grain)
            for v, f in res.fits.items()
        },
        failed={str(v): dict(cost=f.cost, hit=f.hit) for v, f in res.failed.items()},
        evals=res.evals,
        wall_s=res.wall_s,
        counts=res.counts,
        log=res.log,
        scores=[[[int(v), s[0], s[1], bool(s[2])] for v, s in sc.items()] for sc in res.scores],
    )


def refine_pass(
    hooks: Any,
    tracer: Tracer,
    T: G.Truth,
    voxels: List[int],
    starts: Dict[int, np.ndarray],
    start_cost: Dict[int, float],
    start_hit: Dict[int, float],
    entry_best: Dict[int, float],
    cap_s: float,
    trace: bool,
) -> Tuple[List[Dict[str, Any]], bool]:
    """One local fit per voxel in ``voxels`` order. Stops at ``cap_s`` (second return: stopped)."""
    rows: List[Dict[str, Any]] = []
    t0 = time.time()
    tracer.on = trace
    for i, v in enumerate(voxels):
        if time.time() - t0 > cap_s:
            tracer.on = False
            return rows, True
        tracer.reset()
        R0 = starts[v]
        Re, cost, hit, n = hooks.local_fit(v, R0)
        row: Dict[str, Any] = dict(
            voxel=int(v),
            start_cost=float(start_cost[v]),
            start_hit=float(start_hit[v]),
            entry_best_hit=float(entry_best[v]),
            end=np.round(Re, 9).tolist(),
            end_cost=float(cost),
            end_hit=float(hit),
            evals=int(n),
        )
        if trace:
            row["trace"] = tracer.trace(R0, float(start_cost[v]), Re, T, v)
        rows.append(row)
        if (i + 1) % 1000 == 0:
            print(f"  fit {i + 1}/{len(voxels)} {time.time() - t0:.0f}s", flush=True)
    tracer.on = False
    return rows, False


def run_p(a: argparse.Namespace) -> None:
    out = RUNS / "P"
    out.mkdir(parents=True, exist_ok=True)
    label = "not a timing run" if a.no_timing_gate else C.preflight_label(out)
    setup, bfs, t_setup = C.build(
        C.CONFIG, out, sample=G.GRID_MIC, extra="PaintGrains 1", stack_name=G.STACK
    )
    mic = setup.sample.get_mic()
    n = len(mic.voxels)
    p = setup.search_params
    assert bool(p.paint_grains) and p.seed_by_indexing and p.fused_cost
    T = G.Truth()
    assert n == T.n
    order = list(range(n))
    np.random.default_rng(0).shuffle(order)

    rng = np.random.default_rng(0)
    t0 = time.time()
    done = bfs.reconstruct_sample(output_mic=str(out / "recon.mic"), rng=rng)
    wall_p = time.time() - t0
    res, hooks = bfs.paint_result, bfs.paint_hooks
    assert res.log[0]["voxel"] == order[0], "voxel order differs from default_rng(0)"
    (out / "records.json").write_text(
        json.dumps({str(k): vars(v) for k, v in bfs.records.items()}) + "\n"
    )
    G.dump_gz(paint_dict(res), out / "paint.json.gz")
    G.dump_gz(hooks.seed_records, out / "seed_records.json.gz")
    info: Dict[str, Any] = dict(
        name="P",
        n_voxels=n,
        n_processed=len(done),
        wall_total_s=wall_p,
        setup_s=t_setup,
        timing_label=label,
        n_entries=len(res.entries),
        evals=res.evals,
        evals_total=int(sum(res.evals.values())),
        counts=res.counts,
        index_pairing_s=float(bfs.stats.get("index_pairing_s", 0.0)),
        cost_fn_eval_count=int(hooks.cost_fn.eval_count),
        stats=bfs.stats,
        truth_sha256=T.sha,
    )
    (out / "run.json").write_text(json.dumps(info, indent=2, default=float) + "\n")
    # the generator's state after the paint: the refinement stage goes on from it
    (out / "rng_state.json").write_text(json.dumps(hooks.rng.bit_generator.state) + "\n")
    print("P done", wall_p, "s;", info["evals_total"], "evaluations", flush=True)


def load_paint(out: Path) -> Any:
    """A ``PaintResult`` rebuilt from ``paint.json.gz`` (entries, scores, stages, fits, failed)."""
    from icenine import grain_paint as GP

    d = G.load_gz(out / "paint.json.gz")
    res = GP.PaintResult()
    for e in d["entries"]:
        res.entries.append(
            GP.GrainEntry(
                voxel=e["voxel"],
                R=np.array(e["R"]),
                cost=e["cost"],
                hit=e["hit"],
                kind=e["kind"],
                best_hit=e["best_hit"],
            )
        )
    res.scores = [{int(v): (c, h, bool(ok)) for v, c, h, ok in sc} for sc in d["scores"]]
    res.stage = {int(k): s for k, s in d["stage"].items()}
    res.fits = {
        int(k): GP.VoxelFit(np.array(f["R"]), f["cost"], f["hit"], f["from_grain"])
        for k, f in d["fits"].items()
    }
    res.failed = {
        int(k): GP.VoxelFit(np.zeros((3, 3)), f["cost"], f["hit"]) for k, f in d["failed"].items()
    }
    return res


def run_refine(a: argparse.Namespace) -> None:
    """R and R2 from the stored paint of P (``--smoke``: no stored paint; the starts are the
    grains' centre orientations of the first voxels, to test this stage)."""
    from icenine import grain_paint as GP

    out = RUNS / ("smoke" if a.smoke else "P")
    out.mkdir(parents=True, exist_ok=True)
    setup, bfs, _t = C.build(
        C.CONFIG, out, sample=G.GRID_MIC, extra="PaintGrains 1", stack_name=G.STACK
    )
    mic = setup.sample.get_mic()
    n = len(mic.voxels)
    T = G.Truth()
    assert n == T.n
    order = list(range(n))
    np.random.default_rng(0).shuffle(order)
    rng = np.random.default_rng(0)
    hooks = GP.ReconstructorHooks(bfs, mic, rng)
    painted: List[int] = []
    start_cost: Dict[int, float] = {}
    start_hit: Dict[int, float] = {}
    entry_best: Dict[int, float] = {}
    starts: Dict[int, np.ndarray] = {}
    P_map = np.zeros((n, 3, 3))
    if a.smoke:
        for v in order[: a.limit_r or 300]:
            starts[v] = T.R_centre[v]
            start_cost[v], start_hit[v] = hooks.score(starts[v], v)
            entry_best[v] = 1.0
            painted.append(v)
            P_map[v] = starts[v]
    else:
        res = load_paint(RUNS / "P")
        rng.bit_generator.state = json.loads((RUNS / "P" / "rng_state.json").read_text())
        P_map, _have = G.read_map(RUNS / "P" / "recon.mic", n)  # seeds, fits, failed voxels
        for v in range(n):
            _Rv, src = res.final_orientation(v)
            if src == "grain" and res.stage.get(v) != "seed":
                g = res.painted_grain(v)
                assert g is not None
                painted.append(v)
                starts[v] = np.asarray(res.entries[g].R, np.float64)
                P_map[v] = starts[v]
                start_cost[v], start_hit[v] = res.scores[g][v][0], res.scores[g][v][1]
                entry_best[v] = res.entries[g].best_hit
    voxels_r = [v for v in order if v in starts]
    print(f"{len(voxels_r)} painted voxels to refine", flush=True)

    tracer = Tracer()
    tracer.on = False
    if a.limit_r and not a.smoke:
        voxels_r = voxels_r[: a.limit_r]
    # ---- the trace check: the first 200 fits with the trace off, then on, from the same state
    state0 = copy.deepcopy(hooks.rng.bit_generator.state)
    k = min(200, len(voxels_r))
    off, _ = refine_pass(
        hooks, tracer, T, voxels_r[:k], starts, start_cost, start_hit, entry_best, 1e9, False
    )
    hooks.rng.bit_generator.state = copy.deepcopy(state0)
    on, _ = refine_pass(
        hooks, tracer, T, voxels_r[:k], starts, start_cost, start_hit, entry_best, 1e9, True
    )
    same = sum(1 for x, y in zip(off, on) if x["end"] == y["end"] and x["evals"] == y["evals"])
    (out / "trace_check.json").write_text(
        json.dumps(
            dict(
                n_fits=k, same_end_and_evals=same, end_ok=sum(1 for y in on if y["trace"]["end_ok"])
            )
        )
        + "\n"
    )
    print("trace check:", same, "of", k, "identical", flush=True)
    hooks.rng.bit_generator.state = copy.deepcopy(state0)
    nfit_before = len(hooks.fit_records)

    # ---- R
    t1 = time.time()
    rows_r, stopped_r = refine_pass(
        hooks, tracer, T, voxels_r, starts, start_cost, start_hit, entry_best, R_CAP_S, True
    )
    wall_r = time.time() - t1
    G.dump_gz(rows_r, out / "refine_R.json.gz")
    R_map = P_map.copy()
    for r in rows_r:
        R_map[r["voxel"]] = np.array(r["end"])
    print(f"R done {wall_r:.0f}s, {len(rows_r)} fits, stopped={stopped_r}", flush=True)

    # ---- R2: every fourth voxel of the voxel order that R refined
    refined = {r["voxel"]: r for r in rows_r}
    v2 = [v for i, v in enumerate(order) if i % 4 == 0 and v in refined]
    st2 = {v: np.array(refined[v]["end"]) for v in v2}
    sc2 = {v: refined[v]["end_cost"] for v in v2}
    sh2 = {v: refined[v]["end_hit"] for v in v2}
    t2 = time.time()
    rows_r2, stopped_r2 = refine_pass(
        hooks, tracer, T, v2, st2, sc2, sh2, entry_best, R2_CAP_S, True
    )
    wall_r2 = time.time() - t2
    G.dump_gz(rows_r2, out / "refine_R2.json.gz")
    R2_map = R_map.copy()
    for r in rows_r2:
        R2_map[r["voxel"]] = np.array(r["end"])
    np.savez_compressed(
        out / "maps.npz",
        P=P_map.astype(np.float32),
        R=R_map.astype(np.float32),
        R2=R2_map.astype(np.float32),
    )
    rj = out / "run.json"
    info: Dict[str, Any] = json.loads(rj.read_text()) if rj.exists() else dict(name="smoke")
    info.update(
        n_painted=len(painted),
        n_refined_R=len(rows_r),
        n_refined_R2=len(rows_r2),
        wall_R_s=wall_r,
        wall_R2_s=wall_r2,
        stopped_R=stopped_r,
        stopped_R2=stopped_r2,
        evals_R=int(sum(r["evals"] for r in rows_r)),
        evals_R2=int(sum(r["evals"] for r in rows_r2)),
        n_fit_records_before_R=nfit_before,
    )
    (out / "run.json").write_text(json.dumps(info, indent=2, default=float) + "\n")
    print("all done", flush=True)


def run_d(a: argparse.Namespace) -> None:
    out = RUNS / "D"
    out.mkdir(parents=True, exist_ok=True)
    label = "not a timing run" if a.no_timing_gate else C.preflight_label(out)
    setup, bfs, t_setup = C.build(C.CONFIG, out, sample=G.GRID_MIC, stack_name=G.STACK)
    mic = setup.sample.get_mic()
    n = len(mic.voxels)
    p = setup.search_params
    assert not p.paint_grains and p.seed_by_indexing and p.fused_cost
    rows = RL.install_fit_log(bfs, mic)
    rng = np.random.default_rng(0)
    t0 = time.time()
    order = bfs.reconstruct_sample(output_mic=str(out / "recon.mic"), rng=rng)
    wall = time.time() - t0
    (out / "records.json").write_text(
        json.dumps({str(k): vars(v) for k, v in bfs.records.items()}) + "\n"
    )
    G.dump_gz(rows, out / "fit_log.json.gz")
    st = bfs.stats
    info = dict(
        name="D",
        n_voxels=n,
        n_processed=len(order),
        wall_total_s=wall,
        setup_s=t_setup,
        timing_label=label,
        stats=st,
        n_fits_logged=len(rows),
        evals_total=int(
            sum(st[k] for k in ("n_evals_seed", "n_evals_neighbor", "n_evals_revisit"))
            + st.get("n_evals_reseed", 0)
            + st.get("n_evals_restart", 0)
        ),
        index_pairing_s=float(st.get("index_pairing_s", 0.0)),
    )
    (out / "run.json").write_text(json.dumps(info, indent=2, default=float) + "\n")
    print("D done", wall, "s", flush=True)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--name", required=True, choices=["P", "D"])
    ap.add_argument("--limit-r", type=int, default=0, help="refine only the first N (pilot)")
    ap.add_argument("--no-timing-gate", action="store_true")
    ap.add_argument("--stage", default="paint", choices=["paint", "refine"])
    ap.add_argument(
        "--smoke", action="store_true", help="refine stage on centre-orientation starts"
    )
    a = ap.parse_args()
    if a.name == "D":
        run_d(a)
    elif a.stage == "paint":
        run_p(a)
    else:
        run_refine(a)


if __name__ == "__main__":
    main()
