#!/usr/bin/env python3
"""U1 / U2 (a start orientation exists): single-worker, interleaved wall time per case.

Cases: the first 5 sweep voxels x 10 radii x directions x 2 variants (clean, realistic), the cases
of scripts/nn_hybrid (per-voxel images, <= 3 distractor sources; the finisher always sees
experiment B's image and B's MC seed). Pipelines per case, in a rotated order:
  H0   FindOptimal alone from the perturbed nominal
  H1   net x1 -> refine_from_candidates([estimate])
  H3   net x3 -> refine_from_candidates
  MCr  MC told the true radius r (optimizer_sweep.run_mc_adam, box 1.5 r), started at the nominal
  HG   Huber GN x3 -> net x3 -> refine_from_candidates   (radii 1.5, 2, 3, 5 only)
The network stage of a pipeline is run on ONE case (batch 1): its own pass-1 prepare_nominal
(production cost), the pass-1 windows taken from the case's row of the sweep batch (harness: they
stand for the experimental windows), passes 2-3 re-prepared and re-rendered. Rendering is harness
work in every pass; "production" time excludes it, "measured" time includes it (and the pass-1
render share). Everything else (ROI / window spec, window decoding, the forward pass, GN, the
finisher) is charged.

  uv run python scripts/profiling/prof_seeded.py run --workers 1 --repeats 3
  uv run python scripts/profiling/prof_seeded.py run --workers 10 --tag w10 --dirs 2 --repeats 1
"""

import argparse
import json
import sys
import time
from pathlib import Path
from typing import Any, Dict, List, Tuple

import numpy as np
from scipy.spatial.transform import Rotation

sys.path.insert(0, str(Path(__file__).resolve().parent))
import prof_common as PC  # noqa: E402

import findoptimal_sweep as fs  # noqa: E402
import optimizer_sweep as osw  # noqa: E402
import perturbation_sweep as ps  # noqa: E402
import run as NH  # noqa: E402  (scripts/nn_hybrid/run.py)

N_VOX, N_DIRS_MAIN = 5, 4
MODEL = "realistic_s0"
PIPES = ["H0", "H1", "H3", "MCr", "HG"]
HG_RI = set(NH.HG_RADII_IDX)
CACHE = PC.CACHE / "seeded"
_S: Dict[str, Any] = {}


def init_worker(wargs: Dict[str, Any]) -> None:
    fs.init_worker(wargs)
    W = fs._W
    assert W is not None
    paths = dict((n, p) for n, p in W.ctx.args.models)
    net = ps.load_model(paths[MODEL])[0]
    T = PC.ST.StageTimer()
    I = PC.Instrument(T, W.rec.params.max_mc_steps, net=True)  # noqa: E741
    I.install()
    I.mark(W.local_fn, "local")
    I.mark(W.ctx.hard_fn, "mc")
    _S.update(T=T, I=I, net=net, warm=False)


def pipes_for(ri: int) -> List[str]:
    return [p for p in PIPES if p != "HG" or ri in HG_RI]


def _rows(b: Dict[str, Any], j: int) -> Dict[str, Any]:
    return {k: v[j : j + 1] for k, v in b.items()}


def run_pipeline(pipe: str, cs: Dict[str, Any]) -> Dict[str, Any]:
    """One timed pipeline on one case. cs: the case (see run_task)."""
    W = fs._W
    assert W is not None
    T, I, net = _S["T"], _S["I"], _S["net"]  # noqa: E741
    ctx, a = W.ctx, W.ctx.args
    vctx, j, vb = cs["vctx"], cs["j"], cs["vb"]
    R_true, R_nom = vctx.R_true, vb.R_nom0[j]
    lf = W.local_fn
    seed_b, seed_j = cs["seed_b"], vb.seed + 7919 * j
    snap = T.snapshot()
    la = __import__("os").getloadavg()[0]
    n0 = lf.eval_count
    fallback, conv, n_evals_mc = False, False, 0
    t_fin_override = None
    t0 = time.perf_counter()
    if pipe == "H0":
        with T.stage("finisher"):
            R_f, cost_f, conv = NH.refine_fo(R_nom, vctx, seed_b, None)
    elif pipe == "MCr":
        with T.stage("finisher"):
            o = osw.run_mc_adam(
                ctx, vctx, cs["keys"], R_nom, float(a.radii[cs["ri"]]), ["mc"], seed_b
            )["mc"]
        R_f = np.asarray(o["R_final"], dtype=np.float64)
        n_evals_mc = int(o["n_evals"])
        t_fin_override = float(o["seconds"])
        cost_f = float("nan")
    else:
        ctx_p = ctx if pipe != "H1" else PC.clone_ctx_with_passes(ctx, 1)
        dj = [cs["draws"][j]]
        with T.stage("net"):
            p1, _why = ps.prepare_nominal(ctx, cs["vidx"], R_nom)
            ok1 = np.array([p1 is not None and bool(vb.ok1[j])])
            if pipe == "HG":
                with T.stage("gn"):
                    g = osw.run_gn(
                        ctx,
                        cs["vidx"],
                        [p1],
                        _rows(vb.b1, j),
                        ok1,
                        vb.delta0[j : j + 1],
                        R_nom[None],
                        R_true,
                        dj,
                        vctx.sources,
                        vb.variant,
                        seed_j,
                        1.0,
                    )
                e = np.stack([g["err_x"][:, -1], g["err_y"][:, -1], g["err_z"][:, -1]], axis=-1)
                good = np.isfinite(e).all(axis=-1)
                R_gn = R_nom[None].copy()
                if good.any():
                    R_gn[good] = (
                        Rotation.from_rotvec(np.radians(e[good].astype(np.float64))).as_matrix()
                        @ R_true
                    )
                dstart = ps.relative_offset_deg(R_true, R_gn)
                p2, _ = ps.prepare_nominal(ctx, cs["vidx"], R_gn[0]) if ok1[0] else (None, 0)
                with T.stage("render"):
                    b2 = ps.render_batch([p2], dstart, dj, vctx.sources, vb.variant, seed_j, a)
                ok2 = np.array(
                    [p2 is not None and bool(b2["n_present"].numpy()[0] >= a.min_present)]
                )
                Rn, _Ln = NH.net_stage(
                    ctx_p,
                    net,
                    cs["vidx"],
                    [p2],
                    b2,
                    ok2,
                    dstart,
                    R_gn,
                    R_true,
                    dj,
                    vctx.sources,
                    vb.variant,
                    seed_j,
                    timer=T,
                )
                for p_i in range(Rn.shape[1]):
                    miss = np.isnan(Rn[:, p_i, 0, 0]) & ok1
                    Rn[miss, p_i] = R_gn[miss]
            else:
                Rn, _Ln = NH.net_stage(
                    ctx_p,
                    net,
                    cs["vidx"],
                    [p1],
                    _rows(vb.b1, j),
                    ok1,
                    vb.delta0[j : j + 1],
                    R_nom[None],
                    R_true,
                    dj,
                    vctx.sources,
                    vb.variant,
                    seed_j,
                    timer=T,
                )
        est = Rn[0, -1]
        fallback = not (bool(ok1[0]) and np.isfinite(est).all())
        start = R_nom if fallback else est
        cs.setdefault("net_err", {})[pipe] = (
            float("nan") if fallback else PC.reduced_err_deg(est, R_true)
        )
        with T.stage("finisher"):
            R_f, cost_f, conv = NH.refine_fo(start, vctx, seed_b, None)
    wall = time.perf_counter() - t0
    sd = PC.stage_delta(T, snap)
    inc = sd["inclusive"]
    render = inc.get("render", 0.0)
    fin = t_fin_override if t_fin_override is not None else inc.get("finisher", 0.0)
    evals_fo = lf.eval_count - n0 if pipe != "MCr" else 0
    n_have = cs["n_batch"]
    total_meas = (wall if t_fin_override is None else fin) + (
        cs["render1_s"] / n_have if pipe in ("H1", "H3", "HG") else 0.0
    )
    total_prod = (wall - render) if t_fin_override is None else fin
    net_total = inc.get("net", 0.0)
    return dict(
        pipe=pipe,
        wall=wall,
        total_measured=total_meas,
        total_production=total_prod,
        net_measured=net_total + (cs["render1_s"] / n_have if net_total else 0.0),
        net_production=net_total - render,
        finisher=fin,
        evals=int(evals_fo if pipe != "MCr" else n_evals_mc),
        fallback=bool(fallback),
        converged=bool(conv),
        R_final=np.asarray(R_f, dtype=np.float64).tolist(),
        err=PC.reduced_err_deg(R_f, R_true),
        cost_final=float(cost_f),
        net_err=cs.get("net_err", {}).get(pipe, float("nan")),
        loadavg_start=la,
        stages=sd,
        evals_all=I.eval_counts(snap),
    )


def run_task(item: Tuple[Any, ...], record: bool = True) -> List[Dict[str, Any]]:
    vidx, vpos, ri, rep, dirs, base, ref_nroi, ref_fail = item
    W = fs._W
    assert W is not None
    ctx = W.ctx
    out: List[Dict[str, Any]] = []
    counter = base
    for vb in NH.variant_batches(ctx, vidx, vpos, ri, ref_nroi, ref_fail, NH.VARIANTS):
        vctx = vb.vctx
        n_batch = max(int(vb.have.sum()), 1)
        for j in dirs:
            if not vb.have[j]:
                continue
            keys = NH.case_keys(ctx, vb, j)
            fs.attach_images(keys)
            cs = dict(
                vctx=vctx,
                vb=vb,
                j=j,
                ri=ri,
                vidx=vidx,
                keys=keys,
                draws=vb.draws,
                seed_b=NH.b_seed(ctx.args, vpos, ri, j, vb.vi),
                render1_s=vb.render1_s,
                n_batch=n_batch,
            )
            for pipe in PC.rotated(pipes_for(ri), counter):
                r = run_pipeline(pipe, cs)
                r.update(
                    vidx=vidx,
                    vpos=vpos,
                    ri=ri,
                    j=int(j),
                    variant=vb.variant,
                    rep=rep,
                    radius=float(ctx.args.radii[ri]),
                    order_idx=counter,
                )
                if record:
                    out.append(r)
            counter += 1
    return out


def task(item: Tuple[Any, ...]) -> str:
    t0 = time.time()
    *core, path, warm = item
    if not _S["warm"]:
        if warm == "full":
            run_task(tuple(core), record=False)  # untimed warm-up of every pipeline
        _S["warm"] = True
    iso = PC.isolation_record(f"seeded v{item[0]} r#{item[2]} rep{item[3]}")
    recs = run_task(tuple(core))
    Path(path).write_text(json.dumps(dict(records=recs, isolation=iso), default=float))
    return f"seeded v{item[0]} r#{item[2]} rep{item[3]} {time.time() - t0:.0f}s"


def work_items(a: argparse.Namespace) -> List[Tuple[Any, ...]]:
    outdir = CACHE / a.tag
    outdir.mkdir(parents=True, exist_ok=True)
    _cfg, sweep = NH.worker_args([MODEL])
    voxels = [int(v) for v in sweep["voxel_indices"]][:N_VOX]
    its: List[Tuple[Any, ...]] = []
    for rep in range(a.repeats):
        for vpos, v in enumerate(voxels):
            for ri in range(10):
                repeat_task = (vpos + ri) % 5 == 0  # 10 of the 50 (voxel, radius) tasks
                if rep > 0 and not repeat_task:
                    continue
                dirs = list(range(a.dirs if rep == 0 else 2))
                path = outdir / f"v{v}_r{ri}_rep{rep}.json"
                if not path.exists():
                    its.append(
                        (
                            v,
                            vpos,
                            ri,
                            rep,
                            dirs,
                            100 * (vpos * 10 + ri) + 13 * rep,
                            sweep["n_roi"][vpos, ri],
                            sweep["fail_pass1"][vpos, ri],
                            str(path),
                            a.warm,
                        )
                    )
    if a.only:  # re-timing of selected tasks (file stems, e.g. v16905_clean_rep0); default off
        its = [it for it in its if Path(it[-2]).stem in a.only]
        missing = sorted(set(a.only) - {Path(it[-2]).stem for it in its})
        if missing:  # unknown stem, or its output file still exists (move it away first)
            print("warning: --only stems with no pending task:", missing, flush=True)
    return its[: a.limit] if a.limit else its


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    ap.add_argument("cmd", choices=["run"])
    ap.add_argument("--workers", type=int, default=1)
    ap.add_argument("--repeats", type=int, default=3)
    ap.add_argument("--dirs", type=int, default=N_DIRS_MAIN)
    ap.add_argument("--tag", default="w1")
    ap.add_argument("--warm", choices=["full", "light"], default="full")
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument(
        "--only",
        nargs="+",
        default=[],
        help="run only these task file stems (move their existing output files away first)",
    )
    a = ap.parse_args()
    its = work_items(a)
    print(len(its), "tasks", flush=True)
    iso = PC.isolation_record(f"seeded start tag={a.tag} workers={a.workers}")
    PC.write_json(CACHE / a.tag / "isolation_start.json", iso)
    print(json.dumps(iso), flush=True)
    wargs, _ = NH.worker_args([MODEL])
    t0 = time.time()
    PC.run_pool(task, init_worker, wargs, its, a.workers, "seeded")
    wall = time.time() - t0
    PC.write_json(CACHE / a.tag / "isolation_end.json", PC.isolation_record("seeded end"))
    PC.write_json(CACHE / a.tag / "wall.json", dict(wall_s=wall, tasks=len(its), workers=a.workers))
    print(f"seeded done, wall {wall:.0f}s", flush=True)


if __name__ == "__main__":
    main()
