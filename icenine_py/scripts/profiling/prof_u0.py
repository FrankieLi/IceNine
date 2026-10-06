#!/usr/bin/env python3
"""U0 (no start): single-worker, interleaved wall-time of the full reconstruct_voxel pipelines.

Pipelines (all on the E0 images and rng of the earlier tasks, seed 0, so each reproduces its stored
run bit for bit, which is checked per run):
  baseline   reconstruct_voxel                                     (E0)
  F1         baseline answer + the CSL-relative check of f1_run.py (time = baseline run + F1)
  F1b        reconstruct_voxel with the CSL relatives added before each quick MC
  e2         E2 rerank (62-feature GBT as the rank_key)
  proxy      the low-Q Q5 + free Q8 cost regression as the rank_key (row (i) of the proxy study)
  proxy+F1   proxy answer + F1 (time = proxy run + F1)
The four reconstruct runs of one (voxel, variant) are run in a rotated order (case index), each F1
right after its parent run. Instrumentation: scripts/nn_hybrid/stage_timer.py wrappers installed
from prof_common.Instrument (nothing in icenine/ changes).

  uv run python scripts/profiling/prof_u0.py run --workers 1 --repeats 3
  uv run python scripts/profiling/prof_u0.py run --workers 10 --tag w10 --limit 20 --warm light
"""

import argparse
import json
import sys
import time
from pathlib import Path
from typing import Any, Dict, List, Tuple

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import prof_common as PC  # noqa: E402

import base as B  # noqa: E402  (coarse_proxy: also brings findoptimal_robustness common as C)

C, M = B.C, B.M
import csl  # noqa: E402
import endtoend as E  # noqa: E402
import f1_run as F1  # noqa: E402
import features as F  # noqa: E402
import fixes_run as FIX  # noqa: E402
import models as MD  # noqa: E402

RECON = ["baseline", "F1b", "e2", "proxy"]
PIPES = ["baseline", "F1", "F1b", "e2", "proxy", "proxy+F1"]
N_VOX = 20
PROXY_SET, PROXY_TARGET, PROXY_Q = "lowq5+c8", "reg", 5
CACHE = PC.CACHE / "u0"
_S: Dict[str, Any] = {}


def init_worker(wargs: Dict[str, Any]) -> None:
    C.init_worker(wargs)
    W = C.get_worker()
    T = PC.ST.StageTimer()
    I = PC.Instrument(T, W.rec.params.max_mc_steps, net=False)  # noqa: E741
    I.install()
    I.mark(W.local_fn, "local")
    _S.update(T=T, I=I, fe={}, warm=False)


def _fe(q: Any, phase: int) -> Any:
    W = C.get_worker()
    key = (q, phase)
    if key not in _S["fe"]:
        if q is None:
            _S["fe"][key] = F.FeatureExtractor(W.local_fn, W.ctx.geo, phase)
        else:
            _S["fe"][key] = F.FeatureExtractor(W.local_fn, W.ctx.geo, phase, q_max=float(q))
    return _S["fe"][key]


def _models(fold: int) -> Tuple[Any, Any]:
    import joblib

    key = ("m", fold)
    if key not in _S:
        e2 = joblib.load(C.CACHE_DIR / "models" / f"fold{fold}_GBT.joblib")
        px = E._model(str(B.CACHE / "models"), PROXY_SET, PROXY_TARGET, fold)
        _S[key] = (e2, px)
    return _S[key]


def f1_posthoc(
    W: Any, vctx: Any, vpos: int, R_g: np.ndarray, cost_g: float, diameter: float
) -> Tuple[np.ndarray, float, np.ndarray]:
    """f1_run.task for one seed (0), Sigma <= 29 answer (f1_answer): the lowest final local cost
    among the original answer and the refined top-3 (by post-quick-MC cost) relatives. Returns
    (answer, cost, refined orientations in f1_run's order, for the bit-identity check)."""
    T, I = _S["T"], _S["I"]
    vertices, phase = vctx.vertices, vctx.voxel.phase
    rng = np.random.default_rng([20_000 + vpos, 0])
    I.tag = "f1"
    with T.stage("f1"):
        rel, labels = csl.csl_relatives(R_g, max_sigma=29)
        post, cost = [], []
        for Rk in rel:
            Rp, cp = F1.quick_mc(W, rng, vertices, phase, Rk, diameter)
            post.append(Rp)
            cost.append(cp)
        post, cost = np.array(post), np.array(cost)
        sig = np.array([int("".join(ch for ch in lab if ch.isdigit())) for lab in labels])
        pick = set()
        for sel in (sig <= F1.SIG_SMALL, sig <= 29):
            idx = np.nonzero(sel)[0]
            pick.update(idx[np.argsort(cost[idx])[: F1.TOP_K]].tolist())
        ref_idx = sorted(pick)
        ref_R, ref_cost = [], []
        for i in ref_idx:
            R_, c_, _ = F1.refine_one(W, rng, vertices, phase, post[i])
            ref_R.append(R_)
            ref_cost.append(c_)
        idx29 = np.nonzero(sig <= 29)[0]
        top = set(idx29[np.argsort(cost[idx29])[: F1.TOP_K]].tolist())
        best_R, best_c = R_g, cost_g
        for j, i in enumerate(ref_idx):
            if i in top and ref_cost[j] < best_c:
                best_R, best_c = ref_R[j], float(ref_cost[j])
    I.tag = "none"
    return best_R, best_c, np.array(ref_R)


def _stored(pipe: str, vidx: int, variant: str) -> np.ndarray:
    """The stored (earlier task) final answer of a reconstruct pipeline, or F1's refined set."""
    if pipe == "baseline":
        return np.load(C.CACHE_DIR / "e0" / f"v{vidx}_{variant}.npz")["s0_R_final"]
    if pipe == "F1b":
        return np.load(C.CACHE_DIR / "fix" / "F1b" / f"v{vidx}_{variant}_s0.npz")["s0_R_final"]
    if pipe == "e2":
        return np.load(C.CACHE_DIR / "e2_rerank_GBT" / f"v{vidx}_{variant}.npz")["R_final"]
    if pipe == "proxy":
        return np.load(B.CACHE / "e2e" / "p_i" / f"v{vidx}_{variant}.npz")["R_final"]
    if pipe == "F1":
        return np.load(C.CACHE_DIR / "f1" / f"v{vidx}_{variant}.npz")["s0_ref_R"]
    return np.load(B.CACHE / "e2e" / "f1_p_i" / f"v{vidx}_{variant}.npz")["s0_ref_R"]


def run_recon(pipe: str, ctxd: Dict[str, Any]) -> Dict[str, Any]:
    """One timed reconstruct_voxel pipeline (setup of the proxy included in the wall time)."""
    W, vctx, vpos = ctxd["W"], ctxd["vctx"], ctxd["vpos"]
    T, I = _S["T"], _S["I"]
    vertices, phase = vctx.vertices, vctx.voxel.phase
    rec = W.rec
    e2m, pxm = _models(ctxd["fold"])
    fe_full, fe_low = ctxd["fe_full"], ctxd["fe_low"]
    n_scored = [0]

    def key_e2(level: int, cands: List[Any]) -> np.ndarray:
        with T.stage("proxy"):
            with T.stage("proxy_features"):
                X = np.stack([fe_full.features(c.orientation, vertices, phase) for c in cands])
            with T.stage("proxy_predict"):
                s = -M.score(e2m, X, True)
        n_scored[0] += len(cands)
        return s

    def key_proxy(level: int, cands: List[Any]) -> np.ndarray:
        with T.stage("proxy"):
            with T.stage("proxy_features"):
                X = E.proxy_features(fe_low, cands, vertices, phase, True)
            with T.stage("proxy_predict"):
                s = -MD.predict_score(PROXY_TARGET, pxm, X)
        n_scored[0] += len(cands)
        return s

    snap = T.snapshot()
    la = __import__("os").getloadavg()[0]
    t0 = time.perf_counter()
    if pipe == "F1b":
        FIX.apply_knobs(rec, FIX.FIXES["F1b"])
    try:
        if pipe in ("e2", "proxy"):
            with T.stage("proxy"):
                with T.stage("proxy_setimage"):
                    (fe_full if pipe == "e2" else fe_low).set_image(ctxd["keys"])
            rec.rank_key = key_e2 if pipe == "e2" else key_proxy
        with C.quiet():
            res = rec.reconstruct_voxel(vertices, phase, rng=C.run_seed(vpos, 0))
    finally:
        rec.rank_key = None
        if pipe == "F1b":
            FIX.apply_knobs(rec, {})
    wall = time.perf_counter() - t0
    R = np.asarray(res.orientation, dtype=np.float64)
    g, loc, _ = rec.last_eval_counts
    sd = PC.stage_delta(T, snap)
    return dict(
        pipe=pipe,
        wall=wall,
        R_final=R,
        cost_final=float(res.cost),
        err=PC.reduced_err_deg(R, vctx.R_true),
        evals_rec=[int(g), int(loc)],
        evals_all=I.eval_counts(snap),
        n_scored=n_scored[0],
        loadavg_start=la,
        stages=sd,
    )


def run_f1(parent: Dict[str, Any], ctxd: Dict[str, Any]) -> Dict[str, Any]:
    W, vctx, vpos = ctxd["W"], ctxd["vctx"], ctxd["vpos"]
    T, I = _S["T"], _S["I"]
    snap = T.snapshot()
    t0 = time.perf_counter()
    R, c, ref_R = f1_posthoc(
        W, vctx, vpos, parent["R_final"], parent["cost_final"], ctxd["diameter"]
    )
    wall = time.perf_counter() - t0
    name = "F1" if parent["pipe"] == "baseline" else "proxy+F1"
    return dict(
        pipe=name,
        wall=wall,
        wall_parent=parent["wall"],
        wall_total=wall + parent["wall"],
        R_final=R,
        cost_final=c,
        err=PC.reduced_err_deg(R, vctx.R_true),
        ref_R=ref_R,
        evals_all=I.eval_counts(snap),
        stages=PC.stage_delta(T, snap),
    )


def run_case(item: Tuple[Any, ...], record: bool = True) -> List[Dict[str, Any]]:
    """All pipelines of one (voxel, variant): returns the records (empty list if not recorded)."""
    from optimizer_sweep import voxel_context

    vidx, vpos, variant, rep, order_idx = item
    W = C.get_worker()
    keys = np.load(C.CACHE_DIR / "images" / f"v{vidx}_{variant}.npz")["keys"]
    C.attach(keys)
    vctx = voxel_context(W.ctx, vidx)
    phase = vctx.voxel.phase
    e0 = np.load(C.CACHE_DIR / "e0" / f"v{vidx}_{variant}.npz")
    ctxd = dict(
        W=W,
        vctx=vctx,
        vpos=vpos,
        keys=keys,
        fold=int(M.fold_of(np.array([vpos]))[0]),
        fe_full=_fe(None, phase),
        fe_low=_fe(PROXY_Q, phase),
        diameter=float(e0["s0_L3_disc_diameter"]),
    )
    _models(ctxd["fold"])
    out: List[Dict[str, Any]] = []
    for pipe in PC.rotated(RECON, order_idx):
        r = run_recon(pipe, ctxd)
        out.append(r)
        if pipe in ("baseline", "proxy"):
            out.append(run_f1(r, ctxd))
    if not record:
        return []
    for r in out:
        ref = _stored(r["pipe"], vidx, variant)
        got = r["ref_R"] if "ref_R" in r else r["R_final"]
        r["identical_to_stored"] = bool(np.array_equal(got, ref))
        r.pop("ref_R", None)
        r.update(vidx=vidx, vpos=vpos, variant=variant, rep=rep)
        r["R_final"] = np.asarray(r["R_final"]).tolist()
    return out


def task(item: Tuple[Any, ...]) -> str:
    vidx, vpos, variant, rep, order_idx, path, warm = item
    t0 = time.time()
    if not _S["warm"]:
        if warm == "full":
            run_case((vidx, vpos, variant, -1, 0), record=False)
        else:  # light: load models, build the extractors, no reconstruct run
            from optimizer_sweep import voxel_context

            W = C.get_worker()
            C.attach(np.load(C.CACHE_DIR / "images" / f"v{vidx}_{variant}.npz")["keys"])
            _fe(None, voxel_context(W.ctx, vidx).voxel.phase)
        _S["warm"] = True
    iso = PC.isolation_record(f"u0 v{vidx} {variant} rep{rep}")
    recs = run_case((vidx, vpos, variant, rep, order_idx))
    Path(path).write_text(json.dumps(dict(records=recs, isolation=iso), default=float))
    return f"u0 v{vidx} {variant} rep{rep} {time.time() - t0:.0f}s"


def work_items(a: argparse.Namespace) -> List[Tuple[Any, ...]]:
    outdir = CACHE / a.tag
    outdir.mkdir(parents=True, exist_ok=True)
    cases = []
    for vpos, v in enumerate(B.voxels()):
        for var in C.VARIANTS:
            e0 = C.CACHE_DIR / "e0" / f"v{v}_{var}.npz"
            if e0.exists() and "unbuildable" not in np.load(e0).files:
                cases.append((v, vpos, var))
    cases = cases[: 2 * N_VOX]
    its = []
    n_rep_cases = max(1, len(cases) // 10)  # repeats on ~10% of the cases
    rep_cases = cases[:: max(1, len(cases) // n_rep_cases)][:n_rep_cases]
    for rep in range(a.repeats):
        for ci, (v, vpos, var) in enumerate(cases):
            if rep > 0 and (v, vpos, var) not in rep_cases:
                continue
            path = outdir / f"v{v}_{var}_rep{rep}.json"
            if not path.exists():
                its.append((v, vpos, var, rep, ci + 7 * rep, str(path), a.warm))
    return its[: a.limit] if a.limit else its


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    ap.add_argument("cmd", choices=["run"])
    ap.add_argument("--workers", type=int, default=1)
    ap.add_argument("--repeats", type=int, default=3)
    ap.add_argument("--tag", default="w1")
    ap.add_argument("--warm", choices=["full", "light"], default="full")
    ap.add_argument("--limit", type=int, default=0)
    a = ap.parse_args()
    its = work_items(a)
    print(len(its), "tasks", flush=True)
    iso = PC.isolation_record(f"u0 start tag={a.tag} workers={a.workers}")
    PC.write_json(CACHE / a.tag / "isolation_start.json", iso)
    print(json.dumps(iso), flush=True)
    t0 = time.time()
    PC.run_pool(task, init_worker, C.worker_args(), its, a.workers, "u0")
    wall = time.time() - t0
    PC.write_json(CACHE / a.tag / "isolation_end.json", PC.isolation_record("u0 end"))
    PC.write_json(CACHE / a.tag / "wall.json", dict(wall_s=wall, tasks=len(its), workers=a.workers))
    print(f"u0 done, wall {wall:.0f}s", flush=True)


if __name__ == "__main__":
    main()
