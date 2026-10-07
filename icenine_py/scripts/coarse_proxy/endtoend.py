#!/usr/bin/env python3
"""End to end: the proxy as the pruning key of reconstruct_voxel (the `rank_key` hook).

At the end of every level (after the quick MC) the candidates are ordered by the proxy instead of by
the post-quick-MC local cost, then the best `keep_fraction` are kept. Full reconstruct_voxel runs
with the E0 images and rng (`C.run_seed(vpos, seed)`), so a run differs from its E0 twin only
through the key (and keep_fraction). The model of the fold NOT containing the voxel's grain is used.
Per run: R_final, cost_final, runtime, global / local evaluations (the reconstructor's own count,
which excludes the proxy's cost evaluations), the number of proxy-scored candidates, the proxy
seconds and low-Q cost evaluations, and the recorded per-level candidates (for `harvest`).

  uv run python scripts/coarse_proxy/endtoend.py run --tag p5 --set lowq5+c8 --target reg
  uv run python scripts/coarse_proxy/endtoend.py harvest --tag p5     # features of the proxy
        run's own candidates and the offline recall on them (domain-shift check)
"""

import argparse
import sys
import time
from pathlib import Path
from typing import Any, Dict, List, Tuple

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import base as B  # noqa: E402

C, M = B.C, B.M
import features as F  # noqa: E402
import models as MD  # noqa: E402

_STATE: Dict[Any, Any] = {}


def run_dir(tag: str) -> Path:
    return B.CACHE / "e2e" / tag


def run_file(tag: str, v: int, var: str, seed: int) -> Path:
    """Seed 0 uses the plain name (the input format of f1_run.py --source)."""
    return run_dir(tag) / (f"v{v}_{var}.npz" if seed == 0 else f"v{v}_{var}_s{seed}.npz")


def _prep(vidx: int, variant: str, q: int) -> Tuple[Any, Any, Any]:
    from optimizer_sweep import voxel_context

    W = C.get_worker()
    keys = np.load(C.CACHE_DIR / "images" / f"v{vidx}_{variant}.npz")["keys"]
    C.attach(keys)
    vctx = voxel_context(W.ctx, vidx)
    key = ("fe", q, vctx.voxel.phase)
    if key not in _STATE:
        _STATE[key] = F.FeatureExtractor(W.local_fn, W.ctx.geo, vctx.voxel.phase, q_max=float(q))
    fe = _STATE[key]
    fe.set_image(keys)
    return W, vctx, fe


def _model(models_dir: str, name: str, target: str, fold: int) -> Any:
    import joblib

    key = (models_dir, name, target, fold)
    if key not in _STATE:
        p = Path(models_dir) / f"{name}__{target}__fold{fold}.joblib"
        _STATE[key] = joblib.load(p)
    return _STATE[key]


def proxy_features(
    fe: Any, cands: List[Any], vertices: Any, phase: int, c8: bool, batched: bool = False
) -> np.ndarray:
    """Low-Q feature matrix of the candidates. `batched` uses `FeatureExtractor.features_batch`
    (identical values, tested; the per-candidate loop stays the default and the reference)."""
    if batched:
        X = fe.features_batch(np.stack([c.orientation for c in cands]), vertices, phase)
    else:
        X = np.stack([fe.features(c.orientation, vertices, phase) for c in cands])
    if c8:
        X = np.hstack([X, np.array([[float(c.cost)] for c in cands])])  # free Q8 local cost
    return X


def task_run(item: Tuple[Any, ...]) -> str:
    vidx, vpos, variant, seed, path, cfg = item
    t_start = time.time()
    kind, q, c8 = MD.SETS[cfg["set"]]
    assert kind == "lowq", "the end-to-end proxy needs a low-Q feature set"
    W, vctx, fe = _prep(vidx, variant, q)
    fold = int(M.fold_of(np.array([vpos]))[0])
    model = _model(cfg["models_dir"], cfg["set"], cfg["target"], fold)
    vertices, phase = vctx.vertices, vctx.voxel.phase
    st: Dict[str, Any] = dict(n=0, sec=0.0, calls=[])

    def rank_key(level: int, cands: List[Any]) -> np.ndarray:
        t0 = time.perf_counter()
        X = proxy_features(fe, cands, vertices, phase, c8, cfg.get("batched", False))
        key = -MD.predict_score(cfg["target"], model, X)
        st["n"] += len(cands)
        st["calls"].append((level, len(cands)))  # batch size of this rank_key call
        st["sec"] += time.perf_counter() - t0
        return key

    rec = C.Recorder()
    W.rec.rank_key, W.rec.recorder = rank_key, rec
    prev_keep = W.rec.keep_fraction
    W.rec.keep_fraction = cfg["keep"]
    t0 = time.perf_counter()
    try:
        with C.quiet():
            res = W.rec.reconstruct_voxel(vertices, phase, rng=C.run_seed(vpos, seed))
    finally:
        W.rec.rank_key, W.rec.recorder, W.rec.keep_fraction = None, None, prev_keep
    dt = time.perf_counter() - t0
    g, loc, _ = W.rec.last_eval_counts
    arr = rec.to_arrays()
    keep = {
        k: v
        for k, v in arr.items()
        if (k.startswith("L") and (k.endswith("_qmc_R") or k.endswith("_qmc_cost")))
        or k.endswith("_qmc_n_keep")
        or k in ("find_R_out", "find_cost")
    }
    np.savez_compressed(
        path, R_final=np.asarray(res.orientation, float), cost_final=float(res.cost), runtime=dt,
        evals_global=g, evals_local=loc, n_scored=st["n"], proxy_seconds=st["sec"],
        n_proxy_cost_evals=2 * st["n"], call_levels=np.array([c[0] for c in st["calls"]], int),
        call_sizes=np.array([c[1] for c in st["calls"]], int), R_true=vctx.R_true, **keep,
    )  # fmt: skip
    return f"{cfg['tag']} voxel {vidx} {variant} s{seed} {time.time() - t_start:.0f}s"


def task_harvest(item: Tuple[Any, ...]) -> str:
    """Features, errors and scores of the candidates of one proxy run, plus their y_bcost labels."""
    import labels as L

    vidx, vpos, variant, seed, src, path, cfg = item
    kind, q, c8 = MD.SETS[cfg["set"]]
    W, vctx, fe = _prep(vidx, variant, q)
    fold = int(M.fold_of(np.array([vpos]))[0])
    model = _model(cfg["models_dir"], cfg["set"], cfg["target"], fold)
    vertices, phase = vctx.vertices, vctx.voxel.phase
    d = np.load(src)
    lab = dict(np.load(B.CACHE / "labels" / f"v{vidx}_{variant}.npz"))
    Rs, costs, lvl = [], [], []
    for lv in range(4):
        if f"L{lv}_qmc_R" in d.files:
            Rs.append(d[f"L{lv}_qmc_R"])
            costs.append(d[f"L{lv}_qmc_cost"])
            lvl.append(np.full(len(costs[-1]), lv))
    R, cost, level = np.concatenate(Rs), np.concatenate(costs), np.concatenate(lvl)
    cands = [type("Cand", (), dict(orientation=r, cost=c))() for r, c in zip(R, cost)]
    X = proxy_features(fe, cands, vertices, phase, c8)
    err = C.err_deg(R.astype(np.float64), d["R_true"])
    e2 = dict(R=R.astype(np.float64), err=err, source=np.zeros(len(R), int), cost=cost)
    y, cat, _ = L.build_labels(e2, lab)
    np.savez_compressed(
        path, X=X, err=err, level=level, y_bcost=y, category=cat, cost=cost, vpos=vpos,
        variant=C.VARIANTS.index(variant), seed=seed, fold=fold,
        score=MD.predict_score(cfg["target"], model, X),
    )  # fmt: skip
    return f"harvest {vidx} {variant} s{seed}"


def cfg_from(a: argparse.Namespace) -> Dict[str, Any]:
    return dict(
        tag=a.tag, set=a.set, target=a.target, keep=a.keep,
        models_dir=str(B.CACHE / a.models_dir), batched=a.batched,
    )  # fmt: skip


def work_items(a: argparse.Namespace, harvest: bool) -> List[Tuple[Any, ...]]:
    cfg = cfg_from(a)
    out = []
    odir = run_dir(a.tag + ("_harvest" if harvest else ""))
    odir.mkdir(parents=True, exist_ok=True)
    for vpos, v in enumerate(B.voxels()):
        for var in C.VARIANTS:
            if not (C.CACHE_DIR / "e0" / f"v{v}_{var}.npz").exists():
                continue
            for s in a.seeds:
                src = run_file(a.tag, v, var, s)
                if harvest:
                    dst = odir / src.name
                    if src.exists() and not dst.exists():
                        out.append((v, vpos, var, s, str(src), str(dst), cfg))
                elif not src.exists():
                    out.append((v, vpos, var, s, str(src), cfg))
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    ap.add_argument("cmd", choices=["run", "harvest"])
    ap.add_argument("--tag", required=True)
    ap.add_argument("--set", default="lowq5+c8")
    ap.add_argument("--target", default="reg", choices=MD.TARGETS)
    ap.add_argument("--keep", type=float, default=0.25)
    ap.add_argument("--models-dir", default="models")
    ap.add_argument("--seeds", type=int, nargs="+", default=[0])
    ap.add_argument("--workers", type=int, default=10)
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument(
        "--batched", action="store_true", help="batched low-Q feature pass (identical features)"
    )
    a = ap.parse_args()
    its = work_items(a, a.cmd == "harvest")
    if a.limit:
        its = its[: a.limit]
    print(len(its), "tasks", flush=True)
    C.run_pool(task_run if a.cmd == "run" else task_harvest, its, a.workers, a.tag)


if __name__ == "__main__":
    main()
