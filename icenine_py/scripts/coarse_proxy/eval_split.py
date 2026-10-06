#!/usr/bin/env python3
"""Where do the extra local evaluations of the proxy rerank (row (i)) go?

Re-runs reconstruct_voxel (seed 0, E0 images and rng) for the first 20 voxels x 2 variants twice,
baseline and proxy row (i), with the stage timer of scripts/nn_hybrid/stage_timer.py installed
(nothing in icenine/ is changed). Cost evaluations are counted per (stage, kind) with
the enclosing stage: discrete search per level, quick MC per level (MCOptimizer.optimize with
max_mc_steps == 10), FindOptimal (the other MCOptimizer.optimize calls), VarianceMinimizing, the
final evaluation, and the proxy's own low-Q costs. By-name imports are patched where they are
looked up (icenine.reconstructor.run_discrete_search_spaced). The patched runs must reproduce the
unpatched R_final (E0 seed 0; the stored proxy run), which is asserted.

  uv run python scripts/coarse_proxy/eval_split.py run --workers 10
  uv run python scripts/coarse_proxy/eval_split.py aggregate
"""

import argparse
import json
import sys
import time
from collections import defaultdict
from pathlib import Path
from typing import Any, Dict, List, Tuple

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "nn_hybrid"))
import base as B  # noqa: E402

C, M = B.C, B.M
import endtoend as E  # noqa: E402
import models as MD  # noqa: E402
import stage_timer as ST  # noqa: E402

N_VOX = 20
SET, TARGET = "lowq5+c8", "reg"


def one_run(W: Any, vctx: Any, vpos: int, rank_key: Any) -> Tuple[np.ndarray, Dict[str, Any]]:
    import icenine.reconstructor as RC
    from icenine.cost_functions import VoxelCostFunction
    from icenine.orientation_search import MCOptimizer

    timer = ST.StageTimer()
    st = dict(level=-1)
    counts: Dict[str, int] = defaultdict(int)

    def discrete_name(a: tuple, k: dict) -> str:
        st["level"] += 1
        return f"discrete_L{st['level']}"

    def mc_name(a: tuple, k: dict) -> str:
        return f"quick_mc_L{st['level']}" if k.get("max_mc_steps") == 10 else "findoptimal"

    def eval_name(a: tuple, k: dict) -> str:
        parent = timer._stack[-1][0] if timer._stack else "final_eval"
        kind = "local" if a[0].pixel_radius == 0 else "pixel3"
        counts[f"{parent}|{kind}"] += 1
        return "evaluate"

    rec = C.Recorder()
    with timer.installed():
        timer.patch(RC, "run_discrete_search_spaced", discrete_name)
        timer.patch(MCOptimizer, "optimize", mc_name)
        timer.patch(MCOptimizer, "variance_minimizing_optimize", "variance")
        timer.patch(VoxelCostFunction, "evaluate", eval_name)
        W.rec.recorder = rec
        W.rec.rank_key = None if rank_key is None else _timed(timer, rank_key)
        t0 = time.perf_counter()
        try:
            with C.quiet():
                res = W.rec.reconstruct_voxel(
                    vctx.vertices, vctx.voxel.phase, rng=C.run_seed(vpos, 0)
                )
        finally:
            W.rec.rank_key, W.rec.recorder = None, None
        wall = time.perf_counter() - t0
    arr = rec.to_arrays()
    info = dict(
        counts=dict(counts), exclusive={k: v for k, v in timer.exclusive.items()}, wall=wall,
        n_candidates=[len(arr[f"L{lv}_qmc_R"]) for lv in range(4) if f"L{lv}_qmc_R" in arr],
        n_keep=[int(arr[f"L{lv}_qmc_n_keep"]) for lv in range(4) if f"L{lv}_qmc_n_keep" in arr],
        n_findoptimal=int(len(arr.get("find_R_out", []))),
        reconstructor_counts=list(W.rec.last_eval_counts),
    )  # fmt: skip
    return np.asarray(res.orientation, float), info


def _timed(timer: Any, fn: Any) -> Any:
    def wrapped(level: int, cands: List[Any]) -> np.ndarray:
        with timer.stage("proxy"):
            return fn(level, cands)

    return wrapped


def task(item: Tuple[Any, ...]) -> str:
    vidx, vpos, variant, path = item
    W, vctx, fe = E._prep(vidx, variant, 5)
    fold = int(M.fold_of(np.array([vpos]))[0])
    model = E._model(str(B.CACHE / "models"), SET, TARGET, fold)
    vertices, phase = vctx.vertices, vctx.voxel.phase

    def proxy_key(level: int, cands: List[Any]) -> np.ndarray:
        X = E.proxy_features(fe, cands, vertices, phase, True)
        return -MD.predict_score(TARGET, model, X)

    out: Dict[str, Any] = {}
    R0, out["baseline"] = one_run(W, vctx, vpos, None)
    e0 = np.load(C.CACHE_DIR / "e0" / f"v{vidx}_{variant}.npz")
    assert np.array_equal(R0, e0["s0_R_final"]), "patched baseline run differs from E0"
    R1, out["proxy"] = one_run(W, vctx, vpos, proxy_key)
    ref = np.load(E.run_file("p_i", vidx, variant, 0))
    assert np.array_equal(R1, ref["R_final"]), "patched proxy run differs from the stored run"
    Path(path).write_text(json.dumps(out))
    return f"split voxel {vidx} {variant}"


def aggregate() -> None:
    root = B.CACHE / "eval_split"
    res: Dict[str, Any] = {}
    lines: List[str] = []
    for var in C.VARIANTS:
        files = sorted(root.glob(f"v*_{var}.json"))
        data = [json.loads(f.read_text()) for f in files]
        if not data:
            continue
        stats: Dict[str, Dict[str, float]] = {}
        for pipe in ("baseline", "proxy"):
            acc: Dict[str, float] = defaultdict(float)
            for d in data:
                for k, v in d[pipe]["counts"].items():
                    acc[k] += v / len(data)
                acc["wall"] += d[pipe]["wall"] / len(data)
                acc["n_findoptimal"] += d[pipe]["n_findoptimal"] / len(data)
                for lv, n in enumerate(d[pipe]["n_candidates"]):
                    acc[f"n_cand_L{lv}"] += n / len(data)
            stats[pipe] = dict(acc)
        keys = sorted(
            {k for p in stats.values() for k in p if "|" in k},
            key=lambda k: (k.split("|")[0], k.split("|")[1]),
        )
        name = "realistic" if var == "all" else "clean"
        lines.append(
            f"  -- {name}, {len(data)} voxels, mean per run (cost evaluations; proxy's own low-Q "
            f"costs are listed under 'proxy') --"
        )
        lines.append(f"  {'stage|kind':26s} {'baseline':>9s} {'proxy':>9s} {'change':>9s}")
        for k in keys:
            b, p = stats["baseline"].get(k, 0.0), stats["proxy"].get(k, 0.0)
            lines.append(f"  {k:26s} {b:9.0f} {p:9.0f} {p - b:+9.0f}")
        for k in ("n_findoptimal", "n_cand_L0", "n_cand_L1", "n_cand_L2", "n_cand_L3"):
            b, p = stats["baseline"].get(k, 0.0), stats["proxy"].get(k, 0.0)
            lines.append(f"  {k:26s} {b:9.1f} {p:9.1f} {p - b:+9.1f}")
        for kind in ("local", "pixel3"):
            tot = {
                pipe: sum(
                    v for k, v in stats[pipe].items() if k.endswith("|" + kind) and "proxy" not in k
                )
                for pipe in stats
            }
            lines.append(
                f"  {'total ' + kind + ' (reconstructor)':26s} {tot['baseline']:9.0f} "
                f"{tot['proxy']:9.0f} {tot['proxy'] - tot['baseline']:+9.0f}"
            )
        b, p = stats["baseline"]["wall"], stats["proxy"]["wall"]
        lines.append(f"  {'wall (s, 10 workers)':26s} {b:9.2f} {p:9.2f} {p - b:+9.2f}")
        res[var] = stats
    out = dict(
        n_voxels=N_VOX, lines=lines, stats=res,
        note="patched runs reproduce the unpatched R_final (asserted per run)",
    )  # fmt: skip
    C.save_json(B.OUT / "eval_split.json", out)
    print("\n".join(lines))


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    ap.add_argument("cmd", choices=["run", "aggregate"])
    ap.add_argument("--workers", type=int, default=10)
    ap.add_argument("--limit", type=int, default=0)
    a = ap.parse_args()
    if a.cmd == "aggregate":
        return aggregate()
    root = B.CACHE / "eval_split"
    root.mkdir(parents=True, exist_ok=True)
    its = [
        (v, vpos, var, str(root / f"v{v}_{var}.json"))
        for vpos, v in enumerate(B.voxels()[:N_VOX])
        for var in C.VARIANTS
        if not (root / f"v{v}_{var}.json").exists()
    ]
    if a.limit:
        its = its[: a.limit]
    print(len(its), "tasks", flush=True)
    C.run_pool(task, its, a.workers, "split")


if __name__ == "__main__":
    main()
