"""Task E measurements 2 and 3: the smooth costs as the ranking key of levels 0-2.

One no-start search per voxel through the Task C probe (``candidate_probe.Probe``: the existing
recorder hooks, evaluations counted by wrapping the cost entry points). Per level it keeps the
post-quick-MC candidates (in ranked order), their strict cost, the quarter kept and the symmetry-
reduced distance of every candidate to the truth.

* Measurement 2 (``--spec off``, the baseline run): every arm of ``RANK_ARMS`` is evaluated on the
  recorded post-MC candidates of levels 0-2, so its ranking can be compared with the strict one on
  the same list (the evaluations of the search are unchanged; this isolates the ranking).
* Measurement 3 (``--spec <SmoothCostRank spec>``): the search runs with the key on, end to end.

Sets: ``e0`` = the E0 200 voxels (per-voxel clean images, seed 0, ``mc`` arm, min_sin_eta 0);
``full`` = the 30-voxel ray-study set on the clean full-sample images (Q-max 8 configuration,
seed 1, min_sin_eta 0), which contains the six interior voxels. Rules: ``old`` lenient (library
default), ``cpp`` C++ rule (CppOverlapRule 1).

  uv run python scripts/eval_reduction/smooth_rank_probe.py run --kind e0 --rule old --spec off
  uv run python scripts/eval_reduction/smooth_rank_probe.py run --kind e0 --rule old \
      --spec E1+E2 --tag e12
  uv run python scripts/eval_reduction/smooth_rank_probe.py run --kind full --voxels six ...

One npz per voxel under benchmarks/eval_reduction/E/cache/rank/<kind>_<rule>_<tag>/ (gitignored).
Wall times are contended multi-worker and are not timing claims.
"""

import argparse
import json
import multiprocessing as mp
import os
import sys
import time
from pathlib import Path
from typing import Any, Dict, List, Tuple

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")

import numpy as np

HERE = Path(__file__).resolve().parent
ICE = HERE.parents[1]
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(ICE / "scripts" / "findoptimal_robustness"))
sys.path.insert(0, str(ICE / "scripts" / "common"))
sys.path.insert(0, str(ICE / "scripts"))

OUT = ICE / "benchmarks" / "eval_reduction" / "E"
CACHE = OUT / "cache" / "rank"
VOXELS_30 = json.loads(
    (ICE / "benchmarks" / "phase_d_seed_diag" / "ray_overlap" / "voxels.json").read_text()
)
SIX = [15901, 17034, 18558, 19889, 22867, 5242]
RANK_LEVELS = 3
# arms re-scored on the baseline's post-MC lists (measurement 2); widths "auto" = from the level
RANK_ARMS: Dict[str, str] = {
    "S": "S",
    "E1s3": "E1,sigma=3",
    "E1s8": "E1,sigma=8",
    "E1a": "E1",
    "E2a": "E2",
    "E3": "E3",
    "E4t05": "E4,tau=0.05",
    "E4t25": "E4,tau=0.25",
    "E4a": "E4",
    "E12a": "E1+E2",
    "E123a": "E1+E2+E3",
    "E23a": "E2+E3",
    "E34a": "E3+E4",
    "E1234a": "E1+E2+E3+E4",
    # fixed-width twins of the landscape arms (smooth_landscape.ARMS)
    "E1s12": "E1,sigma=12",
    "E2w1": "E2,omega=1",
    "E12s3": "E1+E2,sigma=3,omega=1",
    "E123s3": "E1+E2+E3,sigma=3,omega=1",
    "E23w1": "E2+E3,omega=1",
    "E34t25": "E3+E4,tau=0.25",
}
_W: Dict[str, Any] = {}


# --------------------------------------------------------------------------------------------
def _init(kind: str) -> None:
    import torch

    torch.set_num_threads(1)
    _W["kind"] = kind
    if kind == "e0":
        import common as C

        C.init_worker(C.worker_args())
        _W["C"] = C
    else:
        import parity_common as PC

        _W["ctx"] = PC.PyContext()


def _set_rule(sp: Any, rule: str, spec: str) -> None:
    sp.cpp_overlap_rule = rule == "cpp"
    sp.smooth_rank = spec


def task(item: Tuple[Any, ...]) -> str:
    from candidate_probe import Probe
    from csl import reduced_misorientation_deg
    from icenine.smooth_cost import SmoothCost, SmoothSpec

    kind, vidx, vpos, rule, spec, path = item
    t0 = time.time()
    if kind == "e0":
        from optimizer_sweep import voxel_context

        C = _W["C"]
        W = C.get_worker()
        keys = np.load(C.CACHE_DIR / "images" / f"v{vidx}_clean.npz")["keys"]
        C.attach(keys)
        vctx = voxel_context(W.ctx, vidx)
        rec = W.rec
        Rt = np.asarray(vctx.R_true, dtype=np.float64)
        verts, phase = vctx.vertices, int(vctx.voxel.phase)
        rng = C.run_seed(vpos, 0)
        setup = rec.setup
    else:
        ctx = _W["ctx"]
        rec = ctx.rec
        Rt = ctx.truth(vidx)
        verts, phase = ctx.verts(vidx), int(ctx.mic.voxels[vidx].phase)
        rng = np.random.default_rng(1)
        setup = ctx.setup
    rec.min_sin_eta = 0.0
    sp = rec.setup.search_params
    _set_rule(sp, rule, spec)
    probe = Probe(rec, keep_raw=False)
    try:
        res = probe.run(verts, phase, rng)
    finally:
        _set_rule(sp, "old", "off")
        if kind == "e0":
            rec.min_sin_eta = 0.3
    err = float(reduced_misorientation_deg(np.asarray(res.orientation, dtype=np.float64), Rt))
    ranker = rec.last_smooth_ranker
    out: Dict[str, Any] = dict(
        voxel=np.int64(vidx),
        err_deg=np.float64(err),
        evals=np.int64(probe.total_evals),
        smooth_evals=np.int64(0 if ranker is None else ranker.n_evals),
        R_true=Rt,
        R_final=np.asarray(res.orientation, dtype=np.float64),
        cost_final=np.float64(res.cost),
        n_levels=np.int64(len(probe.levels)),
    )
    sc = None
    if spec == "off":
        sc = SmoothCost(
            setup.simulator,
            setup.detector_list,
            setup.range_map,
            setup.exp_data,
            setup.sample,
            setup.structure_list,
            eta_limit=setup.exp_setup.get_eta_limit(),
            min_sin_eta=0.0,
        )
        specs = [SmoothSpec.parse(s) for s in RANK_ARMS.values()]
    for L in probe.levels:
        lv = int(L["level"])
        Ra = np.asarray(L["R_after_ranked"], dtype=np.float64)
        d = reduced_misorientation_deg(Ra, np.broadcast_to(Rt, Ra.shape).copy())
        pre = f"L{lv}_"
        out[pre + "diameter"] = np.float64(L["diameter"])
        out[pre + "n_keep"] = np.int64(L["n_keep"])
        out[pre + "n_disc"] = np.int64(L["n_disc"])
        out[pre + "cost"] = np.asarray(L["cost_after_ranked"], dtype=np.float64)
        out[pre + "dist"] = d.astype(np.float32)
        out[pre + "screen_recost"] = np.int64(L["screen_recost_evals"])
        out[pre + "quick"] = np.int64(L["quick_evals_total"])
        dd = np.asarray(
            reduced_misorientation_deg(
                np.asarray(L["R_disc"], dtype=np.float64),
                np.broadcast_to(Rt, np.asarray(L["R_disc"]).shape).copy(),
            )
        )
        out[pre + "dist_disc"] = dd.astype(np.float32)
        if sc is not None and lv < RANK_LEVELS:
            keys = sc.evaluate_specs(Ra, verts, phase, specs, diameter_rad=float(L["diameter"]))
            for name, k in zip(RANK_ARMS, keys):
                out[pre + "key_" + name] = k.astype(np.float32)
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(path, **out)
    dt = time.time() - t0
    return f"{kind} v{vidx} {rule} {spec}: err={err:.3f} evals={probe.total_evals} {dt:.0f}s"


def voxel_list(kind: str, which: str) -> List[Tuple[int, int]]:
    if kind == "e0":
        import common as C

        info = dict(np.load(C.OUT_DIR / "voxels.npz"))
        vox = [int(v) for v in info["voxel_indices"]][: C.N_VOXELS]
        return [
            (v, p)
            for p, v in enumerate(vox)
            if (C.CACHE_DIR / "images" / f"v{v}_clean.npz").exists()
            and (C.CACHE_DIR / "e0" / f"v{v}_clean.npz").exists()
        ]
    vox = SIX if which == "six" else VOXELS_30
    return [(v, p) for p, v in enumerate(vox)]


def run_dir(kind: str, rule: str, tag: str) -> Path:
    return CACHE / f"{kind}_{rule}_{tag}"


def run(kind: str, rule: str, spec: str, tag: str, workers: int, limit: int, which: str) -> None:
    root = run_dir(kind, rule, tag)
    root.mkdir(parents=True, exist_ok=True)
    (root / "spec.txt").write_text(spec + "\n")
    items = [
        (kind, v, p, rule, spec, str(root / f"v{v}.npz"))
        for v, p in voxel_list(kind, which)[: limit or None]
        if not (root / f"v{v}.npz").exists()
    ]
    print(f"{root.name} ({spec}): {len(items)} tasks, {workers} workers", flush=True)
    t0 = time.time()
    if workers <= 1:
        _init(kind)
        for i, it in enumerate(items):
            print(f"[{i + 1}/{len(items)} {time.time() - t0:.0f}s] {task(it)}", flush=True)
    else:
        with mp.get_context("spawn").Pool(workers, initializer=_init, initargs=(kind,)) as pool:
            for i, msg in enumerate(pool.imap_unordered(task, items)):
                print(f"[{i + 1}/{len(items)} {time.time() - t0:.0f}s] {msg}", flush=True)
    (root / "run_meta.json").write_text(
        json.dumps(
            dict(
                workers=workers,
                wall_s=time.time() - t0,
                tasks=len(items),
                spec=spec,
                rule=rule,
                note="contended multi-worker wall, not a timing claim",
            )
        )
        + "\n"
    )
    (OUT / "cache" / "logs" / f"rank_{root.name}.done").write_text(f"{time.time() - t0:.0f}s\n")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    ap.add_argument("cmd", choices=["run"])
    ap.add_argument("--kind", default="e0", choices=["e0", "full"])
    ap.add_argument("--rule", default="old", choices=["old", "cpp"])
    ap.add_argument("--spec", default="off")
    ap.add_argument("--tag", default="base")
    ap.add_argument("--workers", type=int, default=6)
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--voxels", default="30", choices=["30", "six"])
    a = ap.parse_args()
    run(a.kind, a.rule, a.spec, a.tag, a.workers, a.limit, a.voxels)


if __name__ == "__main__":
    main()
