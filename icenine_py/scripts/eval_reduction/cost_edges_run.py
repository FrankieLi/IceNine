"""Task B steps 1 and 4: walk rays in orientation space from the truth and attribute every change of
the cost to per-peak events (see cost_edges.py).

Cases: the six interior full-sample voxels (clean full-sample images, Q-max 8 configuration) and the
first 30 voxels of the E0 set (per-voxel clean images). Combinations per orientation (setting x
rule):
  full   (pixel radius 0, all reflections: the complete cost)   x old, cpp, cppE, r2
  screen (pixel radius 3, Q-max 5, centroid: the level-0 screen) x old, cpp, cppE, r2
  old = lenient rule (library default), cpp = C++ overlap rule, cppE = C++ rule + CppEdgeRules,
  r2 = ray rule with the 2-pixel window.
Rays: 20 random directions per case (seed 1000 + voxel); full-sample voxels walk 0 to 5 degrees at
0.0005 degree (10,001 points); E0 voxels walk 0 to 0.5 degree at 0.0005 plus four 0.05-degree
windows at 0.5, 1.5, 3.0 and 4.95 degrees (0.0005 step), and a log-spaced profile (cost only).
Segments are differenced within themselves, never across a gap. (The plan asked for log-spaced
radii; a log grid takes several events per step at large radius, so the walks are dense instead.)

  uv run python scripts/eval_reduction/cost_edges_run.py pilot
  nohup uv run python scripts/eval_reduction/cost_edges_run.py run --kind full --workers 8 &
  nohup uv run python scripts/eval_reduction/cost_edges_run.py run --kind e0 --workers 8 &
  uv run python scripts/eval_reduction/cost_edges_slices.py   # 2-D slices through the truth

Outputs (one npz per task): benchmarks/eval_reduction/B/cache/<kind>/v<voxel>_d<dir>.npz
"""

import argparse
import json
import multiprocessing as mp
import os
import sys
import time
from pathlib import Path
from types import SimpleNamespace
from typing import Any, Dict, List, Tuple

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")

import numpy as np

HERE = Path(__file__).resolve().parent
ICE = HERE.parents[1]
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(ICE / "scripts" / "findoptimal_robustness"))
sys.path.insert(0, str(ICE / "scripts"))
import cost_edges as CE  # noqa: E402

OUT = ICE / "benchmarks" / "eval_reduction" / "B"
CACHE = OUT / "cache"
COMBOS = [
    ("full", "old"),
    ("full", "cpp"),
    ("full", "cppE"),
    ("full", "r2"),
    ("screen", "old"),
    ("screen", "cpp"),
    ("screen", "cppE"),
    ("screen", "r2"),
]
N_DIR = 20
N_E0 = 30
STEP = 0.0005
N_CAUSE = len(CE.EVENT_TYPES) + 2  # + denominator, unattributed
CAUSES = CE.EVENT_TYPES + ["denominator", "unattributed"]
HBINS = np.arange(-7.0, 0.01, 0.25)  # log10 of the per-peak change of the cost
REGIMES = [("fin", 0.0, 0.05), ("mc", 0.05, 0.5), ("grid", 0.5, 5.01)]
_ENV: Dict[str, Any] = {}


def segments(kind: str) -> Dict[str, Tuple[str, np.ndarray]]:
    """name -> (type, radii in degrees); type 'fine' is differenced, 'profile' only evaluated."""
    if kind == "full":
        return {"dense": ("fine", np.linspace(0.0, 5.0, int(round(5.0 / STEP)) + 1))}
    seg: Dict[str, Tuple[str, np.ndarray]] = {
        "dense": ("fine", np.linspace(0.0, 0.5, int(round(0.5 / STEP)) + 1))
    }
    n = int(round(0.05 / STEP)) + 1
    for s0 in (0.5, 1.5, 3.0, 4.95):
        seg[f"win{s0}"] = ("fine", np.linspace(s0, s0 + 0.05, n))
    seg["profile"] = ("profile", CE.log_radii())
    return seg


# --------------------------------------------------------------------------------------------
# environments
# --------------------------------------------------------------------------------------------
def init_full() -> None:
    import parity_common as PC

    ctx = PC.PyContext()
    s = ctx.setup
    _ENV.update(
        kind="full",
        ctx=ctx,
        base=SimpleNamespace(
            simulator=s.simulator,
            detector_list=s.detector_list,
            range_map=s.range_map,
            exp_data=s.exp_data,
            sample=s.sample,
            structure_list=s.structure_list,
            eta_limit=ctx.eta_limit,
            reg_mm=float(s.exp_setup.config_file.detector_dist_deviation),
        ),
    )


def init_e0() -> None:
    import common as C

    C.init_worker(C.worker_args())
    _ENV.update(kind="e0", C=C, current=None)


def case_data(voxel: int) -> Tuple[SimpleNamespace, np.ndarray, Any, int]:
    """(base, R_true, vertices, phase) of a voxel, with its images attached for E0."""
    if _ENV["kind"] == "full":
        ctx = _ENV["ctx"]
        return _ENV["base"], ctx.truth(voxel), ctx.verts(voxel), int(ctx.mic.voxels[voxel].phase)
    from optimizer_sweep import voxel_context

    C = _ENV["C"]
    W = C.get_worker()
    if _ENV["current"] != voxel:
        keys = np.load(C.CACHE_DIR / "images" / f"v{voxel}_clean.npz")["keys"]
        C.attach(keys)
        _ENV["current"] = voxel
    vctx = voxel_context(W.ctx, voxel)
    rs = W.rsetup
    base = SimpleNamespace(
        simulator=rs.simulator,
        detector_list=rs.detector_list,
        range_map=rs.range_map,
        exp_data=rs.exp_data,
        sample=rs.sample,
        structure_list=rs.structure_list,
        eta_limit=rs.exp_setup.get_eta_limit(),
        reg_mm=float(rs.exp_setup.config_file.detector_dist_deviation),
    )
    return base, np.asarray(vctx.R_true, dtype=np.float64), vctx.vertices, int(vctx.voxel.phase)


def e0_voxels() -> List[int]:
    import common as C

    info = dict(np.load(C.OUT_DIR / "voxels.npz"))
    vox = [int(v) for v in info["voxel_indices"]][: C.N_VOXELS]
    ok = [v for v in vox if (C.CACHE_DIR / "e0" / f"v{v}_clean.npz").exists()]
    return ok[:N_E0]


# --------------------------------------------------------------------------------------------
# the walk
# --------------------------------------------------------------------------------------------
def walk(
    cfs: List[Any],
    R_list: np.ndarray,
    verts: Any,
    phase: int,
    dims: Tuple[int, int],
    differenced: bool,
    radii: np.ndarray,
) -> List[Dict[str, Any]]:
    """Trace every orientation with every cost function; difference consecutive traces."""
    T = len(R_list)
    res = []
    for cf in cfs:
        radius = int(cf.pixel_radius)
        cost = np.zeros(T, np.float32)
        nq = np.zeros(T, np.int16)
        npres = np.zeros(T, np.int16)
        pk = np.zeros(T, np.int16)
        r: Dict[str, Any] = dict(cost=cost, nq=nq, npres=npres, pkov=pk, maxdiff=0.0)
        if differenced:
            r.update(
                parts=np.zeros((T - 1, N_CAUSE), np.float32),
                evn=np.zeros((T - 1, len(CE.EVENT_TYPES)), np.int16),
                silent=np.zeros((T - 1, len(CE.EVENT_TYPES)), np.int16),
                flipc=np.zeros((T - 1, len(CE.EVENT_TYPES) + 1), np.int8),
                den_only=np.zeros(T - 1, np.int8),
                dN=np.zeros(T - 1, np.int16),
                hist_reg=np.zeros((len(REGIMES), N_CAUSE, len(HBINS)), np.int32),
            )
        res.append(r)
    prev: List[Any] = [None] * len(cfs)
    for t in range(T):
        Rt = R_list[t].astype(np.float32)
        for c, cf in enumerate(cfs):
            tr = CE.trace_eval(cf, Rt, verts, phase)
            r = res[c]
            r["maxdiff"] = max(r["maxdiff"], CE.check_trace(tr))
            r["cost"][t] = 1.0 - tr["agg"]["quality"]
            r["nq"][t] = tr.get("n_q", 0)
            r["npres"][t] = int((tr["status"] == CE.PRESENT).sum())
            r["pkov"][t] = tr["agg"]["pkov"]
            if differenced and prev[c] is not None:
                a = prev[c]
                ev = CE.classify(a, tr, dims, int(cfs[c].pixel_radius))
                dcost, parts, flags = CE.attribute(a, tr, ev)
                k = t - 1
                for name, v in parts.items():
                    r["parts"][k, CAUSES.index(name)] = v
                for i, typ in enumerate(CE.EVENT_TYPES):
                    sel = ev == i
                    r["evn"][k, i] = int(sel.sum())
                    if sel.any():
                        r["silent"][k, i] = int((tr["q"][sel] == a["q"][sel]).sum())
                r["den_only"][k] = flags["n_den_only"]
                r["dN"][k] = flags["dN"]
                fc = CE.flip_causes(a, tr, ev)
                for name, v in fc.items():
                    r["flipc"][k, CE.EVENT_TYPES.index(name) if name != "unattributed" else -1] = v
                # per-peak jump sizes
                Nb = int(tr.get("n_q", 0))
                rmid = 0.5 * (radii[t - 1] + radii[t])
                reg = next(i for i, (_, lo, hi) in enumerate(REGIMES) if lo <= rmid < hi)
                if Nb:
                    dq = np.abs(tr["q"] - a["q"]) / Nb
                    for i in np.nonzero(dq)[0]:
                        cause = int(ev[i]) if ev[i] >= 0 else len(CE.EVENT_TYPES) + 1
                        b = int(np.clip(np.searchsorted(HBINS, np.log10(dq[i])), 0, len(HBINS) - 1))
                        r["hist_reg"][reg, cause, b] += 1
            prev[c] = tr
    return res


def run_task(item: Tuple[str, int, int]) -> str:
    kind, voxel, d = item
    path = CACHE / kind / f"v{voxel}_d{d}.npz"
    if path.exists():
        return f"{path.name} exists"
    t0 = time.time()
    base, Rt, verts, phase = case_data(voxel)
    cfs = [CE.make_cf(base, s, r) for s, r in COMBOS]
    im = base.exp_data.get_image(0, 0)
    dims = (int(im.num_rows), int(im.num_cols))
    dirs = CE.unit_directions(N_DIR, 1000 + voxel)
    out: Dict[str, Any] = {}
    maxdiff = 0.0
    for name, (typ, radii) in segments(kind).items():
        Rs = CE.rotate_about(Rt, radii, dirs[d])
        res = walk(cfs, Rs, verts, phase, dims, typ == "fine", radii)
        out[f"radii|{name}"] = radii.astype(np.float32)
        for (setting, rule), r in zip(COMBOS, res):
            maxdiff = max(maxdiff, r["maxdiff"])
            for key, val in r.items():
                if key == "maxdiff":
                    continue
                out[f"{setting}|{rule}|{key}|{name}"] = val
    out["dir"] = dirs[d]
    out["maxdiff"] = np.float64(maxdiff)
    path.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(path, **out)
    return f"{path.name} {time.time() - t0:.0f}s maxdiff {maxdiff:.1e}"


def _init(kind: str) -> None:
    import torch

    torch.set_num_threads(1)
    (init_full if kind == "full" else init_e0)()


def cases(kind: str) -> List[int]:
    if kind == "full":
        from parity_common import SIX_VOXELS

        return list(SIX_VOXELS)
    return e0_voxels()


def run(kind: str, workers: int, limit_dirs: int, limit_cases: int) -> None:
    log = CACHE / "logs"
    log.mkdir(parents=True, exist_ok=True)
    vox = cases(kind)[: limit_cases or None]
    items = [(kind, v, d) for d in range(limit_dirs or N_DIR) for v in vox]
    print(f"{kind}: {len(items)} tasks, {workers} workers", flush=True)
    t0 = time.time()
    with mp.get_context("spawn").Pool(workers, initializer=_init, initargs=(kind,)) as pool:
        for i, msg in enumerate(pool.imap_unordered(run_task, items)):
            print(f"[{i + 1}/{len(items)} {time.time() - t0:.0f}s] {msg}", flush=True)
    (log / f"{kind}.done").write_text(f"{time.time() - t0:.0f}s\n")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd", choices=["pilot", "run"])
    ap.add_argument("--kind", default="full", choices=["full", "e0"])
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--dirs", type=int, default=0, help="first N directions only")
    ap.add_argument("--cases", type=int, default=0, help="first N cases only")
    a = ap.parse_args()
    if a.cmd == "pilot":
        _init(a.kind)
        v = cases(a.kind)[0]
        t0 = time.time()
        print(run_task((a.kind, v, 0)), time.time() - t0)
        return
    run(a.kind, a.workers, a.dirs, a.cases)


if __name__ == "__main__":
    main()
