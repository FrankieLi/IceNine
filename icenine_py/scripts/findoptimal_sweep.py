#!/usr/bin/env python3
"""
Multi-level adaptive reconstruction (FindOptimal) on the single-voxel perturbation-sweep cases.

Uses the same 50 voxels and the same per-voxel detector images as scripts/optimizer_sweep.py
(clean and realistic `all`; pixel sets built by its helpers, never a whole-sample simulation) and
the project's AdaptiveVoxelReconstructor with the settings of
Examples/Example2.ThreeVoxels/ConfigFiles/ReconstructQ8.config (SearchParameters.from_config:
5 deg grid radius, local resolution 0..3 = 4 levels, 200 MC steps, 2 successive restarts, MC radius
scale 0.4, max convergence cost 1e-4, 30 discrete candidates; Q_max 8; the geometry files of the
two examples are identical), the ManyGrains physics (Q_max 8) and fundamental-zone file, and the
net's eligibility filter |sin eta| >= 0.3 (AdaptiveVoxelReconstructor(min_sin_eta=...)).

Experiment A (global; no starting guess, so no r): reconstruct_voxel from scratch, once per voxel
  and variant (50 x 2). The image of a voxel is the one the optimizer sweep gave its case
  (smallest radius, direction 0; at r = 0.05 deg the windows sit at the truth). Seeded per voxel.
Experiment B (local; depends on r): AdaptiveVoxelReconstructor.refine_from_candidates, i.e.
  FindOptimal (full MC, convergence check) + VarianceMinimizing + the final evaluation exactly as
  reconstruct_voxel runs them after its levels, with the case's perturbed nominal as the single
  candidate. Its search box is fixed by SearchParameters (max(d/3, 0.2 deg) / 2^min_res with
  d = 5 deg / 1.5^4 -> 0.329 deg) and does NOT depend on r: it is not told r.

Error metric: misorientation angle to the true .mic orientation reduced by cubic symmetry (the
global search can return a symmetry-equivalent orientation); `err_plain` (no reduction) is stored
as well.

Usage (from icenine_py/):
  uv run python scripts/findoptimal_sweep.py pilot
  uv run python scripts/findoptimal_sweep.py run-a --workers 10
  uv run python scripts/findoptimal_sweep.py run-b --workers 10 [--n-dirs 20]
  uv run python scripts/findoptimal_sweep.py summarize
"""

import argparse
import contextlib
import io
import json
import os
import sys
import time
from pathlib import Path
from types import SimpleNamespace
from typing import Any, Dict, Iterator, List, Optional, Sequence, Tuple

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")

import numpy as np
import torch
from scipy.spatial.transform import Rotation

ICENINE_PY = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(Path(__file__).parent))
sys.path.insert(0, str(ICENINE_PY / "benchmarks"))

import optimizer_sweep as osw  # noqa: E402
import perturbation_sweep as ps  # noqa: E402

VARIANTS = ps.VARIANTS
OUT_DIR = ICENINE_PY / "benchmarks" / "toy_orientation_sweep"
CACHE_DIR = ICENINE_PY / "scripts" / "findoptimal_sweep_cache"
RECON_CONFIG = (
    ICENINE_PY.parent / "Examples" / "Example2.ThreeVoxels" / "ConfigFiles" / "ReconstructQ8.config"
)
LOCK_ON_DEG = 2.0  # a level's best candidate "is" the final answer within this angle


# ---------------------------------------------------------------------------
# Symmetry-reduced error
# ---------------------------------------------------------------------------

_CUBIC: Optional[np.ndarray] = None


def cubic_rotations() -> np.ndarray:
    """The 24 proper rotations of the cubic point group, (24, 3, 3), from the project's symmetry
    module (icenine.symmetry.create_cubic_symmetry)."""
    global _CUBIC
    if _CUBIC is None:
        from icenine.symmetry import create_cubic_symmetry

        mats = [np.array(m) for m in create_cubic_symmetry(4.0).get_rotation_matrices()]
        _CUBIC = np.stack([m for m in mats if np.linalg.det(m) > 0])
        assert _CUBIC.shape == (24, 3, 3)
    return _CUBIC


def misorientation_deg(R_est: np.ndarray, R_true: np.ndarray, reduce: bool = True) -> np.ndarray:
    """Angle (deg) of R_est S R_true^T minimised over the 24 cubic operators S (reduce=True; the
    crystal-frame symmetry of orientations, R -> R S), or of R_est R_true^T (reduce=False).
    Batched over leading axes of R_est (R_true broadcastable); the error of the sweep's other
    methods (perturbation_sweep.estimate_error_deg) is the reduce=False value."""
    R_est = np.asarray(R_est, dtype=np.float64)
    R_true = np.asarray(R_true, dtype=np.float64)
    ops = cubic_rotations() if reduce else np.eye(3)[None]
    M = (R_est[..., None, :, :] @ ops) @ np.swapaxes(R_true, -1, -2)[..., None, :, :]
    ang = np.degrees(Rotation.from_matrix(M.reshape(-1, 3, 3)).magnitude()).reshape(M.shape[:-2])
    return ang.min(axis=-1)


# ---------------------------------------------------------------------------
# Worker setup and the cases
# ---------------------------------------------------------------------------

_W: Optional[SimpleNamespace] = None


def init_worker(args_dict: Dict[str, Any]) -> None:
    """The optimizer sweep's physics setup plus an AdaptiveVoxelReconstructor on it."""
    global _W
    osw.init_worker(args_dict)
    ctx = osw._CTX
    assert ctx is not None
    from icenine.config_file import ConfigFile
    from icenine.orientation_search import SearchParameters
    from icenine.reconstructor import AdaptiveVoxelReconstructor, ReconstructionSetup
    from icenine.sampling import load_fundamental_zone_file

    mic, sample, detector_list, range_map, exp_setup, simulator, structure_list, _ = ctx.setup
    config = ConfigFile.from_file(str(RECON_CONFIG))
    fz = load_fundamental_zone_file(str(ctx.example_dir / "DataFiles" / "MyFZ.dat"))
    rsetup = ReconstructionSetup(
        config=config,
        exp_setup=exp_setup,
        exp_data=None,  # type: ignore[arg-type]  # attached per case
        fz_orientations=fz,
        search_params=SearchParameters.from_config(config),
        simulator=simulator,
        detector_list=detector_list,
        range_map=range_map,
        sample=sample,
        structure_list=structure_list,
    )
    rec = AdaptiveVoxelReconstructor(rsetup, min_sin_eta=ctx.args.min_sin_eta)
    _W = SimpleNamespace(ctx=ctx, rsetup=rsetup, rec=rec, local_fn=rec._make_local_cost_fn())


def case_images(
    ctx: SimpleNamespace,
    vidx: int,
    vpos: int,
    ri: int,
    dirs: Sequence[int],
    variants: Sequence[str],
    ref_nroi: np.ndarray,
    ref_fail: np.ndarray,
) -> Iterator[Tuple[int, int, np.ndarray, np.ndarray, np.ndarray]]:
    """Yield (variant index, direction j, pixel keys of the case's detector images, perturbed
    nominal R_nom0[j], R_true), built exactly as osw.sweep_task builds them (asserting the case
    is aligned with the network sweep)."""
    a = ctx.args
    vctx = osw.voxel_context(ctx, vidx)
    D = a.n_dirs
    r = a.radii[ri]
    sigma_comp = a.neighbor_sigma_deg / np.sqrt(3.0)
    delta0, draws = ps.case_draws(
        a.sweep_seed, vidx, ri, r, D, len(vctx.sources), sigma_comp, a.neighbor_p
    )
    R_nom0 = ps.perturbed_nominal(vctx.R_true, delta0)
    prep1: List[Optional[Any]] = []
    reason1 = np.zeros(D, dtype=np.int8)
    for j in range(D):
        p, why = ps.prepare_nominal(ctx, vidx, R_nom0[j])
        prep1.append(p)
        reason1[j] = why
    assert (np.array([p.n if p is not None else 0 for p in prep1]) == ref_nroi).all()
    for vi, variant in enumerate(VARIANTS):
        if variant not in variants:
            continue
        seed = a.realism_seed + 1000003 * vpos + 1009 * ri
        layers: Dict[str, torch.Tensor] = {}
        b1 = ps.render_batch(prep1, delta0, draws, vctx.sources, variant, seed, a, layers=layers)
        have = np.array([p is not None for p in prep1])
        ok1 = have & (b1["n_present"].numpy() >= a.min_present)
        fail = np.where(reason1 > 0, reason1, np.where(ok1, 0, 4)).astype(np.int8)
        assert (fail == ref_fail[:, vi]).all(), "pass-1 failures differ from the sweep"
        for j in dirs:
            if not have[j]:
                continue
            edit = None
            if variant != "clean":
                p = prep1[j]
                n = p.n
                edit = osw.realism_edit(
                    layers["clean"][j, :n].numpy(),
                    layers["dis"][j, :n].numpy(),
                    b1["windows"][j, :n].numpy(),
                    p.spec,
                    p.obs.det_idx.numpy(),
                    a.frame_half_width,
                    ctx.geo,
                )
            keys = osw.case_image_keys(ctx, vctx, variant, draws[j], edit)
            yield vi, j, keys, R_nom0[j], vctx.R_true


def attach_images(keys: np.ndarray) -> None:
    """Make `keys` the experimental data of the reconstructor and of the shared local cost fn."""
    assert _W is not None
    ctx = _W.ctx
    data = osw.PixelSetData(osw.group_pixels(keys, ctx.geo), ctx.geo, ctx.zeros)
    _W.rsetup.exp_data = data
    _W.local_fn.exp_data = data


# ---------------------------------------------------------------------------
# Experiment A: the full multi-level reconstruction
# ---------------------------------------------------------------------------

N_LEVELS = 4


def a_task(item: Tuple[int, int, int, str, np.ndarray, np.ndarray]) -> Tuple[int, int, float]:
    """One (voxel, variant) full reconstruct_voxel run; caches the result."""
    vidx, vpos, vi, path, ref_nroi, ref_fail = item
    assert _W is not None
    ctx, rec = _W.ctx, _W.rec
    t_start = time.time()
    variant = VARIANTS[vi]
    got = None
    for v_i, j, keys, R_nom, R_true in case_images(
        ctx, vidx, vpos, 0, [0], [variant], ref_nroi, ref_fail
    ):
        got = (keys, R_true)
    assert got is not None, "direction 0 of the smallest radius could not be built"
    keys, R_true = got
    attach_images(keys)
    vctx = osw.voxel_context(ctx, vidx)
    voxel, vertices = vctx.voxel, vctx.vertices
    rng = np.random.default_rng(10_000 + vpos)
    buf = io.StringIO()
    t0 = time.perf_counter()
    with contextlib.redirect_stdout(buf):
        res = rec.reconstruct_voxel(vertices, voxel.phase, rng=rng)
    dt = time.perf_counter() - t0
    g, loc, tot = rec.last_eval_counts
    lv = rec.last_level_best
    level_R = np.full((N_LEVELS, 3, 3), np.nan)
    level_cost = np.full(N_LEVELS, np.nan)
    for lvl, R, c in lv:
        level_R[lvl], level_cost[lvl] = R, c
    cost_true = _W.local_fn.evaluate(R_true.astype(np.float32), vertices, voxel.phase).cost
    fo = rec.last_find_optimal
    np.savez(
        path,
        R_final=np.asarray(res.orientation, dtype=np.float64),
        R_true=R_true,
        cost_final=float(res.cost),
        cost_true=float(cost_true),
        runtime=dt,
        evals_global=g,
        evals_local=loc,
        level_R=level_R,
        level_cost=level_cost,
        find_winner=int(fo.get("winner_index", -2)),
        find_evaluated=int(fo.get("n_evaluated", -1)),
        find_candidates=int(fo.get("n_candidates", -1)),
        find_converged=bool(fo.get("converged", False)),
        n_image_pixels=len(keys),
    )
    return vidx, vi, time.time() - t_start


# ---------------------------------------------------------------------------
# Experiment B: FindOptimal + VarianceMinimizing from the perturbed start
# ---------------------------------------------------------------------------


def b_task(item: Tuple[int, int, int, str, np.ndarray, np.ndarray, int]) -> Tuple[int, int, float]:
    """One (voxel, radius): refine_from_candidates from each direction's perturbed nominal, for
    both variants; caches the result."""
    from icenine.orientation_search import MCOptimizer, SearchCandidate

    vidx, vpos, ri, path, ref_nroi, ref_fail, n_dirs = item
    assert _W is not None
    ctx, rec, lf = _W.ctx, _W.rec, _W.local_fn
    t_start = time.time()
    a = ctx.args
    vctx = osw.voxel_context(ctx, vidx)
    voxel, vertices = vctx.voxel, vctx.vertices
    D, V = a.n_dirs, len(VARIANTS)
    out = dict(
        R_final=np.full((D, V, 3, 3), np.nan),
        cost_final=np.full((D, V), np.nan),
        cost_true=np.full((D, V), np.nan),
        runtime=np.full((D, V), np.nan),
        evals=np.zeros((D, V), dtype=np.int32),
        converged=np.zeros((D, V), dtype=bool),
        ran=np.zeros((D, V), dtype=bool),
    )
    for vi, j, keys, R_nom, R_true in case_images(
        ctx, vidx, vpos, ri, list(range(n_dirs)), VARIANTS, ref_nroi, ref_fail
    ):
        attach_images(keys)
        seed = int(a.realism_seed + 1000003 * vpos + 1009 * ri + 7919 * j + 17 * vi)
        mc = MCOptimizer(
            cost_fn=lf, voxel_vertices=vertices, phase_index=voxel.phase,
            rng=np.random.default_rng(seed),
        )  # fmt: skip
        n0 = lf.eval_count
        t0 = time.perf_counter()
        with contextlib.redirect_stdout(io.StringIO()):
            res = rec.refine_from_candidates(
                [SearchCandidate(orientation=R_nom.astype(np.float32), cost=1.0)],
                vertices,
                voxel.phase,
                local_cost_fn=lf,
                mc_optimizer=mc,
            )
        dt = time.perf_counter() - t0
        n_ev = lf.eval_count - n0
        out["R_final"][j, vi] = res.orientation
        out["cost_final"][j, vi] = res.cost
        out["cost_true"][j, vi] = lf.evaluate(R_true.astype(np.float32), vertices, voxel.phase).cost
        out["runtime"][j, vi] = dt
        out["evals"][j, vi] = n_ev
        out["converged"][j, vi] = bool(rec.last_find_optimal.get("converged", False))
        out["ran"][j, vi] = True
    np.savez(path, **out)
    return vidx, ri, time.time() - t_start


# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------


def worker_args(args: argparse.Namespace) -> Tuple[Dict[str, Any], Dict[str, np.ndarray]]:
    cfg, sweep = osw.load_sweep_config(OUT_DIR / "perturbation_sweep_raw.npz")
    wargs = dict(cfg)
    wargs.update(methods=[], opt_dirs=0, mc_steps=0, mc_restarts=0, mc_step_frac=0)
    return wargs, sweep


def run_pool(func: Any, todo: List[Any], workers: int, wargs: Dict[str, Any], label: str) -> None:
    import multiprocessing as mp

    t0 = time.time()
    mpctx = mp.get_context("spawn")
    with mpctx.Pool(workers, initializer=init_worker, initargs=(wargs,)) as pool:
        for k, (vidx, ix, secs) in enumerate(pool.imap_unordered(func, todo), 1):
            print(
                f"  [{label} {k}/{len(todo)}] voxel {vidx} #{ix} done in {secs:.0f}s "
                f"({time.time() - t0:.0f}s total)",
                flush=True,
            )


def a_items(args: argparse.Namespace, sweep: Dict[str, np.ndarray], cache: Path) -> List[Any]:
    voxels = [int(v) for v in sweep["voxel_indices"]][: args.n_voxels or None]
    return [
        (v, vpos, vi, str(cache / f"A_v{v}_{VARIANTS[vi]}.npz"), sweep["n_roi"][vpos, 0],
         sweep["fail_pass1"][vpos, 0])
        for vpos, v in enumerate(voxels)
        for vi in range(len(VARIANTS))
    ]  # fmt: skip


def b_items(args: argparse.Namespace, sweep: Dict[str, np.ndarray], cache: Path) -> List[Any]:
    voxels = [int(v) for v in sweep["voxel_indices"]][: args.n_voxels or None]
    radii_idx = args.radii_idx if args.radii_idx else list(range(len(sweep["radii"])))
    return [
        (v, vpos, ri, str(cache / f"B_v{v}_r{ri}_d{args.n_dirs}.npz"), sweep["n_roi"][vpos, ri],
         sweep["fail_pass1"][vpos, ri], args.n_dirs)
        for vpos, v in enumerate(voxels)
        for ri in radii_idx
    ]  # fmt: skip


def do_run(args: argparse.Namespace, which: str) -> None:
    cache = Path(args.cache_dir).resolve()
    cache.mkdir(parents=True, exist_ok=True)
    wargs, sweep = worker_args(args)
    items = a_items(args, sweep, cache) if which == "A" else b_items(args, sweep, cache)
    todo = [it for it in items if not Path(it[3]).exists()]
    print(
        f"experiment {which}: {len(items)} tasks ({len(items) - len(todo)} cached), "
        f"{args.workers} workers",
        flush=True,
    )
    t0 = time.time()
    run_pool(a_task if which == "A" else b_task, todo, args.workers, wargs, which)
    print(f"experiment {which} finished; wall {time.time() - t0:.0f}s", flush=True)
    assemble(which, items, sweep, wargs, args)


def assemble(
    which: str,
    items: List[Any],
    sweep: Dict[str, np.ndarray],
    wargs: Dict[str, Any],
    args: argparse.Namespace,
) -> None:
    """Stack the cache files of one experiment into <out-dir>/findoptimal_{a,b}_raw.npz."""
    if any(not Path(it[3]).exists() for it in items):
        print("not all tasks cached; nothing assembled")
        return
    parts = [np.load(it[3]) for it in items]
    raw: Dict[str, np.ndarray] = {}
    if which == "A":
        for k in parts[0].files:
            raw[k] = np.stack([p[k] for p in parts])
        raw["voxel_indices"] = np.array([it[0] for it in items])
        raw["variant_index"] = np.array([it[2] for it in items])
    else:
        voxels = sorted({it[0] for it in items})
        ris = sorted({int(Path(it[3]).stem.split("_r")[1].split("_")[0]) for it in items})
        by = {
            (it[0], int(Path(it[3]).stem.split("_r")[1].split("_")[0])): p
            for it, p in zip(items, parts)
        }
        for k in parts[0].files:
            raw[k] = np.stack([np.stack([by[(v, ri)][k] for ri in ris]) for v in voxels])
        raw["voxel_indices"] = np.array(voxels)
        raw["radii"] = np.array([sweep["radii"][i] for i in ris])
        raw["radii_idx"] = np.array(ris)
        raw["n_dirs_run"] = np.array(args.n_dirs)
    raw["variants"] = np.array(VARIANTS)
    raw["config_json"] = np.array(json.dumps(wargs, default=str))
    name = OUT_DIR / f"findoptimal_{which.lower()}_raw.npz"
    np.savez_compressed(name, **raw)
    print(f"saved {name}", flush=True)


def do_pilot(args: argparse.Namespace) -> None:
    """Timing pilot: experiment A on 2 voxels (both variants) and B on 2 voxels at the smallest
    and largest radius (all 20 directions, both variants); prints the projected full wall time."""
    cache = Path(args.cache_dir).resolve() / "pilot"
    cache.mkdir(parents=True, exist_ok=True)
    wargs, sweep = worker_args(args)
    args.n_voxels, args.n_dirs, args.radii_idx = 2, 20, [0, 9]
    ai, bi = a_items(args, sweep, cache), b_items(args, sweep, cache)
    t0 = time.time()
    run_pool(a_task, [it for it in ai if not Path(it[3]).exists()], args.workers, wargs, "A")
    ta = time.time() - t0
    t0 = time.time()
    run_pool(b_task, [it for it in bi if not Path(it[3]).exists()], args.workers, wargs, "B")
    tb = time.time() - t0
    ra = [float(np.load(it[3])["runtime"]) for it in ai]
    rb = np.concatenate([np.load(it[3])["runtime"][np.load(it[3])["ran"]] for it in bi])
    print(
        f"PILOT A: {len(ai)} runs, wall {ta:.0f}s, per-run reconstruct_voxel runtime "
        f"{np.round(ra, 1).tolist()} s",
        flush=True,
    )
    print(
        f"PILOT B: {len(bi)} tasks ({len(rb)} cases), wall {tb:.0f}s, per-case refine runtime "
        f"median {np.median(rb):.2f} s mean {rb.mean():.2f} s max {rb.max():.2f} s",
        flush=True,
    )
    n_full = 50 * 10 * 20 * 2
    per_case = tb / len(rb) * min(args.workers, len(bi))  # wall per case at full occupancy
    print(
        f"PROJECTION B, full 50 voxels x 10 radii x 20 dirs x 2 variants = {n_full} cases: "
        f"~{n_full * rb.mean() / args.workers / 3600:.2f} h refine time on {args.workers} workers "
        f"(+ image preparation; pilot wall/case x workers = {per_case:.2f} s -> "
        f"{n_full * per_case / args.workers / 3600:.2f} h)",
        flush=True,
    )
    print(
        f"PROJECTION A, 100 runs: ~{100 * np.mean(ra) / args.workers / 3600:.2f} h on "
        f"{args.workers} workers",
        flush=True,
    )


def build_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    for name in ("pilot", "run-a", "run-b"):
        r = sub.add_parser(name)
        r.add_argument("--cache-dir", default=str(CACHE_DIR))
        r.add_argument("--workers", type=int, default=10)
        r.add_argument("--n-voxels", type=int, default=0, help="first N sweep voxels (0 = all)")
        if name == "run-b":
            r.add_argument("--n-dirs", type=int, default=20, help="directions per radius")
            r.add_argument("--radii-idx", type=int, nargs="*", default=[])
    s = sub.add_parser("summarize")
    s.add_argument("--out-dir", default=str(OUT_DIR))
    return ap


def main() -> None:
    args = build_parser().parse_args()
    if args.cmd == "summarize":
        from findoptimal_sweep_summary import do_summarize

        do_summarize(Path(args.out_dir).resolve())
    elif args.cmd == "pilot":
        do_pilot(args)
    else:
        do_run(args, "A" if args.cmd == "run-a" else "B")


if __name__ == "__main__":
    main()
