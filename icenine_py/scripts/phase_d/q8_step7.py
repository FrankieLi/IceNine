#!/usr/bin/env python3
"""Q8 step 7 (D-20261014-I8): are the noisy images easier to reconstruct?

Analysis only, no full searches. The 200 E0 voxels (seed 0 of the Step 4 runs), image sets `clean`
and `realistic_nodrop` (= RealismConfig "all_nodrop": the neighbour/twin layer plus jitter, hot
pixels and blobs, no random spot drop), the Step 4 arms at min_sin_eta 0: lenient rule (`off`),
C++ rule (`on`) and ray rule with the 2-pixel window (`rayR2`).

  (a) `level0`: per (voxel, images, arm) the level-0 candidate list of reconstruct_voxel (the code
      path of the Step 4 runs, stopped after level 0; deterministic, no random draw): number of
      candidates (after the screen peak_overlap > 0 and the spacing filter), number passing the
      raw screen, symmetry-reduced misorientation of each candidate to the truth, the closest one,
      how many lie within 1 and 3 degrees, and the closest one's rank in the list (sorted by the
      level-0 cost 1 - confidence; stable sort, ties keep the order of the fundamental-zone
      points). Cross-checked against the n_l0, min_dist_l0, n_within the Step 4 runs recorded.
  (b) `width`: the cost at the truth rotated by r in RADII degrees about the same N_AXES random
      axes (seeded by the voxel only, so identical on both image sets). Two readouts per arm: the
      level-0 cost (1 - confidence of the arm's level-0 local cost function: lenient, C++ rule, or
      the ray rule) and the final cost (1 - quality of the arm's final local cost function:
      lenient, C++ rule, C++ rule for the ray arm, which never uses the ray rule at the end).
      Also the fraction of the rotated truths that pass the level-0 screen (peak_overlap > 0 of
      the level-0 global cost function: pixel radius 3, the level-0 Q-max), which is what
      decides how many level-0 candidates a voxel has.
  (c) `summary`: (a) and (b) overall and per group of the Step 4 pairing of no drop against clean
      (fixed / broken / both right / both wrong; wrong = error > 1 degree), for the lenient and
      the ray 2-pixel arm (and the C++ rule arm).

  uv run python scripts/phase_d/q8_step7.py pilot --workers 2 --limit 4
  uv run python scripts/phase_d/q8_step7.py level0 --workers 10
  uv run python scripts/phase_d/q8_step7.py width --workers 10
  uv run python scripts/phase_d/q8_step7.py summary

Writes benchmarks/q8_nodrop/step7/{summary.json, tables.md, per_voxel.npz}. Not a timing run.
"""

import argparse
import json
import sys
import time
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import numpy as np

HERE = Path(__file__).resolve().parent
ICE = HERE.parents[1]
sys.path.insert(0, str(ICE / "scripts" / "findoptimal_robustness"))
sys.path.insert(0, str(ICE / "scripts" / "common"))
import common as C  # noqa: E402
import doc_tables  # noqa: E402
import stats as S  # noqa: E402

CACHE = HERE / "cache" / "q8_step7"
STEP4 = HERE / "cache"  # <arm>_ms0_q8/v<voxel>_<variant>.npz of the Step 4 runs
OUT = ICE / "benchmarks" / "q8_nodrop" / "step7"
E0_PER_VOXEL = ICE / "benchmarks" / "q8_nodrop" / "e0_per_voxel.npz"

VARIANTS = ("clean", "realistic_nodrop")
VLABEL = {"clean": "clean", "realistic_nodrop": "no drop"}
ARMS = ("off", "on", "rayR2")
ALABEL = {"off": "lenient rule", "on": "C++ rule", "rayR2": "ray rule, 2-pixel window"}
MSIN = 0.0
RADII = (0.03, 0.1, 0.3, 0.5, 1.0, 2.0, 3.0)  # degrees
N_AXES = 50
AXIS_SEED = 7_340_011  # axes are drawn from [AXIS_SEED, voxel]: the same on every image set
GROUPS = ("fixed", "broken", "both right", "both wrong")
HIT_DEG = (1.0, 3.0)


def voxels() -> List[int]:
    info = dict(np.load(C.OUT_DIR / "voxels.npz"))
    return [int(v) for v in info["voxel_indices"]][: C.N_VOXELS]


def axes_of(vidx: int) -> np.ndarray:
    """N_AXES random unit vectors, deterministic in the voxel (not in the image set)."""
    ax = np.random.default_rng([AXIS_SEED, vidx]).normal(size=(N_AXES, 3))
    return ax / np.linalg.norm(ax, axis=1, keepdims=True)


def rotated_truth(R_true: np.ndarray, vidx: int) -> np.ndarray:
    """(1 + len(RADII) * N_AXES, 3, 3): the truth, then R(r, axis) @ truth with r-major order."""
    from scipy.spatial.transform import Rotation

    ax = axes_of(vidx)
    rv = np.radians(np.asarray(RADII))[:, None, None] * ax[None]
    rots = Rotation.from_rotvec(rv.reshape(-1, 3)).as_matrix()
    return np.concatenate([R_true[None], rots @ R_true])


class StopL0(Exception):
    pass


def _configure(rec: Any, arm: str) -> Any:
    """Set the arm on the reconstructor (as e0_cpprule.task: centroid screen, C++ rule for on and
    rayR2) and return its ray rule or None. The screen is left at the library default, the
    batched exact screen (the Step 4 runs forced the per-orientation screen, which is about 14
    times slower); the ray arm falls back to the per-orientation screen by itself. The records of
    the Step 4 runs (n_l0, n_within, min_dist_l0) are compared in `summary` (step4_check)."""
    from icenine.ray_overlap import RayRule

    rec.min_sin_eta = MSIN
    sp = rec.setup.search_params
    sp.cpp_overlap_rule = arm in ("on", "rayR2")
    sp.cpp_center_screen = True
    ray = None
    if arm == "rayR2":
        reg = float(rec.setup.exp_setup.config_file.detector_dist_deviation)
        ray = RayRule(arm="R", radius_px=2, registration_mm=reg)
        sp.ray_rule = ray
    return ray


def _restore(rec: Any) -> None:
    sp = rec.setup.search_params
    sp.cpp_overlap_rule = False
    sp.ray_rule = None
    sp.cpp_center_screen = True
    rec.min_sin_eta = 0.3


def _level0_screen_fn(rec: Any, ray: Any) -> Any:
    """The level-0 global cost function of reconstruct_voxel (pixel radius 3, the level-0 Q-max,
    the arm's rule), built with the same arguments."""
    from icenine.cost_functions import VoxelCostFunction

    setup = rec.setup
    return VoxelCostFunction(
        simulator=setup.simulator,
        detector_list=setup.detector_list,
        range_map=setup.range_map,
        exp_data=setup.exp_data,
        **rec._cost_flags(ray),
        sample=setup.sample,
        structure_list=setup.structure_list,
        mode="hard",
        eta_limit=setup.exp_setup.get_eta_limit(),
        pixel_radius=rec.global_pixel_radius,
        max_q=5.0 + rec.params.min_local_resolution + rec.n_q_start_offset,
        min_sin_eta=rec.min_sin_eta,
        ray_rule=ray,
    )


def level0_task(item: Tuple[Any, ...]) -> str:
    """Level 0 of reconstruct_voxel for one (voxel, images, arm), stopped after the discrete
    search; saves the candidates' distances to the truth and costs."""
    from icenine import reconstructor as RC
    from optimizer_sweep import voxel_context

    vidx, vpos, variant, arm, path = item
    W = C.get_worker()
    t0 = time.time()
    keys = np.load(C.CACHE_DIR / "images" / f"v{vidx}_{variant}.npz")["keys"]
    C.attach(keys)
    vctx = voxel_context(W.ctx, vidx)
    rec = W.rec
    Rt = np.asarray(vctx.R_true, dtype=np.float64)
    _configure(rec, arm)
    got: Dict[str, Any] = {}
    passed = {"n": 0, "n_fz": 0}

    def screen_recorder(name: str, info: Dict[str, Any]) -> None:
        if name == "fz":
            passed["n_fz"] += 1
            passed["n"] += len(info["j"])

    orig = RC.run_discrete_search_spaced

    def wrapped(**kw: Any) -> Any:
        if not got:
            kw["recorder"] = screen_recorder
        out = orig(**kw)
        if not got:
            got["R"] = np.stack([c.orientation for c in out]).astype(np.float64)
            got["cost"] = np.array([c.cost for c in out], dtype=np.float64)
            raise StopL0
        return out

    RC.run_discrete_search_spaced = wrapped
    check = -1
    try:
        with C.quiet():
            rec.reconstruct_voxel(vctx.vertices, vctx.voxel.phase, rng=C.run_seed(vpos, 0))
    except StopL0:
        # _level0_screen_fn is the screen of the run: every returned candidate passes it
        fn_scr = _level0_screen_fn(rec, rec.setup.search_params.ray_rule)
        pick = np.linspace(0, len(got["R"]) - 1, min(len(got["R"]), 12)).astype(int)
        check = sum(
            fn_scr.evaluate(
                got["R"][i].astype(np.float32), vctx.vertices, vctx.voxel.phase
            ).peak_overlap
            > 0
            for i in pick
        ) - len(
            pick
        )  # 0 when all pass
    finally:
        RC.run_discrete_search_spaced = orig
        _restore(rec)
    R, cost = got["R"], got["cost"]
    d = C.err_deg(R, np.broadcast_to(Rt, R.shape).copy())
    np.savez_compressed(
        path, d=d, cost=cost, n_screen_pass=passed["n"], n_fz_pass=passed["n_fz"],
        n_keep=max(1, int(len(cost) * rec.keep_fraction)), R_true=Rt, screen_check=check,
        seconds=time.time() - t0,
    )  # fmt: skip
    return f"voxel {vidx} {variant} {arm} n_l0={len(cost)} {time.time() - t0:.1f}s"


def width_task(item: Tuple[Any, ...]) -> str:
    """The costs at the truth and at the rotated truths for one (voxel, images), all arms."""
    from optimizer_sweep import voxel_context

    vidx, vpos, variant, path = item
    W = C.get_worker()
    t0 = time.time()
    keys = np.load(C.CACHE_DIR / "images" / f"v{vidx}_{variant}.npz")["keys"]
    C.attach(keys)
    vctx = voxel_context(W.ctx, vidx)
    rec = W.rec
    Rs = rotated_truth(np.asarray(vctx.R_true, dtype=np.float64), vidx).astype(np.float32)
    out: Dict[str, Any] = {}
    try:
        for arm in ARMS:
            ray = _configure(rec, arm)
            fn_final = rec._make_local_cost_fn()
            fn_l0 = rec._make_local_cost_fn(ray) if ray is not None else fn_final
            fn_scr = _level0_screen_fn(rec, ray)
            c0 = np.empty(len(Rs))
            cf = np.empty(len(Rs))
            scr = np.empty(len(Rs), dtype=np.int32)
            for i, R in enumerate(Rs):
                info = fn_l0.evaluate(R, vctx.vertices, vctx.voxel.phase)
                c0[i] = 1.0 - info.confidence if info.peak_on_detector > 0 else 1.0
                if fn_l0 is fn_final:
                    cf[i] = info.cost
                else:
                    cf[i] = fn_final.evaluate(R, vctx.vertices, vctx.voxel.phase).cost
                scr[i] = fn_scr.evaluate(R, vctx.vertices, vctx.voxel.phase).peak_overlap
            out[f"l0|{arm}"] = c0
            out[f"fin|{arm}"] = cf
            out[f"scr|{arm}"] = scr
            _restore(rec)
    finally:
        _restore(rec)
    np.savez_compressed(path, seconds=time.time() - t0, **out)
    return f"voxel {vidx} {variant} width {time.time() - t0:.1f}s"


def l0_items(limit: int = 0) -> List[Tuple[Any, ...]]:
    items = [
        (v, vpos, var, arm, str(CACHE / "level0" / f"v{v}_{var}_{arm}.npz"))
        for vpos, v in enumerate(voxels())
        for var in VARIANTS
        for arm in ARMS
    ]
    return items[:limit] if limit else items


def width_items(limit: int = 0) -> List[Tuple[Any, ...]]:
    items = [
        (v, vpos, var, str(CACHE / "width" / f"v{v}_{var}.npz"))
        for vpos, v in enumerate(voxels())
        for var in VARIANTS
    ]
    return items[:limit] if limit else items


def run(kind: str, workers: int, limit: int = 0, pilot: bool = False) -> None:
    (CACHE / "level0").mkdir(parents=True, exist_ok=True)
    (CACHE / "width").mkdir(parents=True, exist_ok=True)
    fn, items = (
        (level0_task, l0_items(limit)) if kind == "level0" else (width_task, width_items(limit))
    )
    todo = [it for it in items if pilot or not Path(it[-1]).exists()]
    print(f"{kind}: {len(items)} tasks, {len(todo)} to run, {workers} workers", flush=True)
    t0 = time.time()
    C.run_pool(fn, todo, workers, f"q8_step7_{kind}")
    print(f"{kind} wall {time.time() - t0:.0f}s (contended multi-worker, not a timing claim)")


# ----------------------------------------------------------------------------------------------
# summary
# ----------------------------------------------------------------------------------------------


def load_level0(vox: List[int]) -> Dict[Tuple[str, str], Dict[str, np.ndarray]]:
    """(images, arm) -> per-voxel arrays (one entry per voxel in vox)."""
    res: Dict[Tuple[str, str], Dict[str, np.ndarray]] = {}
    for var in VARIANTS:
        for arm in ARMS:
            rows = []
            for v in vox:
                p = CACHE / "level0" / f"v{v}_{var}_{arm}.npz"
                if not p.exists():
                    break
                z = np.load(p)
                d, cost = z["d"], z["cost"]
                j = int(np.argmin(d))
                rows.append(
                    dict(
                        n=len(d), n_screen=int(z["n_screen_pass"]), n_keep=int(z["n_keep"]),
                        dmin=float(d[j]), rank=j + 1, rank_ties=int((cost <= cost[j]).sum()),
                        rank_frac=(j + 1) / len(d),
                        in_keep=bool(j < int(z["n_keep"])), cost_closest=float(cost[j]),
                        cost_best=float(cost[0]), d_best=float(d[0]),
                        n1=int((d < HIT_DEG[0]).sum()), n3=int((d < HIT_DEG[1]).sum()),
                    )
                )  # fmt: skip
            else:
                res[(var, arm)] = {k: np.array([r[k] for r in rows]) for k in rows[0]}
    return res


def load_width(vox: List[int]) -> Dict[str, Dict[str, np.ndarray]]:
    """images -> {'l0|arm' / 'fin|arm' / 'scr|arm': (n_vox, 1 + n_r, N_AXES)} with index 0 the
    truth (its value is repeated along the axis dimension). scr = peak_overlap of the level-0
    screen function (a rotated truth passes the screen when it is > 0)."""
    n_r = len(RADII)
    res: Dict[str, Dict[str, np.ndarray]] = {}
    for var in VARIANTS:
        zs = []
        for v in vox:
            p = CACHE / "width" / f"v{v}_{var}.npz"
            if not p.exists():
                break
            zs.append(np.load(p))
        else:
            d: Dict[str, np.ndarray] = {}
            for arm in ARMS:
                for kind in ("l0", "fin", "scr"):
                    a = np.stack([z[f"{kind}|{arm}"] for z in zs]).astype(float)
                    truth = np.repeat(a[:, :1], N_AXES, axis=1)[:, None, :]
                    d[f"{kind}|{arm}"] = np.concatenate(
                        [truth, a[:, 1:].reshape(len(vox), n_r, N_AXES)], axis=1
                    )
            res[var] = d
    return res


def width_of(curve: np.ndarray, radii: np.ndarray) -> Tuple[float, float]:
    """(grid width, interpolated width) of a median cost curve over [0] + RADII: the smallest r
    at which the curve reaches halfway from the truth value (index 0) to the value at the last
    radius, and the same by linear interpolation in log r between the bracketing grid points
    (nan if the curve never reaches it, or the 3 degree value is not above the truth)."""
    c0, c3 = curve[0], curve[-1]
    if not c3 > c0:
        return float("nan"), float("nan")
    half = c0 + 0.5 * (c3 - c0)
    for i in range(1, len(curve)):
        if curve[i] >= half:
            lo_r = radii[i - 1] if i > 1 else radii[0] / 3.0  # r = 0 has no log; start below
            lo_c = curve[i - 1]
            t = (half - lo_c) / (curve[i] - lo_c) if curve[i] > lo_c else 1.0
            interp = float(np.exp(np.log(lo_r) + t * (np.log(radii[i]) - np.log(lo_r))))
            return float(radii[i]), interp
    return float("nan"), float("nan")


def pairing(vox: List[int]) -> Dict[str, Dict[str, np.ndarray]]:
    """arm -> group name -> boolean mask over vox (Step 4 pairing, no drop against clean)."""
    z = np.load(E0_PER_VOXEL)
    order = {int(v): i for i, v in enumerate(z["voxel"])}
    idx = np.array([order[v] for v in vox])
    out: Dict[str, Dict[str, np.ndarray]] = {}
    for arm in ARMS:
        key = f"{arm}_ms0"
        wc = z[f"{key}|clean|err"][idx] > C.WRONG_DEG
        wn = z[f"{key}|realistic_nodrop|err"][idx] > C.WRONG_DEG
        out[arm] = {
            "fixed": wc & ~wn,
            "broken": ~wc & wn,
            "both right": ~wc & ~wn,
            "both wrong": wc & wn,
        }
    return out


def step4_check(vox: List[int], l0: Dict[Tuple[str, str], Dict[str, np.ndarray]]) -> Dict[str, Any]:
    """Level-0 records of the Step 4 runs against this script's: n_l0 and the 3-degree count must
    be equal (deterministic level 0); min_dist_l0 equal to float32 noise."""
    out: Dict[str, Any] = {}
    for (var, arm), d in l0.items():
        root = STEP4 / f"{arm}_ms0_q8"
        files = [root / f"v{v}_{var}.npz" for v in vox]
        if not all(f.exists() for f in files):
            continue
        z = [np.load(f) for f in files]
        n_old = np.array([int(x["n_l0"]) for x in z])
        w_old = np.array([int(x["n_within"]) for x in z])
        m_old = np.array([float(x["min_dist_l0"]) for x in z])
        out[f"{var}|{arm}"] = dict(
            n_l0_equal=int((n_old == d["n"]).sum()), n_within3_equal=int((w_old == d["n3"]).sum()),
            min_dist_max_abs_diff_deg=float(np.abs(m_old - d["dmin"]).max()), n=len(vox),
        )  # fmt: skip
    return out


# table header of a row key (the rows keep their short keys)
HEADER = {
    "arm": "arm", "images": "images", "group": "group", "n_vox": "voxels",
    "n_l0": "median candidates", "n_screen": "median passing the screen",
    "dmin": "median closest (deg)", "hit1": "closest within 1 deg", "hit3": "closest within 3 deg",
    "n1": "median within 1 deg", "n3": "median within 3 deg", "rank": "median rank of closest",
    "top1": "closest ranked first", "in_keep": "closest in top quarter",
    "more_fewer_equal": "candidates: more / fewer / equal (no drop)",
    "ratio": "median ratio of candidates",
    "closer_farther_equal": "closest: nearer / farther / same (no drop)",
    "h1": "within 1 deg: both / only no drop / only clean / neither",
    "h3": "within 3 deg: both / only no drop / only clean / neither",
    "keep": "in top quarter: both / only no drop / only clean / neither",
    "keep_str": "in top quarter: both / only no drop / only clean / neither",
    "hit3_str": "within 3 deg: both / only no drop / only clean / neither",
    "c_n_l0": "median candidates", "c_dmin": "median closest (deg)",
    "c_hit1": "closest within 1 deg", "c_hit3": "closest within 3 deg",
    "c_in_keep": "closest in top quarter", "c_rank": "median rank of closest",
    "c_rank_frac": "median rank / candidates", "wi_l0": "level-0 width (deg)",
    "p1": "screen pass at 1 deg", "p3": "screen pass at 3 deg",
    "width": "width, grid (deg)", "interp": "width, interpolated (deg)",
    "per_voxel": "median per-voxel width, grid (deg)", "cost": "cost", "n": "voxels",
    "wider": "no drop wider", "narrower": "no drop narrower", "tied": "same",
    "median_clean": "median per-voxel width, clean (deg)",
    "median_nodrop": "median per-voxel width, no drop (deg)",
    "pooled_clean": "pooled width, clean (deg)", "pooled_nodrop": "pooled width, no drop (deg)",
    "wratio": "pooled width ratio (no drop / clean)",
    "ratio_ci": "ratio 95% interval (voxel bootstrap)",
    "width_clean": "level-0 width, clean (deg)", "width_nodrop": "level-0 width, no drop (deg)",
    "width_ratio": "width ratio", "width_ratio_ci": "ratio 95% interval (voxel bootstrap)",
    "pass3_clean": "screen pass at 3 deg, clean", "pass3_nodrop": "screen pass at 3 deg, no drop",
}  # fmt: skip


def table(rows: Any, columns: List[str], formats: Optional[Dict[str, str]] = None) -> str:
    """doc_tables.markdown_table with readable column headers (HEADER; unknown keys unchanged)."""
    h = lambda k: HEADER.get(k, k)  # noqa: E731
    rows2 = [None if r is None else {h(k): v for k, v in r.items()} for r in rows]
    fm = {h(k): v for k, v in (formats or {}).items()}
    return doc_tables.markdown_table(rows2, [h(c) for c in columns], fm)


def med(x: np.ndarray) -> float:
    return float(np.median(x)) if len(x) else float("nan")


def count_ci(k: int, n: int) -> str:
    lo, hi = S.wilson(k, n)
    return f"{k} of {n} ({100 * lo:.0f} to {100 * hi:.0f}%)"


def l0_stats(d: Dict[str, np.ndarray], m: Optional[np.ndarray] = None) -> Dict[str, Any]:
    m = np.ones(len(d["n"]), bool) if m is None else m
    n = int(m.sum())
    if n == 0:
        return dict(n=0)
    return dict(
        n=n, n_l0=med(d["n"][m]), n_screen=med(d["n_screen"][m]), dmin=med(d["dmin"][m]),
        n1=med(d["n1"][m]), n3=med(d["n3"][m]), rank=med(d["rank"][m]),
        rank_frac=med(d["rank_frac"][m]),
        hit1=int((d["dmin"][m] < HIT_DEG[0]).sum()), hit3=int((d["dmin"][m] < HIT_DEG[1]).sum()),
        in_keep=int(d["in_keep"][m].sum()), top1=int((d["rank"][m] == 1).sum()),
        best_within1=int((d["d_best"][m] < HIT_DEG[0]).sum()),
    )  # fmt: skip


def curve_stats(a: np.ndarray, m: Optional[np.ndarray] = None) -> Dict[str, Any]:
    """a: (n_vox, 1 + n_r, N_AXES). Pooled median over voxels (and axes) of each radius, the
    widths of that curve, and the per-voxel widths (median over axes, then the same rule)."""
    m = np.ones(len(a), bool) if m is None else m
    if not m.any():
        return dict(n=0)
    sel = a[m]
    radii = np.asarray(RADII)
    curve = np.array([np.median(sel[:, i, :]) for i in range(sel.shape[1])])
    g, ip = width_of(curve, radii)
    pv = np.array([width_of(np.median(x, axis=1), radii)[0] for x in sel])
    return dict(
        n=int(m.sum()),
        curve=curve.tolist(),
        mean_curve=[float(sel[:, i, :].mean()) for i in range(sel.shape[1])],
        width=g,
        width_interp=ip,
        width_voxel_median=float(np.nanmedian(pv)) if np.isfinite(pv).any() else float("nan"),
    )


def cs_mean_at(a: np.ndarray, m: np.ndarray, r: float) -> float:
    """Mean cost at radius r (pooled over the voxels of mask m and the axes)."""
    return float(a[m][:, 1 + RADII.index(r), :].mean())


def pass_curve(a: np.ndarray, m: Optional[np.ndarray] = None) -> List[float]:
    """a: screen peak_overlap (n_vox, 1 + n_r, N_AXES). Fraction of the rotated truths of each
    radius (pooled over voxels and axes) that pass the level-0 screen (peak_overlap > 0)."""
    m = np.ones(len(a), bool) if m is None else m
    return [float((a[m][:, i, :] > 0).mean()) for i in range(a.shape[1])]


def per_voxel_width(a: np.ndarray) -> np.ndarray:
    """Interpolated width of each voxel's own curve (median over its axes)."""
    radii = np.asarray(RADII)
    return np.array([width_of(np.median(x, axis=1), radii)[1] for x in a])


def pooled_width(a: np.ndarray) -> float:
    """Interpolated width of the pooled median curve of a (n_vox, 1 + n_r, N_AXES)."""
    curve = np.array([np.median(a[:, i, :]) for i in range(a.shape[1])])
    return width_of(curve, np.asarray(RADII))[1]


def boot_width_diff(
    a_clean: np.ndarray, a_nodrop: np.ndarray, n_boot: int = 1000, seed: int = 0
) -> Dict[str, float]:
    """Pooled interpolated widths of both image sets and their difference and ratio, with a 95%
    percentile interval over voxels (the same voxels resampled for both sets, n_boot draws)."""
    rng = np.random.default_rng(seed)
    n = len(a_clean)
    diffs, ratios = [], []
    for _ in range(n_boot):
        i = rng.integers(0, n, n)
        wc, wn = pooled_width(a_clean[i]), pooled_width(a_nodrop[i])
        diffs.append(wn - wc)
        ratios.append(wn / wc)
    wc, wn = pooled_width(a_clean), pooled_width(a_nodrop)
    d_lo, d_hi = np.nanpercentile(diffs, [2.5, 97.5])
    r_lo, r_hi = np.nanpercentile(ratios, [2.5, 97.5])
    return dict(
        clean=wc, nodrop=wn, diff=wn - wc, diff_lo=float(d_lo), diff_hi=float(d_hi),
        ratio=wn / wc, ratio_lo=float(r_lo), ratio_hi=float(r_hi), n_boot=n_boot,
    )  # fmt: skip


def summary() -> None:
    vox = voxels()
    l0 = load_level0(vox)
    wd = load_width(vox)
    grp = pairing(vox)
    OUT.mkdir(parents=True, exist_ok=True)
    summ: Dict[str, Any] = dict(
        n_voxels=len(vox), radii_deg=list(RADII), n_axes=N_AXES, min_sin_eta=MSIN,
        wrong_deg=C.WRONG_DEG, hit_deg=list(HIT_DEG),
        groups={arm: {g: int(m.sum()) for g, m in grp[arm].items()} for arm in ARMS},
    )  # fmt: skip
    if l0:
        summ["step4_check"] = step4_check(vox, l0)
    tables: Dict[str, str] = {}

    # ---- (a) candidates ----
    if l0:
        rows, summ["candidates"] = [], {}
        for arm in ARMS:
            for var in VARIANTS:
                d = l0.get((var, arm))
                if d is None:
                    continue
                s = l0_stats(d)
                summ["candidates"][f"{arm}|{var}"] = s
                rows.append(
                    dict(
                        arm=ALABEL[arm], images=VLABEL[var], n_l0=s["n_l0"], n_screen=s["n_screen"],
                        dmin=s["dmin"], hit1=count_ci(s["hit1"], s["n"]),
                        hit3=count_ci(s["hit3"], s["n"]), n1=s["n1"], n3=s["n3"], rank=s["rank"],
                        top1=s["top1"], in_keep=s["in_keep"],
                    )
                )  # fmt: skip
        tables["step7_candidates"] = table(
            rows,
            ["arm", "images", "n_l0", "n_screen", "dmin", "hit1", "hit3", "n1", "n3", "rank", "top1",
             "in_keep"],
            {"n_l0": ".0f", "n_screen": ".0f", "dmin": ".3f", "n1": ".0f", "n3": ".0f",
             "rank": ".0f"},
        )  # fmt: skip
        # paired clean -> no drop per voxel
        prow, summ["candidates_paired"] = [], {}
        for arm in ARMS:
            c, n_ = l0.get(("clean", arm)), l0.get(("realistic_nodrop", arm))
            if c is None or n_ is None:
                continue
            rec: Dict[str, Any] = {}
            for name, th in (("1 deg", HIT_DEG[0]), ("3 deg", HIT_DEG[1])):
                hc, hn = c["dmin"] < th, n_["dmin"] < th
                rec[name] = dict(
                    both=int((hc & hn).sum()), only_nodrop=int((~hc & hn).sum()),
                    only_clean=int((hc & ~hn).sum()), neither=int((~hc & ~hn).sum()),
                )  # fmt: skip
            kc, kn = c["in_keep"], n_["in_keep"]
            rec["keep"] = dict(
                both=int((kc & kn).sum()), only_nodrop=int((~kc & kn).sum()),
                only_clean=int((kc & ~kn).sum()), neither=int((~kc & ~kn).sum()),
            )  # fmt: skip
            rec["more_candidates_nodrop"] = int((n_["n"] > c["n"]).sum())
            rec["fewer_candidates_nodrop"] = int((n_["n"] < c["n"]).sum())
            rec["equal_candidates"] = int((n_["n"] == c["n"]).sum())
            rec["closer_nodrop"] = int((n_["dmin"] < c["dmin"]).sum())
            rec["farther_nodrop"] = int((n_["dmin"] > c["dmin"]).sum())
            rec["equal_dmin"] = int((n_["dmin"] == c["dmin"]).sum())
            rec["median_ratio_n_l0"] = med(n_["n"] / np.maximum(c["n"], 1))
            summ["candidates_paired"][arm] = rec
            prow.append(
                dict(
                    arm=ALABEL[arm],
                    more_fewer_equal=f"{rec['more_candidates_nodrop']} / {rec['fewer_candidates_nodrop']} / {rec['equal_candidates']}",
                    ratio=rec["median_ratio_n_l0"],
                    closer_farther_equal=f"{rec['closer_nodrop']} / {rec['farther_nodrop']} / {rec['equal_dmin']}",
                    h1=f"{rec['1 deg']['both']} / {rec['1 deg']['only_nodrop']} / {rec['1 deg']['only_clean']} / {rec['1 deg']['neither']}",
                    h3=f"{rec['3 deg']['both']} / {rec['3 deg']['only_nodrop']} / {rec['3 deg']['only_clean']} / {rec['3 deg']['neither']}",
                    keep=f"{rec['keep']['both']} / {rec['keep']['only_nodrop']} / {rec['keep']['only_clean']} / {rec['keep']['neither']}",
                )
            )  # fmt: skip
        tables["step7_candidates_paired"] = table(
            prow,
            ["arm", "more_fewer_equal", "ratio", "closer_farther_equal", "h1", "h3", "keep"],
            {"ratio": ".2f"},
        )

    # ---- (b) cost width ----
    if len(wd) == len(VARIANTS):
        wrows, summ["width"] = [], {}
        for arm in ARMS:
            for kind, klabel in (("l0", "level-0 cost"), ("fin", "final cost")):
                for var in VARIANTS:
                    s = curve_stats(wd[var][f"{kind}|{arm}"])
                    summ["width"][f"{arm}|{kind}|{var}"] = s
                    wrows.append(
                        dict(
                            arm=ALABEL[arm], cost=klabel, images=VLABEL[var],
                            **{f"r{i}": s["curve"][i] for i in range(len(s["curve"]))},
                            width=s["width"], interp=s["width_interp"], per_voxel=s["width_voxel_median"],
                        )
                    )  # fmt: skip
        cols = ["arm", "cost", "images"] + [f"r{i}" for i in range(len(RADII) + 1)]
        labels = ["truth"] + [f"{r:g}" for r in RADII]
        tbl = table(
            wrows, cols + ["width", "interp", "per_voxel"],
            {**{f"r{i}": ".4f" for i in range(len(RADII) + 1)},
             "width": "g", "interp": ".3f", "per_voxel": "g"},
        )  # fmt: skip
        head = tbl.split("\n", 1)[0]
        for i, lab in enumerate(["truth"] + [f"{r:g} deg" for r in RADII]):
            head = head.replace(f"r{i}", lab, 1)
        tables["step7_width"] = head + "\n" + tbl.split("\n", 1)[1]
        # the rise of the median cost above its value at the truth, no drop as a fraction of clean
        # (below 1: the cost rises more slowly on the no-drop images)
        rrows, summ["rise_ratio"] = [], {}
        for arm in ARMS:
            for kind, klabel in (("l0", "level-0 cost"), ("fin", "final cost")):
                cc = np.array(summ["width"][f"{arm}|{kind}|clean"]["curve"])
                cn = np.array(summ["width"][f"{arm}|{kind}|realistic_nodrop"]["curve"])
                rr = ((cn - cn[0]) / (cc - cc[0])).tolist()
                summ["rise_ratio"][f"{arm}|{kind}"] = rr
                rrows.append(
                    dict(
                        arm=ALABEL[arm], cost=klabel, **{f"r{i}": x for i, x in enumerate(rr) if i}
                    )
                )
        t = table(
            rrows,
            ["arm", "cost"] + [f"r{i}" for i in range(1, len(RADII) + 1)],
            {f"r{i}": ".3f" for i in range(1, len(RADII) + 1)},
        )
        h = t.split("\n", 1)[0]
        for i, lab in enumerate([f"{r:g} deg" for r in RADII], 1):
            h = h.replace(f"r{i}", lab, 1)
        tables["step7_rise_ratio"] = h + "\n" + t.split("\n", 1)[1]
        # the mean level-0 cost (the median saturates at 1 from 2 degrees) and the screen pass rate
        mrows, srows, summ["screen_pass"] = [], [], {}
        for arm in ARMS:
            for var in VARIANTS:
                mc = summ["width"][f"{arm}|l0|{var}"]["mean_curve"]
                mrows.append(
                    dict(
                        arm=ALABEL[arm],
                        images=VLABEL[var],
                        **{f"r{i}": x for i, x in enumerate(mc)},
                    )
                )
                pc = pass_curve(wd[var][f"scr|{arm}"])
                summ["screen_pass"][f"{arm}|{var}"] = pc
                srows.append(
                    dict(
                        arm=ALABEL[arm],
                        images=VLABEL[var],
                        **{f"r{i}": x for i, x in enumerate(pc)},
                    )
                )
        for key, rows_, fmt in (
            ("step7_width_mean", mrows, ".4f"),
            ("step7_screen_pass", srows, ".3f"),
        ):
            t = table(
                rows_,
                ["arm", "images"] + [f"r{i}" for i in range(len(RADII) + 1)],
                {f"r{i}": fmt for i in range(len(RADII) + 1)},
            )
            h = t.split("\n", 1)[0]
            for i, lab in enumerate(["truth"] + [f"{r:g} deg" for r in RADII]):
                h = h.replace(f"r{i}", lab, 1)
            tables[key] = h + "\n" + t.split("\n", 1)[1]
        # paired per-voxel width, no drop against clean
        prow2, summ["width_paired"] = [], {}
        for arm in ARMS:
            for kind, klabel in (("l0", "level-0 cost"), ("fin", "final cost")):
                ac, an = wd["clean"][f"{kind}|{arm}"], wd["realistic_nodrop"][f"{kind}|{arm}"]
                wc, wn = per_voxel_width(ac), per_voxel_width(an)
                ok = np.isfinite(wc) & np.isfinite(wn)
                bd = boot_width_diff(ac, an)
                rec2 = dict(
                    n=int(ok.sum()), wider=int((wn[ok] > wc[ok]).sum()),
                    narrower=int((wn[ok] < wc[ok]).sum()), tied=int((wn[ok] == wc[ok]).sum()),
                    median_clean=med(wc[ok]), median_nodrop=med(wn[ok]), pooled=bd,
                )  # fmt: skip
                summ["width_paired"][f"{arm}|{kind}"] = rec2
                prow2.append(
                    dict(
                        arm=ALABEL[arm], cost=klabel, n=rec2["n"], wider=rec2["wider"],
                        narrower=rec2["narrower"], tied=rec2["tied"],
                        median_clean=rec2["median_clean"], median_nodrop=rec2["median_nodrop"],
                        pooled_clean=bd["clean"], pooled_nodrop=bd["nodrop"],
                        wratio=bd["ratio"], ratio_ci=f"{bd['ratio_lo']:.2f} to {bd['ratio_hi']:.2f}",
                    )
                )  # fmt: skip
        tables["step7_width_paired"] = table(
            prow2,
            ["arm", "cost", "n", "wider", "narrower", "tied", "median_clean", "median_nodrop",
             "pooled_clean", "pooled_nodrop", "wratio", "ratio_ci"],
            {"median_clean": ".3f", "median_nodrop": ".3f", "pooled_clean": ".3f",
             "pooled_nodrop": ".3f", "wratio": ".3f"},
        )  # fmt: skip

    # ---- (c) per group of the Step 4 pairing ----
    summ["by_group"] = {}
    for arm in ARMS:
        grows = []
        for g in GROUPS:
            m = grp[arm][g]
            for var in VARIANTS:
                rec3: Dict[str, Any] = dict(group=g, images=VLABEL[var], n_vox=int(m.sum()))
                if (var, arm) in l0:
                    s = l0_stats(l0[(var, arm)], m)
                    rec3.update({f"c_{k}": v for k, v in s.items() if k != "n"})
                if var in wd:
                    for kind in ("l0", "fin"):
                        cs = curve_stats(wd[var][f"{kind}|{arm}"], m)
                        rec3[f"w_{kind}"] = cs.get("width", float("nan"))
                        rec3[f"wi_{kind}"] = cs.get("width_interp", float("nan"))
                        rec3[f"ratio3_{kind}"] = (
                            cs["curve"][-1] / cs["curve"][0] if cs.get("curve") and cs["curve"][0] > 0 else float("nan")
                        )  # fmt: skip
                    pc = pass_curve(wd[var][f"scr|{arm}"], m)
                    rec3["p1"] = pc[1 + RADII.index(1.0)]
                    rec3["p3"] = pc[1 + RADII.index(3.0)]
                    rec3["mean_l0_1deg"] = cs_mean_at(wd[var][f"l0|{arm}"], m, 1.0)
                summ["by_group"].setdefault(arm, {}).setdefault(g, {})[var] = rec3
                grows.append(rec3)
        tables[f"step7_groups_{arm}"] = table(
            grows,
            ["group", "images", "n_vox", "c_n_l0", "c_dmin", "c_hit1", "c_hit3", "c_in_keep",
             "c_rank", "c_rank_frac", "wi_l0", "p1", "p3"],
            {"c_n_l0": ".0f", "c_dmin": ".3f", "c_rank": ".0f", "c_rank_frac": ".3f",
             "wi_l0": ".3f", "p1": ".3f", "p3": ".3f"},
        )  # fmt: skip
    # paired within each group: does the closest candidate survive the top-quarter cut / lie within
    # 3 degrees on the clean and the no-drop images; level-0 width ratio with a voxel bootstrap
    if l0 and len(wd) == len(VARIANTS):
        krows, summ["groups_paired"] = [], {}
        for arm in ARMS:
            c, n_ = l0[("clean", arm)], l0[("realistic_nodrop", arm)]
            for g in GROUPS:
                m = grp[arm][g]
                rec4: Dict[str, Any] = dict(arm=ALABEL[arm], group=g, n_vox=int(m.sum()))
                for name, hc, hn in (
                    ("keep", c["in_keep"], n_["in_keep"]),
                    ("hit3", c["dmin"] < HIT_DEG[1], n_["dmin"] < HIT_DEG[1]),
                ):
                    hc, hn = hc[m], hn[m]
                    rec4[name] = dict(
                        both=int((hc & hn).sum()), only_nodrop=int((~hc & hn).sum()),
                        only_clean=int((hc & ~hn).sum()), neither=int((~hc & ~hn).sum()),
                    )  # fmt: skip
                    rec4[f"{name}_str"] = (
                        f"{rec4[name]['both']} / {rec4[name]['only_nodrop']} / "
                        f"{rec4[name]['only_clean']} / {rec4[name]['neither']}"
                    )
                ac, an = wd["clean"][f"l0|{arm}"], wd["realistic_nodrop"][f"l0|{arm}"]
                bd = boot_width_diff(ac[m], an[m], n_boot=1000, seed=1)
                rec4["width_ratio"] = bd["ratio"]
                rec4["width_ratio_ci"] = f"{bd['ratio_lo']:.2f} to {bd['ratio_hi']:.2f}"
                rec4["width_clean"], rec4["width_nodrop"] = bd["clean"], bd["nodrop"]
                pc = pass_curve(wd["clean"][f"scr|{arm}"], m)
                pn = pass_curve(wd["realistic_nodrop"][f"scr|{arm}"], m)
                rec4["pass3_clean"] = pc[1 + RADII.index(3.0)]
                rec4["pass3_nodrop"] = pn[1 + RADII.index(3.0)]
                summ["groups_paired"].setdefault(arm, {})[g] = rec4
                krows.append(rec4)
        tables["step7_groups_paired"] = table(
            krows,
            ["arm", "group", "n_vox", "keep_str", "hit3_str", "width_clean", "width_nodrop",
             "width_ratio", "width_ratio_ci", "pass3_clean", "pass3_nodrop"],
            {"width_clean": ".3f", "width_nodrop": ".3f", "width_ratio": ".2f",
             "pass3_clean": ".3f", "pass3_nodrop": ".3f"},
        )  # fmt: skip
    (OUT / "summary.json").write_text(json.dumps(summ, indent=1, default=float) + "\n")
    doc_tables.write_tables(OUT / "tables.md", tables)
    np.savez_compressed(
        OUT / "per_voxel.npz", voxel=np.array(vox),
        **{f"l0|{arm}|{var}|{k}": v for (var, arm), d in l0.items() for k, v in d.items()},
        **{f"wrong_group|{arm}|{g}": m for arm in ARMS for g, m in grp[arm].items()},
    )  # fmt: skip
    print(json.dumps({k: summ[k] for k in ("n_voxels", "groups")}, indent=1))
    print("wrote", OUT / "summary.json", OUT / "tables.md")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    ap.add_argument("cmd", choices=["pilot", "level0", "width", "summary"])
    ap.add_argument("--workers", type=int, default=10)
    ap.add_argument("--limit", type=int, default=0)
    a = ap.parse_args()
    if a.cmd == "summary":
        summary()
    elif a.cmd == "pilot":
        t0 = time.time()
        run("level0", a.workers, a.limit or 6, pilot=True)
        run("width", a.workers, a.limit or 4, pilot=True)
        print(f"pilot wall {time.time() - t0:.0f}s")
    else:
        run(a.cmd, a.workers, a.limit)


if __name__ == "__main__":
    main()
