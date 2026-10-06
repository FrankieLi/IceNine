#!/usr/bin/env python3
"""
Existing optimizers on the cases of the single-voxel perturbation sweep.

Runs the project's non-learned optimizers on exactly the cases scripts/perturbation_sweep.py
gave the network (same 50 voxels, 10 radii, 20 directions per radius, same true orientation,
same perturbed nominal, same distractor draws and realism seed), so their error versus the
perturbation radius r can be put next to the network's.

Methods (one result per case, variant and pass; error = angle of R_est R_true^T, the sweep's
metric: no crystal-symmetry reduction, irrelevant for r <= 5 deg):

  mc     MCOptimizer on the binary pixel-overlap cost, 3500 steps, 2 restarts, step fraction 0.5,
         search box 1.5 r (the protocol of benchmarks/bench_hp_sweep.py; it is told r).
  adam   Riemannian Adam (geoopt) on the differentiable cost: scale 2 (8x max-pool), omega window 1,
         lr 1e-4, 100 steps (the HP-sweep best at 1 deg).
  gn     CentroidGaussNewton (plain least squares) on the network's frame-coded windows.
  huber  The same with a Huber loss, c = 1.
  GN is run one-shot (pass 1) and iterated x3 with re-centring, exactly like the net.

What the image-based optimizers (mc, adam) see. No full-sample image exists, so each case gets
its own detector images, built from pixel sets:
  clean  every spot of the target voxel at its TRUE orientation (all frames, both detectors),
         lit pixels as the simulator's rasteriser lights them (orientation_eval.lit_pixel_set).
  all    the same, plus (i) every spot of the distractor sources (2 nearest neighbours + a Sigma3
         twin, the sweep's draws) over the WHOLE detector, plus (ii) the sweep's realism edits:
         the pixel-level difference between the net's window before and after
         make_realistic_dataset is applied to the image at the window's detector position (pixels
         dropped / grown / hot / blob; a window zeroed as a missing spot removes all its pixels).
So the optimizers see the whole detector (all of the voxel's spots with |g| <= Q_max, not only the
net's ROI windows), but the same random nuisances inside the windows. Both cost functions use the
net's eligibility: |g| <= Q_max, both detectors, |sin eta| >= min_sin_eta.

Usage (from icenine_py/):
  uv run python scripts/optimizer_sweep.py run --workers 10 \
      --out-dir benchmarks/toy_orientation_sweep
  uv run python scripts/optimizer_sweep.py summarize --out-dir benchmarks/toy_orientation_sweep
"""

import argparse
import contextlib
import json
import os
import sys
import time
from pathlib import Path
from types import SimpleNamespace
from typing import Any, Dict, Iterator, List, NamedTuple, Optional, Sequence, Tuple

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")

import numpy as np
import torch

ICENINE_PY = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(Path(__file__).parent))
sys.path.insert(0, str(ICENINE_PY / "benchmarks"))

import perturbation_sweep as ps  # noqa: E402

METHODS = ["mc", "adam", "gn", "huber"]
VARIANTS = ps.VARIANTS
IMAGE_METHODS = ("mc", "adam")
N_PASS = 3
GN_STATUS = {"step_below_tol": 0, "no_descent": 1, "max_iter": 2, "too_few_spots": 3, "error": 4}
FACTOR = 8  # scale 2 of the HP sweep: 8x max-pool
SCALE_INDEX = 2


# ---------------------------------------------------------------------------
# Pixel-set helpers (no physics): keys, grouping, the window -> image edit
# ---------------------------------------------------------------------------


class Geometry(NamedTuple):
    """Detector image geometry. A pixel is one int64 key = img_id * H * W + row * W + col, with
    img_id = frame * n_det + detector (the flat index of the cost functions' image stacks)."""

    n_det: int
    n_omega: int
    H: int
    W: int


def encode_pixels(
    frame: np.ndarray, det: np.ndarray, row: np.ndarray, col: np.ndarray, geo: Geometry
) -> np.ndarray:
    img = np.asarray(frame, dtype=np.int64) * geo.n_det + np.asarray(det, dtype=np.int64)
    return (img * geo.H + np.asarray(row, dtype=np.int64)) * geo.W + np.asarray(col, dtype=np.int64)


def decode_pixels(keys: np.ndarray, geo: Geometry) -> Tuple[np.ndarray, ...]:
    """(frame, det, row, col) of pixel keys."""
    keys = np.asarray(keys, dtype=np.int64)
    col = keys % geo.W
    rest = keys // geo.W
    row = rest % geo.H
    img = rest // geo.H
    return img // geo.n_det, img % geo.n_det, row, col


def window_pixel_keys(
    windows: np.ndarray,
    spec: Any,
    det_idx: np.ndarray,
    frame_half_width: int,
    geo: Geometry,
) -> np.ndarray:
    """Sorted unique pixel keys of the lit pixels of frame-coded windows (n, W, W) uint8, placed
    at their detector positions (spec.col0 / row0) and frames (spec.frame0 + code - 1 - K).
    Pixels falling outside the detector or the frame range are dropped."""
    m, y, x = np.nonzero(windows)
    code = windows[m, y, x].astype(np.int64)
    f = np.asarray(spec.frame0)[m] + code - 1 - frame_half_width
    r = np.asarray(spec.row0)[m] + y
    c = np.asarray(spec.col0)[m] + x
    ok = (f >= 0) & (f < geo.n_omega) & (r >= 0) & (r < geo.H) & (c >= 0) & (c < geo.W)
    return np.unique(encode_pixels(f[ok], np.asarray(det_idx)[m][ok], r[ok], c[ok], geo))


def realism_edit(
    clean_w: np.ndarray,
    dis_w: np.ndarray,
    real_w: np.ndarray,
    spec: Any,
    det_idx: np.ndarray,
    frame_half_width: int,
    geo: Geometry,
) -> Tuple[np.ndarray, np.ndarray]:
    """Image edit equivalent to the realism layer of one case: (removed, added) pixel keys.

    The windows before realism are the target's own windows overlaid on the distractor layer
    (the target's pixels win); `real_w` is the same case after make_realistic_dataset. Removed =
    pixels lit before and lit in none of the windows after (flipped off or a zeroed window);
    added = pixels lit after and in none before (grown, hot, blob). Taking differences of the
    unions over windows keeps overlapping windows consistent."""
    before = np.where(clean_w > 0, clean_w, dis_w)
    kb = window_pixel_keys(before, spec, det_idx, frame_half_width, geo)
    ka = window_pixel_keys(real_w, spec, det_idx, frame_half_width, geo)
    return np.setdiff1d(kb, ka), np.setdiff1d(ka, kb)


def group_pixels(keys: np.ndarray, geo: Geometry) -> Dict[int, np.ndarray]:
    """img_id -> sorted flat in-image pixel positions (row * W + col) of a sorted unique key
    array."""
    keys = np.asarray(keys, dtype=np.int64)
    hw = geo.H * geo.W
    img, pos = keys // hw, keys % hw
    if len(keys) == 0:
        return {}
    cuts = np.flatnonzero(np.diff(img)) + 1
    return {int(i[0]): p for i, p in zip(np.split(img, cuts), np.split(pos, cuts))}


class _PixelImage:
    """Minimal ImageData stand-in for the hard cost: num_rows / num_cols / get_binary_numpy."""

    def __init__(self, positions: Optional[np.ndarray], geo: Geometry, zeros: np.ndarray):
        self.num_rows, self.num_cols = geo.H, geo.W
        self._pos, self._zeros, self._arr = positions, zeros, None

    def get_binary_numpy(self) -> np.ndarray:
        if self._pos is None:
            return self._zeros
        if self._arr is None:  # zero pages are not touched until written: cheap for 4 MB images
            self._arr = np.zeros((self.num_rows, self.num_cols), dtype=np.uint8)
            self._arr.reshape(-1)[self._pos] = 1
        return self._arr


class PixelSetData:
    """ExperimentalData stand-in for VoxelCostFunction (hard): get_image(frame, det) from sparse
    pixel groups; dark images share one zero array."""

    def __init__(self, groups: Dict[int, np.ndarray], geo: Geometry, zeros: np.ndarray):
        self.n_omega_intervals, self.n_detectors = geo.n_omega, geo.n_det
        self._groups, self._geo, self._zeros = groups, geo, zeros
        self._cache: Dict[int, _PixelImage] = {}

    def get_image(self, omega_idx: int, det_idx: int) -> _PixelImage:
        i = omega_idx * self._geo.n_det + det_idx
        img = self._cache.get(i)
        if img is None:
            img = self._cache[i] = _PixelImage(self._groups.get(i), self._geo, self._zeros)
        return img


class CoarseStack:
    """MultiScaleImageStack stand-in for DifferentiableCostFunction at ONE scale: binary images
    max-pooled by `factor` (a coarse pixel is lit if any of its factor x factor pixels is) and
    omega-blended (max over +-omega_window frames of the same detector, zero padded), densified
    on demand. Equals MultiScaleImageStack(...).get_at_scale(scale_index) for the same data."""

    def __init__(
        self,
        groups: Dict[int, np.ndarray],
        geo: Geometry,
        factor: int = FACTOR,
        omega_window: int = 1,
        scale_index: int = SCALE_INDEX,
    ):
        self.geo, self.factor, self.window, self.scale_index = (
            geo,
            factor,
            omega_window,
            scale_index,
        )
        self.H, self.W = geo.H // factor, geo.W // factor
        self._coarse: Dict[int, np.ndarray] = {}
        for i, pos in groups.items():
            rr, cc = pos // geo.W, pos % geo.W
            self._coarse[i] = np.unique((rr // factor) * self.W + (cc // factor))
        self._dense: Dict[int, torch.Tensor] = {}

    def get_at_scale(self, scale_idx: int) -> "CoarseStack":
        assert scale_idx == self.scale_index, "this stack only holds scale %d" % self.scale_index
        return self

    def _image(self, flat: int) -> torch.Tensor:
        t = self._dense.get(flat)
        if t is None:
            n_det = self.geo.n_det
            f, d = divmod(flat, n_det)
            t = torch.zeros(self.H * self.W, dtype=torch.float32)
            for g in range(max(0, f - self.window), min(self.geo.n_omega, f + self.window + 1)):
                pos = self._coarse.get(g * n_det + d)
                if pos is not None:
                    t[torch.from_numpy(pos)] = 1.0
            t = self._dense[flat] = t.reshape(1, self.H, self.W)
        return t

    def get_images_batch(self, flat_indices: torch.Tensor) -> torch.Tensor:
        return torch.stack([self._image(int(i)) for i in flat_indices.tolist()])


@contextlib.contextmanager
def seeded_default_rng(seed: int) -> Iterator[None]:
    """MCOptimizer draws from an unseeded np.random.default_rng() (the np.random.seed call in
    optimizer_baselines.py has no effect on it); make that generator deterministic for one run."""
    orig = np.random.default_rng

    def patched(*a: Any, **k: Any) -> np.random.Generator:
        return orig(seed) if not a and not k else orig(*a, **k)

    np.random.default_rng = patched  # type: ignore[assignment]
    try:
        yield
    finally:
        np.random.default_rng = orig  # type: ignore[assignment]


# ---------------------------------------------------------------------------
# Physics: target / source spot images, cost functions
# ---------------------------------------------------------------------------

_LIT_CACHE: Dict[Tuple[int, ...], frozenset] = {}


def spot_pixel_keys(obs: Any, delta_deg: np.ndarray, geo: Geometry) -> np.ndarray:
    """Sorted unique pixel keys of all spots of a BatchedObserver that are recorded at the offset
    delta_deg (3,), lit as the simulator's rasteriser lights them."""
    from icenine.orientation_eval import BatchedObserver, lit_pixel_set

    o = obs.observe(torch.as_tensor(np.asarray(delta_deg, dtype=np.float64)[None]))
    keys = BatchedObserver.vertex_keys(o)[0].tolist()
    present, frames = o.present[0].tolist(), o.frame[0].tolist()
    dets = obs.det_idx.tolist()
    ncols, nrows = obs.d_ncols.tolist(), obs.d_nrows.tolist()
    f_l: List[int] = []
    d_l: List[int] = []
    r_l: List[int] = []
    c_l: List[int] = []
    for m in range(obs.M):
        if not present[m] or frames[m] < 0:
            continue
        ck = (dets[m], ncols[m], nrows[m]) + tuple(keys[m])
        pix = _LIT_CACHE.get(ck)
        if pix is None:
            pix = _LIT_CACHE[ck] = lit_pixel_set(tuple(keys[m]), ncols[m], nrows[m])
        for c, r in pix:
            f_l.append(frames[m])
            d_l.append(dets[m])
            r_l.append(r)
            c_l.append(c)
    if len(_LIT_CACHE) > 500_000:
        _LIT_CACHE.clear()
    if not f_l:
        return np.zeros(0, dtype=np.int64)
    return np.unique(encode_pixels(np.array(f_l), np.array(d_l), np.array(r_l), np.array(c_l), geo))


_CTX: Optional[SimpleNamespace] = None


def init_worker(args_dict: Dict[str, Any]) -> None:
    """Per-process setup: physics, the two cost functions (images are attached per case)."""
    global _CTX
    torch.set_num_threads(1)
    from generate_toy_orientation_dataset import example_dir_for, setup_example
    from icenine.cost_functions import VoxelCostFunction
    from icenine.differentiable_cost import DifferentiableCostFunction

    args = SimpleNamespace(**args_dict)
    example_dir = example_dir_for(args.example)
    setup = setup_example(example_dir, max_q=args.max_q)
    mic, sample, detector_list, range_map, exp_setup, simulator, structure_list, get_vertices = (
        setup
    )
    common = dict(
        simulator=simulator,
        detector_list=detector_list,
        range_map=range_map,
        sample=sample,
        structure_list=structure_list,
        max_q=args.max_q,
        min_sin_eta=args.min_sin_eta,
        eta_limit=exp_setup.get_eta_limit(),  # the config's EtaLimit (86 deg), as the net's ROI set
    )
    hard_fn = VoxelCostFunction(exp_data=None, mode="hard", **common)  # type: ignore[arg-type]
    diff_fn = DifferentiableCostFunction(image_stack=None, **common)  # type: ignore[arg-type]
    n_omega = int(max(i for i in range_map.index_list if i is not None)) + 1
    det0 = detector_list[0]
    assert all(
        (d.num_rows, d.num_cols) == (det0.num_rows, det0.num_cols) for d in detector_list
    ), "detectors differ in size"
    geo = Geometry(len(detector_list), n_omega, det0.num_rows, det0.num_cols)
    assert geo.H % FACTOR == 0 and geo.W % FACTOR == 0
    _CTX = SimpleNamespace(
        args=args,
        example_dir=example_dir,
        setup=setup,
        hard_fn=hard_fn,
        diff_fn=diff_fn,
        geo=geo,
        zeros=np.zeros((geo.H, geo.W), dtype=np.uint8),
        voxels={},
    )


def voxel_context(ctx: SimpleNamespace, voxel_index: int) -> SimpleNamespace:
    """Per-voxel data (cached for a few voxels): truth, vertices, distractor sources, and the
    pixel keys of the target's spots at its true orientation."""
    cached = ctx.voxels.get(voxel_index)
    if cached is not None:
        return cached
    from generate_toy_orientation_dataset import build_distractor_sources, build_problem
    from icenine.orientation_eval import BatchedObserver

    a = ctx.args
    mic = ctx.setup[0]
    voxel = mic.voxels[voxel_index]
    vertices = ctx.setup[7](voxel)
    sources, _ = build_distractor_sources(
        ctx.example_dir,
        mic,
        voxel_index,
        ctx.setup,
        SimpleNamespace(neighbors=a.neighbors, twin=True, neighbor_radius_um=a.neighbor_radius_um),
        a.max_q,
        a.detectors,
        a.min_sin_eta,
    )
    pr = build_problem(
        ctx.example_dir,
        voxel_index,
        max_q=a.max_q,
        detectors=a.detectors,
        min_sin_eta=0.0,  # the image holds every spot; the cost applies the |sin eta| filter
        setup=ctx.setup,
    )
    obs = BatchedObserver(
        pr["R_nom"],
        pr["vertices"],
        pr["sample"],
        pr["detector_list"],
        pr["range_map"],
        pr["exp_setup"],
        pr["roi_list"],
    )
    vctx = SimpleNamespace(
        voxel=voxel,
        vertices=vertices,
        R_true=np.asarray(voxel.orientation, dtype=np.float64),
        sources=sources,
        target_keys=spot_pixel_keys(obs, np.zeros(3), ctx.geo),
    )
    if len(ctx.voxels) >= 3:
        ctx.voxels.pop(next(iter(ctx.voxels)))
    ctx.voxels[voxel_index] = vctx
    return vctx


def case_image_keys(
    ctx: SimpleNamespace,
    vctx: SimpleNamespace,
    variant: str,
    draw: Tuple[np.ndarray, np.ndarray],
    edit: Optional[Tuple[np.ndarray, np.ndarray]],
) -> np.ndarray:
    """All lit pixels of one case's detector images (sorted unique keys)."""
    keys = vctx.target_keys
    if variant == "clean":
        return keys
    off, act = draw
    parts = [keys]
    for s, src in enumerate(vctx.sources):
        if act[s]:
            parts.append(spot_pixel_keys(src, off[s], ctx.geo))
    keys = np.unique(np.concatenate(parts))
    assert edit is not None
    removed, added = edit
    return np.union1d(np.setdiff1d(keys, removed), added)


# ---------------------------------------------------------------------------
# Running the optimizers
# ---------------------------------------------------------------------------


def run_mc_adam(
    ctx: SimpleNamespace,
    vctx: SimpleNamespace,
    keys: np.ndarray,
    R_start: np.ndarray,
    r_deg: float,
    methods: Sequence[str],
    seed: int,
) -> Dict[str, Dict[str, float]]:
    """MC and/or Adam on one case's images. Returns method -> dict(R_final, seconds, q_start,
    q_true, q_final, n_evals)."""
    import bench_hp_sweep as hp  # benchmarks/bench_hp_sweep.py

    a, geo = ctx.args, ctx.geo
    groups = group_pixels(keys, geo)
    out: Dict[str, Dict[str, Any]] = {}
    voxel, vertices = vctx.voxel, vctx.vertices
    R0 = R_start.astype(np.float32)
    Rt = vctx.R_true.astype(np.float32)
    if "mc" in methods:
        fn = ctx.hard_fn
        fn.exp_data = PixelSetData(groups, geo, ctx.zeros)
        q_start = 1.0 - fn.evaluate(R0, vertices, voxel.phase).cost
        q_true = 1.0 - fn.evaluate(Rt, vertices, voxel.phase).cost
        with seeded_default_rng(seed):
            t0 = time.perf_counter()
            res = hp.run_one_mc(
                fn,
                voxel,
                vertices,
                R0,
                vctx.R_true,
                a.mc_steps,
                a.mc_restarts,
                a.mc_step_frac,
                np.radians(r_deg),
                record_traj=False,
            )
            dt = time.perf_counter() - t0
        out["mc"] = dict(
            R_final=res["R_final"],
            seconds=dt,
            q_start=q_start,
            q_true=q_true,
            q_final=res["final_quality"],
            n_evals=res["n_evaluations"],
        )
        fn.exp_data = None
    if "adam" in methods:
        fn = ctx.diff_fn
        fn.image_stack = CoarseStack(groups, geo, FACTOR, hp.FIXED_OW, SCALE_INDEX)
        with torch.no_grad():
            q_start = float(
                fn.evaluate(
                    torch.from_numpy(R0), vertices, voxel.phase, scale=hp.FIXED_SCALE
                ).quality
            )
            q_true = float(
                fn.evaluate(
                    torch.from_numpy(Rt), vertices, voxel.phase, scale=hp.FIXED_SCALE
                ).quality
            )
        assert hp.FIXED_SCALE == SCALE_INDEX
        t0 = time.perf_counter()
        res = hp.run_one_riemannian_adam_geoopt(
            fn, voxel, vertices, R0, hp.FIXED_SCALE, a.adam_steps, a.adam_lr, traj_subsample=0
        )
        dt = time.perf_counter() - t0
        out["adam"] = dict(
            R_final=res["R_final"],
            seconds=dt,
            q_start=q_start,
            q_true=q_true,
            q_final=res["final_quality"],
            n_evals=a.adam_steps + 1,
        )
        fn.image_stack = None
    return out


def run_gn(
    ctx: SimpleNamespace,
    voxel_index: int,
    prep1: List[Optional[Any]],
    b1: Dict[str, torch.Tensor],
    ok1: np.ndarray,
    delta0: np.ndarray,
    R_nom0: np.ndarray,
    R_true: np.ndarray,
    draws: List[Tuple[np.ndarray, np.ndarray]],
    sources: List[Any],
    variant: str,
    seed: int,
    huber_c: Optional[float],
) -> Dict[str, np.ndarray]:
    """Gauss-Newton (plain or Huber) with up to N_PASS re-centring passes on the sweep's windows,
    mirroring perturbation_sweep.sweep_voxel's pass loop. Returns arrays (D, N_PASS)."""
    from icenine.orientation_baselines import CentroidGaussNewton, extract_measurements

    a = ctx.args
    D = len(prep1)
    shape = (D, N_PASS)
    nan = lambda: np.full(shape, np.nan, dtype=np.float32)  # noqa: E731
    out = dict(
        err_angle=nan(),
        err_x=nan(),
        err_y=nan(),
        err_z=nan(),
        runtime=nan(),
        status=np.full(shape, -1, dtype=np.int8),
        converged=np.zeros(shape, dtype=bool),
        aux=np.zeros(shape, dtype=np.int32),  # spots used by the last solve
    )
    alive = ok1.copy()
    R_nom, delta_t = R_nom0.copy(), delta0.copy()
    est: Dict[int, Dict[str, Any]] = {}
    spent = np.zeros(D)
    for p_i in range(N_PASS):
        if p_i == 0:
            batch, preps = b1, prep1
        else:
            preps = [None] * D
            for j in np.nonzero(alive)[0]:
                pj, _why = ps.prepare_nominal(ctx, voxel_index, R_nom[j])
                if pj is None:
                    alive[j] = False
                preps[j] = pj
            batch = ps.render_batch(preps, delta_t, draws, sources, variant, seed, a)
            alive &= ~(alive & (batch["n_present"].numpy() < a.min_present))
        rows = np.nonzero(alive)[0]
        dh = np.zeros((len(rows), 3))
        for k, j in enumerate(rows):
            pj = preps[j]
            t0 = time.perf_counter()
            try:
                meas = extract_measurements(batch["windows"][j, : pj.n], pj.spec, pj.obs)
                res = CentroidGaussNewton(pj.obs, huber_c=huber_c).solve(meas)
                delta = np.asarray(res["delta"], dtype=np.float64)
                if not np.all(np.isfinite(delta)):
                    raise FloatingPointError
                status, conv, n_used = (
                    GN_STATUS[str(res["status"])],
                    bool(res["converged"]),
                    res["n_used"],
                )
            except (np.linalg.LinAlgError, FloatingPointError):
                delta, status, conv, n_used = np.zeros(3), GN_STATUS["error"], False, 0
            spent[j] += time.perf_counter() - t0
            dh[k] = delta
            est[j] = dict(status=status, conv=conv, n_used=n_used)
            out["status"][j, p_i], out["converged"][j, p_i] = status, conv
            out["aux"][j, p_i] = n_used
        if len(rows):
            Rn, dt_ = ps.recentre(R_nom[rows], dh, R_true)
            e = ps.estimate_error_deg(Rn, R_true)
            for k, j in enumerate(rows):
                est[j].update(angle=np.linalg.norm(e[k]), x=e[k, 0], y=e[k, 1], z=e[k, 2])
            R_nom[rows], delta_t[rows] = Rn, dt_
        for j, c in est.items():  # carried-forward estimate
            out["err_angle"][j, p_i], out["err_x"][j, p_i] = c["angle"], c["x"]
            out["err_y"][j, p_i], out["err_z"][j, p_i] = c["y"], c["z"]
            out["runtime"][j, p_i] = spent[j]
    return out


def sweep_task(item: Tuple[int, int, int, str, np.ndarray, np.ndarray]) -> Tuple[int, int, float]:
    """One (voxel, radius): every variant and method; cache the result file."""
    vidx, vpos, ri, path, ref_nroi, ref_fail = item
    t_start = time.time()
    ctx = _CTX
    assert ctx is not None
    a = ctx.args
    vctx = voxel_context(ctx, vidx)
    r = a.radii[ri]
    D, P = a.n_dirs, N_PASS
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
    n_roi = np.array([p.n if p is not None else 0 for p in prep1])
    assert (n_roi == ref_nroi).all(), "ROI counts differ from the sweep: cases are not aligned"
    shape = (D, len(VARIANTS), len(METHODS), P)
    nan = lambda: np.full(shape, np.nan, dtype=np.float32)  # noqa: E731
    res = dict(
        err_angle=nan(),
        err_x=nan(),
        err_y=nan(),
        err_z=nan(),
        runtime=nan(),
        status=np.full(shape, -1, dtype=np.int8),
        converged=np.zeros(shape, dtype=bool),
        aux=np.zeros(shape, dtype=np.int32),
        q_start=np.full((D, len(VARIANTS), 2), np.nan, dtype=np.float32),
        q_true=np.full((D, len(VARIANTS), 2), np.nan, dtype=np.float32),
        q_final=np.full((D, len(VARIANTS), 2), np.nan, dtype=np.float32),
        fail_pass1=np.zeros((D, len(VARIANTS)), dtype=np.int8),
    )
    img_methods = [m for m in IMAGE_METHODS if m in a.methods]
    for vi, variant in enumerate(VARIANTS):
        seed = a.realism_seed + 1000003 * vpos + 1009 * ri
        layers: Dict[str, torch.Tensor] = {}
        b1 = ps.render_batch(prep1, delta0, draws, vctx.sources, variant, seed, a, layers=layers)
        have = np.array([p is not None for p in prep1])
        ok1 = have & (b1["n_present"].numpy() >= a.min_present)
        fail = np.where(reason1 > 0, reason1, np.where(ok1, 0, 4)).astype(np.int8)
        assert (fail == ref_fail[:, vi]).all(), "pass-1 failures differ from the sweep"
        res["fail_pass1"][:, vi] = fail
        for mi, name in enumerate(METHODS):
            if name in ("gn", "huber") and name in a.methods:
                g = run_gn(
                    ctx, vidx, prep1, b1, ok1, delta0, R_nom0, vctx.R_true, draws,
                    vctx.sources, variant, seed, None if name == "gn" else a.huber_c,
                )  # fmt: skip
                for k, v in g.items():
                    res[k][:, vi, mi, :] = v
        if img_methods:
            for j in range(min(a.opt_dirs, D)):
                if not have[j]:
                    continue
                edit = None
                if variant != "clean":
                    p = prep1[j]
                    n = p.n
                    edit = realism_edit(
                        layers["clean"][j, :n].numpy(),
                        layers["dis"][j, :n].numpy(),
                        b1["windows"][j, :n].numpy(),
                        p.spec,
                        p.obs.det_idx.numpy(),
                        a.frame_half_width,
                        ctx.geo,
                    )
                keys = case_image_keys(ctx, vctx, variant, draws[j], edit)
                seed_j = int(seed + 7919 * j)
                out = run_mc_adam(ctx, vctx, keys, R_nom0[j], float(r), img_methods, seed_j)
                for name, o in out.items():
                    mi = METHODS.index(name)
                    e = ps.estimate_error_deg(o["R_final"].astype(np.float64), vctx.R_true)
                    res["err_angle"][j, vi, mi, 0] = np.linalg.norm(e)
                    res["err_x"][j, vi, mi, 0], res["err_y"][j, vi, mi, 0] = e[0], e[1]
                    res["err_z"][j, vi, mi, 0] = e[2]
                    res["runtime"][j, vi, mi, 0] = o["seconds"]
                    res["status"][j, vi, mi, 0] = 0
                    res["converged"][j, vi, mi, 0] = True
                    res["aux"][j, vi, mi, 0] = o["n_evals"]
                    k = IMAGE_METHODS.index(name)
                    res["q_start"][j, vi, k] = o["q_start"]
                    res["q_true"][j, vi, k] = o["q_true"]
                    res["q_final"][j, vi, k] = o["q_final"]
    np.savez(path, **res)
    return vidx, ri, time.time() - t_start


# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------


def load_sweep_config(sweep_raw: Path) -> Tuple[Dict[str, Any], Dict[str, np.ndarray]]:
    raw = dict(np.load(sweep_raw, allow_pickle=False))
    cfg = json.loads(str(raw["config_json"]))
    return cfg, raw


def do_run(args: argparse.Namespace) -> None:
    import multiprocessing as mp

    out_dir = Path(args.out_dir).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    cache = Path(args.cache_dir).resolve()
    cache.mkdir(parents=True, exist_ok=True)
    cfg, sweep = load_sweep_config(out_dir / "perturbation_sweep_raw.npz")
    wargs = dict(cfg)
    wargs.update(
        methods=args.methods,
        opt_dirs=args.opt_dirs,
        mc_steps=args.mc_steps,
        mc_restarts=args.mc_restarts,
        mc_step_frac=args.mc_step_frac,
        adam_lr=args.adam_lr,
        adam_steps=args.adam_steps,
        huber_c=args.huber_c,
    )
    voxels = [int(v) for v in sweep["voxel_indices"]][: args.n_voxels or None]
    radii_idx = args.radii_idx if args.radii_idx else list(range(len(cfg["radii"])))
    items = []
    for vpos, v in enumerate(voxels):
        for ri in radii_idx:
            items.append(
                (
                    v,
                    vpos,
                    ri,
                    str(cache / f"v{v}_r{ri}.npz"),
                    sweep["n_roi"][vpos, ri],
                    sweep["fail_pass1"][vpos, ri],
                )
            )
    todo = [it for it in items if not Path(it[3]).exists()]
    print(
        f"{len(items)} tasks ({len(items) - len(todo)} cached), {args.workers} workers; "
        f"methods {args.methods}, optimizers on the first {args.opt_dirs} directions",
        flush=True,
    )
    t0 = time.time()
    ctx = mp.get_context("spawn")
    with ctx.Pool(args.workers, initializer=init_worker, initargs=(wargs,)) as pool:
        for k, (vidx, ri, secs) in enumerate(pool.imap_unordered(sweep_task, todo), 1):
            print(
                f"  [{k}/{len(todo)}] voxel {vidx} r#{ri} done in {secs:.0f}s "
                f"({time.time() - t0:.0f}s total)",
                flush=True,
            )
    if args.no_assemble:
        return
    assemble(out_dir, cache, voxels, radii_idx, cfg, wargs)


def assemble(
    out_dir: Path,
    cache: Path,
    voxels: List[int],
    radii_idx: List[int],
    cfg: Dict[str, Any],
    wargs: Dict[str, Any],
) -> None:
    """Stack the per-(voxel, radius) cache files into optimizer_sweep_raw.npz."""
    parts = {(v, ri): np.load(cache / f"v{v}_r{ri}.npz") for v in voxels for ri in radii_idx}
    first = next(iter(parts.values()))
    raw: Dict[str, np.ndarray] = {}
    for k in first.files:
        raw[k] = np.stack(
            [np.stack([parts[(v, ri)][k] for ri in radii_idx]) for v in voxels]
        )  # (V, R, D, ...)
    raw["voxel_indices"] = np.array(voxels)
    raw["radii"] = np.array([cfg["radii"][i] for i in radii_idx])
    raw["methods"] = np.array(METHODS)
    raw["variants"] = np.array(VARIANTS)
    raw["config_json"] = np.array(json.dumps(wargs, default=str))
    np.savez_compressed(out_dir / "optimizer_sweep_raw.npz", **raw)
    print(f"saved {out_dir / 'optimizer_sweep_raw.npz'}")


def build_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    r = sub.add_parser("run", help="run the optimizers (resumable) and assemble the raw npz")
    r.add_argument("--out-dir", default="benchmarks/toy_orientation_sweep")
    r.add_argument("--cache-dir", default="scripts/optimizer_sweep_cache")
    r.add_argument("--methods", nargs="+", default=METHODS, choices=METHODS)
    r.add_argument("--n-voxels", type=int, default=0, help="first N sweep voxels (0 = all)")
    r.add_argument(
        "--radii-idx", type=int, nargs="*", default=[], help="radius indices (default all)"
    )
    r.add_argument(
        "--opt-dirs",
        type=int,
        default=20,
        help="MC / Adam on the first N of the 20 directions per radius (GN always on all)",
    )
    r.add_argument("--mc-steps", type=int, default=3500)
    r.add_argument("--mc-restarts", type=int, default=2)
    r.add_argument("--mc-step-frac", type=float, default=0.5)
    r.add_argument("--adam-lr", type=float, default=1e-4)
    r.add_argument("--adam-steps", type=int, default=100)
    r.add_argument("--huber-c", type=float, default=1.0)
    r.add_argument("--workers", type=int, default=10)
    r.add_argument("--no-assemble", action="store_true")
    s = sub.add_parser("summarize", help="tables, JSON and plot from the raw npz files")
    s.add_argument("--out-dir", default="benchmarks/toy_orientation_sweep")
    return ap


def main() -> None:
    args = build_parser().parse_args()
    if args.cmd == "summarize":
        from optimizer_sweep_summary import do_summarize

        do_summarize(Path(args.out_dir).resolve())
    else:
        do_run(args)


if __name__ == "__main__":
    main()
