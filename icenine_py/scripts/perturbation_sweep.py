#!/usr/bin/env python3
"""
Single-voxel perturbation sweep for the toy orientation network.

Randomly sampled ManyGrains voxels (not among the 30 used to build the datasets) are each tested on
their own: the nominal orientation is the voxel's true .mic orientation rotated by a random
rotation of angle r, so the truth relative to the nominal is a rotation vector delta of length
exactly r (delta = rotvec(R_true R_nom^T), the dataset convention R(delta) = exp([delta]x) R_nom).
Windows are rendered per voxel exactly as scripts/generate_toy_orientation_dataset.py does (observer
renderer, 32x32 px, K = 4 frames, frame coding, the peak context and sub-pixel nominal offsets
computed at the nominal orientation); nothing renders the whole sample. The `all` variant adds the
generator's distractor layer (2 nearest neighbours + a Sigma3 twin) and make_realistic_dataset.

One-shot = one pass of the net. Iterated = up to --passes passes, each re-centring at the current
estimate (new nominal = previous estimate; ROI set, windows and context are rebuilt there and the
windows of the SAME true orientation re-rendered, with the same distractor draws and the same
realism seed).

Usage (from icenine_py/):
  uv run python scripts/perturbation_sweep.py run --models clean_s0=scripts/..._clean_s0.pt ... \
      --out-dir benchmarks/toy_orientation_sweep
  uv run python scripts/perturbation_sweep.py summarize --out-dir benchmarks/toy_orientation_sweep
"""

import argparse
import json
import sys
import time
from pathlib import Path
from types import SimpleNamespace
from typing import Any, Callable, Dict, List, Optional, Sequence, Tuple

import numpy as np
import torch
from scipy.spatial import cKDTree
from scipy.spatial.transform import Rotation

sys.path.insert(0, str(Path(__file__).parent))

DEG = np.pi / 180.0
RADII_DEG = [0.05, 0.1, 0.25, 0.5, 0.75, 1.0, 1.5, 2.0, 3.0, 5.0]
VARIANTS = ["clean", "all"]
STOP_REASONS = {0: "ok", 1: "no_roi", 2: "few_roi", 3: "spec", 4: "few_present"}
M_PAD = 130  # peak count the datasets are padded to


# ---------------------------------------------------------------------------
# Pure helpers (no physics): voxel sampling, perturbation, re-centring, errors
# ---------------------------------------------------------------------------


def excluded_voxel_indices(dataset_paths: Sequence[str]) -> np.ndarray:
    """Union of the mic voxel indices used by the given multi-voxel dataset files."""
    out: List[int] = []
    for p in dataset_paths:
        d = torch.load(Path(p).resolve(), mmap=True)
        out += [int(i) for i in d["voxel_indices"]]
    return np.unique(np.array(out, dtype=np.int64))


def eligible_candidates(
    r_perp_um: np.ndarray, excluded: Sequence[int], r_max_um: float
) -> np.ndarray:
    """Mic voxel indices with r_perp <= r_max_um that are not in `excluded`."""
    keep = np.asarray(r_perp_um) <= r_max_um
    keep[np.asarray(list(excluded), dtype=np.int64)] = False
    return np.nonzero(keep)[0]


def sample_voxels(
    candidates: np.ndarray, n: int, seed: int, usable: Callable[[int], bool]
) -> List[int]:
    """n voxels drawn uniformly at random (seeded) without replacement from `candidates`: the
    candidates are put in a seeded random order and the first n for which usable(idx) is true are
    kept, so a rejected voxel is replaced by the next one in the same order."""
    rng = np.random.default_rng(seed)
    out: List[int] = []
    for idx in rng.permutation(candidates):
        if usable(int(idx)):
            out.append(int(idx))
            if len(out) == n:
                break
    return out


def grain_ids(orientations: np.ndarray, decimals: int = 5) -> np.ndarray:
    """Grain id per voxel from identical orientation matrices (voxels of a grain share one)."""
    keys = np.round(np.asarray(orientations).reshape(len(orientations), 9), decimals)
    _, inv = np.unique(keys, axis=0, return_inverse=True)
    return inv.reshape(-1).astype(np.int64)


def near_boundary(pos_um: np.ndarray, grain: np.ndarray, radius_um: float) -> np.ndarray:
    """True where another grain has a voxel within radius_um (sample plane)."""
    tree = cKDTree(pos_um)
    out = np.zeros(len(pos_um), dtype=bool)
    for i, nb in enumerate(tree.query_ball_point(pos_um, radius_um)):
        out[i] = bool((grain[nb] != grain[i]).any())
    return out


def random_rotvecs(n: int, angle_deg: float, rng: np.random.Generator) -> np.ndarray:
    """(n, 3) rotation vectors (degrees) of length exactly angle_deg, uniform random axes."""
    d = rng.normal(size=(n, 3))
    d /= np.linalg.norm(d, axis=1, keepdims=True)
    return d * angle_deg


def relative_offset_deg(R_true: np.ndarray, R_nom: np.ndarray) -> np.ndarray:
    """delta (deg) with exp([delta]x) R_nom = R_true, i.e. rotvec(R_true R_nom^T). Batched or not."""
    return Rotation.from_matrix(np.asarray(R_true) @ np.swapaxes(R_nom, -1, -2)).as_rotvec() / DEG


def perturbed_nominal(R_true: np.ndarray, delta_deg: np.ndarray) -> np.ndarray:
    """Nominal orientation such that the truth is delta away from it: R_nom = exp(-[delta]x) R_true,
    so exp([delta]x) R_nom = R_true (the dataset convention, orientation_eval.offsets_to_matrices).
    """
    from icenine.orientation_eval import offsets_to_matrices

    return offsets_to_matrices(-np.asarray(delta_deg, dtype=np.float64), R_true)


def recentre(R_nom: np.ndarray, delta_hat_deg: np.ndarray, R_true: np.ndarray):
    """New nominal = current estimate exp([delta_hat]x) R_nom, and the unchanged truth's offset
    from it (what the re-rendered windows must be drawn at). Returns (R_new, delta_true_new)."""
    from icenine.orientation_eval import offsets_to_matrices

    R_new = offsets_to_matrices(np.asarray(delta_hat_deg, dtype=np.float64), R_nom)
    return R_new, relative_offset_deg(R_true, R_new)


def estimate_error_deg(R_est: np.ndarray, R_true: np.ndarray) -> np.ndarray:
    """(.., 3) rotation vector (deg, sample frame) of R_est R_true^T; its norm is the angle."""
    return Rotation.from_matrix(np.asarray(R_est) @ np.swapaxes(R_true, -1, -2)).as_rotvec() / DEG


def mahalanobis_sq(truth_deg: np.ndarray, pred_deg: np.ndarray, chol: np.ndarray) -> np.ndarray:
    """Squared Mahalanobis distance of the truth under N(pred, chol chol^T), as the trainer."""
    r = (np.asarray(truth_deg, float) - np.asarray(pred_deg, float))[..., None]
    z = np.linalg.solve(np.asarray(chol, float), r)[..., 0]
    return (z**2).sum(-1)


def load_model(path: str) -> Tuple[torch.nn.Module, Dict[str, Any]]:
    """Rebuild a network saved by train_toy_orientation_nn.py --save-model (eval mode, CPU)."""
    from icenine.toy_orientation_model import GNLayerNet

    ck = torch.load(Path(path).resolve(), map_location="cpu", weights_only=False)
    assert ck["arch"] == "gn", "the sweep only supports --arch gn checkpoints"
    net = GNLayerNet(**ck["model_kwargs"])
    net.load_state_dict(ck["state_dict"])
    net.eval()
    return net, ck


# ---------------------------------------------------------------------------
# Physics: one nominal -> observer, window spec, context, nominal offsets
# ---------------------------------------------------------------------------


class Prepared(SimpleNamespace):
    """Everything fixed by one nominal orientation of one voxel."""


def prepare_nominal(
    ctx: SimpleNamespace, voxel_index: int, R_nom: np.ndarray
) -> Tuple[Optional[Prepared], int]:
    """ROI set, observer, window spec, context and nominal offsets at R_nom. Returns
    (prepared, 0) or (None, stop_reason) with reason 1 (no ROI peak), 2 (fewer than min_peaks
    ROI entries) or 3 (window spec failed: a ROI spot is not present at the nominal)."""
    from generate_toy_orientation_dataset import build_problem
    from icenine.orientation_eval import (
        BatchedObserver,
        WindowSpec,
        nominal_offsets,
    )

    a = ctx.args
    try:
        pr = build_problem(
            ctx.example_dir,
            voxel_index,
            max_q=a.max_q,
            detectors=a.detectors,
            min_sin_eta=a.min_sin_eta,
            setup=ctx.setup,
            orientation=R_nom,
        )
    except RuntimeError:
        return None, 1
    if len(pr["roi_list"]) < a.min_peaks:
        return None, 2
    obs = BatchedObserver(
        pr["R_nom"],
        pr["vertices"],
        pr["sample"],
        pr["detector_list"],
        pr["range_map"],
        pr["exp_setup"],
        pr["roi_list"],
    )
    try:
        spec = WindowSpec.from_nominal(obs, a.window_size, a.frame_half_width)
    except ValueError:
        return None, 3
    return (
        Prepared(
            obs=obs,
            spec=spec,
            n=len(pr["roi_list"]),
            context=obs.peak_context().float(),
            nom_off=torch.from_numpy(np.asarray(nominal_offsets(obs))).float(),
            R_nom=np.asarray(R_nom, dtype=np.float64),
        ),
        0,
    )


def render_batch(
    preps: List[Optional[Prepared]],
    delta_true: np.ndarray,
    dis_draws: Optional[List[Tuple[np.ndarray, np.ndarray]]],
    sources: List[Any],
    variant: str,
    seed: int,
    args: SimpleNamespace,
) -> Dict[str, torch.Tensor]:
    """Windows (+ context, nominal offsets, valid mask) for a batch of cases, each with its own
    prepared nominal (None = case not run: zero windows). delta_true (B, 3) is the truth's offset
    from each case's nominal. For variant "all" the distractor layer of each case uses its fixed
    draws dis_draws[j] = (source offsets (S, 3) deg, source active (S,)), and the whole batch goes
    through make_realistic_dataset with a fixed seed and the padding mask."""
    from icenine.orientation_eval import (
        make_realistic_dataset,
        render_distractor_windows,
        render_windows,
    )

    B, W = len(preps), args.window_size
    M = max([M_PAD] + [p.n for p in preps if p is not None])
    win = torch.zeros(B, M, W, W, dtype=torch.uint8)
    dis = torch.zeros(B, M, W, W, dtype=torch.uint8)
    ctx = torch.zeros(B, M, preps[0].context.shape[-1] if preps[0] is not None else 16)
    nom = torch.zeros(B, M, 3)
    valid = torch.zeros(B, M, dtype=torch.bool)
    for j, p in enumerate(preps):
        if p is None:
            continue
        if ctx.shape[-1] != p.context.shape[-1]:
            ctx = torch.zeros(B, M, p.context.shape[-1])
        w, _ = render_windows(p.obs, p.spec, delta_true[j : j + 1])
        win[j, : p.n] = w[0]
        if variant == "all" and sources:
            off, act = dis_draws[j]  # type: ignore[index]
            d = render_distractor_windows(
                p.obs,
                p.spec,
                sources,
                [off[s : s + 1] for s in range(len(sources))],
                source_active=[act[s : s + 1] for s in range(len(sources))],
            )
            dis[j, : p.n] = d[0]
        ctx[j, : p.n] = p.context
        nom[j, : p.n] = p.nom_off
        valid[j, : p.n] = True
    if variant != "clean":
        win = make_realistic_dataset(
            win, dis, variant, args.frame_half_width, seed=seed, valid=valid
        )
    n_present = ((win > 0).flatten(2).any(-1) & valid).sum(1)
    return dict(windows=win, context=ctx, nom_off=nom, valid=valid, n_present=n_present)


def run_net(net: torch.nn.Module, batch: Dict[str, torch.Tensor], rows: np.ndarray, K: int):
    """Network on the selected cases. Returns (delta_hat (n, 3) deg, chol (n, 3, 3)) float64."""
    from icenine.orientation_eval import decode_windows

    with torch.no_grad():
        x = decode_windows(batch["windows"][rows], K)
        m, L = net(x, batch["context"][rows], {"nom_off": batch["nom_off"][rows]})
    return m.double().numpy(), L.double().numpy()


# ---------------------------------------------------------------------------
# Per-voxel sweep
# ---------------------------------------------------------------------------

_CTX: Optional[SimpleNamespace] = None


def init_worker(args_dict: Dict[str, Any]) -> None:
    """Per-process setup: physics, mic, models (CPU, one thread)."""
    global _CTX
    torch.set_num_threads(1)
    from generate_toy_orientation_dataset import example_dir_for, setup_example

    args = SimpleNamespace(**args_dict)
    example_dir = example_dir_for(args.example)
    setup = setup_example(example_dir, max_q=args.max_q)
    models = {name: load_model(path)[0] for name, path in args.models}
    _CTX = SimpleNamespace(args=args, example_dir=example_dir, setup=setup, models=models)


def sweep_voxel(voxel_index: int, voxel_pos: int) -> Dict[str, np.ndarray]:
    """All radii / directions / variants / models / passes for one voxel. voxel_pos is the voxel's
    position in the sampled list (used only for the realism seeds' spacing)."""
    from generate_toy_orientation_dataset import build_distractor_sources

    ctx = _CTX
    assert ctx is not None
    a = ctx.args
    mic = ctx.setup[0]
    R_true = np.asarray(mic.voxels[voxel_index].orientation, dtype=np.float64)
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
    S = len(sources)
    sigma_comp = a.neighbor_sigma_deg / np.sqrt(3.0)
    radii, D, P = a.radii, a.n_dirs, a.passes
    names = [n for n, _ in a.models]
    shape = (len(radii), D, len(VARIANTS), len(names), P)
    f32 = lambda: np.full(shape, np.nan, dtype=np.float32)  # noqa: E731
    out = dict(
        err_angle=f32(),
        err_x=f32(),
        err_y=f32(),
        err_z=f32(),
        maha2=f32(),
        delta_hat_norm=f32(),
        n_present=np.zeros(shape, dtype=np.int16),
        ran=np.zeros(shape, dtype=bool),
        fail_pass1=np.zeros((len(radii), D, len(VARIANTS)), dtype=np.int8),
        stop_reason=np.zeros(shape[:4], dtype=np.int8),
        stop_pass=np.full(shape[:4], -1, dtype=np.int8),
        n_roi=np.zeros((len(radii), D), dtype=np.int16),
    )
    K = a.frame_half_width
    for ri, r in enumerate(radii):
        rng = np.random.default_rng([a.sweep_seed, voxel_index, ri])
        delta0 = random_rotvecs(D, r, rng)
        src_off = rng.normal(0.0, sigma_comp, size=(D, S, 3))
        src_act = rng.random((D, S)) < a.neighbor_p
        draws = [(src_off[j], src_act[j]) for j in range(D)]
        R_nom0 = perturbed_nominal(R_true, delta0)
        prep1: List[Optional[Prepared]] = []
        reason1 = np.zeros(D, dtype=np.int8)
        for j in range(D):
            p, why = prepare_nominal(ctx, voxel_index, R_nom0[j])
            prep1.append(p)
            reason1[j] = why
            out["n_roi"][ri, j] = p.n if p is not None else 0
        for vi, variant in enumerate(VARIANTS):
            seed = a.realism_seed + 1000003 * voxel_pos + 1009 * ri
            b1 = render_batch(prep1, delta0, draws, sources, variant, seed, a)
            ok1 = np.array([p is not None for p in prep1]) & (
                b1["n_present"].numpy() >= a.min_present
            )
            fail = np.where(reason1 > 0, reason1, np.where(ok1, 0, 4)).astype(np.int8)
            out["fail_pass1"][ri, :, vi] = fail
            for mi, name in enumerate(names):
                net = ctx.models[name]
                alive = ok1.copy()
                R_nom = R_nom0.copy()
                delta_t = delta0.copy()
                est_cache: Dict[int, Dict[str, float]] = {}
                for p_i in range(P):
                    if p_i == 0:
                        batch, preps = b1, prep1
                    else:
                        preps = [None] * D
                        for j in np.nonzero(alive)[0]:
                            pj, why = prepare_nominal(ctx, voxel_index, R_nom[j])
                            if pj is None:
                                alive[j] = False
                                out["stop_reason"][ri, j, vi, mi] = why
                                out["stop_pass"][ri, j, vi, mi] = p_i
                            preps[j] = pj
                        batch = render_batch(preps, delta_t, draws, sources, variant, seed, a)
                        low = alive & (batch["n_present"].numpy() < a.min_present)
                        for j in np.nonzero(low)[0]:
                            alive[j] = False
                            out["stop_reason"][ri, j, vi, mi] = 4
                            out["stop_pass"][ri, j, vi, mi] = p_i
                    rows = np.nonzero(alive)[0]
                    if len(rows):
                        dh, L = run_net(net, batch, rows, K)
                        Rn, dt = recentre(R_nom[rows], dh, R_true)
                        e = estimate_error_deg(Rn, R_true)
                        m2 = mahalanobis_sq(delta_t[rows], dh, L)
                        for k, j in enumerate(rows):
                            est_cache[j] = dict(
                                angle=np.linalg.norm(e[k]),
                                x=e[k, 0],
                                y=e[k, 1],
                                z=e[k, 2],
                                maha2=m2[k],
                                dh=np.linalg.norm(dh[k]),
                                npres=int(batch["n_present"][j]),
                            )
                            out["ran"][ri, j, vi, mi, p_i] = True
                        R_nom[rows], delta_t[rows] = Rn, dt
                    # carried-forward estimate for every case that has one
                    for j, c in est_cache.items():
                        out["err_angle"][ri, j, vi, mi, p_i] = c["angle"]
                        out["err_x"][ri, j, vi, mi, p_i] = c["x"]
                        out["err_y"][ri, j, vi, mi, p_i] = c["y"]
                        out["err_z"][ri, j, vi, mi, p_i] = c["z"]
                        out["maha2"][ri, j, vi, mi, p_i] = c["maha2"]
                        out["delta_hat_norm"][ri, j, vi, mi, p_i] = c["dh"]
                        out["n_present"][ri, j, vi, mi, p_i] = c["npres"]
    return out


def work(item: Tuple[int, int, str]) -> Tuple[int, float]:
    """Pool task: sweep one voxel and cache the result file. Returns (voxel_index, seconds)."""
    vidx, vpos, path = item
    t0 = time.time()
    res = sweep_voxel(vidx, vpos)
    np.savez(path, **res)
    return vidx, time.time() - t0


# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------

METRICS = [
    "n_ok",
    "median_angle",
    "rms_angle",
    "rms_z",
    "rms_perp",
    "frac_lt_0p1",
    "frac_improved",
    "mean_maha2",
]


def case_metrics(raw: Dict[str, np.ndarray], sel: np.ndarray, vi: int, mi: int, ri: int, p: int):
    """Metrics over cases sel (V, D bool, pass-1 successes) at radius index ri, pass p."""
    r = raw["radii"][ri]
    ang = raw["err_angle"][:, ri, :, vi, mi, p][sel]
    ex, ey, ez = (raw[k][:, ri, :, vi, mi, p][sel] for k in ("err_x", "err_y", "err_z"))
    m2 = raw["maha2"][:, ri, :, vi, mi, p][sel]
    n = len(ang)
    if n == 0:
        return {k: float("nan") for k in METRICS} | {"n_ok": 0}
    return dict(
        n_ok=int(n),
        median_angle=float(np.median(ang)),
        rms_angle=float(np.sqrt(np.mean(ang.astype(np.float64) ** 2))),
        rms_z=float(np.sqrt(np.mean(ez.astype(np.float64) ** 2))),
        rms_perp=float(np.sqrt(np.mean((ex.astype(np.float64) ** 2 + ey**2) / 2.0))),
        frac_lt_0p1=float(np.mean(ang < 0.1)),
        frac_improved=float(np.mean(ang < r)),
        mean_maha2=float(np.mean(m2)),
    )


def summarize(raw: Dict[str, np.ndarray]) -> Dict[str, Any]:
    """Nested dict summary[group][variant][model][radius][pass] -> metrics, groups all/tercile_k."""
    radii = raw["radii"]
    models = [str(m) for m in raw["models"]]
    r_perp = raw["voxel_r_perp_um"]
    edges = np.quantile(r_perp, [1 / 3, 2 / 3])
    terc = np.digitize(r_perp, edges)
    groups: Dict[str, np.ndarray] = {"all": np.ones(len(r_perp), dtype=bool)}
    for t in range(3):
        groups[f"rperp_tercile{t + 1}"] = terc == t
    out: Dict[str, Any] = {
        "tercile_edges_um": edges.tolist(),
        "tercile_ranges_um": [
            (
                [float(r_perp[terc == t].min()), float(r_perp[terc == t].max())]
                if (terc == t).any()
                else [float("nan")] * 2
            )
            for t in range(3)
        ],
        "tercile_n_voxels": [int((terc == t).sum()) for t in range(3)],
        "failures": {},
        "stops": {},
        "metrics": {},
    }
    for g, gm in groups.items():
        out["metrics"][g] = {}
        for vi, v in enumerate(VARIANTS):
            out["metrics"][g][v] = {}
            for mi, m in enumerate(models):
                out["metrics"][g][v][m] = {}
                for ri, r in enumerate(radii):
                    ok = (raw["fail_pass1"][:, ri, :, vi] == 0) & gm[:, None]
                    out["metrics"][g][v][m][f"{r:g}"] = [
                        case_metrics(raw, ok, vi, mi, ri, p)
                        for p in range(raw["err_angle"].shape[-1])
                    ]
    for vi, v in enumerate(VARIANTS):
        out["failures"][v] = {}
        for ri, r in enumerate(radii):
            f = raw["fail_pass1"][:, ri, :, vi]
            out["failures"][v][f"{r:g}"] = {
                "n_cases": int(f.size),
                "n_fail": int((f > 0).sum()),
                **{STOP_REASONS[k]: int((f == k).sum()) for k in (1, 2, 3, 4)},
            }
        out["stops"][v] = {}
        for mi, m in enumerate(models):
            sr = raw["stop_reason"][:, :, :, vi, mi]
            out["stops"][v][m] = {
                f"{r:g}": {STOP_REASONS[k]: int((sr[:, ri] == k).sum()) for k in (1, 2, 3, 4)}
                for ri, r in enumerate(radii)
            }
    return out


MODEL_TYPES = {"clean": "clean-trained", "realistic": "realism-trained"}


def seed_groups(models: Sequence[str]) -> Dict[str, List[int]]:
    """Model type (name before the seed suffix) -> indices of its seeds."""
    g: Dict[str, List[int]] = {}
    for i, m in enumerate(models):
        g.setdefault(m.rsplit("_s", 1)[0], []).append(i)
    return g


def format_summary(raw: Dict[str, np.ndarray], summ: Dict[str, Any]) -> str:
    """Plain-text tables: per window variant and model type, one row per radius."""
    radii = [f"{r:g}" for r in raw["radii"]]
    models = [str(m) for m in raw["models"]]
    P = raw["err_angle"].shape[-1]
    mt = summ["metrics"]
    n_vox, n_dir = raw["fail_pass1"].shape[0], raw["fail_pass1"].shape[2]
    lines: List[str] = [
        f"Perturbation sweep: {n_vox} voxels x {len(radii)} radii x {n_dir} directions = "
        f"{n_vox * n_dir} cases per radius and variant; models {', '.join(models)}; "
        f"up to {P} passes (pass 1 = one-shot, pass {P} = iterated).",
        "Error = angle (deg) between the estimate and the true orientation, over the cases where "
        "the net could run at pass 1 (n_ok). Per-seed values are given as s0/s1.",
        "",
    ]
    for v in VARIANTS:
        lines.append(f"=== windows: {v} ===")
        fl = summ["failures"][v]
        lines.append("failures at pass 1 (net could not run), counts per radius:")
        lines.append("  r(deg)       " + " ".join(f"{r:>6}" for r in radii))
        lines.append("  n_fail       " + " ".join(f"{fl[r]['n_fail']:>6d}" for r in radii))
        for reason in ("no_roi", "few_roi", "spec", "few_present"):
            if any(fl[r][reason] for r in radii):
                lines.append(f"    {reason:<11}" + " ".join(f"{fl[r][reason]:>6d}" for r in radii))
        for m in models[:1]:
            stops = summ["stops"][v][m]
            tot = sum(sum(stops[r].values()) for r in radii)
            lines.append(f"  later passes stopped early (model {m}): {tot} case-stops")
        lines.append("")
        for tname, idx in seed_groups(models).items():
            head = MODEL_TYPES.get(tname, tname)
            lines.append(f"-- {head} (seeds {len(idx)}); {v} windows")
            lines.append(
                f"{'r':>5} {'n_ok':>5} | {'median one-shot s0/s1 (mean)':>32} | "
                f"{'median pass ' + str(P) + ' s0/s1 (mean)':>32} | {'<0.1deg':>13} | "
                f"{'err<r':>13} | {'Mah^2':>13} | {'rms z/perp (iter)':>17}"
            )
            for r in radii:

                def val(stat, p):
                    return [mt["all"][v][models[i]][r][p][stat] for i in idx]

                def pair(stat, p):
                    x = val(stat, p)
                    return "/".join(f"{y:.4f}" for y in x) + f" ({np.mean(x):.4f})"

                def mean2(stat):
                    return f"{np.mean(val(stat, 0)):.3f}/{np.mean(val(stat, P - 1)):.3f}"

                lines.append(
                    f"{r:>5} {int(np.mean(val('n_ok', 0))):>5} | {pair('median_angle', 0):>32} | "
                    f"{pair('median_angle', P - 1):>32} | {mean2('frac_lt_0p1'):>13} | "
                    f"{mean2('frac_improved'):>13} | "
                    f"{np.mean(val('mean_maha2', 0)):6.3g}/{np.mean(val('mean_maha2', P - 1)):<6.3g} | "
                    f"{np.mean(val('rms_z', P - 1)):.4f}/{np.mean(val('rms_perp', P - 1)):.4f}"
                )
            lines.append("  (<0.1deg, err<r, Mah^2: seed means, one-shot/iterated)")
            lines.append(f"  median error per pass (mean over seeds), by radius:")
            for p in range(P):
                lines.append(
                    f"    pass {p + 1}: "
                    + " ".join(
                        f"{np.mean([mt['all'][v][models[i]][r][p]['median_angle'] for i in idx]):.4f}"
                        for r in radii
                    )
                )
            lines.append("  median error by r_perp tercile (seed mean, one-shot/iterated):")
            for g in ("rperp_tercile1", "rperp_tercile2", "rperp_tercile3"):
                row = []
                for r in radii:
                    a = np.mean([mt[g][v][models[i]][r][0]["median_angle"] for i in idx])
                    b = np.mean([mt[g][v][models[i]][r][P - 1]["median_angle"] for i in idx])
                    row.append(f"{a:.3f}/{b:.3f}")
                lines.append(f"    {g[-8:]:<9}" + " ".join(f"{x:>11}" for x in row))
            lines.append("")
    lines.append(
        f"terciles of r_perp (um): edges {np.round(summ['tercile_edges_um'], 1).tolist()}, ranges "
        f"{[[round(x) for x in rg] for rg in summ['tercile_ranges_um']]}, voxels "
        f"{summ['tercile_n_voxels']}"
    )
    return "\n".join(lines)


def plot_summary(raw: Dict[str, np.ndarray], summ: Dict[str, Any], path: Path) -> None:
    """Median error vs r (log-log): panels = window variant, lines = model type (seed mean, faint
    lines = individual seeds), one-shot solid, iterated dashed, y = x dotted."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    radii = np.asarray(raw["radii"], dtype=float)
    keys = [f"{r:g}" for r in radii]
    models = [str(m) for m in raw["models"]]
    P = raw["err_angle"].shape[-1]
    colors = {"clean": "#1f77b4", "realistic": "#d95f02"}
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.6), sharey=True)
    for ax, v in zip(axes, VARIANTS):
        ax.plot(radii, radii, ":", color="0.4", label="y = x (no correction)")
        for tname, idx in seed_groups(models).items():
            for p, ls in ((0, "-"), (P - 1, "--")):
                ys = np.array(
                    [
                        [summ["metrics"]["all"][v][models[i]][k][p]["median_angle"] for k in keys]
                        for i in idx
                    ]
                )
                for y in ys:
                    ax.plot(radii, y, ls, color=colors.get(tname, "#555555"), alpha=0.25, lw=1)
                ax.plot(
                    radii,
                    ys.mean(0),
                    ls,
                    color=colors.get(tname, "#555555"),
                    lw=2,
                    marker="o",
                    ms=4,
                    label=f"{MODEL_TYPES.get(tname, tname)}, "
                    + ("one-shot" if p == 0 else f"iterated x{P}"),
                )
        ax.axvline(1.0, color="0.8", lw=1)
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlabel("perturbation radius r (deg)")
        ax.set_title("clean windows" if v == "clean" else "realistic windows (all)")
        ax.grid(True, which="both", alpha=0.25)
    axes[0].set_ylabel("median angular error (deg)")
    axes[0].legend(fontsize=7.5, loc="upper left")
    n_v, n_d = raw["fail_pass1"].shape[0], raw["fail_pass1"].shape[2]
    fig.suptitle(
        f"Perturbation sweep: {n_v} voxels x {n_d} directions per radius; faint = single seeds; "
        "training prior ends at r = 1 deg (grey line)",
        fontsize=9,
    )
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)


def load_raw(out_dir: Path) -> Dict[str, np.ndarray]:
    return dict(np.load(out_dir / "perturbation_sweep_raw.npz", allow_pickle=False))


def do_summarize(out_dir: Path) -> None:
    raw = load_raw(out_dir)
    summ = summarize(raw)
    (out_dir / "perturbation_sweep_summary.json").write_text(json.dumps(summ, indent=1))
    text = format_summary(raw, summ)
    (out_dir / "perturbation_sweep_summary.txt").write_text(text + "\n")
    plot_summary(raw, summ, out_dir / "perturbation_sweep.png")
    print(text)


# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------


def select_sweep_voxels(args: argparse.Namespace) -> Dict[str, np.ndarray]:
    """Draw the sweep's voxels (and record their properties) in one process."""
    from generate_toy_orientation_dataset import build_problem, example_dir_for, setup_example
    from icenine.orientation_eval import BatchedObserver, WindowSpec

    example_dir = example_dir_for(args.example)
    setup = setup_example(example_dir, max_q=args.max_q)
    mic = setup[0]
    pos = np.array([v.position for v in mic.voxels], dtype=float)
    r_perp = np.hypot(pos[:, 0], pos[:, 1]) * 1e3
    orients = np.array([v.orientation for v in mic.voxels], dtype=float)
    excluded = excluded_voxel_indices(args.exclude)
    cand = eligible_candidates(r_perp, excluded, args.r_max_um)
    rejected: List[int] = []

    def usable(idx: int) -> bool:
        try:
            pr = build_problem(
                example_dir,
                idx,
                max_q=args.max_q,
                detectors=args.detectors,
                min_sin_eta=args.min_sin_eta,
                setup=setup,
            )
            if len(pr["roi_list"]) < args.min_peaks:
                raise RuntimeError
            obs = BatchedObserver(
                pr["R_nom"],
                pr["vertices"],
                pr["sample"],
                pr["detector_list"],
                pr["range_map"],
                pr["exp_setup"],
                pr["roi_list"],
            )
            WindowSpec.from_nominal(obs, args.window_size, args.frame_half_width)
        except (RuntimeError, ValueError):
            rejected.append(idx)
            return False
        return True

    chosen = sample_voxels(cand, args.n_voxels, args.voxel_seed, usable)
    grain = grain_ids(orients)
    side_um = float(mic.voxels[0].side_length) * 1e3
    nb = near_boundary(pos[:, :2] * 1e3, grain, 1.25 * side_um)
    # the generator's distractor neighbours: 2 nearest voxels within 30 um (not closer than 1 um)
    d_um = np.hypot(*(pos[:, None, :2] - pos[chosen][None, :, :2]).transpose(2, 0, 1)) * 1e3
    cross = np.zeros(len(chosen), dtype=np.int64)
    for k, idx in enumerate(chosen):
        order = np.argsort(d_um[:, k])
        near = [int(i) for i in order if 1.0 < d_um[i, k] <= args.neighbor_radius_um][
            : args.neighbors
        ]
        cross[k] = int(sum(grain[i] != grain[idx] for i in near))
    train_grains = set(grain[excluded].tolist())
    ch = np.array(chosen, dtype=np.int64)
    n_pk = []
    for idx in chosen:
        pr = build_problem(
            example_dir,
            idx,
            max_q=args.max_q,
            detectors=args.detectors,
            min_sin_eta=args.min_sin_eta,
            setup=setup,
        )
        n_pk.append(len(pr["roi_list"]))
    print(
        f"{len(cand)} eligible voxels ({len(excluded)} excluded, r_perp <= {args.r_max_um:g} um); "
        f"drew {len(chosen)}, {len(rejected)} rejected as unusable on the way"
    )
    return dict(
        voxel_indices=ch,
        voxel_position_mm=pos[ch],
        voxel_r_perp_um=r_perp[ch],
        voxel_grain_id=grain[ch],
        voxel_near_boundary=nb[ch],
        voxel_n_cross_grain_neighbours=cross,
        voxel_grain_in_dataset=np.array([grain[i] in train_grains for i in chosen]),
        voxel_n_roi_true=np.array(n_pk),
        excluded_indices=excluded,
        n_candidates=np.array(len(cand)),
        n_rejected=np.array(len(rejected)),
    )


def do_run(args: argparse.Namespace) -> None:
    import multiprocessing as mp

    out_dir = Path(args.out_dir).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    cache = Path(args.cache_dir).resolve()
    cache.mkdir(parents=True, exist_ok=True)
    args.models = [
        (m.split("=", 1)[0], str(Path(m.split("=", 1)[1]).resolve())) for m in args.models
    ]
    args.exclude = [str(Path(p).resolve()) for p in args.exclude]
    t0 = time.time()
    info = select_sweep_voxels(args)
    chosen = [int(i) for i in info["voxel_indices"]]
    print("voxels:", chosen, flush=True)
    items = [(v, k, str(cache / f"voxel_{v}.npz")) for k, v in enumerate(chosen)]
    todo = [it for it in items if not Path(it[2]).exists()]
    print(f"{len(todo)} voxels to run ({len(items) - len(todo)} cached), {args.workers} workers")
    ctx = mp.get_context("spawn")
    # worker setup chdirs into the example, so every path handed to workers is absolute
    with ctx.Pool(args.workers, initializer=init_worker, initargs=(vars(args),)) as pool:
        for vidx, secs in pool.imap_unordered(work, todo):
            print(f"  voxel {vidx} done in {secs:.0f}s ({time.time() - t0:.0f}s total)", flush=True)
    parts = [np.load(p) for _, _, p in items]
    raw: Dict[str, np.ndarray] = {k: np.stack([p[k] for p in parts]) for k in parts[0].files}
    raw.update(info)
    raw["radii"] = np.array(args.radii)
    raw["models"] = np.array([n for n, _ in args.models])
    raw["variants"] = np.array(VARIANTS)
    raw["config_json"] = np.array(json.dumps({k: v for k, v in vars(args).items()}, default=str))
    np.savez_compressed(out_dir / "perturbation_sweep_raw.npz", **raw)
    print(f"saved raw results; wall time {time.time() - t0:.0f}s")
    do_summarize(out_dir)


def build_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    r = sub.add_parser("run", help="run the sweep and summarize")
    r.add_argument("--models", nargs="+", required=True, metavar="NAME=PATH")
    r.add_argument(
        "--exclude",
        nargs="+",
        default=[
            "scripts/toy_orientation_stage3_multi_test.pt",
            "scripts/toy_orientation_arch_dis_test.pt",
        ],
        help="multi-voxel dataset files whose voxels are excluded from the draw",
    )
    r.add_argument("--out-dir", default="benchmarks/toy_orientation_sweep")
    r.add_argument("--cache-dir", default="scripts/perturbation_sweep_cache")
    r.add_argument("--example", default="manygrains")
    r.add_argument("--n-voxels", type=int, default=50)
    r.add_argument("--voxel-seed", type=int, default=0)
    r.add_argument("--sweep-seed", type=int, default=0, help="perturbation directions / draws")
    r.add_argument("--realism-seed", type=int, default=12345)
    r.add_argument("--radii", type=float, nargs="+", default=RADII_DEG)
    r.add_argument("--n-dirs", type=int, default=20)
    r.add_argument("--passes", type=int, default=3)
    r.add_argument("--max-q", type=float, default=8.0)
    r.add_argument("--detectors", default="all")
    r.add_argument("--min-sin-eta", type=float, default=0.3)
    r.add_argument("--min-peaks", type=int, default=40)
    r.add_argument("--r-max-um", type=float, default=500.0)
    r.add_argument("--window-size", type=int, default=32)
    r.add_argument("--frame-half-width", type=int, default=4)
    r.add_argument("--neighbors", type=int, default=2)
    r.add_argument("--neighbor-radius-um", type=float, default=30.0)
    r.add_argument("--neighbor-p", type=float, default=0.5)
    r.add_argument("--neighbor-sigma-deg", type=float, default=0.3)
    r.add_argument(
        "--min-present",
        type=int,
        default=20,
        help="a case needs this many peaks with a lit pixel in their window for the net to run",
    )
    r.add_argument("--workers", type=int, default=8)
    s = sub.add_parser("summarize", help="re-make the summary and plot from the raw npz")
    s.add_argument("--out-dir", default="benchmarks/toy_orientation_sweep")
    return ap


def main() -> None:
    args = build_parser().parse_args()
    if args.cmd == "summarize":
        do_summarize(Path(args.out_dir).resolve())
    else:
        do_run(args)


if __name__ == "__main__":
    main()
