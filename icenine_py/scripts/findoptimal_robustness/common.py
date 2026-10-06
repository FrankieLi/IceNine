"""Shared infrastructure of the FindOptimal robustness study: the voxel set, the per-case detector
images (built exactly as scripts/findoptimal_sweep.py experiment A builds them), the worker setup
and the run-cache helpers."""

import argparse
import contextlib
import io
import json
import os
import sys
import time
from pathlib import Path
from types import SimpleNamespace
from typing import Any, Dict, List, Optional, Sequence, Tuple

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")

import numpy as np
import torch

HERE = Path(__file__).resolve().parent
ICENINE_PY = HERE.parents[1]
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(ICENINE_PY / "scripts"))
sys.path.insert(0, str(ICENINE_PY / "benchmarks"))

import findoptimal_sweep as fs  # noqa: E402
import optimizer_sweep as osw  # noqa: E402
import perturbation_sweep as ps  # noqa: E402

VARIANTS = ["clean", "all"]  # "all" = realistic
OUT_DIR = ICENINE_PY / "benchmarks" / "findoptimal_robustness"
CACHE_DIR = ICENINE_PY / "scripts" / "findoptimal_robustness" / "cache"  # gitignored
SWEEP_RAW = ICENINE_PY / "benchmarks" / "toy_orientation_sweep" / "perturbation_sweep_raw.npz"
N_VOXELS = 200
N_SEEDS = 3
WRONG_DEG = 1.0

_W: Optional[SimpleNamespace] = None


def sweep_config() -> Tuple[Dict[str, Any], Dict[str, np.ndarray]]:
    cfg, sweep = osw.load_sweep_config(SWEEP_RAW)
    return dict(cfg), sweep


def worker_args() -> Dict[str, Any]:
    cfg, _ = sweep_config()
    cfg.update(methods=[], opt_dirs=0, mc_steps=0, mc_restarts=0, mc_step_frac=0)
    return cfg


def select_voxels(n: int) -> Dict[str, np.ndarray]:
    """The study's voxels: the perturbation sweep's selection (same seed, same exclusion of the NN
    dataset voxels) extended to n voxels. Its first 50 are the sweep's 50 voxels (asserted)."""
    cfg, sweep = sweep_config()
    a = argparse.Namespace(**cfg)
    a.n_voxels = n
    a.exclude = [str(Path(p).resolve()) for p in cfg["exclude"]]
    info = ps.select_sweep_voxels(a)
    assert (info["voxel_indices"][:50] == sweep["voxel_indices"]).all()
    return info


def init_worker(args_dict: Dict[str, Any]) -> None:
    """Physics + AdaptiveVoxelReconstructor of findoptimal_sweep (ReconstructQ8.config settings)."""
    global _W
    fs.init_worker(args_dict)
    _W = fs._W


def get_worker() -> SimpleNamespace:
    assert _W is not None
    return _W


def build_case(
    vidx: int, vpos: int, variant: str
) -> Optional[Tuple[np.ndarray, np.ndarray, SimpleNamespace]]:
    """Pixel keys of the detector images of (voxel, variant), the true orientation and the voxel
    context; None if the case cannot be built (the nominal at the truth is not renderable).
    Mirrors findoptimal_sweep.case_images for radius index 0, direction 0 without the sweep
    alignment assertions (the first 50 voxels are checked against the stored experiment A)."""
    from typing import Any as _A

    ctx = get_worker().ctx
    a = ctx.args
    vctx = osw.voxel_context(ctx, vidx)
    D = a.n_dirs
    r = a.radii[0]
    sigma_comp = a.neighbor_sigma_deg / np.sqrt(3.0)
    delta0, draws = ps.case_draws(
        a.sweep_seed, vidx, 0, r, D, len(vctx.sources), sigma_comp, a.neighbor_p
    )
    R_nom0 = ps.perturbed_nominal(vctx.R_true, delta0)
    prep1: List[Optional[_A]] = []
    for j in range(D):
        p, _ = ps.prepare_nominal(ctx, vidx, R_nom0[j])
        prep1.append(p)
    if prep1[0] is None:
        return None
    vi = VARIANTS.index(variant)
    seed = a.realism_seed + 1000003 * vpos + 0
    layers: Dict[str, torch.Tensor] = {}
    b1 = ps.render_batch(prep1, delta0, draws, vctx.sources, variant, seed, a, layers=layers)
    if b1["n_present"].numpy()[0] < a.min_present:
        return None
    edit = None
    if variant != "clean":
        p = prep1[0]
        n = p.n
        edit = osw.realism_edit(
            layers["clean"][0, :n].numpy(),
            layers["dis"][0, :n].numpy(),
            b1["windows"][0, :n].numpy(),
            p.spec,
            p.obs.det_idx.numpy(),
            a.frame_half_width,
            ctx.geo,
        )
    keys = osw.case_image_keys(ctx, vctx, variant, draws[0], edit)
    del vi
    return keys, vctx.R_true, vctx


def attach(keys: np.ndarray) -> None:
    fs.attach_images(keys)


def run_seed(vpos: int, seed: int) -> np.random.Generator:
    """Seed 0 reproduces experiment A of findoptimal_sweep (rng 10_000 + voxel position)."""
    return np.random.default_rng(10_000 + vpos if seed == 0 else [10_000 + vpos, seed])


def err_deg(R: np.ndarray, R_true: np.ndarray) -> np.ndarray:
    from csl import reduced_misorientation_deg

    return reduced_misorientation_deg(R, R_true)


class Recorder:
    """Collects the reconstructor's recorder events into flat arrays."""

    def __init__(self) -> None:
        self.events: List[Tuple[str, Dict[str, Any]]] = []

    def __call__(self, name: str, data: Dict[str, Any]) -> None:
        self.events.append((name, data))

    def to_arrays(self) -> Dict[str, np.ndarray]:
        out: Dict[str, np.ndarray] = {}
        fo: Dict[str, List[Any]] = {"R_in": [], "R_out": [], "cost": [], "index": []}
        for name, d in self.events:
            if name in ("discrete", "quick_mc"):
                tag = f"L{d['level']}_{'disc' if name == 'discrete' else 'qmc'}"
                for k, v in d.items():
                    if k == "level":
                        continue
                    out[f"{tag}_{k}"] = np.asarray(v)
            elif name == "find_candidate":
                for k in fo:
                    fo[k].append(d[k])
            elif name in ("variance", "final"):
                out[f"{name}_R"] = np.asarray(d["R"])
                out[f"{name}_cost"] = np.asarray(d["cost"])
        for k, v in fo.items():
            out[f"find_{k}"] = np.asarray(v)
        return out


def quiet() -> Any:
    return contextlib.redirect_stdout(io.StringIO())


def run_pool(func: Any, todo: List[Any], workers: int, label: str) -> None:
    import multiprocessing as mp

    t0 = time.time()
    mpctx = mp.get_context("spawn")
    with mpctx.Pool(workers, initializer=init_worker, initargs=(worker_args(),)) as pool:
        for k, res in enumerate(pool.imap_unordered(func, todo), 1):
            print(f"  [{label} {k}/{len(todo)}] {res} ({time.time() - t0:.0f}s total)", flush=True)


def save_json(path: Path, obj: Any) -> None:
    path.write_text(json.dumps(obj, indent=1, default=float))


def wilson(k: int, n: int, z: float = 1.96) -> Tuple[float, float, float]:
    """(rate, lo, hi) Wilson 95% interval."""
    if n == 0:
        return float("nan"), float("nan"), float("nan")
    p = k / n
    den = 1 + z * z / n
    c = (p + z * z / (2 * n)) / den
    h = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / den
    return p, max(0.0, c - h), min(1.0, c + h)
