#!/usr/bin/env python3
"""
Generate the windowed local-orientation-refinement datasets for the toy NN.

Picks one voxel from Example2.ThreeVoxels' ground-truth .mic, defines its ROI
peak set at the ground-truth orientation, then renders thresholded windows
around the ROI peaks for many small orientation offsets.

  train: offsets drawn from the prior (uniform in a ball of radius --prior-radius)
  test:  --test-per-bin offsets at each fixed magnitude in --test-magnitudes,
         random directions (the bench_hp_sweep-style protocol)

Windows are stored as uint8 lit / not-lit values (thresholded), so the network
cannot read orientation information from exact simulated intensities.

Usage:
  cd icenine_py
  uv run python scripts/generate_toy_orientation_dataset.py --smoke-test
  uv run python scripts/generate_toy_orientation_dataset.py --n-train 1500 --test-per-bin 30
"""

import argparse
import os
import time
from pathlib import Path

import numpy as np
import torch

project_root = Path(__file__).parent.parent.parent
DEFAULT_EXAMPLE = project_root / "Examples" / "Example2.ThreeVoxels"
EXAMPLES = {
    "threevoxels": DEFAULT_EXAMPLE,
    "manygrains": project_root / "Examples" / "Example2.ManyGrains",
}


def example_dir_for(name: str) -> Path:
    """Example directory for a dataset's meta["example"] (default: ThreeVoxels)."""
    return EXAMPLES[name or "threevoxels"]


def setup_example(example_dir: Path, basename: str = "3Grains.sim", max_q=None):
    """Load physics + ground-truth mic. Same pattern as benchmarks/bench_riemannian_optimization.py.

    max_q (1/Angstrom) overrides the config's MaxQ, which limits the reflection list
    (Example2's config has 16; reconstruction typically uses 8). Note: changes the
    process's working directory to example_dir.
    """
    from icenine.config_file import ConfigFile
    from icenine.experiment_setup import XDMExperimentSetup
    from icenine.mic_file import MicFile
    from icenine.reconstructor import _get_voxel_vertices
    from icenine.sample import Sample
    from icenine.simulation import Simulation

    config_path = example_dir / "ConfigFiles" / "Example2.Simulation.config"
    os.chdir(example_dir)

    config = ConfigFile.from_file(str(config_path))
    config.out_file_basename = basename
    if max_q is not None:
        config.max_q = float(max_q)

    exp_setup = XDMExperimentSetup(config)
    exp_setup.initialize_experiment()
    detector_list = exp_setup.get_detector_list()
    range_map = exp_setup.get_range_to_index_map()

    sample = Sample()
    exp_setup.initialize_sample(sample, detector_list[0])
    simulator = Simulation(exp_setup)
    structure_list = sample.get_structure_list()

    mic_path = Path(config.sample_filename)
    if not mic_path.is_absolute():
        mic_path = example_dir / mic_path
    mic = MicFile.read(str(mic_path))

    return (
        mic,
        sample,
        detector_list,
        range_map,
        exp_setup,
        simulator,
        structure_list,
        _get_voxel_vertices,
    )


def build_problem(
    example_dir: Path,
    voxel_index=None,
    max_q=None,
    detectors: str = "first",
    min_sin_eta: float = 0.0,
    setup=None,
    orientation=None,
):
    """Set up physics and pick a voxel with a non-empty ROI set.

    setup: the tuple returned by setup_example(), to reuse one physics setup for many voxels.

    max_q and detectors are passed to setup_example / define_roi_set; min_sin_eta
    drops near-axis spots (|sin eta| below it at the nominal orientation, decision D3).
    The defaults reproduce Stage 0 (config MaxQ, first detector only, no filter).

    orientation: replace the voxel's orientation (3x3) while keeping its position and vertices,
    e.g. a twin of the voxel (used for distractor sources).

    Returns a dict with everything both this script and the Bayes baseline need.
    """
    from icenine.orientation_nn import define_roi_set

    if setup is None:
        setup = setup_example(example_dir, max_q=max_q)
    mic, sample, detector_list, range_map, exp_setup, simulator, structure_list, get_vertices = (
        setup
    )
    candidates = [voxel_index] if voxel_index is not None else range(len(mic.voxels))
    for idx in candidates:
        voxel = mic.voxels[idx]
        vertices = get_vertices(voxel)
        R_nom = (voxel.orientation if orientation is None else orientation).astype(np.float64)
        roi_list = define_roi_set(
            torch.from_numpy(R_nom).float(),
            vertices,
            sample,
            detector_list,
            range_map,
            exp_setup,
            structure_list,
            simulator,
            phase_index=voxel.phase,
            detectors=detectors,
        )
        if roi_list and min_sin_eta > 0.0:
            from icenine.orientation_eval import BatchedObserver

            obs = BatchedObserver(
                R_nom, vertices, sample, detector_list, range_map, exp_setup, roi_list
            )
            se = obs.sin_eta(torch.zeros(1, 3, dtype=torch.float64))[0].numpy()
            roi_list = [p for p, s in zip(roi_list, se) if s >= min_sin_eta]
        if roi_list:
            return dict(
                voxel_index=idx,
                voxel=voxel,
                vertices=vertices,
                R_nom=R_nom,
                roi_list=roi_list,
                sample=sample,
                detector_list=detector_list,
                range_map=range_map,
                exp_setup=exp_setup,
                simulator=simulator,
            )
    raise RuntimeError("No voxel in the .mic produced a non-empty ROI set")


def render_dataset(problem, offsets_deg: np.ndarray, window_size: int, label: str) -> torch.Tensor:
    """Render thresholded uint8 windows (N, n_peaks, W, W) for the given offsets."""
    from icenine.orientation_eval import offsets_to_matrices
    from icenine.orientation_nn import render_local_windows

    mats = offsets_to_matrices(offsets_deg, problem["R_nom"])
    n_peaks = len(problem["roi_list"])
    windows = torch.zeros(len(offsets_deg), n_peaks, window_size, window_size, dtype=torch.uint8)
    t0 = time.time()
    for i, mat in enumerate(mats):
        w, _missing = render_local_windows(
            torch.from_numpy(mat).float(),
            problem["roi_list"],
            problem["vertices"],
            problem["sample"],
            problem["detector_list"],
            problem["range_map"],
            problem["exp_setup"],
            problem["simulator"],
            window_size=window_size,
        )
        windows[i] = (w > 0).to(torch.uint8)
        if (i + 1) % 100 == 0 or (i + 1) == len(mats):
            print(f"  [{label}] {i + 1}/{len(mats)} ({time.time() - t0:.0f}s)")
    return windows


def render_dataset_observer(
    problem, offsets_deg: np.ndarray, window_size: int, frame_half_width: int
):
    """Frame-coded exact thresholded windows from the batched observer (Stage 1)."""
    from icenine.orientation_eval import BatchedObserver, WindowSpec, render_windows

    obs = BatchedObserver(
        problem["R_nom"],
        problem["vertices"],
        problem["sample"],
        problem["detector_list"],
        problem["range_map"],
        problem["exp_setup"],
        problem["roi_list"],
    )
    spec = WindowSpec.from_nominal(obs, window_size, frame_half_width)
    windows, status = render_windows(obs, spec, offsets_deg)
    return windows, status, obs.peak_context()


def select_voxels(mic, n_voxels: int, r_max_um: float, seed: int):
    """Candidate voxel indices spanning r_perp in [0, r_max_um]: for each of n_voxels evenly
    spaced target radii, up to 60 random voxels within a small tolerance, in random order.
    Returned in order of increasing target radius. The caller accepts the first usable
    candidate per radius and skips candidates whose orientation (grain) was already accepted,
    so no grain is used twice (see main_multi)."""
    rng = np.random.default_rng(seed)
    pos = np.array([v.position for v in mic.voxels], dtype=float)
    r_um = np.hypot(pos[:, 0], pos[:, 1]) * 1e3
    targets = np.linspace(0.0, r_max_um, n_voxels)
    spacing = targets[1] - targets[0] if n_voxels > 1 else r_max_um
    tol = max(5.0, 0.5 * spacing)
    out = []
    for t in targets:
        cand = np.nonzero(np.abs(r_um - t) <= tol)[0]
        cand = cand[rng.permutation(len(cand))]
        out.append([int(c) for c in cand[:60]])
    return targets, out


def accept_voxels(cand_lists, orientation_of, usable, tol: float = 1e-6):
    """Accept at most one voxel per candidate list (one per target radius).

    Candidates whose orientation matches an already ACCEPTED voxel (same grain) are skipped;
    ``usable(idx)`` returns a truthy result for an acceptable voxel and None/False otherwise.
    Returns the list of ``usable`` results for the accepted voxels, in order.
    """
    accepted, used = [], []
    for cands in cand_lists:
        for idx in cands:
            R = orientation_of(idx)
            if any(np.abs(R - u).max() <= tol for u in used):
                continue
            res = usable(idx)
            if res:
                accepted.append(res)
                used.append(R)
                break
    return accepted


def sigma3_matrix() -> np.ndarray:
    """Cubic twin operator: 60 degrees about [111] in the crystal frame."""
    from scipy.spatial.transform import Rotation

    return Rotation.from_rotvec(np.radians(60.0) * np.ones(3) / np.sqrt(3.0)).as_matrix()


def build_distractor_sources(
    example_dir, mic, target_index, setup, args, max_q, detectors, min_sin_eta
):
    """Observers for the spots of the target's neighbours, for distractor rendering.

    Up to --neighbors mic voxels nearest to the target (in the sample plane, excluding ones
    closer than 1 um, within --neighbor-radius-um), each with its own orientation, vertices
    and ROI set; plus, with --twin, a Sigma3 twin (orientation R @ T) at the nearest
    neighbour's position. Returns (observers, descriptions).
    """
    from icenine.orientation_eval import BatchedObserver

    pos = np.array([v.position for v in mic.voxels], dtype=float)
    d_um = np.hypot(*(pos[:, :2] - pos[target_index, :2]).T) * 1e3
    order = np.argsort(d_um)
    near = [int(i) for i in order if 1.0 < d_um[i] <= args.neighbor_radius_um][: args.neighbors]
    specs = [(i, None, f"voxel {i} ({d_um[i]:.0f} um)") for i in near]
    if args.twin and near:
        T = sigma3_matrix()
        specs.append((near[0], mic.voxels[near[0]].orientation @ T, f"twin of voxel {near[0]}"))
    observers, desc = [], []
    for idx, orient, name in specs:
        try:
            pr = build_problem(
                example_dir,
                idx,
                max_q=max_q,
                detectors=detectors,
                min_sin_eta=min_sin_eta,
                setup=setup,
                orientation=orient,
            )
        except RuntimeError:
            continue
        observers.append(
            BatchedObserver(
                pr["R_nom"],
                pr["vertices"],
                pr["sample"],
                pr["detector_list"],
                pr["range_map"],
                pr["exp_setup"],
                pr["roi_list"],
            )
        )
        desc.append(name)
    return observers, desc


def main_multi(args, outdir_abs):
    """Multi-voxel dataset (Stage 3 step 3): padded windows, per-voxel context table."""
    from icenine.orientation_eval import (
        BatchedObserver,
        WindowSpec,
        render_windows,
        sample_fixed_magnitude_offsets,
        sample_prior_offsets,
    )

    example_dir = example_dir_for(args.example)
    setup = setup_example(example_dir, max_q=args.max_q)
    mic = setup[0]
    targets, cand_lists = select_voxels(mic, args.n_voxels, args.r_max_um, args.voxel_seed)
    rng = np.random.default_rng(args.seed)

    def usable(idx):
        try:
            pr = build_problem(
                example_dir,
                idx,
                max_q=args.max_q,
                detectors=args.detectors,
                min_sin_eta=args.min_sin_eta,
                setup=setup,
            )
        except RuntimeError:
            return None
        return pr if len(pr["roi_list"]) >= args.min_peaks else None

    problems = accept_voxels(cand_lists, lambda i: mic.voxels[i].orientation, usable)
    V = len(problems)
    r_perp = np.array(
        [np.hypot(*np.asarray(p["voxel"].position, dtype=float)[:2]) * 1e3 for p in problems]
    )
    n_peaks = np.array([len(p["roi_list"]) for p in problems])
    order = np.argsort(r_perp)
    problems = [problems[i] for i in order]
    r_perp, n_peaks = r_perp[order], n_peaks[order]
    held_out = np.zeros(V, dtype=bool)
    held_out[args.holdout_offset :: args.holdout_every] = True
    M = int(n_peaks.max())
    print(
        f"{V} voxels, r_perp {r_perp.min():.0f}-{r_perp.max():.0f} um, peaks {n_peaks.min()}-{M}, "
        f"held out {int(held_out.sum())}"
    )
    W = args.window_size
    context = None  # (V, M, D), D = 14 + number of detectors, allocated on first voxel
    R_nom = torch.zeros(V, 3, 3)
    train_w, train_off, train_vid = [], [], []
    test_w, test_off, test_vid, test_mag = [], [], [], []
    train_dis, test_dis = [], []
    dis_rng = np.random.default_rng(args.seed + 1000)
    sigma_comp = args.neighbor_sigma_deg / np.sqrt(3.0)
    for v, pr in enumerate(problems):
        obs = BatchedObserver(
            pr["R_nom"],
            pr["vertices"],
            pr["sample"],
            pr["detector_list"],
            pr["range_map"],
            pr["exp_setup"],
            pr["roi_list"],
        )
        spec = WindowSpec.from_nominal(obs, W, args.frame_half_width)
        ctx = obs.peak_context()
        if context is None:
            context = torch.zeros(V, M, ctx.shape[-1])
        context[v, : n_peaks[v]] = ctx
        R_nom[v] = torch.from_numpy(pr["R_nom"])
        n_tr = 0 if held_out[v] else args.per_voxel_train
        sources, src_desc = [], []
        if args.neighbors or args.twin:
            sources, src_desc = build_distractor_sources(
                example_dir,
                mic,
                pr["voxel_index"],
                setup,
                args,
                args.max_q,
                args.detectors,
                args.min_sin_eta,
            )

        def dis_layer(offsets):
            """Distractor layer for these target offsets: each source is perturbed by the
            target's offset plus its own random misorientation (sigma --neighbor-sigma-deg)."""
            if not sources:
                return None
            from icenine.orientation_eval import render_distractor_windows

            dsrc = [offsets + dis_rng.normal(0.0, sigma_comp, size=offsets.shape) for _ in sources]
            active = [dis_rng.random(len(offsets)) < args.neighbor_p for _ in sources]
            d = render_distractor_windows(obs, spec, sources, dsrc, source_active=active)
            pad_d = torch.zeros(len(offsets), M, W, W, dtype=torch.uint8)
            pad_d[:, : n_peaks[v]] = d
            return pad_d

        offs = [sample_prior_offsets(n_tr, args.prior_radius, rng)] if n_tr else []
        if n_tr:
            w, _ = render_windows(obs, spec, offs[0])
            pad = torch.zeros(n_tr, M, W, W, dtype=torch.uint8)
            pad[:, : n_peaks[v]] = w
            dl = dis_layer(offs[0])
            if dl is not None:
                train_dis.append(dl)
            train_w.append(pad)
            train_off.append(torch.from_numpy(offs[0]).float())
            train_vid.append(torch.full((n_tr,), v, dtype=torch.long))
        for mag in args.test_magnitudes:
            o = sample_fixed_magnitude_offsets(args.per_voxel_test, mag, rng)
            w, _ = render_windows(obs, spec, o)
            pad = torch.zeros(len(o), M, W, W, dtype=torch.uint8)
            pad[:, : n_peaks[v]] = w
            dl = dis_layer(o)
            if dl is not None:
                test_dis.append(dl)
            test_w.append(pad)
            test_off.append(torch.from_numpy(o).float())
            test_vid.append(torch.full((len(o),), v, dtype=torch.long))
            test_mag.append(torch.full((len(o),), mag))
        print(
            f"  voxel {v:2d} (mic {pr['voxel_index']:5d}) r_perp {r_perp[v]:5.0f} um  "
            f"{n_peaks[v]:3d} peaks  {'held out' if held_out[v] else 'train'}"
            + (f"  sources: {', '.join(src_desc)}" if src_desc else ""),
            flush=True,
        )
    meta = dict(
        window_size=W,
        n_peaks=M,
        n_peaks_per_voxel=torch.from_numpy(n_peaks),
        voxel_indices=torch.tensor([p["voxel_index"] for p in problems]),
        r_perp_um=torch.from_numpy(r_perp).float(),
        held_out=torch.from_numpy(held_out),
        R_nom=R_nom,
        context=context,
        prior_radius_deg=args.prior_radius,
        max_q=args.max_q if args.max_q is not None else float("nan"),
        detectors=args.detectors,
        min_sin_eta=args.min_sin_eta,
        renderer="observer",
        frame_half_width=args.frame_half_width,
        seed=args.seed,
        example=args.example,
        multi_voxel=True,
    )
    for name, w, off, vid, extra, dis in (
        ("train", train_w, train_off, train_vid, {}, train_dis),
        ("test", test_w, test_off, test_vid, {"magnitudes_deg": torch.cat(test_mag)}, test_dis),
    ):
        windows = torch.cat(w)
        if dis:
            extra = {**extra, "dis_windows": torch.cat(dis)}
        path = outdir_abs / f"toy_orientation_{args.tag}_{name}.pt"
        torch.save(
            {
                "windows": windows,
                "offsets_deg": torch.cat(off),
                "voxel_id": torch.cat(vid),
                **meta,
                **extra,
            },
            path,
        )
        print(f"Saved {path}  ({windows.numel() / 1e9:.2f} GB, {len(windows)} samples)")


def main():
    from icenine.orientation_eval import sample_fixed_magnitude_offsets, sample_prior_offsets

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--example", choices=sorted(EXAMPLES), default="threevoxels")
    parser.add_argument(
        "--voxel-index", type=int, default=None, help="Force a specific voxel; default auto-selects"
    )
    parser.add_argument("--n-train", type=int, default=1500)
    parser.add_argument("--test-per-bin", type=int, default=30)
    parser.add_argument("--test-magnitudes", type=float, nargs="+", default=[0.25, 0.5, 1.0, 2.0])
    parser.add_argument(
        "--prior-radius", type=float, default=2.5, help="Training prior: uniform ball radius (deg)"
    )
    parser.add_argument("--window-size", type=int, default=32)
    parser.add_argument(
        "--max-q", type=float, default=None, help="Override the config MaxQ (1/A); Stage 1 uses 8"
    )
    parser.add_argument(
        "--detectors",
        choices=["first", "all"],
        default="first",
        help="ROI peaks on the first detector only (Stage 0) or on every detector (Stage 1)",
    )
    parser.add_argument(
        "--renderer",
        choices=["simulator", "observer"],
        default="simulator",
        help="simulator: per-peak renderer with intensities (Stage 0); observer: exact "
        "thresholded pixels with the frame offset coded in the pixel value (Stage 1)",
    )
    parser.add_argument(
        "--frame-half-width", type=int, default=4, help="observer renderer: frames kept (+-K)"
    )
    parser.add_argument(
        "--min-sin-eta", type=float, default=0.0, help="drop spots with |sin eta| below this (D3)"
    )
    parser.add_argument(
        "--n-voxels",
        type=int,
        default=0,
        help="multi-voxel mode: this many voxels spanning r_perp in [0, --r-max-um]",
    )
    parser.add_argument("--voxel-seed", type=int, default=0)
    parser.add_argument("--r-max-um", type=float, default=500.0)
    parser.add_argument("--min-peaks", type=int, default=40)
    parser.add_argument("--per-voxel-train", type=int, default=500)
    parser.add_argument("--per-voxel-test", type=int, default=10, help="per magnitude, per voxel")
    parser.add_argument("--holdout-every", type=int, default=5, help="every k-th voxel is held out")
    parser.add_argument("--holdout-offset", type=int, default=2)
    parser.add_argument(
        "--neighbors",
        type=int,
        default=0,
        help="multi-voxel mode: render spots of this many nearest mic voxels as a distractor layer",
    )
    parser.add_argument("--neighbor-radius-um", type=float, default=30.0)
    parser.add_argument(
        "--neighbor-p",
        type=float,
        default=0.5,
        help="probability that each distractor source is present in a given sample",
    )
    parser.add_argument(
        "--neighbor-sigma-deg",
        type=float,
        default=0.3,
        help="random misorientation of each neighbour relative to the target's perturbed orientation",
    )
    parser.add_argument(
        "--twin", action="store_true", help="add a Sigma3 twin of the nearest voxel"
    )
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--outdir", default=str(Path(__file__).parent))
    parser.add_argument("--tag", default="stage0")
    parser.add_argument("--smoke-test", action="store_true", help="Tiny run: 60 train, 5 per bin")
    args = parser.parse_args()
    if args.smoke_test:
        args.n_train, args.test_per_bin, args.tag = 60, 5, "smoke"

    outdir_abs = Path(args.outdir).resolve()  # build_problem() changes directory
    if args.n_voxels:
        return main_multi(args, outdir_abs)
    example_dir = example_dir_for(args.example)
    print(f"Setting up {example_dir} ...")
    problem = build_problem(
        example_dir,
        args.voxel_index,
        max_q=args.max_q,
        detectors=args.detectors,
        min_sin_eta=args.min_sin_eta,
    )
    n_peaks = len(problem["roi_list"])
    pos = np.asarray(problem["voxel"].position, dtype=float)
    r_perp_um = float(np.hypot(pos[0], pos[1]) * 1e3)
    print(
        f"  voxel {problem['voxel_index']}: r_perp {r_perp_um:.1f} um, side "
        f"{problem['voxel'].side_length * 1e3:.2f} um, {n_peaks} ROI peaks"
    )

    rng = np.random.default_rng(args.seed)
    train_offsets = sample_prior_offsets(args.n_train, args.prior_radius, rng)
    test_offsets, test_mags = [], []
    for mag in args.test_magnitudes:
        test_offsets.append(sample_fixed_magnitude_offsets(args.test_per_bin, mag, rng))
        test_mags += [mag] * args.test_per_bin
    test_offsets = np.concatenate(test_offsets)

    meta = dict(
        window_size=args.window_size,
        n_peaks=n_peaks,
        voxel_index=problem["voxel_index"],
        R_nom=torch.from_numpy(problem["R_nom"]),
        prior_radius_deg=args.prior_radius,
        max_q=args.max_q if args.max_q is not None else float("nan"),
        detectors=args.detectors,
        min_sin_eta=args.min_sin_eta,
        renderer=args.renderer,
        frame_half_width=args.frame_half_width,
        seed=args.seed,
        example=args.example,
        r_perp_um=r_perp_um,
        side_um=float(problem["voxel"].side_length * 1e3),
    )
    outdir = outdir_abs
    for name, offsets, extra in (
        ("train", train_offsets, {}),
        ("test", test_offsets, {"magnitudes_deg": torch.tensor(test_mags)}),
    ):
        if args.renderer == "observer":
            windows, status, context = render_dataset_observer(
                problem, offsets, args.window_size, args.frame_half_width
            )
            extra = {**extra, "status": status, "context": context}
            s = status.float()
            print(
                f"  [{name}] spots inside window {(s == 0).float().mean():.1%}, absent "
                f"{(s == 1).float().mean():.1%}, frame outside +-K {(s == 2).float().mean():.1%}, "
                f"left window {(s == 3).float().mean():.1%}"
            )
        else:
            windows = render_dataset(problem, offsets, args.window_size, name)
        path = outdir / f"toy_orientation_{args.tag}_{name}.pt"
        torch.save(
            {"windows": windows, "offsets_deg": torch.from_numpy(offsets).float(), **meta, **extra},
            path,
        )
        print(f"Saved {path}  ({windows.numel() / 1e6:.0f} MB)")


if __name__ == "__main__":
    main()
