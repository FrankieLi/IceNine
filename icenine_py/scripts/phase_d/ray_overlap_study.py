#!/usr/bin/env python3
"""Ray-consistency overlap rule, steps 1-2 of the plan (Phase D task ray-overlap).

Step 1 (peak level, no search) and step 2 (coarse stage levels 0-2, no FindOptimal) on full-sample
images (clean and realistic stacks), about 30 voxels, arms: old rule, C++ rule, and the ray rule
x {R1, R2, R4, R8, A} (arm names f0_*; the f1 variant was dropped 2026-10-12).

  uv run python scripts/phase_d/ray_overlap_study.py pilot            # 2 voxels, timed
  uv run python scripts/phase_d/ray_overlap_study.py run --step 1 --workers 8
  uv run python scripts/phase_d/ray_overlap_study.py run --step 2 --workers 8
  uv run python scripts/phase_d/ray_overlap_study.py summary

Results: benchmarks/phase_d_seed_diag/ray_overlap/{peak,coarse}/<variant>_v<voxel>.json
Every wall time here is "contended" (several workers on a shared laptop).
"""

import argparse
import json
import multiprocessing as mp
import os
import sys
import time
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import numpy as np

HERE = Path(__file__).resolve().parent
ICE = HERE.parents[1]
EX = HERE.parents[2] / "Examples" / "Example2.ManyGrains"
OUT = ICE / "benchmarks" / "phase_d_seed_diag" / "ray_overlap"
sys.path.insert(0, str(ICE / "scripts" / "common"))
sys.path.insert(0, str(ICE / "scripts" / "findoptimal_robustness"))
MIC = "SimInput/rand_500grains_1mm_neworient_s0.mic"
BASIN_DEG = 3.0
VARIANTS = tuple(os.environ.get("RAY_VARIANTS", "clean,realistic").split(","))
TOLS = ("R1", "R2", "R4", "R8", "A")
# arm names keep the "f0_" prefix of the earlier results; f1 was dropped 2026-10-12
# centroid task 2026-10-12: every arm projects the centroid; RAY_VERTEX0=1 = the old vertex-0 arm
CENTER_SCREEN = os.environ.get("RAY_VERTEX0", "") == ""
ARM_NAMES = ["old", "cpp"] + [f"f0_{t}" for t in TOLS]
NEAR_DEG = (0.5, 1.0, 2.0, 5.0)
N_NEAR, N_RANDOM = 10, 200
N_CASES = 8
CASE_NAMES = ("a", "b", "c", "d", "e", "f_dark", "f_bright", "det0_bright")


def arm_spec(name: str, registration_mm: float) -> Tuple[bool, Optional[Any]]:
    """(cpp_overlap_rule, RayRule or None) of an arm name."""
    from icenine.ray_overlap import RayRule

    if name == "old":
        return False, None
    if name == "cpp":
        return True, None
    f, tol = name.split("_")
    if f != "f0":
        raise ValueError(f"arm {name!r}: only f0 exists (f1 dropped 2026-10-12)")
    if tol == "A":
        return True, RayRule(arm="A", registration_mm=registration_mm)
    return True, RayRule(arm="R", radius_px=int(tol[1:]), registration_mm=registration_mm)


def select_voxels() -> List[int]:
    """The 6 seed-diag voxels, the 18 random voxels and 6 more drawn with seed 11 (30 in all)."""
    base = Path(__file__).resolve().parents[2] / "benchmarks" / "phase_d_seed_diag"
    v = [int(d["voxel"]) for d in json.loads((base / "voxels.json").read_text())]
    v += [int(d["voxel"]) for d in json.loads((base / "voxels_random.json").read_text())]
    rng = np.random.default_rng(11)
    n_vox = int(np.loadtxt(EX / MIC, skiprows=1).shape[0])
    extra = [int(x) for x in rng.choice(n_vox, 12, replace=False) if int(x) not in v]
    return v + extra[: 30 - len(v)]


class Ctx:
    """Per-process state: setup and reconstructor of the current image variant (init_worker)."""

    variant: str = ""
    setup: Any = None
    rec: Any = None
    truth: Any = None


def init_worker(variant: str) -> None:
    if Ctx.variant == variant:
        return
    import torch

    torch.set_num_threads(1)
    os.chdir(EX)
    from icenine import reconstructor as RC
    from icenine.config_file import ConfigFile
    from icenine.experimental_data import ExperimentalData
    from icenine.mic_file import MicFile

    cfg = ConfigFile.from_file(str(HERE / "configs" / f"ReconstructPhaseD_mc_{variant}.config"))
    stack = EX / "ScatteringData_PhaseD" / "stacks" / f"{variant}.npy"
    exp = ExperimentalData.from_binary_memmap(stack, 180, 2)
    Ctx.variant = variant
    Ctx.setup = RC.setup_reconstruction(cfg, exp_data=exp)
    Ctx.rec = RC.AdaptiveVoxelReconstructor(Ctx.setup)
    # counts evaluations by wrapping evaluate(), which the batched screen bypasses: pin it off
    Ctx.setup.search_params.batched_screen = False
    Ctx.rec.min_sin_eta = 0.0
    Ctx.truth = MicFile.read(MIC)


def orientation_set(R_true: np.ndarray, seed: int) -> Tuple[np.ndarray, List[str]]:
    """Truth, N_NEAR perturbations at each of NEAR_DEG about random axes, N_RANDOM random."""
    from scipy.spatial.transform import Rotation

    rng = np.random.default_rng(seed)
    Rs = [R_true]
    labels = ["truth"]
    for deg in NEAR_DEG:
        for _ in range(N_NEAR):
            ax = rng.normal(size=3)
            ax /= np.linalg.norm(ax)
            Rs.append(Rotation.from_rotvec(ax * np.radians(deg)).as_matrix() @ R_true)
            labels.append(f"near{deg}")
    for R in Rotation.random(N_RANDOM, random_state=int(rng.integers(1 << 30))).as_matrix():
        Rs.append(R)
        labels.append("random")
    return np.stack(Rs), labels


def peak_task(voxel: int) -> Dict[str, Any]:
    """Step 1 for one voxel: every arm x {screen, local} on the orientation set."""
    from icenine import reconstructor as RC
    from icenine.cost_functions import VoxelCostFunction

    setup, rec = Ctx.setup, Ctx.rec
    v = setup.sample.get_mic().voxels[voxel]
    verts = RC._get_voxel_vertices(v)
    R_true = np.asarray(Ctx.truth.voxels[voxel].orientation, dtype=np.float64)
    Rs, labels = orientation_set(R_true, seed=voxel)
    reg = float(setup.exp_setup.config_file.detector_dist_deviation)
    res: Dict[str, Any] = dict(voxel=voxel, labels=labels)
    t0 = time.perf_counter()
    for fn_name, radius, max_q in (("screen", 3, 5.0), ("local", 0, 5.0)):
        for arm in ARM_NAMES:
            cpp, ray = arm_spec(arm, reg)
            cf = VoxelCostFunction(
                simulator=setup.simulator,
                detector_list=setup.detector_list,
                range_map=setup.range_map,
                exp_data=setup.exp_data,
                sample=setup.sample,
                structure_list=setup.structure_list,
                mode="hard",
                eta_limit=setup.exp_setup.get_eta_limit(),
                pixel_radius=radius,
                max_q=max_q,
                min_sin_eta=0.0,
                cpp_overlap_rule=cpp,
                center_screen=CENTER_SCREEN,
                ray_rule=ray,
            )
            rows = []
            for R in Rs:
                info = cf.evaluate(R.astype(np.float32), verts, v.phase)
                rows.append(
                    [
                        info.peak_overlap,
                        info.peak_on_detector,
                        info.pixel_overlap,
                        info.pixel_on_detector,
                        info.quality,
                    ]
                    + list(info.ray_cases if info.ray_cases is not None else [0] * N_CASES)
                )
            res[f"{fn_name}|{arm}"] = rows
    res["wall_s"] = time.perf_counter() - t0
    return res


class StopCoarse(Exception):
    pass


def coarse_task(voxel: int, arms: Optional[List[str]] = None) -> Dict[str, Any]:
    """Step 2 for one voxel: levels 0-2 of the search for every arm, stopped before FindOptimal."""
    from icenine import cost_functions as CF
    from icenine import reconstructor as RC
    from csl import reduced_misorientation_deg

    setup, rec = Ctx.setup, Ctx.rec
    v = setup.sample.get_mic().voxels[voxel]
    verts = RC._get_voxel_vertices(v)
    R_true = np.asarray(Ctx.truth.voxels[voxel].orientation, dtype=np.float64)
    reg = float(setup.exp_setup.config_file.detector_dist_deviation)
    out: Dict[str, Any] = dict(voxel=voxel, arms={})
    sp = setup.search_params
    orig_eval = CF.VoxelCostFunction.evaluate
    counter = [0]

    def counted(self: Any, *a: Any, **k: Any) -> Any:
        counter[0] += 1
        return orig_eval(self, *a, **k)

    CF.VoxelCostFunction.evaluate = counted  # type: ignore[method-assign]
    try:
        for arm in arms or ARM_NAMES:
            cpp, ray = arm_spec(arm, reg)
            sp.cpp_overlap_rule, sp.ray_rule = cpp, ray
            sp.cpp_center_screen = CENTER_SCREEN
            levels: List[Dict[str, Any]] = []
            counter[0] = 0
            cur: Dict[str, Any] = {}

            def drec(ev: str, d: Dict[str, Any]) -> None:
                if ev == "level":
                    cur.clear()
                    cur.update(n_fz=d["n_fz"], n_local=d["n_local"])
                elif ev == "level_end":
                    cur.update(n_cliques=d["n_cliques"])

            def hook(ev: str, d: Dict[str, Any]) -> None:
                if ev == "discrete":
                    R = d["R"].astype(np.float64)
                    err = reduced_misorientation_deg(R, np.broadcast_to(R_true, R.shape).copy())
                    cur.update(
                        level=d["level"],
                        n_cand=len(R),
                        d_min_disc=float(err.min()) if len(R) else None,
                        n_within_disc=int((err < BASIN_DEG).sum()),
                        evals_discrete=counter[0],
                    )
                    # rank of the closest discrete-stage candidate by the discrete score, before the
                    # quick MC moves it; tie rule: rank = number of candidates with a strictly
                    # smaller score
                    sc = np.asarray(d["score"], dtype=np.float64)
                    if len(R):
                        kc = int(np.argmin(err))
                        cur.update(
                            rank_closest_disc=int((sc < sc[kc]).sum()),
                            n_keep_disc=max(1, int(len(R) * rec.keep_fraction)),
                            disc_idx_closest=kc,
                            n_within_disc_1=int((err < 1.0).sum()),
                            n_within_disc_2=int((err < 2.0).sum()),
                        )
                elif ev == "quick_mc":
                    R = d["R"].astype(np.float64)  # sorted by post-quick-MC cost
                    err = reduced_misorientation_deg(R, np.broadcast_to(R_true, R.shape).copy())
                    n_keep = int(d["n_keep"])
                    k = int(np.argmin(err))
                    within = np.nonzero(err < BASIN_DEG)[0]
                    cur.update(
                        d_min_qmc=float(err[k]),
                        rank_closest=k,
                        n_keep=n_keep,
                        n_within_qmc=int(len(within)),
                        basin_survives=bool(len(within) and within.min() < n_keep),
                        best_in_kept_d=float(err[:n_keep].min()),
                        survives_1=bool((err[:n_keep] < 1.0).any()),
                        survives_2=bool((err[:n_keep] < 2.0).any()),
                        # post-MC rank of the candidate that was closest BEFORE the MC
                        rank_post_of_disc_closest=(
                            int(np.nonzero(np.asarray(d["perm"]) == cur["disc_idx_closest"])[0][0])
                            if "disc_idx_closest" in cur
                            else None
                        ),
                        evals_end=counter[0],
                    )
                    levels.append(dict(cur))
                    if d["level"] >= 2:
                        raise StopCoarse

            rec.recorder = hook
            rec.discrete_recorder = drec
            t0 = time.perf_counter()
            try:
                import io
                from contextlib import redirect_stdout

                with redirect_stdout(io.StringIO()):
                    rec.reconstruct_voxel(verts, v.phase, rng=np.random.default_rng(1000 + voxel))
            except StopCoarse:
                pass
            finally:
                rec.recorder = None
                rec.discrete_recorder = None
            out["arms"][arm] = dict(levels=levels, wall_s=time.perf_counter() - t0)
            print(f"  v{voxel} arm {arm}: {out['arms'][arm]['wall_s']:.0f}s", flush=True)
    finally:
        CF.VoxelCostFunction.evaluate = orig_eval  # type: ignore[method-assign]
        sp.cpp_overlap_rule, sp.ray_rule = False, None
        sp.cpp_center_screen = True
    return out


def _run_task(args: Tuple[str, int, str, int, Optional[str]]) -> str:
    variant, step, outdir, vox, arm = args
    path = Path(outdir) / (
        f"{variant}_v{vox}.json" if arm is None else f"{variant}_v{vox}_{arm}.json"
    )
    if path.exists():
        return f"{path.name} exists"
    init_worker(variant)
    t0 = time.perf_counter()
    res = peak_task(vox) if step == 1 else coarse_task(vox, [arm] if arm else None)
    res["variant"] = variant
    path.write_text(json.dumps(res, default=float) + "\n")
    msg = f"{path.name} {time.perf_counter() - t0:.0f}s"
    print(msg, flush=True)
    return msg


def run(step: int, workers: int, voxels: List[int], arms: List[str]) -> None:
    sub = (
        os.environ.get("RAY_PEAK_DIR", "peak")
        if step == 1
        else os.environ.get("RAY_COARSE_DIR", "coarse")
    )
    outdir = OUT / sub
    outdir.mkdir(parents=True, exist_ok=True)
    jobs = []
    if step == 1:
        jobs = [(var, 1, str(outdir), v, None) for var in VARIANTS for v in voxels]
    else:
        # longest arms first so the pool stays busy to the end
        jobs = [(var, 2, str(outdir), v, a) for a in arms[::-1] for var in VARIANTS for v in voxels]
    with mp.get_context("spawn").Pool(workers) as pool:
        for _ in pool.imap_unordered(_run_task, jobs, chunksize=1):
            pass


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd", choices=["pilot", "run", "summary"])
    ap.add_argument("--step", type=int, default=1)
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--arms", default="")
    a = ap.parse_args()
    voxels = select_voxels()
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "voxels.json").write_text(json.dumps(voxels) + "\n")
    if a.cmd == "pilot":
        init_worker("clean")
        if a.arms:
            ARM_NAMES[:] = a.arms.split(",")
        for vox in voxels[: a.limit or 2]:
            t0 = time.perf_counter()
            r1 = peak_task(vox)
            t1 = time.perf_counter()
            r2 = coarse_task(vox)
            t2 = time.perf_counter()
            print(f"v{vox}: step1 {t1 - t0:.0f}s, step2 {t2 - t1:.0f}s (contended)")
            (OUT / f"pilot_v{vox}.json").write_text(
                json.dumps(dict(peak=r1, coarse=r2), default=float) + "\n"
            )
    elif a.cmd == "run":
        run(
            a.step,
            a.workers,
            voxels[: a.limit] if a.limit else voxels,
            a.arms.split(",") if a.arms else ARM_NAMES,
        )
    else:
        import ray_overlap_summary  # noqa: F401


if __name__ == "__main__":
    main()
