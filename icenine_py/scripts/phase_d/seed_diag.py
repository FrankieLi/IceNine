"""Phase D seed-cost diagnosis: stage-instrumented no-start `reconstruct_voxel` on the full-sample
images (memmap loader) or on an isolated one-voxel render of the same voxel, with optional
option overrides. One process, one voxel per call; results go to
benchmarks/phase_d_seed_diag/runs/<tag>_v<voxel>.json.

Usage (from icenine_py/):
  uv run python scripts/phase_d/seed_diag.py select
  uv run python scripts/phase_d/seed_diag.py run --voxel 6628 --images full --tag base
  uv run python scripts/phase_d/seed_diag.py run --voxel 6628 --images isolated --tag base
Options: --opt NAME=VALUE (search_params field, or rec.<attr>), see OPTS below.
"""

import argparse
import contextlib
import io
import json
import math
import os
import sys
import time
from pathlib import Path
from typing import Any, Dict, List

import numpy as np

HERE = Path(__file__).resolve().parent
ICE = HERE.parents[1]
EX = HERE.parents[2] / "Examples" / "Example2.ManyGrains"
OUT = ICE / "benchmarks" / "phase_d_seed_diag"
sys.path.insert(0, str(ICE / "scripts" / "common"))
sys.path.insert(0, str(ICE / "scripts" / "findoptimal_robustness"))
sys.path.insert(0, str(HERE))
MIC = "SimInput/rand_500grains_1mm_neworient_s0.mic"
SCRATCH = Path(os.environ.get("SEED_DIAG_SCRATCH", "/tmp/seed_diag_stacks"))


def select(n_near: int = 3, n_far: int = 3, seed: int = 7) -> List[Dict[str, Any]]:
    """Grain-interior voxels: every voxel within 2.1 side lengths is in the same grain; the
    grain has >= 40 voxels. Half near the rotation axis (r < 0.12), half far (r > 0.38)."""
    from scipy.spatial import cKDTree

    raw = np.loadtxt(EX / MIC, skiprows=1)
    pos = raw[:, :2]
    grain = np.load(EX / "SimInput" / "rand_500grains_1mm_neworient_s0_grainmap.npy")
    side = 0.009375
    tree = cKDTree(pos)
    sizes = np.bincount(grain)
    r = np.linalg.norm(pos, axis=1)
    ok = np.zeros(len(pos), bool)
    for i, nb in enumerate(tree.query_ball_point(pos, 2.1 * side)):
        ok[i] = bool((grain[nb] == grain[i]).all()) and sizes[grain[i]] >= 40
    rng = np.random.default_rng(seed)
    out = []
    for label, mask, n in (
        ("near", ok & (r < 0.12), n_near),
        ("far", ok & (r > 0.38), n_far),
    ):
        cand = np.nonzero(mask)[0]
        for i in rng.choice(cand, n, replace=False):
            out.append(dict(voxel=int(i), grain=int(grain[i]), r=float(r[i]), where=label))
    return out


def build_isolated_stack(voxel: int) -> Path:
    """Render only `voxel` (clean, same forward path as the full render) into a uint8 stack."""
    path = SCRATCH / f"iso_{voxel}.npy"
    if path.exists():
        return path
    SCRATCH.mkdir(parents=True, exist_ok=True)
    import tempfile
    import render_full as RF

    with tempfile.TemporaryDirectory() as tmp:
        job = dict(
            worker=0,
            workers=1,
            mic=MIC,
            out=tmp,
            voxel_idx=[voxel],
            batch=2000,
            noise_seed=RF.NOISE_SEED,
            max_q=0.0,
            variants=["clean"],
            keep=True,
        )
        with contextlib.redirect_stdout(io.StringIO()):
            kept = RF._worker(job)["kept"]
    st = np.lib.format.open_memmap(str(path), mode="w+", dtype=np.uint8, shape=(360, 2048, 2048))
    for (di, i), (kk, jj, _) in kept.items():
        st[di * 180 + i, kk, jj] = 1
    st.flush()
    del st
    return path


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd", choices=["select", "run"])
    ap.add_argument("--voxel", type=int)
    ap.add_argument("--images", default="full", choices=["full", "isolated"])
    ap.add_argument("--variant", default="clean")
    ap.add_argument("--tag", default="base")
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--opt", action="append", default=[])
    ap.add_argument("--varcap", type=int, default=0, help="hard cap on VarianceMinimizing steps")
    ap.add_argument(
        "--topk0", type=int, default=0, help="keep only the top K discrete candidates at level 0"
    )
    ap.add_argument("--topk", type=int, default=0, help="same, every level >= 1")
    ap.add_argument("--save-cands", action="store_true")
    a = ap.parse_args()
    if a.cmd == "select":
        sel = select()
        (OUT / "voxels.json").write_text(json.dumps(sel, indent=1) + "\n")
        print(json.dumps(sel, indent=1))
        return
    run(a)


# --------------------------------------------------------------------------------------------


def run(a: argparse.Namespace) -> None:
    os.chdir(EX)
    import torch

    torch.set_num_threads(1)
    from csl import reduced_misorientation_deg
    from icenine import cost_functions as CF
    from icenine import orientation_search as OS
    from icenine import reconstructor as RC
    from icenine.config_file import ConfigFile
    from icenine.experimental_data import ExperimentalData
    from icenine.mic_file import MicFile
    import preflight

    cfg = ConfigFile.from_file(str(HERE / "configs" / f"ReconstructPhaseD_mc_{a.variant}.config"))
    stack = (
        EX / "ScatteringData_PhaseD" / "stacks" / f"{a.variant}.npy"
        if a.images == "full"
        else build_isolated_stack(a.voxel)
    )
    exp = ExperimentalData.from_binary_memmap(stack, 180, 2)
    setup = RC.setup_reconstruction(cfg, exp_data=exp)
    rec = RC.AdaptiveVoxelReconstructor(setup)
    for kv in a.opt:
        k, v = kv.split("=")
        tgt, name = (rec, k[4:]) if k.startswith("rec.") else (setup.search_params, k)
        cur = getattr(tgt, name)
        if isinstance(cur, str):
            setattr(tgt, name, v)
        else:
            setattr(tgt, name, bool(int(v)) if isinstance(cur, bool) else type(cur)(float(v)))
    mic = setup.sample.get_mic()
    truth = MicFile.read(MIC)
    v = mic.voxels[a.voxel]
    R_true = np.asarray(truth.voxels[a.voxel].orientation, dtype=np.float64)
    verts = RC._get_voxel_vertices(v)

    # ---- instrumentation ----
    S: Dict[str, Any] = dict(phase="init", level=-1)
    ev: Dict[str, Dict[str, float]] = {}  # key phase|pr -> n, pass, t
    calls: List[Dict[str, Any]] = []
    levels: List[Dict[str, Any]] = []
    orig_eval = CF.VoxelCostFunction.evaluate
    orig_disc = RC.run_discrete_search_spaced
    orig_mc = OS.MCOptimizer.optimize
    orig_var = OS.MCOptimizer.variance_minimizing_optimize

    def n_total() -> int:
        return int(sum(d["n"] for d in ev.values()))

    def eval_w(self, *args, **kw):  # type: ignore[no-untyped-def]
        t0 = time.perf_counter()
        r = orig_eval(self, *args, **kw)
        dt = time.perf_counter() - t0
        d = ev.setdefault(f"{S['phase']}|pr{self.pixel_radius}", dict(n=0, passed=0, t=0.0))
        d["n"] += 1
        d["t"] += dt
        d["passed"] += int(r.peak_overlap > 0)
        return r

    def disc_w(**kw):  # type: ignore[no-untyped-def]
        S["level"] += 1
        S["phase"] = f"discrete_L{S['level']}"
        n0, t0 = n_total(), time.perf_counter()
        out = orig_disc(**kw)
        n_full = len(out)
        cap = a.topk0 if S["level"] == 0 else a.topk
        if cap and len(out) > cap:
            out = out[:cap]  # sorted best first by the discrete-stage cost
        levels.append(
            dict(
                level=S["level"],
                n_fz=int(len(kw["fz_orientations"])),
                n_local=int(len(kw["local_grid"])),
                n_returned=n_full,
                n_after_cap=len(out),
                discrete_evals=n_total() - n0,
                discrete_s=time.perf_counter() - t0,
                n_recip=int(
                    next(iter(kw["global_cost_fn"]._phase_recip_vecs.values()))[0].shape[0]
                ),
                diameter_deg=float(np.degrees(kw["angular_radius"])),
            )
        )
        S["phase"] = f"quick_L{S['level']}"
        return out

    def mc_w(self, **kw):  # type: ignore[no-untyped-def]
        quick = kw["max_mc_steps"] == 10
        if not quick:
            S["phase"] = "find"
        n0, t0 = n_total(), time.perf_counter()
        r = orig_mc(self, **kw)
        n = n_total() - n0
        hit = float(r.overlap_info.hit_ratio) if r.overlap_info is not None else float("nan")
        steps = kw["max_mc_steps"]
        why = (
            "step_cap"
            if n >= steps + 1
            else (
                "max_convergence_cost"
                if r.cost < kw["max_convergence_cost"]
                else "restarts_exhausted"
            )
        )
        calls.append(
            dict(
                phase=S["phase"],
                quick=quick,
                evals=n,
                s=time.perf_counter() - t0,
                cost=float(r.cost),
                hit_ratio=hit,
                ended_by=why,
            )
        )
        return r

    def capped_var(
        self,
        initial_orientation,
        search_box_side,
        max_mc_steps,
        successive_restarts,
        max_convergence_cost,
        convergence_variance,
        hard_cap,
        cost_fn_angular_resolution=math.radians(0.5),
    ):  # type: ignore[no-untyped-def]
        """variance_minimizing_optimize with one added exit: total steps >= hard_cap."""
        from icenine.orientation_search import SearchCandidate as SC
        from icenine.orientation_search import (
            _quat_multiply,
            matrix_to_quaternion,
            quaternion_to_matrix,
        )

        gq = matrix_to_quaternion(initial_orientation)
        cq = gq.copy()
        ginfo = self.cost_fn.evaluate(initial_orientation, self.voxel_vertices, self.phase_index)
        gmin = ginfo.cost
        sub = math.tan(search_box_side) / math.sqrt(48.0)
        taken, mx = 0, max_mc_steps
        while taken < mx and taken < hard_cap:
            n_sub = max(int(math.ceil(sub / cost_fn_angular_resolution) ** 2.7), 10)
            new_o, new_c, var, new_i = self._zero_temp_with_variance(
                quaternion_to_matrix(cq), sub, n_sub
            )
            if var > convergence_variance:
                mx += n_sub
            taken += n_sub
            if new_c >= gmin:
                sub = min(2.0 * sub, search_box_side)
                hb = search_box_side / 2.0
                rq = self._grid_gen.get_near_identity_point(
                    self._rng.uniform(-hb, hb),
                    self._rng.uniform(-hb, hb),
                    self._rng.uniform(-hb, hb),
                )
                cq = _quat_multiply(rq, gq)
            else:
                gmin, gq, ginfo = new_c, matrix_to_quaternion(new_o), new_i
                cq = gq.copy()
                sub *= 0.5
            if gmin < max_convergence_cost and abs(var) < convergence_variance:
                break
        return SC(orientation=quaternion_to_matrix(gq), cost=gmin, overlap_info=ginfo)

    def var_w(self, **kw):  # type: ignore[no-untyped-def]
        S["phase"] = "variance"
        S["err_before_variance"] = float(
            reduced_misorientation_deg(np.asarray(kw["initial_orientation"]), R_true)
        )
        n0, t0 = n_total(), time.perf_counter()
        r = capped_var(self, hard_cap=a.varcap, **kw) if a.varcap else orig_var(self, **kw)
        calls.append(
            dict(
                phase="variance",
                quick=False,
                evals=n_total() - n0,
                s=time.perf_counter() - t0,
                cost=float(r.cost),
                hit_ratio=float("nan"),
                ended_by="variance",
            )
        )
        S["phase"] = "final"
        return r

    events: List[Any] = []
    cands: Dict[str, Any] = {}

    def recorder(name: str, d: Dict[str, Any]) -> None:
        events.append(
            (
                name,
                {
                    k: (np.asarray(x).shape if hasattr(x, "shape") else x)
                    for k, x in d.items()
                    if k in ("level", "n_keep")
                },
            )
        )
        if a.save_cands and d.get("level") == 0:
            for k in ("R", "score", "perm", "cost"):
                if k in d:
                    cands[f"L0_{name}_{k}"] = np.asarray(d[k])

    rec.recorder = recorder
    CF.VoxelCostFunction.evaluate = eval_w  # type: ignore[method-assign]
    RC.run_discrete_search_spaced = disc_w  # type: ignore[assignment]
    OS.MCOptimizer.optimize = mc_w  # type: ignore[method-assign]
    OS.MCOptimizer.variance_minimizing_optimize = var_w  # type: ignore[method-assign]

    # quality of the truth under this image set (uninstrumented-phase evaluate)
    S["phase"] = "truth"
    local = rec._make_local_cost_fn()
    info_t = local.evaluate(R_true.astype(np.float32), verts, v.phase)
    q_true, hit_true = float(1.0 - info_t.cost), float(info_t.hit_ratio)
    ev.clear()

    pf = preflight.preflight()
    try:
        preflight.require_quiet(info=pf)
        quiet = True
    except preflight.MachineBusyError as e:
        quiet = False
        pf = dict(pf, not_quiet=str(e)[:300])
    rng = np.random.default_rng(a.seed)
    S["phase"] = "start"
    t0 = time.perf_counter()
    with contextlib.redirect_stdout(io.StringIO()):
        res = rec.reconstruct_voxel(voxel_vertices=verts, phase_index=v.phase, rng=rng)
    wall = time.perf_counter() - t0
    CF.VoxelCostFunction.evaluate = orig_eval  # type: ignore[method-assign]
    err = float(reduced_misorientation_deg(np.asarray(res.orientation), R_true))
    fo = dict(rec.last_find_optimal)
    for e in events:
        pass
    keep = {}
    for name, d in events:
        if name == "quick_mc":
            keep[d["level"]] = d["n_keep"]
    for lv in levels:
        lv["n_kept"] = int(keep.get(lv["level"], -1))
    by = {}
    for k, d in ev.items():
        by[k] = dict(d, mean_us=1e6 * d["t"] / max(d["n"], 1))
    tot_n = n_total()
    tot_t = sum(d["t"] for d in ev.values())
    out = dict(
        voxel=a.voxel,
        images=a.images,
        variant=a.variant,
        tag=a.tag,
        seed=a.seed,
        opts=a.opt,
        quiet_preflight=quiet,
        preflight=pf,
        label="single-process" + ("" if quiet else ", contended"),
        wall_s=wall,
        evals_total=tot_n,
        eval_time_s=tot_t,
        mean_us_per_eval=1e6 * tot_t / max(tot_n, 1),
        err_deg=err,
        err_before_variance=S.get("err_before_variance"),
        varcap=a.varcap,
        topk0=a.topk0,
        topk=a.topk,
        cost_final=float(res.cost),
        hit_ratio_final=float(res.overlap_info.hit_ratio),
        q_true=q_true,
        hit_true=hit_true,
        find_optimal=fo,
        levels=levels,
        calls=calls,
        by_phase=by,
        n_quick_mc=sum(1 for c in calls if c["quick"]),
        n_find=sum(1 for c in calls if c["phase"] == "find"),
        find_evals=sum(c["evals"] for c in calls if c["phase"] == "find"),
        quick_evals=sum(c["evals"] for c in calls if c["quick"]),
        variance_evals=sum(c["evals"] for c in calls if c["phase"] == "variance"),
        c_ext=bool(getattr(CF, "_HAS_C_RASTERIZE", False)),
    )
    (OUT / "runs").mkdir(parents=True, exist_ok=True)
    if a.save_cands:
        np.savez_compressed(
            OUT / "runs" / f"{a.tag}_{a.images}_v{a.voxel}_cands.npz", R_true=R_true, **cands
        )
    fn = OUT / "runs" / f"{a.tag}_{a.images}_v{a.voxel}.json"
    fn.write_text(json.dumps(out, indent=1, default=float) + "\n")
    print(
        f"{a.tag} {a.images} v{a.voxel}: {wall:.1f}s evals={tot_n} err={err:.4f} "
        f"q_true={q_true:.3f} find {fo.get('n_evaluated')}/{fo.get('n_candidates')} "
        f"conv={fo.get('converged')} -> {fn.name}"
    )


if __name__ == "__main__":
    main()
