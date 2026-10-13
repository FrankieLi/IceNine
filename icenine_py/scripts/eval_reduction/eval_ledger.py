"""Evaluation ledger of a no-start search (Task A step 1): cost evaluations, candidate counts and
wall time per level and per stage, for one voxel, on the Phase D clean full-sample images (Q-max 8).

Stages (per level L = 0 .. MaxLocalResolution, then once):
  screen      pixel-radius-3 cost, FZ points x local grid (run_discrete_search_spaced, phase 1)
  recost      complete cost of every locally passing candidate (phase 2); C++ does one more of these
              (see README of the task) which Python does not
  quick MC    10 steps, 5 restarts per surviving candidate
  ranking     candidates in / kept / dropped (the best quarter is kept)
  FindOptimal MC on the handed-over candidates, one row per candidate
  VarianceMinimizing, final evaluation

No library change: the stage boundaries are the events of the existing hooks
(``AdaptiveVoxelReconstructor.discrete_recorder`` and ``.recorder``); evaluations and their time are
counted by wrapping ``VoxelCostFunction.evaluate`` (as scripts/phase_d/seed_diag.py does). The wrap
adds about 1 microsecond per call; the unwrapped wall time is not measured separately.

  uv run python scripts/eval_reduction/eval_ledger.py run --voxel 15901 --tag cpp \
      --opt cpp_overlap_rule=1
  uv run python scripts/eval_reduction/eval_ledger.py summary \
      --dir benchmarks/eval_reduction/A/ledger

Options (--opt NAME=VALUE, a search_params field; ``ray_rule=R2`` as in seed_diag):
  cpp_overlap_rule=1   C++ CountQualifiedPeaks rule (config key CppOverlapRule)
  cpp_edge_rules=1     C++ edge conventions (config key CppEdgeRules)
  ray_rule=R2          ray-consistency rule, isotropic 2 px (implies cpp_overlap_rule)
"""

import argparse
import contextlib
import io
import json
import os
import sys
import time
from pathlib import Path
from typing import Any, Dict, List, Optional

import numpy as np

HERE = Path(__file__).resolve().parent
ICE = HERE.parents[1]
EX = ICE.parent / "Examples" / "Example2.ManyGrains"
OUT = ICE / "benchmarks" / "eval_reduction" / "A" / "ledger"
CONFIG = ICE / "scripts" / "phase_d" / "configs" / "ReconstructPhaseD_mc_clean.config"
STACK = EX / "ScatteringData_PhaseD" / "stacks" / "clean.npy"
SIX = [15901, 17034, 18558, 19889, 22867, 5242]
sys.path.insert(0, str(ICE / "scripts" / "common"))
sys.path.insert(0, str(ICE / "scripts" / "findoptimal_robustness"))


def apply_opts(setup: Any, opts: List[str]) -> None:
    for kv in opts:
        k, v = kv.split("=")
        if k == "ray_rule":
            from icenine.ray_overlap import RayRule

            reg = float(setup.exp_setup.config_file.detector_dist_deviation)
            setup.search_params.cpp_overlap_rule = True
            setup.search_params.ray_rule = (
                RayRule(arm="A", registration_mm=reg)
                if v == "A"
                else RayRule(arm="R", radius_px=int(v[1:]), registration_mm=reg)
            )
            continue
        cur = getattr(setup.search_params, k)
        if isinstance(cur, bool):
            setattr(setup.search_params, k, bool(int(v)))
        else:
            setattr(setup.search_params, k, type(cur)(float(v)))


def compact_calls(calls: List[Dict[str, Any]]) -> Dict[str, Any]:
    """The quick-MC calls (about 11,500 per search) as a summary, the FindOptimal and
    VarianceMinimizing calls in full."""
    quick = [c for c in calls if c["kind"] == "quick"]
    out: Dict[str, Any] = dict(
        find=[c for c in calls if c["kind"] == "find"],
        variance=[c for c in calls if c["kind"] == "variance"],
    )
    if quick:
        ev = np.array([c["evals"] for c in quick])
        out["quick"] = dict(
            n=len(quick),
            evals_total=int(ev.sum()),
            evals_per_call=dict(min=int(ev.min()), median=float(np.median(ev)), max=int(ev.max())),
            blocks_mean=float(np.mean([c["blocks"] for c in quick])),
            stop_codes={
                str(k): int(sum(1 for c in quick if c["stop"] == k))
                for k in sorted({c["stop"] for c in quick})
            },
        )
    return out


def run_one(voxel: int, tag: str, opts: List[str], seed: int, out_dir: Path) -> Dict[str, Any]:
    os.chdir(EX)
    import torch

    torch.set_num_threads(1)
    from icenine import cost_functions as CF
    from icenine import orientation_search as OS
    from icenine import reconstructor as RC
    from icenine.config_file import ConfigFile
    from icenine.experimental_data import ExperimentalData
    from icenine.mic_file import MicFile
    import preflight

    cfg = ConfigFile.from_file(str(CONFIG))
    exp = ExperimentalData.from_binary_memmap(STACK, 180, 2)
    setup = RC.setup_reconstruction(cfg, exp_data=exp)
    apply_opts(setup, opts)
    rec = RC.AdaptiveVoxelReconstructor(setup)
    v = setup.sample.get_mic().voxels[voxel]
    verts = RC._get_voxel_vertices(v)
    R_true = np.asarray(
        MicFile.read("SimInput/rand_500grains_1mm_neworient_s0.mic").voxels[voxel].orientation,
        dtype=np.float64,
    )

    # ---- counters: evaluations and their time by pixel radius ----
    cnt = {"n3": 0, "n0": 0, "t3": 0.0, "t0": 0.0}
    orig_eval = CF.VoxelCostFunction.evaluate

    def eval_w(self: Any, *args: Any, **kw: Any) -> Any:
        t0 = time.perf_counter()
        r = orig_eval(self, *args, **kw)
        dt = time.perf_counter() - t0
        if self.pixel_radius > 0:
            cnt["n3"] += 1
            cnt["t3"] += dt
        else:
            cnt["n0"] += 1
            cnt["t0"] += dt
        return r

    # the batched / interpolated screen (BatchedScreen, CppOmegaInterp) does not go through
    # evaluate(): count its orientations as pixel-radius-3 evaluations and time the calls
    from icenine import batched_screen as BS

    orig_exact, orig_interp = BS.ScreenEvaluator.overlaps_exact, BS.ScreenEvaluator.overlaps_interp

    def exact_w(self: Any, orientations: Any, *args: Any, **kw: Any) -> Any:
        t0 = time.perf_counter()
        r = orig_exact(self, orientations, *args, **kw)
        cnt["t3"] += time.perf_counter() - t0
        cnt["n3"] += len(orientations)
        return r

    def interp_w(self: Any, centres: Any, local_grid: Any, *args: Any, **kw: Any) -> Any:
        t0 = time.perf_counter()
        r = orig_interp(self, centres, local_grid, *args, **kw)
        cnt["t3"] += time.perf_counter() - t0
        cnt["n3"] += len(centres) * len(local_grid)
        return r

    def snap() -> Dict[str, float]:
        return dict(cnt, wall=time.perf_counter())

    def diff(a: Dict[str, float], b: Dict[str, float]) -> Dict[str, float]:
        return {k: b[k] - a[k] for k in a}

    # ---- MC call records (stop reasons, evaluations per call) ----
    calls: List[Dict[str, Any]] = []
    orig_mc = OS.MCOptimizer.optimize
    orig_var = OS.MCOptimizer.variance_minimizing_optimize

    def mc_w(self: Any, **kw: Any) -> Any:
        s0 = snap()
        r = orig_mc(self, **kw)
        d = diff(s0, snap())
        lr = dict(self.last_run)
        calls.append(
            dict(
                kind="quick" if kw["max_mc_steps"] == 10 else "find",
                evals=int(d["n0"] + d["n3"]),
                cost=float(r.cost),
                blocks=int(lr.get("n_blocks", 0)),
                steps=int(lr.get("steps_run", 0)),
                restarts=int(lr.get("n_restarts", 0)),
                accepts=int(lr.get("n_accept", 0)),
                stop=int(lr.get("stop", -1)),
            )
        )
        return r

    def var_w(self: Any, **kw: Any) -> Any:
        s0 = snap()
        r = orig_var(self, **kw)
        d = diff(s0, snap())
        calls.append(dict(kind="variance", evals=int(d["n0"] + d["n3"]), cost=float(r.cost)))
        return r

    # ---- stage boundaries from the existing hooks ----
    levels: List[Dict[str, Any]] = []
    find: List[Dict[str, Any]] = []
    st: Dict[str, Any] = dict(mark=None, cur=None, var=None, final=None, passed=0)

    def disc_hook(ev: str, d: Dict[str, Any]) -> None:
        if ev == "level":
            st["cur"] = dict(
                level=len(levels),
                fz_in=int(d["n_fz"]),
                n_local=int(d["n_local"]),
                n_recip=d.get("n_recip"),
                diameter_deg=float(np.degrees(d["diameter_rad"])),
                s0=snap(),
            )
            st["passed"] = 0
        elif ev == "fz":
            st["passed"] += len(d["j"])
        elif ev == "level_end":
            c = st["cur"]
            dd = diff(c.pop("s0"), snap())
            c.update(
                candidates=int(d["n_candidates"]),
                cliques=int(d["n_cliques"]),
                passed_screen=int(st["passed"]),
                screen_evals=int(dd["n3"]),
                recost_evals=int(dd["n0"]),
                screen_eval_s=dd["t3"],
                recost_eval_s=dd["t0"],
                discrete_s=dd["wall"],
                kept=0,
                dropped=int(d["n_candidates"]),
                quick_evals=0,
                quick_s=0.0,
            )
            levels.append(c)
            st["mark"] = snap()

    def hook(name: str, d: Dict[str, Any]) -> None:
        now = snap()
        if name == "quick_mc":
            c = levels[-1]
            dd = diff(st["mark"], now)
            c["quick_evals"] = int(dd["n0"] + dd["n3"])
            c["quick_s"] = dd["wall"]
            c["kept"] = int(d["n_keep"])
            c["dropped"] = c["candidates"] - c["kept"]
            c["quick_calls"] = int(len(d["R"]))
            st["mark"] = now
        elif name == "find_candidate":
            dd = diff(st["mark"], now)
            find.append(
                dict(
                    index=int(d["index"]),
                    evals=int(dd["n0"] + dd["n3"]),
                    wall_s=dd["wall"],
                    cost=float(d["cost"]),
                )
            )
            st["mark"] = now
        elif name == "variance":
            dd = diff(st["mark"], now)
            st["var"] = dict(evals=int(dd["n0"] + dd["n3"]), wall_s=dd["wall"])
            st["mark"] = now
        elif name == "final":
            dd = diff(st["mark"], now)
            st["final"] = dict(evals=int(dd["n0"] + dd["n3"]), wall_s=dd["wall"])
            st["mark"] = now

    rec.discrete_recorder = disc_hook
    rec.recorder = hook
    CF.VoxelCostFunction.evaluate = eval_w  # type: ignore[method-assign]
    BS.ScreenEvaluator.overlaps_exact = exact_w  # type: ignore[method-assign]
    BS.ScreenEvaluator.overlaps_interp = interp_w  # type: ignore[method-assign]
    OS.MCOptimizer.optimize = mc_w  # type: ignore[method-assign]
    OS.MCOptimizer.variance_minimizing_optimize = var_w  # type: ignore[method-assign]
    pf = preflight.preflight()
    try:
        preflight.require_quiet(info=pf)
        quiet = True
    except preflight.MachineBusyError as e:
        quiet = False
        pf = dict(pf, not_quiet=str(e)[:300])
    rng = np.random.default_rng(seed)
    t0 = time.perf_counter()
    with contextlib.redirect_stdout(io.StringIO()):
        res = rec.reconstruct_voxel(voxel_vertices=verts, phase_index=v.phase, rng=rng)
    wall = time.perf_counter() - t0
    CF.VoxelCostFunction.evaluate = orig_eval  # type: ignore[method-assign]
    BS.ScreenEvaluator.overlaps_exact = orig_exact  # type: ignore[method-assign]
    BS.ScreenEvaluator.overlaps_interp = orig_interp  # type: ignore[method-assign]
    OS.MCOptimizer.optimize = orig_mc  # type: ignore[method-assign]
    OS.MCOptimizer.variance_minimizing_optimize = orig_var  # type: ignore[method-assign]

    from csl import reduced_misorientation_deg

    err = float(reduced_misorientation_deg(np.asarray(res.orientation), R_true))
    fo = dict(rec.last_find_optimal)
    for c in levels:
        c["stage_evals"] = c["screen_evals"] + c["recost_evals"] + c["quick_evals"]
    total = cnt["n3"] + cnt["n0"]
    totals = dict(
        screen=sum(c["screen_evals"] for c in levels),
        recost=sum(c["recost_evals"] for c in levels),
        quick=sum(c["quick_evals"] for c in levels),
        find=sum(f["evals"] for f in find),
        variance=(st["var"] or {}).get("evals", 0),
        final=(st["final"] or {}).get("evals", 0),
    )
    totals["sum"] = sum(totals.values())
    out = dict(
        voxel=voxel,
        tag=tag,
        opts=opts,
        seed=seed,
        preflight=pf,
        label="single-worker" + ("" if quiet else ", contended"),
        quiet_preflight=quiet,
        wall_s=wall,
        evals_total=int(total),
        evals_by_radius=dict(r3=int(cnt["n3"]), r0=int(cnt["n0"])),
        eval_time_s=cnt["t3"] + cnt["t0"],
        mean_us_per_eval=1e6 * (cnt["t3"] + cnt["t0"]) / max(total, 1),
        library_eval_counts=list(rec.last_eval_counts),
        totals=totals,
        levels=levels,
        find=find,
        find_optimal=fo,
        variance=st["var"],
        final=st["final"],
        calls=compact_calls(calls),
        err_deg=err,
        cost_final=float(res.cost),
        c_ext=bool(getattr(CF, "_HAS_C_RASTERIZE", False)),
    )
    out_dir.mkdir(parents=True, exist_ok=True)
    fn = out_dir / f"{tag}_v{voxel}.json"
    fn.write_text(json.dumps(out, indent=1, default=float) + "\n")
    print(
        f"{tag} v{voxel}: {wall:.1f}s evals={total} (sum of stages {totals['sum']}) "
        f"err={err:.4f} {out['label']} -> {fn.name}",
        flush=True,
    )
    return out


def summary(d: Path) -> None:
    """ledger.json and tables.md from the run files in ``d / 'runs'``."""
    import doc_tables as DT

    runs = [json.loads(p.read_text()) for p in sorted((d / "runs").glob("*_v*.json"))]
    tags = sorted({r["tag"] for r in runs})
    ledger: Dict[str, Any] = {}
    for t in tags:
        ledger[t] = {str(r["voxel"]): r for r in runs if r["tag"] == t}
    (d / "ledger.json").write_text(json.dumps(ledger, indent=1, default=float) + "\n")

    T: Dict[str, str] = {}
    tags = [t for t in tags if not t.startswith("timing_")]  # quiet repeats: step3_summary.py
    runs = [r for r in runs if r["tag"] in tags]
    # per-stage totals by arm and voxel
    rows = []
    for t in tags:
        for r in sorted((r for r in runs if r["tag"] == t), key=lambda r: SIX.index(r["voxel"])):
            tt = r["totals"]
            rows.append(
                dict(
                    arm=t,
                    voxel=r["voxel"],
                    screen=tt["screen"],
                    recost=tt["recost"],
                    quick=tt["quick"],
                    find=tt["find"],
                    variance=tt["variance"],
                    final=tt["final"],
                    total=r["evals_total"],
                    wall_s=r["wall_s"],
                    label=r["label"],
                    err=r["err_deg"],
                )
            )
    cols = [
        "arm",
        "voxel",
        "screen",
        "recost",
        "quick",
        "find",
        "variance",
        "final",
        "total",
        "wall_s",
        "label",
        "err",
    ]
    T["ledger_stages"] = DT.markdown_table(
        rows,
        cols,
        formats={
            "screen": ",d",
            "recost": ",d",
            "quick": ",d",
            "find": ",d",
            "variance": ",d",
            "final": ",d",
            "total": ",d",
            "wall_s": ".0f",
            "err": ".3f",
        },
    )
    # means over the six voxels
    mrows = []
    for t in tags:
        rr = [r for r in runs if r["tag"] == t]
        if len(rr) < 1:
            continue
        mean = lambda f: float(np.mean([f(r) for r in rr]))  # noqa: E731
        mrows.append(
            dict(
                arm=t,
                n=len(rr),
                screen=mean(lambda r: r["totals"]["screen"]),
                recost=mean(lambda r: r["totals"]["recost"]),
                quick=mean(lambda r: r["totals"]["quick"]),
                find=mean(lambda r: r["totals"]["find"]),
                variance=mean(lambda r: r["totals"]["variance"]),
                final=mean(lambda r: r["totals"]["final"]),
                total=mean(lambda r: r["evals_total"]),
                wall_s=mean(lambda r: r["wall_s"]),
                us_per_eval=mean(lambda r: r["mean_us_per_eval"]),
            )
        )
    T["ledger_mean"] = DT.markdown_table(
        mrows,
        [
            "arm",
            "n",
            "screen",
            "recost",
            "quick",
            "find",
            "variance",
            "final",
            "total",
            "wall_s",
            "us_per_eval",
        ],
        formats={
            k: ",.0f"
            for k in (
                "screen",
                "recost",
                "quick",
                "find",
                "variance",
                "final",
                "total",
                "wall_s",
                "us_per_eval",
            )
        }
        | {"n": "d"},
    )
    # per-level table for the first voxel, each arm
    for t in tags:
        r15 = next((r for r in runs if r["tag"] == t and r["voxel"] == SIX[0]), None)
        if r15 is None:
            continue
        lrows = []
        for c in r15["levels"]:
            lrows.append(
                dict(
                    level=c["level"],
                    fz_in=c["fz_in"],
                    local=c["n_local"],
                    screen=c["screen_evals"],
                    passed=c["passed_screen"],
                    recost=c["recost_evals"],
                    cands=c["candidates"],
                    quick=c["quick_evals"],
                    kept=c["kept"],
                    dropped=c["dropped"],
                    screen_s=c["screen_eval_s"],
                    recost_s=c["recost_eval_s"],
                    discrete_s=c["discrete_s"],
                    quick_s=c["quick_s"],
                )
            )
        T[f"levels_{t}_v{SIX[0]}"] = DT.markdown_table(
            lrows,
            [
                "level",
                "fz_in",
                "local",
                "screen",
                "passed",
                "recost",
                "cands",
                "quick",
                "kept",
                "dropped",
                "screen_s",
                "recost_s",
                "discrete_s",
                "quick_s",
            ],
            formats={"screen_s": ".1f", "recost_s": ".1f", "discrete_s": ".1f", "quick_s": ".1f"},
        )
    DT.write_tables(d / "tables.md", T)
    print(f"wrote {d / 'ledger.json'} and tables.md ({len(runs)} runs)")


def shrink(d: Path) -> None:
    """Rewrite stored run files whose `calls` is still the full list (older runs of this task)."""
    for p in sorted((d / "runs").glob("*_v*.json")):
        r = json.loads(p.read_text())
        if isinstance(r.get("calls"), list):
            r["calls"] = compact_calls(r["calls"])
            p.write_text(json.dumps(r, indent=1, default=float) + "\n")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd", choices=["run", "summary", "shrink"])
    ap.add_argument("--voxel", type=int)
    ap.add_argument("--tag", default="base")
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--opt", action="append", default=[])
    ap.add_argument("--dir", default=str(OUT))
    a = ap.parse_args()
    d = Path(a.dir).resolve()  # before run_one changes the working directory
    if a.cmd == "run":
        run_one(a.voxel, a.tag, a.opt, a.seed, d / "runs")
    elif a.cmd == "shrink":
        shrink(d)
    else:
        summary(d)


if __name__ == "__main__":
    main()
