#!/usr/bin/env python3
"""E0 (200 voxels, seed 0, clean + realistic, per-voxel images) with the opt-in C++ overlap rule
(config key CppOverlapRule) on, `mc` arm only, paired per voxel against the stored `mc` result of
the C++-faithful MC rerun (scripts/cma_finisher/cache/e0_run_mcfaithful: 40 / 25 wrong of 200).

Also records, per case, the level-0 candidate set against the truth: number of candidates and the
smallest symmetry-reduced misorientation to the true orientation (truth-basin check). The `off`
arm with --level0-only gives the same level-0 record for the library default (stops after level 0).

  uv run python scripts/phase_d/e0_cpprule.py pilot --workers 2 --limit 4
  uv run python scripts/phase_d/e0_cpprule.py run --arm on --workers 10
  uv run python scripts/phase_d/e0_cpprule.py run --arm off --level0-only --workers 10
  uv run python scripts/phase_d/e0_cpprule.py summary

Q8 (no-drop realistic images, scripts/phase_d/q8_e0_summary.py): `--variant realistic_nodrop`
runs an arm on the E0 images without the random spot drop (build them first with
scripts/findoptimal_robustness/e0_run.py images); `--variant` also takes clean and all. A task
needs only its cached image (cache/images), not the published cache/e0, and `--meta-dir` moves the
run_meta file (default: the Q3 directory) so a Q8 run does not touch the tracked Q3 records.

Evaluation-reduction task A arms (scripts/eval_reduction/e0_interp_summary.py reads them):
  edge    C++ overlap rule + CppEdgeRules
  interp  C++ overlap rule + CppEdgeRules + CppOmegaInterp + BatchedScreen
  lenint  lenient (default) rule + CppEdgeRules + CppOmegaInterp + BatchedScreen
"""

import argparse
import json
import sys
import time
from pathlib import Path
from typing import Any, Dict, List, Tuple

import numpy as np

HERE = Path(__file__).resolve().parent
ICE = HERE.parents[1]
sys.path.insert(0, str(ICE / "scripts" / "findoptimal_robustness"))
sys.path.insert(0, str(ICE / "scripts" / "common"))
import common as C  # noqa: E402
import doc_tables  # noqa: E402
import stats as S  # noqa: E402

CACHE = HERE / "cache"
BASE = ICE / "scripts" / "cma_finisher" / "cache" / "e0_run_mcfaithful"  # stored `mc` arm
OUT = ICE / "benchmarks" / "phase_d_seed_diag" / "level0" / "e0"
BASIN_DEG = 3.0  # truth basin: a level-0 candidate within this symmetry-reduced angle of the truth


class StopL0(Exception):
    pass


def task(item: Tuple[Any, ...]) -> str:
    from icenine import reconstructor as RC
    from optimizer_sweep import voxel_context

    vidx, vpos, variant, path, arm, l0only, msin = item
    W = C.get_worker()
    t0 = time.time()
    keys = np.load(C.CACHE_DIR / "images" / f"v{vidx}_{variant}.npz")["keys"]
    C.attach(keys)
    vctx = voxel_context(W.ctx, vidx)
    rec = W.rec
    Rt = np.asarray(vctx.R_true, dtype=np.float64)
    rec.min_sin_eta = msin
    sp = rec.setup.search_params
    sp.cpp_overlap_rule = arm in ("on", "rule", "rayR2", "rayA", "edge", "interp")
    # evaluation-reduction task A arms (opt-in switches; all off for the other arms)
    sp.cpp_edge_rules = arm in ("edge", "interp", "lenint")
    sp.cpp_omega_interp = arm in ("interp", "lenint")
    sp.batched_screen = arm in ("interp", "lenint")
    if arm in ("rayR2", "rayA"):  # ray rule at levels 0-2 only (owner decision 2026-10-12)
        from icenine.ray_overlap import RayRule

        reg = float(rec.setup.exp_setup.config_file.detector_dist_deviation)
        rec.setup.search_params.ray_rule = (
            RayRule(arm="A", registration_mm=reg)
            if arm == "rayA"
            else RayRule(arm="R", radius_px=2, registration_mm=reg)
        )
    rec.setup.search_params.cpp_center_screen = arm != "rule"
    l0: Dict[str, Any] = {}
    orig = RC.run_discrete_search_spaced

    def wrapped(**kw: Any) -> Any:
        out = orig(**kw)
        if not l0:
            R = np.stack([c.orientation for c in out]).astype(np.float64)
            d = C.err_deg(R, np.broadcast_to(Rt, R.shape).copy())
            l0.update(
                n_l0=len(out), min_dist_l0=float(d.min()), n_within=int((d < BASIN_DEG).sum())
            )
            if l0only:
                raise StopL0
        return out

    RC.run_discrete_search_spaced = wrapped
    res = None
    try:
        t1 = time.perf_counter()
        with C.quiet():
            res = rec.reconstruct_voxel(vctx.vertices, vctx.voxel.phase, rng=C.run_seed(vpos, 0))
        dt = time.perf_counter() - t1
    except StopL0:
        dt = time.perf_counter() - t1
    finally:
        RC.run_discrete_search_spaced = orig
        rec.setup.search_params.cpp_overlap_rule = False
        rec.setup.search_params.cpp_edge_rules = False
        rec.setup.search_params.cpp_omega_interp = False
        rec.setup.search_params.batched_screen = False
        rec.setup.search_params.ray_rule = None
        rec.setup.search_params.cpp_center_screen = True
        rec.min_sin_eta = 0.3
    out: Dict[str, Any] = dict(R_true=Rt, runtime=dt, **l0)
    if res is not None:
        g, loc, _ = rec.last_eval_counts
        out.update(
            R=np.asarray(res.orientation, dtype=np.float64), cost=float(res.cost),
            evals_global=g, evals_local=loc,
        )  # fmt: skip
    np.savez_compressed(path, **out)
    return f"voxel {vidx} {variant} {arm} {time.time() - t0:.0f}s"


def work_items(
    root: Path, arm: str, l0only: bool, limit: int = 0, msin: float = 0.3, only_var: str = ""
) -> List[Tuple[Any, ...]]:
    info = dict(np.load(C.OUT_DIR / "voxels.npz"))
    vox = [int(v) for v in info["voxel_indices"]][: C.N_VOXELS]
    root.mkdir(parents=True, exist_ok=True)
    out = []
    for vpos, v in enumerate(vox):
        for var in [only_var] if only_var else C.VARIANTS:
            # a case that could not be built has no cached image (e0_run.py images / run)
            if not (C.CACHE_DIR / "images" / f"v{v}_{var}.npz").exists():
                continue
            out.append((v, vpos, var, str(root / f"v{v}_{var}.npz"), arm, l0only, msin))
    return out[:limit] if limit else out


def summary(arm: str = "on") -> None:
    items = work_items(CACHE / arm, arm, False)
    rows, summ = [], {}
    for var in C.VARIANTS:
        sel = [it for it in items if it[2] == var and Path(it[3]).exists()]
        if not sel:
            continue
        z = [np.load(it[3]) for it in sel]
        b = [np.load(BASE / Path(it[3]).name) for it in sel]
        Rt = np.stack([x["R_true"] for x in z])
        e_new = C.err_deg(np.stack([x["R"] for x in z]), Rt)
        e_old = C.err_deg(np.stack([x["mc_R"] for x in b]), Rt)
        wn, wo = e_new > C.WRONG_DEG, e_old > C.WRONG_DEG
        only_new, only_old = S.paired_discordant(wn, wo)
        n = len(sel)
        lo_n, hi_n = S.wilson(int(wn.sum()), n)
        lo_o, hi_o = S.wilson(int(wo.sum()), n)
        ev_new = np.array([int(x["evals_global"]) + int(x["evals_local"]) for x in z])
        ev_old = np.array([int(x["mc_evals_global"]) + int(x["mc_evals_local"]) for x in b])
        both = ~wn & ~wo
        label = "realistic" if var == "all" else var
        s: Dict[str, Any] = dict(
            n=n, wrong_before=int(wo.sum()), wrong_after=int(wn.sum()),
            only_after_wrong=only_new, only_before_wrong=only_old,
            both_wrong=int((wn & wo).sum()),
            median_err_right_both_before=float(np.median(e_old[both])),
            median_err_right_both_after=float(np.median(e_new[both])),
            mean_evals_before=float(ev_old.mean()), mean_evals_after=float(ev_new.mean()),
            mean_cost_before=float(np.mean([float(x["mc_cost"]) for x in b])),
            mean_cost_after=float(np.mean([float(x["cost"]) for x in z])),
            median_n_l0_after=float(np.median([int(x["n_l0"]) for x in z])),
            basin_miss_after=int(sum(float(x["min_dist_l0"]) > BASIN_DEG for x in z)),
        )  # fmt: skip
        off = [Path(CACHE / "off" / Path(it[3]).name) for it in sel]
        if all(p.exists() for p in off):
            zo = [np.load(p) for p in off]
            s["median_n_l0_before"] = float(np.median([int(x["n_l0"]) for x in zo]))
            s["basin_miss_before"] = int(sum(float(x["min_dist_l0"]) > BASIN_DEG for x in zo))
            mb = np.array([float(x["min_dist_l0"]) > BASIN_DEG for x in zo])
            ma = np.array([float(x["min_dist_l0"]) > BASIN_DEG for x in z])
            s["basin_only_after_missing"], s["basin_only_before_missing"] = S.paired_discordant(
                ma, mb
            )
        summ[var] = s
        rows.append(
            dict(variant=label, n=n,
                 wrong_before=f"{int(wo.sum())}/{n} [{lo_o:.3f}, {hi_o:.3f}]",
                 wrong_after=f"{int(wn.sum())}/{n} [{lo_n:.3f}, {hi_n:.3f}]",
                 broken=only_new, fixed=only_old, both_wrong=s["both_wrong"],
                 med_both_before=s["median_err_right_both_before"],
                 med_both_after=s["median_err_right_both_after"],
                 evals_before=s["mean_evals_before"], evals_after=s["mean_evals_after"],
                 basin_miss_before=s.get("basin_miss_before"),
                 basin_miss_after=s["basin_miss_after"])
        )  # fmt: skip
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / f"summary_{arm}.json").write_text(json.dumps(summ, indent=1, default=float))
    tbl = doc_tables.markdown_table(
        rows,
        ["variant", "n", "wrong_before", "wrong_after", "broken", "fixed", "both_wrong",
         "med_both_before", "med_both_after", "evals_before", "evals_after",
         "basin_miss_before", "basin_miss_after"],
        formats={"med_both_before": ".4f", "med_both_after": ".4f",
                 "evals_before": ".0f", "evals_after": ".0f"},
    )  # fmt: skip
    doc_tables.write_tables(OUT / f"tables_{arm}.md", {f"cpprule_e0_{arm}": tbl})
    print(json.dumps(summ, indent=1, default=float))


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    ap.add_argument("cmd", choices=["pilot", "run", "summary"])
    ap.add_argument("--workers", type=int, default=10)
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument(
        "--arm",
        default="on",
        choices=["on", "off", "rule", "rayR2", "rayA", "edge", "interp", "lenint"],
    )
    ap.add_argument("--level0-only", action="store_true")
    ap.add_argument(
        "--msin", type=float, default=0.3, help="min_sin_eta of the cost function (E0: 0.3)"
    )
    ap.add_argument("--tag", default="", help="result dir suffix (centroid task: 'centroid')")
    ap.add_argument(
        "--variant",
        default="",
        choices=[""] + C.ALL_VARIANTS,
        help="only this image variant (clean, all, realistic_nodrop; default: clean and all)",
    )
    ap.add_argument("--meta-dir", default=str(OUT), help="where run_meta_*.json is written")
    a = ap.parse_args()
    if a.cmd == "summary":
        summary(a.arm)
        return
    suffix = "" if a.msin == 0.3 else "_ms0"
    tag = f"_{a.tag}" if a.tag else ""
    root = CACHE / ("pilot_" + a.arm if a.cmd == "pilot" else a.arm + suffix + tag)
    its = work_items(
        root, a.arm, a.level0_only, a.limit or (4 if a.cmd == "pilot" else 0), a.msin, a.variant
    )
    its = [it for it in its if not Path(it[3]).exists()]
    print(f"{len(its)} tasks, {a.workers} workers", flush=True)
    t0 = time.time()
    C.run_pool(task, its, a.workers, f"e0_{a.arm}")
    meta = Path(a.meta_dir)
    meta.mkdir(parents=True, exist_ok=True)
    vtag = f"_{a.variant}" if a.variant else ""
    (meta / f"run_meta_{a.cmd}_{a.arm}{suffix}{tag}{vtag}.json").write_text(
        json.dumps(dict(workers=a.workers, wall_s=time.time() - t0, tasks=len(its),
                        level0_only=a.level0_only, arm=a.arm, msin=a.msin, variant=a.variant,
                        note="contended multi-worker wall, not a timing claim"))
    )  # fmt: skip


if __name__ == "__main__":
    main()
