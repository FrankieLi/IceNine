#!/usr/bin/env python3
"""Phase C: no-start reconstruct_voxel on the 200-voxel E0 set, local_optimizer "mc" vs "cma".

Same images (cache of scripts/findoptimal_robustness) and rng as E0 seed 0
(`C.run_seed(vpos, 0)`), so the "mc" run must reproduce the stored E0 `s0_R_final` bit for bit
(checked in `summary`) and the "cma" run differs only through the refinement stage (FindOptimal's
MC per candidate + VarianceMinimizing replaced by one CMAOptimizer run per candidate,
sigma0 0.2 deg, 1000 evaluations, the SearchParameters defaults). The two methods of a
(voxel, variant) run back to back in one worker, in alternating order (voxel position parity).

  uv run python scripts/cma_finisher/e0_cma.py pilot
  uv run python scripts/cma_finisher/e0_cma.py run --workers 10
  uv run python scripts/cma_finisher/e0_cma.py timing           # single worker, preflight-gated
  uv run python scripts/cma_finisher/e0_cma.py summary
"""

import argparse
import json
import sys
import time
from pathlib import Path
from typing import Any, Dict, List, Tuple

import numpy as np

HERE = Path(__file__).resolve().parent
ICENINE_PY = HERE.parents[1]
sys.path.insert(0, str(ICENINE_PY / "scripts" / "findoptimal_robustness"))
sys.path.insert(0, str(ICENINE_PY / "scripts" / "common"))
import common as C  # noqa: E402
import doc_tables  # noqa: E402
import stats as S  # noqa: E402

CACHE = HERE / "cache"
OUT = ICENINE_PY / "benchmarks" / "cma_finisher"
RUN = "e0_run"  # --tag T: cache/e0_run_T and benchmarks/cma_finisher/T/ (rerun, old cache kept)
METHODS = ("mc", "cma")
N_TIMING_VOXELS = 20  # the first 20 of the 200 voxels (U0 convention)


def one_run(W: Any, vctx: Any, vpos: int, method: str) -> Dict[str, Any]:
    rec = W.rec
    rec.params.local_optimizer = method
    try:
        t0 = time.perf_counter()
        with C.quiet():
            res = rec.reconstruct_voxel(vctx.vertices, vctx.voxel.phase, rng=C.run_seed(vpos, 0))
        dt = time.perf_counter() - t0
    finally:
        rec.params.local_optimizer = "mc"
    g, loc, _ = rec.last_eval_counts
    fo = rec.last_find_optimal
    return dict(
        R=np.asarray(res.orientation, dtype=np.float64), cost=float(res.cost), runtime=dt,
        evals_global=g, evals_local=loc, n_find_runs=int(fo["n_evaluated"]),
        converged=bool(fo["converged"]),
    )  # fmt: skip


def task(item: Tuple[Any, ...]) -> str:
    from optimizer_sweep import voxel_context

    vidx, vpos, variant, path = item
    W = C.get_worker()
    t0 = time.time()
    keys = np.load(C.CACHE_DIR / "images" / f"v{vidx}_{variant}.npz")["keys"]
    C.attach(keys)
    vctx = voxel_context(W.ctx, vidx)
    order = METHODS if vpos % 2 == 0 else METHODS[::-1]
    out: Dict[str, Any] = dict(R_true=np.asarray(vctx.R_true, dtype=np.float64), order=list(order))
    for m in order:
        r = one_run(W, vctx, vpos, m)
        for k, v in r.items():
            out[f"{m}_{k}"] = v
    np.savez_compressed(path, **out)
    return f"voxel {vidx} {variant} {time.time() - t0:.0f}s"


def work_items(root: Path, limit: int = 0) -> List[Tuple[Any, ...]]:
    info = dict(np.load(C.OUT_DIR / "voxels.npz"))
    vox = [int(v) for v in info["voxel_indices"]][: C.N_VOXELS]
    root.mkdir(parents=True, exist_ok=True)
    out = []
    for vpos, v in enumerate(vox):
        for var in C.VARIANTS:
            if not (C.CACHE_DIR / "e0" / f"v{v}_{var}.npz").exists():
                continue  # unbuildable
            out.append((v, vpos, var, str(root / f"v{v}_{var}.npz")))
    return out[:limit] if limit else out


def load(root: Path, items: List[Tuple[Any, ...]]) -> Dict[str, Dict[str, List[Any]]]:
    d: Dict[str, Dict[str, List[Any]]] = {}
    for vidx, vpos, var, path in items:
        if not Path(path).exists():
            continue
        z = np.load(path)
        e0 = np.load(C.CACHE_DIR / "e0" / f"v{vidx}_{var}.npz")
        rows = d.setdefault(var, {})
        rows.setdefault("vidx", []).append(vidx)
        rows.setdefault("e0_R", []).append(e0["s0_R_final"])
        rows.setdefault("e0_cost", []).append(float(e0["s0_cost_final"]))
        rows.setdefault("e0_evals", []).append(int(e0["s0_evals_global"] + e0["s0_evals_local"]))
        rows.setdefault("R_true", []).append(z["R_true"])
        for m in METHODS:
            for k in ("R", "cost", "runtime", "evals_global", "evals_local", "n_find_runs"):
                rows.setdefault(f"{m}_{k}", []).append(z[f"{m}_{k}"])
    return d


def summarize(root: Path, items: List[Tuple[Any, ...]]) -> Tuple[Dict[str, Any], Dict[str, str]]:
    data = load(root, items)
    summ: Dict[str, Any] = {}
    rows_tbl: List[Dict[str, Any]] = []
    rows_ev: List[Dict[str, Any]] = []
    for var in C.VARIANTS:
        if var not in data:
            continue
        d = data[var]
        n = len(d["vidx"])
        Rt = np.stack(d["R_true"])
        err = {m: C.err_deg(np.stack(d[f"{m}_R"]), Rt) for m in METHODS}
        wrong = {m: err[m] > C.WRONG_DEG for m in METHODS}
        label = "realistic" if var == "all" else var
        ident = int(
            sum(bool((np.asarray(a) == np.asarray(b)).all()) for a, b in zip(d["mc_R"], d["e0_R"]))
        )
        b, c = S.paired_discordant(wrong["mc"], wrong["cma"])  # b: mc only wrong, c: cma only wrong
        s: Dict[str, Any] = dict(
            n=n, mc_identical_to_stored_e0=ident, mc_only_wrong=b, cma_only_wrong=c,
            mcnemar_p=S.mcnemar_exact(b, c),
        )  # fmt: skip
        both_right = ~wrong["mc"] & ~wrong["cma"]
        s["n_right_in_both"] = int(both_right.sum())
        for m in METHODS:
            k = int(wrong[m].sum())
            lo, hi = S.wilson(k, n)
            right = ~wrong[m]
            ev = np.asarray(d[f"{m}_evals_global"]) + np.asarray(d[f"{m}_evals_local"])
            rt = np.asarray(d[f"{m}_runtime"])
            s[m] = dict(
                n_wrong=k, wrong_rate=k / n, wilson_lo=lo, wilson_hi=hi,
                median_err_right=float(np.median(err[m][right])),
                median_err_right_both=float(np.median(err[m][both_right])),
                mean_evals_total=float(ev.mean()),
                mean_evals_local=float(np.mean(d[f"{m}_evals_local"])),
                mean_evals_global=float(np.mean(d[f"{m}_evals_global"])),
                mean_find_runs=float(np.mean(d[f"{m}_n_find_runs"])),
                mean_runtime_s=float(rt.mean()), median_runtime_s=float(np.median(rt)),
                mean_cost=float(np.mean(d[f"{m}_cost"])),
            )  # fmt: skip
            lt = f"{label}"
            rows_tbl.append(
                dict(variant=lt, method=m, n=n, wrong=f"{k}/{n} = {k / n:.3f} [{lo:.3f}, {hi:.3f}]",
                     med_err=s[m]["median_err_right"], med_both=s[m]["median_err_right_both"],
                     cost=s[m]["mean_cost"])
            )  # fmt: skip
            rows_ev.append(
                dict(variant=lt, method=m, glob=s[m]["mean_evals_global"],
                     loc=s[m]["mean_evals_local"], find=s[m]["mean_find_runs"],
                     wall=s[m]["mean_runtime_s"], wall_med=s[m]["median_runtime_s"])
            )  # fmt: skip
        # win rate of the final cost, paired: cma lower / higher / equal
        dc = np.asarray(d["cma_cost"]) - np.asarray(d["mc_cost"])
        s["evals_total_change_pct"] = 100.0 * (
            s["cma"]["mean_evals_total"] / s["mc"]["mean_evals_total"] - 1.0
        )
        s["cost_cma_lower"], s["cost_cma_higher"] = int((dc < 0).sum()), int((dc > 0).sum())
        s["cost_equal"] = int((dc == 0).sum())
        s["wall_ratio_cma_over_mc_median"] = float(
            np.median(np.asarray(d["cma_runtime"]) / np.asarray(d["mc_runtime"]))
        )
        summ[var] = s
    tables = {
        "e0_accuracy": doc_tables.markdown_table(
            rows_tbl, ["variant", "method", "n", "wrong", "med_err", "med_both", "cost"],
            formats={"med_err": ".4f", "med_both": ".4f", "cost": ".4f"},
        ),
        "e0_evals_wall": doc_tables.markdown_table(
            rows_ev, ["variant", "method", "glob", "loc", "find", "wall", "wall_med"],
            formats={"glob": ".0f", "loc": ".0f", "find": ".2f", "wall": ".2f", "wall_med": ".2f"},
        ),
    }  # fmt: skip
    return summ, tables


def compare(old_root: Path, new_root: Path) -> Tuple[Dict[str, Any], Dict[str, str]]:
    """Before / after the C++-faithful MC and VarianceMinimizing: the same 200 voxels x variants,
    paired by voxel. Wrong = symmetry-reduced misorientation over C.WRONG_DEG."""
    items_old, items_new = work_items(old_root), work_items(new_root)

    def tot(d: Dict[str, List[Any]], m: str) -> np.ndarray:
        return np.asarray(d[f"{m}_evals_global"]) + np.asarray(d[f"{m}_evals_local"])

    d_old, d_new = load(old_root, items_old), load(new_root, items_new)
    out: Dict[str, Any] = {}
    rows: List[Dict[str, Any]] = []
    for var in C.VARIANTS:
        o, n = d_old[var], d_new[var]
        assert o["vidx"] == n["vidx"]
        Rt = np.stack(n["R_true"])
        label = "realistic" if var == "all" else var
        res: Dict[str, Any] = {}
        for m in METHODS:
            eo = C.err_deg(np.stack(o[f"{m}_R"]), Rt)
            en = C.err_deg(np.stack(n[f"{m}_R"]), Rt)
            wo, wn = eo > C.WRONG_DEG, en > C.WRONG_DEG
            b, c = S.paired_discordant(wn, wo)  # b: wrong only after, c: wrong only before
            right_both = ~wo & ~wn
            res[m] = dict(
                wrong_before=int(wo.sum()), wrong_after=int(wn.sum()), only_after_wrong=b,
                only_before_wrong=c, mcnemar_p=S.mcnemar_exact(b, c),
                med_err_right_both_before=float(np.median(eo[right_both])),
                med_err_right_both_after=float(np.median(en[right_both])),
                n_right_both=int(right_both.sum()),
                mean_evals_before=float(np.mean(tot(o, m))),
                mean_evals_after=float(np.mean(tot(n, m))),
                mean_local_before=float(np.mean(o[f"{m}_evals_local"])),
                mean_local_after=float(np.mean(n[f"{m}_evals_local"])),
                mean_cost_before=float(np.mean(o[f"{m}_cost"])),
                mean_cost_after=float(np.mean(n[f"{m}_cost"])),
            )  # fmt: skip
            lo_b, hi_b = S.wilson(res[m]["wrong_before"], len(eo))
            lo_a, hi_a = S.wilson(res[m]["wrong_after"], len(eo))
            rows.append(
                dict(variant=label, method=m, n=len(eo),
                     wrong_before=f"{res[m]['wrong_before']}/{len(eo)} [{lo_b:.3f}, {hi_b:.3f}]",
                     wrong_after=f"{res[m]['wrong_after']}/{len(eo)} [{lo_a:.3f}, {hi_a:.3f}]",
                     only_after=b, only_before=c, p=res[m]["mcnemar_p"],
                     med_both_before=res[m]["med_err_right_both_before"],
                     med_both_after=res[m]["med_err_right_both_after"],
                     evals_before=res[m]["mean_evals_before"],
                     evals_after=res[m]["mean_evals_after"])
            )  # fmt: skip
        out[var] = res
    table = doc_tables.markdown_table(
        rows,
        ["variant", "method", "n", "wrong_before", "wrong_after", "only_after", "only_before", "p",
         "med_both_before", "med_both_after", "evals_before", "evals_after"],
        formats={"p": ".2g", "med_both_before": ".4f", "med_both_after": ".4f",
                 "evals_before": ".0f", "evals_after": ".0f"},
    )  # fmt: skip
    return out, {"e0_before_after": table}


def validate_table(v: Dict[str, Any]) -> str:
    rows = []
    for key, label in (
        ("start", "start (net x3)"),
        ("b3_250", "B3 cma_02, 250 evals"),
        ("lib250", "library CMAOptimizer, 250 evals"),
        ("b3_1000", "B3 cma_02, 1000 evals"),
        ("lib1000", "library CMAOptimizer, 1000 evals"),
        ("refine", "refine_from_candidates, local_optimizer cma (1000 + final)"),
    ):
        r = v[key]
        frac = f"{r['n_below_0p02']}/{v['n_cases']} [{r['wilson_lo']:.2f}, {r['wilson_hi']:.2f}]"
        rows.append(dict(method=label, n=v["n_cases"], med=r["median_err_deg"], under=frac))
    for mode, label in (
        ("mc", "local_optimization, mc (VarianceMinimizing)"),
        ("cma", "local_optimization, cma"),
    ):
        r = v[f"local_optimization_{mode}"]
        frac = f"{r['n_below_0p02']}/{v['n_cases']} [{r['wilson_lo']:.2f}, {r['wilson_hi']:.2f}]"
        rows.append(
            dict(method=label, n=v["n_cases"], med=r["median_err_deg"], under=frac,
                 evals=r["median_evals"])
        )  # fmt: skip
    return (
        doc_tables.markdown_table(
            rows, ["method", "n", "med", "under", "evals"], formats={"med": ".4f", "evals": ".0f"}
        )
        .replace("| med |", "| median error (deg) |")
        .replace("| under |", "| fraction < 0.02 deg (Wilson 95%) |")
        .replace("| evals |", "| median evals (local_optimization rows) |")
    )


def timing_summary(items: List[Tuple[Any, ...]]) -> Tuple[Dict[str, Any], str]:
    out: Dict[str, Any] = {}
    rows = []
    for var in C.VARIANTS:
        sel = [it for it in items if it[2] == var]
        if not sel:
            continue
        z = [np.load(it[3]) for it in sel]
        mc = np.array([float(x["mc_runtime"]) for x in z])
        cm = np.array([float(x["cma_runtime"]) for x in z])
        ev = lambda x, m: int(x[f"{m}_evals_global"]) + int(x[f"{m}_evals_local"])  # noqa: E731
        label = "realistic" if var == "all" else var
        out[var] = dict(
            n=len(sel), mc_mean_s=float(mc.mean()), cma_mean_s=float(cm.mean()),
            ratio_of_means=float(cm.mean() / mc.mean()),
            time_change_pct=float(100.0 * (cm.mean() / mc.mean() - 1.0)),
            median_paired_ratio=float(np.median(cm / mc)),
            mc_mean_evals=float(np.mean([ev(x, "mc") for x in z])),
            cma_mean_evals=float(np.mean([ev(x, "cma") for x in z])),
        )  # fmt: skip
        rows.append(dict(variant=label, **out[var]))
    table = doc_tables.markdown_table(
        rows,
        ["variant", "n", "mc_mean_s", "cma_mean_s", "ratio_of_means", "median_paired_ratio",
         "mc_mean_evals", "cma_mean_evals"],
        formats={"mc_mean_s": ".2f", "cma_mean_s": ".2f", "ratio_of_means": ".3f",
                 "median_paired_ratio": ".3f", "mc_mean_evals": ".0f", "cma_mean_evals": ".0f"},
    )  # fmt: skip
    return out, table


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    ap.add_argument("cmd", choices=["pilot", "run", "timing", "summary", "compare"])
    ap.add_argument("--workers", type=int, default=10)
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--tag", default="")
    a = ap.parse_args()
    global OUT, RUN
    if a.tag:
        OUT, RUN = OUT / a.tag, f"e0_run_{a.tag}"
    OUT.mkdir(parents=True, exist_ok=True)
    if a.cmd == "summary":
        items = work_items(CACHE / RUN)
        summ, tables = summarize(CACHE / RUN, items)
        merged: Dict[str, Any] = dict(e0=summ)
        vl = OUT / "validate_lib.json"
        if vl.exists():
            merged["validate_lib"] = json.loads(vl.read_text())
            tables["validate_lib"] = validate_table(merged["validate_lib"])
        tdir = CACHE / "timing"  # single-worker timing exists for the original run only
        tit = [it for it in work_items(tdir) if Path(it[3]).exists()] if not a.tag else []
        if tit:
            merged["timing_single_worker"], tables["timing_single_worker"] = timing_summary(tit)
        if a.tag:  # distinct marker names: sync_doc_tables must not overwrite the Phase C blocks
            tables = {f"mcf_{k}": v for k, v in tables.items()}
        (OUT / "summary.json").write_text(json.dumps(merged, indent=1))
        doc_tables.write_tables(OUT / "tables.md", tables)
        print(json.dumps(merged, indent=1))
        return
    if a.cmd == "compare":
        res, tbl = compare(CACHE / "e0_run", CACHE / RUN)
        (OUT / "compare.json").write_text(json.dumps(res, indent=1))
        doc_tables.write_tables(OUT / "tables_compare.md", tbl)
        print(json.dumps(res, indent=1))
        return
    if a.cmd == "pilot":
        pdir = CACHE / (f"pilot_{a.tag}" if a.tag else "pilot")
        root, its = pdir, work_items(pdir, a.limit or 8)
        its = [it for it in its if not Path(it[3]).exists()]
        t0 = time.time()
        C.run_pool(task, its, a.workers, "pilot")
        print(f"pilot wall {time.time() - t0:.0f}s, {len(its)} tasks", flush=True)
        summ, _ = summarize(root, work_items(root, a.limit or 8))
        print(json.dumps(summ, indent=1))
        return
    if a.cmd == "run":
        its = [it for it in work_items(CACHE / RUN, a.limit) if not Path(it[3]).exists()]
        print(f"{len(its)} tasks, {a.workers} workers", flush=True)
        t0 = time.time()
        C.run_pool(task, its, a.workers, "e0cma")
        (OUT / "run_meta_e0.json").write_text(
            json.dumps(dict(workers=a.workers, wall_s=time.time() - t0, tasks=len(its),
                            note="contended 10-worker wall, not a timing claim"))
        )  # fmt: skip
        return
    # timing: one process, preflight-gated, mc and cma interleaved (alternating order)
    import preflight as P

    info = P.require_quiet()
    (OUT / "preflight_e0_timing.json").write_text(json.dumps(info, indent=1))
    root = CACHE / "timing"
    its = [it for it in work_items(root) if it[1] < N_TIMING_VOXELS]
    its = [it for it in its if not Path(it[3]).exists()][: a.limit or None]
    C.init_worker(C.worker_args())
    t0 = time.time()
    for it in its:
        print(task(it), f"({time.time() - t0:.0f}s)", flush=True)
    info2 = P.preflight()
    (OUT / "preflight_e0_timing_after.json").write_text(json.dumps(info2, indent=1))


if __name__ == "__main__":
    main()
