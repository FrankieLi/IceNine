"""Summarize the Phase D seed-cost diagnosis: reads benchmarks/phase_d_seed_diag/runs/*.json (and
the
variance traces, rank analysis, seed-count oracle and the C++ reference JSON if present) and writes
summary.json and tables.md (doc_tables blocks).

Usage (from icenine_py/): uv run python scripts/phase_d/seed_diag_summary.py
"""

import glob
import json
import sys
from pathlib import Path
from typing import Any, Dict, List

import numpy as np

HERE = Path(__file__).resolve().parent
ICE = HERE.parents[1]
OUT = ICE / "benchmarks" / "phase_d_seed_diag"
sys.path.insert(0, str(ICE / "scripts" / "common"))
import doc_tables as DT  # noqa: E402
import stats as ST  # noqa: E402


def load(pat: str) -> List[Dict[str, Any]]:
    return [json.load(open(f)) for f in sorted(glob.glob(str(OUT / "runs" / pat)))]


def stage_evals(d: Dict[str, Any]) -> Dict[str, float]:
    bp = d["by_phase"]
    out: Dict[str, float] = dict(
        disc_L0=0, quick_L0=0, disc_L1_3=0, quick_L1_3=0, find=0, variance=0
    )
    for k, v in bp.items():
        ph = k.split("|")[0]
        if ph.startswith("discrete_L"):
            out["disc_L0" if ph.endswith("L0") else "disc_L1_3"] += v["n"]
        elif ph.startswith("quick_L"):
            out["quick_L0" if ph.endswith("L0") else "quick_L1_3"] += v["n"]
        elif ph in ("find", "variance"):
            out[ph] += v["n"]
    return out


def main() -> None:
    summ: Dict[str, Any] = {}
    tables: Dict[str, str] = {}

    base = load("base_full_v*.json")
    iso = load("base_isolated_v*.json")
    rows = []
    for d in base:
        s = stage_evals(d)
        fo = d["find_optimal"]
        rows.append(
            dict(
                voxel=d["voxel"],
                wall_s=d["wall_s"],
                evals=d["evals_total"],
                disc_L0=s["disc_L0"],
                quick_L0=s["quick_L0"],
                rest_L1_3=s["disc_L1_3"] + s["quick_L1_3"],
                find=s["find"],
                variance=s["variance"],
                var_share=s["variance"] / d["evals_total"],
                us=d["mean_us_per_eval"],
                find_n=f"{fo['n_evaluated']}/{fo['n_candidates']}",
                conv=fo["converged"],
                err=d["err_deg"],
                q_true=d["q_true"],
            )
        )
    if rows:
        mean = {
            k: float(np.mean([r[k] for r in rows]))
            for k in (
                "wall_s",
                "evals",
                "disc_L0",
                "quick_L0",
                "rest_L1_3",
                "find",
                "variance",
                "var_share",
                "us",
            )
        }
        summ["base_full_mean"] = mean
        summ["base_full_rows"] = rows
        rows_t = rows + [
            dict(voxel="mean", **mean, find_n="", conv="", err=float("nan"), q_true=float("nan"))
        ]
        tables["seed_diag_stage_breakdown"] = DT.markdown_table(
            rows_t,
            [
                "voxel",
                "wall_s",
                "evals",
                "disc_L0",
                "quick_L0",
                "rest_L1_3",
                "find",
                "variance",
                "var_share",
                "us",
                "find_n",
                "conv",
                "err",
                "q_true",
            ],
            formats=dict(
                wall_s=".0f",
                evals=",.0f",
                disc_L0=",.0f",
                quick_L0=",.0f",
                rest_L1_3=",.0f",
                find=",.0f",
                variance=",.0f",
                var_share=".2f",
                us=".0f",
                err=".3f",
                q_true=".3f",
            ),
        )
        # levels table (mean over voxels)
        lv_rows = []
        for L in range(4):
            sel = [d["levels"][L] for d in base if len(d["levels"]) > L]
            lv_rows.append(
                dict(
                    level=L,
                    n_fz=np.mean([x["n_fz"] for x in sel]),
                    n_returned=np.mean([x["n_returned"] for x in sel]),
                    n_kept=np.mean([x["n_kept"] for x in sel]),
                    discrete_evals=np.mean([x["discrete_evals"] for x in sel]),
                    n_recip=sel[0]["n_recip"],
                    diameter_deg=sel[0]["diameter_deg"],
                )
            )
        tables["seed_diag_levels"] = DT.markdown_table(
            lv_rows,
            ["level", "n_fz", "n_returned", "n_kept", "discrete_evals", "n_recip", "diameter_deg"],
            formats=dict(
                n_fz=",.0f",
                n_returned=",.0f",
                n_kept=",.0f",
                discrete_evals=",.0f",
                diameter_deg=".2f",
            ),
        )
        # pass rate of the peak_overlap>0 screen at level 0
        pr = []
        for d in base + iso:
            v = d["by_phase"]["discrete_L0|pr3"]
            pr.append(
                dict(
                    images=d["images"],
                    voxel=d["voxel"],
                    screened=v["n"],
                    passed=v["passed"],
                    pass_rate=v["passed"] / v["n"],
                    n_returned=d["levels"][0]["n_returned"],
                )
            )
        summ["screen_pass"] = pr
        tables["seed_diag_screen"] = DT.markdown_table(
            pr,
            ["images", "voxel", "screened", "passed", "pass_rate", "n_returned"],
            formats=dict(screened=",d", passed=",d", pass_rate=".3f", n_returned=",d"),
        )
    if iso:
        irows = []
        for d in iso:
            b = next((x for x in base if x["voxel"] == d["voxel"]), None)
            irows.append(
                dict(
                    voxel=d["voxel"],
                    iso_evals=d["evals_total"],
                    iso_wall_s=d["wall_s"],
                    iso_us=d["mean_us_per_eval"],
                    iso_err=d["err_deg"],
                    iso_var_evals=d["variance_evals"],
                    full_evals=b["evals_total"] if b else None,
                    full_wall_s=b["wall_s"] if b else None,
                    full_us=b["mean_us_per_eval"] if b else None,
                    full_err=b["err_deg"] if b else None,
                    evals_ratio=(b["evals_total"] / d["evals_total"]) if b else None,
                )
            )
        summ["isolated_vs_full"] = irows
        tables["seed_diag_isolated_vs_full"] = DT.markdown_table(
            irows,
            [
                "voxel",
                "iso_evals",
                "iso_wall_s",
                "iso_us",
                "iso_err",
                "iso_var_evals",
                "full_evals",
                "full_wall_s",
                "full_us",
                "full_err",
                "evals_ratio",
            ],
            formats=dict(
                iso_evals=",d",
                iso_wall_s=".0f",
                iso_us=".0f",
                iso_err=".3f",
                iso_var_evals=",d",
                full_evals=",d",
                full_wall_s=".0f",
                full_us=".0f",
                full_err=".3f",
                evals_ratio=".1f",
            ),
        )

    # options (interior voxels)
    opt_rows = []
    for tag, label in (
        ("base", "base (mc, as Phase D config)"),
        ("cap2k", "variance cap 2000 steps"),
        ("cma", "CMA-ES finisher (LocalOptimizer cma)"),
        ("lean", "cap 2000 + top 3000 at level 0 + keep 1/8"),
    ):
        ds = load(f"{tag}_full_v*.json")
        ds = [d for d in ds if d["voxel"] in {r["voxel"] for r in rows}]
        if not ds:
            continue
        errs = np.array([d["err_deg"] for d in ds])
        wrong = int((errs > 1.0).sum())
        lo, hi = ST.wilson(wrong, len(ds))
        opt_rows.append(
            dict(
                option=label,
                n=len(ds),
                wall_s=np.mean([d["wall_s"] for d in ds]),
                evals=np.mean([d["evals_total"] for d in ds]),
                err_med=float(np.median(errs)),
                err_max=float(errs.max()),
                wrong=f"{wrong}/{len(ds)} [{lo:.2f}, {hi:.2f}]",
            )
        )
    if opt_rows:
        summ["options_interior"] = opt_rows
        tables["seed_diag_options_interior"] = DT.markdown_table(
            opt_rows,
            ["option", "n", "wall_s", "evals", "err_med", "err_max", "wrong"],
            formats=dict(wall_s=".0f", evals=",.0f", err_med=".3f", err_max=".3f"),
        )
    # random voxels
    rrows = []
    for tag, label in (
        ("rcap2k", "variance cap 2000"),
        ("rlean", "lean (cap 2000, top 3000, keep 1/8)"),
    ):
        ds = load(f"{tag}_full_v*.json")
        if not ds:
            continue
        errs = np.array([d["err_deg"] for d in ds])
        wrong = int((errs > 1.0).sum())
        lo, hi = ST.wilson(wrong, len(ds))
        rrows.append(
            dict(
                option=label,
                n=len(ds),
                wall_s=np.mean([d["wall_s"] for d in ds]),
                wall_max=max(d["wall_s"] for d in ds),
                evals=np.mean([d["evals_total"] for d in ds]),
                err_med=float(np.median(errs)),
                err_max=float(errs.max()),
                wrong=f"{wrong}/{len(ds)} [{lo:.2f}, {hi:.2f}]",
            )
        )
        summ[f"random_{tag}"] = [
            dict(
                voxel=d["voxel"],
                wall_s=d["wall_s"],
                evals=d["evals_total"],
                err=d["err_deg"],
                q_true=d["q_true"],
                hit=d["hit_ratio_final"],
            )
            for d in ds
        ]
    if rrows:
        summ["options_random"] = rrows
        tables["seed_diag_options_random"] = DT.markdown_table(
            rrows,
            ["option", "n", "wall_s", "wall_max", "evals", "err_med", "err_max", "wrong"],
            formats=dict(wall_s=".0f", wall_max=".0f", evals=",.0f", err_med=".3f", err_max=".3f"),
        )
        if len(summ.get("random_rcap2k", [])) and len(summ.get("random_rlean", [])):
            a = {d["voxel"]: d["err"] for d in summ["random_rcap2k"]}
            b = {d["voxel"]: d["err"] for d in summ["random_rlean"]}
            common = sorted(set(a) & set(b))
            wa = [a[v] > 1.0 for v in common]
            wb = [b[v] > 1.0 for v in common]
            bb, cc = ST.paired_discordant(wa, wb)
            summ["paired_cap2k_vs_lean_random"] = dict(
                n=len(common),
                wrong_cap_only=bb,
                wrong_lean_only=cc,
                mcnemar_p=ST.mcnemar_exact(bb, cc),
            )
    # variance traces, ranks, oracle
    vt = []
    for f in sorted(glob.glob(str(OUT / "variance_trace_v*.json"))):
        vt.append(json.load(open(f)))
    if vt:
        summ["variance_trace"] = vt
        vrows = []
        for t in vt:
            for images in ("full", "isolated"):
                for box, r in t[images].items():
                    vrows.append(
                        dict(
                            voxel=t["voxel"],
                            images=images,
                            box_deg=float(box),
                            n_runs=r["n_runs"],
                            ended=r["ended"],
                            p_low=r["frac_runs_var_below_thr"],
                            var_median=r["var_quantiles"][2],
                        )
                    )
        tables["seed_diag_variance_trace"] = DT.markdown_table(
            vrows,
            ["voxel", "images", "box_deg", "ended", "n_runs", "p_low", "var_median"],
            formats=dict(box_deg=".2f", p_low=".4f", var_median=".4f", n_runs=",d"),
        )
    rk = []
    for f in sorted(glob.glob(str(OUT / "rank_analysis_*.json"))):
        rk += [r for r in json.load(open(f)) if r["near_deg"] == 3.0]
    if rk:
        summ["rank_analysis_3deg"] = rk
        tables["seed_diag_ranks"] = DT.markdown_table(
            rk,
            [
                "voxel",
                "n",
                "n_near",
                "best_near_rank_discrete",
                "best_near_rank_after_quickmc",
                "n_near_in_kept_quarter",
            ],
            formats=dict(n=",d"),
        )
    orc = OUT / "seed_count_oracle.json"
    if orc.exists():
        o = json.load(open(orc))
        summ["seed_count_oracle"] = o
        tables["seed_diag_seed_count"] = DT.markdown_table(
            [
                dict(p_fail=k, mean=v["mean"], lo=v["min"], hi=v["max"])
                for k, v in o["oracle_seeds_by_p_fail"].items()
            ],
            ["p_fail", "mean", "lo", "hi"],
            formats=dict(mean=".0f"),
        )
    cpp = OUT / "cpp_reference.json"
    if cpp.exists():
        c = json.load(open(cpp))
        summ["cpp_reference"] = c
        tables["seed_diag_cpp"] = DT.markdown_table(
            c["rows"],
            [
                "config",
                "images",
                "voxel",
                "n_runs",
                "adap_s_mean",
                "adap_evals_mean",
                "us_per_eval",
                "variance_steps_median",
                "variance_steps_max",
            ],
            formats=dict(adap_s_mean=".1f", adap_evals_mean=",.0f", us_per_eval=".1f"),
        )
    est = OUT / "run_estimates.json"
    if est.exists():
        e = json.load(open(est))
        summ["run_estimates"] = e
        tables["seed_diag_run_estimates"] = DT.markdown_table(
            e["rows"],
            [
                "option",
                "per_seed_s",
                "seeds_low",
                "seeds_mid",
                "seeds_high",
                "hours_low",
                "hours_mid",
                "hours_high",
            ],
            formats=dict(per_seed_s=".0f", hours_low=".1f", hours_mid=".1f", hours_high=".1f"),
        )
    tim = OUT / "timing_single_worker.json"
    if tim.exists():
        summ["timing_single_worker"] = json.load(open(tim))
    (OUT / "summary.json").write_text(json.dumps(summ, indent=1, default=float) + "\n")
    DT.write_tables(OUT / "tables.md", tables)
    print("wrote", OUT / "summary.json", OUT / "tables.md", sorted(tables))


if __name__ == "__main__":
    main()
