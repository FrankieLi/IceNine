"""S3 Step 4 summary: timing JSONs (s3_timing.py) and the E0 A/B caches (e0_run.py) to
benchmarks/s3_steps1_5/summary.json and tables.md (doc tables, ``doc_tables.write_tables``).

  uv run python scripts/speed_parity/s3_summary.py --out benchmarks/s3_steps1_5

Inputs found under ``--out``: timing_search_three.json, timing_search_e0.json,
timing_crossover.json. The A/B reads the run caches of scripts/findoptimal_robustness/cache:
``e0`` (key off), ``e0_s3on`` (BatchedQuickMC, CPU) and ``e0_s3mps`` (BatchedQuickMC, SearchDevice
mps); the key-off arm is also compared run by run with the published E0 clean runs
(benchmarks/findoptimal_robustness/e0_runs.npz). No p-values: counts, Wilson intervals and the
paired split (fixed = wrong with the key off and right with it on; broken = the reverse).
"""

import argparse
import json
import statistics
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

import numpy as np

HERE = Path(__file__).resolve().parent
ICENINE_PY = HERE.parents[1]
sys.path.insert(0, str(ICENINE_PY / "scripts" / "common"))
sys.path.insert(0, str(ICENINE_PY / "scripts" / "findoptimal_robustness"))

import doc_tables as DT  # noqa: E402
import stats  # noqa: E402
from csl import reduced_misorientation_deg  # noqa: E402

CACHE = ICENINE_PY / "scripts" / "findoptimal_robustness" / "cache"
PUBLISHED = ICENINE_PY / "benchmarks" / "findoptimal_robustness" / "e0_runs.npz"
WRONG_DEG = 1.0
N_LEVELS = 4


def med(x: Sequence[float]) -> float:
    return float(statistics.median(x)) if len(x) else float("nan")


def table(
    rows: Sequence[Dict[str, Any]], cols: Dict[str, str], fmt: Optional[Dict[str, str]] = None
) -> str:
    """Markdown table of ``rows`` with the columns ``cols`` (row key -> header text)."""
    rel = [{cols[k]: r.get(k) for k in cols} for r in rows]
    return DT.markdown_table(
        rel, list(cols.values()), formats={cols[k]: v for k, v in (fmt or {}).items()}
    )


# ---------------------------------------------------------------------------------------------
# timing
# ---------------------------------------------------------------------------------------------


def summarize_search(path: Path) -> Dict[str, Any]:
    d = json.loads(path.read_text())
    runs = d["runs"]
    out: Dict[str, Any] = dict(
        label=d["label"],
        torch=d["torch"],
        arms=d["arms"],
        repeats=d["repeats"],
        load=d["preflight"].get("loadavg"),
        per_arm={},
        repeat_agreement={},
    )
    arms = d["arms"]
    keyf = lambda r: (r["case"], r["seed"], r["repeat"])  # noqa: E731
    by: Dict[str, Dict[Tuple[str, int, int], Dict[str, Any]]] = {a: {} for a in arms}
    for r in runs:
        by[r["arm"]][keyf(r)] = r
    keys = sorted(by["off"])
    for a in arms:
        rs = [by[a][k] for k in keys if k in by[a]]
        q_us = [r["quick_mc_us_per_eval"] for r in rs]
        out["per_arm"][a] = dict(
            n_runs=len(rs),
            wall_s_median=med([r["wall_s"] for r in rs]),
            wall_s_sum=float(sum(r["wall_s"] for r in rs)),
            quick_mc_s_median=med([r["quick_mc_s"] for r in rs]),
            quick_mc_s_sum=float(sum(r["quick_mc_s"] for r in rs)),
            quick_mc_us_per_eval_median=med(q_us),
            quick_mc_share_of_wall_median=med([r["quick_mc_s"] / r["wall_s"] for r in rs]),
            quick_mc_evals_sum=int(sum(r["quick_mc_evals"] for r in rs)),
            evals_total_median=med([r["evals"]["total"] for r in rs]),
            cand_per_level_median=[
                med([r["levels"][L]["n"] for r in rs if len(r["levels"]) > L])
                for L in range(N_LEVELS)
            ],
            quick_mc_s_per_level_median=[
                med([r["levels"][L]["quick_mc_s"] for r in rs if len(r["levels"]) > L])
                for L in range(N_LEVELS)
            ],
        )
        splits = [r["batched_split"] for r in rs if r.get("batched_split")]
        if splits:
            tot = {k: sum(s["timing"][k] for s in splits) for k in splits[0]["timing"]}
            cnt = {k: sum(s["counts"][k] for s in splits) for k in splits[0]["counts"]}
            out["per_arm"][a]["batched_cost_split_s"] = tot
            out["per_arm"][a]["batched_cost_counts"] = cnt
    off = out["per_arm"]["off"]
    for a in arms:
        if a == "off":
            continue
        pa = out["per_arm"][a]
        pairs = [(by["off"][k], by[a][k]) for k in keys if k in by[a]]
        out["per_arm"][a]["speedup_wall_sum"] = off["wall_s_sum"] / pa["wall_s_sum"]
        out["per_arm"][a]["speedup_wall_median_of_ratios"] = med(
            [o["wall_s"] / n["wall_s"] for o, n in pairs]
        )
        out["per_arm"][a]["speedup_quick_mc_per_eval"] = (
            off["quick_mc_us_per_eval_median"] / pa["quick_mc_us_per_eval_median"]
        )
        out["per_arm"][a]["speedup_quick_mc_sum"] = off["quick_mc_s_sum"] / pa["quick_mc_s_sum"]
    # agreement of the repeats: relative difference of wall and quick-MC time per (case, seed, arm)
    for a in arms:
        diffs_w, diffs_q = [], []
        for (case, seed, rep), r in by[a].items():
            if rep != 0:
                continue
            r1 = by[a].get((case, seed, 1))
            if r1 is None:
                continue
            diffs_w.append(abs(r["wall_s"] - r1["wall_s"]) / (0.5 * (r["wall_s"] + r1["wall_s"])))
            diffs_q.append(
                abs(r["quick_mc_s"] - r1["quick_mc_s"])
                / (0.5 * (r["quick_mc_s"] + r1["quick_mc_s"]))
            )
        out["repeat_agreement"][a] = dict(
            n=len(diffs_w),
            wall_median=med(diffs_w),
            wall_max=max(diffs_w, default=float("nan")),
            quick_mc_median=med(diffs_q),
            quick_mc_max=max(diffs_q, default=float("nan")),
        )
    out["errors_deg_median"] = {
        a: med([r["err_deg"] for r in by[a].values() if r["err_deg"] is not None]) for a in arms
    }
    return out


def search_tables(name: str, s: Dict[str, Any]) -> Dict[str, str]:
    rows = []
    for a in s["arms"]:
        p = s["per_arm"][a]
        rows.append(
            dict(
                arm=a,
                n=p["n_runs"],
                wall=p["wall_s_median"],
                quick=p["quick_mc_s_median"],
                us=p["quick_mc_us_per_eval_median"],
                cand="/".join(f"{x:.0f}" for x in p["cand_per_level_median"]),
                evals=p["evals_total_median"],
                sp_wall=p.get("speedup_wall_sum"),
                sp_q=p.get("speedup_quick_mc_per_eval"),
            )
        )
    t1 = table(
        rows,
        {
            "arm": "arm (off = sequential quick MC)",
            "n": "searches",
            "wall": "wall s, median",
            "quick": "quick MC s, median",
            "us": "quick MC us per evaluation, median",
            "cand": "candidates per level, median",
            "evals": "evaluations per search, median",
            "sp_wall": "wall speedup (sum of searches)",
            "sp_q": "quick-MC us per evaluation speedup",
        },
        {
            "wall": ".2f",
            "quick": ".2f",
            "us": ".0f",
            "evals": ".0f",
            "sp_wall": ".2f",
            "sp_q": ".2f",
        },
    )
    return {f"s3_timing_{name}": t1}


def crossover_tables(d: Dict[str, Any]) -> Dict[str, str]:
    rows = d["rows"]
    t = table(
        rows,
        {
            "B": "batch B",
            "cpu_compact_us_per_orient": "CPU compact us per orientation",
            "cpu_fixed_us_per_orient": "CPU fixed-shape us per orientation",
            "mps_fixed_us_per_orient": "MPS fixed-shape us per orientation",
        },
        {
            c: ".0f"
            for c in (
                "cpu_compact_us_per_orient",
                "cpu_fixed_us_per_orient",
                "mps_fixed_us_per_orient",
            )
        },
    )
    return {"s3_crossover": t}


# ---------------------------------------------------------------------------------------------
# the A/B on E0 clean
# ---------------------------------------------------------------------------------------------


def load_arm(dirname: str, seeds: Sequence[int]) -> Dict[Tuple[int, int], Dict[str, Any]]:
    """{(vpos, seed): run} of cache/<dirname>/v*_clean.npz."""
    runs: Dict[Tuple[int, int], Dict[str, Any]] = {}
    for f in sorted((CACHE / dirname).glob("v*_clean.npz")):
        d = np.load(f)
        if "unbuildable" in d.files:
            continue
        vpos, vidx = int(d["vpos"]), int(d["vidx"])
        R_true = d["R_true"]
        for s in seeds:
            if f"s{s}_R_final" not in d.files:
                continue
            err = float(reduced_misorientation_deg(d[f"s{s}_R_final"], R_true))
            n_cand = []
            for L in range(N_LEVELS):
                k = f"s{s}_L{L}_disc_R"
                n_cand.append(int(d[k].shape[0]) if k in d.files else 0)
            runs[(vpos, s)] = dict(
                vidx=vidx,
                vpos=vpos,
                seed=s,
                err=err,
                wrong=err > WRONG_DEG,
                cost=float(d[f"s{s}_cost_final"]),
                runtime=float(d[f"s{s}_runtime"]),
                evals_global=int(d[f"s{s}_evals_global"]),
                evals_local=int(d[f"s{s}_evals_local"]),
                n_cand=n_cand,
            )
    return runs


def wilson_row(k: int, n: int) -> Dict[str, Any]:
    lo, hi = stats.wilson(k, n)
    return dict(wrong=k, n=n, rate=k / n if n else float("nan"), lo=lo, hi=hi)


def arm_summary(runs: Dict[Tuple[int, int], Dict[str, Any]]) -> Dict[str, Any]:
    seeds = sorted({s for _, s in runs})
    out: Dict[str, Any] = dict(n_runs=len(runs), per_seed={}, voxels_wrong_in_k_seeds={})
    for s in seeds:
        rs = [r for (v, ss), r in runs.items() if ss == s]
        out["per_seed"][str(s)] = wilson_row(sum(r["wrong"] for r in rs), len(rs))
    out["pooled"] = wilson_row(sum(r["wrong"] for r in runs.values()), len(runs))
    vox = sorted({v for v, _ in runs})
    full = [v for v in vox if all((v, s) in runs for s in seeds)]
    hist = {k: 0 for k in range(len(seeds) + 1)}
    for v in full:
        hist[sum(runs[(v, s)]["wrong"] for s in seeds)] += 1
    out["voxels_wrong_in_k_seeds"] = dict(
        n_voxels=len(full), counts={str(k): c for k, c in hist.items()}
    )
    rs = list(runs.values())
    out["evals_local_median"] = med([r["evals_local"] for r in rs])
    out["evals_global_median"] = med([r["evals_global"] for r in rs])
    out["evals_local_mean"] = float(np.mean([r["evals_local"] for r in rs]))
    out["evals_global_mean"] = float(np.mean([r["evals_global"] for r in rs]))
    out["cand_per_level_median"] = [med([r["n_cand"][L] for r in rs]) for L in range(N_LEVELS)]
    out["quick_mc_evals_mean"] = float(np.mean([12 * sum(r["n_cand"]) for r in rs]))
    out["runtime_s_median"] = med([r["runtime"] for r in rs])
    out["runtime_s_mean"] = float(np.mean([r["runtime"] for r in rs]))
    out["cost_median"] = med([r["cost"] for r in rs])
    return out


def paired(
    a: Dict[Tuple[int, int], Dict[str, Any]], b: Dict[Tuple[int, int], Dict[str, Any]]
) -> Dict[str, Any]:
    """Paired split of the runs present in both arms (a = reference, b = the changed arm)."""
    common = sorted(set(a) & set(b))
    wa = np.array([a[k]["wrong"] for k in common])
    wb = np.array([b[k]["wrong"] for k in common])
    broken, fixed = stats.paired_discordant(wb, wa)  # (wrong in b only, wrong in a only)
    out: Dict[str, Any] = dict(
        n_paired=len(common),
        wrong_a=int(wa.sum()),
        wrong_b=int(wb.sum()),
        fixed=fixed,
        broken=broken,
        both_wrong=int((wa & wb).sum()),
        both_right=int((~wa & ~wb).sum()),
        per_seed={},
        flips={},
    )
    for s in sorted({k[1] for k in common}):
        ks = [k for k in common if k[1] == s]
        xa = np.array([a[k]["wrong"] for k in ks])
        xb = np.array([b[k]["wrong"] for k in ks])
        br, fx = stats.paired_discordant(xb, xa)
        out["per_seed"][str(s)] = dict(
            n=len(ks),
            wrong_a=int(xa.sum()),
            wrong_b=int(xb.sum()),
            fixed=fx,
            broken=br,
            both_wrong=int((xa & xb).sum()),
            both_right=int((~xa & ~xb).sum()),
        )
        out["flips"][str(s)] = dict(
            fixed_vpos=[int(k[0]) for k in ks if a[k]["wrong"] and not b[k]["wrong"]],
            broken_vpos=[int(k[0]) for k in ks if not a[k]["wrong"] and b[k]["wrong"]],
        )
    same = [abs(a[k]["err"] - b[k]["err"]) < 1e-3 for k in common]
    out["runs_with_same_answer_1e-3deg"] = int(sum(same))
    out["runs_with_same_eval_counts"] = int(
        sum(
            a[k]["evals_local"] == b[k]["evals_local"]
            and a[k]["evals_global"] == b[k]["evals_global"]
            for k in common
        )
    )
    out["runs_with_same_cost"] = int(sum(a[k]["cost"] == b[k]["cost"] for k in common))
    return out


def published_check(off: Dict[Tuple[int, int], Dict[str, Any]]) -> Optional[Dict[str, Any]]:
    """The key-off arm against the published E0 clean runs, run by run."""
    if not PUBLISHED.exists() or not off:
        return None
    d = np.load(PUBLISHED, allow_pickle=True)
    clean = d["variant"] == 0
    pub = {
        (int(v), int(s)): dict(
            err=float(e), wrong=bool(e > WRONG_DEG), cost=float(c), eg=int(g), el=int(l)
        )
        for v, s, e, c, g, l in zip(
            d["vpos"][clean],
            d["seed"][clean],
            d["final_err"][clean],
            d["cost_final"][clean],
            d["evals_global"][clean],
            d["evals_local"][clean],
        )
    }
    common = sorted(set(off) & set(pub))
    if not common:
        return None
    same_flag = sum(off[k]["wrong"] == pub[k]["wrong"] for k in common)
    same_err = sum(abs(off[k]["err"] - pub[k]["err"]) < 1e-3 for k in common)
    same_evals = sum(
        off[k]["evals_global"] == pub[k]["eg"] and off[k]["evals_local"] == pub[k]["el"]
        for k in common
    )
    return dict(
        n_compared=len(common),
        wrong_here=sum(off[k]["wrong"] for k in common),
        wrong_published=sum(pub[k]["wrong"] for k in common),
        same_wrong_flag=same_flag,
        same_final_error_1e3deg=same_err,
        same_eval_counts=same_evals,
        differing_flag=[list(k) for k in common if off[k]["wrong"] != pub[k]["wrong"]][:20],
    )


ARM_NAMES = {
    "off": "key off",
    "on_cpu": "key on (CPU)",
    "on_mps": "key on (MPS; seed 0, voxels 0-49)",
    "off_2w": "key off (2 workers)",
    "on_cpu_2w": "key on (CPU, 2 workers)",
}


def ab_tables(s: Dict[str, Any]) -> Dict[str, str]:  # s = the "ab" summary
    tabs: Dict[str, str] = {}
    rows, labels = [], []
    for name, a in s["arms"].items():
        for sd, w in a["per_seed"].items():
            rows.append(dict(arm=ARM_NAMES[name], scope=f"seed {sd}", **w))
        rows.append(dict(arm=ARM_NAMES[name], scope="pooled", **a["pooled"]))
    tabs["s3_ab_wrong"] = table(
        rows,
        {
            "arm": "arm",
            "scope": "scope",
            "wrong": "wrong (> 1 deg)",
            "n": "runs",
            "rate": "rate",
            "lo": "Wilson 95% low",
            "hi": "Wilson 95% high",
        },
        {"rate": ".3f", "lo": ".3f", "hi": ".3f"},
    )
    rows = []
    for name, a in s["arms"].items():
        rows.append(
            dict(
                arm=ARM_NAMES[name],
                cand="/".join(f"{x:.0f}" for x in a["cand_per_level_median"]),
                local=a["evals_local_mean"],
                glob=a["evals_global_mean"],
                quick=a["quick_mc_evals_mean"],
                rt=None if name == "on_mps" else a["runtime_s_median"],  # worker counts differ
            )
        )
    tabs["s3_ab_evals"] = table(
        rows,
        {
            "arm": "arm",
            "cand": "candidates per level, median",
            "local": "local evaluations per run, mean",
            "glob": "global evaluations per run, mean",
            "quick": "quick-MC evaluations per run, mean",
            "rt": "run time s, median (10 workers, contended)",
        },
        {"local": ".0f", "glob": ".0f", "quick": ".0f", "rt": ".1f"},
    )
    rows = []
    for name, p in s["paired"].items():
        rows.append(
            dict(
                pair=name,
                scope="all seeds",
                **{
                    k: p[k]
                    for k in (
                        "n_paired",
                        "wrong_a",
                        "wrong_b",
                        "fixed",
                        "broken",
                        "both_wrong",
                        "both_right",
                    )
                },
            )
        )
        for sd, q in p["per_seed"].items():
            rows.append(
                dict(
                    pair=name,
                    scope=f"seed {sd}",
                    n_paired=q["n"],
                    wrong_a=q["wrong_a"],
                    wrong_b=q["wrong_b"],
                    fixed=q["fixed"],
                    broken=q["broken"],
                    both_wrong=q["both_wrong"],
                    both_right=q["both_right"],
                )
            )
    if "timing_2w" in s:
        trows = [
            dict(t, arm=ARM_NAMES.get(n + "_2w", ARM_NAMES[n])) for n, t in s["timing_2w"].items()
        ]
        tabs["s3_ab_timing_2w"] = table(
            trows,
            {
                "arm": "arm",
                "n": "runs",
                "wrong": "wrong (> 1 deg)",
                "runtime_median": "run time s, median",
                "runtime_mean": "run time s, mean",
                "runtime_sum": "run time s, sum",
            },
            {"runtime_median": ".2f", "runtime_mean": ".2f", "runtime_sum": ".1f"},
        )
    if "off_seed_to_seed" in s:
        for name, q in s["off_seed_to_seed"].items():
            rows.append(
                dict(
                    pair=f"key off, {name}",
                    scope="same voxels",
                    n_paired=q["n_paired"],
                    wrong_a=q["wrong_a"],
                    wrong_b=q["wrong_b"],
                    fixed=q["fixed"],
                    broken=q["broken"],
                    both_wrong=q["both_wrong"],
                    both_right=q["both_right"],
                )
            )
    tabs["s3_ab_paired"] = table(
        rows,
        {
            "pair": "pair (a vs b)",
            "scope": "scope",
            "n_paired": "paired runs",
            "wrong_a": "wrong a",
            "wrong_b": "wrong b",
            "fixed": "fixed (wrong a, right b)",
            "broken": "broken (right a, wrong b)",
            "both_wrong": "both wrong",
            "both_right": "both right",
        },
    )
    return tabs


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    ap.add_argument("--out", required=True)
    ap.add_argument("--seeds", type=int, default=3)
    a = ap.parse_args()
    out = Path(a.out)
    summary: Dict[str, Any] = {}
    tables: Dict[str, str] = {}
    for name in ("three", "e0"):
        p = out / f"timing_search_{name}.json"
        if p.exists():
            summary[f"timing_{name}"] = summarize_search(p)
            tables.update(search_tables(name, summary[f"timing_{name}"]))
    p = out / "timing_crossover.json"
    if p.exists():
        d = json.loads(p.read_text())
        summary["crossover"] = dict(
            label=d["label"],
            gpu_min_batch=d["gpu_min_batch"],
            mps_faster_from_B=d["mps_faster_from_B"],
            rows=d["rows"],
        )
        tables.update(crossover_tables(d))
    seeds = list(range(a.seeds))
    arms_dirs = {
        "off": "e0",
        "on_cpu": "e0_s3on",
        "on_mps": "e0_s3mps",
        "off_2w": "e0_s3off2w",
        "on_cpu_2w": "e0_s3cpu2w",
    }
    loaded = {n: load_arm(d, seeds) for n, d in arms_dirs.items() if (CACHE / d).exists()}
    if loaded.get("off") or loaded.get("on_cpu"):
        ab: Dict[str, Any] = dict(
            arms={n: arm_summary(r) for n, r in loaded.items() if r and not n.endswith("_2w")},
            paired={},
        )
        if loaded.get("off") and loaded.get("on_cpu"):
            ab["paired"]["off vs on_cpu"] = paired(loaded["off"], loaded["on_cpu"])
        if loaded.get("on_cpu") and loaded.get("on_mps"):
            ab["paired"]["on_cpu vs on_mps"] = paired(loaded["on_cpu"], loaded["on_mps"])
            n_mps = len(loaded["on_mps"])
            same = {k: v for k, v in loaded["on_cpu"].items() if k in loaded["on_mps"]}
            ab["mps_same_runs_cpu"] = arm_summary(same)
        if loaded.get("on_mps") and loaded.get("on_cpu_2w") and loaded.get("off_2w"):
            # seed 0 of the first 50 voxels, two workers each, run one after the other
            ab["timing_2w"] = {
                name: dict(
                    n=len(r),
                    wrong=sum(x["wrong"] for x in r.values()),
                    runtime_median=med([x["runtime"] for x in r.values()]),
                    runtime_mean=float(np.mean([x["runtime"] for x in r.values()])),
                    runtime_sum=float(sum(x["runtime"] for x in r.values())),
                )
                for name, r in (
                    ("off", loaded["off_2w"]),
                    ("on_cpu", loaded["on_cpu_2w"]),
                    ("on_mps", loaded["on_mps"]),
                )
            }
            ab["paired"]["off_2w vs on_mps"] = paired(loaded["off_2w"], loaded["on_mps"])
            ab["paired"]["on_cpu_2w vs on_mps"] = paired(loaded["on_cpu_2w"], loaded["on_mps"])
        if loaded.get("off"):
            # calibration: how many runs flip between two seeds of the SAME (key-off) code
            ab["off_seed_to_seed"] = {}
            for s0, s1 in ((0, 1), (0, 2), (1, 2)):
                ra = {(v, 0): r for (v, sd), r in loaded["off"].items() if sd == s0}
                rb = {(v, 0): r for (v, sd), r in loaded["off"].items() if sd == s1}
                pr = paired(ra, rb)
                pr.pop("flips", None)
                ab["off_seed_to_seed"][f"seed {s0} vs seed {s1}"] = pr
        chk = published_check(loaded.get("off", {}))
        if chk:
            ab["published_check"] = chk
        summary["ab"] = ab
        tables.update(ab_tables(ab))
    (out / "summary.json").write_text(json.dumps(summary, indent=1, default=float))
    DT.write_tables(out / "tables.md", tables)
    print(f"wrote {out / 'summary.json'} and {out / 'tables.md'}")


if __name__ == "__main__":
    main()
