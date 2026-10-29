"""S2 Step 0 summary: summary.json and tables.md (doc_tables) from the Step 0 results.

Inputs (``--dir``, default benchmarks/s2_step0): parity.json (0d), timing.json (0c),
search_three.json and search_e0.json (0a), fmaf_probe.json (s2_fma_probe.py). Every number in
the tables comes from these files.

  uv run python scripts/speed_parity/s2/s2_summary.py [--dir benchmarks/s2_step0]
"""

import argparse
import json
import statistics
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional

import numpy as np

HERE = Path(__file__).resolve().parent
ICE = HERE.parents[2]
sys.path.insert(0, str(ICE / "scripts" / "common"))

import doc_tables as DT  # noqa: E402
import stats  # noqa: E402

TRIANGLE_MODES = ("local", "cpp_rule", "edge", "edge_cpp_rule", "min_sin_eta")
SCREEN_MODES = ("screen_default", "screen_center", "screen_vertex0", "screen_edge")
CALLERS = ("screen", "recost", "quick_mc", "find_optimal", "variance", "final")
ARM_LABEL = {
    "default": "default keys",
    "bqmc": "BatchedQuickMC 1",
    "fused": "default keys, fused",
    "fused_bqmc": "BatchedQuickMC 1, fused",
}


def load(path: Path) -> Optional[Dict[str, Any]]:
    return json.loads(path.read_text()) if path.exists() else None


def pct(x: Optional[float]) -> Optional[float]:
    return None if x is None else 100.0 * x


def peak_stats(pk: List[Dict[str, Any]]) -> Dict[str, Any]:
    """Per flipped pixel peak (the coordinate with the largest torch margin of each differing
    (row, detector)): which float32 value is nearer the float64 value (tie: equal distance), the
    side of the threshold the float64 value is on, the ratio of the larger to the smaller error,
    and the two largest-margin peaks (by ulp and by px; they can be different peaks)."""
    out: Dict[str, Any] = dict(n_pixel_peaks=len(pk))
    if not pk:
        return out
    et = np.array([abs(x["torch"] - x["ref"]) for x in pk])
    ef = np.array([abs(x["fused"] - x["ref"]) for x in pk])
    n = len(pk)
    side_t = sum(x["ref_side"] == "torch" for x in pk)
    lo, hi = stats.wilson(side_t, n)
    nt = int((et < ef).sum())
    nlo, nhi = stats.wilson(nt, n)
    ratio = np.maximum(et, ef) / np.maximum(np.minimum(et, ef), 1e-300)
    out.update(
        nearer_torch=nt,
        nearer_fused=int((ef < et).sum()),
        nearer_tie=int((et == ef).sum()),
        nearer_torch_wilson_lo_pct=100.0 * nlo,
        nearer_torch_wilson_hi_pct=100.0 * nhi,
        side_torch=side_t,
        side_fused=n - side_t,
        side_torch_wilson_lo_pct=100.0 * lo,
        side_torch_wilson_hi_pct=100.0 * hi,
        torch_err_median_px=float(np.median(et)),
        fused_err_median_px=float(np.median(ef)),
        error_ratio_median=float(np.median(ratio)),
        error_ratio_gt3_frac=float((ratio > 3).mean()),
        fused_max_px_from_edge=max(x["fused_px_from_edge"] for x in pk),
        fused_max_margin_ulp=max(x["fused_margin_ulp"] for x in pk),
    )
    keys = ("dataset", "index", "row", "det", "threshold", "margin_ulp", "torch_px_from_edge",
            "fused_margin_ulp", "fused_px_from_edge")  # fmt: skip
    by_ulp = max(pk, key=lambda x: x["margin_ulp"])
    by_px = max(pk, key=lambda x: x["torch_px_from_edge"])
    by_fused_px = max(pk, key=lambda x: x["fused_px_from_edge"])
    out["max_ulp_peak"] = {k: by_ulp[k] for k in keys}
    out["max_px_peak"] = {k: by_px[k] for k in keys}
    out["max_fused_px_peak"] = {k: by_fused_px[k] for k in keys}
    out["max_ulp_and_max_px_same_peak"] = by_ulp is by_px
    return out


def paired(s: Dict[str, Any]) -> Dict[str, Any]:
    """Per arm pair, by case (voxel), on repeat 0 (the wrong sets of all repeats are compared):
    fixed = wrong under the first arm only, broken = wrong under the second only, both wrong."""
    reps = sorted({r["repeat"] for r in s["runs"]})

    def wrong(arm: str, rep: int) -> Dict[str, bool]:
        return {
            r["case"]: (r["err_deg"] or 0.0) > 1.0
            for r in s["runs"]
            if r["arm"] == arm and r["repeat"] == rep
        }

    def err(arm: str) -> Dict[str, float]:
        return {r["case"]: r["err_deg"] for r in s["runs"] if r["arm"] == arm and r["repeat"] == 0}

    def evals(arm: str) -> Dict[str, int]:
        return {
            r["case"]: r["evaluate_n"] for r in s["runs"] if r["arm"] == arm and r["repeat"] == 0
        }

    out: Dict[str, Any] = {}
    for a, b in (("default", "fused"), ("bqmc", "fused_bqmc"), ("default", "bqmc")):
        if a not in s["arms"] or b not in s["arms"]:
            continue
        wa, wb = wrong(a, 0), wrong(b, 0)
        cases = sorted(wa)
        fixed, broken = stats.paired_discordant([wa[c] for c in cases], [wb[c] for c in cases])
        out[f"{a}_vs_{b}"] = dict(
            n=len(cases),
            wrong_a=sum(wa.values()),
            wrong_b=sum(wb.values()),
            fixed=fixed,
            broken=broken,
            both_wrong=sum(wa[c] and wb[c] for c in cases),
            wrong_cases_a=[c for c in cases if wa[c]],
            wrong_cases_b=[c for c in cases if wb[c]],
            path_changed=[c for c in cases if evals(a)[c] != evals(b)[c]],
            err_deg_a={c: err(a)[c] for c in cases if wa[c] or wb[c]},
            err_deg_b={c: err(b)[c] for c in cases if wa[c] or wb[c]},
            same_wrong_sets_all_repeats=all(
                wrong(x, r) == wrong(x, 0) for x in (a, b) for r in reps
            ),
        )
    return out


def probe(f: Dict[str, Any]) -> Dict[str, Any]:
    rows = f["rows"]
    out: Dict[str, Any] = dict(meta=f["meta"], rows=rows)
    out["table"] = DT.markdown_table(
        rows,
        ["product", "order", "n", "equal", "equal_frac", "max_ulp"],
        formats=dict(equal_frac=".3f"),
    )
    return out


def per_call(t: Dict[str, Any]) -> Dict[str, Any]:
    rows = []
    out: Dict[str, Any] = dict(label=t["label"], n=t["n_orientations"], rounds=t["rounds"])
    for arm, d in t["us_per_call"].items():
        m = {i: d[i]["median"] for i in d}
        rows.append(
            dict(
                arm=arm,
                torch=m["torch"],
                fused=m["fused"],
                fused_acc=m["fused_acc"],
                c_only=m["c_only"],
                speedup=m["torch"] / m["fused"],
                fused_range=f"{d['fused']['min']:.1f}-{d['fused']['max']:.1f}",
            )
        )
        out[arm] = dict(m, speedup=m["torch"] / m["fused"], fused_min=d["fused"]["min"],
                        fused_max=d["fused"]["max"])  # fmt: skip
    out["table"] = DT.markdown_table(
        rows,
        ["arm", "torch", "fused", "fused_range", "fused_acc", "c_only", "speedup"],
        formats=dict(torch=".1f", fused=".1f", fused_acc=".1f", c_only=".1f", speedup=".1f"),
    )
    return out


def search(s: Dict[str, Any]) -> Dict[str, Any]:
    """Per arm: totals over the set's cases (E0: sum over voxels) per repeat, median over
    repeats."""
    arms = s["arms"]
    reps = sorted({r["repeat"] for r in s["runs"]})
    out: Dict[str, Any] = dict(label=s["label"], set=s["set"], keys=s["keys"], arms={})
    for arm in arms:
        per_rep: List[Dict[str, float]] = []
        for rep in reps:
            rr = [r for r in s["runs"] if r["arm"] == arm and r["repeat"] == rep]
            tot: Dict[str, float] = dict(
                wall_s=sum(r["wall_s"] for r in rr),
                evaluate_s=sum(r["evaluate_s"] for r in rr),
                evaluate_n=sum(r["evaluate_n"] for r in rr),
                screen_s=sum(r["screen_s"] for r in rr),
                bqmc_s=sum(r["batched_quick_mc"]["s"] for r in rr),
                bqmc_n=sum(r["batched_quick_mc"]["orientations"] for r in rr),
                wrong=sum((r["err_deg"] or 0.0) > 1.0 for r in rr),
                n_cases=len(rr),
            )
            for c in CALLERS:
                tot[f"{c}_s"] = sum(r["callers"][c]["eval_s"] for r in rr)
                tot[f"{c}_n"] = sum(r["callers"][c]["eval_n"] for r in rr)
            for k in ("quick_mc_s", "find_optimal_s", "variance_s", "discrete_s"):
                tot[f"section_{k}"] = sum(r["sections"][k] for r in rr)
            per_rep.append(tot)
        med = {k: statistics.median(p[k] for p in per_rep) for k in per_rep[0]}
        med["us_per_evaluate"] = 1e6 * med["evaluate_s"] / max(med["evaluate_n"], 1)
        med["evaluate_share"] = med["evaluate_s"] / med["wall_s"]
        med["repeats"] = len(per_rep)
        med["wall_s_repeats"] = [p["wall_s"] for p in per_rep]
        med["wrong_repeat0"] = per_rep[0]["wrong"]
        out["arms"][arm] = med
    base = out["arms"].get("default")
    rows = []
    for arm, m in out["arms"].items():
        rows.append(
            dict(
                arm=ARM_LABEL.get(arm, arm),
                wall=m["wall_s"],
                vs_default=(base["wall_s"] / m["wall_s"]) if base else None,
                evaluate=m["evaluate_s"],
                evaluate_n=int(m["evaluate_n"]),
                us_per_eval=m["us_per_evaluate"],
                screen=m["screen_s"],
                bqmc=m["bqmc_s"],
                wrong=f"{int(m['wrong_repeat0'])} of {int(m['n_cases'])}",
            )
        )
    out["table"] = DT.markdown_table(
        rows,
        ["arm", "wall", "vs_default", "evaluate", "evaluate_n", "us_per_eval", "screen", "bqmc",
         "wrong"],  # fmt: skip
        formats=dict(wall=".2f", vs_default=".2f", evaluate=".2f", us_per_eval=".0f",
                     screen=".2f", bqmc=".2f"),  # fmt: skip
    )
    crow = []
    for arm in ("default", "bqmc"):
        m = out["arms"].get(arm)
        if m is None:
            continue
        for c in CALLERS:
            crow.append(
                dict(
                    arm=ARM_LABEL[arm],
                    caller="screen (BatchedScreen, not in v1)" if c == "screen" else c,
                    calls=int(m[f"{c}_n"]),
                    seconds=m[f"{c}_s"],
                    share=m[f"{c}_s"] / m["wall_s"],
                    us_per_call=(1e6 * m[f"{c}_s"] / m[f"{c}_n"]) if m[f"{c}_n"] else None,
                )
            )
        crow.append(
            dict(arm=ARM_LABEL[arm], caller="batched quick-MC cost", calls=int(m["bqmc_n"]),
                 seconds=m["bqmc_s"], share=m["bqmc_s"] / m["wall_s"],
                 us_per_call=(1e6 * m["bqmc_s"] / m["bqmc_n"]) if m["bqmc_n"] else None)
        )  # fmt: skip
        crow.append(
            dict(arm=ARM_LABEL[arm], caller="covered evaluate (all but screen)",
                 calls=int(m["evaluate_n"] - m["screen_n"]), seconds=m["evaluate_s"],
                 share=m["evaluate_share"], us_per_call=m["us_per_evaluate"])
        )  # fmt: skip
    out["callers_table"] = DT.markdown_table(
        crow,
        ["arm", "caller", "calls", "seconds", "share", "us_per_call"],
        formats=dict(seconds=".2f", share=".3f", us_per_call=".0f"),
    )
    return out


def parity(p: Dict[str, Any]) -> Dict[str, Any]:
    meta = p["meta"]
    tu = meta["trace_ulp"]
    out: Dict[str, Any] = dict(trace_ulp=tu, variants=meta["variants"], datasets={})
    frows = []
    irows = []
    erows = []
    any_frame = 0
    peaks: Dict[Any, List[Dict[str, Any]]] = {}
    for ds, modes in p["datasets"].items():
        out["datasets"][ds] = {}
        for group, names in (("triangle", TRIANGLE_MODES), ("screen", SCREEN_MODES)):
            for v in meta["variants"]:
                for mode in names:
                    if mode not in modes:
                        continue
                    x = modes[mode][v]
                    rows = [r for d in x["differing"] for r in d["rows"]]
                    peaks[(ds, mode, v)] = [
                        dict(r, dataset=ds, index=d["index"])
                        for d in x["differing"]
                        for r in d["rows"]
                        if r["cause"] == "pixel_integer"
                    ]
                    causes = sorted({r["cause"] for r in rows})
                    px = [r["torch_px_from_edge"] for r in rows if "torch_px_from_edge" in r]
                    within = sum(
                        all(r["margin_ulp"] <= tu for r in d["rows"]) and bool(d["rows"])
                        for d in x["differing"]
                    )
                    lo, hi = stats.wilson(x["n_orientations_differing"], x["n"])
                    any_frame += x["frame_rows_differing"]
                    rec = dict(
                        dataset=ds,
                        mode=mode,
                        variant=v,
                        n=x["n"],
                        differing=x["n_orientations_differing"],
                        wilson_lo=lo,
                        wilson_hi=hi,
                        within_trace_ulp=within,
                        beyond_trace_ulp=x["n_orientations_differing"] - within,
                        no_differing_peak=sum(not d["rows"] for d in x["differing"]),
                        max_margin_ulp=x["max_margin_ulp"],
                        max_margin_px=max(px) if px else None,
                        causes=",".join(causes) if causes else "-",
                        ref_side_torch=x["ref_side"]["torch"],
                        ref_side_fused=x["ref_side"]["fused"],
                        fields=",".join(
                            f"{k}:{n}" for k, n in x["field_orientations_differing"].items() if n
                        )
                        or "-",
                        silent=x["silent_input_differences"],
                        field_names=list(x["field_orientations_differing"]),
                    )
                    out["datasets"][ds][f"{mode}/{v}"] = rec
                    frows.append(rec)
            if "local" in modes and group == "triangle":
                for v in meta["variants"]:
                    x = modes["local"][v]
                    for name, s in x["intermediates"].items():
                        if ds == "three0":
                            irows.append(
                                dict(variant=v, intermediate=name, n=s["n"], equal=s["equal"],
                                     equal_frac=s["equal"] / max(s["n"], 1),
                                     max_ulp=s["max_ulp"], max_abs=s["max_abs"])
                            )  # fmt: skip
                    e = x["error_vs_float64"]
                    erows.append(
                        dict(dataset=ds, variant=v,
                             torch_mean_px=e["torch"]["pixel"]["mean"],
                             fused_mean_px=e["fused"]["pixel"]["mean"],
                             torch_max_px=e["torch"]["pixel"]["max"],
                             fused_max_px=e["fused"]["pixel"]["max"])
                    )  # fmt: skip
    out["frame_rows_differing_total"] = any_frame
    # float64 pixel error of the CppEdgeRules mode (reference64 applies the -1 rule since fix 10)
    out["edge_float64_error"] = {
        f"{ds}/{v}": {
            f"{who}_{k}_px": modes["edge"][v]["error_vs_float64"][who]["pixel"][k]
            for who in ("torch", "fused")
            for k in ("mean", "max")
        }
        for ds, modes in p["datasets"].items()
        if "edge" in modes
        for v in meta["variants"]
    }
    # headline: the plan's triangle arm (mode "local") and the screen ("screen_center"), totals
    # over the datasets, per variant
    head: Dict[str, Any] = {}
    for mode in ("local", "screen_center"):
        for v in meta["variants"]:
            recs = [r for r in frows if r["mode"] == mode and r["variant"] == v]
            n = sum(r["n"] for r in recs)
            k = sum(r["differing"] for r in recs)
            lo, hi = stats.wilson(k, n)
            mm = [r["max_margin_ulp"] for r in recs if r["max_margin_ulp"] is not None]
            mp = [r["max_margin_px"] for r in recs if r["max_margin_px"] is not None]
            head[f"{mode}/{v}"] = dict(
                n=n,
                differing=k,
                differing_pct=100.0 * k / n,
                wilson_lo_pct=100.0 * lo,
                wilson_hi_pct=100.0 * hi,
                within_trace_ulp=sum(r["within_trace_ulp"] for r in recs),
                beyond_trace_ulp=sum(r["beyond_trace_ulp"] for r in recs),
                no_differing_peak=sum(r["no_differing_peak"] for r in recs),
                max_margin_ulp=max(mm) if mm else None,
                max_margin_px=max(mp) if mp else None,
                flipped_peaks=sum(r["ref_side_torch"] + r["ref_side_fused"] for r in recs),
                ref_side_torch=sum(r["ref_side_torch"] for r in recs),
                ref_side_fused=sum(r["ref_side_fused"] for r in recs),
            )
            h = head[f"{mode}/{v}"]
            lo, hi = stats.wilson(h["beyond_trace_ulp"], h["differing"]) if k else (None, None)
            h.update(beyond_wilson_lo_pct=pct(lo), beyond_wilson_hi_pct=pct(hi))
            h.update(peak_stats([pk for ds in p["datasets"] for pk in peaks[(ds, mode, v)]]))
            h["fields"] = (
                {
                    f: sum(
                        p["datasets"][ds][mode][v]["field_orientations_differing"][f]
                        for ds in p["datasets"]
                        if mode in p["datasets"][ds]
                    )
                    for f in recs[0]["field_names"]
                }
                if recs
                else {}
            )
    out["headline"] = head
    # one table row per (dataset, variant, mode group): the triangle modes and screen modes
    # differ almost identically within a group, so the group shows its first mode and the range
    trows = []
    for ds in p["datasets"]:
        for v in meta["variants"]:
            for group, names in (("triangle (5 modes)", TRIANGLE_MODES),
                                 ("screen (4 modes)", SCREEN_MODES)):  # fmt: skip
                recs = [r for r in frows if r["dataset"] == ds and r["variant"] == v
                        and r["mode"] in names]  # fmt: skip
                if not recs:
                    continue
                dif = [r["differing"] for r in recs]
                beyond = [r["beyond_trace_ulp"] for r in recs]
                mm = [r["max_margin_ulp"] for r in recs if r["max_margin_ulp"] is not None]
                mp = [r["max_margin_px"] for r in recs if r["max_margin_px"] is not None]
                trows.append(
                    dict(
                        dataset=ds,
                        variant=v,
                        modes=group,
                        differing=f"{min(dif)}-{max(dif)}" if min(dif) != max(dif) else min(dif),
                        beyond_8ulp=(
                            f"{min(beyond)}-{max(beyond)}"
                            if min(beyond) != max(beyond)
                            else min(beyond)
                        ),
                        no_peak=sum(r["no_differing_peak"] for r in recs),
                        causes=",".join(
                            sorted({c for r in recs for c in r["causes"].split(",")} - {"-"})
                        )
                        or "-",  # fmt: skip
                        max_margin_ulp=max(mm) if mm else None,
                        max_margin_px=max(mp) if mp else None,
                    )
                )
    out["fields_table"] = DT.markdown_table(
        trows,
        [
            "dataset",
            "variant",
            "modes",
            "differing",
            "beyond_8ulp",
            "no_peak",
            "causes",
            "max_margin_ulp",
            "max_margin_px",
        ],  # fmt: skip
        formats=dict(max_margin_px=".4f"),
    )
    out["fields_rows"] = frows
    out["intermediates_table"] = DT.markdown_table(
        irows,
        ["variant", "intermediate", "n", "equal", "equal_frac", "max_ulp", "max_abs"],
        formats=dict(equal_frac=".3f", max_abs=".2e"),
    )
    out["error_table"] = DT.markdown_table(
        erows,
        ["dataset", "variant", "torch_mean_px", "fused_mean_px", "torch_max_px", "fused_max_px"],
        formats=dict(torch_mean_px=".2e", fused_mean_px=".2e", torch_max_px=".4f",
                     fused_max_px=".4f"),  # fmt: skip
    )
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--dir", default=str(ICE / "benchmarks" / "s2_step0"))
    a = ap.parse_args()
    d = Path(a.dir)
    summary: Dict[str, Any] = {}
    tables: Dict[str, str] = {}
    t = load(d / "timing.json")
    if t:
        summary["per_call"] = pc = per_call(t)
        tables["s2_step0_per_call"] = pc.pop("table")
    for name in ("three", "e0"):
        s = load(d / f"search_{name}.json")
        if s:
            summary[f"search_{name}"] = ss = search(s)
            ss["paired"] = paired(s)
            tables[f"s2_step0_search_{name}"] = ss.pop("table")
            tables[f"s2_step0_callers_{name}"] = ss.pop("callers_table")
    p = load(d / "parity.json")
    if p:
        summary["parity"] = pp = parity(p)
        tables["s2_step0_parity"] = pp.pop("fields_table")
        tables["s2_step0_intermediates"] = pp.pop("intermediates_table")
        tables["s2_step0_float64_error"] = pp.pop("error_table")
    f = load(d / "fmaf_probe.json")
    if f:
        summary["fmaf_probe"] = fp = probe(f)
        tables["s2_step0_fmaf_probe"] = fp.pop("table")
    (d / "summary.json").write_text(json.dumps(summary, indent=1, default=float) + "\n")
    DT.write_tables(d / "tables.md", tables)
    print(f"wrote summary.json and tables.md ({len(tables)} tables)")


if __name__ == "__main__":
    main()
