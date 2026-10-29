"""S2 Steps 4-6 summary: summary.json and tables.md (doc_tables) from benchmarks/s2_results.

Inputs (``--dir``, default benchmarks/s2_results): parity.json (``s2_lib_parity.py``),
keyoff_compare.txt and keyoff_batched.txt (``s1_parity.py compare`` / ``batched`` against the base
832db36), timing_single.json and timing_batched.json (``s2_lib_timing.py``), search_three.json
and search_e0.json (``s2_search.py --lib``), e0.json (``s2_e0.py``: the E0 A/B records) and
e0_trace.json (``s2_e0_trace.py``); the Step 0 results (benchmarks/s2_step0) for the prototype's
columns. The Task A tables are written by ``taska_summary.py`` into ``taska/tables.md``. Every
number in the tables comes from these files.

  uv run python scripts/speed_parity/s2/s2_results_summary.py [--dir benchmarks/s2_results]
"""

import argparse
import json
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional

HERE = Path(__file__).resolve().parent
ICE = HERE.parents[2]
sys.path.insert(0, str(ICE / "scripts" / "common"))
sys.path.insert(0, str(HERE))

import doc_tables as DT  # noqa: E402
import s2_summary as S0  # noqa: E402
import stats  # noqa: E402

TRIANGLE_MODES = S0.TRIANGLE_MODES
SCREEN_MODES = S0.SCREEN_MODES
STEP0 = ICE / "benchmarks" / "s2_step0"
ARM_LABEL = {
    "default": "default keys",
    "bqmc": "BatchedQuickMC 1",
    "fused": "FusedCost 1",
    "fused_bqmc": "BatchedQuickMC 1, FusedCost 1",
}


def load(path: Path) -> Optional[Dict[str, Any]]:
    return json.loads(path.read_text()) if path.exists() else None


def rng(xs: List[int]) -> str:
    lo, hi = min(xs), max(xs)
    return str(lo) if lo == hi else f"{lo}-{hi}"


# ---------------------------------------------------------------------------------------------
# Step 4: parity
# ---------------------------------------------------------------------------------------------
def parity(p: Dict[str, Any]) -> Dict[str, Any]:
    out: Dict[str, Any] = dict(meta=p["meta"], totals=p["totals"])
    rows = []
    for ds, dd in p["datasets"].items():
        for grp, modes in (
            ("triangle (5 modes)", TRIANGLE_MODES),
            ("screen (4 modes)", SCREEN_MODES),
        ):
            rr = [dd[m] for m in modes if m in dd]
            causes = sorted({c for r in rr for c in r["causes"]})
            rows.append(
                dict(
                    dataset=ds,
                    modes=grp,
                    n=rr[0]["n"],
                    differing=rng([r["n_differing"] for r in rr]),
                    silent=rng([r["n_silent"] for r in rr]),
                    flipped_coords=rng([sum(r["flipped_coords"].values()) for r in rr]),
                    failing=sum(r["n_failing"] for r in rr),
                    causes=", ".join(causes) or "-",
                    max_near_px=max(r["max_near_px"] for r in rr),
                    max_margin_ulp=max(r["max_margin_ulp"] for r in rr),
                    several_coords=sum(r["rows_dets_with_several_coords"] for r in rr),
                    frame_differ=sum(r["frame_differ"] for r in rr),
                )
            )
    out["table"] = DT.markdown_table(
        rows,
        [
            "dataset",
            "modes",
            "n",
            "differing",
            "silent",
            "flipped_coords",
            "failing",
            "causes",
            "max_near_px",
            "max_margin_ulp",
            "several_coords",
            "frame_differ",
        ],  # fmt: skip
        formats=dict(max_near_px=".4f"),
    )
    # the triangle arm ``local`` (as Step 0d): fields, nearer and side of float64, over 4 voxels
    loc = [dd["local"] for dd in p["datasets"].values() if "local" in dd]
    n = sum(r["n"] for r in loc)
    dif = sum(r["n_differing"] for r in loc)
    fields = {
        f: sum(r["field_orientations_differing"][f] for r in loc)
        for f in loc[0]["field_orientations_differing"]
    }
    nearer = {k: sum(r["nearer_f64"][k] for r in loc) for k in ("torch", "fused", "tie")}
    side = {k: sum(r["f64_side"][k] for r in loc) for k in ("torch", "fused")}
    coords = sum(sum(r["flipped_coords"].values()) for r in loc)
    lo, hi = stats.wilson(dif, n)
    nlo, nhi = stats.wilson(nearer["torch"], nearer["torch"] + nearer["fused"])
    slo, shi = stats.wilson(side["torch"], side["torch"] + side["fused"])
    out["local"] = dict(
        n=n, differing=dif, differing_wilson=[lo, hi], fields=fields, flipped_coords=coords,
        nearer=nearer, nearer_torch_wilson=[nlo, nhi], side=side, side_torch_wilson=[slo, shi],
    )  # fmt: skip
    step0 = load(STEP0 / "summary.json")
    frows = [
        dict(path="library (plain g_lab)", n=n, differing=dif, wilson_low=lo, wilson_high=hi,
             pixel_on_detector=fields["pixel_on_detector"], pixel_overlap=fields["pixel_overlap"],
             quality=fields["quality"], peak_overlap=fields["peak_overlap"],
             peak_on_detector=fields["peak_on_detector"], flipped_coords=coords,
             nearer_torch=nearer["torch"], nearer_fused=nearer["fused"],
             side_torch=side["torch"], side_fused=side["fused"])
    ]  # fmt: skip
    out["fields_table"] = DT.markdown_table(
        frows,
        [
            "path",
            "n",
            "differing",
            "wilson_low",
            "wilson_high",
            "pixel_on_detector",
            "pixel_overlap",
            "quality",
            "peak_overlap",
            "peak_on_detector",
            "flipped_coords",
            "nearer_torch",
            "nearer_fused",
            "side_torch",
            "side_fused",
        ],  # fmt: skip
        formats=dict(wilson_low=".3f", wilson_high=".3f"),
    )
    out["step0_present"] = step0 is not None
    # intermediates of three0 local (Step 0d's table, library row)
    irows = []
    for ds in p["datasets"]:
        st = p["datasets"][ds].get("local", {}).get("intermediates", {})
        for k, s in st.items():
            irows.append(
                dict(dataset=ds, intermediate=k, n=s["n"], equal=s["equal"],
                     equal_frac=s["equal"] / s["n"] if s["n"] else None,
                     max_ulp=s["max_ulp"], max_abs=s["max_abs"])
            )  # fmt: skip
    out["intermediates_table"] = DT.markdown_table(
        [r for r in irows if r["dataset"] == "three0"],
        ["dataset", "intermediate", "n", "equal", "equal_frac", "max_ulp", "max_abs"],
        formats=dict(equal_frac=".3f", max_abs=".2e"),
    )
    return out


def keyoff(d: Path) -> Dict[str, Any]:
    def count(path: Path) -> Dict[str, Any]:
        lines = path.read_text().splitlines()
        rows = [ln for ln in lines if ln.strip().endswith("identical") or "differ" in ln.lower()]
        rows = [ln for ln in rows if not ln.startswith("RESULT") and not ln.startswith("dataset")]
        res = [ln for ln in lines if ln.startswith("RESULT")]
        return dict(
            rows=len(rows),
            identical=sum(ln.strip().endswith("identical") for ln in rows),
            result=res[0] if res else None,
        )

    return dict(evaluate=count(d / "keyoff_compare.txt"), batched=count(d / "keyoff_batched.txt"))


# ---------------------------------------------------------------------------------------------
# Step 5: timings
# ---------------------------------------------------------------------------------------------
def per_call(t: Dict[str, Any], t0: Optional[Dict[str, Any]]) -> Dict[str, Any]:
    rows = []
    for arm, d in t["us_per_call"].items():
        tt, ff, cc = d["torch"]["median"], d["fused"]["median"], d["c_only"]["median"]
        p0 = t0["us_per_call"][arm]["fused"]["median"] if t0 else None
        rows.append(
            dict(arm=arm, torch=tt, fused=ff, fused_range=f"{d['fused']['min']:.1f}-{d['fused']['max']:.1f}",
                 c_only=cc, python_share=ff - cc, speedup=tt / ff, step0_prototype=p0)
        )  # fmt: skip
    return dict(
        label=t["meta"]["label"],
        rows=rows,
        table=DT.markdown_table(
            rows,
            ["arm", "torch", "fused", "fused_range", "c_only", "python_share", "speedup",
             "step0_prototype"],  # fmt: skip
            formats=dict(torch=".1f", fused=".1f", c_only=".1f", python_share=".1f",
                         speedup=".1f", step0_prototype=".1f"),  # fmt: skip
        ),
    )


def batched(t: Dict[str, Any]) -> Dict[str, Any]:
    gpu = t.get("gpu")
    rows = []
    out: Dict[str, Any] = dict(label=t["meta"]["label"], gpu=gpu, gpu_min_batch=t["gpu_min_batch"])
    for case in ("three0", "three0_search", "e0_0"):
        c = t.get(case)
        if c is None:
            continue
        u = c["us_per_orientation"]
        ev = u["evaluate"]
        for B in t["batches"]:
            d = u[f"B={B}"]
            fx = d.get(f"fixed_{gpu}", {}).get("median") if gpu else None
            rows.append(
                dict(set=case, K=c["K"], B=B, compact=d["compact"]["median"],
                     fused=d["fused"]["median"], fixed_gpu=fx,
                     gpu_over_fused=(fx / d["fused"]["median"]) if fx else None,
                     evaluate_torch=ev["torch"]["median"], evaluate_fused=ev["fused"]["median"])
            )  # fmt: skip
    out["rows"] = rows
    out["table"] = DT.markdown_table(
        rows,
        ["set", "K", "B", "compact", "fused", "fixed_gpu", "gpu_over_fused", "evaluate_torch",
         "evaluate_fused"],  # fmt: skip
        formats=dict(compact=".1f", fused=".1f", fixed_gpu=".1f", gpu_over_fused=".2f",
                     evaluate_torch=".1f", evaluate_fused=".1f"),  # fmt: skip
    )
    return out


def search(s: Dict[str, Any]) -> Dict[str, Any]:
    ss = S0.search(s)
    rows = []
    base = ss["arms"].get("default")
    for arm, m in ss["arms"].items():
        bq = m["bqmc_n"] > 0
        qmc_n = m["bqmc_n"] if bq else m["quick_mc_n"]
        qmc_s = m["bqmc_s"] if bq else m["quick_mc_s"]
        rows.append(
            dict(
                arm=ARM_LABEL.get(arm, arm), wall=m["wall_s"],
                vs_default=(base["wall_s"] / m["wall_s"]) if base else None,
                evaluate=m["evaluate_s"], evaluate_n=int(m["evaluate_n"]),
                us_per_eval=m["us_per_evaluate"], quick_mc=qmc_s, quick_mc_n=int(qmc_n),
                quick_mc_us=(1e6 * qmc_s / qmc_n) if qmc_n else None, screen=m["screen_s"],
                wrong=f"{int(m['wrong_repeat0'])} of {int(m['n_cases'])}",
            )
        )  # fmt: skip
    ss["table"] = DT.markdown_table(
        rows,
        ["arm", "wall", "vs_default", "evaluate", "evaluate_n", "us_per_eval", "quick_mc",
         "quick_mc_n", "quick_mc_us", "screen", "wrong"],  # fmt: skip
        formats=dict(wall=".2f", vs_default=".2f", evaluate=".2f", us_per_eval=".0f",
                     quick_mc=".2f", quick_mc_us=".0f", screen=".2f"),  # fmt: skip
    )
    ss.pop("callers_table", None)
    ss["paired"] = S0.paired(s)
    # per case and arm, repeat 0: the final error and cost (whether the fused arm took another path)
    ss["final"] = {
        f"{r['case']}/{r['arm']}": dict(
            err_deg=r["err_deg"], cost=r["cost"], evaluate_n=r["evaluate_n"], evals=r["evals"]
        )
        for r in s["runs"]
        if r["repeat"] == 0
    }
    ss["fused_counts"] = {
        r["arm"]: dict(
            evaluate_fused_n=r["evaluate_fused_n"],
            evaluate_n=r["evaluate_n"],
            bqmc_fused=r["batched_quick_mc"].get("orientations_fused"),
            bqmc_n=r["batched_quick_mc"]["orientations"],
        )
        for r in s["runs"]
        if r["repeat"] == 0 and r["case"] == s["runs"][0]["case"]
    }
    return ss


# ---------------------------------------------------------------------------------------------
# Step 6: E0
# ---------------------------------------------------------------------------------------------
def e0(e: Dict[str, Any], tr: Optional[Dict[str, Any]]) -> Dict[str, Any]:
    out: Dict[str, Any] = dict(meta=e["meta"])
    wrows, prows, erows = [], [], []
    for arm, a in e["arms"].items():
        lo, hi = stats.wilson(a["wrong"], a["n"])
        wrows.append(dict(arm=a["label"], wrong=a["wrong"], n=a["n"], rate=a["wrong"] / a["n"],
                          wilson_low=lo, wilson_high=hi))  # fmt: skip
        erows.append(dict(arm=a["label"], local_mean=a["evals_local_mean"],
                          global_mean=a["evals_global_mean"], runtime_median=a["runtime_median"],
                          workers=e["meta"]["workers"]))  # fmt: skip
    for p in e["pairs"]:
        prows.append(dict(pair=p["label"], n=p["n"], wrong_a=p["wrong_a"], wrong_b=p["wrong_b"],
                          fixed=p["fixed"], broken=p["broken"], both_wrong=p["both_wrong"],
                          both_right=p["both_right"], same_final=p["same_final"]))  # fmt: skip
    out["wrong_table"] = DT.markdown_table(
        wrows, ["arm", "wrong", "n", "rate", "wilson_low", "wilson_high"],
        formats=dict(rate=".3f", wilson_low=".3f", wilson_high=".3f"),
    )  # fmt: skip
    out["paired_table"] = DT.markdown_table(
        prows, ["pair", "n", "wrong_a", "wrong_b", "fixed", "broken", "both_wrong", "both_right",
                "same_final"],
    )  # fmt: skip
    out["evals_table"] = DT.markdown_table(
        erows, ["arm", "local_mean", "global_mean", "runtime_median", "workers"],
        formats=dict(local_mean=".1f", global_mean=".1f", runtime_median=".2f"),
    )  # fmt: skip
    if tr is not None:
        trows = []
        for f in tr["flips"]:
            fs = f["first"] or {}
            trows.append(
                dict(pair=f["pair"], voxel=f["vidx"], change=f["change"],
                     reproduced=f["reproduced"], shadowed=f["shadowed"], differing=f["differing"],
                     first_call=fs.get("call"), first_path=fs.get("path"),
                     radius=fs.get("pixel_radius"), fields=", ".join(fs.get("fields", [])) or "-",
                     causes=", ".join(sorted({x["cause"] for x in fs.get("entries", [])})) or "-",
                     within_rule=fs.get("within_rule"), near_px=fs.get("max_near_px"))
            )  # fmt: skip
        out["trace"] = tr["totals"]
        prs = []
        for pair in sorted({f["pair"] for f in tr["flips"]}):
            ff = [f for f in tr["flips"] if f["pair"] == pair]
            prs.append(
                dict(pair=pair, voxels=len(ff),
                     changes=", ".join(f"{c} {sum(f['change'] == c for f in ff)}"
                                       for c in sorted({f["change"] for f in ff})),
                     reproduced=sum(f["reproduced"] for f in ff),
                     shadowed=sum(f["shadowed"] for f in ff),
                     differing=sum(f["differing"] for f in ff),
                     within_rule=sum(f["differing_within_rule"] for f in ff),
                     quality_differs=sum(f["differing_fields"].get("quality", 0) for f in ff),
                     causes=", ".join(sorted({c for f in ff for c in f["differing_causes"]})),
                     max_near_px=max(f["differing_max_near_px"] for f in ff))
            )  # fmt: skip
        out["trace_pairs_table"] = DT.markdown_table(
            prs,
            [
                "pair",
                "voxels",
                "changes",
                "reproduced",
                "shadowed",
                "differing",
                "within_rule",
                "quality_differs",
                "causes",
                "max_near_px",
            ],  # fmt: skip
            formats=dict(max_near_px=".4f"),
        )
        out["trace_table"] = DT.markdown_table(
            trows,
            [
                "pair",
                "voxel",
                "change",
                "reproduced",
                "shadowed",
                "differing",
                "first_call",
                "first_path",
                "radius",
                "fields",
                "causes",
                "within_rule",
                "near_px",
            ],  # fmt: skip
            formats=dict(near_px=".4f"),
        )
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--dir", default=str(ICE / "benchmarks" / "s2_results"))
    a = ap.parse_args()
    d = Path(a.dir).resolve()
    summary: Dict[str, Any] = {}
    tables: Dict[str, str] = {}
    p = load(d / "parity.json")
    if p:
        pp = summary["parity"] = parity(p)
        tables["s2_lib_parity"] = pp.pop("table")
        tables["s2_lib_parity_local"] = pp.pop("fields_table")
        tables["s2_lib_intermediates"] = pp.pop("intermediates_table")
    if (d / "keyoff_compare.txt").exists():
        summary["keyoff"] = keyoff(d)
    t = load(d / "timing_single.json")
    if t:
        pc = summary["per_call"] = per_call(t, load(STEP0 / "timing.json"))
        tables["s2_lib_per_call"] = pc.pop("table")
    t = load(d / "timing_batched.json")
    if t:
        bb = summary["batched"] = batched(t)
        tables["s2_lib_batched"] = bb.pop("table")
    for name in ("three", "e0"):
        s = load(d / f"search_{name}.json")
        if s:
            ss = summary[f"search_{name}"] = search(s)
            tables[f"s2_lib_search_{name}"] = ss.pop("table")
    e = load(d / "e0.json")
    if e:
        ee = summary["e0"] = e0(e, load(d / "e0_trace.json"))
        for k in ("wrong", "paired", "evals", "trace_pairs", "trace"):
            if f"{k}_table" in ee:
                tables[f"s2_e0_{k}"] = ee.pop(f"{k}_table")
    (d / "summary.json").write_text(json.dumps(summary, indent=1, default=float) + "\n")
    DT.write_tables(d / "tables.md", tables)
    for k, tb in tables.items():
        print(f"## {k}\n{tb}\n")


if __name__ == "__main__":
    main()
