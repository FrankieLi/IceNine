"""Task A step 2: the Python evaluation ledger next to the C++ per-stage counts, stage by stage.

  uv run python scripts/eval_reduction/ledger_attribution.py \
      --ledger benchmarks/eval_reduction/A/ledger --cpp benchmarks/eval_reduction/A/step2 \
      --arms cpp cpp_edge cpp_edge_interp --out benchmarks/eval_reduction/A/step2

C++ numbers: cpp_ledger.json (scratch build with counters, see cpp_ledger_patch.py.txt). The stages
and what drives their count:

  screen     9 evaluations per FZ point at every level                 (units: FZ points)
  recost     1 per (FZ point, local-grid point) that passed the screen (units: passed pairs)
  rerecost   C++ only: one more complete-cost call per returned candidate
             (DiscreteAdaptive.tmpl.cpp, the loop after RunDiscreteSearch); the value equals the
             recost, Python does not repeat it
  quick MC   12 per candidate                                          (units: candidates)
  FindOptimal, VarianceMinimizing, final: random (their own distribution), not tied to the units
"""

import argparse
import json
import sys
from pathlib import Path
from typing import Any, Dict, List

import numpy as np

HERE = Path(__file__).resolve().parent
ICE = HERE.parents[1]
sys.path.insert(0, str(ICE / "scripts" / "common"))
import doc_tables as DT  # noqa: E402

SIX = [15901, 17034, 18558, 19889, 22867, 5242]
STAGES = ["screen", "recost", "rerecost", "quick", "find", "variance", "final"]


def py_units(r: Dict[str, Any]) -> Dict[str, Any]:
    lv = r["levels"]
    return dict(
        fz=sum(l["fz_in"] for l in lv),
        passed=sum(l["passed_screen"] for l in lv),
        cands=sum(l["candidates"] for l in lv),
        cands_l0=lv[0]["candidates"],
        kept=[l["kept"] for l in lv],
    )


def cpp_units(c: Dict[str, Any]) -> Dict[str, Any]:
    lv = c["levels"]
    return dict(
        fz=sum(l["n_fz"] for l in lv),
        passed=sum(l["recost"] for l in lv),
        cands=sum(l["cand"] for l in lv),
        cands_l0=lv[0]["cand"],
        kept=[l["kept"] for l in lv],
    )


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--ledger", required=True)
    ap.add_argument("--cpp", required=True)
    ap.add_argument("--arms", nargs="+", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    cpp = json.loads((Path(a.cpp) / "cpp_ledger.json").read_text())
    out = Path(a.out)
    runs: Dict[str, Dict[int, Dict[str, Any]]] = {}
    for arm in a.arms:
        runs[arm] = {}
        for v in SIX:
            p = Path(a.ledger) / "runs" / f"{arm}_v{v}.json"
            if p.exists():
                runs[arm][v] = json.loads(p.read_text())
    S: Dict[str, Any] = {}
    T: Dict[str, str] = {}

    # ---- per arm: totals by stage, mean over the voxels both sides have ----
    for arm, rr in runs.items():
        if not rr:
            continue
        vs = [v for v in SIX if v in rr]
        rows = []
        mean_py = {s: 0.0 for s in STAGES}
        mean_cpp = {s: 0.0 for s in STAGES}
        for s in STAGES:
            py_v = np.array([rr[v]["totals"].get(s, 0) for v in vs], float)
            cp_v = np.array([cpp[str(v)]["totals"][s] for v in vs], float)
            mean_py[s], mean_cpp[s] = float(py_v.mean()), float(cp_v.mean())
        tot_py = np.array([rr[v]["evals_total"] for v in vs], float)
        tot_cpp = np.array([cpp[str(v)]["total_counter"] for v in vs], float)
        gap = float(tot_py.mean() - tot_cpp.mean())
        for s in STAGES:
            d = mean_py[s] - mean_cpp[s]
            rows.append(
                dict(
                    stage=s,
                    python=mean_py[s],
                    cpp=mean_cpp[s],
                    diff=d,
                    share_of_gap=100.0 * d / gap if gap else float("nan"),
                )
            )
        rows.append(
            dict(
                stage="total",
                python=float(tot_py.mean()),
                cpp=float(tot_cpp.mean()),
                diff=gap,
                share_of_gap=100.0,
            )
        )
        S[f"stages_{arm}"] = dict(
            voxels=vs,
            rows=rows,
            ratio_total=float(tot_py.mean() / tot_cpp.mean()),
            # net difference of every stage but the C++-only second re-cost
            diff_net_other_than_rerecost=float(
                sum(r["diff"] for r in rows if r["stage"] not in ("total", "rerecost"))
            ),
            ratio_by_voxel={
                v: float(rr[v]["evals_total"] / cpp[str(v)]["total_counter"]) for v in vs
            },
        )
        T[f"stages_{arm}"] = DT.markdown_table(
            rows,
            ["stage", "python", "cpp", "diff", "share_of_gap"],
            formats={"python": ",.0f", "cpp": ",.0f", "diff": ",.0f", "share_of_gap": ".0f"},
        )

        # units and per-unit costs
        urows = []
        for v in vs:
            pu, cu = py_units(rr[v]), cpp_units(cpp[str(v)])
            pt, ct = rr[v]["totals"], cpp[str(v)]["totals"]
            urows.append(
                dict(
                    voxel=v,
                    fz_py=pu["fz"],
                    fz_cpp=cu["fz"],
                    passed_py=pu["passed"],
                    passed_cpp=cu["passed"],
                    cands_py=pu["cands"],
                    cands_cpp=cu["cands"],
                    l0_py=pu["cands_l0"],
                    l0_cpp=cu["cands_l0"],
                    screen_per_fz_py=pt["screen"] / pu["fz"],
                    screen_per_fz_cpp=ct["screen"] / cu["fz"],
                    quick_per_cand_py=pt["quick"] / pu["cands"],
                    quick_per_cand_cpp=ct["quick"] / cu["cands"],
                    recost_per_pass_py=pt["recost"] / max(pu["passed"], 1),
                    recost_per_pass_cpp=ct["recost"] / max(cu["passed"], 1),
                    ratio_total=rr[v]["evals_total"] / cpp[str(v)]["total_counter"],
                )
            )
        S[f"units_{arm}"] = urows
        T[f"units_{arm}"] = DT.markdown_table(
            urows,
            [
                "voxel",
                "fz_py",
                "fz_cpp",
                "passed_py",
                "passed_cpp",
                "cands_py",
                "cands_cpp",
                "l0_py",
                "l0_cpp",
                "screen_per_fz_py",
                "screen_per_fz_cpp",
                "quick_per_cand_py",
                "quick_per_cand_cpp",
                "recost_per_pass_py",
                "recost_per_pass_cpp",
                "ratio_total",
            ],
            formats={
                "voxel": "d",
                **{
                    k: ",d"
                    for k in (
                        "fz_py",
                        "fz_cpp",
                        "passed_py",
                        "passed_cpp",
                        "cands_py",
                        "cands_cpp",
                        "l0_py",
                        "l0_cpp",
                    )
                },
                "screen_per_fz_py": ".2f",
                "screen_per_fz_cpp": ".2f",
                "quick_per_cand_py": ".2f",
                "quick_per_cand_cpp": ".2f",
                "recost_per_pass_py": ".2f",
                "recost_per_pass_cpp": ".2f",
                "ratio_total": ".3f",
            },
        )

        # FindOptimal / variance: distributions, not exact counts
        frows = []
        for v in vs:
            r = rr[v]
            fo = r["find_optimal"]
            fc = cpp[str(v)]["find"]
            frows.append(
                dict(
                    voxel=v,
                    find_py=r["totals"]["find"],
                    find_cpp=fc["find"],
                    cand_evaluated_py=fo["n_evaluated"],
                    conv_py=int(fo["converged"]),
                    conv_cpp=fc["converged"],
                    var_py=r["totals"]["variance"],
                    var_cpp=fc["var"],
                )
            )
        S[f"random_stages_{arm}"] = frows
        # per-voxel Python / C++ ratio of the random stages: FindOptimal only where both
        # converged (otherwise Python evaluates all 30 candidates), VarianceMinimizing always
        find_ratios = [
            f["find_py"] / f["find_cpp"] for f in frows if f["conv_py"] and f["conv_cpp"]
        ]
        var_ratios = [f["var_py"] / f["var_cpp"] for f in frows]
        both = find_ratios + var_ratios
        S[f"random_ratio_range_{arm}"] = dict(
            find_converged_min=min(find_ratios) if find_ratios else None,
            find_converged_max=max(find_ratios) if find_ratios else None,
            n_find_converged=len(find_ratios),
            variance_min=min(var_ratios),
            variance_max=max(var_ratios),
            min=min(both),
            max=max(both),
        )
        T[f"random_{arm}"] = DT.markdown_table(
            frows,
            [
                "voxel",
                "find_py",
                "find_cpp",
                "cand_evaluated_py",
                "conv_py",
                "conv_cpp",
                "var_py",
                "var_cpp",
            ],
            formats={"voxel": "d"},
        )

        # level table of the first voxel: python vs C++
        v0 = vs[0]
        lrows = []
        cl = cpp[str(v0)]["levels"]
        for i, l in enumerate(rr[v0]["levels"]):
            lrows.append(
                dict(
                    level=i,
                    fz_py=l["fz_in"],
                    fz_cpp=cl[i]["n_fz"],
                    passed_py=l["passed_screen"],
                    passed_cpp=cl[i]["recost"],
                    cands_py=l["candidates"],
                    cands_cpp=cl[i]["cand"],
                    kept_py=l["kept"],
                    kept_cpp=cl[i]["kept"],
                )
            )
        T[f"cmp_levels_{arm}_v{v0}"] = DT.markdown_table(
            lrows, list(lrows[0].keys()), formats={k: ",d" for k in lrows[0]}
        )
    # level-0 candidates of the C++ rule without and with CppEdgeRules (Python, same seed)
    if runs.get("cpp") and runs.get("cpp_edge"):
        l0 = {
            str(v): dict(
                cpp=runs["cpp"][v]["levels"][0]["candidates"],
                cpp_edge=runs["cpp_edge"][v]["levels"][0]["candidates"],
            )
            for v in SIX
            if v in runs["cpp"] and v in runs["cpp_edge"]
        }
        absdiff = [abs(x["cpp_edge"] - x["cpp"]) for x in l0.values()]
        S["l0_candidates_cpp_vs_cpp_edge"] = dict(
            by_voxel=l0, absdiff_min=min(absdiff), absdiff_max=max(absdiff)
        )
    (out / "attribution.json").write_text(json.dumps(S, indent=1, default=float) + "\n")
    DT.write_tables(out / "tables.md", T)
    print(f"wrote attribution.json and tables.md for arms {list(runs)}")


if __name__ == "__main__":
    main()
