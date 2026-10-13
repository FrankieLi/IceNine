#!/usr/bin/env python3
"""Why Python with CppOverlapRule on loses E0 voxels that it gets right with the rule off: Python
(off / on) and C++ side by side, per level (inputs: scripts/phase_d/e0_cpp.py, e0_trace.py).

Per run and level: the closest candidate to the truth (symmetry-reduced angle) among all candidates
after the quick MC, its rank (position in the cost-sorted list, 0 = best), the number kept (n_keep),
and whether it survives; the candidates handed to FindOptimal (top 30 of the last level), the
closest of them and its cost; FindOptimal's result. Classes of a wrong run (final > 1 deg):
  never      no candidate within 3 deg of the truth at level 0
  pruned L   a candidate within 3 deg exists at level L but none among the kept ones (for the last
             level: none among the 30 handed to FindOptimal); L is the level whose ranking cut it
  trap       a candidate within 3 deg was handed to FindOptimal and the result is wrong
  other      anything else (e.g. no basin at level L+1 although one was kept at L)

  uv run python scripts/phase_d/e0_debug.py
"""

import json
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional

import numpy as np

HERE = Path(__file__).resolve().parent
ICE = HERE.parents[1]
sys.path.insert(0, str(ICE / "scripts" / "findoptimal_robustness"))
sys.path.insert(0, str(ICE / "scripts" / "common"))
import common as C  # noqa: E402
import doc_tables as DT  # noqa: E402

CACHE = HERE / "cache"
OUT = ICE / "benchmarks" / "phase_d_seed_diag" / "level0" / "e0"
BASIN = 3.0
WRONG = 1.0
N_FIND = 30


def err(R: np.ndarray, Rt: np.ndarray) -> np.ndarray:
    R = np.asarray(R, dtype=np.float64).reshape(-1, 3, 3)
    return C.err_deg(R, np.broadcast_to(Rt, R.shape).copy())


def levels_from_arrays(
    disc: Dict[int, np.ndarray],
    sorted_R: Dict[int, np.ndarray],
    sorted_cost: Dict[int, np.ndarray],
    n_keep: Dict[int, int], Rt: np.ndarray,
) -> List[Dict[str, Any]]:  # fmt: skip
    out = []
    for lv in sorted(sorted_R):
        e_disc = err(disc[lv], Rt)
        e_s = err(sorted_R[lv], Rt)
        j = int(e_s.argmin())
        nk = n_keep[lv] if lv < 3 else min(N_FIND, len(e_s))
        out.append(
            dict(
                level=lv, n_cand=len(e_s), disc_min=float(e_disc.min()),
                post_min=float(e_s.min()), rank=j, n_keep=int(nk), survives=bool(j < nk),
                kept_min=float(e_s[:nk].min()), cost_of_closest=float(sorted_cost[lv][j]),
                best_cost=float(sorted_cost[lv][0]),
            )
        )  # fmt: skip
    return out


def py_levels(z: Any, Rt: np.ndarray) -> List[Dict[str, Any]]:
    disc, sR, sC, nk = {}, {}, {}, {}
    for lv in range(4):
        if f"quick_mc_{lv}_R" in z.files:
            disc[lv] = z[f"discrete_{lv}_R"]
            sR[lv], sC[lv] = z[f"quick_mc_{lv}_R"], z[f"quick_mc_{lv}_cost"]
            nk[lv] = int(z[f"quick_mc_{lv}_n_keep"])
    return levels_from_arrays(disc, sR, sC, nk, Rt)


def cpp_levels(d: Dict[str, Any], Rt: np.ndarray) -> List[Dict[str, Any]]:
    disc, sR, sC, nk = {}, {}, {}, {}
    for lv_s, rows in d["CQ"].items():
        lv = int(lv_s)
        a = np.array(rows)
        sR[lv], sC[lv] = a[:, 1:].reshape(-1, 3, 3), a[:, 0]
        disc[lv] = np.array(d["CD"][lv_s])[:, 1:].reshape(-1, 3, 3)
        nk[lv] = d["CK"][lv_s]
    return levels_from_arrays(disc, sR, sC, nk, Rt)


def classify(levels: List[Dict[str, Any]], final_err: float, find_min: Optional[float]) -> str:
    if final_err <= WRONG:
        return "right"
    if not levels or levels[0]["disc_min"] > BASIN and levels[0]["post_min"] > BASIN:
        return "never"
    for L in levels:
        if L["kept_min"] > BASIN:
            return f"pruned {L['level']}"
    if find_min is not None and find_min <= BASIN:
        return "trap"
    return "other"


def main() -> None:
    cases = [tuple(x.split(":")) for x in (OUT / "newly_wrong_cases.txt").read_text().split()]
    rows: List[Dict[str, Any]] = []
    detail: Dict[str, Any] = {}
    for v, var in cases:
        Rt = np.load(CACHE / "e0_trace" / f"v{v}_{var}_on.npz")["R_true"]
        run: Dict[str, Dict[str, Any]] = {}
        for arm in ("off", "on"):
            z = np.load(CACHE / "e0_trace" / f"v{v}_{var}_{arm}.npz")
            lv = py_levels(z, Rt)
            fe = float(err(z["R_final"], Rt)[0])
            # FindOptimal candidates (inputs) of the last level = top N_FIND of the sorted list
            top = lv[-1]
            find_in = err(z[f"quick_mc_{lv[-1]['level']}_R"][:N_FIND], Rt)
            run[f"py_{arm}"] = dict(levels=lv, final_err=fe, find_min=float(find_in.min()),
                                    n_find=int(z["n_find"]),
                                    final_cost=float(z["cost_final"]))  # fmt: skip
        cf = CACHE / "e0_cpp" / f"v{v}_{var}.json"
        if cf.exists():
            d = json.loads(cf.read_text())
            if d.get("CV"):
                lv = cpp_levels(d, Rt)
                R = np.array(d["CV"][1:]).reshape(3, 3)
                fe = float(err(R, Rt)[0])
                fm = float(err(np.array(d["CQ"][str(lv[-1]["level"])])[:N_FIND, 1:], Rt).min())
                run["cpp"] = dict(levels=lv, final_err=fe, find_min=fm, n_find=len(d["CF"]),
                                  final_cost=float(d["CV"][0]))  # fmt: skip
        for k, r in run.items():
            r["class"] = classify(r["levels"], r["final_err"], r["find_min"])
        detail[f"{v}_{var}"] = run
        rows.append(dict(case=f"{v} {var}",
                         **{f"{k}_class": r["class"] for k, r in run.items()},
                         **{f"{k}_err": r["final_err"] for k, r in run.items()}))  # fmt: skip
    (OUT / "debug_cases.json").write_text(json.dumps(detail, indent=1, default=float))
    cnt: Dict[str, Dict[str, int]] = {}
    for r in rows:
        for k in ("py_off", "py_on", "cpp"):
            c = r.get(f"{k}_class")
            if c is not None:
                key = c.split(" ")[0] if c.startswith("pruned") else c
                lvl = c if c.startswith("pruned") else c
                cnt.setdefault(k, {})[lvl] = cnt.setdefault(k, {}).get(lvl, 0) + 1
    tbl = []
    classes = sorted({c for d in cnt.values() for c in d})
    for c in classes:
        tbl.append(dict(cls=c, **{k: cnt.get(k, {}).get(c, 0) for k in ("py_off", "py_on", "cpp")}))
    t = DT.markdown_table(tbl, ["cls", "py_off", "py_on", "cpp"])
    DT.write_tables(OUT / "tables_debug.md", {"e0_debug_classes": t})
    (OUT / "debug_summary.json").write_text(json.dumps(cnt, indent=1))
    print(t)


if __name__ == "__main__":
    main()


def per_level_table() -> None:
    """Medians over the newly wrong cases (n = 63): closest candidate per level, Python off / on
    and C++. Also the level-0 candidate-set match of Python on against C++ (same FZ grid, no random
    draws): matched candidates (within 0.01 deg) and the re-cost difference."""
    det = json.loads((OUT / "debug_cases.json").read_text())
    rows = []
    for k in ("py_off", "py_on", "cpp"):
        for lv in range(4):
            sel = [r[k]["levels"][lv] for r in det.values() if k in r and len(r[k]["levels"]) > lv]
            if not sel:
                continue
            f = lambda key: float(np.median([x[key] for x in sel]))  # noqa: E731
            rows.append(
                dict(run=k, level=lv, n=len(sel), n_cand=f("n_cand"), closest_deg=f("post_min"),
                     rank=f("rank"), n_keep=f("n_keep"),
                     survives=sum(x["survives"] for x in sel),
                     cost_closest=f("cost_of_closest"), cost_best=f("best_cost"),
                     within3=sum(x["post_min"] <= BASIN for x in sel))
            )  # fmt: skip
    t1 = DT.markdown_table(
        rows, list(rows[0].keys()),
        formats={"n_cand": ".0f", "closest_deg": ".2f", "rank": ".0f", "n_keep": ".0f",
                 "cost_closest": ".3f", "cost_best": ".3f"},
    )  # fmt: skip
    # level-0 set match
    m_rows = []
    for key in det:
        v, var = key.split("_")
        z = np.load(CACHE / "e0_trace" / f"v{v}_{var}_on.npz")
        zo = np.load(CACHE / "e0_trace" / f"v{v}_{var}_off.npz")
        d = json.loads((CACHE / "e0_cpp" / f"v{v}_{var}.json").read_text())
        Rc = np.array(d["CD"]["0"])
        cc, Rc = Rc[:, 0], Rc[:, 1:].reshape(-1, 3, 3)
        for tag, zz in (("on", z), ("off", zo)):
            Rp, cp = zz["discrete_0_R"], zz["discrete_0_score"]
            hit, dc = 0, []
            for i in range(len(Rc)):
                dd = np.abs(Rp - Rc[i]).reshape(len(Rp), -1).max(axis=1)
                j = int(dd.argmin())
                if dd[j] < 1e-5:
                    hit += 1
                    dc.append(cp[j] - cc[i])
            m_rows.append(dict(case=key, py=tag, n_cpp=len(Rc), n_py=len(Rp), matched=hit,
                               median_cost_diff=float(np.median(dc)) if dc else float("nan"),
                               ))  # fmt: skip
    mt = np.array(
        [[r["matched"] / r["n_cpp"] for r in m_rows if r["py"] == t] for t in ("on", "off")]
    )
    summ = {
        t: dict(
            cpp_cands_matched_frac_median=float(np.median(mt[i])),
            n_py_median=float(np.median([r["n_py"] for r in m_rows if r["py"] == t])),
            n_cpp_median=float(np.median([r["n_cpp"] for r in m_rows if r["py"] == t])),
            median_abs_cost_diff=float(np.nanmedian(
                [abs(r["median_cost_diff"]) for r in m_rows if r["py"] == t])),
        )
        for i, t in enumerate(("on", "off"))
    }  # fmt: skip
    (OUT / "debug_level0_match.json").write_text(
        json.dumps(dict(summary=summ, cases=m_rows), indent=1)
    )
    t2 = DT.markdown_table([dict(py=k, **v) for k, v in summ.items()],
                           ["py", "n_py_median", "n_cpp_median", "cpp_cands_matched_frac_median",
                            "median_abs_cost_diff"],
                           formats={"n_py_median": ".0f", "n_cpp_median": ".0f",
                                    "cpp_cands_matched_frac_median": ".3f",
                                    "median_abs_cost_diff": ".2e"})  # fmt: skip
    DT.write_tables(
        OUT / "tables_debug_levels.md", {"e0_debug_levels": t1, "e0_debug_level0_match": t2}
    )
    print(t1)
    print(t2)


if __name__ == "__main__":
    per_level_table()
