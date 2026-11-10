"""AR2 step 1: summary.json, tables.md, the trace of every wrong voxel, and the figure.

  uv run python scripts/autoresearch/ar2/holdout_summary.py
Reads benchmarks/autoresearch/ar2/holdout/run_part*.json (+ voxels.json, the preflight if present).
"""

import json
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[0] / "ar1"))
import ar1_common as AC  # noqa: E402
import misses as M  # noqa: E402

from icenine import ray_indexing as RI  # noqa: E402

import doc_tables as DT  # noqa: E402
import stats as ST  # noqa: E402

OUT = AC.ICE / "benchmarks" / "autoresearch" / "ar2" / "holdout"
ARMS = ("1", "2", "3")
# bins set after the range of the explained-ray counts was seen (30 to 57); they change no arm
BINS = [(0, 35, "35 or fewer"), (36, 40, "36 to 40"), (41, 45, "41 to 45"), (46, 50, "46 to 50"),
        (51, 10**6, "51 or more")]  # fmt: skip
KS = (1, 5, 20)
TOLS = (0.5, 2.0, 5.0)


def load_rows() -> List[Dict[str, Any]]:
    rows: List[Dict[str, Any]] = []
    for f in sorted(OUT.glob("run_part*.json")):
        rows += json.loads(f.read_text())["rows"]
    rows.sort(key=lambda r: r["voxel"])
    return rows


def recall_at(r: Dict[str, Any], k: int, tol: float) -> bool:
    d = r.get("dist", [])
    return bool(len(d) and min(d[:k]) < tol)


def wil(k: int, n: int) -> str:
    lo, hi = ST.wilson(k, n)
    return f"{k} of {n} ({100 * lo:.1f} to {100 * hi:.1f}%)"


def trace(rows: List[Dict[str, Any]], vox: Dict[str, Any]) -> List[Dict[str, Any]]:
    """For each voxel that is wrong in any arm: where the truth was, what was taken."""
    wrong = [r for r in rows if any(r.get("arms", {}).get(a, {}).get("wrong", True) for a in ARMS)]
    if not wrong:
        return []
    gmap = np.load(M.GRAINMAP)
    vxy, side = M.voxel_vertices_xy(M.MIC)
    AC.init("full")
    ctx = AC.bench().RUN._ENV["ctx"]
    first: Dict[int, int] = {}
    for v, g in enumerate(gmap):
        first.setdefault(int(g), v)
    gids = sorted(first)
    Rg = np.stack([ctx.truth(first[g]) for g in gids])
    e0 = AC.env("full", wrong[0]["voxel"], need_sc=False)
    sym = RI.proper_rotations(AC.symmetry_of(e0))
    out = []
    for r in wrong:
        v = r["voxel"]
        row: Dict[str, Any] = dict(voxel=v, set=r["set"], grain=int(gmap[v]))
        if r.get("no_candidates"):
            row.update(note="no candidates", wrong_arms="1,2,3")
            out.append(row)
            continue
        d = np.array(r["dist"])
        row["wrong_arms"] = ",".join(a for a in ARMS if r["arms"][a]["wrong"])
        hit = np.nonzero(d < 2.0)[0]
        row["truth_rank_in_20"] = int(hit[0]) + 1 if len(hit) else None
        row["truth_in_first5"] = bool(len(hit) and hit[0] < 5)
        fins = r["finishers"]
        refined = [f["cand"] for f in fins]
        row["refined_cands"] = refined  # in order: best 1, best 2, best 3
        tr = [f for f in fins if d[f["cand"]] < 2.0]  # the truth candidate, if refined
        row["truth_refined_in_arms"] = (
            ",".join(str(i + 1) for i, f in enumerate(fins) if d[f["cand"]] < 2.0) or "none"
        )
        if tr:
            row["truth_cost"] = [tr[0]["cost_before"], tr[0]["cost_after"]]
            row["truth_err_after_deg"] = tr[0]["err"]
        a_first = next(a for a in ARMS if r["arms"][a]["wrong"])
        taken = fins[[f["cand"] for f in fins].index(r["arms"][a_first]["taken_cand"])]
        Rf = np.array(taken["R_final"]).reshape(3, 3)
        dg = RI.misorientation_deg(Rf[None], Rg, sym)
        j = int(np.argmin(dg))
        row["arm_traced"] = a_first
        row["taken_cost"] = [taken["cost_before"], taken["cost_after"]]
        row["taken_err_deg"] = taken["err"]
        row["taken_dist_before"] = float(d[taken["cand"]])
        if dg[j] < 1.0:
            g = gids[j]
            row["taken_grain"] = g
            row["taken_touches"] = M.touch(v, g, gmap, vxy, side)["touches"]
        else:
            row["taken_grain"] = None
            row["taken_touches"] = None
        out.append(row)
    return out


def main() -> None:
    vox = json.loads((OUT / "voxels.json").read_text())
    rows = load_rows()
    summ: Dict[str, Any] = dict(label="contended", n_rows=len(rows))
    sets = {"boundary": [r for r in rows if r["set"] == "boundary"],
            "interior": [r for r in rows if r["set"] == "interior"]}  # fmt: skip
    sets["both"] = rows
    # (a) arms
    ta = []
    for name, rs in sets.items():
        for a in ARMS:
            wr = [bool(r.get("arms", {}).get(a, {}).get("wrong", True)) for r in rs]
            ev = np.array([r["arms"][a]["evals"] for r in rs if "arms" in r])
            er = np.array([r["arms"][a]["err"] for r in rs if "arms" in r])
            k = int(sum(wr))
            lo, hi = ST.wilson(k, len(rs))
            summ.setdefault("arms", {}).setdefault(name, {})[a] = dict(
                n=len(rs), wrong=k, wilson=[lo, hi], evals_mean=float(ev.mean()),
                evals_median=float(np.median(ev)), err_median=float(np.median(er)),
            )  # fmt: skip
            ta.append(
                dict(vset=name, arm=f"best {a}", wrong=wil(k, len(rs)), mean_evals=float(ev.mean()),
                     median_evals=float(np.median(ev)), median_err=float(np.median(er)))
            )  # fmt: skip
    # (b) recall
    tb = []
    for name, rs in sets.items():
        for tol in TOLS:
            cells = {f"at{k}": wil(sum(recall_at(r, k, tol) for r in rs), len(rs)) for k in KS}
            tb.append(dict(vset=name, tol=f"{tol:g} degrees", **cells))
            summ.setdefault("recall", {}).setdefault(name, {})[f"{tol:g}"] = {
                str(k): sum(recall_at(r, k, tol) for r in rs) for k in KS
            }
    # (d) recall against explained rays
    td = []
    for lo_, hi_, lab in BINS:
        rs = [r for r in rows if lo_ <= r["n_true_explained"] <= hi_]
        cells = {f"at{k}": f"{sum(recall_at(r, k, 2.0) for r in rs)} of {len(rs)}" for k in KS}
        td.append(dict(bin=lab, n=len(rs), **cells))
        summ.setdefault("recall_by_explained", {})[lab] = dict(
            n=len(rs), **{str(k): sum(recall_at(r, k, 2.0) for r in rs) for k in KS}
        )
    # (e) time
    te = []
    for name, rs in sets.items():
        t = np.array([r["index_s"] for r in rs])
        te.append(dict(vset=name, n=len(t), median_s=float(np.median(t)), mean_s=float(t.mean()),
                       p90_s=float(np.percentile(t, 90)), max_s=float(t.max())))  # fmt: skip
        summ.setdefault("index_s", {})[name] = dict(
            median=float(np.median(t)), mean=float(t.mean()), p90=float(np.percentile(t, 90)),
            max=float(t.max()),
        )  # fmt: skip
    # (c) trace
    tr = trace(rows, vox)
    summ["trace"] = tr
    tc = []
    for t in tr:

        def cc(c: Optional[List[float]]) -> str:
            return "n/a" if not c else f"{c[0]:.3f} to {c[1]:.3f}"

        tg = t.get("taken_grain")
        touch = {True: "touches", False: "does not touch", None: "n/a"}[t.get("taken_touches")]
        tc.append(
            dict(voxel=t["voxel"], vset=t["set"], wrong_in_best=t["wrong_arms"],
                 truth_rank=t.get("truth_rank_in_20") or "not in 20",
                 truth_in_first5="yes" if t.get("truth_in_first5") else "no",
                 truth_refined_in_best=t.get("truth_refined_in_arms", "n/a"),
                 taken=("no grain" if tg is None else f"grain {tg}, {touch}"),
                 taken_err=t.get("taken_err_deg", float("nan")),
                 taken_cost=cc(t.get("taken_cost")), truth_cost=cc(t.get("truth_cost")),
                 truth_err=t.get("truth_err_after_deg", float("nan")))
        )  # fmt: skip
    # the plan's rule
    w2, w3 = summ["arms"]["both"]["2"]["wrong"], summ["arms"]["both"]["3"]["wrong"]
    summ["rule"] = dict(best2_wrong=w2, best3_wrong=w3, main_arm="best 3" if w2 > w3 else "best 2")
    (OUT / "summary.json").write_text(json.dumps(summ, indent=1) + "\n")
    tables = {
        "ar2a_arms": DT.markdown_table(ta, ["vset", "arm", "wrong", "mean_evals", "median_evals",
                                            "median_err"], {"mean_evals": ".0f",
                                            "median_evals": ".0f", "median_err": ".4f"},
                                       ["set", "arm", "wrong (Wilson 95%)", "mean evals",
                                        "median evals", "median error (deg)"]),
        "ar2a_recall": DT.markdown_table(tb, ["vset", "tol", "at1", "at5", "at20"], None,
                                         ["set", "within", "recall at 1", "recall at 5",
                                          "recall at 20"]),
        "ar2a_trace": DT.markdown_table(
            tc, ["voxel", "vset", "wrong_in_best", "truth_rank", "truth_in_first5",
                 "truth_refined_in_best", "taken", "taken_err", "taken_cost", "truth_cost",
                 "truth_err"],
            {"taken_err": ".3f", "truth_err": ".3f"},
            ["voxel", "set", "wrong in best", "truth rank (of 20)", "truth in first 5",
             "truth refined in best", "taken (first wrong arm)", "taken error (deg)",
             "taken cost before to after", "truth cost before to after",
             "truth error after finisher (deg)"],
        ) if tc else "No voxel is wrong in any arm.\n",
        "ar2a_explained": DT.markdown_table(td, ["bin", "n", "at1", "at5", "at20"], None,
                                            ["rays the truth explains", "voxels", "recall at 1",
                                             "recall at 5", "recall at 20"]),
        "ar2a_time": DT.markdown_table(te, ["vset", "n", "median_s", "mean_s", "p90_s", "max_s"],
                                       {"median_s": ".2f", "mean_s": ".2f", "p90_s": ".2f",
                                        "max_s": ".2f"},
                                       ["set", "voxels", "median (s)", "mean (s)", "90th pct (s)",
                                        "max (s)"]),
    }  # fmt: skip
    DT.write_tables(OUT / "tables.md", tables)
    figure(summ)
    print(json.dumps(summ["rule"]), "wrong voxels", len(tr))


def figure(summ: Dict[str, Any]) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(1, 2, figsize=(11, 4))
    w = 0.25
    for i, a in enumerate(ARMS):
        for j, name in enumerate(("boundary", "interior")):
            d = summ["arms"][name][a]
            lo, hi = d["wilson"]
            f = d["wrong"] / d["n"]
            ax[0].bar(j + (i - 1) * w, 100 * f, w * 0.9, color=f"C{i}",
                      label=f"best {a}" if j == 0 else None,
                      yerr=[[100 * (f - lo)], [100 * (hi - f)]], capsize=3)  # fmt: skip
            ax[0].text(j + (i - 1) * w, 100 * hi + 0.3, f"{d['wrong']}", ha="center", fontsize=8)
    ax[0].set_xticks([0, 1], ["boundary (150)", "interior (150)"])
    ax[0].set_ylabel("voxels wrong, % (bar: Wilson 95%)")
    ax[0].set_ylim(0, 9)
    ax[0].set_title("Wrong voxels per finisher arm (count above bar)")
    ax[0].legend()
    labs = list(summ["recall_by_explained"])
    x = np.arange(len(labs))
    for i, k in enumerate(KS):
        fr = [
            100
            * summ["recall_by_explained"][b][str(k)]
            / max(summ["recall_by_explained"][b]["n"], 1)
            for b in labs
        ]
        ax[1].plot(x, fr, "o-", label=f"recall at {k}")
    ax[1].set_xticks(x, [f"{b}\n(n={summ['recall_by_explained'][b]['n']})" for b in labs],
                     fontsize=8)  # fmt: skip
    ax[1].set_xlabel("attached rays the true orientation explains")
    ax[1].set_ylabel("truth within 2 degrees, %")
    ax[1].set_title("Recall against rays explained")
    ax[1].legend()
    fig.tight_layout()
    fig.savefig(OUT / "ar2a_holdout.png", dpi=130)


if __name__ == "__main__":
    main()
