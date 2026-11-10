"""AR2 step 1: summary.json, tables.md, the trace of every wrong voxel, and the figure.

  uv run python scripts/autoresearch/ar2/holdout_summary.py
Reads benchmarks/autoresearch/ar2/holdout/run_part*.json, cand_part*.npz and voxels.json.

The distances and final errors are recomputed here from the stored candidates and final
orientations with the current ``ray_indexing.misorientation_deg`` (review fix 3: the stored values
used an ``arccos`` with a floor of a few thousandths of a degree). The wrong flags and the recall
counts must not change; the script stops if they do. For every traced voxel it evaluates the
strict cost at the exact true orientation and re-evaluates the stored final costs (fix 1). It also
splits the boundary set into voxels that share an edge with another grain and voxels that touch it
only at a corner (fix 4), stores the grain counts per set (fix 9) and the paired splits (fix 5).
"""

import json
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional, Set, Tuple

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
WRONG_DEG = 1.0
GRAIN_TOL_DEG = 1.0  # as select_voxels.py


def load_rows() -> List[Dict[str, Any]]:
    rows: List[Dict[str, Any]] = []
    for f in sorted(OUT.glob("run_part*.json")):
        rows += json.loads(f.read_text())["rows"]
    rows.sort(key=lambda r: r["voxel"])
    return rows


def load_cands() -> Dict[int, np.ndarray]:
    cand: Dict[int, np.ndarray] = {}
    for f in sorted(OUT.glob("cand_part*.npz")):
        z = np.load(f)
        cand.update({int(k[1:]): z[k] for k in z.files})
    return cand


def recall_at(r: Dict[str, Any], k: int, tol: float) -> bool:
    d = r.get("dist", [])
    return bool(len(d) and min(d[:k]) < tol)


def wil(k: int, n: int) -> str:
    lo, hi = ST.wilson(k, n)
    return f"{k} of {n} ({100 * lo:.1f} to {100 * hi:.1f}%)"


def recompute_errors(
    rows: List[Dict[str, Any]], truth: Dict[int, np.ndarray], sym: np.ndarray
) -> Dict[str, Any]:
    """Replace the stored distances and final errors (old ``arccos`` form) by the robust angle,
    in place. Returns the comparison; raises if a wrong flag, a taken candidate or a recall count
    would change."""
    cand = load_cands()
    old_err: List[float] = []
    new_err: List[float] = []
    old_d: List[float] = []
    new_d: List[float] = []
    flips = {f"{t:g}": 0 for t in TOLS}
    for r in rows:
        v = r["voxel"]
        Rt = truth[v]
        d_new = RI.misorientation_deg(cand[v], Rt[None], sym)
        for t in TOLS:
            flips[f"{t:g}"] += int(
                sum(recall_at(dict(dist=r["dist"]), k, t) != recall_at(dict(dist=d_new), k, t)
                    for k in KS)
            )  # fmt: skip
        old_d += r["dist"]
        new_d += d_new.tolist()
        r["dist"] = [float(x) for x in d_new]
        for f in r["finishers"]:
            Rf = np.array(f["R_final"]).reshape(3, 3)
            e = float(RI.misorientation_deg(Rf[None], Rt[None], sym)[0])
            old_err.append(f["err"])
            new_err.append(e)
            f["err"] = e
        for a in ARMS:
            arm = r["arms"][a]
            f = next(f for f in r["finishers"] if f["cand"] == arm["taken_cand"])
            if bool(f["err"] > WRONG_DEG) != arm["wrong"]:
                raise RuntimeError(f"voxel {v} arm {a}: the wrong flag changes")
            arm["err"] = f["err"]
    if any(flips.values()):
        raise RuntimeError(f"recall counts change: {flips}")
    oe, ne = np.array(old_err), np.array(new_err)
    od, nd = np.array(old_d), np.array(new_d)
    small = od < 1.0
    return dict(
        finisher_runs=int(len(oe)),
        err_exact_zero_old=int((oe == 0.0).sum()),
        err_exact_zero_new=int((ne == 0.0).sum()),
        err_max_abs_change=float(np.abs(ne - oe).max()),
        err_median_old=float(np.median(oe[oe < WRONG_DEG])),
        err_median_new=float(np.median(ne[ne < WRONG_DEG])),
        dist_n=int(len(od)),
        dist_max_abs_change_below_1deg=float(np.abs(nd - od)[small].max()),
        dist_max_abs_change_above_1deg=float(np.abs(nd - od)[~small].max()),
        recall_flips=flips,
        wrong_flags_changed=0,
    )


def neighbours(vxy: np.ndarray, side: np.ndarray) -> Tuple[List[Set[int]], List[Set[int]]]:
    """Per voxel: the voxels sharing at least one vertex, and those sharing two (an edge)."""
    q = float(np.median(side)) * 1e-3
    keys = np.rint(vxy / q).astype(np.int64)  # (N, 3, 2)
    owners: Dict[Tuple[int, int], List[int]] = {}
    for v in range(len(vxy)):
        for k in range(3):
            owners.setdefault((int(keys[v, k, 0]), int(keys[v, k, 1])), []).append(v)
    count: List[Dict[int, int]] = [dict() for _ in range(len(vxy))]
    for vs in owners.values():
        for a in vs:
            for b in vs:
                if b != a:
                    count[a][b] = count[a].get(b, 0) + 1
    corner = [set(c) for c in count]
    edge = [set(b for b, n in c.items() if n >= 2) for c in count]
    return corner, edge


def boundary_split(
    vox: Dict[str, Any], gmap: np.ndarray, Rg: Dict[int, np.ndarray], sym: np.ndarray
) -> Tuple[Dict[int, str], Dict[str, Any]]:
    """Kind of each holdout boundary voxel ("edge" or "corner") and the layer counts."""
    vxy, side = M.voxel_vertices_xy(M.MIC)
    corner, edge = neighbours(vxy, side)
    far_cache: Dict[Tuple[int, int], bool] = {}

    def far(a: int, b: int) -> bool:
        if a == b:
            return False
        k = (min(a, b), max(a, b))
        if k not in far_cache:
            far_cache[k] = bool(
                RI.misorientation_deg(Rg[a][None], Rg[b][None], sym)[0] > GRAIN_TOL_DEG
            )
        return far_cache[k]

    def has_far(v: int, nb: List[Set[int]]) -> bool:
        return any(far(int(gmap[v]), int(gmap[u])) for u in nb[v])

    region = np.load(AC.ICE / "benchmarks" / "phase_d_pilot" / "region_index.npy")
    in_region = np.zeros(len(gmap), dtype=bool)
    in_region[region] = True
    elig = np.nonzero(~in_region)[0]
    corner_rule = int(sum(has_far(int(v), corner) for v in elig))
    edge_rule = int(sum(has_far(int(v), edge) for v in elig))
    kind: Dict[int, str] = {}
    for v in vox["sets"]["boundary"]["voxels"]:
        assert has_far(v, corner)
        kind[v] = "edge" if has_far(v, edge) else "corner"
    n_edge = sum(k == "edge" for k in kind.values())
    layer = dict(
        eligible=int(len(elig)),
        boundary_corner_rule=corner_rule,
        boundary_edge_rule=edge_rule,
        holdout_boundary=len(kind),
        holdout_edge=n_edge,
        holdout_corner_only=len(kind) - n_edge,
    )
    return kind, layer


def trace(
    rows: List[Dict[str, Any]], kind: Dict[int, str], gmap: np.ndarray, Rg: Dict[int, np.ndarray],
    sym: np.ndarray,
) -> List[Dict[str, Any]]:  # fmt: skip
    """For each voxel that is wrong in any arm: where the truth was, what was taken, the strict
    cost at the exact truth, and the stored final costs re-evaluated."""
    wrong = [r for r in rows if any(r["arms"][a]["wrong"] for a in ARMS)]
    if not wrong:
        return []
    vxy, side = M.voxel_vertices_xy(M.MIC)
    cand = load_cands()
    gids = sorted(Rg)
    Rg_all = np.stack([Rg[g] for g in gids])
    out = []
    for r in wrong:
        v = r["voxel"]
        e = AC.env("full", v, need_sc=False)
        lf, verts, phase, Rt = e["lf"], e["verts"], e["phase"], e["R_true"]
        row: Dict[str, Any] = dict(voxel=v, set=r["set"], grain=int(gmap[v]))
        row["boundary_kind"] = kind.get(v, "interior")
        d = np.array(r["dist"])
        row["wrong_arms"] = ",".join(a for a in ARMS if r["arms"][a]["wrong"])
        hit = np.nonzero(d < 2.0)[0]
        row["truth_rank_in_20"] = int(hit[0]) + 1 if len(hit) else None
        row["truth_in_first5"] = bool(len(hit) and hit[0] < 5)
        fins = r["finishers"]
        row["refined_cands"] = [f["cand"] for f in fins]  # in order: best 1, best 2, best 3
        # every finisher position (1-based) whose candidate is within 2 degrees of the truth
        row["truth_refined_in_arms"] = (
            ",".join(str(i + 1) for i, f in enumerate(fins) if d[f["cand"]] < 2.0) or "none"
        )
        row["strict_cost_at_truth"] = float(lf.evaluate(Rt, verts, phase).cost)
        tr = [f for f in fins if d[f["cand"]] < 2.0]  # the truth candidate, if refined
        if tr:
            f = tr[0]
            Rf = np.array(f["R_final"]).reshape(3, 3)
            row["truth_cost"] = [f["cost_before"], f["cost_after"]]
            row["truth_cost_after_reeval"] = float(lf.evaluate(Rf, verts, phase).cost)
            row["truth_start_deg"] = float(d[f["cand"]])
            row["truth_err_after_deg"] = f["err"]
            row["truth_moved_deg"] = float(
                RI.misorientation_deg(Rf[None], cand[v][f["cand"]][None], sym)[0]
            )
            row["truth_finisher_evals"] = f["evals"]
        a_first = next(a for a in ARMS if r["arms"][a]["wrong"])
        taken = next(f for f in fins if f["cand"] == r["arms"][a_first]["taken_cand"])
        Rf = np.array(taken["R_final"]).reshape(3, 3)
        dg = RI.misorientation_deg(Rf[None], Rg_all, sym)
        j = int(np.argmin(dg))
        row["arm_traced"] = a_first
        row["taken_cost"] = [taken["cost_before"], taken["cost_after"]]
        row["taken_cost_after_reeval"] = float(lf.evaluate(Rf, verts, phase).cost)
        row["taken_err_deg"] = taken["err"]
        row["taken_dist_before"] = float(d[taken["cand"]])
        if dg[j] < WRONG_DEG:
            g = gids[j]
            row["taken_grain"] = g
            row["taken_touches"] = M.touch(v, g, gmap, vxy, side)["touches"]
        else:
            row["taken_grain"] = None
            row["taken_touches"] = None
        # the stored final costs must reproduce on unchanged code (stop rule of the plan)
        for k in ("taken", "truth"):
            if f"{k}_cost_after_reeval" in row:
                stored, again = row[f"{k}_cost"][1], row[f"{k}_cost_after_reeval"]
                if abs(stored - again) > 1e-9:
                    raise RuntimeError(f"voxel {v}: stored {k} cost {stored} != {again}")
        out.append(row)
    return out


def paired(rows: List[Dict[str, Any]], a: str, b: str) -> Dict[str, int]:
    """Arm ``b`` against arm ``a`` by voxel: fixed (wrong in a only), broken (wrong in b only),
    both wrong, and the voxels where both arms take the same candidate."""
    wa = [r["arms"][a]["wrong"] for r in rows]
    wb = [r["arms"][b]["wrong"] for r in rows]
    fixed, broken = ST.paired_discordant(wa, wb)
    same = sum(r["arms"][a]["taken_cand"] == r["arms"][b]["taken_cand"] for r in rows)
    both = int(sum(x and y for x, y in zip(wa, wb)))
    return dict(n=len(rows), fixed=fixed, broken=broken, both_wrong=both, same_candidate=int(same))


def main() -> None:
    vox = json.loads((OUT / "voxels.json").read_text())
    rows = load_rows()
    no_cand = [r["voxel"] for r in rows if r.get("no_candidates")]
    # every count below assumes each voxel has candidates and arms (true in this run)
    assert not no_cand, f"voxels without candidates: {no_cand}"
    gmap = np.load(M.GRAINMAP)
    AC.init("full")
    ctx = AC.bench().RUN._ENV["ctx"]
    first: Dict[int, int] = {}
    for v, g in enumerate(gmap):
        first.setdefault(int(g), v)
    Rg = {g: ctx.truth(v) for g, v in first.items()}
    e0 = AC.env("full", rows[0]["voxel"], need_sc=False)
    sym = RI.proper_rotations(AC.symmetry_of(e0))
    truth = {r["voxel"]: AC.env("full", r["voxel"], need_sc=False)["R_true"] for r in rows}
    summ: Dict[str, Any] = dict(label="contended", n_rows=len(rows))
    summ["error_metric_fix3"] = recompute_errors(rows, truth, sym)
    summ["grains"] = {
        name: dict(grains=len(set(s["grain"])), in_pilot_region=sum(s["grain_in_pilot_region"]))
        for name, s in vox["sets"].items()
    }
    kind, layer = boundary_split(vox, gmap, Rg, sym)
    summ["boundary_split"] = dict(layer=layer)
    sets = {"boundary": [r for r in rows if r["set"] == "boundary"],
            "interior": [r for r in rows if r["set"] == "interior"]}  # fmt: skip
    sets["both"] = rows
    # (a) arms
    ta = []
    for name, rs in sets.items():
        for a in ARMS:
            wr = [bool(r["arms"][a]["wrong"]) for r in rs]
            ev = np.array([r["arms"][a]["evals"] for r in rs])
            er = np.array([r["arms"][a]["err"] for r in rs])
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
    # (a2) the boundary set split into edge-sharing and corner-only voxels
    tsplit = []
    for sub in ("edge", "corner"):
        rs = [r for r in sets["boundary"] if kind[r["voxel"]] == sub]
        for a in ARMS:
            k = int(sum(r["arms"][a]["wrong"] for r in rs))
            lo, hi = ST.wilson(k, len(rs))
            summ["boundary_split"].setdefault(sub, {})[a] = dict(
                n=len(rs), wrong=k, wilson=[lo, hi]
            )
            tsplit.append(dict(sub=sub, arm=f"best {a}", wrong=wil(k, len(rs))))
        r1 = [r for r in rs if recall_at(r, 1, 2.0)]
        summ["boundary_split"][sub]["recall_at1_2deg"] = len(r1)
    # (a3) paired splits by voxel
    summ["paired"] = {"best1_vs_best2": paired(rows, "1", "2"),
                      "best2_vs_best3": paired(rows, "2", "3")}  # fmt: skip
    tpair = [dict(pair=f"best {b} against best {a}", **summ["paired"][f"best{a}_vs_best{b}"])
             for a, b in (("1", "2"), ("2", "3"))]  # fmt: skip
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
    tr = trace(rows, kind, gmap, Rg, sym)
    summ["trace"] = tr
    summ["boundary_split"]["wrong_any_arm"] = {
        sub: sum(1 for t in tr if t["boundary_kind"] == sub) for sub in ("edge", "corner")
    }
    tc = []
    for t in tr:

        def cc(c: Optional[List[float]]) -> str:
            return "n/a" if not c else f"{c[0]:.3f} to {c[1]:.3f}"

        tg = t.get("taken_grain")
        touch = {True: "touches", False: "does not touch", None: "n/a"}[t.get("taken_touches")]
        tc.append(
            dict(voxel=t["voxel"], vset=t["set"], boundary=t["boundary_kind"],
                 wrong_in_best=t["wrong_arms"],
                 truth_rank=t.get("truth_rank_in_20") or "not in 20",
                 truth_in_first5="yes" if t.get("truth_in_first5") else "no",
                 truth_refined_in_best=t.get("truth_refined_in_arms", "n/a"),
                 taken=("no grain" if tg is None else f"grain {tg}, {touch}"),
                 taken_err=t.get("taken_err_deg", float("nan")),
                 taken_cost=cc(t.get("taken_cost")), truth_cost=cc(t.get("truth_cost")),
                 cost_at_truth=t["strict_cost_at_truth"],
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
        "ar2a_split": DT.markdown_table(tsplit, ["sub", "arm", "wrong"]),
        "ar2a_paired": DT.markdown_table(tpair, ["pair", "n", "fixed", "broken", "both_wrong",
                                                 "same_candidate"]),
        "ar2a_recall": DT.markdown_table(tb, ["vset", "tol", "at1", "at5", "at20"], None,
                                         ["set", "within", "recall at 1", "recall at 5",
                                          "recall at 20"]),
        "ar2a_trace": DT.markdown_table(
            tc, ["voxel", "vset", "boundary", "wrong_in_best", "truth_rank", "truth_in_first5",
                 "truth_refined_in_best", "taken", "taken_err", "taken_cost", "truth_cost",
                 "cost_at_truth", "truth_err"],
            {"taken_err": ".3f", "truth_err": ".4f", "cost_at_truth": ".4f"},
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
    print(json.dumps({k: summ[k] for k in ("error_metric_fix3", "grains", "paired")}))
    print(json.dumps(summ["boundary_split"]))


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
    ax[1].set_xlabel("attached rays the true orientation explains (bins set after the counts)")
    ax[1].set_ylabel("truth within 2 degrees, %")
    ax[1].set_title("Recall against rays explained")
    ax[1].legend()
    fig.tight_layout()
    fig.savefig(OUT / "ar2a_holdout.png", dpi=130)


if __name__ == "__main__":
    main()
