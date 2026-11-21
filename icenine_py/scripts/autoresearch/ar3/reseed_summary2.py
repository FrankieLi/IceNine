"""AR3 step 8: tables and plot for parts 1.2, 1.3, 2 and 3 (called by ``reseed_summary.py``).

Reads only committed files under ``reseed/slim``. Voxels are counted at their first appearance in a
reseed cluster (a voxel can be listed in several clusters, because a cluster that is reseeded
leaves its unsettled voxels to form the next one); the number of cluster entries is given beside.
"""

import gzip
import json
from pathlib import Path
from typing import Any, Dict, List, Tuple

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from doc_tables import markdown_table  # noqa: E402
from stats import wilson  # noqa: E402

HERE = Path(__file__).resolve().parent
OUT = HERE.parents[2] / "benchmarks" / "autoresearch" / "ar3" / "reseed"
SLIM = OUT / "slim"
RIGHT = 1.0
BINS = [0.0, 0.1, 0.2, 0.5, np.inf]


def rj(name: str) -> Any:
    p = SLIM / name
    if name.endswith(".gz"):
        with gzip.open(p, "rt") as f:
            return json.load(f)
    return json.loads(p.read_text())


def med(x: List[float]) -> float:
    return float(np.median(x)) if len(x) else float("nan")


def ci(k: int, n: int, d: int = 1) -> str:
    lo, hi = wilson(k, n)
    return f"{100 * lo:.{d}f} to {100 * hi:.{d}f}"


# ------------------------------------------------------------------------------ parts 1.2, 1.3
def classes() -> Tuple[Dict[str, Any], Dict[str, str]]:
    mem = rj("members_big.json.gz")
    first: Dict[int, Dict[str, Any]] = {}
    for m in mem:
        first.setdefault(m["voxel"], m)
    rows: List[Dict[str, Any]] = []
    out: Dict[str, Any] = dict(entries=len(mem), voxels=len(first))
    out["entries_by_class"] = {c_: sum(1 for m in mem if m["cls"] == c_) for c_ in "ABCN"}
    for L in [1.0, 2.0, 4.0, None]:
        sel = [m for m in first.values() if L is None or m["level"] == L]
        if not sel:
            continue
        n = len(sel)
        a = [m for m in sel if m["cls"] == "A"]
        b = [m for m in sel if m["cls"] == "B"]
        c = [m for m in sel if m["cls"] == "C"]
        nn = [m for m in sel if m["cls"] == "N"]
        bb = [int(np.digitize(m["err_start"], BINS[1:-1])) for m in b]
        d = dict(
            level="all" if L is None else f"{L:g}",
            voxels=n,
            A=len(a),
            B1=bb.count(0),
            B2=bb.count(1),
            B3=bb.count(2),
            B4=bb.count(3),
            C=len(c),
            N=len(nn),
            A_share=100.0 * len(a) / n,
            A_ci=ci(len(a), n),
            AB_share=100.0 * (len(a) + len(b)) / n,
            hit_A=med([m["hit"] for m in a]),
            best_A=med([m["best"] for m in a]),
            hit_B=med([m["hit"] for m in b]),
            hit_C=med([m["hit"] for m in c]),
            err_A=med([m["err_end"] for m in a]),
            start_A=med([m["err_start"] for m in a]),
            dist_A=med([m["dist"] for m in a]),
            dist_C=med([m["dist"] for m in c]),
            start_B=med([m["err_start"] for m in b]),
            err_C=med([m["err_end"] for m in c]),
        )
        rows.append(d)
        out[d["level"]] = d
    tab = markdown_table(
        rows,
        ["level", "voxels", "A", "B1", "B2", "B3", "B4", "C", "N", "A_share", "A_ci", "AB_share"],
        formats=dict(A_share=".1f", AB_share=".1f"),
    )
    tab2 = markdown_table(
        rows,
        [
            "level",
            "hit_A",
            "best_A",
            "err_A",
            "start_A",
            "dist_A",
            "hit_B",
            "start_B",
            "hit_C",
            "err_C",
            "dist_C",
        ],
        formats=dict(
            hit_A=".3f",
            best_A=".3f",
            err_A=".3f",
            start_A=".3f",
            dist_A=".1f",
            hit_B=".3f",
            hit_C=".3f",
            start_B=".1f",
            err_C=".1f",
            dist_C=".1f",
        ),
    )
    # ---- plot
    fig, ax = plt.subplots(figsize=(6.5, 4))
    lv = [d for d in rows if d["level"] != "all"]
    x = np.arange(len(lv))
    bottom = np.zeros(len(lv))
    for key, col, lab in (
        ("A", "#1b9e77", "A: reached, ended right, rejected"),
        ("B", "#d95f02", "B: own grain's start, did not end right"),
        ("C", "#7570b3", "C: start from another grain"),
    ):
        v = np.array(
            [d["B1"] + d["B2"] + d["B3"] + d["B4"] if key == "B" else d[key] for d in lv], float
        )
        ax.bar(x, v, bottom=bottom, color=col, label=lab)
        for xi, (bo, vi) in enumerate(zip(bottom, v)):
            if vi > 0:
                ax.text(xi, bo + vi / 2, f"{int(vi)}", ha="center", va="center", fontsize=8)
        bottom += v
    ax.set_xticks(x)
    ax.set_xticklabels([f"{d['level']} degrees" for d in lv])
    ax.set_ylabel("voxels in reseed clusters (first appearance)")
    ax.set_xlabel("gradient level of the true grain")
    ax.set_title("big, run S5: class of the last rejected fit")
    ax.legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(OUT / "ar3_reseed_classes.png", dpi=110)
    plt.close(fig)
    return {"classes": out}, {"classes_by_level": tab, "classes_numbers": tab2}


def clusters() -> Tuple[Dict[str, Any], Dict[str, str]]:
    c = rj("clusters_big.json.gz")
    tot_seed = sum(x["seed_evals"] for x in c)
    tot_all = sum(x["seed_evals"] + x["flood_evals"] for x in c)
    rows: List[Dict[str, Any]] = []
    out: Dict[str, Any] = dict(total_seed_evals=tot_seed, total_with_flood=tot_all)
    for L in [1.0, 2.0, 4.0, None]:
        s = [x for x in c if L is None or x["level"] == L]
        if not s:
            continue
        sz = [x["size"] for x in s]
        d = dict(
            level="all" if L is None else f"{L:g}",
            clusters=len(s),
            single=sum(1 for z in sz if z == 1),
            size_med=med(sz),
            size_max=max(sz),
            dist_med=med([x["dist"] for x in s]),
            settled_med=med([x["settled"] for x in s]),
            members=sum(sz),
            settled_sum=sum(x["settled"] for x in s),
            flood_acc_med=med([x["flood_accepted"] for x in s]),
            seed_evals=sum(x["seed_evals"] for x in s),
            flood_evals=sum(x["flood_evals"] for x in s),
            share=100.0 * sum(x["seed_evals"] + x["flood_evals"] for x in s) / tot_all,
            two_grains=sum(1 for x in s if x["n_grains"] > 1),
        )
        rows.append(d)
        out[d["level"]] = d
    tab = markdown_table(
        rows,
        [
            "level",
            "clusters",
            "single",
            "size_med",
            "size_max",
            "dist_med",
            "settled_med",
            "members",
            "settled_sum",
            "flood_acc_med",
            "seed_evals",
            "flood_evals",
            "share",
            "two_grains",
        ],
        formats=dict(
            size_med=".0f",
            dist_med=".1f",
            settled_med=".0f",
            flood_acc_med=".0f",
            seed_evals="d",
            flood_evals="d",
            share=".1f",
        ),
    )
    # cost share of a class: the reseed cost (seed step plus flood) of the clusters in which that
    # class has the most voxels (judged on the cluster's member rows; ties to the lower letter).
    # Chosen after the classes were seen, to read the plan's "cost share" for part 2.
    from collections import Counter

    mem = rj("members_big.json.gz")
    by: Dict[int, List[str]] = {}
    for m in mem:
        by.setdefault(m["k"], []).append(m["cls"])
    plu: Dict[str, List[int]] = {c_: [0, 0] for c_ in ("A", "B", "C")}
    for x in c:
        cnt = Counter(by[x["k"]])
        top = max(cnt.values())
        cls = sorted(k for k, v in cnt.items() if v == top)[0]
        plu[cls][0] += 1
        plu[cls][1] += x["seed_evals"] + x["flood_evals"]
    out["plurality"] = {
        k: dict(clusters=v[0], evals=v[1], share=100.0 * v[1] / tot_all) for k, v in plu.items()
    }
    rows_p = [dict(cls=k, **out["plurality"][k]) for k in ("A", "B", "C")]
    tab_p = markdown_table(
        rows_p, ["cls", "clusters", "evals", "share"], formats=dict(evals="d", share=".1f")
    )
    return {"clusters": out}, {"clusters_by_level": tab, "clusters_by_class": tab_p}


# --------------------------------------------------------------------------------------- part 2
def remedies() -> Tuple[Dict[str, Any], Dict[str, str]]:
    e1 = rj("../remedies_e1.json")
    pk = rj("../remedies_pick.json")
    e2 = rj("../remedies_e2.json")
    out: Dict[str, Any] = dict(e1=e1, pick=pk, e2=e2)
    rows = []
    for t, d in pk["pick_pooled_S5_both_samples"].items():
        rows.append(
            dict(
                tau=float(t),
                n=d["n"],
                stored=d["evals_stored"],
                e1=d["evals_e1"],
                cut=100 * d["cut"],
                fixed=d["fixed"],
                broken=d["broken"],
                chosen="yes" if float(t) == pk["chosen_tau"] else "",
            )
        )
    tab_pick = markdown_table(
        rows,
        ["tau", "n", "stored", "e1", "cut", "fixed", "broken", "chosen"],
        formats=dict(tau=".2f", stored="d", e1="d", cut=".1f"),
    )
    names = {"big/S5L": "big S5", "big/D4": "big D4", "grad/S5": "grad S5", "grad/D4": "grad D4"}
    rows = []
    for key, nm in names.items():
        for tag, lab in (("tau0.80", "0.80"), ("tau_chosen", f"{pk['chosen_tau']:.2f}")):
            for st in ("pick", "test"):
                d = e1["results"][f"{key}/{tag}/{st}"]
                rows.append(
                    dict(
                        run=nm,
                        tau=lab,
                        set=st,
                        n=d["n"],
                        stored=d["evals_stored"],
                        e1=d["evals_e1"],
                        cut=100 * d["cut"],
                        fixed=d["fixed"],
                        broken=d["broken"],
                        both_wrong=d["both_wrong"],
                        switched=d["switched"],
                    )
                )
    tab_e1 = markdown_table(
        rows,
        [
            "run",
            "tau",
            "set",
            "n",
            "stored",
            "e1",
            "cut",
            "fixed",
            "broken",
            "both_wrong",
            "switched",
        ],
        formats=dict(stored="d", e1="d", cut=".1f"),
    )
    rows = []
    for n_, d in e2["pick"].items():
        rows.append(
            dict(
                test=n_,
                cases=d["cases"],
                kept=d["kept"],
                wrong=d["settled_wrong"],
                fixed=d["fixed"],
                broken=d["broken"],
                gone=d["clusters_gone"],
                clusters=d["clusters"],
                cut=100 * d["cut"],
            )
        )
    tab_e2 = markdown_table(
        rows,
        ["test", "cases", "kept", "wrong", "fixed", "broken", "gone", "clusters", "cut"],
        formats=dict(cut=".1f"),
    )
    return out, {"e1_pick": tab_pick, "e1_sets": tab_e1, "e2_pick": tab_e2}


# --------------------------------------------------------------------------------------- part 3
def part3() -> Tuple[Dict[str, Any], Dict[str, str]]:
    out: Dict[str, Any] = {}
    tabs: Dict[str, str] = {}
    stage_rows: List[Dict[str, Any]] = []
    wrong_rows: List[Dict[str, Any]] = []
    pair_rows: List[Dict[str, Any]] = []
    for smp in ("big", "grad"):
        if not (SLIM / f"s8_{smp}.json").exists():
            continue
        s8 = rj(f"s8_{smp}.json")
        base = rj(f"baselines_{smp}.json")
        z = np.load(SLIM / f"errors_{smp}.npz")
        runs = {"D4": base["D4"], "S5": base["S5"], "S8": s8}
        res: Dict[str, Any] = {}
        for k, r in runs.items():
            s = r["stats"]
            pe = int(r["pass"]["evals"])
            st = dict(
                seed_step=int(s["index_evals"]),
                grid_fallback=int(s.get("index_fallback_grid_evals", 0)),
                neighbour=int(s["n_evals_neighbor"]),
                revisit=int(s["n_evals_revisit"]),
                pass_=pe,
            )
            tot = sum(st.values())
            assert tot - pe == int(r["evals_total"]), (smp, k, tot - pe, r["evals_total"])
            before = z[f"{k}_before"] > RIGHT
            after = z[f"{k}_after"] > RIGHT
            n = len(before)
            res[k] = dict(
                stages=st,
                total=tot,
                seeds=int(s["n_seeds"]),
                reseeds=int(s["n_reseeds"]),
                fallbacks=int(
                    s.get("index_fallback_seed_test", 0) + s.get("index_fallback_no_candidate", 0)
                ),
                wrong_before=int(before.sum()),
                wrong_after=int(after.sum()),
                n=n,
                ci_after=wilson(int(after.sum()), n),
                ci_before=wilson(int(before.sum()), n),
            )
            stage_rows.append(
                dict(
                    sample=smp,
                    run=k + " + pass",
                    seeds=res[k]["seeds"],
                    reseeds=res[k]["reseeds"],
                    grid=res[k]["fallbacks"],
                    **st,
                    total=tot,
                )
            )
            wrong_rows.append(
                dict(
                    sample=smp,
                    run=k,
                    before=f"{int(before.sum())}",
                    before_ci=ci(int(before.sum()), n, 2),
                    after=f"{int(after.sum())}",
                    after_ci=ci(int(after.sum()), n, 2),
                )
            )
        for k in ("S5", "D4"):
            for when in ("before", "after"):
                a = z[f"S8_{when}"] > RIGHT
                b = z[f"{k}_{when}"] > RIGHT
                pair_rows.append(
                    dict(
                        sample=smp,
                        vs=k,
                        when=when,
                        wrong_S8=int(a.sum()),
                        wrong_base=int(b.sum()),
                        fixed=int((b & ~a).sum()),
                        broken=int((a & ~b).sum()),
                        both_wrong=int((a & b).sum()),
                    )
                )
        res["e1"] = s8["e1"]
        res["s8_wall_s"] = s8["wall_total_s"]
        res["timing_label"] = s8["timing_label"]
        res["pass_wall_s"] = s8["pass"].get("wall_tries_s")
        out[smp] = res
        # by level (after the pass)
        lv = z["level"]
        rows = []
        for L in sorted(set(lv.tolist())):
            m = lv == L
            d = dict(sample=smp, level=f"{L:g}", voxels=int(m.sum()))
            for k in ("D4", "S5", "S8"):
                d[k] = int((z[f"{k}_after"][m] > RIGHT).sum())
            rows.append(d)
        tabs[f"s8_wrong_by_level_{smp}"] = markdown_table(
            rows, ["sample", "level", "voxels", "D4", "S5", "S8"]
        )
    if stage_rows:
        tabs["s8_stages"] = markdown_table(
            stage_rows,
            [
                "sample",
                "run",
                "seeds",
                "reseeds",
                "grid",
                "seed_step",
                "grid_fallback",
                "neighbour",
                "revisit",
                "pass_",
                "total",
            ],
            formats={
                k: "d"
                for k in ("seed_step", "grid_fallback", "neighbour", "revisit", "pass_", "total")
            },
        )
        tabs["s8_wrong"] = markdown_table(
            wrong_rows, ["sample", "run", "before", "before_ci", "after", "after_ci"]
        )
        tabs["s8_paired"] = markdown_table(
            pair_rows,
            ["sample", "vs", "when", "wrong_S8", "wrong_base", "fixed", "broken", "both_wrong"],
        )
    if out:
        fig, axes = plt.subplots(1, len(out), figsize=(4.5 * len(out), 4), sharey=False)
        axes = np.atleast_1d(axes)
        cols = ["#1f78b4", "#ff7f00", "#33a02c", "#b2df8a", "#6a3d9a"]
        stages = ["seed_step", "grid_fallback", "neighbour", "revisit", "pass_"]
        labs = ["seed step", "grid fallbacks", "neighbour", "revisit", "boundary pass"]
        for ax, (smp, res) in zip(axes, out.items()):
            for i, k in enumerate(("D4", "S5", "S8")):
                bottom = 0.0
                for st, col, lab in zip(stages, cols, labs):
                    v = res[k]["stages"][st] / 1e6
                    ax.bar(i, v, bottom=bottom, color=col, label=lab if i == 0 else None)
                    bottom += v
                ax.text(i, bottom, f"{bottom:.2f}", ha="center", va="bottom", fontsize=8)
            ax.set_xticks(range(3))
            ax.set_xticklabels(["D4", "S5", "S8"])
            ax.set_title(f"{smp}: evaluations (millions), with the pass")
        axes[0].legend(fontsize=7)
        fig.tight_layout()
        fig.savefig(OUT / "ar3_reseed_s8.png", dpi=110)
        plt.close(fig)
    return {"s8": out}, tabs


def build() -> Tuple[Dict[str, Any], Dict[str, str]]:
    summ: Dict[str, Any] = {}
    tabs: Dict[str, str] = {}
    for fn in (classes, clusters, remedies, part3):
        s, t = fn()
        summ.update(s)
        tabs.update(t)
    return summ, tabs
