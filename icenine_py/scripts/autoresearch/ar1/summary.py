"""AR1 summary: summary.json, tables.md (doc tables) and the figure from the saved results.

  uv run python scripts/autoresearch/ar1/summary.py
Reads benchmarks/autoresearch/ar1/{sanity_*,recall_*,e2e_*,misses_pilot_base,sparse_pilot,
timing_pilot}.json; writes summary.json, tables.md and ar1_recall.png beside them.
"""

import json
import sys
from pathlib import Path
from typing import Any, Dict, List

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "common"))
import ar1_common as AC  # noqa: E402
import doc_tables  # noqa: E402
from stats import wilson  # noqa: E402

OUT = AC.OUT
SETS = [("e0", "E0 (200 voxels)"), ("pilot", "Pilot full searches (69 voxels)")]
KS = (1, 5, 20)
WRONG_DEG = 1.0


def frac(k: int, n: int) -> str:
    lo, hi = wilson(k, n)
    return f"{k}/{n} ({100 * k / n:.1f}%, {100 * lo:.1f}-{100 * hi:.1f})"


def load_e2e(kind: str) -> List[Dict[str, Any]]:
    rows: List[Dict[str, Any]] = []
    for f in sorted(OUT.glob(f"e2e_{kind}_base_part*.json")):
        rows += json.loads(f.read_text())["rows"]
    rows.sort(key=lambda r: r["voxel"])
    return rows


def arm(rows: List[Dict[str, Any]], n_fin: int, k_strict: int = 5) -> Dict[str, Any]:
    """Finisher on the best ``n_fin`` of the first ``k_strict`` candidates by strict cost; the
    lowest final cost wins."""
    err, ev = [], []
    for r in rows:
        f = r["finishers"][:n_fin]
        err.append(min(f, key=lambda x: x["cost"])["err"])
        ev.append(k_strict + sum(x["evals"] for x in f))
    return dict(err=np.array(err), evals=np.array(ev))


def main() -> None:
    S: Dict[str, Any] = {}
    tables: Dict[str, str] = {}
    # ---- sanity
    san_rows = []
    for kind, label in SETS:
        fj = sorted(OUT.glob(f"sanity_{kind}.json")) or sorted(OUT.glob(f"sanity_{kind}_p*.json"))
        if not fj:
            continue
        srows: List[Dict[str, Any]] = []
        for f in fj:
            srows += json.loads(f.read_text())
        zz = [np.load(f.with_suffix(".npz")) for f in fj]
        ang = np.concatenate([z["ang"] for z in zz])
        rel = np.concatenate([z["rel"] for z in zz])
        natt = sum(r["n_att"] for r in srows)
        tr = {t: sum(r["truth"][t] for r in srows if r["n_att"]) for t in ("0.25", "0.5", "1", "2")}
        co = {
            t: sum(r["control"][t] for r in srows if r["n_att"]) for t in ("0.25", "0.5", "1", "2")
        }
        S[f"sanity_{kind}"] = dict(
            n_voxels=len(srows),
            n_attached=natt,
            truth_within=tr,
            control_within=co,
            ang_p50=float(np.median(ang)),
            ang_p95=float(np.percentile(ang, 95)),
            ang_p99=float(np.percentile(ang, 99)),
            ang_max=float(ang.max()),
            len_p95=float(np.percentile(rel, 95)),
            len_max=float(rel.max()),
        )
        san_rows.append(
            dict(
                set=label,
                voxels=len(srows),
                attached=natt,
                t05=f"{tr['0.5']} ({100 * tr['0.5'] / natt:.1f}%)",
                c05=f"{co['0.5']} ({100 * co['0.5'] / natt:.2f}%)",
                t1=f"{tr['1']}",
                c1=f"{co['1']}",
                ang50=float(np.median(ang)),
                ang95=float(np.percentile(ang, 95)),
                angmax=float(ang.max()),
                len95=float(np.percentile(rel, 95)),
            )
        )
    tables["ar1_sanity"] = doc_tables.markdown_table(
        san_rows,
        [
            "set",
            "voxels",
            "attached",
            "t05",
            "c05",
            "t1",
            "c1",
            "ang50",
            "ang95",
            "angmax",
            "len95",
        ],
        {"ang50": ".3f", "ang95": ".3f", "angmax": ".3f", "len95": ".4f"},
    )
    # ---- recall
    rec_rows, thr_rows, dist_rows, time_rows = [], [], [], []
    fig_data: Dict[str, Any] = {}
    for kind, label in SETS:
        rows, _ = AC.load_recall(kind)
        if not rows:
            continue
        n = len(rows)
        D = [np.array(r["dist"]) for r in rows]

        def hit(K: int, t: float) -> int:
            return sum(1 for d in D if len(d) and d[:K].min() < t)

        rec_rows.append(dict(set=label, **{f"K{K}": frac(hit(K, 2.0), n) for K in KS}))
        for t in (0.5, 2.0, 5.0):
            thr_rows.append(
                dict(set=label, within=f"{t:g} deg", **{f"K{K}": f"{hit(K, t)}/{n}" for K in KS})
            )
        best20 = np.array([d[:20].min() if len(d) else np.inf for d in D])
        top1 = np.array([d[0] if len(d) else np.inf for d in D])
        ok20 = best20 < 2.0
        dist_rows.append(
            dict(
                set=label,
                best_med=float(np.median(best20[ok20])),
                best_p90=float(np.percentile(best20[ok20], 90)),
                best_max=float(best20[ok20].max()),
                top1_med=float(np.median(top1[top1 < 2.0])) if (top1 < 2.0).any() else None,
                n_cand_med=float(np.median([len(d) for d in D])),
            )
        )
        sec = np.array([r["seconds"].get("total", np.nan) for r in rows])
        time_rows.append(
            dict(
                set=label,
                med_ms=1e3 * float(np.nanmedian(sec)),
                p90_ms=1e3 * float(np.nanpercentile(sec, 90)),
                med_votes=float(np.median([r["n_votes"] for r in rows])),
                med_att=float(np.median([r["n_att"] for r in rows])),
                med_true=float(np.median([r["n_true_explained"] for r in rows])),
            )
        )
        fig_data[kind] = [
            sum(1 for d in D if len(d) and d[:K].min() < 2.0) / n for K in range(1, 21)
        ]
        S[f"recall_{kind}"] = dict(
            n=n,
            **{f"k{K}_{t:g}": hit(K, t) for K in KS for t in (0.5, 2.0, 5.0)},
            best_med=dist_rows[-1]["best_med"],
            best_max=dist_rows[-1]["best_max"],
            median_index_ms_contended=time_rows[-1]["med_ms"],
        )
    tables["ar1_recall"] = doc_tables.markdown_table(rec_rows, ["set", "K1", "K5", "K20"])
    tables["ar1_recall_thr"] = doc_tables.markdown_table(
        thr_rows, ["set", "within", "K1", "K5", "K20"]
    )
    tables["ar1_dist"] = doc_tables.markdown_table(
        dist_rows,
        ["set", "best_med", "best_p90", "best_max", "top1_med", "n_cand_med"],
        {
            "best_med": ".3f",
            "best_p90": ".3f",
            "best_max": ".3f",
            "top1_med": ".3f",
            "n_cand_med": ".0f",
        },
    )
    tables["ar1_index_time"] = doc_tables.markdown_table(
        time_rows,
        ["set", "med_ms", "p90_ms", "med_votes", "med_att", "med_true"],
        {"med_ms": ".0f", "p90_ms": ".0f", "med_votes": ".0f", "med_att": ".0f", "med_true": ".0f"},
    )
    # ---- misses (pilot)
    mp = OUT / "misses_pilot_base.json"
    if mp.exists():
        M = json.loads(mp.read_text())
        miss = [r for r in M if "reason" in r]
        reasons = {
            k: sum(r["reason"] == k for r in miss)
            for k in ("too_few_rays", "neighbour", "no_grain")
        }
        kinds = {
            k: float(np.median([r[k] for r in M]))
            for k in ("n_att", "n_own_pure", "n_own_mixed", "n_other", "n_true_explained")
        }
        gain = [r["rank1_explained"] - r["n_true_explained"] for r in miss]
        rank = {k: sum(1 for r in M if r["rank_of_truth"] == k) for k in (1, 2, 3)}
        S["misses_pilot"] = dict(
            n=len(M),
            rank1_misses=len(miss),
            reasons=reasons,
            medians=kinds,
            truth_rank=rank,
            true_explained_mixed_total=int(sum(r["true_explained_mixed"] for r in M)),
            true_explained_other_total=int(sum(r["true_explained_other"] for r in M)),
            true_explained_total=int(sum(r["n_true_explained"] for r in M)),
            rank1_extra_rays_median=float(np.median(gain)) if gain else None,
            rank1_extra_rays_min=int(min(gain)) if gain else None,
        )
        rows_m = [
            dict(what="rank-1 candidate is not the truth", n=f"{len(miss)}/{len(M)}"),
            dict(what="  reason: too few rays (truth explains < 8)", n=reasons["too_few_rays"]),
            dict(what="  reason: another real grain explains more rays", n=reasons["neighbour"]),
            dict(what="  reason: matches no grain (coincidence of rays)", n=reasons["no_grain"]),
            dict(
                what="truth is candidate 1 / 2 / 3 (of 20)", n=f"{rank[1]} / {rank[2]} / {rank[3]}"
            ),
        ]
        tables["ar1_misses"] = doc_tables.markdown_table(rows_m, ["what", "n"])
        tables["ar1_ray_kinds"] = doc_tables.markdown_table(
            [
                dict(
                    quantity="attached rays",
                    median=kinds["n_att"],
                ),
                dict(quantity="drawn only by the voxel's own grain", median=kinds["n_own_pure"]),
                dict(
                    quantity="drawn by the own grain and another grain", median=kinds["n_own_mixed"]
                ),
                dict(quantity="drawn by other grains only", median=kinds["n_other"]),
                dict(
                    quantity="explained by the true orientation", median=kinds["n_true_explained"]
                ),
            ],
            ["quantity", "median"],
            {"median": ".0f"},
        )
    # ---- end to end
    e2e_rows = []
    fig_ev: Dict[str, Any] = {}
    for kind, label in SETS:
        rows = load_e2e(kind)
        if not rows:
            continue
        n = len(rows)
        stored_err = np.array([r["stored_err"] for r in rows])
        stored_ev = np.array([r["stored_evals"] for r in rows])
        sw = stored_err > WRONG_DEG

        def add(name: str, err: np.ndarray, ev: np.ndarray) -> None:
            w = err > WRONG_DEG
            e2e_rows.append(
                dict(
                    set=label,
                    arm=name,
                    wrong=frac(int(w.sum()), n),
                    fixed=int((sw & ~w).sum()),
                    broken=int((~sw & w).sum()),
                    both=int((sw & w).sum()),
                    mean_ev=float(np.mean(ev)),
                    med_ev=float(np.median(ev)),
                    med_err=float(np.median(err)),
                )
            )
            S.setdefault(f"e2e_{kind}", {})[name] = dict(
                n=n,
                wrong=int(w.sum()),
                fixed=int((sw & ~w).sum()),
                broken=int((~sw & w).sum()),
                both_wrong=int((sw & w).sum()),
                mean_evals=float(np.mean(ev)),
                median_evals=float(np.median(ev)),
                median_err_deg=float(np.median(err)),
            )

        add("stored full search", stored_err, stored_ev)
        for K in (1, 5, 20):
            err = np.array([r["err_pick"][str(K)] for r in rows])
            add(f"strict cost on first {K}, no finisher", err, np.full(n, float(K)))
        a1, a3 = arm(rows, 1), arm(rows, 3)
        add("AR1: strict on 5, finisher on best 1", a1["err"], a1["evals"])
        add("AR1: strict on 5, finisher on best 3", a3["err"], a3["evals"])
        fig_ev[kind] = dict(stored=float(np.median(stored_ev)), ar1=float(np.median(a1["evals"])))
        S.setdefault(f"e2e_{kind}", {})["stored_wrong_voxels"] = [
            int(r["voxel"]) for r, w in zip(rows, sw) if w
        ]
    tables["ar1_e2e"] = doc_tables.markdown_table(
        e2e_rows,
        ["set", "arm", "wrong", "fixed", "broken", "both", "mean_ev", "med_ev", "med_err"],
        {"mean_ev": ",.0f", "med_ev": ",.0f", "med_err": ".3f"},
    )
    # ---- sparse
    sp = OUT / "sparse_pilot.json"
    if sp.exists():
        d = json.loads(sp.read_text())
        S["sparse_pilot"] = d
        tables["ar1_sparse"] = doc_tables.markdown_table(
            [
                dict(
                    case=k,
                    rays=v["median_true_rays"],
                    K1=frac(v["k1"], v["n"]),
                    K5=frac(v["k5"], v["n"]),
                    K20=frac(v["k20"], v["n"]),
                )
                for k, v in d.items()
            ],
            ["case", "rays", "K1", "K5", "K20"],
            {"rays": ".0f"},
        )
    # ---- timing
    tp = OUT / "timing_pilot.json"
    if tp.exists():
        t = json.loads(tp.read_text())
        t.pop("rows", None)
        S["timing_pilot"] = t
        tables["ar1_timing"] = doc_tables.markdown_table(
            [
                dict(
                    quantity=f"one-time pairing of all spots ({t['n_rays_layer']:,} rays)",
                    seconds=t["observed_rays_build_s"],
                ),
                dict(quantity="attach rays to a voxel (median)", seconds=t["attach_s_median"]),
                dict(quantity="index a voxel (median)", seconds=t["index_s_median"]),
                dict(quantity="index a voxel (90th percentile)", seconds=t["index_s_p90"]),
                dict(quantity="  of which voting (median)", seconds=t["vote_s_median"]),
                dict(quantity="  of which clustering (median)", seconds=t["cluster_s_median"]),
                dict(quantity="  of which refining (median)", seconds=t["refine_s_median"]),
            ],
            ["quantity", "seconds"],
            {"seconds": ".3f"},
        )
    # ---- files
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "summary.json").write_text(json.dumps(S, indent=1, default=float) + "\n")
    doc_tables.write_tables(OUT / "tables.md", tables)
    # ---- figure
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        fig, ax = plt.subplots(1, 2, figsize=(10, 3.8))
        for kind, label in SETS:
            if kind in fig_data:
                ax[0].plot(range(1, 21), 100 * np.array(fig_data[kind]), "o-", ms=4, label=label)
        ax[0].set_xlabel("K (candidates kept, ranked by rays explained)")
        ax[0].set_ylabel("voxels with the truth within 2 deg (%)")
        ax[0].set_ylim(0, 102)
        ax[0].set_xticks([1, 5, 10, 15, 20])
        ax[0].legend(loc="lower right", fontsize=8)
        ax[0].set_title("Recall at K")
        names, vals = [], []
        for kind, label in SETS:
            if kind in fig_ev:
                names += [f"{kind}\nstored", f"{kind}\nAR1"]
                vals += [fig_ev[kind]["stored"], fig_ev[kind]["ar1"]]
        if names:
            bars = ax[1].bar(range(len(vals)), vals, color=["0.6", "C0"] * (len(vals) // 2))
            ax[1].set_yscale("log")
            ax[1].set_xticks(range(len(vals)))
            ax[1].set_xticklabels(names, fontsize=8)
            ax[1].set_ylabel("evaluations per voxel (median)")
            for b, v in zip(bars, vals):
                ax[1].text(
                    b.get_x() + b.get_width() / 2, v * 1.1, f"{v:,.0f}", ha="center", fontsize=8
                )
            ax[1].set_title("Evaluations of one full search")
        fig.tight_layout()
        fig.savefig(OUT / "ar1_recall.png", dpi=130)
    except ImportError:
        print("matplotlib missing; no figure")
    print(
        json.dumps(
            {k: v for k, v in S.items() if k not in ("sparse_pilot",)}, indent=1, default=float
        )[:3000]
    )


if __name__ == "__main__":
    main()
