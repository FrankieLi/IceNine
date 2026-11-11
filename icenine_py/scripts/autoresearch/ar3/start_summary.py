"""AR3 step 9: the summary (reads only committed files under ``start/`` and ``reseed/slim/``).

  uv run python scripts/autoresearch/ar3/start_summary.py

Writes ``start/summary.json``, ``start/tables.md`` (marker blocks for ``sync_doc_tables.py``) and
the plots ``start/ar3_start_evals.png`` (evaluations by stage with the pass: D4, S5, S8, S9 on the
three samples) and ``start/ar3_start_share.png`` (the share of neighbour fits whose start passes,
by gradient level).

Evaluations are raw counts of strict-cost calls (the stages of step 8's page: seed step, grid
fallbacks, neighbour fits, revisits, boundary pass). In S9 the start tests of F1 and the extra
fits of F2 are inside the neighbour and revisit stages (their counts are given beside). Intervals
are Wilson 95% intervals on voxels, fits or seed steps taken as independent, although they
cluster in grains. Errors are in degrees; "right" is within 1 degree of the truth.
"""

import json
from pathlib import Path
from typing import Any, Dict, List, Optional

import ar3_common  # noqa: F401  (sets sys.path for doc_tables and stats)
import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from doc_tables import markdown_table, write_tables  # noqa: E402
from stats import wilson  # noqa: E402

HERE = Path(__file__).resolve().parent
ICE = HERE.parents[2]
START = ICE / "benchmarks" / "autoresearch" / "ar3" / "start"
SLIM8 = ICE / "benchmarks" / "autoresearch" / "ar3" / "reseed" / "slim"
RIGHT = 1.0
SAMPLES = ("big", "grad", "mid")
STAGES = ["seed_step", "grid_fallback", "neighbour", "revisit", "pass_"]


def rj(p: Path) -> Any:
    return json.loads(p.read_text())


def ci(k: int, n: int, d: int = 1) -> str:
    lo, hi = wilson(k, n)
    return f"{100 * lo:.{d}f} to {100 * hi:.{d}f}"


def pc(a: np.ndarray, q: float) -> float:
    return float(np.percentile(a, q))


def stage_row(r: Dict[str, Any]) -> Dict[str, Any]:
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
    assert tot - pe == int(r["evals_total"]), (tot - pe, r["evals_total"])
    return dict(
        stages=st,
        total=tot,
        seeds=int(s["n_seeds"]),
        reseeds=int(s["n_reseeds"]),
        grid=int(s.get("index_fallback_seed_test", 0) + s.get("index_fallback_no_candidate", 0)),
        neighbour_fits=int(s["n_neighbor_fits"]),
        revisits=int(s["n_revisit_attempted"]),
    )


def load_runs() -> Dict[str, Dict[str, Any]]:
    """{sample: {run: run info with pass}} from step 8's slim files and this task's."""
    out: Dict[str, Dict[str, Any]] = {}
    for smp in ("big", "grad"):
        s8 = rj(SLIM8 / f"s8_{smp}.json")
        base = rj(SLIM8 / f"baselines_{smp}.json")
        out[smp] = {"D4": base["D4"], "S5": base["S5"], "S8": s8}
    for smp in SAMPLES:
        p = START / f"runs_{smp}.json"
        if p.exists():
            out.setdefault(smp, {}).update(rj(p))
    return out


def load_errors(smp: str) -> Dict[str, np.ndarray]:
    z = dict(np.load(SLIM8 / f"errors_{smp}.npz")) if smp != "mid" else {}
    p = START / f"errors_{smp}.npz"
    if p.exists():
        z.update(dict(np.load(p)))
    return z


# ------------------------------------------------------------------------------------- part 1
def part1() -> Any:
    out: Dict[str, Any] = {}
    tabs: Dict[str, str] = {}
    rows: List[Dict[str, Any]] = []
    mrows: List[Dict[str, Any]] = []
    cls_rows: List[Dict[str, Any]] = []
    for smp in ("big", "grad"):
        p = START / f"part1_{smp}.json"
        if not p.exists():
            continue
        d = rj(p)
        out[smp] = d
        for ph in ("neighbour", "revisit"):
            keys = [k for k in d["by_level"] if k.startswith(ph + "/")]
            for k in keys:
                L = k.split("/")[1]
                r = d["by_level"][k]
                n = r["fits"]
                if n == 0:
                    continue
                rows.append(
                    dict(
                        sample=smp,
                        phase=ph,
                        level=L,
                        fits=n,
                        evals=r["evals"],
                        s90=f"{r['start_ge_0.90']} ({100 * r['start_ge_0.90'] / n:.1f}%)",
                        s95=f"{r['start_ge_0.95']} ({100 * r['start_ge_0.95'] / n:.1f}%)",
                        s98=f"{r['start_ge_0.98']} ({100 * r['start_ge_0.98'] / n:.1f}%)",
                        ci95=ci(r["start_ge_0.95"], n),
                    )
                )

                def mv(t: str) -> str:
                    m50, m90 = r[f"move_median_{t}"], r[f"move_p90_{t}"]
                    return "n/a" if m50 != m50 else f"{m50:.3f} / {m90:.3f}"

                mrows.append(
                    dict(
                        sample=smp,
                        phase=ph,
                        level=L,
                        m90=mv("0.90"),
                        m95=mv("0.95"),
                        m98=mv("0.98"),
                        w90=r["wrong_start_0.90"],
                        w95=r["wrong_start_0.95"],
                        w98=r["wrong_start_0.98"],
                    )
                )
        if "class_a" in d:
            c = d["class_a"]
            cls_rows.append(
                dict(
                    sample=smp,
                    voxels=c["voxels_in_clusters_first_appearance"],
                    class_a=c["class_a"],
                    truth_pass=c["truth_passes_relative_test"],
                    truth_fail=c["truth_fails_relative_test"],
                    after_pass=c["after_extra_fit_passes"],
                    after_ci=ci(c["after_extra_fit_passes"], max(c["class_a"], 1)),
                    hit_stored=c["median_hit_stored"],
                    hit_truth=c["median_hit_truth"],
                    best=c["median_best"],
                    hit_after=c["median_hit_after"],
                    err_end=c["median_err_end"],
                    err_after=c["median_err_after"],
                    extra_evals=c["evals_extra_total"],
                )
            )
    if rows:
        tabs["p1_start_share"] = markdown_table(
            rows,
            [
                "sample",
                "phase",
                "level",
                "fits",
                "evals",
                "s90",
                "s95",
                "s98",
                "ci95",
            ],
            formats={"evals": "d"},
        )
        tabs["p1_move"] = markdown_table(
            mrows,
            ["sample", "phase", "level", "m90", "m95", "m98", "w90", "w95", "w98"],
        )
    if cls_rows:
        tabs["p1_class_a"] = markdown_table(
            cls_rows,
            [
                "sample",
                "voxels",
                "class_a",
                "truth_pass",
                "truth_fail",
                "after_pass",
                "after_ci",
                "hit_stored",
                "hit_truth",
                "best",
                "hit_after",
                "err_end",
                "err_after",
                "extra_evals",
            ],
            formats={k: ".3f" for k in ("hit_stored", "hit_truth", "best", "hit_after", "err_end")}
            | {"err_after": ".3f", "extra_evals": "d"},
        )
    # start and end of the neighbour fits and revisits (medians), from the committed features
    se_rows: List[Dict[str, Any]] = []
    for smp in ("big", "grad"):
        p = START / f"feat_{smp}.npz"
        if not p.exists():
            continue
        f = dict(np.load(p))
        for ph, name in ((0, "neighbour"), (1, "revisit")):
            for L in sorted(set(f["level"].tolist())):
                sel = (f["phase"] == ph) & (f["level"] == L)
                if not sel.any():
                    continue
                se_rows.append(
                    dict(
                        sample=smp,
                        phase=name,
                        level=f"{L:g}",
                        fits=int(sel.sum()),
                        hit_start=float(np.median(f["start_hit"][sel])),
                        hit_end=float(np.median(f["hit_end"][sel])),
                        err_start=float(np.median(f["err_start"][sel])),
                        err_end=float(np.median(f["err_end"][sel])),
                        evals_per_fit=float(f["evals"][sel].mean()),
                    )
                )
    if se_rows:
        tabs["p1_start_end"] = markdown_table(
            se_rows,
            ["sample", "phase", "level", "fits", "hit_start", "hit_end", "err_start"]
            + ["err_end", "evals_per_fit"],
            formats={k: ".3f" for k in ("hit_start", "hit_end", "err_start", "err_end")}
            | {"evals_per_fit": ".1f"},
        )
    # predictions 1 and 2
    pred: Dict[str, Any] = {}
    if "big" in out:
        bl = out["big"]["by_level"]
        low = [bl[f"neighbour/{L}"] for L in ("0.0", "0.2", "0.5", "1.0") if f"neighbour/{L}" in bl]
        n_low = sum(r["fits"] for r in low)
        k_low = sum(r["start_ge_0.95"] for r in low)
        hi = bl["neighbour/4.0"]
        pred["p2_big_le1deg"] = dict(
            n=n_low, k=k_low, share=k_low / n_low, ci=ci(k_low, n_low), prediction=">0.5"
        )
        pred["p2_big_4deg"] = dict(
            n=hi["fits"],
            k=hi["start_ge_0.95"],
            share=hi["start_ge_0.95"] / hi["fits"],
            ci=ci(hi["start_ge_0.95"], hi["fits"]),
            prediction="<0.25",
        )
    if cls_rows:
        c = out["big"]["class_a"]
        pred["p1"] = dict(
            class_a=c["class_a"],
            truth_passes=c["truth_passes_relative_test"],
            share=c["truth_passes_relative_test"] / c["class_a"],
            prediction="most",
        )
    return out, tabs, pred


# ------------------------------------------------------------------------------------- part 2
def part2() -> Any:
    out: Dict[str, Any] = {}
    tabs: Dict[str, str] = {}
    pk_p, vd_p = START / "f_pick.json", START / "part2_verdict.json"
    if not pk_p.exists():
        return out, tabs
    pk = rj(pk_p)
    out["pick"] = pk
    rows = []
    for t, r in ((k, v) for k, v in pk["F1"]["pooled"].items() if isinstance(v, dict)):
        rows.append(
            dict(
                set="pick",
                remedy="F1",
                setting=f"hit >= {t}",
                kept=r["kept"],
                evals=r["evals_f1"],
                cut=100 * r["cut"],
                fixed=r["fixed"],
                broken=r["broken"],
                both_wrong=r["both_wrong"],
            )
        )
    if "pooled" in pk.get("F2", {}):
        for s, r in ((k, v) for k, v in pk["F2"]["pooled"].items() if isinstance(v, dict)):
            rows.append(
                dict(
                    set="pick",
                    remedy="F2",
                    setting=f"{s} deg",
                    kept=r["clusters_vanish"],
                    evals=r["f2_evals"],
                    cut=100 * r["cut"],
                    fixed=r["fixed"],
                    broken=r["broken"],
                    both_wrong=r["both_wrong"],
                )
            )
    if vd_p.exists():
        vd = rj(vd_p)
        out["verdict"] = vd
        for t, r in ((k, v) for k, v in vd["F1"]["pooled"].items() if isinstance(v, dict)):
            rows.append(
                dict(
                    set="test",
                    remedy="F1",
                    setting=f"hit >= {t}",
                    kept=r["kept"],
                    evals=r["evals_f1"],
                    cut=100 * r["cut"],
                    fixed=r["fixed"],
                    broken=r["broken"],
                    both_wrong=r["both_wrong"],
                )
            )
        if "pooled" in vd.get("F2", {}):
            for s, r in ((k, v) for k, v in vd["F2"]["pooled"].items() if isinstance(v, dict)):
                rows.append(
                    dict(
                        set="test",
                        remedy="F2",
                        setting=f"{s} deg",
                        kept=r["clusters_vanish"],
                        evals=r["f2_evals"],
                        cut=100 * r["cut"],
                        fixed=r["fixed"],
                        broken=r["broken"],
                        both_wrong=r["both_wrong"],
                    )
                )
    tabs["p2_remedies"] = markdown_table(
        rows,
        ["set", "remedy", "setting", "kept", "evals", "cut", "fixed", "broken", "both_wrong"],
        formats={"cut": ".1f", "evals": "d"},
    )
    return out, tabs


# ------------------------------------------------------------------------------------- part 3
def part3() -> Any:
    runs = load_runs()
    out: Dict[str, Any] = {}
    tabs: Dict[str, str] = {}
    stage_rows: List[Dict[str, Any]] = []
    wrong_rows: List[Dict[str, Any]] = []
    err_rows: List[Dict[str, Any]] = []
    pair_rows: List[Dict[str, Any]] = []
    s9_rows: List[Dict[str, Any]] = []
    for smp in SAMPLES:
        if smp not in runs:
            continue
        z = load_errors(smp)
        res: Dict[str, Any] = {}
        for k in ("D4", "S5", "S8", "S9"):
            if k not in runs[smp] or "pass" not in runs[smp][k]:
                continue
            sr = stage_row(runs[smp][k])
            before = z[f"{k}_before"]
            after = z[f"{k}_after"]
            n = len(after)
            res[k] = dict(
                **sr,
                wrong_before=int((before > RIGHT).sum()),
                wrong_after=int((after > RIGHT).sum()),
                n=n,
                err_after_median=pc(after, 50),
                err_after_p90=pc(after, 90),
                err_after_p99=pc(after, 99),
                err_before_median=pc(before, 50),
            )
            stage_rows.append(
                dict(sample=smp, run=k + " + pass", seeds=sr["seeds"], reseeds=sr["reseeds"])
                | dict(grid=sr["grid"], **sr["stages"], total=sr["total"])
            )
            wb, wa = res[k]["wrong_before"], res[k]["wrong_after"]
            wrong_rows.append(
                dict(
                    sample=smp,
                    run=k,
                    before=wb,
                    before_ci=ci(wb, n, 2),
                    after=wa,
                    after_ci=ci(wa, n, 2),
                )
            )
            err_rows.append(
                dict(
                    sample=smp,
                    run=k,
                    med=res[k]["err_after_median"],
                    p90=res[k]["err_after_p90"],
                    p99=res[k]["err_after_p99"],
                )
            )
            if k == "S9":
                c = runs[smp][k].get("s9", {}).get("counters", {})
                s9_rows.append(
                    dict(
                        sample=smp,
                        f1=runs[smp][k]["s9"]["settings"]["f1"],
                        f2=runs[smp][k]["s9"]["settings"]["f2"],
                        f1_tests=c.get("f1_tests"),
                        f1_kept=c.get("f1_kept"),
                        f2_runs=c.get("f2_runs"),
                        f2_evals=c.get("f2_evals"),
                        f2_replaced=c.get("f2_replaced"),
                        f2_would_pass=c.get("f2_would_pass"),
                    )
                )
        pairs = [("S8", "S5"), ("S8", "D4"), ("S5", "D4"), ("S9", "S8"), ("S9", "S5"), ("S9", "D4")]
        for a, b in pairs:
            if a not in res or b not in res:
                continue
            for when in ("before", "after"):
                wa = z[f"{a}_{when}"] > RIGHT
                wb = z[f"{b}_{when}"] > RIGHT
                pair_rows.append(
                    dict(
                        sample=smp,
                        run=a,
                        vs=b,
                        when=when,
                        wrong_run=int(wa.sum()),
                        wrong_vs=int(wb.sum()),
                        fixed=int((wb & ~wa).sum()),
                        broken=int((wa & ~wb).sum()),
                        both_wrong=int((wa & wb).sum()),
                    )
                )
        out[smp] = res
    tabs["p3_stages"] = markdown_table(
        stage_rows,
        ["sample", "run", "seeds", "reseeds", "grid"] + STAGES + ["total"],
        formats={k: "d" for k in STAGES + ["total"]},
    )
    tabs["p3_wrong"] = markdown_table(
        wrong_rows, ["sample", "run", "before", "before_ci", "after", "after_ci"]
    )
    tabs["p3_error"] = markdown_table(
        err_rows,
        ["sample", "run", "med", "p90", "p99"],
        formats={k: ".4f" for k in "med p90 p99".split()},
    )
    tabs["p3_paired"] = markdown_table(
        pair_rows,
        ["sample", "run", "vs", "when", "wrong_run", "wrong_vs", "fixed", "broken", "both_wrong"],
    )
    if s9_rows:
        tabs["p3_s9"] = markdown_table(
            s9_rows,
            ["sample", "f1", "f2", "f1_tests", "f1_kept", "f2_runs", "f2_evals", "f2_replaced"],
        )
    # the effect of F1 by stage (S9 against S8) and by gradient level
    eff_rows: List[Dict[str, Any]] = []
    lvl_rows: List[Dict[str, Any]] = []
    for smp, res in out.items():
        if "S8" not in res or "S9" not in res:
            continue
        a, b = res["S8"], res["S9"]
        nr8 = a["stages"]["neighbour"] + a["stages"]["revisit"]
        nr9 = b["stages"]["neighbour"] + b["stages"]["revisit"]
        eff_rows.append(
            dict(
                sample=smp,
                nr_s8=nr8,
                nr_s9=nr9,
                nr_change=100.0 * (nr9 / nr8 - 1),
                total_s8=a["total"],
                total_s9=b["total"],
                total_change=100.0 * (b["total"] / a["total"] - 1),
                reseeds_s8=a["reseeds"],
                reseeds_s9=b["reseeds"],
            )
        )
        z = load_errors(smp)
        for L in sorted(set(z["level"].tolist())):
            m = z["level"] == L
            lvl_rows.append(
                dict(
                    sample=smp,
                    level=f"{L:g}",
                    voxels=int(m.sum()),
                    wrong_s8=int((z["S8_after"][m] > RIGHT).sum()),
                    wrong_s9=int((z["S9_after"][m] > RIGHT).sum()),
                    med_s8=pc(z["S8_after"][m], 50),
                    med_s9=pc(z["S9_after"][m], 50),
                    p99_s8=pc(z["S8_after"][m], 99),
                    p99_s9=pc(z["S9_after"][m], 99),
                )
            )
    if eff_rows:
        tabs["p3_f1_effect"] = markdown_table(
            eff_rows,
            [
                "sample",
                "nr_s8",
                "nr_s9",
                "nr_change",
                "total_s8",
                "total_s9",
                "total_change",
                "reseeds_s8",
                "reseeds_s9",
            ],
            formats={"nr_s8": "d", "nr_s9": "d", "total_s8": "d", "total_s9": "d"}
            | {"nr_change": ".1f", "total_change": ".1f"},
        )
        tabs["p3_by_level"] = markdown_table(
            lvl_rows,
            ["sample", "level", "voxels", "wrong_s8", "wrong_s9"]
            + ["med_s8", "med_s9", "p99_s8", "p99_s9"],
            formats={k: ".4f" for k in ("med_s8", "med_s9", "p99_s8", "p99_s9")},
        )
    # the plot
    fig, axes = plt.subplots(1, len(out), figsize=(4.6 * len(out), 4.2), sharey=False)
    axes = np.atleast_1d(axes)
    cols = ["#1f78b4", "#ff7f00", "#33a02c", "#b2df8a", "#6a3d9a"]
    labs = ["seed step", "grid fallbacks", "neighbour", "revisit", "boundary pass"]
    for ax, (smp, res) in zip(axes, out.items()):
        names = [k for k in ("D4", "S5", "S8", "S9") if k in res]
        for i, k in enumerate(names):
            bottom = 0.0
            for st, col, lab in zip(STAGES, cols, labs):
                v = res[k]["stages"][st] / 1e6
                ax.bar(i, v, bottom=bottom, color=col, label=lab if i == 0 else None)
                bottom += v
            ax.text(i, bottom, f"{bottom:.2f}", ha="center", va="bottom", fontsize=8)
        ax.set_xticks(range(len(names)))
        ax.set_xticklabels(names)
        ax.set_title(f"{smp}: evaluations (millions), with the pass")
    axes[0].legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(START / "ar3_start_evals.png", dpi=110)
    plt.close(fig)
    return out, tabs


def plot_share(p1: Dict[str, Any]) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(10, 3.8), sharey=True)
    for ax, smp in zip(axes, ("big", "grad")):
        if smp not in p1:
            ax.set_visible(False)
            continue
        bl = p1[smp]["by_level"]
        levels = [
            k.split("/")[1] for k in bl if k.startswith("neighbour/") and k != "neighbour/all"
        ]
        x = np.arange(len(levels))
        for j, (t, col) in enumerate(
            zip(("0.90", "0.95", "0.98"), ("#a6cee3", "#1f78b4", "#08306b"))
        ):
            y = [
                100 * bl[f"neighbour/{L}"][f"start_ge_{t}"] / max(bl[f"neighbour/{L}"]["fits"], 1)
                for L in levels
            ]
            ax.bar(x + (j - 1) * 0.27, y, 0.27, color=col, label=f"start hit ratio >= {t}")
        ax.set_xticks(x)
        ax.set_xticklabels(
            [f"{float(L):g}\n(n={bl[f'neighbour/{L}']['fits']})" for L in levels], fontsize=8
        )
        ax.set_xlabel("gradient level of the true grain (degrees)")
        ax.set_title(f"{smp}: neighbour fits whose start already passes")
    axes[0].set_ylabel("share of neighbour fits (%)")
    axes[0].legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(START / "ar3_start_share.png", dpi=110)
    plt.close(fig)


def main() -> None:
    p1, t1, pred = part1()
    p2, t2 = part2()
    p3, t3 = part3()
    plot_share(p1)
    summ = dict(part1=p1, part2=p2, part3=p3, predictions=pred)
    (START / "summary.json").write_text(json.dumps(summ, indent=1, default=float) + "\n")
    write_tables(START / "tables.md", {**t1, **t2, **t3})
    print(json.dumps(pred, indent=1, default=float))


if __name__ == "__main__":
    main()
