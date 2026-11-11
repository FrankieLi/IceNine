"""AR3 step 7: the summary of the four runs on ``big`` against the gradient sample.

  uv run python scripts/autoresearch/ar3/big_summary.py

Reads (all committed): ``big/errors.npz`` (per-voxel errors of each map before and after the pass),
``big/{gradient,fine,seed}/runs/*/run.json`` and the pass folders' ``pass.json``, ``big/gradient/
sample/sample.json``; and, for the gradient sample, the committed summaries of steps 4, 5 and 6
(nothing is rerun there). Writes ``big/summary.json``, ``big/tables.md`` and the two plots.

Stages of the evaluations (strict-equivalent units, as step 6): *seed step* (all index runs),
*grid fallbacks* (full grid searches after a failed seed step), *neighbour* (the BFS's fits from
a neighbour; for the paint run the fits of the painting loop), *revisit* (the BFS's revisit fits;
for the paint run the R4 fits), *scoring* (the paint run only), *pass* (step 5's pass, rule hit).
"""

import json
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import ar3_common as C
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from doc_tables import markdown_table, write_tables  # noqa: E402
from stats import wilson  # noqa: E402

BIG = C.ICE / "benchmarks" / "autoresearch" / "ar3" / "big"
GRAD = C.ICE / "benchmarks" / "autoresearch" / "ar3"
WRONG = 1.0
PIPES = ("D", "D4", "S5", "P4R4")
NAMES = {
    "D": "D (stock BFS, grid search)",
    "D4": "D4 (index seeds, CMA fit)",
    "S5": "S5 (D4 with step 6's seed step)",
    "P4R4": "P4R4 (paint, CMA fit)",
}
LEVELS = (0.0, 0.2, 0.5, 1.0, 2.0, 4.0)
# matched levels: big level -> gradient level with about the same voxel-to-voxel change
MATCH = {0.0: 0.0, 0.5: 0.1, 1.0: 0.2, 2.0: 0.5, 4.0: 1.0}
STAGES = ("seed_step", "grid_fallback", "neighbour", "revisit", "scoring", "pass_")
N = 24570


def jl(p: Path) -> Any:
    return json.loads(p.read_text())


def lv(x: float) -> str:
    return f"{x:g}"


def ci(k: int, n: int) -> str:
    lo, hi = wilson(k, n)
    d = 2 if hi < 0.1 else 1
    return f"{k} of {n} ({100 * lo:.{d}f} to {100 * hi:.{d}f}%)"


def ci_d(k: int, n: int) -> Dict[str, Any]:
    lo, hi = wilson(k, n)
    return dict(k=int(k), n=int(n), lo=lo, hi=hi)


# ---------------------------------------------------------------------------------- big side
def big_stats(pipe: str) -> Dict[str, Any]:
    """Counters and stages of one run on ``big`` (run.json and pass.json)."""
    r = jl(BIG / "slim" / f"{pipe}_run.json")
    pj = jl(BIG / "passes" / f"{pipe}_pass.json")
    return stages_of(pipe, r, pj["evals"], pj)


def make_slim() -> None:
    """Copy each run's ``run.json`` without its non-scalar stats (index logs and the like) to
    ``big/slim/<pipe>_run.json`` and the pass folder's ``pass.json`` to ``big/passes/``. Run on the
    full results; the summary reads the slim copies only."""
    src = {
        "D": BIG / "gradient" / "runs" / "D" / "run.json",
        "D4": BIG / "fine" / "runs" / "D4" / "run.json",
        "S5": BIG / "seed" / "runs" / "S5" / "run.json",
        "P4R4": BIG / "fine" / "runs" / "P4" / "run.json",
    }
    (BIG / "slim").mkdir(exist_ok=True)
    (BIG / "passes").mkdir(exist_ok=True)
    for pipe, p in src.items():
        r = jl(p)
        if "stats" in r:
            r["stats"] = {k: v for k, v in r["stats"].items() if isinstance(v, (int, float))}
        (BIG / "slim" / f"{pipe}_run.json").write_text(json.dumps(r, indent=1) + "\n")
        pj = BIG / "boundary" / "cache" / "passes" / f"D4_{pipe}" / "pass.json"
        (BIG / "passes" / f"{pipe}_pass.json").write_text(pj.read_text())
    import shutil

    shutil.copyfile(
        BIG / "seed" / "runs" / "S5" / "seed_log.json.gz", BIG / "slim" / "S5_seed_log.json.gz"
    )


def stages_of(pipe: str, r: Dict[str, Any], pass_evals: int, pj: Dict[str, Any]) -> Dict[str, Any]:
    if pipe == "P4R4":
        e = r["evals"]
        st = dict(
            seed_step=e["seed_index"] + e["index_rejected"],
            grid_fallback=e["seed_grid"] + e["grid_rejected"],
            neighbour=e["fit_accepted"] + e["fit_rejected"],
            revisit=r["evals_r4"],
            scoring=e["score"],
            pass_=pass_evals,
        )
        assert sum(st.values()) - pass_evals == r["evals_total"] + r["evals_r4"]
        out = dict(
            seeds=r["counts"]["seeds"],
            reseeds=None,
            grid_fallbacks=None,
            wall_run_s=r["wall_p_s"] + r["wall_r4_s"],
        )
    else:
        s = r["stats"]
        st = dict(
            seed_step=int(s["index_evals"]),
            grid_fallback=int(s.get("index_fallback_grid_evals", 0)),
            neighbour=int(s["n_evals_neighbor"]),
            revisit=int(s["n_evals_revisit"]),
            scoring=0,
            pass_=pass_evals,
        )
        assert (
            st["seed_step"] + st["grid_fallback"] == s["n_evals_seed"] + s["n_evals_reseed"]
        ), "the seed step and grid fallbacks are the seed and reseed buckets"
        out = dict(
            seeds=int(s["n_seeds"]),
            reseeds=int(s["n_reseeds"]),
            grid_fallbacks=int(
                s.get("index_fallback_seed_test", 0) + s.get("index_fallback_no_candidate", 0)
            ),
            wall_run_s=float(r["wall_total_s"]),
        )
        assert sum(st.values()) - pass_evals == r["evals_total"]
    out.update(stages=st, total=int(sum(st.values())), wall_pass_s=pj.get("wall_tries_s"))
    return out


# ------------------------------------------------------------------------------ gradient side
def grad_stats() -> Dict[str, Dict[str, Any]]:
    sp = jl(GRAD / "seed" / "summary_part2.json")
    bs = jl(GRAD / "boundary" / "summary.json")
    out: Dict[str, Dict[str, Any]] = {}
    for pipe, key in (("D", "D + pass"), ("D4", "D4 + pass"), ("S5", "S5 + pass")):
        row = next(x for x in sp["stages"] if x["run"] == key)
        st = dict(
            seed_step=row["seed_step"],
            grid_fallback=row["grid_fallback"],
            neighbour=row["neighbour"],
            revisit=row["revisit"],
            scoring=0,
            pass_=row["pass_"],
        )
        assert sum(st.values()) == row["total"]
        out[pipe] = dict(
            seeds=row["seeds"],
            reseeds=row["reseeds"],
            grid_fallbacks=row["grid_fallbacks"],
            wall_run_s=60.0 * row["wall_min"],
            stages=st,
            total=int(row["total"]),
        )
    p4 = jl(GRAD / "fine" / "runs" / "P4" / "run.json")
    pe = bs["maps"]["P4R4"]["evals_pass"]
    out["P4R4"] = stages_of("P4R4", p4, pe, dict(wall_tries_s=bs["maps"]["P4R4"]["wall_tries_s"]))
    assert out["P4R4"]["total"] == 54563676  # boundary/summary.json pipeline row
    return out


def grad_wrong(pipe: str, when: str, level: str) -> Tuple[int, int]:
    """(wrong, n) on the gradient sample; ``level`` 'all' or a level string such as '0.1'."""
    if pipe == "P4R4":
        bs = jl(GRAD / "boundary" / "summary.json")["maps"]["P4R4"]["err"]
        key = ("all" if level == "all" else f"level {level}") + "/" + when
        w = bs[key]["wrong"]
        return int(w["k"]), int(w["n"])
    sp = jl(GRAD / "seed" / "summary_part2.json")
    key = f"{pipe} {when}"
    lab = "all" if level == "all" else f"level {level}"
    row = next(x for x in sp["levels"] if x["map"] == key and x["level"] == lab)
    return int(row["wrong_k"]), int(row["n"])


def paired(wa: np.ndarray, wb: np.ndarray) -> Dict[str, int]:
    """A against baseline B, by voxel: fixed = wrong in B and right in A."""
    return dict(
        wrong_a=int(wa.sum()),
        wrong_b=int(wb.sum()),
        fixed=int((wb & ~wa).sum()),
        broken=int((~wb & wa).sum()),
        both_wrong=int((wa & wb).sum()),
    )


# ------------------------------------------------------------------------------------- main
def main() -> None:
    import sys

    if "--slim" in sys.argv:  # needs the full run folders and the pass cache (not committed)
        make_slim()
    z = np.load(BIG / "errors.npz")
    level = z["level"].astype(float)
    err = {f"{p}_{w}": z[f"{p}_{w}"].astype(float) for p in PIPES for w in ("before", "after")}
    wrong = {k: v > WRONG for k, v in err.items()}
    info = jl(BIG / "gradient" / "sample" / "sample.json")
    n_grains = info["n_grains"]
    grains_by_level = {lv(L): int(info["n_by_level"][str(L)]) for L in LEVELS}
    vox_by_level = {lv(L): int((np.isclose(level, L)).sum()) for L in LEVELS}
    bigS = {p: big_stats(p) for p in PIPES}
    grS = grad_stats()
    S: Dict[str, Any] = dict(
        n_voxels=N,
        n_grains_big=n_grains,
        n_grains_gradient=497,
        grains_by_level=grains_by_level,
        voxels_by_level=vox_by_level,
        big=bigS,
        gradient=grS,
    )

    # ---- wrong by level
    lev: List[Dict[str, Any]] = []
    for p in PIPES:
        for when in ("before", "after"):
            w = wrong[f"{p}_{when}"]
            for L in ("all",) + tuple(lv(x) for x in LEVELS):
                sel = np.ones(N, bool) if L == "all" else np.isclose(level, float(L))
                row = dict(pipe=p, when=when, level=L, **ci_d(int(w[sel].sum()), int(sel.sum())))
                row["median_err_deg"] = float(np.median(err[f"{p}_{when}"][sel]))
                lev.append(row)
    S["levels_big"] = lev
    gl: List[Dict[str, Any]] = []
    for p in PIPES:
        for when in ("before", "after"):
            for L in ("all", "0", "0.1", "0.2", "0.5", "1"):
                k, n = grad_wrong(p, when, L)
                gl.append(dict(pipe=p, when=when, level=L, **ci_d(k, n)))
    S["levels_gradient"] = gl

    # ---- matched levels (prediction 3)
    match_rows: List[Dict[str, Any]] = []
    for p in PIPES:
        for bL, gL in MATCH.items():
            b = next(
                x for x in lev if x["pipe"] == p and x["when"] == "after" and x["level"] == lv(bL)
            )
            k, n = grad_wrong(p, "after", lv(gL))
            g = ci_d(k, n)
            status = (
                "higher (intervals do not overlap)"
                if b["lo"] > g["hi"]
                else "lower (intervals do not overlap)" if b["hi"] < g["lo"] else "overlap"
            )
            match_rows.append(
                dict(
                    pipe=p,
                    big_level=lv(bL),
                    grad_level=lv(gL),
                    big=b,
                    grad=g,
                    point="higher" if b["k"] / b["n"] > g["k"] / g["n"] else "not higher",
                    status=status,
                )
            )
    S["matched"] = match_rows

    # ---- paired (all with the pass)
    pr = {}
    for a, b in (("S5", "D4"), ("S5", "D"), ("P4R4", "S5")):
        for when in ("after", "before"):
            pr[f"{a} vs {b} {when}"] = paired(wrong[f"{a}_{when}"], wrong[f"{b}_{when}"])
    S["paired"] = pr

    # ---- predictions
    s5 = bigS["S5"]
    rest = s5["stages"]["neighbour"] + s5["stages"]["revisit"] + s5["stages"]["pass_"]
    seed_grid = s5["stages"]["seed_step"] + s5["stages"]["grid_fallback"]
    g5 = grS["S5"]
    g_rest = g5["stages"]["neighbour"] + g5["stages"]["revisit"] + g5["stages"]["pass_"]
    S["prediction1"] = dict(
        reseeds=s5["reseeds"],
        grains=n_grains,
        reseeds_limit=3 * n_grains,
        reseeds_per_grain=s5["reseeds"] / n_grains,
        gradient_reseeds_per_grain=g5["reseeds"] / 497,
        seed_plus_grid=seed_grid,
        total=s5["total"],
        share=seed_grid / s5["total"],
        agrees=bool(s5["reseeds"] <= 3 * n_grains and seed_grid / s5["total"] < 0.1),
    )
    S["prediction2"] = dict(
        total_per_voxel=s5["total"] / N,
        local_plus_pass_per_voxel=rest / N,
        ratio=s5["total"] / rest,
        gradient_total_per_voxel=g5["total"] / N,
        gradient_local_plus_pass_per_voxel=g_rest / N,
        agrees=bool(s5["total"] / rest <= 1.10),
    )
    s5m = [m for m in match_rows if m["pipe"] == "S5"]
    S["prediction3"] = dict(
        rows=[
            dict(
                big_level=m["big_level"],
                grad_level=m["grad_level"],
                big=m["big"],
                grad=m["grad"],
                point=m["point"],
                status=m["status"],
            )
            for m in s5m
        ],
        any_higher=bool(any(m["status"].startswith("higher") for m in s5m)),
        agrees=bool(not any(m["status"].startswith("higher") for m in s5m)),
    )

    # ---- what the result decides
    w5, w4 = wrong["S5_after"].sum(), wrong["D4_after"].sum()
    lo5, hi5 = wilson(int(w5), N)
    lo4, hi4 = wilson(int(w4), N)
    # ---- added: seed steps of S5 by the true level of the voxel's grain (a reading of the log)
    import gzip

    with gzip.open(BIG / "slim" / "S5_seed_log.json.gz", "rt") as f:
        slog = json.load(f)
    sl: List[Dict[str, Any]] = []
    for L in LEVELS:
        sel = np.isclose(level, L)
        steps = [e for e in slog if sel[e["voxel"]]]
        sl.append(
            dict(
                level=lv(L),
                grains=grains_by_level[lv(L)],
                voxels=int(sel.sum()),
                seed_steps=len(steps),
                reseeds=int(sum(1 for e in steps if e["reseed"])),
                evals=int(sum(e["evals"] for e in steps)),
            )
        )
    assert sum(x["seed_steps"] for x in sl) == len(slog)
    S["seed_steps_by_level"] = sl
    # the same reading on the gradient sample, from its committed log and sample record
    with gzip.open(GRAD / "seed" / "runs" / "S5" / "seed_log.json.gz", "rt") as f:
        glog = json.load(f)
    ginfo = jl(GRAD / "gradient" / "sample" / "sample.json")
    gmap = np.load(C.EX / "SimInput" / "rand_500grains_1mm_neworient_s0_grainmap.npy")
    glev_of = {g["id"]: g["level_deg"] for g in ginfo["grains"]}
    glevel = np.array([glev_of[int(g)] for g in gmap])
    gl2: List[Dict[str, Any]] = []
    for L in (0.0, 0.1, 0.2, 0.5, 1.0):
        sel = np.isclose(glevel, L)
        steps = [e for e in glog if sel[e["voxel"]]]
        gl2.append(
            dict(
                level=lv(L),
                voxels=int(sel.sum()),
                seed_steps=len(steps),
                reseeds=int(sum(1 for e in steps if e["reseed"])),
                evals=int(sum(e["evals"] for e in steps)),
            )
        )
    assert sum(x["seed_steps"] for x in gl2) == len(glog)
    S["seed_steps_by_level_gradient"] = gl2

    S["decision"] = dict(
        s5_wrong=int(w5),
        d4_wrong=int(w4),
        s5_evals=s5["total"],
        d4_evals=bigS["D4"]["total"],
        intervals_overlap=bool(not (lo5 > hi4 or hi5 < lo4)),
        s5_more_wrong=bool(lo5 > hi4),
        s5_fewer_or_equal_evals=bool(s5["total"] <= bigS["D4"]["total"]),
    )

    (BIG / "summary.json").write_text(json.dumps(S, indent=1, default=float) + "\n")
    write_tables(BIG / "tables.md", tables(S))
    plots(S)
    print(json.dumps({k: S[k] for k in ("prediction1", "prediction2", "decision")}, indent=1))


# ---------------------------------------------------------------------------------- tables
def tables(S: Dict[str, Any]) -> Dict[str, str]:
    T: Dict[str, str] = {}
    # evaluations per voxel by stage
    rows: List[Optional[Dict[str, Any]]] = []
    for p in PIPES:
        for side, src in (("big", S["big"]), ("gradient", S["gradient"])):
            d = src[p]
            row: Dict[str, Any] = dict(pipeline=p, sample=side)
            for s in STAGES:
                row[s] = d["stages"][s] / N
            row["total"] = d["total"] / N
            row["total_M"] = d["total"] / 1e6
            rows.append(row)
    cols = ["pipeline", "sample", *STAGES, "total", "total_M"]
    T["evals_per_voxel"] = markdown_table(
        rows, cols, {s: ".1f" for s in (*STAGES, "total")} | {"total_M": ".2f"}
    )
    # counts
    rows = []
    for p in PIPES:
        for side, src in (("big", S["big"]), ("gradient", S["gradient"])):
            d = src[p]
            rows.append(
                dict(
                    pipeline=p,
                    sample=side,
                    seeds=d["seeds"],
                    reseeds=d["reseeds"],
                    grid_fallbacks=d["grid_fallbacks"],
                    wall_run_min=d["wall_run_s"] / 60.0,
                )
            )
    T["counts"] = markdown_table(
        rows,
        ["pipeline", "sample", "seeds", "reseeds", "grid_fallbacks", "wall_run_min"],
        {"wall_run_min": ".1f"},
    )
    # wrong, all voxels, before and after
    rows = []
    for p in PIPES:
        b = {(x["when"]): x for x in S["levels_big"] if x["pipe"] == p and x["level"] == "all"}
        g = {(x["when"]): x for x in S["levels_gradient"] if x["pipe"] == p and x["level"] == "all"}
        rows.append(
            dict(
                pipeline=p,
                big_before=ci(b["before"]["k"], b["before"]["n"]),
                big_after=ci(b["after"]["k"], b["after"]["n"]),
                gradient_before=ci(g["before"]["k"], g["before"]["n"]),
                gradient_after=ci(g["after"]["k"], g["after"]["n"]),
            )
        )
    T["wrong_all"] = markdown_table(
        rows, ["pipeline", "big_before", "big_after", "gradient_before", "gradient_after"]
    )
    # wrong by level on big (after the pass), and before
    rows = []
    for p in PIPES:
        for L in ("0", "0.2", "0.5", "1", "2", "4"):
            b = next(
                x
                for x in S["levels_big"]
                if x["pipe"] == p and x["when"] == "before" and x["level"] == L
            )
            a = next(
                x
                for x in S["levels_big"]
                if x["pipe"] == p and x["when"] == "after" and x["level"] == L
            )
            rows.append(
                dict(
                    pipeline=p,
                    level=L,
                    grains=S["grains_by_level"][L],
                    voxels=a["n"],
                    before=ci(b["k"], b["n"]),
                    after=ci(a["k"], a["n"]),
                    median_after_deg=a["median_err_deg"],
                )
            )
    T["wrong_level_big"] = markdown_table(
        rows,
        ["pipeline", "level", "grains", "voxels", "before", "after", "median_after_deg"],
        {"median_after_deg": ".4f"},
    )
    # matched levels
    rows = []
    for m in S["matched"]:
        rows.append(
            dict(
                pipeline=m["pipe"],
                big_level=m["big_level"],
                big_after=ci(m["big"]["k"], m["big"]["n"]),
                gradient_level=m["grad_level"],
                gradient_after=ci(m["grad"]["k"], m["grad"]["n"]),
                reading=m["status"],
            )
        )
    T["matched"] = markdown_table(
        rows,
        ["pipeline", "big_level", "big_after", "gradient_level", "gradient_after", "reading"],
    )
    # paired
    rows = []
    for k, v in S["paired"].items():
        rows.append(dict(pair=k, **v))
    T["paired"] = markdown_table(
        rows, ["pair", "wrong_a", "wrong_b", "fixed", "broken", "both_wrong"]
    )
    # S5 split
    s5 = S["big"]["S5"]
    rows = []
    for s in STAGES[:3] + ("revisit", "pass_"):
        rows.append(
            dict(stage=s, evals=s5["stages"][s], share_pct=100.0 * s5["stages"][s] / s5["total"])
        )
    rows.append(dict(stage="total", evals=s5["total"], share_pct=100.0))
    T["s5_split"] = markdown_table(rows, ["stage", "evals", "share_pct"], {"share_pct": ".1f"})
    rows = []
    for x in S["seed_steps_by_level"]:
        rows.append(dict(x, per_1000_voxels=1000.0 * x["seed_steps"] / x["voxels"]))
    T["seed_by_level"] = markdown_table(
        rows,
        ["level", "grains", "voxels", "seed_steps", "reseeds", "per_1000_voxels", "evals"],
        {"per_1000_voxels": ".1f"},
    )
    rows = []
    for x in S["seed_steps_by_level_gradient"]:
        rows.append(dict(x, per_1000_voxels=1000.0 * x["seed_steps"] / x["voxels"]))
    T["seed_by_level_gradient"] = markdown_table(
        rows,
        ["level", "voxels", "seed_steps", "reseeds", "per_1000_voxels", "evals"],
        {"per_1000_voxels": ".1f"},
    )
    return T


# ------------------------------------------------------------------------------------- plots
def plots(S: Dict[str, Any]) -> None:
    colors = ["#4c72b0", "#dd8452", "#55a868", "#c44e52", "#8172b3", "#937860"]
    labs = ["seed step", "grid fallbacks", "neighbour", "revisit / R4", "scoring", "pass"]
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2), sharey=True)
    for ax, side in zip(axes, ("gradient", "big")):
        x = np.arange(len(PIPES))
        bottom = np.zeros(len(PIPES))
        for s, c, lab in zip(STAGES, colors, labs):
            v = np.array([S[side][p]["stages"][s] / N for p in PIPES])
            ax.bar(x, v, bottom=bottom, color=c, label=lab)
            bottom += v
        ax.set_xticks(x)
        ax.set_xticklabels(PIPES)
        ax.set_title(("gradient sample (497 grains)" if side == "gradient" else "big (32 grains)"))
        ax.set_ylabel("evaluations per voxel")
    axes[1].legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(BIG / "ar3_big_evals.png", dpi=110)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(7.5, 4.2))
    w = 0.38
    for j, side in enumerate(("gradient", "big")):
        for i, p in enumerate(PIPES):
            src = S["levels_gradient"] if side == "gradient" else S["levels_big"]
            r = next(
                x for x in src if x["pipe"] == p and x["when"] == "after" and x["level"] == "all"
            )
            share = r["k"] / r["n"]
            ax.bar(
                i + (j - 0.5) * w,
                100 * share,
                w,
                color=["#8c8c8c", "#4c72b0"][j],
                label=side if i == 0 else None,
            )
            ax.errorbar(
                i + (j - 0.5) * w,
                100 * share,
                yerr=[[100 * (share - r["lo"])], [100 * (r["hi"] - share)]],
                color="k",
                capsize=3,
            )
    ax.set_xticks(range(len(PIPES)))
    ax.set_xticklabels(PIPES)
    ax.set_ylabel("wrong voxels after the pass (%)")
    ax.legend()
    fig.tight_layout()
    fig.savefig(BIG / "ar3_big_wrong.png", dpi=110)
    plt.close(fig)


if __name__ == "__main__":
    main()
