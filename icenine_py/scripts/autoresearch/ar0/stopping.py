"""AR0: what the stored results say about stopping a search early (no new searches).

Reads, read-only, from the main checkout (--eval-root = benchmarks/eval_reduction):
  E/cache/rank/e0_{old,cpp}_base/v*.npz   E0, 200 hard voxels, one no-start search each, written by
                                          scripts/eval_reduction/smooth_rank_probe.py (spec "off").
  G/g2/bfs/runs/{keyoff_rerun,outward_only,integrated}, G/g3/bfs/runs/smooth_only
                                          2,000-voxel BFS pilot runs (records.json, recon.mic).

Sign convention: the cost is LOWER = BETTER (0 is a perfect match). ``L{k}_cost`` is the strict cost
of the candidates after the quick Monte-Carlo at level k, best first; ``cost_final`` is the cost of
the returned orientation after the finisher. Evaluations per level = ``L{k}_screen_recost`` +
``L{k}_quick``; the rest of ``evals`` is the finisher.

  uv run python scripts/autoresearch/ar0/stopping.py \
      --eval-root <main>/icenine_py/benchmarks/eval_reduction

Writes summary.json, tables.md, fig_stopping.png, fig_pilot_hit_ratio.png to --out
(default benchmarks/autoresearch/ar0). Counts only; no p-values.
"""

import argparse
import glob
import json
import shutil
import sys
import tempfile
from pathlib import Path
from typing import Any, Dict, List

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(ROOT / "scripts" / "common"))
sys.path.insert(0, str(ROOT / "scripts" / "phase_d"))
import doc_tables as DT  # noqa: E402

WRONG_DEG = 1.0
FINISH_DEG = 2.0  # a start within 2 degrees is enough for the finisher (plan, Task F)
LEVELS = 4
C_GRID = [0.2, 0.3, 0.4, 0.5, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95]


def load_e0(root: Path, rule: str) -> Dict[str, np.ndarray]:
    fs = sorted(glob.glob(str(root / "E" / "cache" / "rank" / f"e0_{rule}_base" / "v*.npz")))
    assert len(fs) == 200, (rule, len(fs))
    out: Dict[str, List[Any]] = {k: [] for k in ["voxel", "err", "evals", "cost_final"]}
    for k in range(LEVELS):
        for s in [
            "best_cost",
            "best_dist",
            "kept_min_dist",
            "cum_evals",
            "argmin_ok",
            "screen",
            "quick",
            "n_disc",
            "first_basin",
        ]:
            out[f"L{k}_{s}"] = []
    for f in fs:
        d = np.load(f)
        assert int(d["n_levels"]) == LEVELS
        out["voxel"].append(int(d["voxel"]))
        out["err"].append(float(d["err_deg"]))
        out["evals"].append(int(d["evals"]))
        out["cost_final"].append(float(d["cost_final"]))
        cum = 0
        for k in range(LEVELS):
            cost, dist = d[f"L{k}_cost"], d[f"L{k}_dist"]
            nk = int(d[f"L{k}_n_keep"])
            cum += int(d[f"L{k}_screen_recost"]) + int(d[f"L{k}_quick"])
            out[f"L{k}_best_cost"].append(float(cost[0]))
            out[f"L{k}_best_dist"].append(float(dist[0]))
            out[f"L{k}_kept_min_dist"].append(float(dist[:nk].min()))
            out[f"L{k}_cum_evals"].append(cum)
            out[f"L{k}_screen"].append(int(d[f"L{k}_screen_recost"]))
            out[f"L{k}_quick"].append(int(d[f"L{k}_quick"]))
            nd = int(d[f"L{k}_n_disc"])
            out[f"L{k}_n_disc"].append(nd)
            dd = d[f"L{k}_dist_disc"]
            hit = np.nonzero(dd <= np.degrees(float(d[f"L{k}_diameter"])))[0]
            out[f"L{k}_first_basin"].append(int(hit[0]) if len(hit) else -1)
            out[f"L{k}_argmin_ok"].append(bool(np.argmin(cost) == 0 or cost[0] == cost.min()))
    return {k: np.array(v) for k, v in out.items()}


def separation(cost: np.ndarray, wrong: np.ndarray) -> Dict[str, Any]:
    """Counts either side of thresholds; pure thresholds (no mistakes) and the gap."""
    right_max, wrong_min = float(cost[~wrong].max()), float(cost[wrong].min())
    rows = []
    for c in C_GRID:
        rows.append(
            dict(
                c=c,
                right_le=int((~wrong & (cost <= c)).sum()),
                wrong_le=int((wrong & (cost <= c)).sum()),
                right_gt=int((~wrong & (cost > c)).sum()),
                wrong_gt=int((wrong & (cost > c)).sum()),
            )
        )
    return dict(
        right_max=right_max,
        wrong_min=wrong_min,
        n_right=int((~wrong).sum()),
        n_wrong=int(wrong.sum()),
        separable=bool(right_max < wrong_min),
        rows=rows,
    )


def level_table(E: Dict[str, np.ndarray]) -> List[Dict[str, Any]]:
    wrong = E["err"] > WRONG_DEG
    tot = float(E["evals"].sum())
    rows = []
    for k in range(LEVELS):
        bd, kd = E[f"L{k}_best_dist"], E[f"L{k}_kept_min_dist"]
        after = E["evals"] - E[f"L{k}_cum_evals"]
        rows.append(
            dict(
                level=k,
                best_le1_right=int((~wrong & (bd <= WRONG_DEG)).sum()),
                best_le2_right=int((~wrong & (bd <= FINISH_DEG)).sum()),
                best_le2_wrong=int((wrong & (bd <= FINISH_DEG)).sum()),
                kept_le2_right=int((~wrong & (kd <= FINISH_DEG)).sum()),
                kept_le2_wrong=int((wrong & (kd <= FINISH_DEG)).sum()),
                evals_after_pct=100.0 * float(after.sum()) / tot,
                evals_after_right_mean=float(after[~wrong].mean()),
            )
        )
    return rows


def stage_table(E: Dict[str, np.ndarray]) -> List[Dict[str, Any]]:
    """Where the evaluations go (all 200 voxels) and where the true basin sits in the screened list.

    first_basin = position (0 = best screen score) of the first discrete candidate within the
    level's grid diameter of the truth; reported for the baseline-right voxels that have one."""
    wrong = E["err"] > WRONG_DEG
    tot = float(E["evals"].sum())
    rows = []
    for k in range(LEVELS):
        fb, nd = E[f"L{k}_first_basin"], E[f"L{k}_n_disc"]
        m = ~wrong & (fb >= 0)
        frac = fb[m] / nd[m]
        rows.append(
            dict(
                level=k,
                screen_pct=100.0 * float(E[f"L{k}_screen"].sum()) / tot,
                quick_pct=100.0 * float(E[f"L{k}_quick"].sum()) / tot,
                n_disc_median=float(np.median(nd)),
                right_with_basin=int(m.sum()),
                right_no_basin=int((~wrong & (fb < 0)).sum()),
                wrong_with_basin=int((wrong & (fb >= 0)).sum()),
                basin_pos_median=float(np.median(fb[m])),
                basin_frac_median=float(np.median(frac)),
                basin_frac_p90=float(np.percentile(frac, 90)),
            )
        )
    return rows


def stop_grid(E: Dict[str, np.ndarray], tol: float, cs: List[float]) -> List[Dict[str, Any]]:
    """Rule: after level k, stop if best strict cost <= c; answer = that level's best candidate.

    broken = baseline right but the stopped candidate is farther than ``tol`` from the truth.
    saved = evaluations after level k (the finisher of the stopped voxel is not charged: gross)."""
    wrong = E["err"] > WRONG_DEG
    tot = float(E["evals"].sum())
    rows = []
    for k in range(LEVELS - 1):
        cost, dist = E[f"L{k}_best_cost"], E[f"L{k}_best_dist"]
        saved_v = E["evals"] - E[f"L{k}_cum_evals"]
        for c in cs:
            st = cost <= c
            rows.append(
                dict(
                    level=k,
                    c=c,
                    n_stop=int(st.sum()),
                    stop_ok=int((st & ~wrong & (dist <= tol)).sum()),
                    broken=int((st & ~wrong & (dist > tol)).sum()),
                    stop_wrong_base=int((st & wrong).sum()),
                    stop_wrong_base_but_near=int((st & wrong & (dist <= tol)).sum()),
                    saved_pct=100.0 * float(saved_v[st].sum()) / tot,
                )
            )
    return rows


def best_safe(rows: List[Dict[str, Any]], max_broken: int) -> Dict[int, Dict[str, Any]]:
    out: Dict[int, Dict[str, Any]] = {}
    for k in range(LEVELS - 1):
        ok = [r for r in rows if r["level"] == k and r["broken"] <= max_broken]
        if ok:
            out[k] = max(ok, key=lambda r: (r["saved_pct"], -r["c"]))
    return out


# ----------------------------------------------------------------------------- pilot


def pilot(runs: Dict[str, Path], tmp: Path) -> Dict[str, Any]:
    import score_bfs as SB

    res: Dict[str, Any] = {}
    for name, p in runs.items():
        d = tmp / name
        d.mkdir(parents=True)
        for f in ["records.json", "recon.mic", "run.json"]:
            shutil.copy(p / f, d / f)
        s = SB.score_run(d)
        v = np.load(d / "voxels.npz")
        rec = {int(k): r for k, r in json.loads((d / "records.json").read_text()).items()}
        n = len(v["err_deg"])
        hit = np.array([rec[i]["hit_ratio"] if i in rec else np.nan for i in range(n)])
        res[name] = dict(
            err=v["err_deg"],
            src=v["source"],
            hit=hit,
            nev=v["n_evals"],
            reseeded=np.array([bool(rec[i]["reseeded"]) if i in rec else False for i in range(n)]),
            wrong_total=s["wrong_total"],
        )
    return res


def pilot_tables(P: Dict[str, Any]) -> Dict[str, Any]:
    out: Dict[str, Any] = {}
    for name, r in P.items():
        wrong = r["err"] > WRONG_DEG
        ok = ~np.isnan(r["hit"])
        full = (r["nev"] > 50000) & ok  # a full search (seed or reseed); local fits use <= ~2,000
        groups = dict(all=ok, full_search=full, local_fit=ok & ~full)
        g_out: Dict[str, Any] = {}
        for gname, m in groups.items():
            rows = []
            for h in [0.90, 0.93, 0.95, 0.96, 0.97, 0.98, 0.99]:
                lo = m & (r["hit"] < h)
                hi = m & (r["hit"] >= h)
                rows.append(
                    dict(
                        h=h,
                        right_ge=int((hi & ~wrong).sum()),
                        wrong_ge=int((hi & wrong).sum()),
                        right_lt=int((lo & ~wrong).sum()),
                        wrong_lt=int((lo & wrong).sum()),
                        evals_in_wrong_pct=100.0
                        * float(r["nev"][m & wrong].sum())
                        / max(float(r["nev"][m].sum()), 1.0),
                    )
                )
            g_out[gname] = dict(
                n=int(m.sum()),
                n_wrong=int((m & wrong).sum()),
                hit_right_min=float(r["hit"][m & ~wrong].min()) if (m & ~wrong).any() else None,
                hit_wrong_max=float(r["hit"][m & wrong].max()) if (m & wrong).any() else None,
                hit_wrong_median=(
                    float(np.median(r["hit"][m & wrong])) if (m & wrong).any() else None
                ),
                hit_right_median=float(np.median(r["hit"][m & ~wrong])),
                rows=rows,
            )
        out[name] = dict(wrong_total=r["wrong_total"], groups=g_out)
    return out


# ----------------------------------------------------------------------------- figures


def figures(
    E: Dict[str, Dict[str, np.ndarray]], grid: List[Dict[str, Any]], P: Dict[str, Any], out: Path
) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    e = E["old"]
    wrong = e["err"] > WRONG_DEG
    rng = np.random.default_rng(0)  # plotting jitter only
    fig, ax = plt.subplots(1, 3, figsize=(16, 4.6))
    for j, (nm, lab) in enumerate(
        [("L0", "level 0"), ("L1", "level 1"), ("L2", "level 2"), ("L3", "level 3")]
    ):
        pass
    # (1) cost vs final error: final cost and the best cost after each level
    names = ["L0", "L1", "L2", "L3", "final"]
    for i, nm in enumerate(names):
        c = e["cost_final"] if nm == "final" else e[f"{nm}_best_cost"]
        x = i + rng.uniform(-0.3, 0.3, len(c))
        ax[0].scatter(
            x[~wrong], c[~wrong], s=9, c="tab:blue", label="right (<=1 deg)" if i == 0 else None
        )
        ax[0].scatter(
            x[wrong], c[wrong], s=9, c="tab:red", label="wrong (>1 deg)" if i == 0 else None
        )
    ax[0].set_xticks(range(5))
    ax[0].set_xticklabels(["after L0", "after L1", "after L2", "after L3", "final\n(finisher)"])
    ax[0].set_ylabel("best strict cost (lower = better)")
    ax[0].set_title("E0, old rule, 200 voxels: cost of the best candidate")
    ax[0].legend(fontsize=8)
    # (2) distance of best candidate to truth vs its cost, per level
    for k, col in zip(range(3), ["tab:orange", "tab:green", "tab:purple"]):
        ax[1].scatter(
            e[f"L{k}_best_cost"],
            np.maximum(e[f"L{k}_best_dist"], 0.02),
            s=8,
            c=col,
            label=f"level {k}",
        )
    ax[1].set_yscale("log")
    ax[1].axhline(FINISH_DEG, c="k", ls=":", lw=1)
    ax[1].axhline(WRONG_DEG, c="k", ls="--", lw=1)
    ax[1].set_xlabel("best strict cost after the level")
    ax[1].set_ylabel("distance of that candidate to truth (deg)")
    ax[1].set_title("Is the best candidate in the right basin?")
    ax[1].legend(fontsize=8)
    # (3) saved vs broken
    for k, col in zip(range(3), ["tab:orange", "tab:green", "tab:purple"]):
        rr = [r for r in grid if r["level"] == k]
        ax[2].plot(
            [r["broken"] for r in rr],
            [r["saved_pct"] for r in rr],
            "o-",
            c=col,
            label=f"stop after level {k}",
        )
    ax[2].set_xlabel("voxels broken (of 159 right in baseline; candidate > 2 deg off)")
    ax[2].set_ylabel("evaluations saved (% of baseline total)")
    ax[2].set_title("Stop rule: best cost <= c (c swept 0.2 to 0.95)")
    ax[2].legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(out / "fig_stopping.png", dpi=130)
    plt.close(fig)

    # pilot hit ratio
    names = [n for n in P if n.startswith("keyoff")] + [n for n in P if not n.startswith("keyoff")]
    fig, ax = plt.subplots(1, 2, figsize=(12, 4.4))
    r = P["keyoff_rerun"]
    wr = r["err"] > WRONG_DEG
    ok = ~np.isnan(r["hit"])
    full = (r["nev"] > 50000) & ok
    for a, m, t in [
        (ax[0], full, "full searches (seed and reseed)"),
        (ax[1], ok & ~full, "local fits"),
    ]:
        bins = np.linspace(0.5, 1.0, 51)
        a.hist(
            r["hit"][m & ~wr],
            bins=bins,
            alpha=0.7,
            label=f"right ({int((m & ~wr).sum())})",
            color="tab:blue",
        )
        a.hist(
            r["hit"][m & wr],
            bins=bins,
            alpha=0.7,
            label=f"wrong ({int((m & wr).sum())})",
            color="tab:red",
        )
        a.set_yscale("log")
        a.set_xlabel("hit ratio of the returned orientation (higher = better)")
        a.set_title(f"Pilot key-off rerun: {t}")
        a.legend()
    fig.tight_layout()
    fig.savefig(out / "fig_pilot_hit_ratio.png", dpi=130)
    plt.close(fig)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    ap.add_argument("--eval-root", required=True, help="main checkout benchmarks/eval_reduction")
    ap.add_argument("--out", default=str(ROOT / "benchmarks" / "autoresearch" / "ar0"))
    a = ap.parse_args()
    root, out = Path(a.eval_root), Path(a.out)
    out.mkdir(parents=True, exist_ok=True)

    E = {rule: load_e0(root, rule) for rule in ["old", "cpp"]}
    summ: Dict[str, Any] = {
        "sign": "cost lower = better",
        "wrong_deg": WRONG_DEG,
        "finish_deg": FINISH_DEG,
    }
    for rule, e in E.items():
        wrong = e["err"] > WRONG_DEG
        d: Dict[str, Any] = dict(
            n=int(len(wrong)),
            n_wrong=int(wrong.sum()),
            evals_total=int(e["evals"].sum()),
            final=separation(e["cost_final"], wrong),
            levels=level_table(e),
            stages=stage_table(e),
            finisher_pct=100.0
            * float((e["evals"] - e[f"L{LEVELS-1}_cum_evals"]).sum())
            / float(e["evals"].sum()),
            argmin_ok_all=bool(all(e[f"L{k}_argmin_ok"].all() for k in range(LEVELS))),
            grid_tol2=stop_grid(e, FINISH_DEG, C_GRID),
            grid_tol1=stop_grid(e, WRONG_DEG, C_GRID),
        )
        d["levels_sep"] = {f"L{k}": separation(e[f"L{k}_best_cost"], wrong) for k in range(LEVELS)}
        d["safe0_tol2"] = {str(k): v for k, v in best_safe(d["grid_tol2"], 0).items()}
        d["safe1_tol2"] = {str(k): v for k, v in best_safe(d["grid_tol2"], 1).items()}
        summ[rule] = d
    # threshold picked on the old-rule voxels, applied to the C++-rule run of the same voxels
    cross = []
    for k, r in summ["old"]["safe0_tol2"].items():
        e = E["cpp"]
        wrong = e["err"] > WRONG_DEG
        st = e[f"L{k}_best_cost"] <= r["c"]
        cross.append(
            dict(
                level=int(k),
                c=r["c"],
                n_stop=int(st.sum()),
                broken=int((st & ~wrong & (e[f"L{k}_best_dist"] > FINISH_DEG)).sum()),
                saved_pct=100.0
                * float((e["evals"] - e[f"L{k}_cum_evals"])[st].sum())
                / float(e["evals"].sum()),
            )
        )
    summ["cross_old_to_cpp"] = cross

    root_g = root / "G"
    runs = {
        "keyoff_rerun": root_g / "g2/bfs/runs/keyoff_rerun",
        "outward_only": root_g / "g2/bfs/runs/outward_only",
        "integrated": root_g / "g2/bfs/runs/integrated",
        "smooth_only": root_g / "g3/bfs/runs/smooth_only",
    }
    with tempfile.TemporaryDirectory() as td:
        P = pilot(runs, Path(td))
    summ["pilot"] = pilot_tables(P)
    figures(E, summ["old"]["grid_tol2"], P, out)
    (out / "summary.json").write_text(json.dumps(summ, indent=1))

    # ---- tables
    T: Dict[str, str] = {}
    for rule in ["old", "cpp"]:
        s = summ[rule]
        T[f"final_sep_{rule}"] = DT.markdown_table(
            s["final"]["rows"],
            ["c", "right_le", "wrong_le", "right_gt", "wrong_gt"],
            {"c": ".2f"},
        )
        T[f"levels_{rule}"] = DT.markdown_table(
            s["levels"],
            [
                "level",
                "best_le1_right",
                "best_le2_right",
                "best_le2_wrong",
                "kept_le2_right",
                "kept_le2_wrong",
                "evals_after_pct",
                "evals_after_right_mean",
            ],
            {"evals_after_pct": ".1f", "evals_after_right_mean": ".0f"},
        )
        T[f"stages_{rule}"] = DT.markdown_table(
            s["stages"],
            [
                "level",
                "screen_pct",
                "quick_pct",
                "n_disc_median",
                "right_with_basin",
                "right_no_basin",
                "wrong_with_basin",
                "basin_pos_median",
                "basin_frac_median",
                "basin_frac_p90",
            ],
            {
                "screen_pct": ".1f",
                "quick_pct": ".1f",
                "n_disc_median": ".0f",
                "basin_pos_median": ".0f",
                "basin_frac_median": ".2f",
                "basin_frac_p90": ".2f",
            },
        )
        sel = [r for r in s["grid_tol2"] if r["c"] in (0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)]
        T[f"stop_grid_{rule}"] = DT.markdown_table(
            sel,
            ["level", "c", "n_stop", "stop_ok", "broken", "stop_wrong_base", "saved_pct"],
            {"c": ".2f", "saved_pct": ".1f"},
        )
        T[f"sep_levels_{rule}"] = DT.markdown_table(
            [
                dict(
                    level=k,
                    **{
                        x: v[x]
                        for x in ["right_max", "wrong_min", "n_right", "n_wrong", "separable"]
                    },
                )
                for k, v in s["levels_sep"].items()
            ],
            ["level", "right_max", "wrong_min", "n_right", "n_wrong", "separable"],
            {"right_max": ".3f", "wrong_min": ".3f"},
        )
        safe = []
        for mb, key in [(0, "safe0_tol2"), (1, "safe1_tol2")]:
            for k, r in s[key].items():
                safe.append(
                    dict(
                        max_broken=mb,
                        level=int(k),
                        c=r["c"],
                        n_stop=r["n_stop"],
                        broken=r["broken"],
                        saved_pct=r["saved_pct"],
                    )
                )
        T[f"safe_{rule}"] = DT.markdown_table(
            safe,
            ["max_broken", "level", "c", "n_stop", "broken", "saved_pct"],
            {"c": ".2f", "saved_pct": ".1f"},
        )
    T["cross_old_to_cpp"] = DT.markdown_table(
        cross, ["level", "c", "n_stop", "broken", "saved_pct"], {"c": ".2f", "saved_pct": ".1f"}
    )
    prow = []
    for name, p in summ["pilot"].items():
        for g, v in p["groups"].items():
            if g == "all":
                continue
            prow.append(
                dict(
                    run=name,
                    group=g,
                    n=v["n"],
                    n_wrong=v["n_wrong"],
                    hit_right_min=v["hit_right_min"],
                    hit_right_median=v["hit_right_median"],
                    hit_wrong_median=v["hit_wrong_median"],
                    hit_wrong_max=v["hit_wrong_max"],
                )
            )
    T["pilot_hit"] = DT.markdown_table(
        prow,
        [
            "run",
            "group",
            "n",
            "n_wrong",
            "hit_right_min",
            "hit_right_median",
            "hit_wrong_median",
            "hit_wrong_max",
        ],
        {
            k: ".3f"
            for k in ["hit_right_min", "hit_right_median", "hit_wrong_median", "hit_wrong_max"]
        },
    )
    for gname in ["full_search", "local_fit"]:
        rows = summ["pilot"]["keyoff_rerun"]["groups"][gname]["rows"]
        T[f"pilot_thr_{gname}"] = DT.markdown_table(
            rows,
            ["h", "right_ge", "wrong_ge", "right_lt", "wrong_lt", "evals_in_wrong_pct"],
            {"h": ".2f", "evals_in_wrong_pct": ".1f"},
        )
    DT.write_tables(out / "tables.md", T)
    print("wrote", out)


if __name__ == "__main__":
    main()
