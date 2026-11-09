"""AR0: what the stored results say about stopping a search early (no new searches).

Reads, read-only, from the main checkout (--eval-root = benchmarks/eval_reduction):
  E/cache/rank/e0_{old,cpp}_base/v*.npz   E0, 200 hard voxels, one no-start search each, written by
                                          scripts/eval_reduction/smooth_rank_probe.py (spec "off").
  G/g2/bfs/runs/{keyoff_rerun,outward_only,integrated}, G/g3/bfs/runs/smooth_only
                                          2,000-voxel BFS pilot runs (records.json, recon.mic,
                                          run.json).
  G/g2/ledger/runs/g2_base_v*.json        six full-sample voxels, one no-start search each, with
                                          evaluations and wall time per stage (Task G2 ledger).

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
            "diameter_deg",
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
            # "screened list" = R_disc: the survivors of the spacing filter, ordered by screen
            # score, before the quick MC (every member goes through the quick MC). "kept" = the
            # first n_keep entries after the quick MC, re-ordered by strict cost. True basin at
            # level k = the first screened candidate within the level's grid diameter of the truth.
            nd = int(d[f"L{k}_n_disc"])
            out[f"L{k}_n_disc"].append(nd)
            dd = d[f"L{k}_dist_disc"]
            diam = float(np.degrees(float(d[f"L{k}_diameter"])))
            out[f"L{k}_diameter_deg"].append(diam)
            hit = np.nonzero(dd <= diam)[0]
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
        bc = E[f"L{k}_best_cost"]
        after = E["evals"] - E[f"L{k}_cum_evals"]
        far = bd > FINISH_DEG
        rows.append(
            dict(
                level=k,
                # lowest best-candidate cost among best candidates more than 2 degrees off (all
                # 200 voxels): below it every best candidate is within 2 degrees
                far_cost_min=float(bc[far].min()) if far.any() else None,
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
    """Where the evaluations go (all 200 voxels, raw counts) and where the true basin sits in the
    screened list.

    first_basin = position (0 = best screen score) of the first screened candidate (R_disc, before
    the quick MC) within the level's grid diameter of the truth; the position medians are over the
    baseline-right (resp. baseline-wrong) voxels that have one."""
    wrong = E["err"] > WRONG_DEG
    tot = float(E["evals"].sum())
    rows = []
    for k in range(LEVELS):
        fb, nd = E[f"L{k}_first_basin"], E[f"L{k}_n_disc"]
        m = ~wrong & (fb >= 0)
        mw = wrong & (fb >= 0)
        frac = fb[m] / nd[m]
        rows.append(
            dict(
                level=k,
                diameter_deg=float(E[f"L{k}_diameter_deg"][0]),
                screen_pct=100.0 * float(E[f"L{k}_screen"].sum()) / tot,
                quick_pct=100.0 * float(E[f"L{k}_quick"].sum()) / tot,
                n_disc_median=float(np.median(nd)),
                right_with_basin=int(m.sum()),
                right_no_basin=int((~wrong & (fb < 0)).sum()),
                wrong_with_basin=int((wrong & (fb >= 0)).sum()),
                basin_pos_median=float(np.median(fb[m])),
                wrong_basin_pos_median=float(np.median(fb[mw])) if mw.any() else None,
                basin_frac_median=float(np.median(frac)),
                basin_frac_p90=float(np.percentile(frac, 90)),
            )
        )
    return rows


def stop_grid(E: Dict[str, np.ndarray], tol: float, cs: List[float]) -> List[Dict[str, Any]]:
    """Rule: after level k, stop if best strict cost <= c; answer = that level's best candidate.

    Raw evaluation counts (the level-0 screen counts interpolated evaluations; not
    strict-equivalent). Two accountings of the saving:
      saved_pct / saved_gross_pct: every evaluation after level k, the finisher included; i.e. the
        stopped voxel returns the level-k best candidate with no finisher. Under this accounting
        "right" means within 1 degree, so read ``broken`` from the tol = 1 grid.
      saved_net_pct: the finisher still runs on the stopped voxel, at its baseline cost; only the
        later levels are saved (L3_cum - Lk_cum). ``broken`` at tol = 2 assumes the finisher reaches
        the truth from within 2 degrees, which is untested for the MC finisher.
    broken = baseline right but the stopped candidate is farther than ``tol`` from the truth."""
    wrong = E["err"] > WRONG_DEG
    tot = float(E["evals"].sum())
    last = E[f"L{LEVELS - 1}_cum_evals"]
    rows = []
    for k in range(LEVELS - 1):
        cost, dist = E[f"L{k}_best_cost"], E[f"L{k}_best_dist"]
        saved_v = E["evals"] - E[f"L{k}_cum_evals"]
        net_v = last - E[f"L{k}_cum_evals"]
        for c in cs:
            st = cost <= c
            gross = 100.0 * float(saved_v[st].sum()) / tot
            rows.append(
                dict(
                    level=k,
                    c=c,
                    n_stop=int(st.sum()),
                    stop_ok=int((st & ~wrong & (dist <= tol)).sum()),
                    broken=int((st & ~wrong & (dist > tol)).sum()),
                    stop_wrong_base=int((st & wrong).sum()),
                    stop_wrong_base_but_near=int((st & wrong & (dist <= tol)).sum()),
                    saved_pct=gross,
                    saved_gross_pct=gross,
                    saved_net_pct=100.0 * float(net_v[st].sum()) / tot,
                )
            )
    return rows


def merge_grids(g1: List[Dict[str, Any]], g2: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
    """One row per (level, c) with broken counted at 1 and at 2 degrees."""
    out = []
    for a, b in zip(g1, g2):
        assert (a["level"], a["c"]) == (b["level"], b["c"])
        out.append(
            dict(
                level=a["level"],
                c=a["c"],
                n_stop=a["n_stop"],
                broken_1deg=a["broken"],
                broken_2deg=b["broken"],
                stop_wrong_base=a["stop_wrong_base"],
                saved_gross_pct=a["saved_gross_pct"],
                saved_net_pct=a["saved_net_pct"],
            )
        )
    return out


def best_safe(rows: List[Dict[str, Any]], max_broken: int) -> Dict[int, Dict[str, Any]]:
    """Per level, the row with at most ``max_broken`` broken voxels at 1 degree (the stricter
    count, which is never below the 2-degree count) and the largest gross saving; on a tie in the
    saving, the smaller c. Rows come from ``merge_grids``."""
    out: Dict[int, Dict[str, Any]] = {}
    for k in range(LEVELS - 1):
        ok = [r for r in rows if r["level"] == k and r["broken_1deg"] <= max_broken]
        if ok:
            out[k] = max(ok, key=lambda r: (r["saved_gross_pct"], -r["c"]))
    return out


STAGES = [
    "level 0 screen and re-cost",
    "level 0 quick MC",
    "level 1",
    "level 2",
    "level 3",
    "finisher",
]


def e0_split(E: Dict[str, np.ndarray]) -> Dict[str, float]:
    """E0's raw evaluations by stage, in % of all its evaluations (all 200 voxels)."""
    tot = float(E["evals"].sum())
    v = [float(E["L0_screen"].sum()), float(E["L0_quick"].sum())]
    v += [float((E[f"L{k}_screen"] + E[f"L{k}_quick"]).sum()) for k in range(1, LEVELS)]
    v.append(tot - float(E[f"L{LEVELS - 1}_cum_evals"].sum()))
    return {s: 100.0 * x / tot for s, x in zip(STAGES, v)}


def full_sample_split(ledger: Path) -> Dict[str, Any]:
    """The six full-sample voxels of the G2 ledger (baseline, no options): raw evaluations and
    wall time by stage, summed over the six voxels.

    Source: <eval-root>/G/g2/ledger/runs/g2_base_v*.json (one no-start search per voxel, written by
    the Task G2 ledger). Wall times are as labelled in the files (single-worker, contended; the
    preflight was not quiet), so the time shares are indicative only. The level-0 "screen" counts
    interpolated evaluations (one per orientation, summed with strict ones), so the raw counts are
    not strict-equivalent."""
    fs = sorted(glob.glob(str(ledger / "g2_base_v*.json")))
    assert len(fs) == 6, len(fs)
    ev = np.zeros(len(STAGES))
    ts = np.zeros(len(STAGES))
    wall = 0.0
    screen_ev = recost_ev = 0
    screen_s = 0.0
    labels = set()
    voxels = []
    for f in fs:
        d = json.loads(Path(f).read_text())
        voxels.append(int(d["voxel"]))
        labels.add(str(d["label"]))
        lv = {int(L["level"]): L for L in d["levels"]}
        L0 = lv[0]
        e = [L0["screen_evals"] + L0["recost_evals"], L0["quick_evals"]]
        e += [lv[k]["stage_evals"] for k in range(1, LEVELS)]
        fin = sum(int(x["evals"]) for x in d["find"]) + d["variance"]["evals"] + d["final"]["evals"]
        e.append(fin)
        assert sum(e) == d["evals_total"] == d["totals"]["sum"], f
        t = [L0["discrete_s"], L0["quick_s"]]
        t += [lv[k]["discrete_s"] + lv[k]["quick_s"] for k in range(1, LEVELS)]
        t.append(
            sum(float(x["wall_s"]) for x in d["find"])
            + d["variance"]["wall_s"]
            + d["final"]["wall_s"]
        )
        ev += np.array(e, dtype=float)
        ts += np.array(t)
        wall += float(d["wall_s"])
        screen_ev += int(L0["screen_evals"])
        recost_ev += int(L0["recost_evals"])
        screen_s += float(L0["screen_eval_s"])
    tot = float(ev.sum())
    n = len(fs)
    return dict(
        source="G/g2/ledger/runs/g2_base_v*.json (six full-sample voxels, baseline)",
        voxels=voxels,
        n=n,
        timing=" / ".join(sorted(labels)),
        evals_total=int(tot),
        evals_per_voxel_mean=int(round(tot / n)),
        evals_pct={s: 100.0 * x / tot for s, x in zip(STAGES, ev)},
        evals_per_voxel={s: int(round(x / n)) for s, x in zip(STAGES, ev)},
        level0_screen_per_voxel=int(round(screen_ev / n)),
        level0_recost_per_voxel=int(round(recost_ev / n)),
        level0_screen_evals_pct=100.0 * screen_ev / tot,
        level0_screen_time_pct=100.0 * screen_s / wall,
        wall_s_total=wall,
        time_pct={s: 100.0 * x / wall for s, x in zip(STAGES, ts)},
        time_other_pct=100.0 * (wall - float(ts.sum())) / wall,
    )


def split_tables(e0: Dict[str, float], full: Dict[str, Any]) -> Dict[str, Any]:
    rows = [
        dict(
            stage=s,
            e0_evals_pct=e0[s],
            full_evals_pct=full["evals_pct"][s],
            full_time_pct=full["time_pct"][s],
        )
        for s in STAGES
    ]
    rows.append(
        dict(
            stage="not in a stage (time only)",
            e0_evals_pct=None,
            full_evals_pct=None,
            full_time_pct=full["time_other_pct"],
        )
    )
    after = []
    for k in range(LEVELS):
        tail = STAGES[2 + k :]  # stages after level k (later levels and the finisher)
        after.append(
            dict(
                after_level=k,
                e0_evals_pct=sum(e0[s] for s in tail),
                full_evals_pct=sum(full["evals_pct"][s] for s in tail),
                full_time_pct=sum(full["time_pct"][s] for s in tail),
            )
        )
    return dict(stages=rows, after=after)


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

        def flag(key: str) -> np.ndarray:
            return np.array([bool(rec[i][key]) if i in rec else False for i in range(n)])

        stats = json.loads((d / "run.json").read_text())["stats"]
        res[name] = dict(
            err=v["err_deg"],
            hit=hit,
            nev=v["n_evals"],
            # a seed or reseed voxel: _fit_from_seed sets source "seed"; a rejected seed keeps
            # seed_rejected after a later revisit overwrites its source; a reseed sets reseeded.
            # In OutwardSearch runs a reseed can be answered by the policy without a full search.
            full=(v["source"] == "seed") | flag("seed_rejected") | flag("reseeded"),
            wrong_total=s["wrong_total"],
            n_full_run=int(stats["n_seeds"]) + int(stats["n_reseeds"]),
            evals_full_run=int(stats["n_evals_seed"]) + int(stats["n_evals_reseed"]),
        )
    return res


def pilot_tables(P: Dict[str, Any]) -> Dict[str, Any]:
    out: Dict[str, Any] = {}
    for name, r in P.items():
        wrong = r["err"] > WRONG_DEG
        ok = ~np.isnan(r["hit"])
        full = r["full"] & ok
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
                    )
                )
            g_out[gname] = dict(
                n=int(m.sum()),
                n_wrong=int((m & wrong).sum()),
                # does not depend on the cut h: one value per group
                evals_in_wrong_pct=100.0
                * float(r["nev"][m & wrong].sum())
                / max(float(r["nev"][m].sum()), 1.0),
                hit_right_min=float(r["hit"][m & ~wrong].min()) if (m & ~wrong).any() else None,
                hit_wrong_max=float(r["hit"][m & wrong].max()) if (m & wrong).any() else None,
                hit_wrong_median=(
                    float(np.median(r["hit"][m & wrong])) if (m & wrong).any() else None
                ),
                hit_right_median=float(np.median(r["hit"][m & ~wrong])),
                rows=rows,
            )
        out[name] = dict(
            wrong_total=r["wrong_total"],
            groups=g_out,
            # from run.json stats: seeds + reseeds and their evaluations (the full searches)
            n_full_run=r["n_full_run"],
            evals_full_run=r["evals_full_run"],
            evals_per_full_search=int(round(r["evals_full_run"] / r["n_full_run"])),
        )
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
    for k, col in zip(range(LEVELS), ["tab:orange", "tab:green", "tab:purple", "tab:brown"]):
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
    # (3) saved vs broken: gross (no finisher; broken at 1 degree) solid, net (finisher charged)
    # dashed at the same broken count
    for k, col in zip(range(3), ["tab:orange", "tab:green", "tab:purple"]):
        rr = [r for r in grid if r["level"] == k]
        x = [r["broken_1deg"] for r in rr]
        ax[2].plot(x, [r["saved_gross_pct"] for r in rr], "o-", c=col, label=f"level {k}, gross")
        ax[2].plot(x, [r["saved_net_pct"] for r in rr], "s--", c=col, label=f"level {k}, net")
    n_right = int((~wrong).sum())
    ax[2].set_xlabel(f"voxels broken (of {n_right} right in baseline; candidate > 1 deg off)")
    ax[2].set_ylabel("raw evaluations saved (% of baseline total)")
    ax[2].set_title("Stop rule: best cost <= c (c swept 0.2 to 0.95)")
    ax[2].legend(fontsize=7, ncol=2)
    fig.tight_layout()
    fig.savefig(out / "fig_stopping.png", dpi=130)
    plt.close(fig)

    # pilot hit ratio
    fig, ax = plt.subplots(1, 2, figsize=(12, 4.4))
    r = P["keyoff_rerun"]
    wr = r["err"] > WRONG_DEG
    ok = ~np.isnan(r["hit"])
    full = r["full"] & ok
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


PILOT_CUT = 0.93


def zero_counts(summ: Dict[str, Any]) -> List[Dict[str, Any]]:
    """Wilson 95% intervals (scripts/common/stats.py) for the counts the page reports, in %."""
    import stats as ST

    rows: List[Dict[str, Any]] = []

    def add(name: str, k: int, n: int) -> None:
        lo, hi = ST.wilson(k, n)
        rows.append(dict(count=name, k=k, n=n, lo_pct=100.0 * lo, hi_pct=100.0 * hi))

    for rule, lab in [("old", "old"), ("cpp", "C++")]:
        s = summ[rule]
        mis = min(r["right_gt"] + r["wrong_le"] for r in s["final"]["rows"])
        add(f"E0 {lab} rule: voxels misclassified by the final cost (best c)", mis, s["n"])
    s = summ["old"]
    n_right = s["n"] - s["n_wrong"]
    for k, r in s["safe0"].items():
        add(
            f"E0 old rule: broken at 1 degree, safe rule after level {k}", r["broken_1deg"], n_right
        )
    c = summ["cpp"]
    for r in summ["cross_old_to_cpp"]:
        lv = r["level"]
        add(
            f"E0 C++ rule: broken at 1 degree, old-rule c after level {lv}",
            r["broken_1deg"],
            c["n"] - c["n_wrong"],
        )
    g = summ["pilot"]["keyoff_rerun"]["groups"]
    cut = [
        [x for x in g[nm]["rows"] if abs(x["h"] - PILOT_CUT) < 1e-9][0]
        for nm in ["full_search", "local_fit"]
    ]
    n_wrong = sum(g[nm]["n_wrong"] for nm in ["full_search", "local_fit"])
    add("pilot key-off: wrong voxels at or above the cut", sum(x["wrong_ge"] for x in cut), n_wrong)
    n_ge = sum(x["wrong_ge"] + x["right_ge"] for x in cut)
    add(
        "pilot key-off: answers at or above the cut that are wrong",
        sum(x["wrong_ge"] for x in cut),
        n_ge,
    )
    return rows


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    ap.add_argument("--eval-root", required=True, help="main checkout benchmarks/eval_reduction")
    ap.add_argument("--out", default=str(ROOT / "benchmarks" / "autoresearch" / "ar0"))
    a = ap.parse_args()
    root, out = Path(a.eval_root), Path(a.out)
    out.mkdir(parents=True, exist_ok=True)

    E = {rule: load_e0(root, rule) for rule in ["old", "cpp"]}
    for e in E.values():
        for k in range(LEVELS):
            assert np.allclose(e[f"L{k}_diameter_deg"], e[f"L{k}_diameter_deg"][0])
    summ: Dict[str, Any] = {
        "sign": "cost lower = better",
        "units": "raw evaluation counts; the level-0 screen counts interpolated evaluations "
        "(one per orientation, summed with strict ones), so the counts are not strict-equivalent",
        "wrong_deg": WRONG_DEG,
        "finish_deg": FINISH_DEG,
        "diameters_deg": [float(E["old"][f"L{k}_diameter_deg"][0]) for k in range(LEVELS)],
        "true_basin": "first screened candidate (R_disc: spacing-filter survivors ordered by "
        "screen score, before the quick MC) within the level's grid diameter of the truth",
        "safe_rule": "per level, at most max_broken broken at 1 degree; largest gross saving; "
        "on a tie, the smaller c",
    }
    for rule, e in E.items():
        wrong = e["err"] > WRONG_DEG
        d: Dict[str, Any] = dict(
            n=int(len(wrong)),
            n_wrong=int(wrong.sum()),
            evals_total=int(e["evals"].sum()),
            evals_per_voxel_mean=int(round(float(e["evals"].mean()))),
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
        d["grid"] = merge_grids(d["grid_tol1"], d["grid_tol2"])
        d["levels_sep"] = {f"L{k}": separation(e[f"L{k}_best_cost"], wrong) for k in range(LEVELS)}
        d["safe0"] = {str(k): v for k, v in best_safe(d["grid"], 0).items()}
        d["safe1"] = {str(k): v for k, v in best_safe(d["grid"], 1).items()}
        d["split"] = e0_split(e)
        summ[rule] = d
    # threshold picked on the old-rule voxels, applied to the C++-rule run of the SAME 200 voxels
    # and images (a different overlap rule, not new voxels)
    cross = []
    for k, r in summ["old"]["safe0"].items():
        g = [x for x in summ["cpp"]["grid"] if x["level"] == int(k) and x["c"] == r["c"]]
        assert len(g) == 1
        cross.append(g[0])
    summ["cross_old_to_cpp"] = cross

    # full-sample voxels: stage split from the six-voxel G2 ledger
    full = full_sample_split(root / "G" / "g2" / "ledger" / "runs")
    summ["full_sample"] = full
    summ["split_old"] = split_tables(summ["old"]["split"], full)

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
    summ["zero_counts"] = zero_counts(summ)
    figures(E, summ["old"]["grid"], P, out)
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
                "far_cost_min",
                "evals_after_pct",
                "evals_after_right_mean",
            ],
            {"far_cost_min": ".3f", "evals_after_pct": ".1f", "evals_after_right_mean": ".0f"},
        )
        T[f"stages_{rule}"] = DT.markdown_table(
            s["stages"],
            [
                "level",
                "diameter_deg",
                "screen_pct",
                "quick_pct",
                "n_disc_median",
                "right_with_basin",
                "right_no_basin",
                "wrong_with_basin",
                "basin_pos_median",
                "wrong_basin_pos_median",
                "basin_frac_median",
                "basin_frac_p90",
            ],
            {
                "diameter_deg": ".2f",
                "screen_pct": ".1f",
                "quick_pct": ".1f",
                "n_disc_median": ".0f",
                "basin_pos_median": ".0f",
                "wrong_basin_pos_median": ".0f",
                "basin_frac_median": ".2f",
                "basin_frac_p90": ".2f",
            },
        )
        sel = [r for r in s["grid"] if r["c"] in (0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)]
        T[f"stop_grid_{rule}"] = DT.markdown_table(
            sel,
            [
                "level",
                "c",
                "n_stop",
                "broken_1deg",
                "broken_2deg",
                "stop_wrong_base",
                "saved_gross_pct",
                "saved_net_pct",
            ],
            {"c": ".2f", "saved_gross_pct": ".2f", "saved_net_pct": ".2f"},
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
        for mb, key in [(0, "safe0"), (1, "safe1")]:
            for k, r in s[key].items():
                safe.append(dict(max_broken=mb, **r))
        T[f"safe_{rule}"] = DT.markdown_table(
            safe,
            [
                "max_broken",
                "level",
                "c",
                "n_stop",
                "broken_1deg",
                "broken_2deg",
                "saved_gross_pct",
                "saved_net_pct",
            ],
            {"c": ".2f", "saved_gross_pct": ".2f", "saved_net_pct": ".2f"},
        )
    T["cross_old_to_cpp"] = DT.markdown_table(
        cross,
        ["level", "c", "n_stop", "broken_1deg", "broken_2deg", "saved_gross_pct", "saved_net_pct"],
        {"c": ".2f", "saved_gross_pct": ".2f", "saved_net_pct": ".2f"},
    )
    sp = summ["split_old"]
    T["split_stages"] = DT.markdown_table(
        sp["stages"],
        ["stage", "e0_evals_pct", "full_evals_pct", "full_time_pct"],
        {"e0_evals_pct": ".1f", "full_evals_pct": ".1f", "full_time_pct": ".1f"},
    )
    T["split_after"] = DT.markdown_table(
        sp["after"],
        ["after_level", "e0_evals_pct", "full_evals_pct", "full_time_pct"],
        {"e0_evals_pct": ".1f", "full_evals_pct": ".1f", "full_time_pct": ".1f"},
    )
    T["zero_counts"] = DT.markdown_table(
        summ["zero_counts"],
        ["count", "k", "n", "lo_pct", "hi_pct"],
        {"lo_pct": ".1f", "hi_pct": ".1f"},
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
                    evals_in_wrong_pct=v["evals_in_wrong_pct"],
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
            "evals_in_wrong_pct",
        ],
        {
            **{
                k: ".3f"
                for k in ["hit_right_min", "hit_right_median", "hit_wrong_median", "hit_wrong_max"]
            },
            "evals_in_wrong_pct": ".1f",
        },
    )
    for gname in ["full_search", "local_fit"]:
        rows = summ["pilot"]["keyoff_rerun"]["groups"][gname]["rows"]
        T[f"pilot_thr_{gname}"] = DT.markdown_table(
            rows, ["h", "right_ge", "wrong_ge", "right_lt", "wrong_lt"], {"h": ".2f"}
        )
    DT.write_tables(out / "tables.md", T)
    print("wrote", out)


if __name__ == "__main__":
    main()
