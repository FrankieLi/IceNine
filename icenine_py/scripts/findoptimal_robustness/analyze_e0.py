#!/usr/bin/env python3
"""E0 analysis: where is the truth lost? Reads cache/e0/*.npz (e0_run.py), writes
benchmarks/findoptimal_robustness/e0_runs.csv-like npz, e0_summary.{txt,json} and the where-lost
plot.

Definitions (per run, error = cubic-symmetry-reduced misorientation to the truth; wrong > 1 deg):
  basin(t)   a candidate within t deg of the truth (t = 1 and 3 deg)
  present[L] a level-L candidate after the quick MC is in the basin
  kept[L]    a basin candidate is among the n_keep best (by post-MC local cost) of level L
  hand-off   the last level's sorted candidates (what FindOptimal receives)
  S1  wrong, no basin(3) candidate at level 0 after the quick MC          (sampling / discrete)
  S2  wrong, a basin(3) candidate exists at some level but is not kept     (ranking: pruned)
  S2b wrong, a basin candidate was kept at level L but level L+1 has none  (not regenerated)
  S3  wrong, a basin(3) candidate is in the hand-off, FindOptimal returns something else
"""

import json
import sys
from collections import Counter
from pathlib import Path
from typing import Any, Dict, List

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import common as C  # noqa: E402
import csl  # noqa: E402

N_LEVELS = 4
BASIN = (1.0, 3.0)


def run_metrics(d: Any, s: int, R_true: np.ndarray) -> Dict[str, Any]:
    p = f"s{s}_"
    final_err = float(C.err_deg(d[p + "R_final"], R_true))
    m: Dict[str, Any] = dict(final_err=final_err, wrong=final_err > C.WRONG_DEG)
    m["cost_final"] = float(d[p + "cost_final"])
    m["runtime"] = float(d[p + "runtime"])
    m["evals_global"] = int(d[p + "evals_global"])
    m["evals_local"] = int(d[p + "evals_local"])
    present = {t: [] for t in BASIN}
    kept = {t: [] for t in BASIN}
    minerr, minerr_pre, rank, ncand, nkeep = [], [], {t: [] for t in BASIN}, [], []
    rank_disc = {t: [] for t in BASIN}
    for L in range(N_LEVELS):
        if p + f"L{L}_qmc_R" not in d.files:
            for t in BASIN:
                present[t].append(False)
                kept[t].append(False)
                rank[t].append(-1)
                rank_disc[t].append(-1)
            minerr.append(np.nan)
            minerr_pre.append(np.nan)
            ncand.append(0)
            nkeep.append(0)
            continue
        R = d[p + f"L{L}_qmc_R"]
        cost = d[p + f"L{L}_qmc_cost"]
        perm = d[p + f"L{L}_qmc_perm"]
        nk = int(d[p + f"L{L}_qmc_n_keep"])
        e = C.err_deg(R, R_true)  # sorted order (post quick MC)
        pre = C.err_deg(d[p + f"L{L}_disc_R"], R_true)
        score = d[p + f"L{L}_disc_score"]
        disc_rank = np.empty(len(score), dtype=int)
        disc_rank[np.argsort(score, kind="stable")] = np.arange(len(score))
        minerr.append(float(e.min()))
        minerr_pre.append(float(pre.min()))
        ncand.append(len(e))
        nkeep.append(nk)
        for t in BASIN:
            idx = np.nonzero(e < t)[0]
            present[t].append(bool(len(idx)))
            kept[t].append(bool(len(idx) and idx.min() < nk))
            rank[t].append(int(idx.min()) if len(idx) else -1)
            rank_disc[t].append(int(disc_rank[perm[idx]].min()) if len(idx) else -1)
    m.update(
        present=present, kept=kept, min_err=minerr, min_err_pre=minerr_pre, rank=rank,
        rank_disc=rank_disc, n_cand=ncand, n_keep=nkeep,
    )  # fmt: skip
    ho = d[p + f"L{N_LEVELS - 1}_qmc_R"] if p + f"L{N_LEVELS - 1}_qmc_R" in d.files else None
    m["handoff_n"] = 0 if ho is None else len(ho)
    m["handoff_min_err"] = float("nan") if ho is None else float(C.err_deg(ho, R_true).min())
    # FindOptimal's per-candidate results
    fe = C.err_deg(d[p + "find_R_out"], R_true) if len(d[p + "find_R_out"]) else np.array([])
    m["find_err"] = fe.tolist()
    m["find_cost"] = d[p + "find_cost"].tolist()
    fin = d[p + "find_R_in"]
    m["find_in_err"] = C.err_deg(fin, R_true).tolist() if len(fin) else []
    # the stage the truth is lost in
    cls = "right"
    if m["wrong"]:
        t = 3.0
        if not present[t][0]:
            cls = "S1"
        elif m["handoff_n"] and m["handoff_min_err"] < t:
            cls = "S3"
        else:
            cls = "S2?"
            for L in range(N_LEVELS):
                if present[t][L] and not kept[t][L]:
                    cls = f"S2@{L}"
                    break
                if L + 1 < N_LEVELS and kept[t][L] and not present[t][L + 1]:
                    cls = f"S2b@{L}"
                    break
            if cls == "S2?":
                cls = "S2c"  # present and kept everywhere yet not in the hand-off (final level)
    m["class"] = cls
    m["csl"] = csl.csl_classify(R_true, d[p + "R_final"]) if m["wrong"] else None
    return m


def load_all(cache: Path, info: Dict[str, np.ndarray]) -> List[Dict[str, Any]]:
    rows: List[Dict[str, Any]] = []
    vox = [int(v) for v in info["voxel_indices"]]
    for vpos, v in enumerate(vox):
        for var in C.VARIANTS:
            f = cache / f"v{v}_{var}.npz"
            if not f.exists():
                continue
            d = np.load(f)
            if "unbuildable" in d.files:
                continue
            R_true = d["R_true"]
            for s in range(C.N_SEEDS):
                if f"s{s}_R_final" not in d.files:
                    continue
                m = run_metrics(d, s, R_true)
                m.update(vidx=v, vpos=vpos, variant=var, seed=s)
                m["r_perp"] = float(info["voxel_r_perp_um"][vpos])
                m["near_boundary"] = bool(info["voxel_near_boundary"][vpos])
                m["n_roi_true"] = int(info["voxel_n_roi_true"][vpos])
                rows.append(m)
    return rows


def main() -> None:
    info = dict(np.load(C.OUT_DIR / "voxels.npz"))
    rows = load_all(C.CACHE_DIR / "e0", info)
    # keep the first N_VOXELS buildable voxels
    built = sorted({r["vpos"] for r in rows})[: C.N_VOXELS]
    rows = [r for r in rows if r["vpos"] in set(built)]
    print(len(rows), "runs", len(built), "voxels")
    summ: Dict[str, Any] = dict(n_voxels=len(built))
    lines: List[str] = []
    for var in C.VARIANTS:
        R = [r for r in rows if r["variant"] == var]
        if not R:
            continue
        n = len(R)
        k = sum(r["wrong"] for r in R)
        p, lo, hi = C.wilson(k, n)
        s: Dict[str, Any] = dict(n_runs=n, n_wrong=int(k), wrong_rate=p, ci=[lo, hi])
        lines.append(f"[{var}] runs {n}  wrong {k} = {p:.3f} (95% CI {lo:.3f}-{hi:.3f})")
        per_seed = {}
        for sd in range(C.N_SEEDS):
            Rs = [r for r in R if r["seed"] == sd]
            ks = sum(r["wrong"] for r in Rs)
            pp = C.wilson(ks, len(Rs))
            per_seed[sd] = dict(n=len(Rs), wrong=int(ks), rate=pp[0], ci=[pp[1], pp[2]])
            lines.append(f"   seed {sd}: {ks}/{len(Rs)} = {pp[0]:.3f} ({pp[1]:.3f}-{pp[2]:.3f})")
        s["per_seed"] = per_seed
        # voxels wrong in all seeds vs some
        by_v: Dict[int, List[bool]] = {}
        for r in R:
            by_v.setdefault(r["vpos"], []).append(r["wrong"])
        full = [v for v, w in by_v.items() if len(w) == C.N_SEEDS]
        n_all = sum(all(by_v[v]) for v in full)
        n_some = sum(any(by_v[v]) and not all(by_v[v]) for v in full)
        n_none = sum(not any(by_v[v]) for v in full)
        s["voxels_all_seeds_wrong"], s["voxels_some_seeds_wrong"] = n_all, n_some
        s["voxels_never_wrong"], s["n_voxels_full"] = n_none, len(full)
        # expected under independence with the pooled rate
        pw = k / n
        s["indep_expected"] = dict(
            all=len(full) * pw**3, none=len(full) * (1 - pw) ** 3,
            some=len(full) * (1 - pw**3 - (1 - pw) ** 3),
        )  # fmt: skip
        lines.append(
            f"   voxels (all {C.N_SEEDS} seeds run, n={len(full)}): wrong in all seeds {n_all}, "
            f"in some {n_some}, never {n_none}; expected if seeds were independent with the pooled "
            f"rate: all {s['indep_expected']['all']:.1f}, some {s['indep_expected']['some']:.1f}, "
            f"never {s['indep_expected']['none']:.1f}"
        )
        # where lost
        wrongs = [r for r in R if r["wrong"]]
        cls = Counter(r["class"].split("@")[0] for r in wrongs)
        cls_full = Counter(r["class"] for r in wrongs)
        s["where_lost"] = dict(cls)
        s["where_lost_detail"] = dict(cls_full)
        lines.append(
            f"   where lost (wrong runs {len(wrongs)}): {dict(cls)}  detail {dict(cls_full)}"
        )
        # S2: would the discrete-stage score have kept it
        s2 = [r for r in wrongs if r["class"].startswith("S2")]
        if s2:
            ranks = []
            for r in s2:
                L = int(r["class"].split("@")[1]) if "@" in r["class"] else 0
                ranks.append(
                    (r["rank"][3.0][L], r["n_keep"][L], r["n_cand"][L], r["rank_disc"][3.0][L])
                )
            ranks_a = np.array(ranks)
            would_keep_disc = float(np.mean(ranks_a[:, 3] < ranks_a[:, 1]))
            s["S2_rank_summary"] = dict(
                median_rank=float(np.median(ranks_a[:, 0])),
                median_n_keep=float(np.median(ranks_a[:, 1])),
                median_n_cand=float(np.median(ranks_a[:, 2])),
                frac_rank_lt_2x_keep=float(np.mean(ranks_a[:, 0] < 2 * ranks_a[:, 1])),
                frac_kept_if_ranked_by_discrete_score=would_keep_disc,
            )
            lines.append(f"   S2 rank of best basin(3) candidate: {s['S2_rank_summary']}")
        # S3 detail
        s3 = [r for r in wrongs if r["class"] == "S3"]
        if s3:
            frac_best_in = float(
                np.mean([r["find_in_err"] and min(r["find_in_err"]) < 3 for r in s3])
            )
            fo_basin_cost_higher = np.nanmean(
                [
                    (
                        min(c for c, e in zip(r["find_cost"], r["find_err"]) if e < 1.0)
                        > r["cost_final"]
                        if any(e < 1.0 for e in r["find_err"])
                        else np.nan
                    )
                    for r in s3
                ]
            )
            lost_by_findoptimal = sum(any(e < 1.0 for e in r["find_err"]) for r in s3)
            s["S3_summary"] = dict(
                n=len(s3), frac_basin_in_input=frac_best_in,
                n_findoptimal_refined_basin_cand_below_1deg=int(lost_by_findoptimal),
                frac_basin_result_cost_above_winner_cost=float(fo_basin_cost_higher),
            )  # fmt: skip
            lines.append(f"   S3 detail: {s['S3_summary']}")
        # CSL classes
        cc = Counter((f"Sigma{r['csl']['label']}" if r["csl"]["sigma"] else "none") for r in wrongs)
        s["wrong_csl"] = dict(cc)
        lines.append(f"   CSL class of wrong answers: {dict(cc)}")
        # error of right answers
        right = [r["final_err"] for r in R if not r["wrong"]]
        s["right_median_err"] = float(np.median(right))
        s["wrong_final_cost_above_truth_cost"] = None
        # level progress: fraction of runs with basin(3) present / kept per level
        s["present3_by_level"] = [
            float(np.mean([r["present"][3.0][L] for r in R])) for L in range(N_LEVELS)
        ]
        s["kept3_by_level"] = [
            float(np.mean([r["kept"][3.0][L] for r in R])) for L in range(N_LEVELS)
        ]
        s["median_n_cand_by_level"] = [
            float(np.median([r["n_cand"][L] for r in R])) for L in range(N_LEVELS)
        ]
        lines.append(
            f"   basin(3) present by level {np.round(s['present3_by_level'], 3).tolist()}, "
            f"kept {np.round(s['kept3_by_level'], 3).tolist()}, median #candidates "
            f"{s['median_n_cand_by_level']}"
        )
        lines.append(
            f"   right answers: median error {s['right_median_err']:.4f} deg; runtime median "
            f"{np.median([r['runtime'] for r in R]):.1f} s; evals global median "
            f"{np.median([r['evals_global'] for r in R]):.0f}, local {np.median([r['evals_local'] for r in R]):.0f}"
        )
        # dependences
        for name, key in (("r_perp", "r_perp"),):
            rp = np.array([r[key] for r in R])
            w = np.array([r["wrong"] for r in R])
            q = np.quantile(rp, [0, 0.25, 0.5, 0.75, 1.0])
            tab = []
            for a, b in zip(q[:-1], q[1:]):
                sel = (rp >= a) & (rp <= b)
                tab.append((float(a), float(b), int(sel.sum()), float(w[sel].mean())))
            s["wrong_by_r_perp_quartile"] = tab
            lines.append(
                "   wrong rate by r_perp quartile (um): "
                + ", ".join(f"{a:.0f}-{b:.0f}: {r_:.2f}" for a, b, n_, r_ in tab)
            )
        for flag in (True, False):
            sel = [r for r in R if r["near_boundary"] == flag]
            if sel:
                s[f"wrong_near_boundary_{flag}"] = [
                    len(sel),
                    float(np.mean([r["wrong"] for r in sel])),
                ]
                lines.append(
                    f"   near grain boundary = {flag}: n={len(sel)}, wrong {np.mean([r['wrong'] for r in sel]):.3f}"
                )
        summ[var] = s
    (C.OUT_DIR / "e0_summary.txt").write_text("\n".join(lines) + "\n")
    C.save_json(C.OUT_DIR / "e0_summary.json", summ)
    print("\n".join(lines))
    # compact per-run table for later stages
    keys = [
        "vidx",
        "vpos",
        "seed",
        "final_err",
        "cost_final",
        "runtime",
        "evals_global",
        "evals_local",
    ]
    tab = {k: np.array([r[k] for r in rows]) for k in keys}
    tab["variant"] = np.array([C.VARIANTS.index(r["variant"]) for r in rows])
    tab["class"] = np.array([r["class"] for r in rows])
    tab["csl_label"] = np.array([r["csl"]["label"] if r["csl"] else "" for r in rows])
    tab["r_perp"] = np.array([r["r_perp"] for r in rows])
    tab["near_boundary"] = np.array([r["near_boundary"] for r in rows])
    np.savez_compressed(C.OUT_DIR / "e0_runs.npz", **tab)


if __name__ == "__main__":
    main()
