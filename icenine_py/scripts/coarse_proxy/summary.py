#!/usr/bin/env python3
"""Summary of the Q_max-8 cost proxy study: q-level check, label validation, offline metrics,
decision points D1-D3, end-to-end rows paired with the existing baseline / F1 / F1b / E2 rows,
cost accounting and every success criterion. Anything not run prints "not evaluated".
Writes benchmarks/coarse_proxy/summary.{txt,json}.

  uv run python scripts/coarse_proxy/summary.py
"""

import argparse
import json
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import base as B  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "common"))
import stats as shared_stats  # noqa: E402

C, M = B.C, B.M
import summarize_fixes as SF  # noqa: E402
import models as MD  # noqa: E402

NE = "not evaluated"


def jload(p: Path) -> Optional[Dict[str, Any]]:
    return json.loads(p.read_text()) if p.exists() else None


# -- per-run loaders ----------------------------------------------------------------------------
def runs_e0(var: str, seeds: List[int]) -> Dict[Any, Dict[str, Any]]:
    out = {}
    for vpos, v in enumerate(B.voxels()):
        f = C.CACHE_DIR / "e0" / f"v{v}_{var}.npz"
        if not f.exists() or "unbuildable" in np.load(f).files:
            continue
        d = np.load(f)
        for s in seeds:
            out[(vpos, s)] = dict(
                R=d[f"s{s}_R_final"], R_true=d["R_true"], cost=float(d[f"s{s}_cost_final"]),
                eg=int(d[f"s{s}_evals_global"]), el=int(d[f"s{s}_evals_local"]),
                rt=float(d[f"s{s}_runtime"]), n_scored=0, v=v,
            )  # fmt: skip
    return out


def runs_files(root: Path, var: str, seeds: List[int], seeded_names: bool) -> Dict[Any, Any]:
    out = {}
    for vpos, v in enumerate(B.voxels()):
        for s in seeds:
            name = f"v{v}_{var}.npz" if (s == 0 or not seeded_names) else f"v{v}_{var}_s{s}.npz"
            f = root / name
            if not f.exists():
                continue
            d = np.load(f)
            out[(vpos, s)] = dict(
                R=d["R_final"], R_true=d["R_true"], cost=float(d["cost_final"]),
                eg=int(d["evals_global"]), el=int(d["evals_local"]), rt=float(d["runtime"]),
                n_scored=int(d["n_scored"]) if "n_scored" in d.files else 0,
                proxy_seconds=float(d["proxy_seconds"]) if "proxy_seconds" in d.files else 0.0,
                v=v,
            )  # fmt: skip
    return out


def runs_f1b(var: str, seeds: List[int]) -> Dict[Any, Any]:
    out = {}
    for vpos, v in enumerate(B.voxels()):
        f = C.CACHE_DIR / "fix" / "F1b" / f"v{v}_{var}_s0.npz"
        if f.exists() and 0 in seeds:
            d = np.load(f)
            out[(vpos, 0)] = dict(
                R=d["s0_R_final"], cost=float(d["s0_cost_final"]), eg=int(d["s0_evals_global"]),
                el=int(d["s0_evals_local"]), rt=float(d["s0_runtime"]), n_scored=0, v=v,
            )  # fmt: skip
    return out


def with_f1(base: Dict[Any, Any], f1_dir: Path, var: str, e0: Dict[Any, Any]) -> Dict[Any, Any]:
    """F1 (Sigma <= 29) on top of `base` runs: answer = lowest local cost among the run's answer and
    the refined relatives; extra evaluations as in summarize_fixes (time scaled as there)."""
    out = {}
    for (vpos, s), r in base.items():
        # baseline F1 cache holds s0..s2 in one file; a proxy run's F1 holds seed 0 in v*.npz and
        # seeds 1-2 in v*_s12.npz (f1_run.py --source-seeds 1 2)
        name = f"v{r['v']}_{var}" + ("_s12" if (s > 0 and f1_dir.name != "f1") else "")
        f = f1_dir / f"{name}.npz"
        if not f.exists():
            continue
        R, c, extra = SF.f1_answer(
            np.load(f), s, 29, dict(R=r["R"], cost=r["cost"]), e0[(vpos, s)]["R_true"]
        )
        tot = r["eg"] + r["el"]
        out[(vpos, s)] = dict(
            r, R=R, cost=c, el=r["el"] + extra, rt=r["rt"] + extra * r["rt"] / max(tot, 1)
        )
    return out


def row(name: str, runs: Dict[Any, Any], e0: Dict[Any, Any], eq: float = 0.0) -> Optional[Dict]:
    """Row stats paired with the E0 baseline on the same (voxel, seed) keys. `eq` = evaluation
    equivalents of one proxy-scored candidate (0 = report the reconstructor's count only)."""
    keys = sorted(k for k in runs if k in e0)
    if not keys:
        return None
    err = np.array([float(C.err_deg(runs[k]["R"], e0[k]["R_true"])) for k in keys])
    base = np.array([float(C.err_deg(e0[k]["R"], e0[k]["R_true"])) for k in keys])
    n = len(keys)
    p, lo, hi = C.wilson(int((err > 1).sum()), n)
    right = err[err <= 1]
    g = np.array([runs[k]["eg"] for k in keys], float)
    lc = np.array([runs[k]["el"] for k in keys], float)
    sc = np.array([runs[k]["n_scored"] for k in keys], float)
    return dict(
        _errs=dict(zip(keys, err.tolist())),
        name=name, n=n, wrong=int((err > 1).sum()), rate=p, lo=lo, hi=hi,
        fixed=int(((err <= 1) & (base > 1)).sum()), broken=int(((err > 1) & (base <= 1)).sum()),
        base_wrong=int((base > 1).sum()), median_err_right=float(np.median(right)),
        evals_global=float(g.mean()), evals_local=float(lc.mean()), scored=float(sc.mean()),
        evals_total=float((g + lc).mean()), evals_with_proxy_eq=float((g + lc + eq * sc).mean()),
        base_evals=float(np.mean([e0[k]["eg"] + e0[k]["el"] for k in keys])),
        wall=float(np.mean([runs[k]["rt"] for k in keys])),
        base_wall=float(np.mean([e0[k]["rt"] for k in keys])),
        proxy_sec=float(np.mean([runs[k].get("proxy_seconds", 0.0) for k in keys])),
    )  # fmt: skip


def mcnemar(a: Optional[Dict], b: Optional[Dict]) -> Optional[Dict[str, Any]]:
    """Exact (two-sided binomial) McNemar test of two rows on their common runs: wrong = error
    above 1 deg. Returns the discordant counts (a wrong / b right, a right / b wrong) and p."""
    if not a or not b:
        return None
    keys = sorted(set(a["_errs"]) & set(b["_errs"]))
    wa = np.array([a["_errs"][k] > C.WRONG_DEG for k in keys])
    wb = np.array([b["_errs"][k] > C.WRONG_DEG for k in keys])
    n_a, n_b = shared_stats.paired_discordant(wa, wb)
    p = shared_stats.mcnemar_exact(n_a, n_b)
    return dict(a=a["name"], b=b["name"], n=len(keys), a_only_wrong=n_a, b_only_wrong=n_b, p=p)


def fmt_mc(m: Optional[Dict[str, Any]]) -> str:
    if not m:
        return NE
    return (
        f"{m['a']} vs {m['b']}: discordant {m['a_only_wrong']} (only the first wrong) vs "
        f"{m['b_only_wrong']} (only the second wrong), exact McNemar p = {m['p']:.3g} (n={m['n']})"
    )


def fmt_row(r: Dict[str, Any]) -> str:
    return (
        f"{r['name']:34s} n={r['n']:3d} wrong {r['wrong']:3d} = {r['rate']:.3f} "
        f"[{r['lo']:.3f},{r['hi']:.3f}]  fixed {r['fixed']:3d} broken {r['broken']:2d}  "
        f"med.err(right) {r['median_err_right']:.4f}  evals global {r['evals_global']:.0f} "
        f"local {r['evals_local']:.0f} (total {r['evals_total']:.0f}, incl. proxy (1-worker eq) "
        f"{r['evals_with_proxy_eq']:.0f}, base {r['base_evals']:.0f}) scored {r['scored']:.0f}  "
        f"wall {r['wall']:.1f}s (base {r['base_wall']:.1f}s, proxy {r['proxy_sec']:.2f}s)"
    )


def _tags() -> List[str]:
    root = B.CACHE / "e2e"
    names = [p.name for p in root.glob("*") if p.is_dir()] if root.exists() else []
    return sorted(n for n in names if "harvest" not in n and not n.startswith("f1_"))


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--rows", nargs="*", default=None, help="proxy row tags (default: all found)")
    a = ap.parse_args()
    out: Dict[str, Any] = {}
    L: List[str] = ["=== Q_max-8 cost proxy: summary (caveat: per-voxel images, <= 3 distractor "
                    "sources; 10-worker wall times are noisy) ==="]  # fmt: skip
    # 1. q levels
    ql = jload(B.OUT / "qlevels.json")
    L.append("\n[1] |q| families (Cu) and the cut-offs")
    if ql:
        L.append(
            f"  q_levels {ql['q_levels']}  reflections per family {ql['reflections_per_family']}"
        )
        for q, v in ql["per_cutoff"].items():
            L.append(f"  Q<={q}: families {v['families']}  ({v['n_reflections']} reflections)")
        L.append("  Q3 contains no reflection and is dropped; Q4 = {111},{200}; Q5 adds {220}.")
    else:
        L.append("  " + NE)
    out["qlevels"] = ql
    # 2. labels
    lv = jload(B.OUT / "label_validation.json")
    L.append("\n[2] y_bcost label validation (censored branch, refined candidates)")
    if lv:
        L.append(
            f"  n={lv['n']}  Spearman(y_bcost, refined cost) pooled {lv['spearman_pooled']:.3f}, "
            f"mean within case {lv['spearman_within_case_mean']:.3f}  -> {lv['verdict']}"
        )
        L.append(f"  by level {lv['by_level']}")
    else:
        L.append("  " + NE)
    out["label_validation"] = lv
    # 3. offline
    off = jload(B.OUT / "offline_metrics.json")
    L.append("\n[3] offline metrics (grain-disjoint held-out folds; trained on clean+realistic)")
    if off:
        for tv in ("both", "clean", "all"):
            L.append(f"  -- test on {'realistic' if tv == 'all' else tv}: see offline_metrics.txt")
        L.append(open(B.OUT / "offline_metrics.txt").read().rstrip())
        d = off["decisions"]
        L.append("\n[4] decision points (pooled held-out, recall L0-2 at keep 1/4 unless noted)")
        d1 = d["D1"]
        L.append(
            f"  D1: F-cache Q5 clf3 recall {d1['cache5_clf3_recall']:.3f} "
            f"(reg {d1['cache5_reg_recall']:.3f}; E2 full {d1['e2full_clf3_recall']:.3f}) "
            f"-> low-Q signal mostly lost: "
            f"{d1['low_q_signal_mostly_lost']}. Best deployable {d1['best_deployable']} "
            f"{d1['best_deployable_recall']:.3f}, gap to E2 {d1['e2_gap']:.3f} -> adds nothing "
            f"beyond E2 (gap > 0.01): {d1['adds_nothing_beyond_e2']}"
        )
        L.append(
            f"  D2: {d['D2']['model']} recall at keep 1/8 = {d['D2']['recall_1_8']:.3f} "
            f"(>= 0.98: {d['D2']['holds']})"
        )
        for fs, v in d["D3"].items():
            L.append(
                f"  D3 {fs}: reg {v['reg']:.3f} vs best classifier {v['best_clf']:.3f} "
                f"-> carry {v['carry']}"
            )
        out["offline_decisions"] = d
    else:
        L.append("  " + NE)
    # 5. timing
    tm = jload(B.OUT / "timing.json")
    L.append(
        "\n[5] cost accounting, one worker, median of %s candidates"
        % (tm or {}).get("n_candidates", "?")
    )
    eq_proxy = eq_e2 = 0.0
    if tm:
        s_, e_ = tm["seconds_median"], tm["equivalents_of_local_q8"]
        for k in s_:
            L.append(f"  {k:30s} {s_[k] * 1e3:9.3f} ms  = {e_[k]:7.2f} local Q8 evaluations")
        eq_proxy = e_["lowq5_full"] + e_["model_predict_per_candidate"]
        eq_e2 = e_["e2_full"] + e_["model_predict_per_candidate"]
        L.append(
            f"  per scored candidate: proxy (lowq5 pass + prediction) {eq_proxy:.2f} eq; "
            f"E2 (full pass + "
            f"prediction) {eq_e2:.2f} eq (model timed: {tm['set']}|{tm['target']} for both)"
        )
    else:
        L.append("  " + NE)
    out["timing"] = tm
    # 6. end to end
    tags = a.rows
    if tags is None:
        tags = _tags()
    L.append("\n[6] end to end, 200 voxels x 2 variants, seed 0 (paired with E0 seed 0)")
    e2e: Dict[str, Any] = {}
    for var in C.VARIANTS:
        e0 = runs_e0(var, [0, 1, 2])
        e0s0 = {k: v for k, v in e0.items() if k[1] == 0}
        rows: List[Dict[str, Any]] = []
        nm = "realistic" if var == "all" else "clean"
        L.append(f"--- {nm} ---")
        rows.append(row("E0 baseline", e0s0, e0))
        base_f1 = with_f1(e0s0, C.CACHE_DIR / "f1", var, e0)
        rows.append(row("F1 (Sigma<=29)", base_f1, e0))
        rows.append(row("F1b", runs_f1b(var, [0]), e0))
        e2r = runs_files(C.CACHE_DIR / "e2_rerank_GBT", var, [0], False)
        rows.append(row("E2 rerank (GBT, full pass)", e2r, e0, eq_e2))
        rows.append(
            row(
                "E2 rerank + F1", with_f1(e2r, C.CACHE_DIR / "f1_e2_rerank_GBT", var, e0), e0, eq_e2
            )
        )
        for tag in tags:
            pr = runs_files(B.CACHE / "e2e" / tag, var, [0], True)
            rows.append(row(f"proxy {tag}", pr, e0, eq_proxy))
            f1d = B.CACHE / "e2e" / f"f1_{tag}"
            if f1d.exists():
                rows.append(row(f"proxy {tag} + F1", with_f1(pr, f1d, var, e0), e0, eq_proxy))
        rows = [r for r in rows if r]
        for r in rows:
            L.append("  " + fmt_row(r))
        b0 = rows[0]
        L.append("  per-run change vs the baseline (mean over the paired runs):")
        for r in rows[1:]:
            if r["scored"] == 0 and "rerank" not in r["name"] and "proxy" not in r["name"]:
                continue
            eq = r["evals_with_proxy_eq"] - r["evals_total"]
            L.append(
                f"    {r['name']:30s} d global {r['evals_global'] - b0['evals_global']:+7.0f}  "
                f"d local {r['evals_local'] - b0['evals_local']:+7.0f}  scoring (equivalents) "
                f"{eq:+6.0f}  d total {r['evals_with_proxy_eq'] - b0['evals_total']:+7.0f} "
                f"({(r['evals_with_proxy_eq'] / b0['evals_total'] - 1) * 100:+.1f}%)  "
                f"d wall (10 workers) {r['wall'] - b0['wall']:+.1f}s"
            )
        by_name = {r["name"]: r for r in rows}
        tests = [
            mcnemar(by_name.get("proxy p_i"), by_name.get("E2 rerank (GBT, full pass)")),
            mcnemar(by_name.get("proxy p_i + F1"), by_name.get("F1 (Sigma<=29)")),
            mcnemar(by_name.get("proxy p_i + F1"), by_name.get("E2 rerank + F1")),
            mcnemar(by_name.get("proxy p_i + F1"), by_name.get("F1b")),
        ]
        out.setdefault("paired_tests", {})[var] = [t for t in tests if t]
        L.append("  paired exact McNemar tests (seed 0, same runs):")
        L += ["    " + fmt_mc(t) for t in tests if t]
        for r in rows:
            r.pop("_errs", None)
        e2e[var] = rows
    out["end_to_end"] = e2e
    # 7. success criteria
    L.append("\n[7] success criteria")
    L += criteria(out, e2e, tm)
    # 8. seeds 1-2 of the proxy rows
    ms: Dict[str, Any] = {}
    L.append(
        "\n[8] proxy rows over seeds 0-2 (paired with E0 seeds 0-2; F1b exists for seed 0 only)"
    )
    for tag in tags:
        for var in C.VARIANTS:
            e0 = runs_e0(var, [0, 1, 2])
            pr = runs_files(B.CACHE / "e2e" / tag, var, [0, 1, 2], True)
            if not any(k[1] > 0 for k in pr):
                continue
            rows8 = [
                row("E0 baseline (3 seeds)", e0, e0),
                row("F1 (Sigma<=29, 3 seeds)", with_f1(e0, C.CACHE_DIR / "f1", var, e0), e0),
                row(f"proxy {tag} (3 seeds)", pr, e0, eq_proxy),
            ]
            f1d = B.CACHE / "e2e" / f"f1_{tag}"
            if f1d.exists():
                rows8.append(
                    row(f"proxy {tag} + F1 (3 seeds)", with_f1(pr, f1d, var, e0), e0, eq_proxy)
                )
            rows8 = [r for r in rows8 if r]
            L.append(f"  -- {'realistic' if var == 'all' else 'clean'} --")
            L += ["  " + fmt_row(r) for r in rows8]
            nm = {r["name"]: r for r in rows8}
            mc = mcnemar(nm.get(f"proxy {tag} + F1 (3 seeds)"), nm.get("F1 (Sigma<=29, 3 seeds)"))
            if mc:
                L.append("    " + fmt_mc(mc))
            ms[f"{tag}/{var}"] = dict(
                rows=[{k: v for k, v in r.items() if k != "_errs"} for r in rows8], mcnemar=mc
            )
    if not ms:
        L.append("  " + NE)
    out["multi_seed"] = ms
    # 8b. where the extra local evaluations go (scripts/coarse_proxy/eval_split.py)
    L.append("\n[8b] evaluation split, baseline vs proxy row (i) (20 voxels x 2 variants, seed 0)")
    sp = jload(B.OUT / "eval_split.json")
    if sp:
        L += sp["lines"]
    else:
        L.append("  " + NE)
    out["eval_split"] = sp
    # 9. domain shift
    L.append("\n[9] domain shift: recall on the proxy runs' own candidates")
    ds = domain_shift(tags)
    L += ds["lines"] or ["  " + NE]
    out["domain_shift"] = ds["data"]
    (B.OUT / "summary.txt").write_text("\n".join(L) + "\n")
    C.save_json(B.OUT / "summary.json", out)
    print("\n".join(L))


def verdict(ok: bool) -> str:
    return "met" if ok else "NOT met"


def criteria(out: Dict[str, Any], e2e: Dict[str, Any], tm: Optional[Dict[str, Any]]) -> List[str]:
    L: List[str] = []
    lv = out.get("label_validation")
    if lv:
        v = (
            f"pooled {verdict(lv['spearman_pooled'] >= 0.8)} ({lv['spearman_pooled']:.3f}); "
            f"within-case mean {lv['spearman_within_case_mean']:.3f} is borderline BELOW 0.8 "
            f"(L2 {lv['by_level']['2']['spearman']:.3f}); ranking happens within a case and "
            f"level, but the classifiers tie the regression offline (D3), so a fallback would "
            f"not change the offline recall"
        )
    else:
        v = NE
    L.append(f"  label: Spearman of the censored y_bcost >= 0.8 -> {v}")
    d = (out.get("offline_decisions") or {}).get("D1")
    off = jload(B.OUT / "offline_metrics.json")
    if d and off:
        k, rec = "recall_L012_1/4", {}
        for tv in ("both", "clean", "all"):
            best = max((m for m in off[tv] if m.startswith("e2full")), key=lambda m: off[tv][m][k])
            rec[tv] = (off[tv][d["best_deployable"]][k], best, off[tv][best][k])
        v = (
            f"{verdict(not d['adds_nothing_beyond_e2'])} against e2full|clf3 "
            f"(gap {d['e2_gap']:.3f}); against the best E2 variant the gap is "
            f"{rec['both'][2] - rec['both'][0]:.3f} pooled "
            f"({rec['both'][1]} {rec['both'][2]:.3f}), "
            f"{rec['clean'][2] - rec['clean'][0]:.3f} clean, "
            f"{rec['all'][2] - rec['all'][0]:.3f} realistic; the deployable model was picked as "
            f"the best of 6 held-out results (mildly optimistic)"
        )
    else:
        v = NE
    L.append(
        f"  informational D1 sub-check (NOT a plan success criterion; borderline): best deployable "
        f"proxy recall within 0.01 of E2 -> {v}"
    )
    by = {var: {r["name"]: r for r in rows} for var, rows in e2e.items()}
    found = False
    for var in C.VARIANTS:
        nm = "realistic" if var == "all" else "clean"
        b = by[var]
        e2, f1, f1b, base = (
            b.get("E2 rerank (GBT, full pass)"), b.get("F1 (Sigma<=29)"), b.get("F1b"),
            b.get("E0 baseline"),
        )  # fmt: skip
        for name, r in b.items():
            if not name.startswith("proxy "):
                continue
            found = True
            if name.endswith("+ F1"):
                if f1 and f1b:
                    ok = r["wrong"] <= f1["wrong"] and r["evals_with_proxy_eq"] < f1b["evals_total"]
                    mcs = " ; ".join(
                        fmt_mc(t)
                        for t in out.get("paired_tests", {}).get(var, [])
                        if t["a"] == name and t["b"] in ("F1 (Sigma<=29)", "E2 rerank + F1")
                    )
                    L.append(f"  [{nm}] CAVEAT (iii): wrong counts are small; paired tests: {mcs}")
                    L.append(
                        f"  [{nm}] (iii) {name}: wrong {r['wrong']} <= F1 {f1['wrong']} and evals "
                        f"{r['evals_with_proxy_eq']:.0f} < F1b {f1b['evals_total']:.0f} -> "
                        f"{verdict(ok)}"
                    )
                continue
            if e2:
                msg = ""
                if tm:
                    pe = tm["seconds_median"]
                    t_p = r["scored"] * (pe["lowq5_full"] + pe["model_predict_per_candidate"])
                    t_e = e2["scored"] * (pe["e2_full"] + pe["model_predict_per_candidate"])
                    t_d = e2["scored"] * (
                        pe["e2_full"] - pe["local_q8"] + pe["model_predict_per_candidate"]
                    )
                    msg = (
                        f"; proxy time per run {t_p:.2f}s ({r['scored']:.0f} scored) vs E2 "
                        f"{t_e:.2f}s as timed ({e2['scored']:.0f} scored) -> {verdict(t_p < t_e)}; "
                        f"vs a deployable E2 without the free Q8 cost {t_d:.2f}s -> "
                        f"{verdict(t_p < t_d)}"
                    )
                L.append(
                    f"  [{nm}] rerank {name}: wrong {r['wrong']}/{r['n']} <= E2 rerank "
                    f"{e2['wrong']}/{e2['n']} -> {verdict(r['wrong'] <= e2['wrong'])}{msg}"
                )
            if base:
                ok = r["rate"] <= 0.05 and r["evals_with_proxy_eq"] <= base["evals_total"]
                L.append(
                    f"  [{nm}] (ii)-criterion {name}: wrong {r['rate']:.3f} <= 0.05 and evals "
                    f"incl. proxy {r['evals_with_proxy_eq']:.0f} <= baseline "
                    f"{base['evals_total']:.0f} -> {verdict(ok)}"
                )
    if not found:
        L.append("  end-to-end criteria -> " + NE)
    return L


def offline_seed0(name: str) -> Dict[str, Dict[str, float]]:
    """Held-out recall (L0-2, keep 1/4 and 1/8) of the offline scores of model `name` on the
    baseline (E0 seed 0) candidates only, per variant: the reference of the domain-shift check."""
    sp = B.CACHE / "scores.npz"
    if not sp.exists():
        return {}
    D = MD.load_all()
    S = {name: np.load(sp)[name.replace("|", "__")]}
    out = {}
    for var in C.VARIANTS:
        mask = (D["variant"] == C.VARIANTS.index(var)) & ((D["seed"] == 0) | (D["source"] >= 2))
        r = MD.evaluate(D, S, mask)[name]
        out[var] = {"1/4": r["recall_L012_1/4"], "1/8": r["recall_L012_1/8"], "n": r["n_L012"]}
    return out


def domain_shift(tags: List[str]) -> Dict[str, Any]:
    lines: List[str] = []
    data: Dict[str, Any] = {}
    ref = offline_seed0("lowq5+c8|reg") if tags else {}
    for tag in tags:
        d = B.CACHE / "e2e" / f"{tag}_harvest"
        if not d.exists():
            continue
        files = sorted(d.glob("*.npz"))
        for var in C.VARIANTS:
            ok = {0.25: [0, 0], 0.125: [0, 0]}
            for fpath in files:
                z = np.load(fpath)
                if int(z["variant"]) != C.VARIANTS.index(var) or int(z["seed"]) != 0:
                    continue
                for lvl in range(3):
                    m = z["level"] == lvl
                    if not m.any() or not (z["err"][m] < 3).any():
                        continue
                    s, e = z["score"][m], z["err"][m]
                    rk = np.empty(m.sum(), int)
                    rk[np.argsort(-s, kind="stable")] = np.arange(m.sum())
                    for fr in ok:
                        ok[fr][0] += int(rk[e < 3].min() < max(1, int(m.sum() * fr)))
                        ok[fr][1] += 1
            if ok[0.25][1]:
                data[f"{tag}/{var}"] = {str(k): v[0] / v[1] for k, v in ok.items()} | dict(
                    n=ok[0.25][1]
                )
                if var in ref:
                    data[f"{tag}/{var}"]["offline_baseline_seed0"] = ref[var]
                    lines.append(
                        f"  {tag} {var} offline reference (baseline candidates, seed 0, held-out): "
                        f"keep 1/4 {ref[var]['1/4']:.3f}, 1/8 {ref[var]['1/8']:.3f} "
                        f"(n={ref[var]['n']})"
                    )
                lines.append(
                    f"  {tag} {var}: recall L0-2 on the run's own candidates keep 1/4 "
                    f"{ok[0.25][0] / ok[0.25][1]:.3f}, 1/8 {ok[0.125][0] / ok[0.125][1]:.3f} "
                    f"(n={ok[0.25][1]} level groups)"
                )
    return dict(lines=lines, data=data)


if __name__ == "__main__":
    main()
