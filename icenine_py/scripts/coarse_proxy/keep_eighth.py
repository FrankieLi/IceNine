#!/usr/bin/env python3
"""T4: proxy rerank (symmetry-agnostic) with keep_fraction 1/8 and 1/6, end to end.

summary : rows paired with the stored E0 baseline and the stored keep-1/4 proxy row (p_i), wrong
          rate with Wilson CI, exact McNemar, median error of right answers, evaluations
          (reconstructor count and + proxy in evaluation equivalents at the BATCHED cost of the
          batch size of every rank_key call), the "useful" verdict.
          Writes benchmarks/coarse_proxy/keep_eighth.{txt,json} and keep_eighth_tables.md.
timing  : single-worker timing of a subset, arms interleaved per case, preflight gate.
          Writes keep_eighth_timing.json, keep_eighth_preflight.json and the timing table.

uv run python scripts/coarse_proxy/keep_eighth.py summary
uv run python scripts/coarse_proxy/keep_eighth.py timing --arms base p_i:0.25 p_k8:0.125
"""

import argparse
import json
import sys
import time
from pathlib import Path
from typing import Any, Dict, List, Tuple

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import base as B  # noqa: E402

C, M = B.C, B.M
import summary as S  # noqa: E402

sys.path.insert(0, str(B.ICENINE_PY / "scripts" / "common"))
import doc_tables  # noqa: E402
import stats as shared_stats  # noqa: E402

# rows: (tag, label). p_i is the stored keep-1/4 proxy (per-candidate pass, same R_final as p_k4b)
ROWS = [("p_i", "proxy keep 1/4"), ("p_k6", "proxy keep 1/6"), ("p_k8", "proxy keep 1/8")]
SIZE_SOURCE = {"p_i": "p_k4b"}  # rank_key call sizes of the batched re-run of the same runs


def eq_curve() -> Tuple[np.ndarray, np.ndarray]:
    """Proxy cost (evaluation equivalents per candidate, scoring + GBT) against the batch size,
    from the T3 single-worker timing (batch 1, 50, 200; log-linear interpolation, clipped)."""
    t = json.loads((B.OUT / "timing_batch.json").read_text())["proxy_equivalents"]
    return np.array([1.0, 50.0, 200.0]), np.array([t["batch1"], t["batch50"], t["batch200"]])


def batched_eq(sizes: np.ndarray) -> float:
    x, y = eq_curve()
    return float(sum(n * np.interp(np.log(n), np.log(x), y) for n in sizes))


def per_run_extra(tag: str, var: str) -> Dict[int, Tuple[float, int]]:
    """vpos -> (proxy equivalents at the batched cost, number of rank_key calls)."""
    src = SIZE_SOURCE.get(tag, tag)
    out = {}
    for vpos, v in enumerate(B.voxels()):
        f = B.CACHE / "e2e" / src / f"v{v}_{var}.npz"
        if f.exists():
            d = np.load(f)
            out[vpos] = (batched_eq(d["call_sizes"]), len(d["call_sizes"]))
    return out


def summary() -> None:
    eq_ref = json.loads((B.OUT / "timing_batch.json").read_text())["proxy_equivalents"]["ref"]
    ys = eq_curve()[1]
    batched_txt = ", ".join(f"{k}: {v:.2f}" for k, v in zip(("1", "50", "200"), ys))
    L: List[str] = [
        "=== T4: proxy rerank at keep 1/8 and 1/6 (200 voxels x 2 variants, seed 0, paired) ===",
        "proxy cost per candidate (T3, single-worker): per-candidate pass "
        f"{eq_ref:.2f} eq; batched {batched_txt} "
        "(interpolated in log batch size at each rank_key call)",
    ]
    out: Dict[str, Any] = dict(variants={})
    L.append("candidates per rank_key call by level (mean [min-max] over the 400 runs):")
    out["call_sizes"] = {}
    for tag in ("p_k4b", "p_k6", "p_k8"):
        per: Dict[int, List[int]] = {lv: [] for lv in range(4)}
        for f in (B.CACHE / "e2e" / tag).glob("v*_*.npz"):
            if "_s1" in f.name or "_s2" in f.name:
                continue
            d = np.load(f)
            for lv, n in zip(d["call_levels"], d["call_sizes"]):
                per[int(lv)].append(int(n))
        L.append(
            f"  {tag}: "
            + "; ".join(f"L{lv} {np.mean(v):.1f} [{min(v)}-{max(v)}]" for lv, v in per.items() if v)
        )
        out["call_sizes"][tag] = {lv: [float(np.mean(v)), min(v), max(v)] for lv, v in per.items()}
    pooled: Dict[str, Dict[str, Any]] = {}
    for var in C.VARIANTS:
        nm = "realistic" if var == "all" else "clean"
        e0 = S.runs_e0(var, [0])
        rows: Dict[str, Dict[str, Any]] = {"E0 baseline": S.row("E0 baseline", e0, e0)}
        runs: Dict[str, Any] = {}
        for tag, label in ROWS:
            pr = S.runs_files(B.CACHE / "e2e" / tag, var, [0], True)
            if not pr:
                continue
            runs[tag] = pr
            ex = per_run_extra(tag, var)
            r = S.row(label, pr, e0, eq_ref)
            keys = r["_errs"].keys()
            sc_eq = np.array([ex[k[0]][0] for k in keys if k[0] in ex])
            tot = np.array([pr[k]["eg"] + pr[k]["el"] for k in keys if k[0] in ex])
            r["proxy_eq_batched"] = float(sc_eq.mean()) if len(sc_eq) else float("nan")
            r["evals_batched_proxy"] = float((tot + sc_eq).mean()) if len(sc_eq) else float("nan")
            r["n_calls"] = float(np.mean([ex[k[0]][1] for k in keys if k[0] in ex]))
            rows[label] = r
        b = rows["E0 baseline"]
        ref = rows.get("proxy keep 1/4")
        L.append(f"--- {nm} ---")
        for r in rows.values():
            L.append("  " + S.fmt_row(r))
        for lab, r in rows.items():
            if "proxy_eq_batched" in r:
                L.append(
                    f"    {lab}: proxy {r['proxy_eq_batched']:.0f} eq per run at the batched cost "
                    f"({r['n_calls']:.1f} rank_key calls per run) -> total "
                    f"{r['evals_batched_proxy']:.0f} vs baseline {b['evals_total']:.0f} "
                    f"({(r['evals_batched_proxy'] / b['evals_total'] - 1) * 100:+.1f}%)"
                )
        for lab, r in rows.items():
            if lab != "E0 baseline" and ref is not None:
                ok_w, ok_e = (
                    r["wrong"] / r["n"] <= ref["hi"],
                    r["evals_batched_proxy"] < b["evals_total"],
                )
                r["useful"] = bool(ok_w and ok_e)
                L.append(
                    f"    {lab}: useful (this variant) = {r['useful']} (wrong {r['rate']:.3f} <= "
                    f"keep-1/4 upper bound {ref['hi']:.3f}: {ok_w}; total evals below baseline: "
                    f"{ok_e})"
                )
        for lab, r in rows.items():
            if lab == "E0 baseline":
                continue
            for other in (b, ref):
                if other is not None and other is not r:
                    L.append("    " + S.fmt_mc(S.mcnemar(r, other)))
        for lab, r in rows.items():
            if lab == "E0 baseline":
                continue
            pool = pooled.setdefault(lab, dict(wrong=0, n=0, errs={}, base={}, tot=[], base_tot=[]))
            pool["wrong"] += r["wrong"]
            pool["n"] += r["n"]
            pool["errs"].update({(var, k): v for k, v in r["_errs"].items()})
            pool["tot"].append(r["evals_batched_proxy"] * r["n"])
            pool["base_tot"].append(b["evals_total"] * r["n"])
            pool["wr"] = pool.get("wr", {})
            pool["wr"][var] = r
        pooled.setdefault("E0 baseline", dict(errs={}))["errs"].update(
            {(var, k): v for k, v in b["_errs"].items()}
        )
        out["variants"][var] = {
            k: {kk: vv for kk, vv in v.items() if kk != "_errs"} for k, v in rows.items()
        }
    # pooled verdict
    L.append("--- pooled over both variants (400 runs) ---")
    ref = pooled.get("proxy keep 1/4")
    ref_ci = shared_stats.wilson(ref["wrong"], ref["n"]) if ref else None
    verdicts = {}
    for lab, pool in pooled.items():
        if lab == "E0 baseline":
            continue
        w, n = pool["wrong"], pool["n"]
        lo, hi = shared_stats.wilson(w, n)
        bw = sum(r["base_wrong"] for r in pool["wr"].values())

        def mc(o: Dict[Any, float]) -> str:
            ks = sorted(set(pool["errs"]) & set(o))
            wa = np.array([pool["errs"][k] > C.WRONG_DEG for k in ks])
            wb = np.array([o[k] > C.WRONG_DEG for k in ks])
            na, nb = shared_stats.paired_discordant(wa, wb)
            return (
                f"discordant {na}/{nb}, exact McNemar p = {shared_stats.mcnemar_exact(na, nb):.3g}"
            )

        tot, btot = sum(pool["tot"]) / n, sum(pool["base_tot"]) / n
        useful = bool(ref_ci and w / n <= ref_ci[1] and tot < btot)
        verdicts[lab] = dict(wrong=w, n=n, lo=lo, hi=hi, total=tot, base_total=btot, useful=useful)
        L.append(
            f"  {lab}: wrong {w}/{n} = {w / n:.3f} [{lo:.3f},{hi:.3f}] (baseline {bw}); "
            f"total {tot:.0f} vs baseline {btot:.0f} ({(tot / btot - 1) * 100:+.1f}%); "
            f"vs baseline: {mc(pooled['E0 baseline']['errs'])}; "
            + (f"vs keep 1/4: {mc(ref['errs'])}; " if ref and pool is not ref else "")
            + f"USEFUL (wrong <= keep-1/4 upper bound {ref_ci[1]:.3f}, fewer evals than baseline): "
            f"{useful}"
        )
    out["pooled"] = verdicts
    L.append(
        "useful = per-row wrong rate (pooled over the 400 runs) <= the Wilson upper bound of the "
        "keep-1/4 proxy row, and total evaluations (reconstructor + proxy at the batched cost) "
        "below the baseline's. Per-variant numbers above are given for reference."
    )
    (B.OUT / "keep_eighth.txt").write_text("\n".join(L) + "\n")
    C.save_json(B.OUT / "keep_eighth.json", out)
    print("\n".join(L))
    write_table(out)


def write_table(out: Dict[str, Any]) -> None:
    rows = []
    for var, d in out["variants"].items():
        nm = "realistic" if var == "all" else "clean"
        for lab, r in d.items():
            rows.append(
                dict(
                    variant=nm, row=lab, wrong=r["wrong"], rate=r["rate"],
                    ci=f"[{r['lo']:.3f}, {r['hi']:.3f}]", med=r["median_err_right"],
                    glob=r["evals_global"], loc=r["evals_local"], scored=r["scored"],
                    peq=r.get("proxy_eq_batched"),
                    tot=r.get("evals_batched_proxy", r["evals_total"]),
                )
            )  # fmt: skip
    tab = doc_tables.markdown_table(
        rows,
        ["variant", "row", "wrong", "rate", "ci", "med", "glob", "loc", "scored", "peq", "tot"],
        formats=dict(
            rate=".3f", med=".4f", glob=".0f", loc=".0f", scored=".0f", peq=".0f", tot=".0f"
        ),
        labels=["variant", "row", "wrong /200", "rate", "Wilson 95%", "median err right (deg)",
                "global evals", "local evals", "scored", "proxy eq (batched)", "total evals"],
    )  # fmt: skip
    doc_tables.write_tables(B.OUT / "keep_eighth_tables.md", {"t4_keep_eighth": tab})


# -- single-worker timing --------------------------------------------------------------------
def timing(a: argparse.Namespace) -> None:
    import torch
    import endtoend as E
    import models as MD
    from icenine.reconstructor import AdaptiveVoxelReconstructor  # noqa: F401

    sys.path.insert(0, str(B.ICENINE_PY / "scripts" / "common"))
    import preflight

    torch.set_num_threads(1)
    B.OUT.mkdir(parents=True, exist_ok=True)
    info = preflight.require_quiet()
    (B.OUT / "keep_eighth_preflight.json").write_text(json.dumps(info, indent=1, default=str))
    C.init_worker(C.worker_args())
    W = C.get_worker()
    arms: List[Tuple[str, float]] = []
    for s in a.arms:
        name, _, k = s.partition(":")
        arms.append((name, float(k) if k else 0.0))
    kind, q, c8 = MD.SETS["lowq5+c8"]
    items = [(v, vpos, var) for vpos, v in enumerate(B.voxels()[: a.n_vox]) for var in C.VARIANTS]
    res: List[Dict[str, Any]] = []
    for i, (vidx, vpos, var) in enumerate(items):
        order = arms[i % len(arms) :] + arms[: i % len(arms)]  # rotate so no arm is always first
        W_, vctx, fe = E._prep(vidx, var, q)
        fold = int(M.fold_of(np.array([vpos]))[0])
        model = E._model(str(B.CACHE / "models"), "lowq5+c8", "reg", fold)
        vertices, phase = vctx.vertices, vctx.voxel.phase
        for name, keep in order:
            st: Dict[str, Any] = dict(sec=0.0, n=0)

            def rank_key(level: int, cands: List[Any]) -> np.ndarray:
                t0 = time.perf_counter()
                X = E.proxy_features(fe, cands, vertices, phase, c8, True)
                k = -MD.predict_score("reg", model, X)
                st["sec"] += time.perf_counter() - t0
                st["n"] += len(cands)
                return k

            W.rec.rank_key = None if name == "base" else rank_key
            prev = W.rec.keep_fraction
            if name != "base":
                W.rec.keep_fraction = keep
            t0 = time.perf_counter()
            try:
                with C.quiet():
                    r = W.rec.reconstruct_voxel(vertices, phase, rng=C.run_seed(vpos, 0))
            finally:
                W.rec.rank_key, W.rec.keep_fraction = None, prev
            dt = time.perf_counter() - t0
            g, loc, _ = W.rec.last_eval_counts
            res.append(dict(voxel=vidx, variant=var, arm=name, wall=dt, proxy_sec=st["sec"],
                            n_scored=st["n"], evals_global=g, evals_local=loc,
                            R_final=np.asarray(r.orientation, float).tolist()))  # fmt: skip
        print(f"timed {i + 1}/{len(items)}", flush=True)
    info2 = preflight.preflight()
    C.save_json(
        B.OUT / "keep_eighth_timing.json",
        dict(label="single-worker", arms=a.arms, n_vox=a.n_vox, runs=res, preflight_end=info2),
    )


def timing_summary() -> None:
    d = json.loads((B.OUT / "keep_eighth_timing.json").read_text())
    arms = [s.partition(":")[0] for s in d["arms"]]
    keys = sorted({(r["voxel"], r["variant"]) for r in d["runs"]})
    by = {(r["voxel"], r["variant"], r["arm"]): r for r in d["runs"]}
    rows, L = [], [
        f"=== T4 timing, {d['label']}, {len(keys)} cases (voxel x variant), interleaved ==="
    ]
    base_wall = np.array([by[k + ("base",)]["wall"] for k in keys])
    base_ev = np.array(
        [by[k + ("base",)]["evals_global"] + by[k + ("base",)]["evals_local"] for k in keys]
    )
    for arm in arms:
        w = np.array([by[k + (arm,)]["wall"] for k in keys])
        ev = np.array(
            [by[k + (arm,)]["evals_global"] + by[k + (arm,)]["evals_local"] for k in keys]
        )
        px = np.array([by[k + (arm,)]["proxy_sec"] for k in keys])
        same = sum(by[k + (arm,)]["R_final"] == by[k + ("base",)]["R_final"] for k in keys)
        stored = 0
        if arm != "base":
            for k in keys:
                f = B.CACHE / "e2e" / ("p_k4b" if arm == "p_i" else arm) / f"v{k[0]}_{k[1]}.npz"
                stored += int(np.array_equal(np.load(f)["R_final"], by[k + (arm,)]["R_final"]))
        rows.append(
            dict(
                arm=arm, wall=w.mean(), speed=base_wall.sum() / w.sum(), ev=ev.mean(),
                dev=(ev.mean() / base_ev.mean() - 1) * 100, proxy=px.mean(),
                msev=w.sum() / ev.sum() * 1e3, same=same,
            )
        )  # fmt: skip
        L.append(
            f"  {arm}: wall {w.mean():.2f}s (base {base_wall.mean():.2f}s), time ratio to base "
            f"{w.sum() / base_wall.sum():.3f}, evals {ev.mean():.0f} (base {base_ev.mean():.0f}), "
            f"proxy {px.mean():.2f}s, ms per eval {w.sum() / ev.sum() * 1e3:.3f} "
            f"(base {base_wall.sum() / base_ev.sum() * 1e3:.3f}), "
            f"R_final equal to base in {same}/{len(keys)}, "
            f"equal to the stored 10-worker run in {stored}/{len(keys)}"
        )
    tab = doc_tables.markdown_table(
        rows, ["arm", "wall", "speed", "ev", "dev", "proxy", "msev", "same"],
        formats=dict(wall=".2f", speed=".3f", ev=".0f", dev="+.1f", proxy=".2f", msev=".3f"),
        labels=["arm", "wall s/run", "speed vs base", "evals/run", "evals vs base %", "proxy s/run",
                "ms per eval", "R_final same as base"],
    )  # fmt: skip
    doc_tables.write_tables(B.OUT / "keep_eighth_timing_tables.md", {"t4_keep_eighth_timing": tab})
    (B.OUT / "keep_eighth_timing.txt").write_text("\n".join(L) + "\n")
    print("\n".join(L))


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    ap.add_argument("cmd", choices=["summary", "timing", "timing-summary"])
    ap.add_argument("--arms", nargs="+", default=["base", "p_i:0.25", "p_k8:0.125"])
    ap.add_argument("--n-vox", type=int, default=20)
    a = ap.parse_args()
    {"summary": summary, "timing": lambda: timing(a), "timing-summary": timing_summary}[a.cmd]()


if __name__ == "__main__":
    main()
