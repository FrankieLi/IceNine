#!/usr/bin/env python3
"""E2 end-to-end summary: classifier as pruning reranker (rerank), classifier as final chooser
among the F1 candidates (final), each against its non-learned counterpart; writes
benchmarks/findoptimal_robustness/e2_endtoend.{txt,json}."""

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import common as C  # noqa: E402
import summarize_fixes as SF  # noqa: E402

FEATURE_EVALS = 3  # one feature pass = local cost + pixel-radius-3 cost + one geometry pass


def main(model="GBT"):
    info = dict(np.load(C.OUT_DIR / "voxels.npz"))
    E0 = SF.load_e0(info)
    vox = [int(v) for v in info["voxel_indices"]][: C.N_VOXELS]
    out, lines = {}, []
    for var in C.VARIANTS:
        rows = []
        k0 = [(i, var, 0) for i in range(len(vox)) if (i, var, 0) in E0]
        err0 = np.array([C.err_deg(E0[k]["R"], E0[k]["R_true"]) for k in k0])
        ev0 = np.array([E0[k]["evals"] for k in k0], float)
        rt0 = np.array([E0[k]["rt"] for k in k0])
        rows.append(SF.row("E0 baseline, seed 0", err0, ev0, rt0, ev0.mean(), rt0.mean()))
        # (a) rerank
        errs, evs, rts, ids = [], [], [], []
        for vpos, v in enumerate(vox):
            f = C.CACHE_DIR / f"e2_rerank_{model}" / f"v{v}_{var}.npz"
            if f.exists() and (vpos, var, 0) in E0:
                d = np.load(f)
                errs.append(float(C.err_deg(d["R_final"], d["R_true"])))
                evs.append(
                    int(d["evals_global"] + d["evals_local"]) + FEATURE_EVALS * int(d["n_scored"])
                )
                rts.append(float(d["runtime"]))
                ids.append(vpos)
        if errs:
            pb = np.array([C.err_deg(E0[(i, var, 0)]["R"], E0[(i, var, 0)]["R_true"]) for i in ids])
            errs = np.array(errs)
            rows.append(
                SF.row(
                    f"C-a classifier rerank ({model})",
                    errs,
                    np.array(evs),
                    np.array(rts),
                    ev0.mean(),
                    rt0.mean(),
                    dict(
                        fixed=int(((errs <= 1) & (pb > 1)).sum()),
                        broken=int(((errs > 1) & (pb <= 1)).sum()),
                        baseline_wrong_same_cases=int((pb > 1).sum()),
                        n_paired=len(pb),
                    ),
                )
            )
            # rerank + F1
            e2, v2 = [], []
            for vpos in ids:
                v = vox[vpos]
                fp = C.CACHE_DIR / f"f1_e2_rerank_{model}" / f"v{v}_{var}.npz"
                if fp.exists():
                    d = np.load(C.CACHE_DIR / f"e2_rerank_{model}" / f"v{v}_{var}.npz")
                    base = dict(R=d["R_final"], cost=float(d["cost_final"]), R_true=d["R_true"])
                    R, c, extra = SF.f1_answer(np.load(fp), 0, 29, base, d["R_true"])
                    e2.append(float(C.err_deg(R, d["R_true"])))
                    v2.append(
                        int(d["evals_global"] + d["evals_local"])
                        + FEATURE_EVALS * int(d["n_scored"])
                        + extra
                    )
            if e2:
                e2 = np.array(e2)
                rows.append(
                    SF.row(
                        "C-a + F1 (Sigma<=29)",
                        e2,
                        np.array(v2),
                        np.array(v2) * rt0.mean() / ev0.mean(),
                        ev0.mean(),
                        rt0.mean(),
                    )
                )
        # (b) final choice by classifier among F1 candidates (all seeds) vs by cost
        ec, evc, e_cost, ev_cost = [], [], [], []
        for vpos, v in enumerate(vox):
            f = C.CACHE_DIR / f"e2_final_{model}" / f"v{v}_{var}.npz"
            fp = C.CACHE_DIR / "f1" / f"v{v}_{var}.npz"
            if not f.exists():
                continue
            d, f1 = np.load(f), np.load(fp)
            for s in range(C.N_SEEDS):
                R = d[f"s{s}_R"]
                sc = d[f"s{s}_score"]
                ec.append(float(C.err_deg(R[int(np.argmax(sc))], d["R_true"])))
                e_cost.append(float(C.err_deg(R[int(np.argmin(d[f"s{s}_cost"]))], d["R_true"])))
                k = (vpos, var, s)
                _, _, extra = SF.f1_answer(f1, s, 29, E0[k], d["R_true"])
                evc.append(E0[k]["evals"] + extra + FEATURE_EVALS * len(R))
                ev_cost.append(E0[k]["evals"] + extra)
        if ec:
            ec, e_cost = np.array(ec), np.array(e_cost)
            rows.append(
                SF.row(
                    "F1 final = lowest cost among {answer, refined relatives} (same candidates)",
                    e_cost,
                    np.array(ev_cost),
                    np.array(ev_cost) * rt0.mean() / ev0.mean(),
                    ev0.mean(),
                    rt0.mean(),
                )
            )
            rows.append(
                SF.row(
                    f"C-b final = classifier ({model}) among the same candidates",
                    ec,
                    np.array(evc),
                    np.array(evc) * rt0.mean() / ev0.mean(),
                    ev0.mean(),
                    rt0.mean(),
                    dict(
                        fixed=int(((ec <= 1) & (e_cost > 1)).sum()),
                        broken=int(((ec > 1) & (e_cost <= 1)).sum()),
                    ),
                )
            )
        out[var] = rows
        lines.append(f"=== {var} ===")
        for r in rows:
            note = (
                f"  vs its counterpart: fixed {r['fixed']}, broken {r['broken']}"
                if "fixed" in r
                else ""
            )
            lines.append(
                f"{r['name']:80s} n={r['n']:4d} wrong {r['wrong']:4d} = {r['rate']:.3f} "
                f"[{r['lo']:.3f},{r['hi']:.3f}]  evals {r['evals_mean']:.0f} (extra "
                f"{r['extra_evals']:.0f})  time {r['runtime_mean']:.1f}s{note}"
            )
    (C.OUT_DIR / "e2_endtoend.txt").write_text("\n".join(lines) + "\n")
    C.save_json(C.OUT_DIR / "e2_endtoend.json", out)
    print("\n".join(lines))


if __name__ == "__main__":
    main()
