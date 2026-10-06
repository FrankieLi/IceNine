#!/usr/bin/env python3
"""E1 summary: wrong rate (with Wilson 95% CI), median error of the right answers, extra cost
evaluations and runtime for the baseline and every fix, per variant.

Reads cache/e0 (E0), cache/f1 (F1, post hoc), cache/fix/<F>/ (F1b, F2, F3 re-runs) and, if present,
benchmarks/findoptimal_robustness/e2_endtoend.json (classifier rows, written by e2_models.py);
writes fixes_summary.{txt,json}.
"""

import json
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import common as C  # noqa: E402

SIGS = (11, 29)


def row(
    name: str,
    err: np.ndarray,
    evals: np.ndarray,
    rt: np.ndarray,
    base_evals: float,
    base_rt: float,
    extra: Optional[Dict[str, Any]] = None,
):
    n = len(err)
    k = int((err > C.WRONG_DEG).sum())
    p, lo, hi = C.wilson(k, n)
    right = err[err <= C.WRONG_DEG]
    r = dict(
        name=name, n=n, wrong=k, rate=p, lo=lo, hi=hi,
        median_err_right=float(np.median(right)) if len(right) else float("nan"),
        evals_mean=float(np.mean(evals)), extra_evals=float(np.mean(evals) - base_evals),
        runtime_mean=float(np.mean(rt)), extra_runtime=float(np.mean(rt) - base_rt),
    )  # fmt: skip
    if extra:
        r.update(extra)
    return r


def load_e0(info):
    vox = [int(v) for v in info["voxel_indices"]][: C.N_VOXELS]
    D: Dict[Any, Dict[str, Any]] = {}
    for vpos, v in enumerate(vox):
        for var in C.VARIANTS:
            f = C.CACHE_DIR / "e0" / f"v{v}_{var}.npz"
            if not f.exists():
                continue
            d = np.load(f)
            if "unbuildable" in d.files:
                continue
            for s in range(C.N_SEEDS):
                D[(vpos, var, s)] = dict(
                    v=v, R=d[f"s{s}_R_final"], R_true=d["R_true"], cost=float(d[f"s{s}_cost_final"]),
                    evals=int(d[f"s{s}_evals_global"] + d[f"s{s}_evals_local"]),
                    rt=float(d[f"s{s}_runtime"]),
                )  # fmt: skip
    return D


def f1_answer(f1, s, S, base, R_true):
    """F1 answer for run s: the lowest final local cost among the original answer and the refined
    relatives with Sigma <= S in the top-3 (by post-quick-MC cost) of that Sigma subset. Returns
    (R, cost, extra local evals)."""
    sig = f1[f"s{s}_rel_sigma"]
    cost = f1[f"s{s}_rel_cost_post"]
    ev = f1[f"s{s}_rel_evals"]
    idx = np.nonzero(sig <= S)[0]
    top = set(idx[np.argsort(cost[idx])[:3]].tolist())
    ref_idx = f1[f"s{s}_ref_idx"].tolist()
    best_R, best_c = base["R"], base["cost"]
    extra = int(ev[idx].sum())
    for j, i in enumerate(ref_idx):
        if i in top:
            extra += int(f1[f"s{s}_ref_evals"][j])
            if f1[f"s{s}_ref_cost"][j] < best_c:
                best_R, best_c = f1[f"s{s}_ref_R"][j], float(f1[f"s{s}_ref_cost"][j])
    return best_R, best_c, extra


def main():
    info = dict(np.load(C.OUT_DIR / "voxels.npz"))
    E0 = load_e0(info)
    vox = [int(v) for v in info["voxel_indices"]][: C.N_VOXELS]
    out: Dict[str, Any] = {}
    lines: List[str] = []
    for var in C.VARIANTS:
        rows: List[Dict[str, Any]] = []
        keys0 = sorted(k for k in E0 if k[1] == var)
        if not keys0:
            continue
        err0 = np.array([C.err_deg(E0[k]["R"], E0[k]["R_true"]) for k in keys0])
        ev0 = np.array([E0[k]["evals"] for k in keys0], float)
        rt0 = np.array([E0[k]["rt"] for k in keys0])
        rows.append(row("E0 baseline (all seeds)", err0, ev0, rt0, ev0.mean(), rt0.mean()))
        # seed-0 subset of the baseline (the fix re-runs are seed 0)
        k_s0 = [i for i, k in enumerate(keys0) if k[2] == 0]
        base_s0 = row(
            "E0 baseline, seed 0",
            err0[k_s0],
            ev0[k_s0],
            rt0[k_s0],
            ev0[k_s0].mean(),
            rt0[k_s0].mean(),
        )
        rows.append(base_s0)
        # F1 (all seeds)
        for S in SIGS:
            errs, evs, rts, nfix, nbreak = [], [], [], 0, 0
            for k in keys0:
                vpos, _, s = k
                fp = C.CACHE_DIR / "f1" / f"v{E0[k]['v']}_{var}.npz"
                if not fp.exists():
                    continue
                f1 = np.load(fp)
                R, c, extra = f1_answer(f1, s, S, E0[k], E0[k]["R_true"])
                e = float(C.err_deg(R, E0[k]["R_true"]))
                e_base = float(C.err_deg(E0[k]["R"], E0[k]["R_true"]))
                errs.append(e)
                evs.append(E0[k]["evals"] + extra)
                rts.append(E0[k]["rt"] + extra * E0[k]["rt"] / max(E0[k]["evals"], 1))  # approx.
                nfix += (e <= 1.0) and (e_base > 1.0)
                nbreak += (e > 1.0) and (e_base <= 1.0)
            if errs:
                rows.append(
                    row(
                        f"F1 CSL relatives Sigma<={S}",
                        np.array(errs),
                        np.array(evs),
                        np.array(rts),
                        ev0.mean(),
                        rt0.mean(),
                        dict(fixed=int(nfix), broken=int(nbreak)),
                    )
                )
        # F4 best of N seeds by final local cost (per voxel-variant)
        for N in (2, 3):
            errs, evs, rts = [], [], []
            for vpos in range(len(vox)):
                runs = [E0.get((vpos, var, s)) for s in range(N)]
                if any(r is None for r in runs):
                    continue
                b = min(runs, key=lambda r: r["cost"])
                errs.append(float(C.err_deg(b["R"], b["R_true"])))
                evs.append(sum(r["evals"] for r in runs))
                rts.append(sum(r["rt"] for r in runs))
            if errs:
                rows.append(
                    row(
                        f"F4 best of {N} seeds (by local cost)",
                        np.array(errs),
                        np.array(evs),
                        np.array(rts),
                        ev0.mean(),
                        rt0.mean(),
                    )
                )
        # F4 + F1: best of N F1 answers
        for N in (2, 3):
            errs, evs = [], []
            for vpos in range(len(vox)):
                cand = []
                ok = True
                for s in range(N):
                    k = (vpos, var, s)
                    if k not in E0:
                        ok = False
                        break
                    fp = C.CACHE_DIR / "f1" / f"v{E0[k]['v']}_{var}.npz"
                    if not fp.exists():
                        ok = False
                        break
                    R, c, extra = f1_answer(np.load(fp), s, 29, E0[k], E0[k]["R_true"])
                    cand.append((c, R, E0[k]["evals"] + extra, E0[k]["R_true"]))
                if ok and cand:
                    b = min(cand, key=lambda t: t[0])
                    errs.append(float(C.err_deg(b[1], b[3])))
                    evs.append(sum(t[2] for t in cand))
            if errs:
                rows.append(
                    row(
                        f"F4+F1 best of {N} F1 answers",
                        np.array(errs),
                        np.array(evs),
                        np.array(evs) * rt0.mean() / ev0.mean(),
                        ev0.mean(),
                        rt0.mean(),
                    )
                )
        # re-run fixes (seed 0)
        for fix in ("F1b", "F2a", "F2b", "F3a", "F3b", "F3c") + tuple(
            sorted(
                p.name
                for p in (C.CACHE_DIR / "fix").glob("*")
                if p.name not in ("F1b", "F2a", "F2b", "F3a", "F3b", "F3c")
            )
        ):
            d_fix = C.CACHE_DIR / "fix" / fix
            if not d_fix.exists():
                continue
            errs, evs, rts, pair_base = [], [], [], []
            for vpos, v in enumerate(vox):
                fs_ = sorted(d_fix.glob(f"v{v}_{var}_s*.npz"))
                if not fs_ or (vpos, var, 0) not in E0:
                    continue
                d = np.load(fs_[0])
                if "s0_R_final" not in d.files:
                    continue
                k = (vpos, var, 0)
                errs.append(float(C.err_deg(d["s0_R_final"], E0[k]["R_true"])))
                evs.append(int(d["s0_evals_global"] + d["s0_evals_local"]))
                rts.append(float(d["s0_runtime"]))
                pair_base.append(float(C.err_deg(E0[k]["R"], E0[k]["R_true"])))
            if errs:
                errs, pb = np.array(errs), np.array(pair_base)
                b0 = row("base", pb, np.array([0.0]), np.array([0.0]), 0.0, 0.0)
                nfix = int(((errs <= 1) & (pb > 1)).sum())
                nbr = int(((errs > 1) & (pb <= 1)).sum())
                rows.append(
                    row(
                        f"{fix} (seed 0)",
                        errs,
                        np.array(evs),
                        np.array(rts),
                        np.mean(
                            [
                                E0[(vpos, var, 0)]["evals"]
                                for vpos in range(len(vox))
                                if (vpos, var, 0) in E0
                            ]
                        ),
                        base_s0["runtime_mean"],
                        dict(
                            fixed=nfix,
                            broken=nbr,
                            baseline_wrong_same_cases=int((pb > 1).sum()),
                            n_paired=len(pb),
                        ),
                    )
                )
        out[var] = rows
        lines.append(f"=== {var} ===")
        lines.append(
            f"{'method':42s} {'n':>4s} {'wrong':>5s} {'rate [95% CI]':>22s} {'med err right':>13s} {'evals/run':>10s} {'extra':>9s} {'time/run':>9s} {'extra t':>8s}  note"
        )
        for r in rows:
            note = f"fixed {r['fixed']}, broken {r['broken']}" if "fixed" in r else ""
            if "baseline_wrong_same_cases" in r:
                note += f" (baseline wrong on these cases: {r['baseline_wrong_same_cases']}/{r['n_paired']})"
            lines.append(
                f"{r['name']:42s} {r['n']:4d} {r['wrong']:5d} {r['rate']:7.3f} [{r['lo']:.3f},{r['hi']:.3f}] "
                f"{r['median_err_right']:13.4f} {r['evals_mean']:10.0f} {r['extra_evals']:9.0f} {r['runtime_mean']:9.1f} {r['extra_runtime']:8.1f}  {note}"
            )
    txt = "\n".join(lines)
    (C.OUT_DIR / "fixes_summary.txt").write_text(txt + "\n")
    C.save_json(C.OUT_DIR / "fixes_summary.json", out)
    print(txt)


if __name__ == "__main__":
    main()
