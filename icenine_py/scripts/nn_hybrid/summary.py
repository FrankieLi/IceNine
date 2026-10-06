#!/usr/bin/env python3
"""
Metrics, success criteria and decision points of the hybrid network -> FindOptimal study.

Reads benchmarks/nn_hybrid/<PIPE>_<model>_raw.npz (scripts/nn_hybrid/run.py) and the sweep raws
of benchmarks/toy_orientation_sweep/, writes <out-dir>/summary.txt and summary.json.

All arrays are re-indexed to the sweep's voxel order (findoptimal_b_raw and the hybrid raws are
sorted by voxel index; the perturbation / optimizer raws are in sweep order).

Errors: cubic-symmetry-reduced misorientation to the truth (FindOptimal can return an equivalent
orientation); "unreduced" is the plain angle. "Wrong" = reduced error > 1 deg. Metrics are over
the sweep's pass-1 success mask (fail_pass1 == 0); cases whose pass 1 failed (the net has no
estimate) are fall-backs (FindOptimal from the perturbed nominal) and are counted separately.
Paired win = the pipeline's error is smaller than the comparator's by more than 0.002 deg; a tie
is |difference| < 0.002 deg (ties count as half a win in "win" and are also reported).

Usage (from icenine_py/):  uv run python scripts/nn_hybrid/summary.py [--out-dir DIR] [--model M]
"""

import argparse
import json
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import numpy as np

HERE = Path(__file__).resolve().parent
ICENINE_PY = HERE.parents[1]
sys.path.insert(0, str(HERE.parent))
sys.path.insert(0, str(ICENINE_PY / "benchmarks"))

SWEEP_DIR = ICENINE_PY / "benchmarks" / "toy_orientation_sweep"
OUT_DIR = ICENINE_PY / "benchmarks" / "nn_hybrid"
VARIANT_LABEL = {0: "clean", 1: "realistic"}
TIE_DEG = 0.002
PIPES = ["H0", "N1", "N3", "H1", "H3", "H3c", "H3m", "HG"]
TAU_GRID = [0.1, 0.2, 0.3, 0.5, 1.0]


def wilson(k: int, n: int, z: float = 1.96) -> Tuple[float, float]:
    if n == 0:
        return float("nan"), float("nan")
    p = k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return float(c - h), float(c + h)


def reorder(arr: np.ndarray, src: np.ndarray, dst: np.ndarray) -> np.ndarray:
    """arr's first axis is in voxel order `src`; return it in order `dst` (voxels absent from
    `src`, as in a pilot on a subset, are NaN / zero)."""
    pos = {int(v): i for i, v in enumerate(src)}
    out = np.full((len(dst),) + arr.shape[1:], np.nan if arr.dtype.kind == "f" else 0, dtype=arr.dtype)
    for k, v in enumerate(dst):
        if int(v) in pos:
            out[k] = arr[pos[int(v)]]
    return out


class Study:
    """The sweep reference data and the pipelines' raws, in sweep voxel order."""

    def __init__(self, out_dir: Path, model: str = "realistic_s0", n_voxels: int = 0) -> None:
        import findoptimal_sweep as fs

        self.fs = fs
        self.out_dir, self.model = out_dir, model
        sw = np.load(SWEEP_DIR / "perturbation_sweep_raw.npz")
        nv = n_voxels or len(sw["voxel_indices"])  # n_voxels > 0: first n sweep voxels (pilot)
        self.vox = sw["voxel_indices"][:nv]
        self.radii = sw["radii"]
        self.mi = [str(m) for m in sw["models"]].index(model)
        self.mask = (sw["fail_pass1"] == 0)[:nv]  # (V, R, D, 2) pass-1 successes
        self.sweep = {k: sw[k][:nv] if k in ("err_angle", "n_roi", "fail_pass1") else sw[k] for k in sw.files}
        a = np.load(SWEEP_DIR / "findoptimal_a_raw.npz")
        self.R_true = np.stack(
            [a["R_true"][np.nonzero((a["voxel_indices"] == v) & (a["variant_index"] == 0))[0][0]]
             for v in self.vox]
        )  # fmt: skip
        b = np.load(SWEEP_DIR / "findoptimal_b_raw.npz")
        self.b = {k: reorder(b[k], b["voxel_indices"], self.vox)
                  for k in ("R_final", "cost_final", "cost_true", "runtime", "evals", "ran")}  # fmt: skip
        self.opt = np.load(SWEEP_DIR / "optimizer_sweep_raw.npz")
        self.opt = {k: self.opt[k][:nv] for k in ("err_angle", "voxel_indices")}
        assert (self.opt["voxel_indices"] == self.vox).all()
        self.raws: Dict[str, Dict[str, np.ndarray]] = {}

    def load(self, pipe: str, model: Optional[str] = None, variant: Optional[str] = None) -> bool:
        model = model or self.model
        name = f"{pipe}_{model}" + (f"_{variant}" if variant else "") + "_raw.npz"
        p = self.out_dir / name
        if not p.exists():
            return False
        r = np.load(p)
        ris = [int(i) for i in r["radii_idx"]]
        d = {}
        for k in r.files:
            if k in ("voxel_indices", "radii_idx", "radii", "variants"):
                continue
            d[k] = reorder(r[k], r["voxel_indices"], self.vox)
        d["radii_idx"] = np.array(ris)
        key = pipe
        self.raws[key] = d
        return True

    # -- errors (voxel, radius, direction, variant) ----------------------------------------------

    def err(self, R: np.ndarray, reduce: bool = True) -> np.ndarray:
        """Angle (deg) of R (V, r, D, 2, 3, 3) to the truth, NaN where R is NaN."""
        Rt = self.R_true[:, None, None, None]
        out = np.full(R.shape[:4], np.nan)
        ok = np.isfinite(R).all(axis=(-1, -2))
        Rb = np.broadcast_to(Rt, R.shape)
        out[ok] = self.fs.misorientation_deg(R[ok], Rb[ok], reduce=reduce)
        return out

    def pipe_err(self, pipe: str, reduce: bool = True) -> Optional[np.ndarray]:
        """Error array (voxel, radius, direction, variant) over all 10 radii, NaN where the
        pipeline was not run. H0 is the live re-run when present, else the stored experiment B."""
        if pipe in ("N1", "N3"):
            src = next((self.raws[k] for k in ("H3", "H1", "H0", "H3c") if k in self.raws), None)
            if src is None:
                return None
            R, ris = src["R_x1" if pipe == "N1" else "R_x3"], src["radii_idx"]
        elif pipe in self.raws:
            R, ris = self.raws[pipe]["R_final"], self.raws[pipe]["radii_idx"]
        elif pipe == "H0":
            return self.err(self.b["R_final"], reduce)
        else:
            return None
        full = np.full((len(self.vox), 10) + R.shape[2:], np.nan)
        full[:, ris] = R
        return self.err(full, reduce)

    def baseline_err(self, name: str) -> np.ndarray:
        """Unreduced angle of the optimizer sweep's `mc` (one-shot) or `huber` (x3)."""
        mi, p = {"MC": (0, 0), "Huber3": (3, 2)}[name]
        return self.opt["err_angle"][:, :, :, :, mi, p].astype(np.float64)

    def net_err_sweep(self, p: int) -> np.ndarray:
        return self.sweep["err_angle"][:, :, :, :, self.mi, p].astype(np.float64)


# ---------------------------------------------------------------------------
# Metrics
# ---------------------------------------------------------------------------


def stats(e: np.ndarray, e_unred: Optional[np.ndarray]) -> Dict[str, float]:
    e = e[np.isfinite(e)]
    n = len(e)
    if n == 0:
        return dict(n=0, median=np.nan, rms=np.nan, rms_unred=np.nan, f01=np.nan, wrong=np.nan,
                    wrong_lo=np.nan, wrong_hi=np.nan)  # fmt: skip
    k = int((e > 1.0).sum())
    lo, hi = wilson(k, n)
    ru = np.nan
    if e_unred is not None:
        eu = e_unred[np.isfinite(e_unred)]
        ru = float(np.sqrt(np.mean(eu**2))) if len(eu) else np.nan
    return dict(n=n, median=float(np.median(e)), rms=float(np.sqrt(np.mean(e**2))), rms_unred=ru,
                f01=float(np.mean(e < 0.1)), wrong=k / n, wrong_lo=lo, wrong_hi=hi)  # fmt: skip


def win_rate(a: np.ndarray, b: np.ndarray) -> Tuple[float, float, int]:
    """(wins of a over b incl. half ties, tie fraction, n) over cases where both are finite."""
    ok = np.isfinite(a) & np.isfinite(b)
    if ok.sum() == 0:
        return float("nan"), float("nan"), 0
    d = a[ok] - b[ok]
    tie = np.abs(d) < TIE_DEG
    win = (d < -TIE_DEG).sum() + 0.5 * tie.sum()
    return float(win / ok.sum()), float(tie.mean()), int(ok.sum())


def pick(arr: np.ndarray, ri: int, vi: int, mask: np.ndarray) -> np.ndarray:
    """Cases of radius ri, variant vi on the pass-1 success mask, as a 1-D array."""
    return arr[:, ri, :, vi][mask[:, ri, :, vi]]


def build(st: Study) -> Dict[str, Any]:
    """All metrics: res[variant][ri][pipe] -> dict."""
    errs = {p: st.pipe_err(p) for p in PIPES}
    errs_u = {p: st.pipe_err(p, reduce=False) for p in PIPES}
    base = {"MC": st.baseline_err("MC"), "Huber3": st.baseline_err("Huber3")}
    res: Dict[str, Any] = {}
    for vi, vname in VARIANT_LABEL.items():
        res[vname] = {}
        for ri in range(10):
            row: Dict[str, Any] = {}
            m = st.mask
            for p in PIPES:
                if errs[p] is None:
                    continue
                e = pick(errs[p], ri, vi, m)
                if np.isfinite(e).sum() == 0:
                    continue
                d = stats(e, pick(errs_u[p], ri, vi, m))
                for q in ("N3", "H0", "H1"):
                    if q != p and errs.get(q) is not None:
                        w, t, n = win_rate(e, pick(errs[q], ri, vi, m))
                        d[f"win_vs_{q}"], d[f"tie_vs_{q}"] = w, t
                for q in ("MC", "Huber3"):
                    w, t, n = win_rate(e, pick(base[q], ri, vi, m))
                    d[f"win_vs_{q}"], d[f"tie_vs_{q}"] = w, t
                ev = None
                if p in st.raws and ri in st.raws[p]["radii_idx"]:
                    k = list(st.raws[p]["radii_idx"]).index(ri)
                    ev = st.raws[p]["evals"][:, k, :, vi][m[:, ri, :, vi]]
                elif p == "H0":
                    ev = st.b["evals"][:, ri, :, vi][m[:, ri, :, vi]]
                d["evals_mean"] = float(np.mean(ev)) if ev is not None and len(ev) else np.nan
                row[p] = d
            # fall-back counts
            nfb = int((~m[:, ri, :, vi] & (st.sweep["n_roi"][:, ri] > 0)).sum())
            row["_n_fallback"] = nfb
            row["_n_success"] = int(m[:, ri, :, vi].sum())
            res[vname][ri] = row
    return res


# ---------------------------------------------------------------------------
# Tables
# ---------------------------------------------------------------------------


def fmt(x: float, p: int = 3) -> str:
    return "  nan" if x is None or not np.isfinite(x) else f"{x:.{p}g}"


def table(res: Dict[str, Any], radii: np.ndarray, vname: str, key: str, label: str, pct: bool = False) -> str:
    pipes = [p for p in PIPES if any(p in res[vname][ri] for ri in range(10))]
    lines = [f"{label} ({vname})", "  r(deg)  " + "".join(f"{p:>10s}" for p in pipes)]
    for ri in range(10):
        cells = []
        for p in pipes:
            d = res[vname][ri].get(p)
            v = d[key] * (100 if pct else 1) if d else np.nan
            cells.append(f"{fmt(v):>10s}")
        lines.append(f"  {radii[ri]:<7g} " + "".join(cells))
    return "\n".join(lines)


def wrong_table(res: Dict[str, Any], radii: np.ndarray, vname: str) -> str:
    pipes = [p for p in PIPES if any(p in res[vname][ri] for ri in range(10))]
    lines = [f"wrong rate % [95% Wilson] ({vname})", "  r(deg)  " + "".join(f"{p:>20s}" for p in pipes)]
    for ri in range(10):
        cells = []
        for p in pipes:
            d = res[vname][ri].get(p)
            c = (f"{100*d['wrong']:.2f} [{100*d['wrong_lo']:.1f},{100*d['wrong_hi']:.1f}]"
                 if d else "")  # fmt: skip
            cells.append(f"{c:>20s}")
        lines.append(f"  {radii[ri]:<7g} " + "".join(cells))
    return "\n".join(lines)


def win_table(res: Dict[str, Any], radii: np.ndarray, vname: str, pipe: str) -> str:
    cols = [c for c in ("N3", "H0", "H1", "MC", "Huber3") if c != pipe]
    lines = [f"paired win rate of {pipe} vs ... (ties = half; tie fraction in brackets) ({vname})",
             "  r(deg)  " + "".join(f"{c:>16s}" for c in cols)]  # fmt: skip
    for ri in range(10):
        d = res[vname][ri].get(pipe)
        cells = []
        for c in cols:
            if d and f"win_vs_{c}" in d and np.isfinite(d[f"win_vs_{c}"]):
                cells.append(f"{d[f'win_vs_{c}']:.2f} [{d[f'tie_vs_{c}']:.2f}]".rjust(16))
            else:
                cells.append(" " * 16)
        lines.append(f"  {radii[ri]:<7g} " + "".join(cells))
    return "\n".join(lines)


def time_table(st: Study, vi: int) -> str:
    """Mean seconds per case (pass-1 successes) by stage, per pipeline (all radii run)."""
    keys = [("t_prep", "prepare"), ("t_render", "render"), ("t_forward", "forward"), ("t_gn", "GN"),
            ("t_fo", "FindOpt"), ("t_vm", "VarMin"), ("t_eval", "eval(all)"), ("t_finish", "finisher")]  # fmt: skip
    lines = [f"mean seconds per case by stage, 10 workers (contended) ({VARIANT_LABEL[vi]})",
             "  pipe  " + "".join(f"{k[1]:>10s}" for k in keys) + "   evals"]  # fmt: skip
    for p, d in st.raws.items():
        if "t_finish" not in d:
            continue
        m = np.isfinite(d["t_finish"][..., vi]) & (d["fail_pass1"][..., vi] == 0)
        cells = [f"{np.mean(d[k][..., vi][m]) if m.any() else np.nan:>10.3f}" for k, _ in keys]
        ev = np.mean(d["evals"][..., vi][m]) if m.any() else np.nan
        lines.append(f"  {p:<5s} " + "".join(cells) + f"{ev:>9.0f}")
    lines.append("  H0 (stored findoptimal_b_raw): mean finisher seconds "
                 f"{np.nanmean(st.b['runtime'][..., vi][st.mask[..., vi]]):.3f}, "
                 f"evals {np.mean(st.b['evals'][..., vi][st.mask[..., vi]]):.0f}")  # fmt: skip
    return "\n".join(lines)


# ---------------------------------------------------------------------------
# Criteria and decisions
# ---------------------------------------------------------------------------


def med(res: Dict[str, Any], v: str, ri: int, p: str, k: str = "median") -> float:
    d = res[v][ri].get(p)
    return float(d[k]) if d else float("nan")


def criteria(st: Study, res: Dict[str, Any], vname: str) -> List[str]:
    """The plan's success criteria for variant vname (realistic is the primary)."""
    radii = st.radii
    le3 = [ri for ri in range(10) if radii[ri] <= 3.0]
    ge075 = [ri for ri in range(10) if radii[ri] >= 0.75]
    le05 = [ri for ri in range(10) if radii[ri] <= 0.5]
    out = [f"Success criteria, {vname} data (H3 = net x3 -> FindOptimal):"]

    def line(name: str, ok: bool, detail: str) -> None:
        out.append(f"  [{'MET    ' if ok else 'NOT MET'}] {name}: {detail}")

    ms = [med(res, vname, ri, "H3") for ri in le3]
    fs_ = [med(res, vname, ri, "H3", "f01") for ri in le3]
    ok = all(m <= 0.035 for m in ms) and all(f >= 0.9 for f in fs_)
    line("C1 H3 median <= 0.035 deg and fraction < 0.1 deg >= 0.9 for r <= 3", ok,
         "medians " + ", ".join(f"{m:.3g}" for m in ms) + "; frac<0.1 " + ", ".join(f"{f:.2f}" for f in fs_))  # fmt: skip
    w = [med(res, vname, ri, "H3", "wrong") for ri in le3]
    line("C2 H3 wrong <= 1% for r <= 3", all(x <= 0.01 for x in w),
         ", ".join(f"{100*x:.2f}%" for x in w))  # fmt: skip
    w = [med(res, vname, ri, "H3", "win_vs_N3") for ri in le3]
    line("C3 H3 beats N3 in >= 70% paired for r <= 3", all(x >= 0.7 for x in w),
         ", ".join(f"{x:.2f}" for x in w))  # fmt: skip
    w = [med(res, vname, ri, "H3", "win_vs_H0") for ri in ge075]
    line("C4 H3 beats H0 in >= 70% for r >= 0.75", all(x >= 0.7 for x in w),
         "r=" + ",".join(f"{radii[ri]:g}" for ri in ge075) + ": " + ", ".join(f"{x:.2f}" for x in w))  # fmt: skip
    dm = [med(res, vname, ri, "H3") - med(res, vname, ri, "H0") for ri in le05]
    w = [med(res, vname, ri, "H3", "win_vs_H0") for ri in le05]
    ok = all(d <= 0.005 for d in dm) and all(x >= 0.45 for x in w)
    line("C5 r <= 0.5: H3 median within 0.005 deg of H0 and win >= 45%", ok,
         "median diff " + ", ".join(f"{d:+.3f}" for d in dm) + "; win " + ", ".join(f"{x:.2f}" for x in w))  # fmt: skip
    ri5 = 9
    mg = med(res, vname, ri5, "HG")
    line("C6 HG median <= 0.1 deg at r = 5", np.isfinite(mg) and mg <= 0.1, f"HG median {mg:.3g}")
    return out


def decisions(st: Study, res: Dict[str, Any], vname: str, vi: int) -> List[str]:
    radii = st.radii
    out = [f"Decision points, {vname} data:"]
    le2 = [ri for ri in range(10) if radii[ri] <= 2.0]
    # D1
    if all(p in res[vname][ri] for ri in le2 for p in ("H1", "H3")):
        dd = [med(res, vname, ri, "H1") - med(res, vname, ri, "H3") for ri in le2]
        ok = all(d <= 0.005 for d in dd)
        out.append(f"  D1 H1 within 0.005 deg of H3 at r <= 2 (median H1 - H3): "
                   + ", ".join(f"{d:+.4f}" for d in dd) + f" -> {'recommend x1' if ok else 'keep x3'}")  # fmt: skip
    else:
        out.append("  D1 not evaluated (H1 or H3 missing)")
    # D2: floor
    if "H3" in st.raws:
        d = st.raws["H3"]
        ri_ = list(d["radii_idx"])
        rows = []
        for ri in range(10):
            if ri not in ri_ or radii[ri] > 3:
                continue
            k = ri_.index(ri)
            m = st.mask[:, ri, :, vi]
            e3 = st.pipe_err("H3")[:, ri, :, vi][m]
            cf = d["cost_final"][:, k, :, vi][m]
            ct = d["cost_true"][:, k, :, vi][m]
            c0 = st.b["cost_final"][:, ri, :, vi][m]
            floor = np.isfinite(e3) & (e3 > 0.035)
            better_true = float(np.mean(ct[floor] < cf[floor] - 1e-12)) if floor.any() else np.nan
            e0 = st.pipe_err("H0")[:, ri, :, vi][m]
            mc = np.where(c0 < cf, e0, e3)
            rows.append(f"r={radii[ri]:g}: H3>0.035 deg in {100*floor.mean():.0f}% of cases; of those "
                        f"cost(truth)<cost(H3) in {100*better_true:.0f}%; min-cost(H0,H3) median "
                        f"{np.nanmedian(mc):.3g} vs H3 {np.nanmedian(e3):.3g}, wrong {100*np.mean(mc[np.isfinite(mc)]>1):.2f}% vs {100*np.mean(e3[np.isfinite(e3)]>1):.2f}%")  # fmt: skip
        out.append("  D2 H3 above 0.035 deg floor; cost at truth vs at result; min-cost rule (H0,H3):")
        out += ["      " + r for r in rows]
    # D3 H3c vs H3 where the covariance box exceeds the default
    if "H3c" in st.raws and "H3" in st.raws:
        c, h = st.raws["H3c"], st.raws["H3"]
        ec, eh = st.pipe_err("H3c"), st.pipe_err("H3")
        sel_all = []
        for ri in range(10):
            if ri not in c["radii_idx"]:
                continue
            k = list(c["radii_idx"]).index(ri)
            big = (c["box_deg"][:, k, :, vi] > 0) & ~c["copied"][:, k, :, vi] & st.mask[:, ri, :, vi]
            if big.sum() == 0:
                continue
            a_, b_ = ec[:, ri, :, vi][big], eh[:, ri, :, vi][big]
            w, t, n = win_rate(a_, b_)
            out.append(f"  D3 r={radii[ri]:g}: {int(big.sum())} cases with b > default; median H3c "
                       f"{np.nanmedian(a_):.3g} vs H3 {np.nanmedian(b_):.3g}; wrong "
                       f"{100*np.mean(a_>1):.2f}% vs {100*np.mean(b_>1):.2f}%; H3c win {w:.2f}")  # fmt: skip
            sel_all.append((a_, b_))
        if sel_all:
            a_ = np.concatenate([s[0] for s in sel_all])
            b_ = np.concatenate([s[1] for s in sel_all])
            w, t, n = win_rate(a_, b_)
            better = np.nanmedian(a_) < np.nanmedian(b_) - 0.002 and w > 0.5
            out.append(f"  D3 pooled ({n} cases): median H3c {np.nanmedian(a_):.3g} vs H3 "
                       f"{np.nanmedian(b_):.3g}, win {w:.2f} -> "
                       f"{'keep the covariance box' if better else 'drop the covariance box (not better)'}")  # fmt: skip
        else:
            out.append("  D3: no case has b > default box")
    else:
        out.append("  D3 not evaluated")
    # D4 H3m vs H3, time
    if "H3m" in st.raws and "H3" in st.raws:
        mm, hh = st.raws["H3m"], st.raws["H3"]
        em, eh = st.pipe_err("H3m"), st.pipe_err("H3")
        for ri in range(10):
            if ri not in mm["radii_idx"] or radii[ri] > 3:
                continue
            k = list(mm["radii_idx"]).index(ri)
            m = st.mask[:, ri, :, vi]
            tm = np.nanmean(mm["t_finish"][:, k, :, vi][m])
            th = np.nanmean(hh["t_finish"][:, k, :, vi][m])
            dmed = np.nanmedian(em[:, ri, :, vi][m]) - np.nanmedian(eh[:, ri, :, vi][m])
            out.append(f"  D4 r={radii[ri]:g}: median H3m - H3 = {dmed:+.4f} deg; finisher time "
                       f"{tm:.2f}s vs {th:.2f}s ({tm/th:.2f}x) -> "
                       f"{'H3m recommended' if abs(dmed) <= 0.01 and tm < 0.5 * th else 'no'}")  # fmt: skip
    else:
        out.append("  D4 not evaluated")
    # post hoc rules
    if "H3" in st.raws:
        out.append("  Post hoc (a): H3 falling back to HG where pass-1 sigma_max > tau (radii with HG):")
        eh, eg = st.pipe_err("H3"), st.pipe_err("HG")
        ri_h3 = list(st.raws["H3"]["radii_idx"])
        if eg is not None:
            sig1 = st.raws["H3"]["sig1"]
            for tau in TAU_GRID:
                for ri in st.raws["HG"]["radii_idx"]:
                    k = ri_h3.index(ri)
                    m = st.mask[:, ri, :, vi]
                    use = sig1[:, k, :, vi][m] > tau
                    e = np.where(use, eg[:, ri, :, vi][m], eh[:, ri, :, vi][m])
                    out.append(f"      tau={tau:g} r={radii[ri]:g}: HG used in {100*use.mean():.0f}%; "
                               f"median {np.nanmedian(e):.3g} (H3 {np.nanmedian(eh[:, ri, :, vi][m]):.3g}, "
                               f"HG {np.nanmedian(eg[:, ri, :, vi][m]):.3g}); wrong "
                               f"{100*np.mean(e>1):.2f}% (H3 {100*np.mean(eh[:, ri, :, vi][m]>1):.2f}%, "
                               f"HG {100*np.mean(eg[:, ri, :, vi][m]>1):.2f}%)")  # fmt: skip
        out.append("  Post hoc (b): min cost of H0 and H3 (all radii):")
        for ri in range(10):
            if ri not in ri_h3:
                continue
            k = ri_h3.index(ri)
            m = st.mask[:, ri, :, vi]
            cf = st.raws["H3"]["cost_final"][:, k, :, vi][m]
            c0 = st.b["cost_final"][:, ri, :, vi][m]
            e0 = st.pipe_err("H0")[:, ri, :, vi][m]
            e3 = eh[:, ri, :, vi][m]
            e = np.where(c0 < cf, e0, e3)
            out.append(f"      r={radii[ri]:g}: median {np.nanmedian(e):.3g} (H3 {np.nanmedian(e3):.3g}, "
                       f"H0 {np.nanmedian(e0):.3g}); wrong {100*np.mean(e>1):.2f}% "
                       f"(H3 {100*np.mean(e3>1):.2f}%, H0 {100*np.mean(e0>1):.2f}%); picks H0 {100*np.mean(c0<cf):.0f}%")  # fmt: skip
    return out


def checks(st: Study) -> List[str]:
    """Consistency of the live runs with the stored sweeps."""
    out = ["Consistency checks:"]
    for pipe in ("H0", "H1", "H3", "H3c", "H3m", "HG"):
        d = st.raws.get(pipe)
        if d is None:
            continue
        ri_ = list(d["radii_idx"])
        if pipe == "HG":
            out.append("  HG: nets start at the Huber GN estimate (not comparable to the sweep's N1/N3)")
            continue
        # net estimates vs the stored perturbation sweep
        for tag, key, p in (("N1", "R_x1", 0), ("N3", "R_x3", 2)):
            R = np.full((len(st.vox), 10) + d[key].shape[2:], np.nan)
            R[:, ri_] = d[key]
            e = st.err(R, reduce=False)
            ref = st.net_err_sweep(p)
            ok = np.isfinite(e) & np.isfinite(ref)
            both_nan = (np.isnan(e) == np.isnan(ref)) | ~np.isfinite(R[..., 0, 0])
            sel = np.zeros_like(ok)
            sel[:, ri_] = True
            diff = np.abs(e - ref)[ok & sel]
            out.append(f"  {pipe}: {tag} vs sweep err_angle: n={len(diff)}, max |diff| = {diff.max():.2e} deg")
        if pipe == "H0":
            R, Rb = d["R_final"], st.b["R_final"][:, ri_]
            ran = d["ran"]
            nneq = int((np.abs(R - Rb).max(axis=(-1, -2))[ran] > 0).sum())
            out.append(f"  H0: R_final identical to findoptimal_b_raw in {int(ran.sum()) - nneq}/{int(ran.sum())} cases")
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out-dir", default=str(OUT_DIR))
    ap.add_argument("--model", default="realistic_s0")
    ap.add_argument("--name", default="summary")
    ap.add_argument("--n-voxels", type=int, default=0)
    ap.add_argument("--variant", default=None, help="raw files restricted to one variant (s1 / clean_s0)")
    ap.add_argument("--pipes", nargs="*", default=["H0", "H1", "H3", "H3c", "H3m", "HG"])
    args = ap.parse_args()
    out_dir = Path(args.out_dir)
    st = Study(out_dir, args.model, args.n_voxels)
    for p in args.pipes:
        if not st.load(p, variant=args.variant) and args.variant:
            print(f"(no {p} raw for {args.model}/{args.variant})")
    res = build(st)
    lines: List[str] = [
        "Hybrid network -> FindOptimal: per-voxel images (<= 3 distractor sources), not full-sample",
        "renders. Noise pairing: H1 is the strictly paired hybrid (net pass 1 and the finisher see",
        "the same image); passes 2-3 and HG re-render. 'realistic' = the sweep variant 'all'.",
        f"net model: {args.model}; voxels {len(st.vox)}; radii {list(map(float, st.radii))}",
        "",
    ]
    lines += checks(st) + [""]
    for vi, vname in VARIANT_LABEL.items():
        n_ok = [res[vname][ri]["_n_success"] for ri in range(10)]
        n_fb = [res[vname][ri]["_n_fallback"] for ri in range(10)]
        lines += [f"=== {vname} ===", f"pass-1 successes per radius: {n_ok}; fall-backs (no net estimate): {n_fb}", ""]
        lines += [table(res, st.radii, vname, "median", "median error (deg, cubic-reduced)"), ""]
        lines += [table(res, st.radii, vname, "rms", "RMS error (deg, reduced)"), ""]
        lines += [table(res, st.radii, vname, "rms_unred", "RMS error (deg, unreduced)"), ""]
        lines += [table(res, st.radii, vname, "f01", "fraction < 0.1 deg"), ""]
        lines += [wrong_table(res, st.radii, vname), ""]
        for p in ("H3", "H1", "H3c", "H3m", "HG"):
            if any(p in res[vname][ri] for ri in range(10)):
                lines += [win_table(res, st.radii, vname, p), ""]
        lines += [table(res, st.radii, vname, "evals_mean", "mean cost evaluations per case (finisher)"), ""]
        lines += [time_table(st, vi), ""]
        lines += criteria(st, res, vname) + [""]
        lines += decisions(st, res, vname, vi) + [""]
    txt = "\n".join(lines)
    (out_dir / f"{args.name}.txt").write_text(txt)
    (out_dir / f"{args.name}.json").write_text(json.dumps(res, indent=1, default=float))
    print(txt)


if __name__ == "__main__":
    main()
