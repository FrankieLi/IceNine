"""Tables, JSON and plot for scripts/findoptimal_sweep.py: the multi-level reconstruction
(FindOptimal) next to the network, MC, Adam and GN on the perturbation sweep's cases."""

import json
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
from scipy.spatial.transform import Rotation

import findoptimal_sweep as fs
import optimizer_sweep_summary as oss
import perturbation_sweep as ps

VARIANTS = ps.VARIANTS
WRONG_DEG = 1.0
B_COLOR = "#111111"
A_COLOR = "#56B4E9"
COMPARE = ["MC", "Adam", "GN x3", "net realistic x3", "net realistic one-shot"]


def _stats(err: np.ndarray, r: float) -> Dict[str, float]:
    return oss._stats(err, r)


# ---------------------------------------------------------------------------
# Coincidence-site-lattice classification of the wrong solutions
# ---------------------------------------------------------------------------

_CSL: Optional[List[Tuple[int, float, Tuple[int, int, int]]]] = None
BRANDON_DEG = 15.0
FAMILIES = {"<100>": (1, 0, 0), "<110>": (1, 1, 0), "<111>": (1, 1, 1)}


def csl_table(max_sigma: int = 29) -> List[Tuple[int, float, Tuple[int, int, int]]]:
    """Cubic CSL misorientations with odd Sigma <= max_sigma as (Sigma, angle deg, axis uvw): for an
    axis [uvw] and integer m, tan(theta/2) = sqrt(N)/m with N = u^2+v^2+w^2 gives Sigma = m^2 + N
    (divided by 2 until odd). Axes with indices <= 3, one entry per (Sigma, reduced angle)."""
    global _CSL
    if _CSL is None:
        seen: Dict[Tuple[int, float], Tuple[int, float, Tuple[int, int, int]]] = {}
        for u in range(4):
            for v in range(u + 1):
                for w in range(v + 1):
                    if (u, v, w) == (0, 0, 0) or np.gcd.reduce([u, v, w]) != 1:
                        continue
                    N = u * u + v * v + w * w
                    for m in range(1, 60):
                        sig = m * m + N
                        while sig % 2 == 0:
                            sig //= 2
                        if sig == 1 or sig > max_sigma:
                            continue
                        th = 2 * np.degrees(np.arctan(np.sqrt(N) / m))
                        Rc = Rotation.from_rotvec(np.radians(th) * np.array([u, v, w]) / np.sqrt(N))
                        red = float(fs.misorientation_deg(Rc.as_matrix(), np.eye(3)))
                        key = (sig, round(red, 2))
                        if key not in seen:
                            seen[key] = (sig, float(th), (u, v, w))
        _CSL = sorted(seen.values())
    return _CSL


def classify_csl(R_true: np.ndarray, R_final: np.ndarray) -> Dict[str, Any]:
    """Crystal-frame misorientation M = R_true^T R_final S minimised in angle over the 24 cubic S:
    its angle, the axis folded into the cubic fundamental family (|components| sorted descending),
    the nearest of <100>/<110>/<111> and its angular distance, and the lowest-Sigma CSL
    (Sigma <= 29) within the Brandon criterion 15 deg / sqrt(Sigma) (deviation = smallest angle
    between M and the ideal CSL rotation over the cubic symmetry on both sides); Sigma 0 = none."""
    ops = fs.cubic_rotations()
    Mk = np.asarray(R_true).T @ np.asarray(R_final) @ ops
    mags = Rotation.from_matrix(Mk).magnitude()
    k = int(np.argmin(mags))
    rv = Rotation.from_matrix(Mk[k]).as_rotvec()
    axis = np.sort(np.abs(rv / np.linalg.norm(rv)))[::-1]
    dist = {
        n: float(np.degrees(np.arccos(np.clip(axis @ (np.array(f) / np.linalg.norm(f)), -1, 1))))
        for n, f in FAMILIES.items()
    }
    near = min(dist, key=lambda n: dist[n])
    both = ops[:, None] @ Mk[k] @ ops[None]  # (24, 24, 3, 3): S1 M S2
    best = None
    for sig, th, ax in csl_table():
        Rc = Rotation.from_rotvec(np.radians(th) * np.array(ax) / np.linalg.norm(ax)).as_matrix()
        dev = float(
            np.degrees(Rotation.from_matrix(both.reshape(-1, 3, 3) @ Rc.T).magnitude().min())
        )
        if dev <= BRANDON_DEG / np.sqrt(sig):
            best = (sig, th, ax, dev)
            break  # table is sorted by Sigma
    return dict(
        angle=float(np.degrees(mags[k])),
        axis=[float(x) for x in axis],
        nearest_family=near,
        nearest_family_deg=dist[near],
        sigma=best[0] if best else 0,
        sigma_angle=best[1] if best else float("nan"),
        sigma_axis=list(best[2]) if best else [],
        deviation_deg=best[3] if best else float("nan"),
        brandon_limit_deg=float(BRANDON_DEG / np.sqrt(best[0])) if best else float("nan"),
    )


# ---------------------------------------------------------------------------
# Experiment A
# ---------------------------------------------------------------------------


def _csl_counts(rows: List[Dict[str, Any]]) -> Dict[str, int]:
    out: Dict[str, int] = {}
    for r_ in rows:
        if r_["csl"] is None:
            continue
        key = f"Sigma{r_['csl']['sigma']}" if r_["csl"]["sigma"] else "no low-Sigma CSL"
        out[key] = out.get(key, 0) + 1
    return dict(sorted(out.items()))


def summarize_a(a: Dict[str, np.ndarray]) -> Dict[str, Any]:
    out: Dict[str, Any] = {}
    vox = a["voxel_indices"]
    for vi, v in enumerate(VARIANTS):
        sel = np.nonzero(a["variant_index"] == vi)[0]
        err = fs.misorientation_deg(a["R_final"][sel], a["R_true"][sel])
        plain = fs.misorientation_deg(a["R_final"][sel], a["R_true"][sel], reduce=False)
        wrong = err > WRONG_DEG
        rows = []
        lock = []
        for k, s in enumerate(sel):
            lev = [
                (
                    float(fs.misorientation_deg(R, a["R_final"][s]))
                    if np.all(np.isfinite(R))
                    else float("nan")
                )
                for R in a["level_R"][s]
            ]
            lev_err = [
                (
                    float(fs.misorientation_deg(R, a["R_true"][s]))
                    if np.all(np.isfinite(R))
                    else float("nan")
                )
                for R in a["level_R"][s]
            ]
            hit = [i for i, d in enumerate(lev) if d == d and d < fs.LOCK_ON_DEG]
            lock.append(hit[0] if hit else -1)
            csl = classify_csl(a["R_true"][s], a["R_final"][s]) if err[k] > WRONG_DEG else None
            rows.append(
                dict(
                    csl=csl,
                    voxel=int(vox[s]),
                    err=float(err[k]),
                    err_plain=float(plain[k]),
                    cost_final=float(a["cost_final"][s]),
                    cost_true=float(a["cost_true"][s]),
                    level_best_err=lev_err,
                    level_best_to_final=lev,
                    find_winner=int(a["find_winner"][s]),
                    find_evaluated=int(a["find_evaluated"][s]),
                    find_converged=bool(a["find_converged"][s]),
                    runtime=float(a["runtime"][s]),
                )
            )
        good = err[~wrong]
        out[v] = dict(
            n=int(len(sel)),
            median=float(np.median(err)),
            rms=float(np.sqrt(np.mean(err**2))),
            max=float(err.max()),
            p25=float(np.percentile(err, 25)),
            p75=float(np.percentile(err, 75)),
            frac_lt_0p1=float(np.mean(err < 0.1)),
            frac_wrong=float(np.mean(wrong)),
            n_wrong=int(wrong.sum()),
            wrong_voxels=[int(vox[s]) for s, w in zip(sel, wrong) if w],
            wrong_errors=[float(e) for e in err[wrong]],
            wrong_search_miss=int(
                sum(a["cost_final"][s] > a["cost_true"][s] + 1e-9 for s, w in zip(sel, wrong) if w)
            ),
            wrong_near_60=int(sum(1 for e in err[wrong] if abs(e - 60.0) < 2.0)),
            wrong_near_37_38=int(sum(1 for e in err[wrong] if 36.0 <= e <= 39.0)),
            wrong_other=int(
                sum(1 for e in err[wrong] if not (abs(e - 60.0) < 2.0 or 36.0 <= e <= 39.0))
            ),
            wrong_lost_after_near_truth_level=int(
                sum(
                    1
                    for r_ in rows
                    if r_["err"] > WRONG_DEG
                    and np.nanmin(r_["level_best_err"]) < 10.0
                    and r_["level_best_err"][-1] > 10.0
                )
            ),
            csl_counts=_csl_counts(rows),
            n_plain_differs_from_reduced=int((np.abs(err - plain) > 1e-6).sum()),
            median_right=float(np.median(good)) if len(good) else float("nan"),
            rms_right=float(np.sqrt(np.mean(good**2))) if len(good) else float("nan"),
            max_right=float(good.max()) if len(good) else float("nan"),
            sym_equals_plain_max_diff=float(np.abs(err - plain).max()),
            runtime_median=float(np.median(a["runtime"][sel])),
            runtime_mean=float(np.mean(a["runtime"][sel])),
            evals_global_median=float(np.median(a["evals_global"][sel])),
            evals_local_median=float(np.median(a["evals_local"][sel])),
            find_converged=int(a["find_converged"][sel].sum()),
            find_winner_rank0=int((a["find_winner"][sel] == 0).sum()),
            lock_on_level_counts={
                str(k): int(sum(1 for x in lock if x == k)) for k in range(-1, 4)
            },
            per_voxel=rows,
        )
    return out


# ---------------------------------------------------------------------------
# Experiment B
# ---------------------------------------------------------------------------


def b_errors(b: Dict[str, np.ndarray], sweep: Dict[str, np.ndarray], R_true_by_voxel: Any):
    """(err_sym, err_plain) with shape (V, R, D, 2) of the B results (NaN where not run)."""
    vox = b["voxel_indices"]
    Rt = np.stack([R_true_by_voxel(int(v)) for v in vox])[:, None, None, None]  # (V,1,1,1,3,3)
    R = b["R_final"]  # (V, R, D, 2, 3, 3)
    ok = b["ran"] & np.isfinite(R).all((-1, -2))
    Rf = np.where(ok[..., None, None], R, np.eye(3))
    e_sym = fs.misorientation_deg(Rf, Rt)
    e_pl = fs.misorientation_deg(Rf, Rt, reduce=False)
    return np.where(ok, e_sym, np.nan), np.where(ok, e_pl, np.nan)


def summarize_b(
    b: Dict[str, np.ndarray],
    sweep: Dict[str, np.ndarray],
    e_sym: np.ndarray,
    e_pl: np.ndarray,
) -> Dict[str, Any]:
    radii = [float(r) for r in b["radii"]]
    ris = [int(i) for i in b["radii_idx"]]
    vpos = {int(v): k for k, v in enumerate(sweep["voxel_indices"])}
    sel_v = [vpos[int(v)] for v in b["voxel_indices"]]
    models = [str(m) for m in sweep["models"]]
    groups = ps.seed_groups(models)
    out: Dict[str, Any] = {"radii": radii, "variants": {}}
    for vi, v in enumerate(VARIANTS):
        vo: Dict[str, Any] = {}
        for k, (ri, r) in enumerate(zip(ris, radii)):
            net_ok = (sweep["fail_pass1"][sel_v][:, ri, :, vi] == 0)[:, : e_sym.shape[2]]
            err = e_sym[:, k, :, vi]
            ran = net_ok & np.isfinite(err)
            s = _stats(err[ran], r)
            s["median_runtime_s"] = float(np.median(b["runtime"][:, k, :, vi][ran]))
            s["median_evals"] = float(np.median(b["evals"][:, k, :, vi][ran]))
            s["frac_converged_hit_ratio"] = float(np.mean(b["converged"][:, k, :, vi][ran]))
            s["frac_cost_ge_1_identity"] = float(np.mean(b["cost_final"][:, k, :, vi][ran] >= 1.0))
            s["frac_cost_below_truth_cost"] = float(
                np.mean(b["cost_final"][:, k, :, vi][ran] < b["cost_true"][:, k, :, vi][ran])
            )
            s["frac_worse_than_start"] = float(np.mean(err[ran] > r))
            fb = np.where(b["cost_final"][:, k, :, vi] >= 1.0, r, err)  # identity -> keep the input
            s["median_with_fallback"] = float(np.median(fb[ran]))
            s["rms_with_fallback"] = float(np.sqrt(np.mean(fb[ran] ** 2)))
            s["frac_lt_0p1_with_fallback"] = float(np.mean(fb[ran] < 0.1))
            s["median_plain"] = float(np.median(e_pl[:, k, :, vi][ran]))
            s["max_sym_minus_plain_abs"] = float(np.max(np.abs(err[ran] - e_pl[:, k, :, vi][ran])))
            paired: Dict[str, Any] = {}
            for tname, idx in groups.items():
                for p, tag in ((0, "one-shot"), (2, "x3")):
                    wins, ns = [], []
                    for i in idx:
                        ne = sweep["err_angle"][sel_v][:, ri, :, vi, i, p][:, : e_sym.shape[2]]
                        okp = ran & np.isfinite(ne)
                        ns.append(int(okp.sum()))
                        wins.append(float(np.mean(ne[okp] < err[okp])) if okp.any() else np.nan)
                    paired[f"net {tname} {tag} < FindOptimal"] = dict(
                        n=ns[0], win_per_seed=wins, win_mean=float(np.mean(wins))
                    )
            s["paired"] = paired
            vo[f"{r:g}"] = s
        out["variants"][v] = vo
    return out


# ---------------------------------------------------------------------------
# Symmetry check on the other methods
# ---------------------------------------------------------------------------


def symmetry_check(
    sweep: Dict[str, np.ndarray], opt: Dict[str, np.ndarray], R_true_by_voxel: Any, n: int = 20000
) -> Dict[str, Any]:
    """Symmetry-reduced vs plain error of the net / MC / Adam / GN / Huber results, on random cases
    (R_est = exp(error vector) R_true, error vectors are stored in the sweeps)."""
    rng = np.random.default_rng(0)
    vox = sweep["voxel_indices"]
    Rt = np.stack([R_true_by_voxel(int(v)) for v in vox])
    res: Dict[str, Any] = {}

    def check(name: str, ex: np.ndarray, ey: np.ndarray, ez: np.ndarray, ang: np.ndarray) -> None:
        fin = np.argwhere(np.isfinite(ang) & np.isfinite(ex))
        pick = fin[rng.choice(len(fin), size=min(n, len(fin)), replace=False)]
        e = np.stack([ex[tuple(pick.T)], ey[tuple(pick.T)], ez[tuple(pick.T)]], -1)
        R = Rotation.from_rotvec(np.radians(e)).as_matrix() @ Rt[pick[:, 0]]
        sym = fs.misorientation_deg(R, Rt[pick[:, 0]])
        pl = fs.misorientation_deg(R, Rt[pick[:, 0]], reduce=False)
        res[name] = dict(
            n_checked=int(len(pick)),
            max_plain_error_all_cases_deg=float(np.nanmax(ang)),
            max_abs_sym_minus_plain=float(np.abs(sym - pl).max()),
            n_sym_smaller=int((sym < pl - 1e-6).sum()),
        )

    models = [str(m) for m in sweep["models"]]
    for i, m in enumerate(models):
        ex, ey, ez, an = (
            sweep[k][:, :, :, :, i, :] for k in ("err_x", "err_y", "err_z", "err_angle")
        )
        # (V, R, D, 2, P) -> index 0 is the voxel axis, as Rt[pick[:, 0]] expects
        check(f"net {m}", ex, ey, ez, an)
    for mi, m in enumerate(str(x) for x in opt["methods"]):
        ex, ey, ez, an = (
            opt[k][:, :, :, :, mi, :] for k in ("err_x", "err_y", "err_z", "err_angle")
        )
        check(f"opt {m}", ex, ey, ez, an)
    return res


# ---------------------------------------------------------------------------
# Output
# ---------------------------------------------------------------------------


def _f(x: float, f: str = ".4f") -> str:
    return "   nan" if x != x else format(x, f)


def format_text(
    a: Dict[str, Any], b: Dict[str, Any], opt_summ: Dict[str, Any], symc: Dict[str, Any],
    n_vox: int, n_dirs: int, cfg: Dict[str, Any],
) -> str:  # fmt: skip
    L: List[str] = []
    L.append(
        "FindOptimal (AdaptiveVoxelReconstructor) on the perturbation-sweep cases. Settings: "
        "ReconstructQ8.config SearchParameters (5 deg grid radius, 4 levels, 200 MC steps, 2 "
        "restarts, MC radius scale 0.4, convergence cost 1e-4, 30 candidates), Q_max 8, "
        f"|sin eta| >= {cfg['min_sin_eta']:g}. Error = cubic-symmetry-reduced misorientation (deg)."
    )
    L.append("")
    L.append(
        f"=== A: full multi-level reconstruction from scratch ({n_vox} voxels per variant) ==="
    )
    for v in VARIANTS:
        s = a[v]
        L.append(
            f"[{v}] n={s['n']}  median {s['median']:.4f}  RMS {s['rms']:.4f}  max {s['max']:.3f}  "
            f"IQR {s['p25']:.4f}-{s['p75']:.4f}  <0.1deg {s['frac_lt_0p1']:.3f}  "
            f"wrong (>{WRONG_DEG:g} deg) {s['n_wrong']}/{s['n']}"
        )
        L.append(
            f"    excluding the wrong ones: median {s['median_right']:.4f}  RMS {s['rms_right']:.4f}  "
            f"max {s['max_right']:.4f}"
        )
        L.append(
            f"    runtime/voxel: median {s['runtime_median']:.1f} s, mean {s['runtime_mean']:.1f} s; "
            f"cost evals: global (pixel_radius 3) median {s['evals_global_median']:.0f}, "
            f"local median {s['evals_local_median']:.0f}; FindOptimal converged (hit ratio) "
            f"{s['find_converged']}/{s['n']}, winner = best candidate {s['find_winner_rank0']}/{s['n']}"
        )
        L.append(
            f"    first level whose best candidate is within {fs.LOCK_ON_DEG:g} deg of the final "
            f"answer (-1 = none): {s['lock_on_level_counts']}"
        )
        L.append(
            f"    wrong solutions: ~60 deg (Sigma3-twin-like) {s['wrong_near_60']}, 36-39 deg "
            f"{s['wrong_near_37_38']}, other {s['wrong_other']}; final cost above the cost at the "
            f"truth (search miss, not a cost preference): {s['wrong_search_miss']}/{s['n_wrong']}; "
            f"an earlier level's best candidate was within 10 deg of the truth but the last was "
            f"not: {s['wrong_lost_after_near_truth_level']}; plain error differs from the "
            f"symmetry-reduced one (equivalent orientation returned): "
            f"{s['n_plain_differs_from_reduced']}/{s['n']}"
        )
        L.append(f"    CSL classification of the wrong solutions: {s['csl_counts']}")
        for row in s["per_voxel"]:
            if row["err"] > WRONG_DEG:
                c = row["csl"]
                sig = (
                    f"Sigma{c['sigma']} (dev {c['deviation_deg']:.2f} deg <= Brandon "
                    f"{c['brandon_limit_deg']:.2f})"
                    if c["sigma"]
                    else "no CSL Sigma<=29 within Brandon"
                )
                L.append(
                    f"      wrong: voxel {row['voxel']} err {row['err']:.2f} deg (plain "
                    f"{row['err_plain']:.2f}) cost {row['cost_final']:.3f} vs truth "
                    f"{row['cost_true']:.3f}; axis {np.round(c['axis'], 3).tolist()} is "
                    f"{c['nearest_family_deg']:.1f} deg from {c['nearest_family']}; {sig}; "
                    f"level-best errors {[round(x, 2) for x in row['level_best_err']]}"
                )
    L.append("")
    L.append(
        f"=== B: FindOptimal + VarianceMinimizing from the perturbed start ({n_dirs} directions "
        "per radius; cases where the net could run) ==="
    )
    radii = [f"{r:g}" for r in b["radii"]]
    for v in VARIANTS:
        L.append(f"--- {v} data ---")
        vo, oo = b["variants"][v], opt_summ["variants"][v]
        blocks = [
            ("median error (deg)", "median", ".4f"),
            ("RMS error (deg)", "rms", ".4f"),
            ("fraction < 0.1 deg", "frac_lt_0p1", ".3f"),
            ("fraction improved (error < r)", "frac_improved", ".3f"),
        ]
        for title, key, fm in blocks:
            L.append(f"-- {title}")
            L.append(
                f"{'r':>5} {'n':>5} {'FindOptimal':>12} " + " ".join(f"{c:>22}" for c in COMPARE)
            )
            for r in radii:
                L.append(
                    f"{r:>5} {vo[r]['n']:>5} {_f(vo[r][key], fm):>12} "
                    + " ".join(f"{_f(oo[r][c][key], fm):>22}" for c in COMPARE)
                )
            L.append("")
        L.append(
            "-- FindOptimal if a result with cost >= 1 (no overlap; reconstruct_voxel then returns "
            "the identity matrix) is replaced by the unchanged input (error = r)"
        )
        L.append(f"{'r':>5} {'median':>9} {'RMS':>9} {'<0.1deg':>9}")
        for r in radii:
            s = vo[r]
            L.append(
                f"{r:>5} {s['median_with_fallback']:>9.4f} {s['rms_with_fallback']:>9.4f} "
                f"{s['frac_lt_0p1_with_fallback']:>9.3f}"
            )
        L.append("")
        L.append("-- FindOptimal details per radius")
        L.append(
            f"{'r':>5} {'runtime(s)':>10} {'evals':>7} {'hitratio-conv':>14} {'cost>=1(identity)':>18} "
            f"{'cost<truth':>11} {'worse>r':>8} {'max|sym-plain|':>15}"
        )
        for r in radii:
            s = vo[r]
            L.append(
                f"{r:>5} {s['median_runtime_s']:>10.3f} {s['median_evals']:>7.0f} "
                f"{s['frac_converged_hit_ratio']:>14.3f} {s['frac_cost_ge_1_identity']:>18.3f} "
                f"{s['frac_cost_below_truth_cost']:>11.3f} {s['frac_worse_than_start']:>8.3f} "
                f"{s['max_sym_minus_plain_abs']:>15.2g}"
            )
        L.append(
            "  (runtime and evals: median per case, refine_from_candidates only; MC median 2.4 s)"
        )
        L.append("")
        for tag in ("x3", "one-shot"):
            L.append(
                f"-- paired: fraction of identical cases where the realistic-trained net ({tag}) "
                "has the smaller error than FindOptimal (seed mean)"
            )
            L.append("   " + " ".join(f"{r:>6}" for r in radii))
            L.append(
                "   "
                + " ".join(
                    f"{_f(vo[r]['paired'][f'net realistic {tag} < FindOptimal']['win_mean'], '.3f'):>6}"
                    for r in radii
                )
            )
            L.append("")
    L.append(
        "=== symmetry reduction on the sweep's other methods (random cases, R_est = exp(err) R_true) ==="
    )
    for k, s in symc.items():
        L.append(
            f"  {k:>24}: n {s['n_checked']:>6}  max plain error over all cases "
            f"{s['max_plain_error_all_cases_deg']:.3f} deg  max |sym - plain| "
            f"{s['max_abs_sym_minus_plain']:.2g}  cases where reduction lowers the error "
            f"{s['n_sym_smaller']}"
        )
    return "\n".join(L)


def plot(
    b: Dict[str, Any], a: Dict[str, Any], opt_summ: Dict[str, Any], path: Path, n_vox: int
) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    radii = np.array(b["radii"])
    keys = [f"{r:g}" for r in radii]
    fig, axes = plt.subplots(1, 2, figsize=(11.5, 5.0), sharey=True)
    grid, ink = "#d9d9d9", "#333333"
    for ax, v in zip(axes, VARIANTS):
        oo = opt_summ["variants"][v]
        ax.plot(radii, radii, ":", color="0.45", lw=1.2, label="y = x (no correction)")
        ax.axvline(1.0, color="0.7", lw=1.2)
        for label, (c, mk, ls) in oss.PLOT_STYLE.items():
            ax.plot(radii, [oo[k][label]["median"] for k in keys], ls, color=c, marker=mk,
                    ms=4, lw=1.4, alpha=0.8, label=label)  # fmt: skip
        for tag, ls in (("x3", "-"), ("one-shot", "--")):
            lab = "net, realistic-trained, " + ("iterated x3" if tag == "x3" else "one-shot")
            ax.plot(radii, [oo[k][f"net realistic {tag}"]["median"] for k in keys], ls,
                    color=oss.NET_COLOR, lw=2.0, marker="o", ms=3.5, label=lab)  # fmt: skip
        ax.plot(radii, [b["variants"][v][k]["median"] for k in keys], "-", color=B_COLOR, lw=2.4,
                marker="D", ms=5, label="FindOptimal from the perturbed start (B)")  # fmt: skip
        sa = a[v]
        ax.axhline(sa["median"], color=A_COLOR, lw=2.2, ls="-",
                   label="full multi-level reconstruction (A), median over all voxels")  # fmt: skip
        ax.axhline(sa["median_right"], color=A_COLOR, lw=1.6, ls="--",
                   label="A, median over the voxels it got right (>1 deg = wrong)")  # fmt: skip
        ax.text(0.98, 0.03, f"A: {sa['n_wrong']}/{sa['n']} voxels wrong (> {WRONG_DEG:g} deg)",
                transform=ax.transAxes, ha="right", va="bottom", fontsize=8, color=ink)  # fmt: skip
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlabel("perturbation radius r (deg)")
        ax.set_title("clean data" if v == "clean" else "realistic data (all)", color=ink)
        ax.grid(True, which="both", color=grid, lw=0.6)
        ax.set_axisbelow(True)
    axes[0].set_ylabel("median angular error (deg; symmetry-reduced for FindOptimal)")
    axes[0].legend(fontsize=7, loc="upper left", frameon=False)
    fig.suptitle(
        f"Multi-level reconstruction (FindOptimal) vs the network and optimizers: {n_vox} voxels; "
        "A has no starting guess (flat in r)",
        fontsize=9,
    )
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)


def do_summarize(out_dir: Path) -> None:
    from icenine.config_file import ConfigFile
    from icenine.mic_file import MicFile
    from generate_toy_orientation_dataset import example_dir_for

    a_raw = dict(np.load(out_dir / "findoptimal_a_raw.npz", allow_pickle=False))
    b_raw = dict(np.load(out_dir / "findoptimal_b_raw.npz", allow_pickle=False))
    sweep = dict(np.load(out_dir / "perturbation_sweep_raw.npz", allow_pickle=False))
    opt = dict(np.load(out_dir / "optimizer_sweep_raw.npz", allow_pickle=False))
    opt_summ = json.loads((out_dir / "optimizer_sweep_summary.json").read_text())
    cfg = json.loads(str(a_raw["config_json"]))
    ex = example_dir_for(cfg["example"])
    mic = MicFile.read(
        str(
            ex
            / ConfigFile.from_file(
                str(ex / "ConfigFiles/Example2.Simulation.config")
            ).sample_filename
        )
    )

    def R_true(v: int) -> np.ndarray:
        return np.asarray(mic.voxels[v].orientation, dtype=np.float64)

    a = summarize_a(a_raw)
    e_sym, e_pl = b_errors(b_raw, sweep, R_true)
    b = summarize_b(b_raw, sweep, e_sym, e_pl)
    symc = symmetry_check(sweep, opt, R_true)
    n_dirs = int(b_raw["n_dirs_run"])
    n_vox = len(b_raw["voxel_indices"])
    summ = dict(A=a, B=b, symmetry_check=symc, n_dirs_B=n_dirs)
    (out_dir / "findoptimal_sweep_summary.json").write_text(json.dumps(summ, indent=1))
    text = format_text(a, b, opt_summ, symc, n_vox, n_dirs, cfg)
    (out_dir / "findoptimal_sweep_summary.txt").write_text(text + "\n")
    plot(b, a, opt_summ, out_dir / "perturbation_sweep_vs_findoptimal.png", n_vox)
    print(text)
