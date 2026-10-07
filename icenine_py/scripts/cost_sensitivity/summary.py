#!/usr/bin/env python3
"""
Phase A summary: A1 landscape (landscape.py cache), A2 resolution (resolution.py cache) and A3 (T5
cache + A1) -> benchmarks/cost_sensitivity/{summary.json,tables.md}.

Usage (from icenine_py/):
  uv run python scripts/cost_sensitivity/summary.py
"""

import glob
import json
import math
import sys
from pathlib import Path
from typing import Any, Dict, List, Tuple

import numpy as np
from scipy import stats as sps

HERE = Path(__file__).resolve().parent
ICENINE_PY = HERE.parents[1]
sys.path.insert(0, str(HERE.parent / "common"))

from doc_tables import markdown_table, write_tables  # noqa: E402
from stats import mcnemar_exact, wilson  # noqa: E402

OUT_DIR = ICENINE_PY / "benchmarks" / "cost_sensitivity"
T5_CACHE = ICENINE_PY / "scripts" / "finisher_diagnosis" / "cache" / "full"
A1_CACHE = HERE / "cache" / "run"
VARIANT_NAMES = ["clean", "realistic"]
TOL = 1e-9
MIN_PTS = 20  # plateau points needed for a shape estimate
R_MAX_DEG = 0.1  # plateau region: radii on the grid up to this
GN_FLOOR_DEG = 0.0127  # centroid Gauss-Newton, clean (October perturbation sweep)
LOWEST_K = 10
# Median MC trial rotation angle / step, measured by drawing 4000 trials with
# QuaternionGrid.get_near_identity_point at radius tan(step)/sqrt(12) (Phase A session): the median
# angle was 0.1286 deg at a step of 0.1317 deg and 0.0096 deg at 0.01 deg.
TRIAL_PER_STEP = 0.976


def q3(x: Any, ratio: bool = False) -> List[float]:
    if ratio:  # x: list of ascending principal sigmas -> largest / smallest
        x = [v[-1] / v[0] for v in x]
    x = np.asarray(x, dtype=float)
    x = x[np.isfinite(x)]
    return [float(v) for v in np.quantile(x, [0.25, 0.5, 0.75])] if x.size else [math.nan] * 3


def fq(q: List[float], nd: int) -> str:
    """q25/q50/q75 as text."""
    return "/".join(f"{x:.{nd}f}" for x in q)


def frac(k: int, n: int) -> Dict[str, Any]:
    lo, hi = wilson(int(k), int(n))
    return dict(k=int(k), n=int(n), p=k / max(n, 1), lo=lo, hi=hi)


def fmt_frac(f: Dict[str, Any]) -> str:
    return f"{f['k']}/{f['n']} ({100 * f['p']:.0f}%, {100 * f['lo']:.0f}-{100 * f['hi']:.0f})"


def plateau_stats(
    cost: np.ndarray, dirs: np.ndarray, radii: np.ndarray, cost_true: float, eps: float
) -> Dict[str, Any]:
    """Plateau {cost <= cost_true + eps} of one case. cost (n_rad, n_dir) over radii (deg, <= R_MAX
    only), dirs (n_dir, 3) unit. Returns fraction of directions inside per radius, r50 (smallest
    grid radius where < half the directions are inside; inf if none), r_any (largest radius with
    any member), contiguous extent per direction, and the anisotropy of the plateau points."""
    inside = cost <= cost_true + eps + TOL  # (n_rad, n_dir)
    frac_r = inside.mean(axis=1)
    below = np.nonzero(frac_r < 0.5)[0]
    r50 = float(radii[below[0]]) if below.size else math.inf
    anyr = np.nonzero(inside.any(axis=1))[0]
    r_any = float(radii[anyr[-1]]) if anyr.size else 0.0
    # contiguous extent per direction: largest radius with every smaller grid radius also inside
    run = np.cumprod(inside, axis=0).astype(bool)
    ext = np.where(run.any(axis=0), radii[np.maximum(run.sum(axis=0) - 1, 0)], 0.0)
    pts = (radii[:, None, None] * dirs[None])[inside]  # (n, 3), deg
    aniso = math.nan
    if len(pts) >= MIN_PTS:
        M = pts.T @ pts / len(pts)  # second-moment matrix of the plateau points about the truth
        w = np.linalg.eigvalsh(M)
        aniso = float(math.sqrt(w[-1] / max(w[0], 1e-300)))
    return dict(frac_r=frac_r, r50=r50, r_any=r_any, ext_med=float(np.median(ext)), aniso=aniso,
                n_pts=int(len(pts)))  # fmt: skip


def min_set_stats(
    flat: np.ndarray, rad: np.ndarray, pts: np.ndarray, cost_true: np.ndarray
) -> Dict[str, Any]:
    """Tie-safe statistics of the sampled minimum. flat, rad (case, n) sample costs and radii
    (deg), pts (case, n, 3) sample offsets (deg), cost_true (case,). Costs are quantised, so many
    samples tie at the minimum; everything uses the whole set of samples at the minimum cost."""
    n_case = flat.shape[0]
    cmin = flat.min(axis=1)
    lower = cmin < cost_true - TOL
    at_min = flat <= cmin[:, None] + TOL
    cen = np.stack([pts[c][at_min[c]].mean(axis=0) for c in range(n_case)])
    rmin = np.array([rad[c][at_min[c]].min() for c in range(n_case)])
    in_pl = flat <= cost_true[:, None] + TOL
    return dict(
        n_lower=int(lower.sum()),
        min_drop_among_lower_quartiles=q3((cost_true - cmin)[lower]),
        n_tied_at_min_quartiles=q3(at_min.sum(axis=1)),
        min_radius_smallest_tied_quartiles=q3(rmin),
        plateau_n_samples_quartiles=q3(in_pl.sum(axis=1)),
        centroid_offset=np.linalg.norm(cen, axis=1),
    )


def monotone_rays(cost_true: np.ndarray, cost: np.ndarray, grid: np.ndarray) -> Dict[str, Any]:
    """Per-ray monotonicity. cost (case, n_rad, n_dir) over the ascending radii grid; a ray is
    monotone if its cost never decreases from the truth (radius 0) outward, or (from_0p005) from the
    0.005 deg sample outward."""
    n_case = cost.shape[0]
    full = np.concatenate([np.repeat(cost_true[:, None, None], cost.shape[2], 2), cost], axis=1)
    mono = (np.diff(full, axis=1) >= -TOL).all(axis=1)  # (case, dir)
    i5 = int(np.searchsorted(grid, 0.005 - 1e-12))
    mono5 = (np.diff(cost[:, i5:, :], axis=1) >= -TOL).all(axis=1)
    return dict(
        all_radii=float(mono.mean()), from_0p005=float(mono5.mean()),
        cases_all_rays_all_radii=frac(int(mono.all(axis=1).sum()), n_case),
        cases_all_rays_from_0p005=frac(int(mono5.all(axis=1).sum()), n_case),
    )  # fmt: skip


def lower_than_result(
    flat: np.ndarray, rad: np.ndarray, cost_res: np.ndarray, ang_res: np.ndarray
) -> Dict[str, Any]:
    """Among sampled points with cost below the finisher result's, how many lie farther from the
    truth than the result. Only cases whose result is within R_MAX_DEG of the truth are usable."""
    usable = ang_res < R_MAX_DEG
    n_any = n_far_cases = tot = far = 0
    per_case = []
    for c in np.nonzero(usable)[0]:
        lo = flat[c] < cost_res[c] - TOL
        if not lo.any():
            continue
        f = int((lo & (rad[c] > ang_res[c])).sum())
        n_any += 1
        tot += int(lo.sum())
        far += f
        n_far_cases += int(f > 0)
        per_case.append(f / int(lo.sum()))
    return dict(
        n_usable=int(usable.sum()), cases_with_lower=n_any, points=tot, points_farther=far,
        cases_with_farther=n_far_cases, per_case_fraction_farther_quartiles=q3(np.array(per_case)),
    )  # fmt: skip


def load_a1() -> Dict[str, Dict[str, np.ndarray]]:
    """Per pipe: arrays stacked over cases; first axis variant for variant-dependent entries."""
    out: Dict[str, Dict[str, Any]] = {}
    for pipe in ("H3", "H0"):
        parts = []
        for f in sorted(glob.glob(str(A1_CACHE / f"{pipe}_v*_r*.npz"))):
            d = np.load(f)
            t5 = np.load(T5_CACHE / Path(f).name)
            assert (t5["dirs"] == d["dirs"]).all()
            parts.append((d, t5))

        def cat(k: str, src: int = 0) -> np.ndarray:
            """Concatenate over cases: axis 1 for the landscape arrays (axis 0 is the variant)."""
            ax = 1 if src == 0 and k in ("cost", "dir", "cost_true", "cost_res", "ang_res") else 0
            return np.concatenate([p[src][k] for p in parts], axis=ax)

        rec = dict(
            cost=cat("cost"), dir=cat("dir"), cost_true=cat("cost_true"), cost_res=cat("cost_res"),
            ang_res=cat("ang_res"), radii=parts[0][0]["radii"],
            step_deg=float(parts[0][0]["step_deg"]),
            vidx=np.concatenate([np.full(len(p[0]["dirs"]), int(p[0]["vidx"])) for p in parts]),
            ri=np.concatenate([np.full(len(p[0]["dirs"]), int(p[0]["ri"])) for p in parts]),
        )  # fmt: skip
        keys = ("cost_true", "cost_res", "ang_res_true", "c_vm_smallbox_cost", "c_vm_smallbox_ang")
        keys += ("c_mc_smallstep_cost", "c_mc_smallstep_ang", "ang_start_true", "fo_log")
        keys += ("c_vm_smallbox_vm",)
        for k in keys:
            rec["t5_" + k] = cat(k, 1)
        out[pipe] = rec
    return out


def summarise_a1(rec: Dict[str, Any]) -> Dict[str, Any]:
    radii_all = rec["radii"]
    grid = radii_all[:-1]  # the 0.0005..0.1 deg grid; the last column is the initial MC step
    assert grid.max() <= R_MAX_DEG + 1e-12
    n_case = rec["cost"].shape[1]
    res: Dict[str, Any] = dict(
        n=int(n_case), step_deg=rec["step_deg"], radii=[float(r) for r in grid]
    )
    # sanity: the realistic cost at the truth equals the T5 value
    res["cost_true_vs_t5_maxabs"] = float(np.abs(rec["cost_true"][1] - rec["t5_cost_true"]).max())
    res["cost_res_vs_t5_maxabs"] = float(np.abs(rec["cost_res"][1] - rec["t5_cost_res"]).max())
    for vi, vname in enumerate(VARIANT_NAMES):
        cost = rec["cost"][vi]  # (case, n_rad+1, n_dir)
        ct = rec["cost_true"][vi]
        dirs = rec["dir"][vi]
        d_step = cost[:, -1, :] - ct[:, None]  # change at the finisher's step radius
        eps_init = np.median(np.abs(d_step), axis=1)  # per case, at the initial MC step
        # at the step the finisher's MC had reached when it stopped (trial = TRIAL_PER_STEP x step)
        dc_med = np.median(np.abs(cost[:, :-1, :] - ct[:, None, None]), axis=2)  # (case, nr)
        th_f = np.maximum(rec["t5_fo_log"][:, 5] * TRIAL_PER_STEP, 1e-6)
        eps_final = np.array(
            [np.exp(np.interp(np.log(th_f[c]), np.log(grid), np.log(np.maximum(dc_med[c], 1e-9))))
             for c in range(n_case)]
        )  # fmt: skip
        acc = rec["t5_fo_log"][:, 2] > 0  # MC accepted >= once, so its step was reduced
        v: Dict[str, Any] = dict(
            n_mc_accepted=int(acc.sum()),
            eps_init_quartiles=q3(eps_init),
            eps_final_quartiles=q3(eps_final[acc]),
            final_step_trial_deg_quartiles=q3(th_f[acc]),
            final_step_below_grid=frac(int((th_f[acc] < grid[0]).sum()), int(acc.sum())),
            step_change_signed_med=float(np.median(d_step)),
            step_frac_not_above=float(np.mean(d_step <= TOL)),
        )
        # median |cost change| and plateau fraction by radius
        dc = cost[:, :-1, :] - ct[:, None, None]
        v["abs_change_med_by_radius"] = [
            float(np.median(np.abs(dc[:, k, :]))) for k in range(len(grid))
        ]
        v["frac_le_truth_by_radius"] = [
            float(np.mean(dc[:, k, :] <= TOL)) for k in range(len(grid))
        ]
        for ename in ("eps0", "epsfinal", "epsinit"):
            ps_list = []
            for c in range(n_case):
                eps = {"eps0": 0.0, "epsfinal": eps_final[c], "epsinit": eps_init[c]}[ename]
                ps_list.append(plateau_stats(cost[c, :-1], dirs[c], grid, float(ct[c]), eps))
            keep = acc if ename == "epsfinal" else np.ones(n_case, dtype=bool)
            ps_list = [p for p, k in zip(ps_list, keep) if k]
            r50 = np.array([p["r50"] for p in ps_list])
            r_any = np.array([p["r_any"] for p in ps_list])
            ext = np.array([p["ext_med"] for p in ps_list])
            an = np.array([p["aniso"] for p in ps_list])
            v[ename] = dict(
                n_used=int(keep.sum()),
                r50_quartiles=q3(np.where(np.isfinite(r50), r50, R_MAX_DEG * 1.5)),
                r50_ge_max=int(np.sum(~np.isfinite(r50))),
                r_any_quartiles=q3(r_any),
                ext_med_quartiles=q3(ext),
                no_plateau_member=int(np.sum([p["n_pts"] == 0 for p in ps_list])),
                n_aniso=int(np.isfinite(an).sum()),
                aniso_quartiles=q3(an),
                r50=r50.tolist(), r_any=r_any.tolist(), ext_med=ext.tolist(), aniso=an.tolist(),
            )  # fmt: skip
        # sampled minimum (all grid radii <= R_MAX). Costs are quantised, so many samples tie at the
        # minimum: use tie-safe statistics (the set of ALL samples at the minimum cost).
        flat = cost[:, :-1, :].reshape(n_case, -1)
        rad = np.broadcast_to(grid[None, :, None], cost[:, :-1, :].shape).reshape(n_case, -1)
        pts = (grid[None, :, None, None] * dirs[:, None, :, :]).reshape(n_case, -1, 3)  # deg
        ms = min_set_stats(flat, rad, pts, ct)
        v["min_below_truth"] = frac(ms["n_lower"], n_case)
        v.update({k: x for k, x in ms.items() if k != "centroid_offset"})
        v["min_set_centroid_offset"] = ms["centroid_offset"].tolist()
        v["min_set_centroid_offset_quartiles"] = q3(ms["centroid_offset"])
        v["monotone_rays"] = monotone_rays(ct, cost[:, :-1, :], grid)
        # A3: within-case Spearman of cost vs angle
        rho_all, rho_fine = [], []
        for c in range(n_case):
            rho_all.append(sps.spearmanr(rad[c], flat[c])[0])
            m = rad[c] <= 0.03 + 1e-12
            rho_fine.append(sps.spearmanr(rad[c][m], flat[c][m])[0])
        v["spearman_quartiles"] = q3(np.array(rho_all))
        v["spearman_le0p03_quartiles"] = q3(np.array(rho_fine))
        v["spearman_lt0_count"] = int(np.sum(np.array(rho_all) < 0))
        # points with lower cost than the T5 finisher result: how many are farther than it?
        v["lower_than_result"] = lower_than_result(
            flat, rad, rec["cost_res"][vi], rec["ang_res"][vi]
        )
        res[vname] = v
    # paired realistic - clean on the tie-safe minimum-set centroid offset
    a_ = np.array(res["clean"]["min_set_centroid_offset"])
    b_ = np.array(res["realistic"]["min_set_centroid_offset"])
    nb, nc = int((b_ > a_ + TOL).sum()), int((b_ < a_ - TOL).sum())
    res["paired_min_set_centroid_offset"] = dict(
        real_gt_clean=nb, clean_gt_real=nc, ties=int(n_case - nb - nc),
        p_sign=mcnemar_exact(nb, nc), median_diff=float(np.median(b_ - a_)),
    )  # fmt: skip
    res["n_voxels"] = int(len(np.unique(rec["vidx"])))
    return res


def summarise_a3_t5(rec: Dict[str, Any]) -> Dict[str, Any]:
    ar = rec["t5_ang_res_true"]
    gap = rec["t5_cost_res"] - rec["t5_cost_true"]
    out: Dict[str, Any] = dict(n=int(len(ar)), ang_res_quartiles=q3(ar), gap_quartiles=q3(gap))
    out["n_voxels"] = int(len(np.unique(rec["vidx"])))
    for key, nm in (("c_vm_smallbox", "vm_smallbox"), ("c_mc_smallstep", "mc_smallstep")):
        ang2 = rec[f"t5_{key}_ang"]
        cost2 = rec[f"t5_{key}_cost"]
        closed = (rec["t5_cost_res"] - cost2) / np.where(gap > TOL, gap, np.nan)
        better = ang2 < ar - 1e-6
        worse = ang2 > ar + 1e-6
        w = sps.wilcoxon(ar, ang2) if np.any(ar != ang2) else None
        rho = sps.spearmanr(closed[np.isfinite(closed)], (ar - ang2)[np.isfinite(closed)])[0]
        out[nm] = dict(
            ang_after_quartiles=q3(ang2), closed_median=float(np.nanmedian(closed)),
            error_reduced=frac(int(better.sum()), len(ar)),
            error_increased=frac(int(worse.sum()), len(ar)),
            wilcoxon_p=float(w.pvalue) if w is not None else math.nan,
            spearman_closed_vs_error_drop=float(rho),
            ratio_med=float(np.median(ang2 / ar)),
        )  # fmt: skip
        if nm == "vm_smallbox":
            out[nm]["step_cap_hit"] = frac(int(rec["t5_c_vm_smallbox_vm"][:, 10].sum()), len(ar))
        # voxel-clustered sign test on the error drop: per-voxel median of (before - after)
        drop = ar - ang2
        vd = np.array([np.median(drop[rec["vidx"] == u]) for u in np.unique(rec["vidx"])])
        out[nm]["voxel_clustered"] = dict(
            n_voxels=int(len(vd)), voxels_improved=int((vd > 0).sum()),
            voxels_worse=int((vd < 0).sum()),
            p_sign=mcnemar_exact(int((vd > 0).sum()), int((vd < 0).sum())),
        )  # fmt: skip
    return out


def summarise_a2(a1: Dict[str, Dict[str, Any]]) -> Dict[str, Any]:
    res = json.loads((HERE / "cache" / "resolution.json").read_text())
    by_v = {r["net_filter"]["vidx"]: r for r in res}
    out: Dict[str, Any] = dict(n_voxels=len(res))
    for tag in ("net_filter", "all_peaks"):
        rows = [r[tag] for r in res]
        out[tag] = dict(
            n_peaks_quartiles=q3([r["n_peaks"] for r in rows]),
            rms3_quartiles=q3([r["both"]["rms3_deg"] for r in rows]),
            frame_only_quartiles=q3([r["frame_only"]["rms3_deg"] for r in rows]),
            pixel_only_quartiles=q3([r["pixel_only"]["rms3_deg"] for r in rows]),
            principal_min_quartiles=q3([r["both"]["principal_deg"][0] for r in rows]),
            principal_max_quartiles=q3([r["both"]["principal_deg"][-1] for r in rows]),
            aniso_quartiles=q3([r["both"]["principal_deg"] for r in rows], ratio=True),
            sin_eta_med=q3([r["median_abs_sin_eta"] for r in rows]),
            two_theta_med_deg=q3([r["median_two_theta_deg"] for r in rows]),
        )  # fmt: skip
    r0 = res[0]["net_filter"]
    out["frame_width_deg"] = r0["frame_width_deg"]
    out["pixel_size"] = r0["pixel_size"]
    out["det_dist"] = r0["det_dist"]
    for pipe, rec in a1.items():
        ratio = np.array([rec["t5_ang_res_true"][i] / by_v[int(v)]["net_filter"]["both"]["rms3_deg"]
                          for i, v in enumerate(rec["vidx"])])  # fmt: skip
        out[f"{pipe}_err_over_bound_quartiles"] = q3(ratio)
        out[f"{pipe}_err_below_bound"] = frac(int((ratio < 1).sum()), len(ratio))
    return out


def tables(S: Dict[str, Any]) -> Dict[str, str]:
    t: Dict[str, str] = {}
    rows = []
    for pipe in ("H3", "H0"):
        for vn in VARIANT_NAMES:
            v = S["a1"][pipe][vn]
            for en, el in (("eps0", "0"), ("epsfinal", "final step"), ("epsinit", "initial step")):
                e = v[en]
                rows.append(dict(
                    set=pipe, variant=vn, eps=el, n=e["n_used"],
                    r50=f"{fq(e['r50_quartiles'], 4)} ({e['r50_ge_max']} ≥ 0.1)",
                    r_any=f"{e['r_any_quartiles'][1]:.4f}", ext=f"{e['ext_med_quartiles'][1]:.4f}",
                    no_member=e["no_plateau_member"],
                    aniso=f"{e['aniso_quartiles'][1]:.2f} (n={e['n_aniso']})",
                ))  # fmt: skip
    t["a1_plateau"] = markdown_table(
        rows,
        ["set", "variant", "eps", "n", "r50", "r_any", "ext", "no_member", "aniso"],
        labels=None,
    )
    rows = []
    for pipe in ("H3", "H0"):
        for vn in VARIANT_NAMES:
            v = S["a1"][pipe][vn]
            rows.append(dict(
                set=pipe, variant=vn, eps_final=f"{v['eps_final_quartiles'][1]:.4f}",
                eps_init=f"{v['eps_init_quartiles'][1]:.3f}",
                clamped=fmt_frac(v["final_step_below_grid"]),
                below=fmt_frac(v["min_below_truth"]),
                drop=f"{v['min_drop_among_lower_quartiles'][1]:.5f} (n={v['n_lower']})",
                tied=f"{v['n_tied_at_min_quartiles'][1]:.0f}",
                rad=fq(v["min_radius_smallest_tied_quartiles"], 4),
                cen=fq(v["min_set_centroid_offset_quartiles"], 4),
            ))  # fmt: skip
    t["a1_minimum"] = markdown_table(
        rows,
        [
            "set",
            "variant",
            "eps_final",
            "eps_init",
            "clamped",
            "below",
            "drop",
            "tied",
            "rad",
            "cen",
        ],
    )
    rows = []
    for pipe in ("H3", "H0"):
        for vn in VARIANT_NAMES:
            m = S["a1"][pipe][vn]["monotone_rays"]
            rows.append(dict(
                set=pipe, variant=vn, rays=f"{100 * m['all_radii']:.1f}%",
                rays5=f"{100 * m['from_0p005']:.1f}%",
                cases=fmt_frac(m["cases_all_rays_all_radii"]),
                cases5=fmt_frac(m["cases_all_rays_from_0p005"]),
            ))  # fmt: skip
    t["a1_monotone"] = markdown_table(rows, ["set", "variant", "rays", "rays5", "cases", "cases5"])
    # change and plateau fraction by radius (H3)
    rows = []
    for pipe in ("H3", "H0"):
        for vn in VARIANT_NAMES:
            v = S["a1"][pipe][vn]
            row: Dict[str, Any] = dict(set=pipe, variant=vn)
            for k, r in enumerate(S["a1"][pipe]["radii"]):
                row[f"{r:g}"] = (
                    f"{v['abs_change_med_by_radius'][k]:.3f} / "
                    f"{v['frac_le_truth_by_radius'][k]:.2f}"
                )
            rows.append(row)
    cols = ["set", "variant"] + [f"{r:g}" for r in S["a1"]["H3"]["radii"]]
    t["a1_by_radius"] = markdown_table(rows, cols)
    # A2
    a2 = S["a2"]
    rows = []
    for tag, lab in (
        ("net_filter", "|sin eta| >= 0.3 (net set)"),
        ("all_peaks", "all recorded peaks"),
    ):
        x = a2[tag]
        rows.append(dict(
            peaks=lab, n_peaks=f"{x['n_peaks_quartiles'][1]:.0f}",
            both=fq(x["rms3_quartiles"], 4),
            frame=f"{x['frame_only_quartiles'][1]:.4f}",
            pixel=f"{x['pixel_only_quartiles'][1]:.4f}",
            aniso=f"{x['aniso_quartiles'][1]:.2f}",
        ))  # fmt: skip
    t["a2_quant_scale"] = markdown_table(
        rows, ["peaks", "n_peaks", "both", "frame", "pixel", "aniso"]
    )
    rows = []
    for pipe in ("H3", "H0"):
        rows.append(dict(
            set=pipe, ratio=fq(a2[pipe + "_err_over_bound_quartiles"], 2),
            below=fmt_frac(a2[pipe + "_err_below_bound"]),
        ))  # fmt: skip
    t["a2_error_vs_scale"] = markdown_table(rows, ["set", "ratio", "below"])
    # A3
    rows = []
    for pipe in ("H3", "H0"):
        for vn in VARIANT_NAMES:
            v = S["a1"][pipe][vn]
            lt = v["lower_than_result"]
            rows.append(dict(
                set=pipe, variant=vn,
                rho=fq(v["spearman_quartiles"], 2),
                rho_fine=f"{v['spearman_le0p03_quartiles'][1]:.2f}",
                neg=f"{v['spearman_lt0_count']}/{S['a1'][pipe]['n']}",
                lower=f"{lt['cases_with_lower']}/{lt['n_usable']}",
                far=f"{lt['points_farther']}/{lt['points']}",
                cases_far=f"{lt['cases_with_farther']}/{lt['cases_with_lower']}",
                case_frac=f"{lt['per_case_fraction_farther_quartiles'][1]:.3f}",
            ))  # fmt: skip
    t["a3_correlation"] = markdown_table(
        rows, ["set", "variant", "rho", "rho_fine", "neg", "lower", "far", "cases_far", "case_frac"]
    )
    rows = []
    for pipe in ("H3", "H0"):
        a = S["a3_t5"][pipe]
        for nm in ("vm_smallbox", "mc_smallstep"):
            x = a[nm]
            rows.append(dict(
                set=pipe, cont=nm, before=f"{a['ang_res_quartiles'][1]:.4f}",
                after=f"{x['ang_after_quartiles'][1]:.4f}", closed=f"{x['closed_median']:.2f}",
                reduced=fmt_frac(x["error_reduced"]), p=f"{x['wilcoxon_p']:.2g}",
                rho=f"{x['spearman_closed_vs_error_drop']:.2f}",
                voxels=(
                    f"{x['voxel_clustered']['voxels_improved']}/{x['voxel_clustered']['n_voxels']}"
                    f" (p={x['voxel_clustered']['p_sign']:.2g})"
                ),
            ))  # fmt: skip
    t["a3_gap_closure"] = markdown_table(
        rows, ["set", "cont", "before", "after", "closed", "reduced", "p", "rho", "voxels"]
    )
    return t


def main() -> None:
    a1 = load_a1()
    S: Dict[str, Any] = dict(a1={}, a3_t5={})
    for pipe, rec in a1.items():
        S["a1"][pipe] = summarise_a1(rec)
        S["a3_t5"][pipe] = summarise_a3_t5(rec)
    S["a2"] = summarise_a2(a1)
    S["gn_floor_deg"] = GN_FLOOR_DEG
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    (OUT_DIR / "summary.json").write_text(json.dumps(S, indent=1))
    write_tables(OUT_DIR / "tables.md", tables(S))
    print((OUT_DIR / "tables.md").read_text())


if __name__ == "__main__":
    main()
