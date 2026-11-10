"""Task D steps 1 and 2: summaries of the provenance (truth) and of the geometric pairing.

Reads the provenance written by scripts/phase_d/render_provenance.py and the pairs written by
pairing_run.py; writes benchmarks/eval_reduction/D/step1/{summary.json,tables.md} (what the
spots are) and .../step2/{summary.json,tables.md} (pairing precision / recall, origin error,
ambiguity). Figures are made by pairing_figures.py.

Truth: two components (one per detector) of a frame are a true pair when they share a
(voxel, reflection) source. Classes of a component: single-source (one voxel and reflection),
one grain with several voxels, several grains. A pair belongs to the larger class of its two
components. Rates carry Wilson intervals (scripts/common/stats.py); no significance tests.

  uv run python scripts/eval_reduction/pairing_summary.py
"""

import json
import sys
from pathlib import Path
from typing import Any, Dict, List

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[0] / "common"))
import doc_tables as DT  # noqa: E402
import pairing_common as pc  # noqa: E402
import pairing_run as PR  # noqa: E402
import stats as ST  # noqa: E402

from icenine import spot_pairing as sp  # noqa: E402

MAIN = "ring_tight"
OUT1 = pc.BENCH / "step1"
OUT2 = pc.BENCH / "step2"
CLS = pc.CLASS_NAMES


def pct(x: np.ndarray, q: List[float]) -> List[float]:
    return [float(v) for v in np.percentile(x, q)] if len(x) else [float("nan")] * len(q)


def rate(k: int, n: int) -> Dict[str, Any]:
    lo, hi = ST.wilson(int(k), int(n)) if n else (float("nan"), float("nan"))
    return {
        "k": int(k),
        "n": int(n),
        "rate": float(k / n) if n else float("nan"),
        "lo": lo,
        "hi": hi,
    }


def fmt_rate(r: Dict[str, Any]) -> str:
    if r["n"] == 0:
        return "n/a (0)"
    return f"{r['rate']:.3f} [{r['lo']:.3f}, {r['hi']:.3f}] ({r['k']}/{r['n']})"


def step1(prov: Dict, ct: Dict, cls: np.ndarray, a: np.ndarray, b: np.ndarray, npk: np.ndarray):
    N = len(cls)
    det = ct["det"]
    linked = np.zeros(N, dtype=bool)
    linked[a] = True
    linked[b] = True
    rows = []
    out: Dict[str, Any] = {"n_components": int(N), "n_triangles": int(len(prov["tri_comp"]))}
    for d in (0, 1):
        for k in range(3):
            m = (det == d) & (cls == k)
            rows.append(
                {
                    "detector": d,
                    "class": CLS[k],
                    "components": int(m.sum()),
                    "with a true partner": int((m & linked).sum()),
                    "median pixels": float(np.median(ct["npix"][m])) if m.any() else None,
                    "p90 pixels": float(np.percentile(ct["npix"][m], 90)) if m.any() else None,
                    "median sources": float(np.median(ct["n_src"][m])) if m.any() else None,
                    "max sources": int(ct["n_src"][m].max()) if m.any() else None,
                }
            )
            out[f"det{d}_class{k}"] = rows[-1]
    lc = np.maximum(cls[a], cls[b])
    link_rows = []
    for k in range(3):
        m = lc == k
        link_rows.append(
            {
                "class of pair": CLS[k],
                "true links": int(m.sum()),
                "peaks (voxel, reflection)": int(npk[m].sum()),
            }
        )
    out["links"] = link_rows
    out["n_links"] = int(len(a))
    out["components_per_detector"] = {f"det{d}": int((det == d).sum()) for d in (0, 1)}
    out["several_grain_components_both_detectors"] = int((cls == 2).sum())
    out["single_source_components_both_detectors"] = int((cls == 0).sum())
    out["n_peaks_both_detectors"] = int(npk.sum())
    out["unpaired_components"] = {f"det{d}": int(((det == d) & ~linked).sum()) for d in (0, 1)}
    out["unpaired_single_source"] = int(((cls == 0) & ~linked).sum())
    out["n_single_source_total"] = int((cls == 0).sum())
    out["edge_touching_components"] = {f"det{d}": int(ct["edge"][det == d].sum()) for d in (0, 1)}
    # partners per component
    pa = np.bincount(a, minlength=N)
    pb = np.bincount(b, minlength=N)
    out["partners_per_linked_det0"] = {
        "max": int(pa.max()),
        "share_more_than_one": float((pa[det == 0] > 1).sum() / max((pa[det == 0] > 0).sum(), 1)),
    }
    out["partners_per_linked_det1"] = {
        "max": int(pb.max()),
        "share_more_than_one": float((pb[det == 1] > 1).sum() / max((pb[det == 1] > 0).sum(), 1)),
    }
    tbl1 = DT.markdown_table(
        rows,
        [
            "detector",
            "class",
            "components",
            "with a true partner",
            "median pixels",
            "p90 pixels",
            "median sources",
            "max sources",
        ],
        {"median pixels": ".0f", "p90 pixels": ".0f", "median sources": ".0f"},
    )
    tbl2 = DT.markdown_table(
        link_rows, ["class of pair", "true links", "peaks (voxel, reflection)"]
    )
    return out, {"step1_components": tbl1, "step1_links": tbl2}


def grain_spot_counts(
    prov: Dict, ct: Dict, cls: np.ndarray, geo: Dict[str, Any], a: np.ndarray
) -> Dict[str, Any]:
    """Per grain: one-grain detector-0 spots that have a true partner (what an indexing fit over a
    grain's spots could use), and the grain size in voxels."""
    comp = prov["tri_comp"].astype(np.int64)
    first = np.unique(comp, return_index=True)[1]
    grain_of = np.full(len(cls), -1, dtype=np.int64)
    grain_of[comp[first]] = geo["voxel_grain"][prov["tri_voxel"][first].astype(np.int64)]
    linked = np.zeros(len(cls), dtype=bool)
    linked[a] = True
    sel = (ct["det"] == 0) & (cls == 1) & linked
    n_grains = int(geo["voxel_grain"].max()) + 1
    per_grain = np.bincount(grain_of[sel], minlength=n_grains)
    size = np.bincount(geo["voxel_grain"], minlength=n_grains)
    return {
        "n_grains": n_grains,
        "one_grain_det0_spots_with_partner_per_grain_median_p10_min": [
            float(np.median(per_grain)),
            float(np.percentile(per_grain, 10)),
            int(per_grain.min()),
        ],
        "grain_size_voxels_median_p10_min": [
            float(np.median(size)),
            float(np.percentile(size, 10)),
            int(size.min()),
        ],
    }


def partner_size_ratio(ct: Dict, cls: np.ndarray, a: np.ndarray, b: np.ndarray) -> Dict[str, Any]:
    """Size of a detector-1 spot against its true detector-0 partner (one-grain pairs, neither
    spot touching a detector border). A grain's peak is a parallel projection, so the two spots
    should have about the same size; a scaling by L1/L0 would give an area ratio of about 2.6."""
    m = (cls[a] == 1) & (cls[b] == 1) & ~ct["edge"][a] & ~ct["edge"][b]
    i0, i1 = a[m], b[m]
    out: Dict[str, Any] = {"n_pairs": int(m.sum())}
    for key, x0, x1 in (
        ("npix", ct["npix"][i0], ct["npix"][i1]),
        ("width", ct["cmax"][i0] - ct["cmin"][i0] + 1, ct["cmax"][i1] - ct["cmin"][i1] + 1),
        ("height", ct["rmax"][i0] - ct["rmin"][i0] + 1, ct["rmax"][i1] - ct["rmin"][i1] + 1),
    ):
        out[f"{key}_ratio_det1_over_det0_pct10_50_90"] = pct(
            x1.astype(float) / x0.astype(float), [10, 50, 90]
        )
    return out


def evaluate(
    name: str,
    suffix: str,
    ct: Dict,
    cls: np.ndarray,
    a: np.ndarray,
    b: np.ndarray,
    npk: np.ndarray,
) -> Dict[str, Any]:
    N = len(cls)
    r = np.load(pc.CACHE / f"pairs_full_{name}{suffix}.npz")
    ch = r["chosen"].astype(bool)
    c0, c1 = r["c0"][ch], r["c1"][ch]
    tk = pc.pair_keys(a, b, N)
    pk = pc.pair_keys(c0, c1, N)
    correct = np.isin(pk, tk)
    found = np.isin(tk, pk)
    pcl = np.maximum(cls[c0], cls[c1])
    lcl = np.maximum(cls[a], cls[b])
    linked = np.zeros(N, dtype=bool)
    linked[a] = True
    linked[b] = True
    res: Dict[str, Any] = {
        "predicted": int(ch.sum()),
        "candidates_before_assignment": int(len(r["c0"])),
        "precision": rate(correct.sum(), len(correct)),
        "link_recall": rate(found.sum(), len(found)),
        "peak_recall": rate(int(npk[found].sum()), int(npk.sum())),
        "truth_in_candidates": rate(
            int(np.isin(tk, pc.pair_keys(r["c0"], r["c1"], N)).sum()), len(tk)
        ),
    }
    for k in range(3):
        res[f"precision_class{k}"] = rate(correct[pcl == k].sum(), (pcl == k).sum())
        res[f"link_recall_class{k}"] = rate(found[lcl == k].sum(), (lcl == k).sum())
    # what the false pairs are
    wrong = ~correct
    both_orphan = wrong & ~linked[c0] & ~linked[c1]
    one_orphan = wrong & (linked[c0] ^ linked[c1])
    wrong_partner = wrong & linked[c0] & linked[c1]
    res["false_pairs"] = {
        "total": int(wrong.sum()),
        "both spots have no true partner": int(both_orphan.sum()),
        "one spot has no true partner": int(one_orphan.sum()),
        "both have a true partner, but another": int(wrong_partner.sum()),
    }
    # ambiguity: candidates per detector-0 component that has a true link
    ncand = r["n_cand_c0"]
    have = np.zeros(N, dtype=bool)
    have[a] = True
    d0 = have & (ct["det"] == 0)
    res["candidates_per_linked_det0"] = {
        "median": float(np.median(ncand[d0])),
        "p90": float(np.percentile(ncand[d0], 90)),
        "share_with_more_than_one": float((ncand[d0] > 1).mean()),
        "share_with_zero": float((ncand[d0] == 0).mean()),
    }
    return res


def origin_error(
    name: str, suffix: str, ct: Dict, cls: np.ndarray, a: np.ndarray, b: np.ndarray, geo
):
    """Origin error of the correct, chosen pairs of single-grain components (truth = mean of the
    source voxel centroids of the detector-0 component)."""
    N = len(cls)
    r = np.load(pc.CACHE / f"pairs_full_{name}{suffix}.npz")
    ch = r["chosen"].astype(bool)
    c0, c1 = r["c0"][ch], r["c1"][ch]
    ok = np.isin(pc.pair_keys(c0, c1, N), pc.pair_keys(a, b, N)) & (cls[c0] == 1) & (cls[c1] == 1)
    c0, c1 = c0[ok], c1[ok]
    q_true = ct["truth_xy"][c0]
    om_true = ct["omega_mean"][c0]
    lab_true = sp.sample_to_lab_xy(q_true, om_true)
    o_lab = r["origin_lab"][ch][ok][:, :2]
    o_smp = r["origin_sample"][ch][ok]
    err_lab = np.linalg.norm(o_lab - lab_true, axis=1) * 1000.0
    err_smp = np.linalg.norm(o_smp - q_true, axis=1) * 1000.0
    cond = r["cond"][ch][ok]
    sig = r["sigma_mm"][ch][ok] * 1000.0
    omc = 0.5 * (geo["omega_lo"] + geo["omega_hi"])[ct["frame"][c0]]
    # arc part: the position error from using the frame-centre omega instead of the exact one
    q_at_centre = sp.lab_to_sample_xy(lab_true, omc)
    err_arc = np.linalg.norm(q_at_centre - q_true, axis=1) * 1000.0
    out: Dict[str, Any] = {
        "n_pairs": int(len(c0)),
        "lab_um_pct10_50_90_99": pct(err_lab, [10, 50, 90, 99]),
        "sample_um_pct10_50_90_99": pct(err_smp, [10, 50, 90, 99]),
        "arc_term_um_pct50_90_max": [
            float(np.median(err_arc)),
            float(np.percentile(err_arc, 90)),
            float(err_arc.max()),
        ],
        "share_sample_err_below_um": {
            str(t): float((err_smp <= t).mean()) for t in (5, 10, 20, 50, 100)
        },
        "sigma_um_pct10_50_90": pct(sig, [10, 50, 90]),
        # predicted sigma (tight tolerances of the preset) over the measured sample-frame error
        "sigma_over_sample_error_medians_overall": float(np.median(sig) / np.median(err_smp)),
    }
    bins = [(0.0, 0.02), (0.02, 0.05), (0.05, 0.1), (0.1, 1.01)]
    rows = []
    for lo, hi in bins:
        m = (cond >= lo) & (cond < hi)
        rows.append(
            {
                "ray z-slope (abs u_z)": f"{lo:g} to {min(hi, 1):g}",
                "pairs": int(m.sum()),
                "lab-frame error median (um)": float(np.median(err_lab[m])) if m.any() else None,
                "lab-frame error p90 (um)": (
                    float(np.percentile(err_lab[m], 90)) if m.any() else None
                ),
                "sample-frame error median (um)": float(np.median(err_smp[m])) if m.any() else None,
                "sample-frame error p90 (um)": (
                    float(np.percentile(err_smp[m], 90)) if m.any() else None
                ),
                "predicted sigma median (um)": float(np.median(sig[m])) if m.any() else None,
            }
        )
    out["by_cond"] = rows
    out["sigma_over_sample_error_medians_by_bin"] = [
        (
            r["predicted sigma median (um)"] / r["sample-frame error median (um)"]
            if r["pairs"]
            else None
        )
        for r in rows
    ]
    tbl = DT.markdown_table(
        rows,
        list(rows[0].keys()),
        {k: ".1f" for k in rows[0] if "um" in k},
    )
    arrays = dict(
        err_lab=err_lab,
        err_smp=err_smp,
        cond=cond,
        sigma=sig,
        truth_xy=q_true,
        om_true=om_true,
        c0=c0,
    )
    return out, tbl, arrays


def derived_geometry(geo: Dict[str, Any], name: str) -> Dict[str, Any]:
    """Geometry numbers quoted in the prose, and the per-spot error budget of a scattering
    vector from one matched pair (the "Proposal, not tested" block of the findings doc).

    One pixel moves a ray direction by pixel / gap on one detector (pixel / gap x sqrt(2) for
    both, in quadrature). An azimuthal ray error delta maps to about delta / sin(2 theta) on the
    scattering vector. The frames are omega_step wide, so the frame-centre omega is off by up to
    omega_step / 2; that moves the sample-frame scattering vector by that angle x sin(g, z)."""
    gs = pc.detector_geoms(geo)
    L0, L1 = gs[0].distance_x(), gs[1].distance_x()
    pw = float(geo["det0_pw"])
    side = float(geo["voxel_side"])
    gap = L1 - L0
    one = np.degrees(pw / gap)
    rings = np.degrees(np.asarray(geo["rings_two_theta"]))
    om_step = float(np.degrees(np.mean(geo["omega_hi"] - geo["omega_lo"])))
    r = np.load(pc.CACHE / f"pairs_full_{name}.npz")
    dn = r["dirn"][r["chosen"].astype(bool)]
    g = dn - sp.BEAM_DIR  # scattering vector direction (k_out - k_in), lab frame
    sin_gz = np.linalg.norm(g[:, :2], axis=1) / np.linalg.norm(g, axis=1)
    return {
        "pixel_um": pw * 1000.0,
        "voxel_side_um": side * 1000.0,
        "L0_mm": L0,
        "L1_mm": L1,
        "gap_mm": gap,
        "L1_over_L0": L1 / L0,
        "voxel_term_um": side * (L1 / L0 - 1.0) * 1000.0,
        "ray_dir_err_1px_one_detector_deg": one,
        "ray_dir_err_1px_both_quadrature_deg": one * np.sqrt(2.0),
        "ring_two_theta_deg_min_max": [float(rings.min()), float(rings.max())],
        "n_rings": int(len(rings)),
        "g_err_1px_deg_at_min_max_two_theta": [
            float(one / np.sin(np.radians(rings.min()))),
            float(one / np.sin(np.radians(rings.max()))),
        ],
        "g_err_1px_both_quadrature_deg_at_min_max_two_theta": [
            float(one * np.sqrt(2.0) / np.sin(np.radians(rings.min()))),
            float(one * np.sqrt(2.0) / np.sin(np.radians(rings.max()))),
        ],
        "omega_frame_width_deg": om_step,
        "omega_half_width_deg": om_step / 2.0,
        "median_sin_g_z_chosen_pairs": float(np.median(sin_gz)),
        "g_err_omega_half_width_median_deg": float(om_step / 2.0 * np.median(sin_gz)),
    }


def measured_ray_error(
    geo: Dict[str, Any], name: str, ct: Dict, cls: np.ndarray, a: np.ndarray, b: np.ndarray
) -> Dict[str, Any]:
    """Measured radial error of the ray direction: 2 theta of the line through the two centroids
    minus the nearest ring angle, over the correct one-grain pairs of `name` (the origin-error
    subset). The line uses both centroids, so this includes both detectors' errors. The azimuthal
    error (a turn about the beam) is not measured; the scattering-vector figure assumes it is
    about as large as the radial one."""
    N = len(cls)
    r = np.load(pc.CACHE / f"pairs_full_{name}.npz")
    ch = r["chosen"].astype(bool)
    c0, c1 = r["c0"][ch], r["c1"][ch]
    ok = np.isin(pc.pair_keys(c0, c1, N), pc.pair_keys(a, b, N)) & (cls[c0] == 1) & (cls[c1] == 1)
    dn = r["dirn"][ch][ok]
    tt = np.degrees(np.arccos(np.clip(dn @ sp.BEAM_DIR, -1.0, 1.0)))
    rings = np.degrees(np.asarray(geo["rings_two_theta"]))
    d = np.abs(tt[:, None] - rings).min(axis=1)
    med = float(np.median(d))
    p90 = float(np.percentile(d, 90))
    return {
        "measured_radial_ray_err_subset": (
            f"correct chosen one-grain pairs of {name}; selected by the ring gate, which cuts the"
            " tail of this same error"
        ),
        "measured_radial_ray_err_n_pairs": int(ok.sum()),
        "measured_radial_ray_err_deg_pct50_90_99_max": pct(d, [50, 90, 99]) + [float(d.max())],
        "g_err_measured_median_deg_at_min_max_two_theta": [
            float(med / np.sin(np.radians(rings.min()))),
            float(med / np.sin(np.radians(rings.max()))),
        ],
        "g_err_measured_p90_deg_at_min_max_two_theta": [
            float(p90 / np.sin(np.radians(rings.min()))),
            float(p90 / np.sin(np.radians(rings.max()))),
        ],
    }


def measured_ray_error_all_links(
    geo: Dict[str, Any], ct: Dict, cls: np.ndarray, a: np.ndarray, b: np.ndarray
) -> Dict[str, Any]:
    """The same radial ray-direction error as :func:`measured_ray_error`, but over ALL true links
    whose two components are both one-grain, with no gate and no assignment: the line through the
    two true partner centroids (pixel-centre shift +0.5, as in ``pair_frame``). Also gives the
    ring-gate tolerance of the tight presets at these links, (1 px x pw0 + 1 px x pw1) /
    separation (``spot_pairing.pair_frame``), and how many links that gate would cut."""
    gs = pc.detector_geoms(geo)
    g0, g1 = gs[0], gs[1]
    m = (cls[a] == 1) & (cls[b] == 1)
    i0, i1 = a[m], b[m]
    c0, r0 = ct["i_col"][i0] + 0.5, ct["i_row"][i0] + 0.5
    c1, r1 = ct["i_col"][i1] + 0.5, ct["i_row"][i1] + 0.5
    lo = sp.line_origin(g0, g1, c0, r0, c1, r1, 1.0, 1.0)
    dn = lo["dirn"]
    tt = np.degrees(np.arccos(np.clip(dn @ sp.BEAM_DIR, -1.0, 1.0)))
    rings = np.degrees(np.asarray(geo["rings_two_theta"]))
    d = np.abs(tt[:, None] - rings).min(axis=1)
    sep = np.linalg.norm(sp.pixel_to_lab(g1, c1, r1) - sp.pixel_to_lab(g0, c0, r0), axis=-1)
    tol = np.degrees((g0.pw + g1.pw) / sep)
    one = float(np.degrees(float(geo["det0_pw"]) / (g1.distance_x() - g0.distance_x())))
    lo_s, hi_s = np.sin(np.radians(rings.min())), np.sin(np.radians(rings.max()))
    q = pct(d, [50, 90, 99])
    return {
        "all_links_radial_ray_err_subset": (
            "all true links with both components one-grain; no gate, no assignment"
        ),
        "all_links_radial_ray_err_n_links": int(m.sum()),
        "all_links_radial_ray_err_deg_pct50_90_99_max": q + [float(d.max())],
        "all_links_radial_ray_err_share_above_1px_one_detector": float((d > one).mean()),
        "all_links_radial_ray_err_n_above_1px_one_detector": int((d > one).sum()),
        "all_links_g_err_median_deg_at_min_max_two_theta": [q[0] / lo_s, q[0] / hi_s],
        "all_links_g_err_p90_deg_at_min_max_two_theta": [q[1] / lo_s, q[1] / hi_s],
        "all_links_g_err_p99_deg_at_min_max_two_theta": [q[2] / lo_s, q[2] / hi_s],
        "ring_gate_tol_deg_tight_pct50_min_max": [
            float(np.median(tol)),
            float(tol.min()),
            float(tol.max()),
        ],
        "all_links_n_cut_by_ring_gate_tight": int((d > tol).sum()),
    }


def true_link_rejections(geo, ct, cls, a, b):
    """First failing gate (ray, hull, sigma, ring, size) of the true links under MAIN."""
    gs = pc.detector_geoms(geo)
    spec = PR.PRESETS[MAIN]
    params = PR.make_params(geo, spec)
    hull = sp.SampleHull(geo["voxel_tri"].reshape(-1, 2))
    om_c = 0.5 * (geo["omega_lo"] + geo["omega_hi"])
    frame, det = ct["frame"], ct["det"]
    names = ["gate_ray", "gate_hull", "gate_sigma", "gate_ring", "gate_size"]
    fails = {n: 0 for n in names}
    fails["passes all gates"] = 0
    by_class: Dict[str, Dict[str, int]] = {}
    lcl = np.maximum(cls[a], cls[b])
    ck, rk = "i_col", "i_row"
    for f in range(len(om_c)):
        s0 = np.nonzero((frame == f) & (det == 0))[0]
        s1 = np.nonzero((frame == f) & (det == 1))[0]
        loc = {}
        comps = []
        for sel in (s0, s1):
            comps.append(
                {
                    "c_col": ct[ck][sel],
                    "c_row": ct[rk][sel],
                    "npix": ct["npix"][sel].astype(float),
                    "edge": ct["edge"][sel],
                    "w": (ct["cmax"][sel] - ct["cmin"][sel]).astype(float),
                    "h": (ct["rmax"][sel] - ct["rmin"][sel]).astype(float),
                }
            )
        r = sp.pair_frame(gs, hull, float(om_c[f]), comps[0], comps[1], params, return_gates=True)
        pos0 = {int(g): i for i, g in enumerate(s0)}
        pos1 = {int(g): i for i, g in enumerate(s1)}
        m = frame[a] == f
        for ia, ib, k in zip(a[m], b[m], lcl[m]):
            i0, i1 = pos0[int(ia)], pos1[int(ib)]
            why = "passes all gates"
            for n in names:
                if not r[n][i0, i1]:
                    why = n
                    break
            fails[why] += 1
            by_class.setdefault(CLS[k], {}).setdefault(why, 0)
            by_class[CLS[k]][why] += 1
    return fails, by_class


def main() -> None:
    OUT1.mkdir(parents=True, exist_ok=True)
    OUT2.mkdir(parents=True, exist_ok=True)
    geo = pc.load_geometry()
    prov = pc.load_provenance()
    ct = pc.component_table(prov, geo)
    cls = pc.component_class(ct)
    a, b, npk = pc.truth_links(prov)

    s1, t1 = step1(prov, ct, cls, a, b, npk)
    s1["per_grain"] = grain_spot_counts(prov, ct, cls, geo, a)
    s1["partner_size_ratio"] = partner_size_ratio(ct, cls, a, b)
    (OUT1 / "summary.json").write_text(json.dumps(s1, indent=2) + "\n")
    DT.write_tables(OUT1 / "tables.md", t1)

    summary: Dict[str, Any] = {"main_preset": MAIN, "presets": {}}
    rows = []
    names = list(PR.PRESETS)
    todo = [(n, "") for n in names] + [("ring_tight", "_binary"), ("shape_tight", "_binary")]
    for n, sfx in todo:
        if not (pc.CACHE / f"pairs_full_{n}{sfx}.npz").exists():
            continue
        e = evaluate(n, sfx, ct, cls, a, b, npk)
        key = n + sfx
        summary["presets"][key] = e
        rows.append(
            {
                "preset": key,
                "predicted pairs": e["predicted"],
                "precision": fmt_rate(e["precision"]),
                "link recall": fmt_rate(e["link_recall"]),
                "peak recall": f"{e['peak_recall']['rate']:.3f}",
                "recall, one grain": fmt_rate(e["link_recall_class1"]),
                "recall, several grains": fmt_rate(e["link_recall_class2"]),
            }
        )
    tbl_cols = list(rows[0].keys())
    t2: Dict[str, str] = {"step2_presets": DT.markdown_table(rows, tbl_cols)}
    # main preset: false pairs, ambiguity, origin error, rejections
    em = summary["presets"][MAIN]
    fp = em["false_pairs"]
    t2["step2_false_pairs"] = DT.markdown_table(
        [{"cause": k, "pairs": v} for k, v in fp.items()], ["cause", "pairs"]
    )
    oe, tbl_oe, arrays = origin_error(MAIN, "", ct, cls, a, b, geo)
    summary["origin_error"] = oe
    t2["step2_origin_error"] = tbl_oe
    np.savez_compressed(pc.CACHE / "origin_error_arrays.npz", **arrays)
    summary["derived_geometry"] = derived_geometry(geo, MAIN)
    summary["derived_geometry"].update(measured_ray_error(geo, MAIN, ct, cls, a, b))
    summary["derived_geometry"].update(measured_ray_error_all_links(geo, ct, cls, a, b))
    fails, by_class = true_link_rejections(geo, ct, cls, a, b)
    summary["true_link_first_failing_gate"] = fails
    summary["true_link_first_failing_gate_by_class"] = by_class
    t2["step2_rejections"] = DT.markdown_table(
        [{"first failing gate": k, "true links": v} for k, v in fails.items()],
        ["first failing gate", "true links"],
    )
    # centroid variant
    # decision rule (plan): single-source recall below 0.95 runs Step 4
    r_single = em["link_recall_class0"]
    r_one = em["link_recall_class1"]
    summary["decision_step4"] = {
        "rule": "run Step 4 only if recall on single-source components is below 0.95",
        "single_source_links": r_single["n"],
        "single_source_recall": r_single["rate"] if r_single["n"] else None,
        "one_grain_recall": r_one,
        "run_step4": bool(
            (r_single["n"] > 0 and r_single["rate"] < 0.95)
            or (r_single["n"] == 0 and r_one["rate"] < 0.95)
        ),
    }
    (OUT2 / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    DT.write_tables(OUT2 / "tables.md", t2)
    print(json.dumps(summary["decision_step4"], indent=2))
    print(t2["step2_presets"])


if __name__ == "__main__":
    main()
