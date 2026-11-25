"""tex1: the hard gate before the freeze, on the ``.mic`` of ``tex24`` read back.

  uv run python scripts/autoresearch/tex1/tex1_sample_check.py

Checks (the tolerance ``FILE_ROUNDING_DEG`` is the ``.mic`` Euler rounding: phi1 is written to
1e-4 degree):
  1. sharp voxels (grains of no class) equal ``R_g``;
  2. every noise voxel is within 1 degree of the smooth field (``R_g`` or the gradient field);
  3. edge pairs of gradient-only grains are at most 0.5 degree apart;
  4. edge and touching pairs of noise-only grains are at most 2 degrees apart; edge pairs of grains
     with both at most 2.5 (2 of noise plus 0.5 of gradient: the one place the owner's 2 is
     exceeded; their touching pairs are only reported, a corner pair has about twice the step);
  5. the sample exists as planned: 250 grains, the voxel count, class counts agree with the file.
Reports per class: grains, voxels, edge-pair and touching-pair max, median and p99; per gradient
grain its strength ``t_g`` and its spread; per grain boundary the least misorientation across it
(counts under 5, 2 and 1 degree). All go to ``tex1/sample/sample_check.json``
and ``sample_tables.md``. Exit 1 if the gate fails.
"""

import json
import sys
from pathlib import Path
from typing import Any, Dict, List

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1] / "common"))
import doc_tables as DT  # noqa: E402
import tex1_lib as F  # noqa: E402
import tex1_paths as P  # noqa: E402

CHECK_JSON = P.BENCH / "sample" / "sample_check.json"
TABLES = P.BENCH / "sample" / "sample_tables.md"
TOL = F.FILE_ROUNDING_DEG


def boundary_minima(
    ids: np.ndarray, R: np.ndarray, pairs: np.ndarray, sym: np.ndarray
) -> Dict[Any, float]:
    """Least cubic misorientation (degrees) across every grain boundary, over the touching voxel
    pairs of different grains: {(grain a, grain b): minimum}."""
    from icenine.ray_indexing import misorientation_deg

    cross = ids[pairs[:, 0]] != ids[pairs[:, 1]]
    pr = pairs[cross]
    mis = np.asarray(misorientation_deg(R[pr[:, 0]], R[pr[:, 1]], sym), dtype=float)
    ga = np.minimum(ids[pr[:, 0]], ids[pr[:, 1]])
    gb = np.maximum(ids[pr[:, 0]], ids[pr[:, 1]])
    best: Dict[Any, float] = {}
    for a, b, m in zip(ga.tolist(), gb.tolist(), mis.tolist()):
        k = (a, b)
        if k not in best or m < best[k]:
            best[k] = m
    return best


def check(
    R: np.ndarray,
    R_smooth: np.ndarray,
    ids: np.ndarray,
    cls_of_grain: Dict[int, int],
    R_g: Dict[int, np.ndarray],
    edge: np.ndarray,
    touch: np.ndarray,
    sym: np.ndarray,
) -> Dict[str, Any]:
    """The numbers and the gate, from arrays (so the tests can run it on small samples)."""
    cls = np.array([cls_of_grain[int(g)] for g in ids])
    Rg_v = np.stack([R_g[int(g)] for g in ids])
    out: Dict[str, Any] = dict(tolerance_deg=TOL)
    ang_g = F.relative_angle_deg(R, Rg_v)
    sharp = cls == F.SHARP
    out["sharp_max_deg_from_R_g"] = float(ang_g[sharp].max()) if sharp.any() else 0.0
    noisy = (cls == F.NOISE_ONLY) | (cls == F.BOTH)
    ang_s = F.relative_angle_deg(R, R_smooth)
    out["noise_max_deg_from_smooth_field"] = float(ang_s[noisy].max()) if noisy.any() else 0.0
    out["smooth_vs_R_g_max_deg_noise_only_grains"] = float(
        ang_g[cls == F.NOISE_ONLY].max() if (cls == F.NOISE_ONLY).any() else 0.0
    )
    pc = {}
    for name, pairs in (("edge", edge), ("touching", touch)):
        same = ids[pairs[:, 0]] == ids[pairs[:, 1]]
        pr = pairs[same]
        d = F.relative_angle_deg(R[pr[:, 0]], R[pr[:, 1]])
        pcls = cls[pr[:, 0]]
        pc[name] = {F.CLASS_NAMES[c]: F.summary(d[pcls == c]) for c in range(4)}
        pc[name]["all"] = F.summary(d)
    out["pairs"] = pc
    # per class grains and voxels
    out["classes"] = {
        F.CLASS_NAMES[c]: dict(
            n_grains=int(sum(1 for g, k in cls_of_grain.items() if k == c)),
            n_voxels=int((cls == c).sum()),
        )
        for c in range(4)
    }
    best = boundary_minima(ids, R, touch, sym)
    v = np.array(list(best.values()))
    out["boundaries"] = dict(
        n=int(len(v)),
        under_5=int((v < 5.0).sum()),
        under_2=int((v < 2.0).sum()),
        under_1=int((v < 1.0).sum()),
        least_deg=float(v.min()) if len(v) else None,
        median_least_deg=float(np.median(v)) if len(v) else None,
        under_5_share=float((v < 5.0).mean()) if len(v) else None,
    )
    out["boundary_minima"] = {f"{a}-{b}": m for (a, b), m in sorted(best.items())}
    e_go = pc["edge"]["gradient"]["max"]
    e_no = pc["edge"]["noise"]["max"]
    e_bo = pc["edge"]["both"]["max"]
    t_no = pc["touching"]["noise"]["max"]
    out["gate"] = F.gate(
        {
            "sharp_voxels_equal_R_g": out["sharp_max_deg_from_R_g"] <= TOL,
            "noise_within_1_degree_of_smooth_field": (
                out["noise_max_deg_from_smooth_field"] <= F.NOISE_BALL_DEG + TOL
            ),
            "gradient_only_edge_pairs_at_most_0.5": e_go is None
            or e_go <= F.GRADIENT_CAP_DEG + TOL,
            "noise_only_pairs_at_most_2": all(
                x is None or x <= 2 * F.NOISE_BALL_DEG + TOL for x in (e_no, t_no)
            ),
            "both_edge_pairs_at_most_2.5": e_bo is None
            or e_bo <= 2 * F.NOISE_BALL_DEG + F.GRADIENT_CAP_DEG + TOL,
            "has_pairs": len(edge) > 0 and len(touch) > 0,
        }
    )
    return out


def tables(res: Dict[str, Any], grads: List[Dict[str, Any]]) -> Dict[str, str]:
    rows = []
    for c in F.CLASS_NAMES:
        e, t = res["pairs"]["edge"][c], res["pairs"]["touching"][c]

        def f(x: Any) -> str:
            return "-" if x is None else f"{x:.4f}"

        rows.append(
            dict(
                cls=c,
                grains=res["classes"][c]["n_grains"],
                voxels=res["classes"][c]["n_voxels"],
                edge_max=f(e["max"]),
                edge_median=f(e["median"]),
                edge_p99=f(e["p99"]),
                touch_max=f(t["max"]),
                touch_median=f(t["median"]),
                touch_p99=f(t["p99"]),
            )
        )
    t1 = DT.markdown_table(
        rows,
        [
            "cls",
            "grains",
            "voxels",
            "edge_max",
            "edge_median",
            "edge_p99",
            "touch_max",
            "touch_median",
            "touch_p99",
        ],
        formats={"grains": ",d", "voxels": ",d"},
    )
    b = res["boundaries"]
    t2 = DT.markdown_table(
        [
            dict(
                boundaries=b["n"],
                under_5=b["under_5"],
                under_2=b["under_2"],
                under_1=b["under_1"],
                least=b["least_deg"],
                median_least=b["median_least_deg"],
            )
        ],
        ["boundaries", "under_5", "under_2", "under_1", "least", "median_least"],
        formats={"least": ".3f", "median_least": ".2f"},
    )
    sp = np.array([g["measured_spread_deg"] for g in grads])
    t3 = DT.markdown_table(
        [
            dict(
                gradient_grains=len(grads),
                spread_median=float(np.median(sp)) if len(sp) else 0.0,
                spread_p90=float(np.percentile(sp, 90)) if len(sp) else 0.0,
                spread_max=float(sp.max()) if len(sp) else 0.0,
            )
        ],
        ["gradient_grains", "spread_median", "spread_p90", "spread_max"],
        formats={k: ".2f" for k in ("spread_median", "spread_p90", "spread_max")},
    )
    return dict(tex1_sample_pairs=t1, tex1_sample_boundaries=t2, tex1_sample_spreads=t3)


def main() -> int:
    from icenine.mic_file import MicFile
    from icenine.ray_indexing import proper_rotations
    from icenine.symmetry import create_cubic_symmetry

    sys.path.insert(0, str(P.ICE / "scripts" / "autoresearch" / "ar1"))
    import misses as M  # noqa: E402

    info = json.loads(P.SAMPLE_JSON.read_text())
    assert P.sha256(P.TRUTH) == info["truth_sha256"], "the truth .mic differs from sample.json"
    back = MicFile.read(str(P.TRUTH))
    R = np.array([v.orientation for v in back.voxels], dtype=np.float64)
    ids = np.load(str(P.GRAINMAP))["grain"]
    assert len(R) == len(ids) == P.N_FINE
    fv, side = M.voxel_vertices_xy(P.GRID)
    fc = fv.mean(axis=1)
    s = float(side.max())
    edge, touch = F.edge_pairs(fc, s), F.touching_pairs(fc, s)
    cls_of = {int(g["id"]): int(g["cls"]) for g in info["grains"]}
    R_g = {int(g["id"]): np.array(g["R_centre"]) for g in info["grains"]}
    noise_g = [g for g, c in cls_of.items() if c in (F.NOISE_ONLY, F.BOTH)]
    grad_g = [g for g, c in cls_of.items() if c in (F.GRADIENT_ONLY, F.BOTH)]
    params = {
        int(r["id"]): dict(M1=r["M1"], M2=r["M2"], t_g_deg=r["t_g_deg"], scale=r["scale"])
        for r in info["gradient_grains"]
    }
    # the smooth field, rebuilt from the recorded field parameters
    R2, R_smooth, _ = F.textured_orientations(
        fc, ids, R_g, noise_g, grad_g, edge, P.SEED, field_params=params
    )
    sym = proper_rotations(create_cubic_symmetry(3.6))
    res = check(R, R_smooth, ids, cls_of, R_g, edge, touch, sym)
    res["rebuilt_truth_vs_file_max_deg"] = float(F.relative_angle_deg(R, R2).max())
    res["n_grains"] = len(cls_of)
    res["n_voxels"] = int(len(R))
    res["gate"]["checks"]["file_equals_rebuild_to_rounding"] = bool(
        res["rebuilt_truth_vs_file_max_deg"] <= TOL
    )
    res["gate"]["checks"]["n_grains_250"] = res["n_grains"] == P.N_CELLS
    res["gate"]["passed"] = bool(all(res["gate"]["checks"].values()))
    # per gradient grain: strength and the spread measured on the file
    ang_g = F.relative_angle_deg(R_smooth, np.stack([R_g[int(g)] for g in ids]))
    grads = []
    for r in info["gradient_grains"]:
        g = int(r["id"])
        sel = ids == g
        grads.append(
            dict(
                id=g,
                cls=F.CLASS_NAMES[cls_of[g]],
                n_voxels=int(sel.sum()),
                t_g_deg=r["t_g_deg"],
                measured_spread_deg=float(ang_g[sel].max()),
                edge_max_smooth_deg=r["edge_max_deg"],
            )
        )
    res["gradient_grains"] = grads
    CHECK_JSON.write_text(json.dumps(res, indent=1) + "\n")
    DT.write_tables(TABLES, tables(res, grads))
    print(json.dumps(res["gate"], indent=1))
    print(TABLES.read_text())
    return 0 if res["gate"]["passed"] else 1


if __name__ == "__main__":
    sys.exit(main())
