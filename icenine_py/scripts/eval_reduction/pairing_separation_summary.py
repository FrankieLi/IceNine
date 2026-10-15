"""Task D step 3 (b): summarize pairing_attribution.py separate (truth vs wrong orientations).

Scores per orientation (see pairing_attribution.py): `pixel` (1 - cost, strict C++ rule, radius 0)
and the completeness of the attributed peak list (A: pairs with origin within 15.7 um, about the
90th percentile of the origin error, B: ray through both spot boxes, A30: 30 um, oracle: the
voxel's own true pairs), ray-direction tolerance 0.05 degrees (``TOL``). Everything is normalised
by the voxel's own score at the true orientation, so the columns are comparable. Counts are over
the 30 voxels; no tests.

  uv run python scripts/eval_reduction/pairing_separation_summary.py
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

OUT = pc.BENCH / "step3"
TOL = 0.05  # default direction tolerance (degrees), frame window 0
FW = 0
COSTS = [
    ("pixel", "pixel overlap (C++ rule)"),
    (f"A_c_{TOL}_{FW}", "attributed, A r=15.7 um"),
    (f"A30_c_{TOL}_{FW}", "attributed, A r=30 um"),
    (f"B_c_{TOL}_{FW}", "attributed, B footprint"),
    (f"oracle_c_{TOL}_{FW}", "attributed, oracle (true pairs)"),
]


def half_radius(angles: np.ndarray, vals: np.ndarray) -> float:
    """Angle (degrees, log-linear interpolation) where the normalised median score first falls to
    0.5; NaN when it never does within the sweep."""
    for i in range(1, len(angles)):
        if vals[i] <= 0.5 < vals[i - 1]:
            f = (vals[i - 1] - 0.5) / (vals[i - 1] - vals[i])
            return float(
                np.exp(np.log(angles[i - 1]) + f * (np.log(angles[i]) - np.log(angles[i - 1])))
            )
    return float("nan")


def main() -> None:
    rec = json.loads((pc.CACHE / "separate_records.json").read_text())
    vox = sorted({r["voxel"] for r in rec})
    by_v: Dict[int, List[Dict[str, Any]]] = {v: [r for r in rec if r["voxel"] == v] for v in vox}
    truth = {v: next(r for r in by_v[v] if r["group"] == "truth") for v in vox}
    bflag = {v: truth[v]["boundary"] for v in vox}
    summary: Dict[str, Any] = {
        "n_voxels": len(vox),
        "n_boundary": int(sum(bflag.values())),
        "tol_deg": TOL,
    }
    summary["truth_score_median"] = {
        k: float(np.median([truth[v][k] for v in vox])) for k, _ in COSTS
    }
    summary["truth_score_min"] = {k: float(np.min([truth[v][k] for v in vox])) for k, _ in COSTS}
    angles = sorted({r["angle"] for r in rec if r["group"] == "random axis"})
    norm_rows = []
    curves: Dict[str, List[float]] = {k: [] for k, _ in COSTS}
    for a in angles:
        row: Dict[str, Any] = {"misorientation (deg)": a}
        for k, name in COSTS:
            per_v = []
            for v in vox:
                t = truth[v][k]
                m = [r[k] for r in by_v[v] if r["group"] == "random axis" and r["angle"] == a]
                per_v.append(np.mean(m) / t if t > 0 else np.nan)
            row[name] = float(np.nanmedian(per_v))
            curves[k].append(row[name])
        norm_rows.append(row)
    summary["normalised_median_by_angle"] = {
        k: dict(zip(map(str, angles), v)) for k, v in curves.items()
    }
    summary["half_score_radius_deg"] = {
        k: half_radius(np.array(angles), np.array(curves[k])) for k, _ in COSTS
    }
    # separation: best wrong orientation at >= 1 degree (random axis), lattice-axis rotations
    sep_rows = []
    sep: Dict[str, Any] = {}
    for gname, sel in (
        (
            "random axis, 1 degree and more",
            lambda r: r["group"] == "random axis" and r["angle"] >= 1.0,
        ),
        (
            "random axis, 0.25 to 0.5 degree",
            lambda r: r["group"] == "random axis" and 0.25 <= r["angle"] <= 0.5,
        ),
        ("rotations about <100>, <110>, <111>", lambda r: r["group"].startswith("about")),
    ):
        row = {"wrong orientations": gname}
        for k, name in COSTS:
            ratio = np.array(
                [
                    (
                        max(r[k] for r in by_v[v] if sel(r)) / truth[v][k]
                        if truth[v][k] > 0
                        else np.nan
                    )
                    for v in vox
                ]
            )
            row[name + ": median of best wrong / truth"] = float(np.nanmedian(ratio))
            row[name + ": voxels with best wrong >= 0.5 truth"] = int((ratio >= 0.5).sum())
            row[name + ": voxels with best wrong >= truth"] = int((ratio >= 1.0).sum())
            sep.setdefault(gname, {})[k] = {
                "median_ratio": float(np.nanmedian(ratio)),
                "n_ge_half": int((ratio >= 0.5).sum()),
                "n_ge_truth": int((ratio >= 1.0).sum()),
                "max_ratio": float(np.nanmax(ratio)),
            }
        sep_rows.append(row)
    summary["separation"] = sep
    tables = {}
    tables["step3_sharpness"] = DT.markdown_table(
        norm_rows,
        list(norm_rows[0].keys()),
        {k: ".3f" for k in norm_rows[0] if k != "misorientation (deg)"},
    )
    # compact separation table (one line per group x cost)
    lines = []
    for gname in sep:
        for k, name in COSTS:
            s = sep[gname][k]
            lines.append(
                {
                    "wrong orientations": gname,
                    "score": name,
                    "median best wrong / truth": s["median_ratio"],
                    "max over voxels": s["max_ratio"],
                    "voxels >= 0.5 truth": s["n_ge_half"],
                    "voxels >= truth": s["n_ge_truth"],
                }
            )
    tables["step3_separation"] = DT.markdown_table(
        lines,
        list(lines[0].keys()),
        {"median best wrong / truth": ".3f", "max over voxels": ".3f"},
    )
    trows = [
        {
            "score": name,
            "median at the truth": summary["truth_score_median"][k],
            "minimum at the truth": summary["truth_score_min"][k],
            "half-score radius (deg)": summary["half_score_radius_deg"][k],
        }
        for k, name in COSTS
    ]
    tables["step3_truth_scores"] = DT.markdown_table(
        trows,
        list(trows[0].keys()),
        {
            "median at the truth": ".3f",
            "minimum at the truth": ".3f",
            "half-score radius (deg)": ".3f",
        },
    )
    # sensitivity to the direction tolerance and the frame window (B footprint list)
    tol_rows = []
    for fw in (0, 1):
        for tol in (0.02, 0.05, 0.1, 0.2):
            k = f"B_c_{tol}_{fw}"
            ko = f"B_p_{tol}_{fw}"
            ks = f"oracle_p_{tol}_{fw}"
            curve = []
            for a in angles:
                per_v = []
                for v in vox:
                    t = truth[v][k]
                    m = [r[k] for r in by_v[v] if r["group"] == "random axis" and r["angle"] == a]
                    per_v.append(np.mean(m) / t if t > 0 else np.nan)
                curve.append(float(np.nanmedian(per_v)))
            tol_rows.append(
                {
                    "frame window": fw,
                    "direction tolerance (deg)": tol,
                    "B completeness at truth (median)": float(
                        np.median([truth[v][k] for v in vox])
                    ),
                    "oracle purity at truth (median)": float(
                        np.median([truth[v][ks] for v in vox])
                    ),
                    "B half-score radius (deg)": half_radius(np.array(angles), np.array(curve)),
                }
            )
    summary["tolerance_sensitivity"] = tol_rows
    tables["step3_tolerance"] = DT.markdown_table(
        tol_rows,
        list(tol_rows[0].keys()),
        {
            "B completeness at truth (median)": ".3f",
            "oracle purity at truth (median)": ".3f",
            "B half-score radius (deg)": ".3f",
        },
    )
    (OUT / "separate_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    DT.write_tables(OUT / "tables_separate.md", tables)
    for t in tables.values():
        print(t)
        print()


if __name__ == "__main__":
    main()
