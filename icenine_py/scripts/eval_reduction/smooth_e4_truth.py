"""Task E review fix 7: what makes the E4 cost high at the true orientation.

E4's cost is the weighted mean of a per-peak residual ``r^2 / (r^2 + tau^2)``; a predicted peak
with no attached observed ray within one frame scores exactly 1. At the truth of each voxel of
Measurement 1 (the six full-sample voxels and the 30 E0 voxels, Task B's sets), this script splits
the E4 cost (tau 0.05 deg, the landscape's E4t05) into the part from peaks with no attached ray
nearby and the part from peaks that have one, and counts the predicted peaks of each kind. The
split is arithmetic (the two parts sum to the cost); it does not say why a peak has no ray. For the
peaks with a ray it also gives the residual (angle plus omega distance) to the nearest one.

  uv run python scripts/eval_reduction/smooth_e4_truth.py
Writes benchmarks/eval_reduction/E/e4_truth/{summary.json, tables.md}.
"""

import json
import math
import os
import sys
from pathlib import Path
from typing import Any, Dict, List

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")

import numpy as np

HERE = Path(__file__).resolve().parent
ICE = HERE.parents[1]
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(ICE / "scripts" / "common"))
import cost_edges_run as RUN  # noqa: E402
import doc_tables as DT  # noqa: E402

OUT = ICE / "benchmarks" / "eval_reduction" / "E" / "e4_truth"
TAU_DEG = 0.05


def split(sc: Any, Rt: np.ndarray, verts: Any, phase: int) -> Dict[str, float]:
    """E4 cost at ``Rt`` and its two parts (peaks with no attached ray nearby / with one)."""
    from icenine.smooth_cost import SmoothSpec

    spec = SmoothSpec.parse(f"E4,tau={TAU_DEG}")
    pr = sc.predict(Rt[None], verts, phase)
    W, rho = sc._e4_terms(pr, verts, math.radians(TAU_DEG), spec)
    W, rho = W.detach().numpy(), rho.detach().numpy()
    sw = float(W.sum())
    none = rho == 1.0  # no attached ray within one frame (a finite residual scores below 1)
    cost = float((W * rho).sum() / sw)
    # residual of the peaks with a ray, from rho = r^2 / (r^2 + tau^2)
    has = (W > 0) & ~none
    r_deg = TAU_DEG * np.sqrt(rho[has] / np.clip(1.0 - rho[has], 1e-300, None))
    check = float(sc.evaluate_many(Rt[None], verts, phase, spec)[0])
    assert abs(cost - check) < 1e-9, (cost, check)
    return dict(
        cost=cost,
        part_no_ray=float((W * none).sum() / sw),
        part_with_ray=float((W * rho * ~none).sum() / sw),
        peaks=int((W > 0).sum()),
        peaks_no_ray=int(((W > 0) & none).sum()),
        attached=int(len(sc.attached_rays(verts))),
        median_residual_deg=float(np.median(r_deg)) if has.any() else float("nan"),
    )


def run(kind: str) -> List[Dict[str, Any]]:
    import torch

    from icenine.smooth_cost import SmoothCost

    torch.set_num_threads(1)
    (RUN.init_full if kind == "full" else RUN.init_e0)()
    rows = []
    sc = None
    for v in RUN.cases(kind):
        base, Rt, verts, phase = RUN.case_data(v)
        if sc is None or kind == "e0":  # E0: every voxel has its own images
            sc = SmoothCost(
                base.simulator,
                base.detector_list,
                base.range_map,
                base.exp_data,
                base.sample,
                base.structure_list,
                eta_limit=base.eta_limit,
            )
        x = split(sc, np.asarray(Rt, dtype=np.float64), verts, phase)
        x["voxel"] = int(v)
        rows.append(x)
        print(kind, json.dumps(x), flush=True)
    return rows


def summary(rows: List[Dict[str, Any]]) -> Dict[str, Any]:
    def med(k: str) -> float:
        return float(np.median([r[k] for r in rows]))

    share = [r["peaks_no_ray"] / r["peaks"] for r in rows]
    return dict(
        n_voxels=len(rows),
        median_cost=med("cost"),
        median_part_no_ray=med("part_no_ray"),
        median_part_with_ray=med("part_with_ray"),
        median_peaks=med("peaks"),
        median_share_peaks_no_ray=float(np.median(share)),
        min_share_peaks_no_ray=float(np.min(share)),
        max_share_peaks_no_ray=float(np.max(share)),
        median_attached=med("attached"),
        median_residual_deg=med("median_residual_deg"),
        voxels=rows,
    )


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    S = {kind: summary(run(kind)) for kind in ("full", "e0")}
    (OUT / "summary.json").write_text(json.dumps(S, indent=1) + "\n")
    names = {"full": "six full-sample voxels", "e0": "30 E0 voxels"}
    trows = []
    for kind, s in S.items():
        trows.append(
            {
                "voxel set": names[kind],
                "n": s["n_voxels"],
                "E4 cost at the truth, tau 0.05 deg (median)": s["median_cost"],
                "part from peaks with no observed ray nearby (median)": s["median_part_no_ray"],
                "part from peaks with a ray (median)": s["median_part_with_ray"],
                "predicted peaks (median)": s["median_peaks"],
                "share of peaks with no ray nearby, median (range)": (
                    f"{s['median_share_peaks_no_ray']:.1%} ({s['min_share_peaks_no_ray']:.1%}"
                    f"-{s['max_share_peaks_no_ray']:.1%})"
                ),
                "residual to the nearest ray, peaks with a ray (median over voxels of the "
                "voxel median, deg)": s["median_residual_deg"],
                "observed rays attached (median)": s["median_attached"],
            }
        )
    fm = {k: ".4f" for k in trows[0] if "(median)" in k and "peaks (median)" not in k}
    fm.update({"predicted peaks (median)": ".1f", "observed rays attached (median)": ".1f"})
    fm[[k for k in trows[0] if k.startswith("residual")][0]] = ".4f"
    T = {"evalE_e4_truth": DT.markdown_table(trows, list(trows[0].keys()), fm)}
    DT.write_tables(OUT / "tables.md", T)
    print(T["evalE_e4_truth"])


if __name__ == "__main__":
    main()
