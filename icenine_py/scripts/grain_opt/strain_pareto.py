"""Build the strain cost-against-accuracy table for the tracking page.

Reads the per-grain rows of the strain arms (``benchmarks/grain_opt/go5*/...json``) listed in
ARMS and writes ``benchmarks/grain_opt/strain_pareto.json``: per arm the mean objective
evaluations per grain (forward predictions of the grain's spots, every cost call and every
Gauss-Newton residual evaluation counted), mean wall time per grain, and per strain level (1e-4,
5e-4, 1e-3 tensor RMS) the median strain error and median voxel orientation error. Arms whose
result file does not exist yet are skipped, so a new arm is one line here.

Usage (from icenine_py/): uv run python scripts/grain_opt/strain_pareto.py
"""

import json
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
B = HERE.parents[1] / "benchmarks" / "grain_opt"
LEVELS = (1e-4, 5e-4, 1e-3)

# (result file, label, family, noisy images)
# families: ref (strain only at the ground-truth orientation), coverage, centroid, hybrid, gn, fast
ARMS = [
    ("go5b/E_GT.json", "coverage, strain only at ground-truth orientation", "ref", False),
    ("go5f/E_GT_cen.json", "centroid, strain only at ground-truth orientation", "ref", False),
    ("go5g/E_GT_hyb2.json", "hybrid lambda 2, strain only at ground-truth orientation", "ref", False),
    ("go5c/J_P.json", "coverage, one long joint fit", "coverage", False),
    ("go5h/S_P_S36.json", "coverage, staged", "coverage", False),
    ("go5f/S_P_cen.json", "centroid, staged", "centroid", False),
    ("go5h/S_P_hyb8_S36.json", "hybrid lambda 8, staged (CMA)", "hybrid", False),
    ("go5h/S_P_hyb16_S36.json", "hybrid lambda 16, staged (CMA)", "hybrid", False),
    ("go5i/GN_P_hyb8_S36.json", "Gauss-Newton on centroids", "gn", False),
    ("go5i/GNH_P_hyb8_S36.json", "Gauss-Newton + hybrid polish 300", "gn", False),
    ("go5j/GA2_P_hyb8_S36.json", "Gauss-Newton + 2 gradient refits", "gn", False),
    ("go5k/GFF_P_hyb8_S36.json", "Gauss-Newton, frozen weights + fallback", "gn", False),
    ("go5l/GP4_P_hyb8_S36.json", "fast fit GP4@P", "fast", False),
    ("go5l/GP6_P_hyb8_S36.json", "fast fit GP6@P (polish 600)", "fast", False),
    ("go5o_clean/GP4_P_hyb8_S36_intclean.json", "fast fit, intclean centroids", "fast", False),
    ("go5m/S_P_hyb8_S36.json", "hybrid lambda 8, staged (CMA), noisy", "hybrid", True),
    ("go5m/GP4_P_hyb8_S36.json", "fast fit GP4@P, noisy", "fast", True),
    ("go5n/GP4_P_hyb8_S36_bin.json", "fast fit, binary centroids, noisy", "fast", True),
    ("go5o/GP4_P_hyb8_S36_intclean.json", "fast fit, intclean centroids, noisy", "fast", True),
    ("go7b/strain_psf10/GP4_P_hyb8_S36.json", "fast fit, 1 px PSF images", "fast", True),
    ("go7c/GP4_P_hyb8_S36_psf1.json", "fast fit + PSF-aware coverage, 1 px PSF images", "fast", True),
]


def main() -> None:
    pts = []
    for path, label, fam, noisy in ARMS:
        p = B / path
        if not p.exists():
            continue
        rows = json.loads(p.read_text())["rows"]
        lv = [[r for r in rows if abs(r["level"] - L) < 1e-12] for L in LEVELS]
        err = [float(np.median([r["err_rms"] for r in x])) for x in lv]
        mis = [float(np.median([r["mis_fit_med"] for r in x])) for x in lv]
        pts.append(
            dict(
                label=label,
                family=fam,
                noisy=noisy,
                file=path,
                n_grains=len(rows),
                evals=float(np.mean([r["evals"] for r in rows])),
                t_s=float(np.mean([r["t_s"] for r in rows])),
                err=err,
                mis=mis,
                orient_fail=fam != "ref" and max(mis) > 0.02,
            )
        )
    out = B / "strain_pareto.json"
    out.write_text(json.dumps(dict(levels=LEVELS, points=pts), indent=1) + "\n")
    print(out, len(pts), "arms")


if __name__ == "__main__":
    main()
