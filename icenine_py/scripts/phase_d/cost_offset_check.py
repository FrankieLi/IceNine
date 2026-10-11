"""Does the opt-in C++ overlap rule remove the Python-minus-C++ cost offset (-0.013 on average)?

Re-evaluates the 160 start orientations whose C++ cost was printed by the debug build
(benchmarks/phase_d_seed_diag/variance_parity/cpp_trace2.log, voxel 15901, clean full-sample
images) with the local cost function, rule off and on. Writes cost_offset_check.json.

  uv run python scripts/phase_d/cost_offset_check.py
"""

import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ICE = HERE.parents[1]
sys.path.insert(0, str(HERE))
import seed_diag_variance_parity as P  # noqa: E402

LOG = ICE / "benchmarks" / "phase_d_seed_diag" / "variance_parity" / "cpp_trace2.log"
OUT = ICE / "benchmarks" / "phase_d_seed_diag" / "level0" / "cost_offset_check.json"


def stats(py: np.ndarray, cpp: np.ndarray) -> dict:
    d = py - cpp
    return dict(
        n=len(d),
        n_within_1e6=int((np.abs(d) < 1e-6).sum()),
        n_within_1e3=int((np.abs(d) < 1e-3).sum()),
        mean_signed=float(d.mean()),
        median_signed=float(np.median(d)),
        min_signed=float(d.min()),
        max_signed=float(d.max()),
        n_python_lower=int((d < -1e-6).sum()),
        n_python_higher=int((d > 1e-6).sum()),
        corr=float(np.corrcoef(py, cpp)[0, 1]),
    )


def main() -> None:
    rec, local_off, verts, phase, truth = P.setup(15901)
    vm, _ = P.parse_cpp(LOG)
    cpp = np.array([c for _, c in vm])
    rec.setup.search_params.cpp_overlap_rule = True
    local_on = rec._make_local_cost_fn()
    off = np.array([local_off.evaluate(R.astype(np.float32), verts, phase).cost for R, _ in vm])
    on = np.array([local_on.evaluate(R.astype(np.float32), verts, phase).cost for R, _ in vm])
    res = dict(voxel=15901, rule_off=stats(off, cpp), rule_on=stats(on, cpp))
    OUT.write_text(json.dumps(res, indent=1) + "\n")
    print(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()
