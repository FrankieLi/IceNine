"""Score a C++ pilot run against the truth (cpptime).

  uv run python scripts/autoresearch/cpptime/score_cpp.py RUN_DIR [--python-run DIR] [--out F]
      [--region-grid FILE]          (default: the pilot region grid)

Matches the voxels of ``RUN_DIR/cpp_pilot.mic.LBFS`` to the rows of the region grid by centroid and
up/down within 1e-3 mm (every output voxel must match exactly one row; region rows without an
output voxel are "missing", counted apart and not wrong). Wrong = more than 1 degree from the truth
(``rand_500grains_1mm_neworient_s0.mic`` at ``region_index.npy``), cubic symmetry. With
``--python-run`` the stored Python ``recon.mic`` is scored the same way and the paired split
(both wrong, only C++, only Python, neither) is given over the voxels C++ reached.
"""

import argparse
import json
import sys
from pathlib import Path
from typing import Any, Dict, Optional, Tuple

import numpy as np

ICE = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ICE / "scripts" / "phase_d"))
sys.path.insert(0, str(ICE / "scripts" / "common"))

from stats import wilson  # noqa: E402

MATCH_TOL_MM = 1e-3
WRONG_DEG = 1.0
BENCH = ICE / "benchmarks" / "phase_d_pilot"
TRUTH = ICE.parent / "Examples" / "Example2.ManyGrains" / "SimInput"
MIC_NAME = "cpp_pilot.mic.LBFS"


def read_positions(path: Path) -> np.ndarray:
    """(N, 3) array [x, y, direction] of the voxel rows of a .mic file (header line skipped)."""
    rows = []
    for ln in path.read_text().splitlines()[1:]:
        t = ln.split()
        if len(t) >= 4:
            rows.append([float(t[0]), float(t[1]), float(t[3])])
    return np.array(rows, dtype=float).reshape(-1, 3)


def match_voxels(
    region: np.ndarray, out: np.ndarray, tol: float = MATCH_TOL_MM
) -> Tuple[np.ndarray, np.ndarray]:
    """Region row of each output voxel, and the region rows with no output voxel.

    Both arrays are (N, 3) [x, y, direction]. Asserts that every output voxel matches exactly one
    unused region row."""
    from scipy.spatial import cKDTree

    tree = cKDTree(region[:, :2])
    idx = np.full(len(out), -1, dtype=int)
    for i, (x, y, d) in enumerate(out):
        cand = [j for j in tree.query_ball_point([x, y], tol) if region[j, 2] == d]
        assert len(cand) == 1, f"output voxel {i} at ({x}, {y}) matches {len(cand)} region rows"
        idx[i] = cand[0]
    assert len(set(idx.tolist())) == len(idx), "two output voxels match the same region row"
    missing = np.setdiff1d(np.arange(len(region)), idx)
    return idx, missing


def paired_wrong(py_wrong: np.ndarray, cpp_wrong: np.ndarray) -> Dict[str, int]:
    """Counts of voxels wrong in both, only in C++, only in Python, in neither."""
    py, cp = np.asarray(py_wrong, bool), np.asarray(cpp_wrong, bool)
    return dict(
        both=int((py & cp).sum()),
        only_cpp=int((~py & cp).sum()),
        only_python=int((py & ~cp).sum()),
        neither=int((~py & ~cp).sum()),
    )


def _errors(mic_path: Path, truth_q: np.ndarray, rows: np.ndarray) -> np.ndarray:
    """Misorientation (deg) of each voxel of ``mic_path`` to ``truth_q[rows]``."""
    from icenine.mic_file import MicFile
    from icenine.orientation_search import get_symmetry_quaternions
    from icenine.symmetry import create_cubic_symmetry
    from score_bfs import mats_to_quats, pair_misorientation_deg

    sym = get_symmetry_quaternions(create_cubic_symmetry(1.0))
    mic = MicFile.read(str(mic_path))
    q = mats_to_quats(np.array([v.orientation for v in mic.voxels]))
    return pair_misorientation_deg(q, truth_q[rows], sym)


def score_run(
    run: Path, python_run: Optional[Path] = None, region_grid: Path = BENCH / "region_grid.mic"
) -> Dict[str, Any]:
    from icenine.mic_file import MicFile
    from score_bfs import mats_to_quats

    region = read_positions(region_grid)
    out = read_positions(run / MIC_NAME)
    idx, missing = match_voxels(region, out)
    truth_idx = np.load(BENCH / "region_index.npy")
    assert len(truth_idx) == len(region)
    truth = MicFile.read(str(TRUTH / "rand_500grains_1mm_neworient_s0.mic"))
    q_all = mats_to_quats(np.array([truth.voxels[i].orientation for i in truth_idx]))
    err = _errors(run / MIC_NAME, q_all, idx)
    wrong = err > WRONG_DEG
    lo, hi = wilson(int(wrong.sum()), len(err))
    res: Dict[str, Any] = dict(
        n_region=int(len(region)),
        n_output=int(len(out)),
        n_missing=int(len(missing)),
        wrong=int(wrong.sum()),
        wrong_frac=float(wrong.mean()) if len(err) else float("nan"),
        wrong_wilson=[lo, hi],
        median_err_deg=float(np.median(err)) if len(err) else float("nan"),
    )
    if python_run is not None:
        perr = _errors(python_run / "recon.mic", q_all, np.arange(len(region)))
        pw = perr > WRONG_DEG
        res["python_wrong_all"] = int(pw.sum())
        res["python_wrong_on_cpp_voxels"] = int(pw[idx].sum())
        res["paired_on_cpp_voxels"] = paired_wrong(pw[idx], wrong)
    return res


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("run", type=Path)
    ap.add_argument("--python-run", type=Path, default=None)
    ap.add_argument("--region-grid", type=Path, default=BENCH / "region_grid.mic")
    ap.add_argument("--out", type=Path, default=None)
    a = ap.parse_args()
    res = score_run(a.run, a.python_run, a.region_grid)
    text = json.dumps(res, indent=1) + "\n"
    if a.out:
        a.out.write_text(text)
    print(text)
    return 0


if __name__ == "__main__":
    sys.exit(main())
