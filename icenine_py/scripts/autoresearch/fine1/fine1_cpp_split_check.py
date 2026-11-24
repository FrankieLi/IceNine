"""fine1: the fine grid against C++'s own split (the grid check of the plan, run on the real code).

  uv run python scripts/autoresearch/fine1/fine1_cpp_split_check.py --binary <cpp_split_check>

``cpp_split_check.cpp`` (built as its header says; links the stock build's
``ReconstructionStrategies.cpp.o``) reads a grid and calls C++'s ``SetMinimumResolution`` with
``MinSideLength`` 0.0046875, then writes the voxel list. Two runs: (1) on the generation-6 source
grid (all 24,570 voxels, the same call the C++ reconstruction makes): its output must equal
``fine7_grid.mic`` row by row (position within 1e-6 mm, direction, generation); (2) on
``fine7_grid.mic`` itself: nothing is split again (98,280 voxels out, the same rows). Writes
``benchmarks/autoresearch/fine1/cpp_split_check.json``; exits 1 on any difference. No
reconstruction runs; the binary only reads and writes a grid.
"""

import argparse
import json
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import Any, Dict

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import fine1_paths as P  # noqa: E402


def rows(path: Path) -> np.ndarray:
    """(n, 5): x, y, direction, generation, phase of every row of a .mic file."""
    a = np.loadtxt(str(path), skiprows=1, ndmin=2)
    return a[:, [0, 1, 3, 4, 5]]


def compare(a: np.ndarray, b: np.ndarray) -> Dict[str, Any]:
    out: Dict[str, Any] = dict(n_a=len(a), n_b=len(b))
    if len(a) != len(b):
        out["ok"] = False
        return out
    out["max_abs_position_diff_mm"] = float(np.abs(a[:, :2] - b[:, :2]).max())
    out["n_direction_differs"] = int((a[:, 2] != b[:, 2]).sum())
    out["n_generation_differs"] = int((a[:, 3] != b[:, 3]).sum())
    out["ok"] = bool(
        out["max_abs_position_diff_mm"] < 1e-6
        and out["n_direction_differs"] == 0
        and out["n_generation_differs"] == 0
    )
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--binary", type=Path, required=True)
    a = ap.parse_args()
    mine = rows(P.GRID)
    res: Dict[str, Any] = dict(min_side_length=P.SIDE_FINE_MM)
    with tempfile.TemporaryDirectory() as td:
        for name, src in (("source_grid_split_by_cpp", P.SRC_GRID), ("fine_grid_again", P.GRID)):
            out = Path(td) / f"{name}.mic"
            r = subprocess.run(
                [str(a.binary), str(src), str(P.SIDE_FINE_MM), str(out)],
                capture_output=True,
                text=True,
            )
            res[name] = dict(
                rc=r.returncode,
                stdout=r.stdout.strip().splitlines()[-2:],
                **(compare(rows(out), mine) if r.returncode == 0 else dict(ok=False)),
            )
    res["binary_sha256"] = P.sha256(a.binary)
    res["ok"] = bool(res["source_grid_split_by_cpp"]["ok"] and res["fine_grid_again"]["ok"])
    (P.BENCH / "cpp_split_check.json").write_text(json.dumps(res, indent=1) + "\n")
    print(json.dumps(res, indent=1))
    return 0 if res["ok"] else 1


if __name__ == "__main__":
    sys.exit(main())
