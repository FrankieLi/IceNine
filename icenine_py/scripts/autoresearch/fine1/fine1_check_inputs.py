"""fine1: do the C++ arms and the Python arms read the same inputs?

  uv run python scripts/autoresearch/fine1/fine1_check_inputs.py

The same questions as part 1's input check (``cppchk/check_inputs.py``, on the unmerged branch
``feature/eval-reduction-cppchk``, so it is rewritten here for the fine layer, not imported):
for the whole layer and the slice, (1) the C++ config (arm CS) against the Python config of arm
PS, key by key after paths are made absolute: the listed differences must be exactly the three
Python-only keys, the output basename and the save interval (the plan: only the grid, the images,
the two side lengths and the output name differ from the base); (2) the sha256 of every file
named by a path key, read at the path the Python run resolves (relative to
``Examples/Example2.ManyGrains``) and at the path in the C++ config; (3) the grid file and its
voxel count; (4) the image frames C++ reads against the stack Python reads (``stack_check.json``).
Writes ``benchmarks/autoresearch/fine1/inputs_check.json``; exits 1 if a difference is not on
the expected list.
"""

import json
import sys
from pathlib import Path
from typing import Any, Dict

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import fine1_config as CF  # noqa: E402
import fine1_paths as P  # noqa: E402

EXPECTED = {
    "BFSRevisitRefit",
    "BFSReseed",
    "FusedCost",  # the three Python-only keys
    "OutStructureBasename",
    "SecondsBetweenSave",
}


def one(grid: Path, n_expect: int) -> Dict[str, Any]:
    frames = str(P.FRAMES / "500Grains.sim")
    py = CF.keys_of(CF.python_config("PS", grid, frames, "x"))
    cpp_text = CF.cpp_config(grid, frames, "y")
    cpp = CF.keys_of(cpp_text)
    diff: Dict[str, Any] = {}
    for k in sorted(set(py) | set(cpp)):
        a, b = py.get(k), cpp.get(k)
        if a is None or b is None:
            diff[k] = dict(python=a, cpp=b)
        elif a != b:
            # a path key: same file if the resolved paths agree
            pa = Path(a) if Path(a).is_absolute() else P.EX / a
            if pa.resolve() != Path(b).resolve():
                diff[k] = dict(python=a, cpp=b)
    files: Dict[str, Any] = {}
    for k in CF.MC.PATH_KEYS:
        if k == "InfileBasename":
            continue
        pa = Path(py[k]) if Path(py[k]).is_absolute() else P.EX / py[k]
        pb = Path(cpp[k])
        files[k] = dict(
            python_sha256=P.sha256(pa), cpp_sha256=P.sha256(pb), same=P.sha256(pa) == P.sha256(pb)
        )
    rows = [ln for ln in grid.read_text().splitlines()[1:] if ln.strip()]
    return dict(
        grid=str(grid.relative_to(P.BENCH)),
        grid_sha256=P.sha256(grid),
        n_voxels=len(rows),
        n_voxels_expected=n_expect,
        sample_filename_same=py["SampleFilename"] == cpp["SampleFilename"],
        infile_basename_same=py["InfileBasename"] == cpp["InfileBasename"],
        side_lengths=[cpp["MinSideLength"], cpp["MaxInitSideLength"]],
        differences=diff,
        differences_unexpected=sorted(set(diff) - EXPECTED),
        expected_missing=sorted(EXPECTED - set(diff)),
        path_key_files=files,
        path_key_files_all_same=all(v["same"] for v in files.values()),
    )


def match_uniqueness(grid: Path, tol_mm: float = 1e-3) -> Dict[str, Any]:
    """``score_cpp.py`` matches an output voxel to the grid row of the same direction within
    ``tol_mm`` of its left vertex (and needs exactly one). On the fine grid: the nearest other
    row of the same direction, and the number of rows with a second candidate in the tolerance."""
    import numpy as np
    from scipy.spatial import cKDTree

    a = np.loadtxt(str(grid), skiprows=1, ndmin=2)
    xy, d = a[:, :2], a[:, 3]
    n_amb = 0
    nearest = np.inf
    for code in (1, 2):
        pts = xy[d == code]
        dist, _ = cKDTree(pts).query(pts, k=2)
        nearest = min(nearest, float(dist[:, 1].min()))
        n_amb += int((dist[:, 1] <= tol_mm).sum())
    return dict(
        tolerance_mm=tol_mm,
        nearest_same_direction_mm=nearest,
        n_rows_with_second_candidate=n_amb,
        centre_distance_of_edge_neighbours_mm=P.SIDE_FINE_MM / 3**0.5,
        ok=bool(n_amb == 0),
    )


def main() -> int:
    stack = json.loads(P.STACK_JSON.read_text())
    out: Dict[str, Any] = dict(
        whole=one(P.GRID, P.N_FINE),
        slice=one(P.SLICE, P.SLICE_N),
        frames_equal_stack=bool(stack["ok"]),
        lit_pixels=stack["lit_pixels_ascii"],
    )
    out["score_match"] = match_uniqueness(P.GRID)
    ok = bool(out["frames_equal_stack"] and out["score_match"]["ok"])
    for part in ("whole", "slice"):
        r = out[part]
        ok = ok and (
            r["n_voxels"] == r["n_voxels_expected"]
            and r["sample_filename_same"]
            and r["infile_basename_same"]
            and not r["differences_unexpected"]
            and r["path_key_files_all_same"]
        )
    out["ok"] = ok
    (P.BENCH / "inputs_check.json").write_text(json.dumps(out, indent=1) + "\n")
    print(
        json.dumps({k: (v if k in ("ok", "frames_equal_stack") else "...") for k, v in out.items()})
    )
    for part in ("whole", "slice"):
        r = out[part]
        print(part, r["n_voxels"], "differences:", sorted(r["differences"]))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
