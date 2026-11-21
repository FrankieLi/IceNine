"""Trace the C++ client log against the truth, voxel by voxel (cpptime, review fixes 8 and 9).

  uv run python scripts/autoresearch/cpptime/trace_cpp_log.py RUN_DIR [--out FILE]

Every ``QueueSize`` record of ``ParallelReconLog_1.log`` is one local optimization of one voxel:
the log gives the voxel centre, ``fBestConf`` and the voxel's hit ratio, the seed orientation
(``Original``) and the orientation after the fit (``Optimized``). The C++ client marks the voxel
REFIT when hit ratio / fBestConf is at most 0.9 (``BreadthFirstReconstructor.tmpl.cpp:163-172``),
otherwise FITTED. The log centres are matched to the pilot grid by centroid; voxels are scored
with ``score_cpp``. Counts are for this run only; no cause is tested.
"""

import argparse
import json
import re
import sys
from pathlib import Path
from typing import Any, Dict, List

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import score_cpp as SC  # noqa: E402

REC = re.compile(
    r"\|=\|\s*(\S+)\s+(\S+)\s+Original :\s*(\S+\s+\S+\s+\S+)\s+Optimized:\s*(\S+\s+\S+\s+\S+)"
    r"\s+VoxelLocation\s+(\S+)\s+(\S+)\s+(\S+)"
)
FIRST = re.compile(r"First Voxel Required")
ACCEPT_RATIO = 0.9


def parse_records(text: str) -> List[Dict[str, Any]]:
    """The local-optimization records of a client log, with the index of their seed."""
    out: List[Dict[str, Any]] = []
    seed = -1
    for ln in text.splitlines():
        if FIRST.search(ln):
            seed += 1
            continue
        m = REC.search(ln)
        if not m:
            continue
        best, hit = float(m.group(1)), float(m.group(2))
        out.append(
            dict(
                seed=seed,
                best=best,
                hit=hit,
                ratio=hit / best,
                refit=hit / best <= ACCEPT_RATIO,
                original=m.group(3).split(),
                optimized=m.group(4).split(),
                x=float(m.group(5)),
                y=float(m.group(6)),
            )
        )
    return out


def summarize(recs: List[Dict[str, Any]], wrong_by_centroid: Dict[int, bool], row_of: List[int]):
    """Counts of the refit and accepted records among the wrong voxels (pure helper)."""
    refit = [r["refit"] for r in recs]
    wrong = [wrong_by_centroid[row] for row in row_of]
    return dict(
        records=len(recs),
        refit=int(sum(refit)),
        refit_wrong=int(sum(w for r, w in zip(refit, wrong) if r)),
        accepted_wrong=int(sum(w for r, w in zip(refit, wrong) if not r)),
        unchanged=int(sum(r["original"] == r["optimized"] for r in recs)),
    )


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("run", type=Path)
    ap.add_argument("--out", type=Path, default=None)
    a = ap.parse_args()
    from scipy.spatial import cKDTree

    from icenine.orientation_search import get_symmetry_quaternions
    from icenine.symmetry import create_cubic_symmetry
    from make_region import centroids_from_lines
    from score_bfs import mats_to_quats, pair_misorientation_deg

    recs = parse_records((a.run / "ParallelReconLog_1.log").read_text())
    grid_lines = (SC.BENCH / "region_grid.mic").read_text().splitlines()
    raw = np.array([[float(t) for t in ln.split()[:6]] for ln in grid_lines[1:] if ln.strip()])
    cent, _ = centroids_from_lines(raw, float(grid_lines[0].split()[0]))
    tree = cKDTree(cent)
    xy = np.array([[r["x"], r["y"]] for r in recs])
    dist, row_of = tree.query(xy)
    assert dist.max() < 1e-4, f"log centre not on a grid centroid (max {dist.max()})"

    # scored output voxels (rows of the pilot grid) and their errors
    full = SC.read_positions(SC.BENCH / "region_grid.mic")
    out = SC.read_positions(a.run / SC.MIC_NAME)
    idx, _ = SC.match_voxels(full, out)
    q_true = SC.truth_quats()
    err = np.full(len(full), np.nan)
    err[idx] = SC._errors(a.run / SC.MIC_NAME, q_true, idx)
    wrong_row = {int(i): bool(err[i] > SC.WRONG_DEG) for i in range(len(full))}
    row_list = [int(r) for r in row_of]
    res: Dict[str, Any] = summarize(recs, wrong_row, row_list)
    res["distinct_record_voxels"] = len(set(row_list))
    res["seeds_in_log"] = (max(r["seed"] for r in recs) + 1) if recs else 0
    seed_rows = sorted(set(range(len(full))) - set(row_list))
    res["voxels_not_in_records"] = len(seed_rows)
    res["wrong_among_voxels_not_in_records"] = int(sum(wrong_row[i] for i in seed_rows))
    res["wrong_total"] = int(np.nansum(err > SC.WRONG_DEG))

    # every wrong voxel: is its orientation within 1 degree of another grain's truth?
    sym = get_symmetry_quaternions(create_cubic_symmetry(1.0))
    from icenine.mic_file import MicFile

    mic = MicFile.read(str(a.run / SC.MIC_NAME))
    q_out = mats_to_quats(np.array([v.orientation for v in mic.voxels]))
    grain = np.load(SC.TRUTH / "rand_500grains_1mm_neworient_s0_grainmap.npy")[
        np.load(SC.BENCH / "region_index.npy")
    ]
    near_other = 0
    n_wrong = 0
    for k, row in enumerate(idx):
        if err[row] <= SC.WRONG_DEG:
            continue
        n_wrong += 1
        e = pair_misorientation_deg(np.repeat(q_out[k : k + 1], len(full), axis=0), q_true, sym)
        near_other += int(bool(((e <= SC.WRONG_DEG) & (grain != grain[row])).any()))
    res["wrong_near_another_grain_truth"] = near_other
    res["wrong_checked"] = n_wrong
    cols = np.array([ln.split()[6:9] for ln in (a.run / SC.MIC_NAME).read_text().splitlines()[1:]])
    res["distinct_output_orientations"] = int(len({tuple(c) for c in cols.tolist()}))
    text = json.dumps(res, indent=1) + "\n"
    if a.out:
        a.out.write_text(text)
    print(text)
    return 0


if __name__ == "__main__":
    sys.exit(main())
