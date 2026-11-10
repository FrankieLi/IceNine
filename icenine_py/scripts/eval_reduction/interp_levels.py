"""Task B (extra): pass-set parity of the Python screen against C++ at levels 0 to 3.

Task A measured the omega interpolation (CppOmegaInterp) against C++ at level 0 only (4,885 of 4,886
FZ points with identical pass sets for voxel 15901). Here a scratch C++ build prints, for every
level of the discrete search, each FZ point's search centre, the local grid and the peak_overlap of
each of its nine local-grid points (benchmarks/eval_reduction/B/levels/cpp_levels_patch.py.txt, env
L0TRACE=1). The Python screen is evaluated on the same centres and local grid, with the level's
Q-max (5, 6, 7, 8), pixel radius 3, the C++ overlap rule and the centroid screen, in four ways:
exact or interpolated omegas, with and without CppEdgeRules. A pair passes when peak_overlap > 0;
an FZ point's pass set is the set of passing local-grid indices. "overlap_equal_pairs" counts the
pairs whose peak_overlap count is equal in C++ and Python (stricter than equal pass/fail).

  uv run python scripts/eval_reduction/interp_levels.py --log DIR/log_15901.txt --voxel 15901 \
      --out DIR
"""

import argparse
import json
import sys
from pathlib import Path
from typing import Any, Dict, List

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import parity_common as PC  # noqa: E402


def parse_log(path: Path) -> Dict[int, Dict[str, Any]]:
    levels: Dict[int, Dict[str, Any]] = {}
    with open(path) as fh:
        for line in fh:
            t = line.split()
            if not t or not t[0].startswith("L0"):
                continue
            if t[0] == "L0H":
                c = int(t[2])
                levels[c] = dict(
                    n_fz=int(t[4]),
                    n_local=int(t[6]),
                    q_max=float(t[8]),
                    diameter_rad=float(t[10]),
                    G={},
                    C={},
                    O={},
                    nobs={},
                )
            elif t[0] == "L0G":
                levels[int(t[1])]["G"][int(t[2])] = [float.fromhex(x) for x in t[3:12]]
            elif t[0] == "L0C":
                levels[int(t[1])]["C"][int(t[2])] = [float.fromhex(x) for x in t[3:12]]
            elif t[0] == "L0O":
                c, i = int(t[1]), int(t[2])
                levels[c]["nobs"][i] = int(t[3])
                levels[c]["O"][i] = [int(x) for x in t[4:]]
    return levels


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--log", required=True)
    ap.add_argument("--voxel", type=int, required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    out = Path(a.out).resolve()
    log = Path(a.log).resolve()
    levels = parse_log(log)
    ctx = PC.PyContext()
    from icenine.batched_screen import ScreenEvaluator
    from icenine.cost_functions import VoxelCostFunction
    from icenine.sampling import generate_local_grid_multi_level

    s = ctx.setup
    verts, phase = ctx.verts(a.voxel), int(ctx.mic.voxels[a.voxel].phase)
    res: Dict[str, Any] = dict(voxel=a.voxel, levels=[])
    for c in sorted(levels):
        L = levels[c]
        n_fz, J = L["n_fz"], L["n_local"]
        assert len(L["C"]) == n_fz == len(L["O"]), (len(L["C"]), n_fz, len(L["O"]))
        cen = np.array([L["C"][i] for i in range(n_fz)], np.float64).reshape(n_fz, 3, 3)
        lg = np.array([L["G"][j] for j in range(J)], np.float64).reshape(J, 3, 3)
        ov = np.array([L["O"][i] for i in range(n_fz)], np.int64)
        py_lg = np.asarray(generate_local_grid_multi_level(L["diameter_rad"], 0, 1), np.float64)
        rec: Dict[str, Any] = dict(
            level=c,
            n_fz=n_fz,
            q_max=L["q_max"],
            cpp_pass_pairs=int((ov > 0).sum()),
            cpp_fz_with_pass=int((ov > 0).any(axis=1).sum()),
            local_grid_max_abs_diff=(
                float(np.abs(py_lg - lg).max()) if py_lg.shape == lg.shape else None
            ),
            arms={},
        )
        for edge in (False, True):
            cf = VoxelCostFunction(
                simulator=s.simulator,
                detector_list=s.detector_list,
                range_map=s.range_map,
                exp_data=s.exp_data,
                sample=s.sample,
                structure_list=s.structure_list,
                mode="hard",
                eta_limit=ctx.eta_limit,
                pixel_radius=3,
                max_q=L["q_max"],
                cpp_overlap_rule=True,
                center_screen=True,
                cpp_edge_rules=edge,
            )
            ev = ScreenEvaluator(cf)
            for kind in ("exact", "interp"):
                m = np.zeros((n_fz, J), np.int64)
                for st in range(0, n_fz, 32):
                    ch = cen[st : st + 32]
                    if kind == "interp":
                        m[st : st + len(ch)] = ev.overlaps_interp(ch, lg, verts, phase)[0]
                    else:
                        prod = np.matmul(lg[None], ch[:, None])
                        m[st : st + len(ch)] = ev.overlaps_exact(
                            prod.reshape(-1, 3, 3), verts, phase
                        )[0].reshape(len(ch), J)
                same = ((m > 0) == (ov > 0)).all(axis=1)
                rec["arms"][f"{kind}|edge{int(edge)}"] = dict(
                    fz_identical_pass_set=int(same.sum()),
                    pairs_pass_differs=int(((m > 0) != (ov > 0)).sum()),
                    pairs_cpp_only=int(((ov > 0) & (m == 0)).sum()),
                    pairs_py_only=int(((m > 0) & (ov == 0)).sum()),
                    overlap_equal_pairs=int((m == ov).sum()),
                    n_pairs=int(m.size),
                )
        res["levels"].append(rec)
        print(json.dumps(rec), flush=True)
    out.mkdir(parents=True, exist_ok=True)
    (out / f"levels_v{a.voxel}.json").write_text(json.dumps(res, indent=1) + "\n")


if __name__ == "__main__":
    main()
