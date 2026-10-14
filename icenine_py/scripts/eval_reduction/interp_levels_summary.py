"""Task B (extra): tables of the per-level pass-set parity (interp_levels.py output, six voxels).

uv run python scripts/eval_reduction/interp_levels_summary.py --dir
benchmarks/eval_reduction/B/levels
"""

import argparse
import json
import sys
from pathlib import Path
from typing import Any, Dict, List

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[0] / "common"))
import doc_tables as DT  # noqa: E402

ARMS = [
    ("interp|edge1", "Python interpolated omegas + CppEdgeRules (both switches)"),
    ("interp|edge0", "Python interpolated omegas, edge rules off"),
    ("exact|edge1", "Python exact omegas + CppEdgeRules"),
    ("exact|edge0", "Python exact omegas, edge rules off (default)"),
]


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--dir", required=True)
    a = ap.parse_args()
    d = Path(a.dir).resolve()
    files = sorted(d.glob("levels_v*.json"))
    data = [json.loads(f.read_text()) for f in files]
    voxels = [x["voxel"] for x in data]
    summ: Dict[str, Any] = dict(voxels=voxels, levels={}, per_voxel={})
    rows: List[Dict[str, Any]] = []
    rows_v: List[Dict[str, Any]] = []
    for lvl in range(4):
        recs = [next(r for r in x["levels"] if r["level"] == lvl) for x in data]
        n_fz = sum(r["n_fz"] for r in recs)
        n_pairs = sum(r["arms"]["interp|edge1"]["n_pairs"] for r in recs)
        cpp_pass = sum(r["cpp_pass_pairs"] for r in recs)
        summ["levels"][lvl] = dict(
            n_fz=n_fz, n_pairs=n_pairs, cpp_pass_pairs=cpp_pass, q_max=recs[0]["q_max"], arms={}
        )
        grid = max(r["local_grid_max_abs_diff"] for r in recs)
        summ["levels"][lvl]["local_grid_max_abs_diff"] = grid
        for key, label in ARMS:
            s = {
                k: sum(r["arms"][key][k] for r in recs)
                for k in (
                    "fz_identical_pass_set",
                    "pairs_pass_differs",
                    "pairs_cpp_only",
                    "pairs_py_only",
                    "overlap_equal_pairs",
                )
            }
            summ["levels"][lvl]["arms"][key] = s
            rows.append(
                dict(
                    level=lvl,
                    q_max=recs[0]["q_max"],
                    arm=label,
                    fz=f"{s['fz_identical_pass_set']} / {n_fz}",
                    pairs_differ=s["pairs_pass_differs"],
                    cpp_only=s["pairs_cpp_only"],
                    py_only=s["pairs_py_only"],
                    equal=f"{s['overlap_equal_pairs']} / {n_pairs}",
                )
            )
    for x in data:
        for key, label in ARMS[:1] + ARMS[3:]:
            row = dict(voxel=x["voxel"], arm=label)
            for r in x["levels"]:
                s = r["arms"][key]
                row[f"L{r['level']}"] = f"{s['fz_identical_pass_set']} / {r['n_fz']}"
            rows_v.append(row)
    t1 = DT.markdown_table(
        rows,
        ["level", "q_max", "arm", "fz", "pairs_differ", "cpp_only", "py_only", "equal"],
        labels=None,
    )
    t2 = DT.markdown_table(rows_v, ["voxel", "arm", "L0", "L1", "L2", "L3"])
    DT.write_tables(
        d / "tables.md", {"evalB_levels_parity": t1, "evalB_levels_parity_per_voxel": t2}
    )
    (d / "summary.json").write_text(json.dumps(summ, indent=1) + "\n")
    print(t1)
    print(t2)


if __name__ == "__main__":
    main()
