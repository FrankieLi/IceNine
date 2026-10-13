"""Task A step 3: the Python omega interpolation (ScreenEvaluator.overlaps_interp) against the C++
interpolated screen on the same (FZ centre, local rotation) pairs of the fixed-orientation test
(step0: results.json.gz has the C++ INT arm per pair), and the interpolation error against the exact
omegas as a function of the local-rotation angle.

  uv run python scripts/eval_reduction/interp_vs_cpp.py --dir DIR \
      --step0 benchmarks/eval_reduction/A/step0 --out benchmarks/eval_reduction/A/step3
"""

import argparse
import gzip
import json
import sys
from pathlib import Path
from typing import Any, Dict, List

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import parity_common as PC  # noqa: E402

sys.path.insert(0, str(PC.ICE / "scripts" / "common"))
import doc_tables as DT  # noqa: E402
from stats import wilson  # noqa: E402


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--dir", required=True, help="in_<voxel>.json of the step-0 inputs")
    ap.add_argument("--step0", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    d, out = Path(a.dir).resolve(), Path(a.out).resolve()
    step0 = Path(a.step0).resolve()
    out.mkdir(parents=True, exist_ok=True)
    rows = json.load(gzip.open(step0 / "results.json.gz"))["rows"]
    ctx = PC.PyContext()
    from icenine.batched_screen import ScreenEvaluator

    cf_e = ctx.cost_fn("screen", cpp_edge_rules=True)
    ev = ScreenEvaluator(cf_e)
    res: List[Dict[str, Any]] = []
    for v in PC.PARITY_VOXELS:
        meta = {m["id"]: m for m in json.load(open(d / f"in_{v}.json"))}
        verts, phase = ctx.verts(v), ctx.mic.voxels[v].phase
        for r in [x for x in rows if x["voxel"] == v and x["kind"] == 1]:
            m = meta[r["id"]]
            c = np.asarray(m["centre"], dtype=np.float64)[None]
            lo = np.asarray(m["local"], dtype=np.float64)[None]
            po_i, pon_i = ev.overlaps_interp(c, lo, verts, phase)
            # exact omegas of the C++ float32 product (the C++ driver's own product)
            prod = (lo[0].astype(np.float32) @ c[0].astype(np.float32)).astype(np.float64)
            po_e, _ = ev.overlaps_exact(prod[None], verts, phase)
            res.append(
                dict(
                    voxel=v,
                    id=r["id"],
                    label=r["label"],
                    theta=float(m.get("theta", np.nan)),
                    cpp_int_pkov=r["INT"]["pkov"],
                    cpp_int_pkon=r["INT"]["pkon"],
                    py_int_pkov=int(po_i[0, 0]),
                    py_int_pkon=int(pon_i[0, 0]),
                    cpp_exact_pkov=r["EXS"]["cpp"]["pkov"],
                    py_exact_pkov=int(po_e[0]),
                )
            )
    (out / "interp_vs_cpp.json").write_text(json.dumps(res))
    T: Dict[str, str] = {}
    tab = []
    for g in ["pair_random", "sweep", "pair_differ_pass", "pair_differ_count", "ALL"]:
        rr = [x for x in res if g == "ALL" or x["label"] == g]
        n = len(rr)
        same_cnt = sum(
            (x["cpp_int_pkov"], x["cpp_int_pkon"]) == (x["py_int_pkov"], x["py_int_pkon"])
            for x in rr
        )
        same_pass = sum((x["cpp_int_pkov"] > 0) == (x["py_int_pkov"] > 0) for x in rr)
        py_exact_same_pass = sum((x["cpp_int_pkov"] > 0) == (x["py_exact_pkov"] > 0) for x in rr)
        lo, hi = wilson(same_pass, n)
        tab.append(
            dict(
                group=g,
                n=n,
                py_interp_pass_matches_cpp=f"{same_pass} [{lo:.4f}, {hi:.4f}]",
                py_interp_counts_match_cpp=same_cnt,
                py_exact_pass_matches_cpp_interp=py_exact_same_pass,
            )
        )
    T["interp_vs_cpp"] = DT.markdown_table(
        tab,
        [
            "group",
            "n",
            "py_interp_pass_matches_cpp",
            "py_interp_counts_match_cpp",
            "py_exact_pass_matches_cpp_interp",
        ],
    )
    DT.write_tables(out / "tables_interp_vs_cpp.md", T)
    print(T["interp_vs_cpp"])


if __name__ == "__main__":
    main()
