"""Task A: why does FindOptimal fail its hit-ratio-1 convergence test on some voxels without
CppEdgeRules? Hit ratio of the complete cost at the true orientation of the six voxels with the
default, with each of the three edge rules alone (by switching the helper functions) and with all.

  uv run python scripts/eval_reduction/hit_ratio_check.py --out benchmarks/eval_reduction/A/step2
"""

import argparse
import json
import sys
from pathlib import Path
from typing import Any, Dict

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import parity_common as PC  # noqa: E402


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    out = Path(a.out).resolve()
    ctx = PC.PyContext()
    from icenine import cost_functions as CF

    o_bins, o_pix = CF._omega_bins, CF._pixel_coordinate
    res: Dict[str, Any] = {}
    # each rule alone: the frame and pixel rules through the module helpers
    # (the vertex rule has no variant of its own here)
    variants = {
        "default": dict(edge=False, frame=False, pixel=False),
        "frame rule only": dict(edge=False, frame=True, pixel=False),
        "pixel -1 rule only": dict(edge=False, frame=False, pixel=True),
        "all three (CppEdgeRules)": dict(edge=True, frame=False, pixel=False),
    }
    for name, v in variants.items():
        cf = ctx.cost_fn("full", cpp_edge_rules=v["edge"])
        CF._omega_bins = (lambda om, lo, w, e: o_bins(om, lo, w, True)) if v["frame"] else o_bins
        CF._pixel_coordinate = (lambda d, h, s, e: o_pix(d, h, s, True)) if v["pixel"] else o_pix
        rows = {}
        for vox in PC.SIX_VOXELS:
            info = cf.evaluate(
                ctx.truth(vox).astype(np.float32), ctx.verts(vox), ctx.mic.voxels[vox].phase
            )
            rows[str(vox)] = dict(
                pixel_overlap=int(info.pixel_overlap),
                pixel_on=int(info.pixel_on_detector),
                hit_ratio=float(info.hit_ratio),
                converged=bool(info.pixel_overlap >= info.pixel_on_detector),
                peak_overlap=int(info.peak_overlap),
                peak_on=int(info.peak_on_detector),
            )
        res[name] = rows
    CF._omega_bins, CF._pixel_coordinate = o_bins, o_pix
    (out / "hit_ratio_at_truth.json").write_text(json.dumps(res, indent=1) + "\n")
    for name, rows in res.items():
        print(
            name,
            {
                k: (r["pixel_overlap"], r["pixel_on"], round(r["hit_ratio"], 4))
                for k, r in rows.items()
            },
        )


if __name__ == "__main__":
    main()
