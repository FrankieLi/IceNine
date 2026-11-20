"""GO7f: warm-evaluate time of the strict cost with the rocking window on and off.

The rocking window ORs several binary frames per peak; ``FrameTable`` caches those unions as
lit-pixel indices. This script measures repeated evaluates (single-worker, preflight-gated) and the
distinct union keys one evaluate needs. Usage: ``uv run python scripts/grain_opt/go7f_union_timing.py
<out_dir>``; writes ``preflight.json`` and ``timing.json`` there.
"""

import json
import math
import sys
import time
from pathlib import Path
from typing import Any, Dict

import numpy as np

HERE = Path(__file__).resolve()
sys.path.insert(0, str(HERE.parents[2] / "tests"))
sys.path.insert(0, str(HERE.parents[1] / "common"))

import preflight  # noqa: E402
from test_center_screen import three_voxels  # noqa: E402

from icenine.cost_functions import VoxelCostFunction  # noqa: E402


def main(out: Path) -> None:
    out.mkdir(parents=True, exist_ok=True)
    info = preflight.require_quiet()
    (out / "preflight.json").write_text(json.dumps(info, indent=1, default=str))
    c = three_voxels.__wrapped__()
    phase, verts = c["voxel"].phase, c["vertices"]
    R = np.asarray(c["voxel"].orientation, dtype=np.float64)
    width = math.degrees(float(c["kwargs"]["range_map"].width))
    res: Dict[str, Any] = {"label": "single-worker", "rows": []}
    for frames in (0.0, 0.5, 1.0, 2.0):
        cf = VoxelCostFunction(
            pixel_radius=0, mode="hard", rocking_deg=frames * width, **c["kwargs"]
        )
        t0 = time.perf_counter()
        cf.evaluate(R, verts, phase)
        first = time.perf_counter() - t0
        times = []
        for _ in range(10):
            t0 = time.perf_counter()
            cf.evaluate(R, verts, phase)
            times.append(time.perf_counter() - t0)
        table = cf.frame_table
        keys = len(table.unions) if table is not None else 0
        res["rows"].append(
            {
                "rocking_frames": frames,
                "first_evaluate_s": first,
                "warm_evaluate_median_s": float(np.median(times)),
                "warm_evaluate_max_s": float(max(times)),
                "union_keys_cached": keys,
            }
        )
        print(res["rows"][-1], flush=True)
    (out / "timing.json").write_text(json.dumps(res, indent=1))


if __name__ == "__main__":
    main(Path(sys.argv[1]))
