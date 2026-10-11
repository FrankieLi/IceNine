#!/usr/bin/env python3
"""Step 4b: quality / confidence of C++ (debug build, patch 4 `CEVAL`) and Python at the true
orientation and at 3 nearby orientations of ~10 E0 cases (frame conversion and setup check).

Nearby orientations: the truth rotated by 0.1, 0.3 and 1.0 deg about fixed random axes (seed 0).
Python: the E0 local cost function (radius 0) with the images attached, for min_sin_eta 0.3 (E0 as
run) and 0 (C++ has no such filter), and the C++ overlap rule off and on.

  E0_CPP_BIN=<scratch>/IceNine E0_CPP_SCRATCH=<dir> uv run python scripts/phase_d/e0_cost_check.py
"""

import json
import os
import subprocess
import sys
from pathlib import Path
from typing import Any, Dict, List

import numpy as np

HERE = Path(__file__).resolve().parent
ICE = HERE.parents[1]
sys.path.insert(0, str(ICE / "scripts" / "findoptimal_robustness"))
sys.path.insert(0, str(ICE / "scripts" / "common"))
sys.path.insert(0, str(HERE))
import common as C  # noqa: E402
import e0_cpp as E  # noqa: E402

OUT = ICE / "benchmarks" / "phase_d_seed_diag" / "level0" / "e0"
ANGLES = [0.0, 0.1, 0.3, 1.0]


def pick() -> List[str]:
    cases = json.loads((OUT / "cpp_e0_summary.json").read_text())["cases"]
    wrong = lambda e: e > C.WRONG_DEG  # noqa: E731
    groups = {
        ("clean", True, True): 2, ("clean", False, True): 2, ("clean", True, False): 1,
        ("all", True, True): 2, ("all", False, True): 2, ("all", True, False): 1,
    }  # fmt: skip  # (variant, cpp_right, py_on_wrong): count
    out: List[str] = []
    for (var, cpp_right, on_wrong), n in groups.items():
        sel = [k for k in sorted(cases) if k.endswith("_" + var)
               and (not wrong(cases[k]["cpp_err"])) == cpp_right
               and wrong(cases[k]["py_on_err"]) == on_wrong][:n]  # fmt: skip
        out += sel
    return out


def rot(axis: np.ndarray, deg: float) -> np.ndarray:
    a = axis / np.linalg.norm(axis)
    K = np.array([[0, -a[2], a[1]], [a[2], 0, -a[0]], [-a[1], a[0], 0]])
    t = np.radians(deg)
    return np.eye(3) + np.sin(t) * K + (1 - np.cos(t)) * (K @ K)


def main() -> None:
    names = pick()
    rng = np.random.default_rng(0)
    axes = rng.normal(size=(len(names), 3))
    C.init_worker(C.worker_args())
    W = C.get_worker()
    from optimizer_sweep import voxel_context

    rows: List[Dict[str, Any]] = []
    for ci, name in enumerate(names):
        v, var = name[1:].split("_", 1)
        vidx = int(v)
        keys = np.load(C.CACHE_DIR / "images" / f"v{vidx}_{var}.npz")["keys"]
        C.attach(keys)
        vctx = voxel_context(W.ctx, vidx)
        Rt = np.asarray(vctx.R_true, dtype=np.float64)
        Rs = [rot(axes[ci], d) @ Rt for d in ANGLES]
        # C++
        root = E.scratch() / f"chk_{name}"
        cfg = E.write_case(vidx, var, root)
        mf = root / "orient.txt"
        mf.write_text("\n".join(" ".join(f"{x:.12g}" for x in R.reshape(-1)) for R in Rs) + "\n")
        p = subprocess.run(
            [os.environ["E0_CPP_BIN"], "r", str(cfg)], cwd=E.EX,
            env=dict(os.environ, CEVAL=str(mf)), capture_output=True, text=True,
        )  # fmt: skip
        cpp = [ln.split() for ln in p.stdout.splitlines() if ln.startswith("CE ")]
        import shutil

        shutil.rmtree(root, ignore_errors=True)
        for msin in (0.3, 0.0):
            for rule in (False, True):
                W.rec.min_sin_eta = msin
                W.rec.setup.search_params.cpp_overlap_rule = rule
                fn = W.rec._make_local_cost_fn()
                for k, (ang, R) in enumerate(zip(ANGLES, Rs)):
                    info = fn.evaluate(R, vctx.vertices, vctx.voxel.phase)
                    row = dict(case=name, angle_deg=ang, py_min_sin_eta=msin, py_rule=rule,
                               py_quality=float(info.quality), py_confidence=float(info.confidence),
                               py_hit_ratio=float(info.hit_ratio),
                               py_peak_overlap=int(info.peak_overlap),
                               py_peak_on_det=int(info.peak_on_detector),
                               py_pixel_overlap=int(info.pixel_overlap),
                               py_pixel_on_det=int(info.pixel_on_detector))  # fmt: skip
                    if k < len(cpp):
                        t = cpp[k]
                        row.update(cpp_quality=float(t[2]), cpp_confidence=float(t[3]),
                                   cpp_hit_ratio=float(t[4]), cpp_peak_overlap=int(t[5]),
                                   cpp_peak_on_det=int(t[6]), cpp_pixel_overlap=int(t[7]),
                                   cpp_pixel_on_det=int(t[8]))  # fmt: skip
                    rows.append(row)
        W.rec.setup.search_params.cpp_overlap_rule = False
        W.rec.min_sin_eta = 0.3
        print(name, "done", flush=True)
    (OUT / "cost_check.json").write_text(json.dumps(rows, indent=1))


if __name__ == "__main__":
    main()
