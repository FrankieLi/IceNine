"""Trace of VarianceMinimizing (variance_minimizing_optimize) started at the truth, on the full
images vs the isolated one-voxel render of the same voxel: per subregion run the radius, the steps,
the cost variance and the cost mean; capped at --max-total steps. Diagnoses why the budget
extension (variance > 0.02^2 adds the subregion steps to the budget) does not terminate.

Usage (from icenine_py/): uv run python scripts/phase_d/seed_diag_variance.py --voxel 15901
"""

import argparse
import json
import math
import os
import sys
from pathlib import Path
from typing import Any, Dict, List

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import seed_diag as SD  # noqa: E402


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--voxel", type=int, required=True)
    ap.add_argument("--max-total", type=int, default=20000)
    ap.add_argument("--boxes", type=float, nargs="+", default=[0.3333, 0.22, 0.2])
    a = ap.parse_args()
    os.chdir(SD.EX)
    from scipy.spatial.transform import Rotation
    from icenine import orientation_search as OS
    from icenine import reconstructor as RC
    from icenine.config_file import ConfigFile
    from icenine.experimental_data import ExperimentalData
    from icenine.mic_file import MicFile

    cfg = ConfigFile.from_file(str(SD.HERE / "configs" / "ReconstructPhaseD_mc_clean.config"))
    truth = MicFile.read(SD.MIC)
    R_true = np.asarray(truth.voxels[a.voxel].orientation, dtype=np.float64)
    out: Dict[str, Any] = {"voxel": a.voxel, "max_total": a.max_total}
    for images in ("full", "isolated"):
        stack = (
            SD.EX / "ScatteringData_PhaseD" / "stacks" / "clean.npy"
            if images == "full"
            else SD.build_isolated_stack(a.voxel)
        )
        setup = RC.setup_reconstruction(
            cfg, exp_data=ExperimentalData.from_binary_memmap(stack, 180, 2)
        )
        rec = RC.AdaptiveVoxelReconstructor(setup)
        mic = setup.sample.get_mic()
        v = mic.voxels[a.voxel]
        verts = RC._get_voxel_vertices(v)
        local = rec._make_local_cost_fn()
        # cost vs random misorientation of fixed magnitude (the landscape the variance sees)
        prof = {}
        prng = np.random.default_rng(5)
        for mag in (0.05, 0.1, 0.2, 0.33, 0.5, 1.0):
            cs = []
            for _ in range(200):
                ax = prng.normal(size=3)
                ax /= np.linalg.norm(ax)
                Rm = Rotation.from_rotvec(np.radians(mag) * ax).as_matrix() @ R_true
                cs.append(local.evaluate(Rm.astype(np.float32), verts, v.phase).cost)
            prof[str(mag)] = dict(
                mean=float(np.mean(cs)),
                var=float(np.var(cs)),
                q10=float(np.quantile(cs, 0.1)),
                q90=float(np.quantile(cs, 0.9)),
            )
        out.setdefault("profile", {})[images] = prof
        print(
            images,
            "profile (cost at fixed misorientation):",
            {k: (round(x["mean"], 3), round(x["var"], 5)) for k, x in prof.items()},
            flush=True,
        )
        for box_deg in a.boxes:
            log: List[Dict[str, float]] = []
            state = {"total": 0}
            mc, _ = rec._make_optimizers(local, verts, v.phase, np.random.default_rng(3))
            orig = mc._zero_temp_with_variance

            def traced(init, radius, n_steps, _o=orig, _l=log, _s=state):  # type: ignore
                r = _o(init, radius, n_steps)
                _s["total"] += n_steps
                _l.append(
                    dict(
                        radius_deg=math.degrees(radius),
                        steps=n_steps,
                        var=float(r[2]),
                        best=float(r[1]),
                    )
                )
                if _s["total"] >= a.max_total:
                    raise StopIteration
                return r

            mc._zero_temp_with_variance = traced  # type: ignore[method-assign]
            try:
                mc.variance_minimizing_optimize(
                    R_true.astype(np.float32), math.radians(box_deg), 200, 2, 0.0, 0.02**2
                )
                ended = "terminated"
            except StopIteration:
                ended = "cap"
            vs = np.array([l["var"] for l in log])
            out.setdefault(images, {})[str(box_deg)] = dict(
                ended=ended,
                n_runs=len(log),
                total_steps=state["total"],
                frac_runs_var_below_thr=float((vs <= 0.02**2).mean()),
                var_quantiles=[float(x) for x in np.quantile(vs, [0.05, 0.25, 0.5, 0.75, 0.95])],
            )
            print(images, box_deg, out[images][str(box_deg)], flush=True)
    p = SD.OUT / f"variance_trace_v{a.voxel}.json"
    p.write_text(json.dumps(out, indent=1) + "\n")


if __name__ == "__main__":
    main()
