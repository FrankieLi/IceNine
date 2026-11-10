"""AR1 timing: the indexer's time per voxel on the pilot's full-search voxels, and the one-time cost
of the observed rays (pairing the spots of all frames), single process.

Gated by the timing preflight: if the machine is not quiet the preflight is still saved and the
run is labelled "contended". Writes benchmarks/autoresearch/ar1/timing_pilot.json and
preflight_timing.json.

  uv run python scripts/autoresearch/ar1/timing.py
"""

import json
import sys
import time
from pathlib import Path
from typing import Any, Dict, List

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "common"))
import ar1_common as AC  # noqa: E402
import preflight  # noqa: E402

from icenine import ray_indexing as RI  # noqa: E402


def main() -> None:
    AC.OUT.mkdir(parents=True, exist_ok=True)
    pf = preflight.preflight()
    try:
        preflight.require_quiet()
        label = "single-worker"
    except preflight.MachineBusyError as exc:
        label = "contended"
        pf["require_quiet_refused"] = str(exc)[:300]
    (AC.OUT / "preflight_timing.json").write_text(json.dumps(pf, indent=1, default=str) + "\n")
    stored = AC.pilot_stored()
    AC.init("full")
    p = RI.IndexParams()
    voxels = sorted(stored)
    # one-time: observed rays of a fresh SmoothCost
    e0 = AC.env("full", voxels[0])
    from icenine.smooth_cost import SmoothCost

    b = AC.bench()._E["env"]["base"]
    sc = SmoothCost(
        b.simulator,
        b.detector_list,
        b.range_map,
        b.exp_data,
        b.sample,
        b.structure_list,
        eta_limit=b.eta_limit,
    )
    t = time.perf_counter()
    rays = sc.observed_rays()
    t_rays = time.perf_counter() - t
    crystal = RI.crystal_from_smooth_cost(sc, AC.symmetry_of(e0), e0["phase"], p.q_max)
    rows: List[Dict[str, float]] = []
    for v in voxels:
        e = AC.env("full", v)
        t0 = time.perf_counter()
        att = sc.attached_rays(e["verts"])
        t1 = time.perf_counter()
        g = RI.scattering_vectors(sc, att)
        res = RI.index_rays(g, crystal, p, 20)
        t2 = time.perf_counter()
        rows.append(dict(voxel=v, attach_s=t1 - t0, index_s=t2 - t1, n_att=len(att), **res.seconds))
    a = lambda k: np.array([r[k] for r in rows])  # noqa: E731
    out: Dict[str, Any] = dict(
        label=label,
        n_voxels=len(rows),
        n_rays_layer=int(len(rays.frame)),
        observed_rays_build_s=t_rays,
        attach_s_median=float(np.median(a("attach_s"))),
        index_s_median=float(np.median(a("index_s"))),
        index_s_p90=float(np.percentile(a("index_s"), 90)),
        index_s_max=float(a("index_s").max()),
        vote_s_median=float(np.median(a("vote"))),
        cluster_s_median=float(np.median(a("cluster"))),
        refine_s_median=float(np.median(a("refine"))),
        per_voxel_total_s_median=float(np.median(a("attach_s") + a("index_s"))),
        rows=rows,
    )
    (AC.OUT / "timing_pilot.json").write_text(json.dumps(out, indent=1) + "\n")
    print({k: v for k, v in out.items() if k != "rows"})


if __name__ == "__main__":
    main()
