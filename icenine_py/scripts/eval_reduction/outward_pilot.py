"""Task F pilots: the settings of F4 and F5, chosen on four E0 voxels the benchmark does not use.

  schedule  F4 with several robust-scale schedules (tau, deg), 4 voxels x 4 offsets x 8 directions
  window    F4 with frame windows of 1 to 12 frames (how far a peak looks for a ray), 4 voxels
            x offsets 2, 3, 5 deg x 20 directions
  spacing   F5 with several shell spacings and acceptance costs, 4 voxels x offsets 2 and 5 deg x
            4 directions

  uv run python scripts/eval_reduction/outward_pilot.py schedule
  uv run python scripts/eval_reduction/outward_pilot.py spacing
Writes benchmarks/eval_reduction/F/pilot/{f4_schedule,f4_window,f5_spacing}.json and .txt.
"""

import json
import sys
import time
from pathlib import Path
from typing import Any, Dict, List

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import outward_bench as B  # noqa: E402

OUT = B.OUT / "pilot"
PILOT_VOXELS = json.loads((B.OUT / "voxel_sets.json").read_text())["pilot"]
SCHEDULES = {
    "0.25": (0.25,),
    "1, 0.5, 0.25": (1.0, 0.5, 0.25),
    "2, 1, 0.5, 0.25": (2.0, 1.0, 0.5, 0.25),
    "0.5, 0.25": (0.5, 0.25),
    "4, 2, 1, 0.5, 0.25": (4.0, 2.0, 1.0, 0.5, 0.25),
    "1, 0.25": (1.0, 0.25),
}


def schedule() -> None:
    from icenine.outward_search import angle_deg

    B.init_worker("e0")
    res: Dict[str, Dict[float, List[List[float]]]] = {k: {} for k in SCHEDULES}
    for v in PILOT_VOXELS:
        env = B.voxel_env(v)
        for off in (0.5, 2.0, 3.0, 5.0):
            starts = B.starts_for(env, off, 8)
            with B.Tracer(env["lf"], env["R_true"], env["verts"], env["phase"]) as tr:
                for k, taus in SCHEDULES.items():
                    gn = B.make_gn(env, tr, taus)
                    for d in range(8):
                        r = gn.solve(starts[d])
                        res[k].setdefault(off, []).append(
                            [angle_deg(r.R, env["R_true"]), r.n_jac, r.n_val]
                        )
    out: Dict[str, Any] = {}
    lines = [
        "F4 robust-scale schedules (deg): share within 1 deg / within 0.1 deg, "
        "median E4 Jacobian evaluations"
    ]
    for k, per in res.items():
        out[k] = {}
        for off, v in per.items():
            a = np.array(v)
            out[k][f"{off:g}"] = dict(
                n=len(a),
                within_1=float((a[:, 0] < 1).mean()),
                within_01=float((a[:, 0] < 0.1).mean()),
                median_jac=float(np.median(a[:, 1])),
            )
        lines.append(
            f"tau {k}: "
            + "; ".join(
                f"{o} deg {d['within_1']:.2f} / {d['within_01']:.2f} / {d['median_jac']:.0f}"
                for o, d in out[k].items()
            )
        )
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "f4_schedule.json").write_text(json.dumps(out, indent=1) + "\n")
    (OUT / "f4_schedule.txt").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


def window() -> None:
    from icenine.outward_search import RayGaussNewton, angle_deg

    B.init_worker("e0")
    wins = (1, 2, 3, 5, 8, 12)
    res: Dict[int, Dict[float, List[List[float]]]] = {w: {} for w in wins}
    for v in PILOT_VOXELS:
        env = B.voxel_env(v)
        for off in (2.0, 3.0, 5.0):
            starts = B.starts_for(env, off, 20)
            with B.Tracer(env["lf"], env["R_true"], env["verts"], env["phase"]) as tr:
                sc = B.lazy(env, "sc")
                for w in wins:
                    gn = RayGaussNewton(sc, env["verts"], tr.strict, env["phase"],
                                        tau_schedule_deg=B.TAU_SCHEDULE, accept_cost=B.ACCEPT,
                                        frame_window=w)  # fmt: skip
                    for R0 in starts:
                        r = gn.solve(R0)
                        res[w].setdefault(off, []).append(
                            [angle_deg(r.R, env["R_true"]), r.n_jac, r.n_val, float(r.accepted)]
                        )
    out: Dict[str, Any] = {}
    lines = [
        "F4 frame window (frames): share within 1 deg / within 0.1 deg / accepted, "
        "median Jacobian evaluations"
    ]
    for w, per in res.items():
        out[str(w)] = {}
        for off, v in per.items():
            a = np.array(v)
            out[str(w)][f"{off:g}"] = dict(n=len(a), within_1=float((a[:, 0] < 1).mean()),
                                           within_01=float((a[:, 0] < 0.1).mean()),
                                           accepted=float(a[:, 3].mean()),
                                           median_jac=float(np.median(a[:, 1])))  # fmt: skip
        lines.append(
            f"window {w}: "
            + "; ".join(
                f"{o} deg {d['within_1']:.2f} / {d['within_01']:.2f} / {d['accepted']:.2f} / "
                f"{d['median_jac']:.0f}"
                for o, d in out[str(w)].items()
            )
        )
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "f4_window.json").write_text(json.dumps(out, indent=1) + "\n")
    (OUT / "f4_window.txt").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


def spacing() -> None:
    from icenine.outward_search import (
        ShellSearch,
        angle_deg,
        screen_fn_from,
        smooth_fn_from,
    )

    B.init_worker("e0")
    rows: Dict[str, List[List[float]]] = {}
    for v in PILOT_VOXELS:
        env = B.voxel_env(v)
        for off in (2.0, 5.0):
            starts = B.starts_for(env, off, 4)
            for sp in (0.15, 0.3, 0.45, 0.6, 0.9):
                for acc in (0.5, 0.8):
                    with B.Tracer(env["lf"], env["R_true"], env["verts"], env["phase"]) as tr:
                        scr = screen_fn_from(B.lazy(env, "screen_cf"), env["verts"], env["phase"])
                        smo = smooth_fn_from(B.lazy(env, "sc"), env["verts"], env["phase"])
                        ss = ShellSearch(scr, tr.strict, smo, spacing_deg=sp, r_max_deg=B.R_MAX_DEG,
                                         accept_cost=acc)  # fmt: skip
                        for d in range(4):
                            t0 = time.perf_counter()
                            r = ss.search(starts[d])
                            rows.setdefault(f"{off:g}|{sp:g}|{acc:g}", []).append(
                                [angle_deg(r.R, env["R_true"]), r.n_evals, float(r.found),
                                 time.perf_counter() - t0]
                            )  # fmt: skip
    out: Dict[str, Any] = {}
    lines = [
        "F5 spacing and acceptance: offset | spacing | accept: share within 1 deg, share found, "
        "median final error (deg), median evaluations, median seconds"
    ]
    for k, v in sorted(rows.items(), key=lambda kv: tuple(float(x) for x in kv[0].split("|"))):
        a = np.array(v)
        out[k] = dict(n=len(a), within_1=float((a[:, 0] < 1).mean()), found=float(a[:, 2].mean()),
                      median_err=float(np.median(a[:, 0])), median_evals=float(np.median(a[:, 1])),
                      median_secs=float(np.median(a[:, 3])))  # fmt: skip
        d = out[k]
        lines.append(
            f"{k}: {d['within_1']:.2f}, {d['found']:.2f}, {d['median_err']:.3f}, "
            f"{d['median_evals']:.0f}, {d['median_secs']:.2f}"
        )
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "f5_spacing.json").write_text(json.dumps(out, indent=1) + "\n")
    (OUT / "f5_spacing.txt").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    {"schedule": schedule, "spacing": spacing, "window": window}[sys.argv[1]]()
