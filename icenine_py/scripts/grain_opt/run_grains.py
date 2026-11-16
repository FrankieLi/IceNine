"""Work items (2) and (3): stage A (CMA on grain rotation + linear field) and stage B (Adam on free
per-voxel offsets with a smoothness penalty) on picked grains of ``grad``, from a start map.

Grains are the flood-fill groups of the START map (neighbour misorientation < 2 deg). Picks: per
gradient level, true grains at size quantiles 0.1 .. 0.95; the group fitted for a true grain is
the start-map group holding most of its voxels. Every voxel is scored against its own truth.

Acceptance (per stage): the update is kept only if the objective and the mean strict cost of the
group's voxels do not get worse; otherwise the previous orientations stay.

Writes ``benchmarks/grain_opt/runs/<name>/``: ``run.json`` (settings, per-grain records),
``voxels.json.gz`` (per-voxel rows).

Usage (from icenine_py/):
  uv run python scripts/grain_opt/run_grains.py --name A_sink --stages A --objective sinkhorn
  uv run python scripts/grain_opt/run_grains.py --name AB_sink --stages AB
"""

import argparse
import json
import time
from pathlib import Path
from typing import Any, Dict, List

import common as CM  # noqa: I001
import numpy as np
import torch

from icenine import grain_fit as GF


def pick_true_grains(T: Any, n: int) -> List[int]:
    gm = np.asarray(T.gmap)
    levels = sorted(set(T.level_of.values()))
    per = max(n // len(levels), 1)
    picks: List[int] = []
    for lv in levels:
        gids = [g for g in T.gids if T.level_of[g] == lv]
        sizes = {g: int((gm == g).sum()) for g in gids}
        gids = sorted(gids, key=lambda g: (sizes[g], g))
        for q in np.linspace(0.1, 0.95, per):
            picks.append(gids[int(q * (len(gids) - 1))])
    return list(dict.fromkeys(picks))[:n]


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--name", required=True)
    ap.add_argument("--stages", default="A", help="A, B (from the start) or AB")
    ap.add_argument("--objective", default="sinkhorn", choices=["sinkhorn", "ls"])
    ap.add_argument("--grains", type=int, default=20)
    ap.add_argument("--true-grains", type=int, nargs="*", default=None)
    ap.add_argument("--start", default="paint")
    ap.add_argument("--order", type=int, default=1)
    ap.add_argument("--n-sub", type=int, default=2)
    ap.add_argument("--sigma0", type=float, default=0.05)
    ap.add_argument("--max-evals", type=int, default=200)
    ap.add_argument("--b-steps", type=int, default=60)
    ap.add_argument("--b-lr", type=float, default=0.01)
    ap.add_argument("--b-mu", type=float, default=1.0)
    ap.add_argument("--reach", type=float, default=3.0)
    ap.add_argument("--blur", type=float, default=1.0)
    a = ap.parse_args()
    out = CM.OUT / "runs" / a.name
    S = CM.Sample(out, start=a.start)
    T = S.T
    gm = np.asarray(T.gmap)
    t0 = time.time()
    lab = GF.group_grains(S.R0, S.nb, S.sym, 2.0)
    t_group = time.time() - t0
    picks = a.true_grains or pick_true_grains(T, a.grains)
    spec = GF.ObjectiveSpec(kind=a.objective, reach_px=a.reach or None, blur_px=a.blur)
    C = {r["voxel"]: r for r in CM.G.load_gz(CM.C.OUT / "fine" / "runs" / "C" / "refine.json.gz")}
    R_cur = S.R0.copy()
    grains: List[Dict[str, Any]] = []
    vox_rows: List[Dict[str, Any]] = []
    for g in picks:
        tv = np.nonzero(gm == g)[0]
        labs, cnt = np.unique(lab[tv], return_counts=True)
        grp = int(labs[np.argmax(cnt)])
        members = [int(v) for v in np.nonzero(lab == grp)[0]]
        shell = GF.shell_of(members, S.nb)
        bnd = set(GF.boundary_members(members, S.nb))
        t1 = time.time()
        prob = GF.GrainProblem(
            S.model, S.images, R_cur, S.V, S.xy, members, shell, S.sym,
            n_sub=a.n_sub, order=a.order, spec=spec,
        )
        t_region = time.time() - t1
        R_start = prob.R_start
        strict0 = S.strict(R_start, members)
        f0 = prob.value(R_start)
        rec: Dict[str, Any] = dict(
            true_grain=int(g), level=T.level_of[g], group=grp, n_vox=len(members),
            n_true=int(len(tv)), n_in_true=int((gm[members] == g).sum()), n_shell=len(shell),
            n_boundary=len(bnd), n_spots=prob.n_spots, t_region_s=t_region, f_start=f0,
            strict_start=float(strict0.mean()),
        )
        R_now, f_now, s_now = R_start, f0, strict0
        stages: Dict[str, np.ndarray] = {"start": R_start}
        strict_of: Dict[str, np.ndarray] = {"start": strict0}
        if "A" in a.stages:
            t2 = time.time()
            th, R_A, info = prob.stage_a(sigma0_deg=a.sigma0, max_evals=a.max_evals)
            sA = S.strict(R_A, members)
            fA = prob.value(R_A)
            acc = fA <= f_now and sA.mean() <= s_now.mean()
            info.update(t_s=time.time() - t2, f=fA, strict=float(sA.mean()), accepted=bool(acc),
                        theta0_deg=np.degrees(prob.theta0).tolist(), theta_deg=np.degrees(th).tolist())
            rec["A"] = info
            stages["A_raw"], strict_of["A_raw"] = R_A, sA
            if acc:
                R_now, f_now, s_now = R_A, fA, sA
            stages["A"], strict_of["A"] = R_now, s_now
        if "B" in a.stages:
            t3 = time.time()
            R_B, info = prob.stage_b(R_now, steps=a.b_steps, lr_deg=a.b_lr, mu=a.b_mu,
                                     neighbours=S.nb)
            sB = S.strict(R_B, members)
            fB = prob.value(R_B)
            acc = fB <= f_now and sB.mean() <= s_now.mean()
            info.update(t_s=time.time() - t3, f=fB, strict=float(sB.mean()), accepted=bool(acc))
            rec["B"] = info
            stages["B_raw"], strict_of["B_raw"] = R_B, sB
            if acc:
                R_now, f_now, s_now = R_B, fB, sB
            stages["B"], strict_of["B"] = R_now, s_now
        errs = {k: T.err(v, members) for k, v in stages.items()}
        for k in errs:
            rec[f"err_{k}"] = CM.err_stats(errs[k])
        for i, v in enumerate(members):
            row = dict(voxel=v, true_grain=int(gm[v]), level=T.level_of[int(gm[v])],
                       group=grp, pick=int(g), boundary=v in bnd)
            for k in errs:
                row[f"err_{k}"] = float(errs[k][i])
                row[f"strict_{k}"] = float(strict_of[k][i])
            c = C.get(v)
            if c is not None:
                row["err_C"] = float(T.err(np.asarray(c["end"])[None], [v])[0])
                row["strict_C"] = float(c["end_cost"])
            vox_rows.append(row)
        R_cur[members] = R_now
        grains.append(rec)
        print(json.dumps({k: rec[k] for k in ("true_grain", "level", "n_vox", "n_spots")}
                         | {k: (rec[k]["median"], rec[k]["n_gt01"], rec[k]["n_wrong"])
                            for k in rec if k.startswith("err_")}
                         | {s: (rec[s]["f_start"], rec[s]["f_end"], rec[s]["accepted"], round(rec[s]["t_s"], 1))
                            for s in ("A", "B") if s in rec}), flush=True)
    run = dict(args=vars(a), t_group_s=t_group, grains=grains, timing=CM.C.preflight_label(out))
    (out / "run.json").write_text(json.dumps(run, indent=1, default=float) + "\n")
    CM.G.dump_gz(vox_rows, out / "voxels.json.gz")


if __name__ == "__main__":
    main()
