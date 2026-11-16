"""GO2 onwards: run a pipeline of fitting stages on benchmark grains and record per-voxel results.

A pipeline is a comma list of stages, applied in order to each grain (largest first in the pick
list order), each stage gated by the acceptance rule (objective of the stage and mean strict cost
of the group's voxels do not get worse):

* ``A<obj>/<evals>[/o2]``: CMA over the grain field (rotation + linear field; ``o2`` quadratic)
  on objective ``obj`` in {sink, ls, strict}, at most ``evals`` evaluations.
* ``B<obj>/<steps>``: Adam on free per-voxel offsets with smoothness (``obj`` sink or ls).
* ``V/<evals>``: per-voxel CMA on that voxel's strict cost (3 numbers, sigma 0.02 deg), at most
  ``evals`` strict evaluations per voxel (AR3 C used 61 from the painted start).

Cost accounting: ``renders_per_voxel`` = grain-objective evaluations of every A/B stage (each
renders every voxel of the group once) + per-voxel strict evaluations of V stages. Wall seconds
per stage are recorded too.

Writes ``benchmarks/grain_opt/runs/<name>/``: ``run.json`` and ``voxels.json.gz``.

Usage (from icenine_py/):
  uv run python scripts/grain_opt/run_pipeline.py --name s_A100 --pipeline Astrict/100
  uv run python scripts/grain_opt/run_pipeline.py --name hyb --pipeline Asink/150,Astrict/100,V/20
"""

import argparse
import json
import math
import time
from typing import Any, Dict, List, Tuple

import common as CM  # noqa: I001
import numpy as np
import torch

from icenine import grain_fit as GF
from run_grains import pick_true_grains

OBJ = {"sink": "sinkhorn", "ls": "ls", "strict": "strict"}


def parse(p: str) -> List[Tuple[str, str, int, int]]:
    out = []
    for st in p.split(","):
        parts = st.split("/")
        head, n = parts[0], int(parts[1])
        order = 2 if len(parts) > 2 and parts[2] == "o2" else 1
        kind, obj = head[0], head[1:]
        assert kind in "ABV", st
        if kind != "V":
            assert obj in OBJ, st
        out.append((kind, obj, n, order))
    return out


def voxel_cma(S: Any, R0: np.ndarray, v: int, evals: int, sigma_deg: float = 0.02) -> Tuple[np.ndarray, int]:
    import cma

    R0 = np.asarray(R0, float)

    def f(w: np.ndarray) -> float:
        R = GF.rotvec_to_matrix(torch.as_tensor(np.asarray(w, float))).numpy() @ R0
        return float(S.strict(R[None], [v])[0])

    f0 = f(np.zeros(3))
    es = cma.CMAEvolutionStrategy([0.0, 0.0, 0.0], math.radians(sigma_deg),
                                  {"seed": 1, "maxfevals": evals, "verbose": -9})
    best, fb, n = np.zeros(3), f0, 1
    while not es.stop():
        X = es.ask()
        F = [f(x) for x in X]
        n += len(X)
        es.tell(X, F)
        i = int(np.argmin(F))
        if F[i] < fb:
            fb, best = F[i], np.asarray(X[i], float)
    return GF.rotvec_to_matrix(torch.as_tensor(best)).numpy() @ R0, n


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--name", required=True)
    ap.add_argument("--pipeline", required=True)
    ap.add_argument("--sample", default="grad")
    ap.add_argument("--start", default="paint")
    ap.add_argument("--grains", type=int, default=20)
    ap.add_argument("--true-grains", type=int, nargs="*", default=None)
    ap.add_argument("--group", default="flood", choices=["flood", "true"],
                    help="grain membership: flood fill of the start map (2 deg) or true grain ids")
    ap.add_argument("--sigma0", type=float, default=0.05)
    ap.add_argument("--b-lr", type=float, default=0.01)
    ap.add_argument("--b-mu", type=float, default=1.0)
    a = ap.parse_args()
    stages = parse(a.pipeline)
    out = CM.OUT / "runs" / a.name
    S = CM.Sample(out, start=a.start, sample=a.sample)
    T = S.T
    gm = np.asarray(T.gmap)
    lab = GF.group_grains(S.R0, S.nb, S.sym, 2.0) if a.group == "flood" else gm.copy()
    picks = a.true_grains or pick_true_grains(T, a.grains)
    Cref: Dict[int, Any] = {}
    if a.sample == "grad" and a.start == "paint":
        Cref = {r["voxel"]: r for r in CM.G.load_gz(CM.C.OUT / "fine" / "runs" / "C" / "refine.json.gz")}
    R_cur = S.R0.copy()
    grains, vox_rows = [], []
    t_all = time.time()
    for g in picks:
        tv = np.nonzero(gm == g)[0]
        labs, cnt = np.unique(lab[tv], return_counts=True)
        grp = int(labs[np.argmax(cnt)])
        members = [int(v) for v in np.nonzero(lab == grp)[0]]
        shell = GF.shell_of(members, S.nb)
        bnd = set(GF.boundary_members(members, S.nb))
        t1 = time.time()
        prob = GF.GrainProblem(S.model, S.images, R_cur, S.V, S.xy, members, shell, S.sym,
                               strict_fn=S.strict)
        rec: Dict[str, Any] = dict(true_grain=int(g), level=T.level_of[g], group=grp,
                                   n_vox=len(members), n_in_true=int((gm[members] == g).sum()),
                                   n_spots=prob.n_spots, t_region_s=time.time() - t1, stages=[])
        R_now = prob.R_start
        s_now = S.strict(R_now, members)
        renders = 0
        snaps = {"start": (R_now, s_now)}
        for si, (kind, obj, n, order) in enumerate(stages):
            t2 = time.time()
            info: Dict[str, Any] = dict(stage=f"{kind}{obj}/{n}")
            if kind == "A":
                spec = GF.ObjectiveSpec(kind=OBJ[obj])
                prob.spec = spec
                if order != 1 or prob.field.feats.shape[1] != 2 * order + (order == 2):
                    prob.field = GF.FieldA.build(prob.field.R_ref, S.xy[members], order)
                th0 = prob.field.fit_theta(R_now, S.sym)
                e0 = prob.n_evals
                th, R_new, inf = prob.stage_a(sigma0_deg=a.sigma0, max_evals=n, theta0=th0)
                f_old, f_new = prob.value(R_now), prob.value(R_new)
                used = prob.n_evals - e0
                info.update(inf, f_old=f_old, f_new=f_new)
            elif kind == "B":
                prob.spec = GF.ObjectiveSpec(kind=OBJ[obj])
                e0 = prob.n_evals
                R_new, inf = prob.stage_b(R_now, steps=n, lr_deg=a.b_lr, mu=a.b_mu, neighbours=S.nb)
                f_old, f_new = prob.value(R_now), prob.value(R_new)
                used = prob.n_evals - e0
                info.update(inf, f_old=f_old, f_new=f_new)
            else:  # per-voxel strict CMA
                R_new = R_now.copy()
                used_v = []
                for i, v in enumerate(members):
                    R_new[i], k = voxel_cma(S, R_now[i], v, n)
                    used_v.append(k)
                used = float(np.mean(used_v))
                f_old = f_new = float("nan")
            s_new = S.strict(R_new, members)
            if kind == "V":  # per voxel: keep a voxel's fit only if its strict cost did not rise
                keep = s_new <= s_now
                R_new = np.where(keep[:, None, None], R_new, R_now)
                s_new = np.where(keep, s_new, s_now)
                acc = True
                info["voxels_kept"] = int(keep.sum())
            else:
                acc = (f_new <= f_old) and (s_new.mean() <= s_now.mean())
            info.update(accepted=bool(acc), strict_old=float(s_now.mean()),
                        strict_new=float(s_new.mean()), renders_per_voxel=float(used),
                        t_s=time.time() - t2)
            renders += used
            if acc:
                R_now, s_now = R_new, s_new
            snaps[f"s{si}"] = (R_now, s_now)
            rec["stages"].append(info)
        rec["renders_per_voxel"] = renders
        rec["strict_n"] = prob.n_strict
        errs = {k: T.err(r, members) for k, (r, _s) in snaps.items()}
        for k in errs:
            rec[f"err_{k}"] = CM.err_stats(errs[k])
        for i, v in enumerate(members):
            row = dict(voxel=v, true_grain=int(gm[v]), level=T.level_of[int(gm[v])], pick=int(g),
                       boundary=v in bnd)
            for k, (r, sc) in snaps.items():
                row[f"err_{k}"] = float(errs[k][i])
                row[f"strict_{k}"] = float(sc[i])
            c = Cref.get(v)
            if c is not None:
                row["err_C"] = float(T.err(np.asarray(c["end"])[None], [v])[0])
                row["strict_C"] = float(c["end_cost"])
            vox_rows.append(row)
        R_cur[members] = R_now
        grains.append(rec)
        last = f"s{len(stages) - 1}"
        print(json.dumps(dict(g=int(g), lv=rec["level"], n=len(members),
                              start=round(rec["err_start"]["median"], 4),
                              end=round(rec[f"err_{last}"]["median"], 4),
                              gt01=(rec["err_start"]["n_gt01"], rec[f"err_{last}"]["n_gt01"]),
                              acc=[s["accepted"] for s in rec["stages"]],
                              r=round(renders, 1), t=round(sum(s["t_s"] for s in rec["stages"]), 1))),
              flush=True)
    run = dict(args=vars(a), stages=a.pipeline, grains=grains, t_total_s=time.time() - t_all,
               timing=CM.C.preflight_label(out))
    (out / "run.json").write_text(json.dumps(run, indent=1, default=float) + "\n")
    CM.G.dump_gz(vox_rows, out / "voxels.json.gz")


if __name__ == "__main__":
    main()
