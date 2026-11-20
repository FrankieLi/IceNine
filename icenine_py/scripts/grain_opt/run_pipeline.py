"""GO2 onwards: run a pipeline of fitting stages on benchmark grains and record per-voxel results.

A pipeline is a comma list of stages, applied in order to each grain (largest first in the pick
list order), each stage gated by the acceptance rule (objective of the stage and mean strict cost
of the group's voxels do not get worse):

* ``A<obj>/<evals>[/o2]``: CMA over the grain field (rotation + linear field; ``o2`` quadratic)
  on objective ``obj`` in {sink, sinki (intensity-weighted Sinkhorn), ls, strict}, at most
  ``evals`` evaluations.
* ``B<obj>/<steps>``: Adam on free per-voxel offsets with smoothness (``obj`` sink or ls).
* ``N/<evals>[/<sigma>]``: neighbour continuation (two breadth-first sweeps from the best 10%;
  candidates = own + fitted in-grain neighbours' orientations; then per-voxel CMA polish).
* ``L/<evals>[/3]``: per-voxel CMA step ladder 1.0 / 0.3 / 0.05 deg (or 3 / 1 / 0.3 / 0.05),
  ``evals`` split evenly over the legs.
* ``G/<n>[/<threshold>]``: gated global search: voxels whose strict cost exceeds the threshold
  (default n/100) run the reconstructor's ray-indexing seed path; the seed is kept if its strict
  cost is lower. Its evaluations count per voxel of the group.
* ``V/<evals>[/<sigma>]``: per-voxel CMA on that voxel's strict cost (3 numbers, sigma 0.02 deg
  unless given), at most
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

OBJ = {"sink": "sinkhorn", "sinki": "sinkhorn", "ls": "ls", "strict": "strict"}


def parse(p: str) -> List[Tuple[str, str, int, int]]:
    out = []
    for st in p.split(","):
        parts = st.split("/")
        head, n = parts[0], int(parts[1])
        order = 2 if len(parts) > 2 and parts[2] == "o2" else 1
        kind, obj = head[0], head[1:]
        if head == "P2":  # GO2c: boundary pass, polish candidates first
            kind, obj = "P", ""
        if kind in "VNLGWXQP" and len(parts) > 2:  # V/<evals>/<sigma deg>
            obj = parts[2]
        assert kind in "ABVNLGWXQP", st
        if kind in "NL" and len(parts) > 2:
            obj = parts[2]
        if kind in "AB":
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


LADDERS = {1.0: (1.0, 0.3, 0.05), 3.0: (3.0, 1.0, 0.3, 0.05)}


def voxel_ladder(
    S: Any, R0: np.ndarray, v: int, evals: int, first: float = 1.0
) -> Tuple[np.ndarray, int]:
    """Per-voxel CMA in legs (steps 1.0, 0.3, 0.05 deg, or 3, 1, 0.3, 0.05 when ``first`` is 3),
    ``evals`` split evenly over the legs, each leg from the best so far."""
    legs = LADDERS[first]
    R, n = np.asarray(R0, float), 0
    fb = float(S.strict(R[None], [v])[0])
    n += 1
    for sig in legs:
        Rn, k = voxel_cma(S, R, v, max(evals // len(legs), 7), sig)
        n += k
        fn = float(S.strict(Rn[None], [v])[0])
        n += 1
        if fn <= fb:
            R, fb = Rn, fn
    return R, n


def neighbour_continuation(
    S: Any, R_now: np.ndarray, s_now: np.ndarray, members: List[int], evals: int, sigma: float,
    sweeps: int = 2,
) -> Tuple[np.ndarray, np.ndarray, float]:
    """Breadth-first from the best-fitting 10% of voxels: each voxel scores its own and its fitted
    in-grain neighbours' orientations (one strict evaluation each), keeps the best, then a
    per-voxel CMA polish. Returns (orientations, strict costs, evaluations per voxel)."""
    from collections import deque

    pos = {v: i for i, v in enumerate(members)}
    R = R_now.copy()
    sc = s_now.copy()
    n_ev = 0
    for _sweep in range(sweeps):
        order = np.argsort(sc)
        seeds = list(order[: max(1, len(order) // 10)])
        done = set(int(i) for i in seeds)
        q = deque(int(i) for i in seeds)
        while q:
            i = q.popleft()
            for u in S.nb[members[i]]:
                j = pos.get(int(u))
                if j is None or j in done:
                    continue
                done.add(j)
                cands = [R[j]] + [R[pos[int(w)]] for w in S.nb[members[j]]
                                  if int(w) in pos and pos[int(w)] in done and pos[int(w)] != j]
                cs = S.strict(np.stack(cands), [members[j]] * len(cands))
                n_ev += len(cands)
                k = int(np.argmin(cs))
                Rb, fb = cands[k], float(cs[k])
                Rp, m = voxel_cma(S, Rb, members[j], evals, sigma)
                n_ev += m
                fp = float(S.strict(Rp[None], [members[j]])[0])
                n_ev += 1
                if fp <= fb:
                    Rb, fb = Rp, fp
                if fb <= sc[j]:
                    R[j], sc[j] = Rb, fb
                q.append(j)
    return R, sc, n_ev / len(members)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--name", required=True)
    ap.add_argument("--pipeline", required=True)
    ap.add_argument("--sample", default="grad")
    ap.add_argument("--stack", default=None, help="GO7: image stack in place of the sample's clean one")
    ap.add_argument("--start", default="paint")
    ap.add_argument("--grains", type=int, default=20)
    ap.add_argument("--true-grains", type=int, nargs="*", default=None)
    ap.add_argument("--group", default="flood", choices=["flood", "true"],
                    help="grain membership: flood fill of the start map (2 deg) or true grain ids")
    ap.add_argument("--sigma0", type=float, default=0.05)
    ap.add_argument("--b-lr", type=float, default=0.01)
    ap.add_argument("--b-mu", type=float, default=1.0)
    ap.add_argument("--rocking_deg", type=float, default=0.0, help="GO7f: strict-cost rocking window (deg)")
    ap.add_argument("--gate", default="strict", choices=["strict", "pois", "either"],
                    help="GO7g: A-stage acceptance on the strict cost, the Poisson deviance, or either")
    ap.add_argument("--model_rocking", type=float, default=0.0, help="GO7g: SpotModel rocking width (deg)")
    ap.add_argument("--store", default=None, help="GO7g: intensity store tag in place of the stack's")
    ap.add_argument("--ratio_off", action="store_true", help="GO7g: region without the area-ratio filter")
    a = ap.parse_args()
    stages = parse(a.pipeline)
    out = CM.OUT / "runs" / a.name
    S = CM.Sample(out, start=a.start, sample=a.sample, stack=a.stack, rocking_deg=a.rocking_deg,
                  model_rocking=a.model_rocking, store=a.store)
    gain = 1244.0
    if a.gate != "strict":  # GO7g: Poisson deviance (GO9, b = 0.5) on the stored intensities
        from probe_poisson import Poisson

        gain = float(np.percentile(S.images.store["val"], 5)) / 50.0
    T = S.T
    gm = np.asarray(T.gmap)
    lab = GF.group_grains(S.R0, S.nb, S.sym, 2.0) if a.group == "flood" else gm.copy()
    picks = a.true_grains or pick_true_grains(T, a.grains)
    Cref: Dict[int, Any] = {}
    if a.sample == "grad" and a.start == "paint":
        Cref = {r["voxel"]: r for r in CM.G.load_gz(CM.C.OUT / "fine" / "runs" / "C" / "refine.json.gz")}
    R_cur = S.R0.copy()
    hooks = None
    nnref = None
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
                               strict_fn=S.strict,
                               region_kw=dict(ratio=None) if a.ratio_off else None)
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
                spec = GF.ObjectiveSpec(kind=OBJ[obj], intensity=(obj == "sinki"))
                if spec.intensity != prob.spec.intensity:  # observed mass weights change
                    prob.y, prob.b, _ = GF.pack_observed(prob.region, spec.intensity, prob.bin_px)
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
            elif kind == "N":
                R_new, _sn, used = neighbour_continuation(S, R_now, s_now, members, n,
                                                          float(obj) if obj else 0.3)
                f_old = f_new = float("nan")
            elif kind == "G":  # gated global search (ray-indexing seed path) for high-cost voxels
                thr = float(obj) if obj else float(n) / 100.0
                if hooks is None:
                    import ray_cache as RC  # GO4d: paired rays from the per-sample cache

                    from icenine import grain_paint as GP

                    _sc, t_rc, loaded = RC.ensure(S, a.sample)
                    print(json.dumps(dict(ray_cache_s=t_rc, loaded=loaded)), flush=True)
                    hooks = GP.ReconstructorHooks(S.bfs, S.mic, np.random.default_rng(0))
                R_new = R_now.copy()
                ev_tot, n_gate, n_found = 0, 0, 0
                for i, v in enumerate(members):
                    if s_now[i] <= thr:
                        continue
                    n_gate += 1
                    Rs, _c, _h, ev = hooks.index_seed(v)
                    ev_tot += ev + 1
                    if Rs is not None:
                        Rs = GF.align_to(Rs[None], R_now[i], S.sym)[0]
                        if float(S.strict(Rs[None], [v])[0]) < s_now[i]:
                            R_new[i] = Rs
                            n_found += 1
                used = ev_tot / len(members)
                info.update(gated=n_gate, improved=n_found, threshold=thr)
                f_old = f_new = float("nan")
            elif kind == "P":  # GO2b: boundary pass, neighbours' orientations from other groups
                polish_first = "P2/" in a.pipeline  # GO2c
                thr = float(obj) if obj else 0.12
                R_new = R_now.copy()
                pos = {v: i for i, v in enumerate(members)}
                n_gate, n_found, ev_tot = 0, 0, 0
                for i, v in enumerate(members):
                    if s_now[i] <= thr:
                        continue
                    n_gate += 1
                    cands = [S.R0[u] for u in S.nb[v] if lab[u] != grp]
                    cands += [R_now[pos[u]] for u in S.nb[v] if u in pos and u != v]
                    if not cands:
                        continue
                    C_ = np.stack(cands)
                    sc = np.array([float(S.strict(C_[j : j + 1], [v])[0]) for j in range(len(C_))])
                    ev_tot += len(C_)
                    if polish_first:  # GO2c: up to 3 distinct candidates, ladder each, take best
                        from icenine.ray_indexing import misorientation_deg as _mis

                        keep: list = []
                        for j in np.argsort(sc):
                            if all(float(_mis(C_[j][None], C_[q][None], S.sym)[0]) > 0.5 for q in keep):
                                keep.append(int(j))
                            if len(keep) == 3:
                                break
                        best_R, best_s = None, s_now[i]
                        for j in keep:
                            Rl, kl = voxel_ladder(S, C_[j], v, n, 1.0)
                            ev_tot += kl
                            sl = float(S.strict(Rl[None], [v])[0])
                            if sl < best_s:
                                best_R, best_s = Rl, sl
                        if best_R is not None:
                            R_new[i] = best_R
                            n_found += 1
                        continue
                    j = int(np.argmin(sc))
                    if sc[j] < s_now[i]:
                        Rb, kb = voxel_cma(S, C_[j], v, n)
                        ev_tot += kb
                        sb = float(S.strict(Rb[None], [v])[0])
                        R_new[i] = Rb if sb <= sc[j] else C_[j]
                        n_found += 1
                used = ev_tot / len(members)
                info.update(gated=n_gate, improved=n_found, threshold=thr)
                f_old = f_new = float("nan")
            elif kind == "Q":  # GO8b: GNLayerNet refiner, rounds kept only if the strict cost drops
                thr = float(obj) if obj else 0.0
                if nnref is None:
                    import nn_refine as NR

                    nnref = NR.NNRefiner(S, str(CM.ICE / "scripts" / "toy_orientation_sweep_model_realistic_s0.pt"))
                R_new = R_now.copy()
                s_cur = s_now.copy()
                n_gate, n_found, calls0 = 0, 0, nnref.n_calls
                for i, v in enumerate(members):
                    if s_now[i] <= thr:
                        continue
                    n_gate += 1
                    for _r in range(n):
                        Rq = nnref.propose(v, R_new[i])
                        if Rq is None:
                            break
                        sq = float(S.strict(Rq[None], [v])[0])
                        if sq < s_cur[i]:
                            R_new[i], s_cur[i] = Rq, sq
                        else:
                            break
                    if s_cur[i] < s_now[i]:
                        n_found += 1
                used = (nnref.n_calls - calls0) / len(members)
                info.update(gated=n_gate, improved=n_found, threshold=thr)
                f_old = f_new = float("nan")
            elif kind in "WX":  # GO4g: wider seed path (W) or grid search (X) above the gate
                thr = float(obj) if obj else 0.15
                if hooks is None:
                    import ray_cache as RC

                    from icenine import grain_paint as GP

                    RC.ensure(S, a.sample)
                    hooks = GP.ReconstructorHooks(S.bfs, S.mic, np.random.default_rng(0))
                R_new = R_now.copy()
                ev_tot, n_gate, n_found = 0, 0, 0
                for i, v in enumerate(members):
                    if s_now[i] <= thr:
                        continue
                    n_gate += 1
                    if kind == "W":
                        sres = S.bfs.reconstructor.index_seed(
                            hooks.vertices(v), S.mic.voxels[v].phase, n_cost=20, n_refine=5
                        )
                        Rs, ev = sres.orientation, int(sres.evals)
                    else:
                        Rs, _c, _h, ev = hooks.grid_seed(v)
                    ev_tot += ev + 1
                    if Rs is not None:
                        Rs = GF.align_to(np.asarray(Rs, np.float64)[None], R_now[i], S.sym)[0]
                        if float(S.strict(Rs[None], [v])[0]) < s_now[i]:
                            R_new[i] = Rs
                            n_found += 1
                used = ev_tot / len(members)
                info.update(gated=n_gate, improved=n_found, threshold=thr)
                f_old = f_new = float("nan")
            elif kind == "L":
                R_new = R_now.copy()
                used_v = []
                for i, v in enumerate(members):
                    R_new[i], k = voxel_ladder(S, R_now[i], v, n, float(obj) if obj else 1.0)
                    used_v.append(k)
                used = float(np.mean(used_v))
                f_old = f_new = float("nan")
            else:  # per-voxel strict CMA
                R_new = R_now.copy()
                used_v = []
                for i, v in enumerate(members):
                    R_new[i], k = voxel_cma(S, R_now[i], v, n, float(obj) if obj else 0.02)
                    used_v.append(k)
                used = float(np.mean(used_v))
                f_old = f_new = float("nan")
            s_new = S.strict(R_new, members)
            if kind in "VNLGWXQP":  # per voxel: keep a voxel's fit only if its strict cost did not rise
                keep = s_new <= s_now
                R_new = np.where(keep[:, None, None], R_new, R_now)
                s_new = np.where(keep, s_new, s_now)
                acc = True
                info["voxels_kept"] = int(keep.sum())
            else:
                acc = (f_new <= f_old) and (s_new.mean() <= s_now.mean())
                if a.gate != "strict":  # GO7g: Poisson deviance in place of, or alongside, the strict gate
                    pois = Poisson(prob, gain, 0.5)
                    with torch.no_grad():
                        p_old = pois(torch.as_tensor(R_now, dtype=GF.DT))
                        p_new = pois(torch.as_tensor(R_new, dtype=GF.DT))
                    acc_p = p_new <= p_old
                    info.update(pois_old=float(p_old), pois_new=float(p_new), acc_strict=bool(acc))
                    acc = acc_p if a.gate == "pois" else (acc or acc_p)
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
