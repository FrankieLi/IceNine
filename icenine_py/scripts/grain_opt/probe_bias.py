"""Bias of the objectives: on grains of level 0 (one true orientation per grain), CMA over the grain
rotation only (3 numbers) starting AT the truth; an unbiased objective stays at the truth. Reports
the optimum's distance to the truth, objective at truth and optimum, and the mean strict cost.

Usage (from icenine_py/): uv run python scripts/grain_opt/probe_bias.py --grains 381 302 41
"""

import argparse
import json
import time

import common as CM  # noqa: I001
import numpy as np
import torch

from icenine import grain_fit as GF

VARIANTS = {
    "sink_r3_b1": dict(kind="sinkhorn", reach_px=3.0, blur_px=1.0),
    "sink_bal_b1": dict(kind="sinkhorn", reach_px=None, blur_px=1.0),
    "sink_r3_b05": dict(kind="sinkhorn", reach_px=3.0, blur_px=0.5),
    "sink_r3_b1_dil": dict(kind="sinkhorn", reach_px=3.0, blur_px=1.0, dilate=0.29),
    "sink_r3_b05_dil": dict(kind="sinkhorn", reach_px=3.0, blur_px=0.5, dilate=0.29),
    "sink_r3_b1_int": dict(kind="sinkhorn", reach_px=3.0, blur_px=1.0, intensity=True),
    "ls": dict(kind="ls"),
}


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--grains", type=int, nargs="+", default=[381, 302, 41])
    ap.add_argument("--variants", nargs="*", default=list(VARIANTS))
    ap.add_argument("--evals", type=int, default=150)
    ap.add_argument("--tag", default="bias")
    a = ap.parse_args()
    S = CM.Sample(CM.OUT / "probe")
    T = S.T
    gm = np.asarray(T.gmap)
    out = []
    for g in a.grains:
        members = [int(v) for v in np.nonzero(gm == g)[0]]
        shell = GF.shell_of(members, S.nb)
        Rt = T.R[members]
        st_truth = float(S.strict(Rt, members).mean())
        for name in a.variants:
            spec = GF.ObjectiveSpec(**VARIANTS[name])
            prob = GF.GrainProblem(S.model, S.images, S.R0, S.V, S.xy, members, shell, S.sym,
                                   spec=spec)
            Rt_al = GF.align_to(Rt, prob.field.R_ref, S.sym)
            w_t = GF.matrix_to_rotvec(Rt_al[:1] @ prob.field.R_ref.T)[0]
            th0 = np.zeros(prob.field.n_params)
            th0[:3] = w_t
            # rotation only: wrap the field so CMA moves w0 alone
            import cma

            def f(w):
                th = th0.copy()
                th[:3] = w
                with torch.no_grad():
                    return float(prob.objective(prob.field.orientations(torch.as_tensor(th))))

            t0 = time.time()
            es = cma.CMAEvolutionStrategy(list(w_t), np.radians(0.03),
                                          {"seed": 1, "maxfevals": a.evals, "verbose": -9})
            es.optimize(f)
            w = np.asarray(es.result.xbest)
            th = th0.copy()
            th[:3] = w
            with torch.no_grad():
                R_opt = prob.field.orientations(torch.as_tensor(th)).numpy()
            off = float(np.median(T.err(R_opt, members)))
            row = dict(grain=g, variant=name, n_spots=prob.n_spots, bin_px=prob.bin_px,
                       n_sub=prob.n_sub, offset_deg=off, f_truth=f(w_t), f_opt=f(w),
                       strict_truth=st_truth, strict_opt=float(S.strict(R_opt, members).mean()),
                       t_s=time.time() - t0)
            print(json.dumps(row), flush=True)
            with open(CM.OUT / "probe" / f"{a.tag}.jsonl", "a") as fh:
                fh.write(json.dumps(row) + "\n")
            out.append(row)
    (CM.OUT / "probe" / f"{a.tag}.json").write_text(json.dumps(out, indent=1) + "\n")


if __name__ == "__main__":
    main()
