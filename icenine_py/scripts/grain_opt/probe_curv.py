"""GO1b (I4): noise floor against curvature of each objective near the ground truth.

Grains of ``grad``, grain rotation only. Along 3 orthonormal rotation axes, each objective at 21
offsets in [-0.05, 0.05] deg around the ground truth; per axis a least-squares parabola gives the
curvature k, the vertex d* and the residual rms sigma; r = sqrt(sigma / k) is the noise-limited
resolution. Writes ``benchmarks/grain_opt/go1b/curv.json`` (per axis rows and the scans).

Usage (from icenine_py/): uv run python scripts/grain_opt/probe_curv.py [--grains 381 302]
"""

import argparse
import json
import time

import common as CM  # noqa: I001
import numpy as np
import torch
from probe_poisson import Poisson

from icenine import grain_fit as GF

VARIANTS = {
    "strict": dict(kind="strict"),
    "ls": dict(kind="ls"),
    "sink_r3_b1_int": dict(kind="sinkhorn", reach_px=3.0, blur_px=1.0, intensity=True),
}


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--grains", type=int, nargs="+", default=[381, 302])
    ap.add_argument("--n", type=int, default=21)
    ap.add_argument("--span", type=float, default=0.05)
    a = ap.parse_args()
    S = CM.Sample(CM.OUT / "go1b" / "run")
    T = S.T
    gm = np.asarray(T.gmap)
    allv = S.images.store["val"] if getattr(S.images, "store", None) else None
    gain = float(np.percentile(allv, 5)) / 50.0 if allv is not None else 1244.0
    d = np.linspace(-a.span, a.span, a.n)
    rows = []
    for g in a.grains:
        members = [int(v) for v in np.nonzero(gm == g)[0]]
        shell = GF.shell_of(members, S.nb)
        Rt = T.R[members]
        for name in list(VARIANTS) + ["pois"]:
            spec = GF.ObjectiveSpec(**VARIANTS.get(name, dict(kind="ls")))
            prob = GF.GrainProblem(S.model, S.images, S.R0, S.V, S.xy, members, shell, S.sym,
                                   spec=spec, strict_fn=S.strict)
            Rt_al = GF.align_to(Rt, prob.field.R_ref, S.sym)
            w_t = GF.matrix_to_rotvec(Rt_al[:1] @ prob.field.R_ref.T)[0]
            th0 = np.zeros(prob.field.n_params)
            pois = Poisson(prob, gain, 0.5) if name == "pois" else None

            def f(w: np.ndarray) -> float:
                th = th0.copy()
                th[:3] = w
                with torch.no_grad():
                    R = prob.field.orientations(torch.as_tensor(th))
                    return float(pois(R)) if pois is not None else float(prob.objective(R))

            t0 = time.time()
            for ax in range(3):
                e = np.zeros(3)
                e[ax] = 1.0
                ys = np.array([f(w_t + np.radians(x) * e) for x in d])
                c2, c1, c0 = np.polyfit(d, ys, 2)
                res = ys - np.polyval([c2, c1, c0], d)
                sig = float(res.std())
                k = float(c2)
                vertex = float(-c1 / (2 * c2)) if c2 > 0 else float("nan")
                r = float(np.sqrt(sig / k)) if k > 0 else float("inf")
                row = dict(grain=g, objective=name, axis=ax, k=k, vertex_deg=vertex, sigma=sig,
                           r_deg=r, f_truth=float(ys[a.n // 2]), scan=ys.tolist())
                rows.append(row)
                print(json.dumps({k_: v for k_, v in row.items() if k_ != "scan"}), flush=True)
            print(json.dumps(dict(grain=g, objective=name, t_s=time.time() - t0)), flush=True)
    out = CM.OUT / "go1b"
    out.mkdir(exist_ok=True)
    (out / "curv.json").write_text(json.dumps(dict(offsets_deg=d.tolist(), rows=rows), indent=1) + "\n")


if __name__ == "__main__":
    main()
