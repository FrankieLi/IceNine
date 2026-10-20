"""Task F diagnostic: why F4 fails from far starts (a check of one explanation, not a proof).

For F4 runs from 5 deg starts (and 3 deg), on the full sample (six voxels x 20 directions) and
on the
E0 voxels of the benchmark (a subset), compare: (a) the share of predicted peaks whose nearest
observed ray at the START is the ray they have at the TRUTH (the matching F4's first step uses), and
(b) the E4 cost (tau 0.25 deg) at the F4 end point against the E4 cost at the truth (a lower value
at the end point means the end point is a deeper minimum of E4 than the truth; a higher value means
F4 stopped in a local minimum or with a wrong matching). Runs are split into captured (error under
1 deg) and not captured.

  uv run python scripts/eval_reduction/outward_diag.py full
  uv run python scripts/eval_reduction/outward_diag.py e0 --n-voxels 8
Writes benchmarks/eval_reduction/F/diag/<kind>_w<window>.json (frame windows 1 and 5).
"""

import argparse
import json
import math
import sys
from pathlib import Path
from typing import Any, Dict, List

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import outward_bench as B  # noqa: E402


def peak_ray_map(gn: Any, R: np.ndarray) -> Dict[int, int]:
    """{peak key (reflection * 2 + branch): nearest attached ray} for the weighted peaks of R."""
    import torch

    with torch.no_grad():
        pr = gn.sc.predict(torch.as_tensor(R, dtype=torch.float64)[None], gn.verts, gn.phase)
        W = gn._weights(pr)
        arg, _ = gn._match(pr)
    key = (pr.refl * 2 + pr.branch).numpy()
    use = (W > 0).numpy()
    return {int(k): int(a) for k, a, u in zip(key, arg.numpy(), use) if u}


def run(kind: str, n_voxels: int, offsets: List[int], window: int) -> Dict[str, Any]:
    from icenine.outward_search import angle_deg

    B.init_worker(kind)
    if kind == "full":
        from parity_common import SIX_VOXELS

        vox = list(SIX_VOXELS)
    else:
        vs = json.loads((B.OUT / "voxel_sets.json").read_text())
        vox = vs["near"][: n_voxels // 2] + vs["far"][: n_voxels - n_voxels // 2]
    rows: Dict[str, List[List[float]]] = {}
    att: List[List[int]] = []  # per voxel: attached rays, weighted peaks at the truth
    tau = math.radians(0.25)
    for v in vox:
        env = B.voxel_env(v)
        with B.Tracer(env["lf"], env["R_true"], env["verts"], env["phase"]) as tr:
            gn = B.make_gn(env, tr, window=window)
            gn._tables()
            att.append([gn._tables()["n"], len(peak_ray_map(gn, env["R_true"]))])
            truth_map = peak_ray_map(gn, env["R_true"])
            e4_truth = gn.e4_value(env["R_true"], tau)
            for oi in offsets:
                starts = B.starts_for(env, B.OFFSETS[oi], 20)
                for R0 in starts:
                    res = gn.solve(R0)
                    m0 = peak_ray_map(gn, R0)
                    common = [k for k in m0 if k in truth_map]
                    same = (
                        float(np.mean([m0[k] == truth_map[k] for k in common]))
                        if common
                        else math.nan
                    )
                    err = angle_deg(res.R, env["R_true"])
                    rows.setdefault(f"{B.OFFSETS[oi]:g}", []).append(
                        [err, same, gn.e4_value(res.R, tau) - e4_truth, float(res.accepted), v]
                    )
    out: Dict[str, Any] = {"kind": kind, "tau_deg": 0.25, "frame_window": window}
    out["median_attached_rays"] = float(np.median([a[0] for a in att]))
    out["median_weighted_peaks_at_truth"] = float(np.median([a[1] for a in att]))
    for off, r in rows.items():
        a = np.array(r)
        ok = a[:, 0] < 1.0 - B.CAPTURE_TOL
        d: Dict[str, Any] = dict(n=len(a), captured=int(ok.sum()))
        for name, sel in (("captured", ok), ("not_captured", ~ok)):
            if sel.any():
                d[name] = dict(
                    n=int(sel.sum()),
                    median_share_same_ray_at_start=float(np.nanmedian(a[sel, 1])),
                    median_e4_end_minus_truth=float(np.median(a[sel, 2])),
                    share_end_below_truth=float((a[sel, 2] < -1e-9).mean()),
                    share_accepted=float(a[sel, 3].mean()),
                )
        out[off] = d
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("kind", choices=["e0", "full"])
    ap.add_argument("--n-voxels", type=int, default=8)
    a = ap.parse_args()
    d = B.OUT / "diag"
    d.mkdir(exist_ok=True)
    for w in (1, 5):
        out = run(a.kind, a.n_voxels, [4, 5], w)
        (d / f"{a.kind}_w{w}.json").write_text(json.dumps(out, indent=1) + "\n")
        print(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
