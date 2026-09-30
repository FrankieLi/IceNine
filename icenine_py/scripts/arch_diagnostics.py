#!/usr/bin/env python3
"""Step 1 diagnostics for the parallax architecture plan (no training).

Per voxel of a test set: detector/pair counts, Gauss-Newton (GN) restricted to detector 0 /
detector 1 / both (error per axis), one undamped GN step (the linear model a GN layer computes),
and the conditioning of J^T W J at the nominal orientation (sigma per axis, condition number,
stage-axis share of the weakest eigenvector).

Usage:
  cd icenine_py
  uv run python scripts/arch_diagnostics.py --test scripts/toy_orientation_stage3_multi_test.pt \
      --out benchmarks/toy_orientation_arch/diag_multi.json [--max-samples 400]
"""

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
from dataset_problems import iter_voxels, load_meta  # noqa: E402


def axis_stats(pred, truth):
    from icenine.orientation_eval import error_summary

    ok = np.isfinite(pred).all(axis=1)
    if ok.sum() == 0:
        return None
    s = error_summary(pred[ok], truth[ok])
    return dict(
        n=int(ok.sum()),
        median_angle=float(s["median_angle"]),
        rms_z=float(s["rms_z"]),
        rms_perp=float(s["rms_perp"]),
    )


def main():
    from icenine.orientation_baselines import CentroidGaussNewton, extract_measurements
    from icenine.orientation_eval import pair_index

    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--test", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--max-samples", type=int, default=0, help="per voxel (0 = all)")
    args = ap.parse_args()
    out_path = Path(args.out).resolve()
    data = load_meta(args.test)
    truth_all = data["offsets_deg"].double().numpy()
    configs = {"both": None, "det0": [0], "det1": [1]}
    rows = []
    t0 = time.time()
    for vox in iter_voxels(data):
        obs, spec, n_pk = vox["obs"], vox["spec"], vox["n_pk"]
        pidx = pair_index(vox["problem"]["roi_list"])
        det = obs.det_idx.numpy()
        gn = CentroidGaussNewton(obs)
        idx = vox["sample_idx"]
        if args.max_samples and len(idx) > args.max_samples:
            idx = idx[:: len(idx) // args.max_samples][: args.max_samples]
        preds = {k: np.full((len(idx), 3), np.nan) for k in configs}
        preds["linear_both"] = np.full((len(idx), 3), np.nan)
        both_present = []
        used_n = {k: [] for k in configs}
        for j, n in enumerate(idx):
            w = data["windows"][n][:n_pk]
            for name, dets in configs.items():
                meas = extract_measurements(w, spec, obs, detectors=dets)
                used_n[name].append(int(meas.used.sum()))
                r = gn.solve(meas)
                preds[name][j] = r["delta"]
                if name == "both":
                    try:
                        preds["linear_both"][j] = gn.solve_linear(meas)
                    except np.linalg.LinAlgError:
                        pass
                    u = meas.used
                    has_p = pidx >= 0
                    both_present.append(int((u & has_p & u[np.clip(pidx, 0, None)]).sum()))
        truth = truth_all[idx]
        # conditioning of J^T W J at nominal for each detector config (nominal windows: the
        # nominal orientation's own spots, all present by construction)
        from icenine.orientation_baselines import Measurements
        from icenine.orientation_baselines import frame_center_omega

        nom = obs.observe(__import__("torch").zeros(1, 3, dtype=obs.dtype))
        cent = nom.verts[0].mean(dim=1).numpy()
        omg = frame_center_omega(obs, nom.frame[0].numpy())
        cond = {}
        for name, dets in configs.items():
            mask = np.ones(n_pk, dtype=bool) if dets is None else np.isin(det, dets)
            meas = Measurements(used=mask & np.isfinite(omg), omega=omg, centroid=cent)
            A = gn.information(meas)
            ev, evec = np.linalg.eigh(A)
            cov = np.linalg.inv(A)
            cond[name] = dict(
                sigma_xyz_deg=np.sqrt(np.diag(cov)).tolist(),
                eig_min=float(ev[0]),
                eig_max=float(ev[-1]),
                cond=float(ev[-1] / ev[0]),
                weakest_eigvec_z_share=float(evec[2, 0] ** 2),
            )
        row = dict(
            voxel=vox["v"],
            r_perp_um=vox["r_perp_um"],
            held_out=bool(data["held_out"][vox["v"]]) if "held_out" in data else False,
            n_entries=n_pk,
            n_det=[int((det == d).sum()) for d in (0, 1)],
            n_pair_entries=int((pidx >= 0).sum()),
            n_pairs=int((pidx >= 0).sum() // 2),
            mean_both_present_entries=float(np.mean(both_present)),
            mean_used=dict((k, float(np.mean(v))) for k, v in used_n.items()),
            gn={k: axis_stats(preds[k], truth) for k in preds},
            info=cond,
            n_samples=len(idx),
        )
        rows.append(row)
        print(
            f"voxel {row['voxel']:2d} r_perp {row['r_perp_um']:5.0f} entries {n_pk} pairs {row['n_pairs']} "
            + " | ".join(
                f"{k}: med {row['gn'][k]['median_angle']:.4f} z {row['gn'][k]['rms_z']:.4f} "
                f"perp {row['gn'][k]['rms_perp']:.4f}"
                for k in row["gn"]
            )
            + f"  cond(both) {cond['both']['cond']:.1e} ({time.time() - t0:.0f}s)",
            flush=True,
        )
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(rows, indent=1))
    print(f"saved {out_path}")


if __name__ == "__main__":
    main()
