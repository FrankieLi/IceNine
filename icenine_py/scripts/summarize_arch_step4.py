#!/usr/bin/env python3
"""Step 4 table: Gauss-Newton (plain / Huber) and GNLayerNet runs on clean and corrupted test sets.

Reads the saved predictions (npz) rather than the result jsons, so every method is scored by
the same code on the same samples. Net predictions: <dir>/<run>_s<seed>.npz for the clean
variant and <dir>/<run>_s<seed>_<variant>.npz otherwise (what train_toy_orientation_nn.py
--save-predictions writes). GN predictions: <dir>/dis_pred_<gn|huberC>_<variant>.npz.

Usage:
  cd icenine_py
  uv run python scripts/summarize_arch_step4.py --test scripts/toy_orientation_arch_dis_test.pt \
      --huber 1 --runs "GNLayerNet (corrupted-trained)=dis_res_gn_k3_corr" \
      "paired=dis_res_gnpair_k3_corr" --seeds 0 1 \
      --out benchmarks/toy_orientation_arch/step4_summary.json

Per method, variant and group (in-dist / held-out voxels): median misorientation angle pooled over
the four |delta| bins and per bin (0.1/0.25/0.5/1.0), the fraction of cases below 0.1 deg, z and
perpendicular RMS, mean Mahalanobis^2 (nets), and over the 30 voxels the correlation of the
per-voxel median error with r_perp. Nets are averaged over seeds; per-seed pooled medians are
listed.
"""

import argparse
import json
from pathlib import Path
from typing import Any, Dict, Optional

import numpy as np
import torch

from icenine.orientation_eval import error_summary

VARIANTS = ["clean", "neighbours", "noise", "all"]
MAGS = [0.1, 0.25, 0.5, 1.0]


def stats(
    pred: np.ndarray,
    truth: np.ndarray,
    mags: np.ndarray,
    vid: np.ndarray,
    held: np.ndarray,
    r_perp: np.ndarray,
    chol: Optional[np.ndarray] = None,
) -> Dict[str, Any]:
    ok = np.isfinite(pred).all(axis=1)
    out = {"n_nan": int((~ok).sum())}
    ev = np.full(len(r_perp), np.nan)
    for v in range(len(r_perp)):
        m = (vid == v) & ok
        if m.any():
            ev[v] = error_summary(pred[m], truth[m])["median_angle"]
    hv = held
    out["corr_err_rperp"] = float(np.corrcoef(r_perp, ev)[0, 1])
    for name, gsel in (("in-dist", ~hv), ("held-out", hv)):
        sel = gsel[vid] & ok
        s = error_summary(pred[sel], truth[sel])
        g = {
            "median_angle": s["median_angle"],
            "success_0p1": s["success_0p1"],
            "rms_z": s["rms_z"],
            "rms_perp": s["rms_perp"],
            "median_voxel_err": float(np.nanmedian(ev[gsel])),
            "per_mag_median": [],
        }
        for mg in MAGS:
            ms = sel & np.isclose(mags, mg, atol=1e-3)
            g["per_mag_median"].append(error_summary(pred[ms], truth[ms])["median_angle"])
        if chol is not None:
            c = chol[sel]
            r = (truth[sel] - pred[sel])[:, :, None]
            z = np.linalg.solve(c, r)[:, :, 0]
            g["mahalanobis_sq"] = float((z**2).sum(axis=1).mean())
        out[name] = g
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--test", required=True)
    ap.add_argument("--dir", default="benchmarks/toy_orientation_arch")
    ap.add_argument("--huber", default="1")
    ap.add_argument(
        "--runs", nargs="+", required=True, help="LABEL=run-name-prefix (without _s<seed>)"
    )
    ap.add_argument("--seeds", nargs="+", type=int, default=[0, 1])
    ap.add_argument("--out", default=None)
    args = ap.parse_args()
    d = Path(args.dir)

    meta_path = d / "dis_test_meta.npz"
    test_file = Path(args.test).resolve()
    test_key = f"{test_file}:{test_file.stat().st_size}"  # the cache is only valid for this file
    meta = None
    if meta_path.exists():
        cached = np.load(meta_path)
        if "test_key" in cached.files and str(cached["test_key"]) == test_key:
            meta = cached
    if meta is not None:
        vid, held, r_perp = meta["voxel_id"], meta["held_out"], meta["r_perp_um"]
    else:
        te = torch.load(args.test)
        vid, held, r_perp = te["voxel_id"].numpy(), te["held_out"].numpy(), te["r_perp_um"].numpy()
        np.savez(meta_path, voxel_id=vid, held_out=held, r_perp_um=r_perp, test_key=test_key)

    results = {}
    for variant in VARIANTS:
        results[variant] = {}
        for label, tag in (
            ("plain GN", "gn"),
            (f"Huber GN (c={args.huber})", f"huber{args.huber}"),
        ):
            z = np.load(d / f"dis_pred_{tag}_{variant}.npz")
            results[variant][label] = stats(
                z["pred_deg"], z["truth_deg"], z["magnitudes_deg"], vid, held, r_perp
            )
        for spec in args.runs:
            label, prefix = spec.split("=", 1)
            per_seed = []
            for seed in args.seeds:
                f = d / (
                    f"{prefix}_s{seed}" + ("" if variant == "clean" else f"_{variant}") + ".npz"
                )
                if not f.exists():
                    continue
                z = np.load(f)
                per_seed.append(
                    stats(
                        z["pred_deg"],
                        z["truth_deg"],
                        z["magnitudes_deg"],
                        vid,
                        held,
                        r_perp,
                        z["chol"],
                    )
                )
            if not per_seed:
                continue
            results[variant][label] = {"seeds": per_seed, "n_seeds": len(per_seed)}

    def avg(vals):
        return float(np.mean(vals))

    print(
        f"{'variant':<11}{'method':<40}{'group':<9}{'median':>7}{'<0.1':>6}{'zRMS':>7}{'pRMS':>7}"
        f"{'Mah2':>8}  per-bin median (0.1/0.25/0.5/1.0)"
    )
    for variant in VARIANTS:
        for label, r in results[variant].items():
            for g in ("in-dist", "held-out"):
                rows = r["seeds"] if "seeds" in r else [r]
                m = lambda k: avg([x[g][k] for x in rows])  # noqa: E731
                pm = np.mean([x[g]["per_mag_median"] for x in rows], axis=0)
                mah = f"{m('mahalanobis_sq'):8.1f}" if "mahalanobis_sq" in rows[0][g] else " " * 8
                seeds = ""
                if "seeds" in r:
                    seeds = "  seeds: " + " ".join(f"{x[g]['median_angle']:.4f}" for x in rows)
                name = label + (f" (n={len(rows)})" if "seeds" in r else "")
                print(
                    f"{variant:<11}{name:<40}{g:<9}{m('median_angle'):7.4f}{m('success_0p1'):6.2f}"
                    f"{m('rms_z'):7.3f}{m('rms_perp'):7.3f}"
                    f"{mah}  {'/'.join(f'{x:.3f}' for x in pm)}{seeds}"
                )
            rows = r["seeds"] if "seeds" in r else [r]
            cs = [x["corr_err_rperp"] for x in rows]
            print(
                f"{'':<11}{'':<40}corr(voxel median error, r_perp) = "
                + ", ".join(f"{c:+.2f}" for c in cs)
                + (" (per seed)" if len(cs) > 1 else "")
                + f"; median voxel error in-dist "
                f"{avg([x['in-dist']['median_voxel_err'] for x in rows]):.4f}"
                f" held-out {avg([x['held-out']['median_voxel_err'] for x in rows]):.4f}"
            )
    if args.out:
        Path(args.out).write_text(json.dumps(results, indent=1))
        print(f"saved {args.out}")


if __name__ == "__main__":
    main()
