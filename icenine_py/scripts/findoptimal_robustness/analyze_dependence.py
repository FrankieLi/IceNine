#!/usr/bin/env python3
"""E0 dependence of the wrong rate on how many of the truth's reflections are shared with its
Sigma relatives (features.py: S<sigma>_shared_frac = fraction of the truth's eligible (peak,
detector) pairs whose reflection is also produced by at least one Sigma relative), plus the number
of eligible pairs. Needs cache/e2 (truth_features) and benchmarks/.../e0_runs.npz."""

import sys
from pathlib import Path

import numpy as np
from scipy.stats import spearmanr

sys.path.insert(0, str(Path(__file__).resolve().parent))
import common as C  # noqa: E402
import features as F  # noqa: E402


def main():
    runs = dict(np.load(C.OUT_DIR / "e0_runs.npz"))
    tf = {}
    for v in np.unique(runs["vidx"]):
        for vi, var in enumerate(C.VARIANTS):
            f = C.CACHE_DIR / "e2" / f"v{v}_{var}.npz"
            if f.exists():
                tf[(int(v), vi)] = np.load(f)["truth_features"]
    nm = F.feature_names_for_width(len(next(iter(tf.values()))))
    lines = []
    out = {}
    for vi, var in enumerate(C.VARIANTS):
        sel = runs["variant"] == vi
        vv = runs["vidx"][sel]
        w = (runs["final_err"][sel] > C.WRONG_DEG).astype(float)
        cls = runs["csl_label"][sel]
        lines.append(f"[{var}]")
        for col in ["log_n_pairs"] + [f"S{s}_shared_frac" for s in F.SIGMAS]:
            xs = np.array(
                [tf[(int(v), vi)][nm.index(col)] if (int(v), vi) in tf else np.nan for v in vv]
            )
            ok = np.isfinite(xs)
            rho, p = spearmanr(xs[ok], w[ok])
            q = np.quantile(xs[ok], [0, 1 / 3, 2 / 3, 1])
            tab = []
            for a, b in zip(q[:-1], q[1:]):
                m = ok & (xs >= a) & (xs <= b)
                tab.append((float(a), float(b), int(m.sum()), float(w[m].mean())))
            out[f"{var}/{col}"] = dict(spearman=float(rho), p=float(p), tertiles=tab)
            lines.append(
                f"  {col:18s} Spearman(wrong, x) = {rho:+.3f} (p={p:.2g}, n={ok.sum()} runs); "
                "wrong rate by tertile: "
                + ", ".join(f"[{a:.2f},{b:.2f}] {r:.2f}" for a, b, n, r in tab)
            )
        # Sigma3-wrong runs vs Sigma3 shared frac: only wrong runs classed Sigma3
        xs3 = np.array(
            [
                tf[(int(v), vi)][nm.index("S3_shared_frac")] if (int(v), vi) in tf else np.nan
                for v in vv
            ]
        )
        s3 = cls == "3"
        lines.append(
            f"  mean S3_shared_frac: Sigma3-wrong runs {np.nanmean(xs3[s3]):.3f} (n={s3.sum()}), "
            f"all other runs {np.nanmean(xs3[~s3]):.3f}"
        )
    (C.OUT_DIR / "e0_dependence.txt").write_text("\n".join(lines) + "\n")
    C.save_json(C.OUT_DIR / "e0_dependence.json", out)
    print("\n".join(lines))


if __name__ == "__main__":
    main()
