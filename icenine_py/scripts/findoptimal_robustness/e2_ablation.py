#!/usr/bin/env python3
"""Which feature groups carry the classifier? GBT (4 voxel-disjoint folds) with feature groups
removed: CSL-aware features, per-family/detector breakdown, pixel-radius-3 features."""

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import common as C  # noqa: E402
import e2_models as M  # noqa: E402
import features as F  # noqa: E402


def main():
    D = M.load_dataset()
    nm = ["log_n_pairs", "hit0", "hit_any", "hit1", "hit3"]
    nm += [f"fam{f}_{k}" for f in range(8) for k in ("hit0", "hit3", "frac")]
    nm += [f"det{d}_{k}" for d in range(2) for k in ("hit0", "hit3", "frac")]
    for sg in F.SIGMAS:
        nm += [
            f"S{sg}_{k}"
            for k in (
                "shared_frac",
                "hit0_shared",
                "hit0_nonshared_min",
                "hit3_nonshared_min",
                "hit0_nonshared_mean",
            )
        ]
    nm += ["cost_local", "cost_global3"]
    assert len(nm) == D["X"].shape[1]
    groups = {
        "all features": [],
        "no cost features": ["cost_local", "cost_global3"],
        "no CSL-aware features": [n for n in nm if n.startswith("S")],
        "no CSL-aware, no cost": [n for n in nm if n.startswith("S")]
        + ["cost_local", "cost_global3"],
        "no pixel-radius-3 features": [
            n for n in nm if "hit3" in n or n in ("hit1", "cost_global3")
        ],
        "only overall hit fractions (hit0, hit_any, hit1, hit3)": [
            n for n in nm if n not in ("log_n_pairs", "hit0", "hit_any", "hit1", "hit3")
        ],
    }
    fold = M.fold_of(D["vpos"])
    lines = []
    res = {}
    for gname, drop in groups.items():
        keep = [i for i, n in enumerate(nm) if n not in drop]
        sc = np.full(len(D["err"]), np.nan)
        for f in range(M.FOLDS):
            tr, te = fold != f, fold == f
            m = M.make_model("GBT")
            m.fit(D["X"][tr][:, keep], (D["err"][tr] < 3.0).astype(int))
            sc[te] = m.predict_proba(D["X"][te][:, keep])[:, 1]
        r = M.metrics(D, {"GBT": sc})["GBT"]
        res[gname] = r
        lines.append(
            f"{gname:58s} AUC A {r['auc_A']:.4f}  B {r['auc_B']:.3f}  C {r['auc_C']:.4f}  prune recall {r['prune_recall']:.3f}  final precision {r['final_precision']:.3f}  ({len(keep)} features)"
        )
    (C.OUT_DIR / "e2_ablation.txt").write_text("\n".join(lines) + "\n")
    C.save_json(C.OUT_DIR / "e2_ablation.json", res)
    print("\n".join(lines))


if __name__ == "__main__":
    main()
