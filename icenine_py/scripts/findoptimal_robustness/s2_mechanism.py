#!/usr/bin/env python3
"""Why is the truth basin outranked at pruning (S2)? For runs lost at level L (S2@L): the post-MC
local cost of the best basin(3) candidate against the cost of the candidate that ranked first, the
error of that first-ranked candidate, and how much of the cost gap a perfect final alignment would
close (the cost at the truth, E0 cost_true)."""

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import analyze_e0 as A  # noqa: E402
import common as C  # noqa: E402
import csl  # noqa: E402


def main():
    info = dict(np.load(C.OUT_DIR / "voxels.npz"))
    rows = A.load_all(C.CACHE_DIR / "e0", info)
    lines = []
    for var in C.VARIANTS:
        stats = {0: [], 1: [], 2: []}
        for r in rows:
            if r["variant"] != var or not r["class"].startswith("S2@"):
                continue
            L = int(r["class"].split("@")[1])
            d = np.load(C.CACHE_DIR / "e0" / f"v{r['vidx']}_{var}.npz")
            p = f"s{r['seed']}_L{L}_qmc_"
            e = C.err_deg(d[p + "R"], d["R_true"])
            cost = d[p + "cost"]
            if not (e < 3.0).any():
                continue
            bb = int(np.where(e < 3.0)[0][0])
            c0 = csl.csl_classify(d["R_true"], d[p + "R"][0])
            stats[L].append((cost[bb], cost[0], e[bb], e[0], c0["label"], float(d["cost_true"])))
        for L, v in stats.items():
            if not v:
                continue
            a = np.array([(x[0], x[1], x[2], x[3], x[5]) for x in v], float)
            labs = [x[4] for x in v]
            from collections import Counter

            lines.append(
                f"[{var}] S2@{L} n={len(v)}: median cost best basin(3) cand "
                f"{np.median(a[:, 0]):.3f} vs first-ranked {np.median(a[:, 1]):.3f} "
                f"(basin cand cost higher in {np.mean(a[:, 0] > a[:, 1]):.2f}); median error basin "
                f"cand {np.median(a[:, 2]):.2f} deg, first-ranked {np.median(a[:, 3]):.1f} deg; "
                f"cost at the exact truth {np.median(a[:, 4]):.3f}; first-ranked is CSL: "
                f"{dict(Counter(labs).most_common(5))}"
            )
    (C.OUT_DIR / "s2_mechanism.txt").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
