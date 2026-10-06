#!/usr/bin/env python3
"""What F1 (Sigma<=29) does not fix: for the runs still wrong afterwards, the E0 stage class, the
CSL class of the original wrong answer, whether any relative's quick-MC result came within 3 deg
of the truth, and the error of the final F1 answer."""

import sys
from collections import Counter
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import common as C  # noqa: E402
import summarize_fixes as SF  # noqa: E402


def main():
    info = dict(np.load(C.OUT_DIR / "voxels.npz"))
    E0 = SF.load_e0(info)
    runs = dict(np.load(C.OUT_DIR / "e0_runs.npz"))
    cls = {
        (int(v), int(vi), int(s)): (c, l)
        for v, vi, s, c, l in zip(
            runs["vidx"], runs["variant"], runs["seed"], runs["class"], runs["csl_label"]
        )
    }
    lines = []
    for var in C.VARIANTS:
        rem = []
        for k, e in E0.items():
            if k[1] != var:
                continue
            fp = C.CACHE_DIR / "f1" / f"v{e['v']}_{var}.npz"
            f1 = np.load(fp)
            s = k[2]
            R, c, extra = SF.f1_answer(f1, s, 29, e, e["R_true"])
            err = float(C.err_deg(R, e["R_true"]))
            if err > 1.0:
                rel_err = C.err_deg(f1[f"s{s}_rel_R_post"], e["R_true"])
                rem.append(
                    (
                        cls[(e["v"], C.VARIANTS.index(var), s)],
                        err,
                        float(rel_err.min()),
                        int(f1[f"s{s}_rel_sigma"][np.argmin(rel_err)]),
                    )
                )
        lines.append(f"[{var}] runs still wrong after F1: {len(rem)}")
        lines.append("  E0 class: " + str(dict(Counter(c[0].split("@")[0] for c, _, _, _ in rem))))
        lines.append(
            "  CSL class of the original wrong answer: "
            + str(dict(Counter(str(c[1]) for c, _, _, _ in rem)))
        )
        near = sum(m < 3.0 for _, _, m, _ in rem)
        lines.append(
            f"  a relative's quick-MC result within 3 deg of the truth in {near}/{len(rem)} of "
            "them (the relatives covered the truth but the cost ranking did not pick it)"
        )
        lines.append(
            "  final errors (deg): "
            + ", ".join(f"{e:.1f}" for _, e, _, _ in sorted(rem, key=lambda t: t[1]))
        )
    (C.OUT_DIR / "f1_remaining.txt").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
