"""GO7b (I24): emulate a detector point-spread function on a rendered image set.

Each frame's clean intensities (``<tag>/intensity_clean.npz``) are blurred with a Gaussian PSF
of ``--sigma`` pixels and thresholded at the 5th percentile of the clean intensities (a fixed
detector threshold). Writes ``stacks/<out>_clean.npy`` and ``<out>/intensity_clean.npz`` under
PHASE_D in this worktree.

Usage (from icenine_py/):
  uv run python scripts/grain_opt/psf_stack.py --tag grad --sigma 1.0 --out grad_psf10
"""

import argparse
import json

import common as CM  # noqa: I001
import numpy as np
from scipy import ndimage

N_OMEGA, N_DET, NROW, NCOL = 180, 2, 2048, 2048


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--tag", required=True)
    ap.add_argument("--sigma", type=float, required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    z = np.load(CM.PHASE_D / a.tag / "intensity_clean.npz")
    offs, rows, cols, vals = z["frame"], z["row"], z["col"], z["val"]
    thr = float(np.percentile(vals, 5))
    stack_path = CM.PHASE_D / "stacks" / f"{a.out}_clean.npy"
    (CM.PHASE_D / a.out).mkdir(parents=True, exist_ok=True)
    stack = np.lib.format.open_memmap(
        str(stack_path), mode="w+", dtype=np.uint8, shape=(N_OMEGA * N_DET, NROW, NCOL)
    )
    o_offs, o_rows, o_cols, o_vals = [0], [], [], []
    for f in range(N_OMEGA * N_DET):
        img = np.zeros((NROW, NCOL), dtype=np.float32)
        s = slice(offs[f], offs[f + 1])
        img[rows[s].astype(np.int64), cols[s].astype(np.int64)] = vals[s]
        if offs[f + 1] > offs[f]:
            img = ndimage.gaussian_filter(img, a.sigma, mode="constant", truncate=4.0)
        lit = img >= thr
        stack[f] = lit.astype(np.uint8)
        kk, jj = np.nonzero(lit)
        o_rows.append(kk.astype(np.uint16))
        o_cols.append(jj.astype(np.uint16))
        o_vals.append(img[kk, jj].astype(np.float32))
        o_offs.append(o_offs[-1] + len(kk))
    stack.flush()
    del stack
    np.savez(CM.PHASE_D / a.out / "intensity_clean.npz", frame=np.asarray(o_offs, np.int64),
             row=np.concatenate(o_rows), col=np.concatenate(o_cols), val=np.concatenate(o_vals))
    info = dict(tag=a.tag, sigma_px=a.sigma, threshold=thr, lit_clean=int(len(vals)),
                lit_psf=int(o_offs[-1]), stack=str(stack_path))
    out = CM.OUT / "go7b"
    out.mkdir(parents=True, exist_ok=True)
    (out / f"psf_{a.out}.json").write_text(json.dumps(info, indent=1) + "\n")
    print(json.dumps(info))


if __name__ == "__main__":
    main()
