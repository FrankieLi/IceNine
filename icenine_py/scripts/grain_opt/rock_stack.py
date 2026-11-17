"""GO7 (I23 part 1): emulate a rocking curve on a rendered image set.

Each lit pixel's clean intensity is spread over the omega frames that the window
[omega_c - D/2, omega_c + D/2] overlaps, D = alpha / |sin eta| (docs/omega_peak_width_derivation.md;
eta = azimuth of the pixel's lab direction from the sample origin, measured from the rotation
axis z around the beam x), omega_c the frame centre (the render keeps no sub-frame omega), D
capped at ``--cap`` frames; each frame gets the overlap fraction of the intensity. A pixel is lit
where the spread intensity reaches the threshold (the 5th percentile of the clean intensities: a
fixed detector threshold). Frames past either end of the omega range are dropped.

Writes ``stacks/<out>_clean.npy`` (uint8) and ``<out>/intensity_clean.npz`` (spread intensities
of the lit pixels) under PHASE_D in this worktree, and prints counts.

Usage (from icenine_py/):
  uv run python scripts/grain_opt/rock_stack.py --tag grad --alpha 0.3 --out grad_rock03
"""

import argparse
import json
import math

import common as CM  # noqa: I001
import numpy as np

N_OMEGA, N_DET, NROW, NCOL = 180, 2, 2048, 2048


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--tag", default="grad")
    ap.add_argument("--alpha", type=float, required=True, help="Bragg window, deg")
    ap.add_argument("--out", required=True)
    ap.add_argument("--cap", type=float, default=8.0, help="max window, frames")
    a = ap.parse_args()
    from icenine import spot_pairing as sp

    S = CM.Sample(CM.OUT / "go7" / f"setup_{a.out}", start="paint", intensity=False, sample=a.tag)
    geoms = sp.geometry_from_detectors(S.setup.exp_setup.get_detector_list())
    width_deg = math.degrees(float(S.sc.width))
    z = np.load(CM.PHASE_D / a.tag / "intensity_clean.npz")
    offs, rows, cols, vals = z["frame"], z["row"], z["col"], z["val"]
    thr = float(np.percentile(vals, 5))
    stack_path = CM.PHASE_D / "stacks" / f"{a.out}_clean.npy"
    store_dir = CM.PHASE_D / a.out
    store_dir.mkdir(parents=True, exist_ok=True)
    stack = np.lib.format.open_memmap(
        str(stack_path), mode="w+", dtype=np.uint8, shape=(N_OMEGA * N_DET, NROW, NCOL)
    )
    o_offs, o_rows, o_cols, o_vals = [0], [], [], []
    n_lit_clean = int(len(vals))
    n_spread = []
    for d in range(N_DET):
        g = geoms[d]
        contrib = {o: [] for o in range(N_OMEGA)}
        for o in range(N_OMEGA):
            f = d * N_OMEGA + o
            r = rows[offs[f] : offs[f + 1]].astype(np.int64)
            c = cols[offs[f] : offs[f + 1]].astype(np.int64)
            v = vals[offs[f] : offs[f + 1]].astype(np.float64)
            if not len(v):
                continue
            p = sp.pixel_to_lab(g, c + 0.5, r + 0.5)
            s_eta = np.abs(p[:, 1]) / np.maximum(np.hypot(p[:, 1], p[:, 2]), 1e-12)
            D = np.minimum(a.alpha / width_deg / np.maximum(s_eta, 1e-9), a.cap)  # frames
            lo, hi = o + 0.5 - D / 2, o + 0.5 + D / 2
            n_spread.append(float(np.mean(D > 1.0)))
            kmax = int(math.ceil(a.cap / 2)) + 1
            for k in range(-kmax, kmax + 1):
                t = o + k
                if t < 0 or t >= N_OMEGA:
                    continue
                w = np.clip(np.minimum(hi, t + 1) - np.maximum(lo, t), 0.0, None) / D
                m = w > 0
                if m.any():
                    contrib[t].append((r[m], c[m], v[m] * w[m]))
        for t in range(N_OMEGA):
            img = np.zeros((NROW, NCOL), dtype=np.float64)
            for r, c, v in contrib[t]:
                np.add.at(img, (r, c), v)
            lit = img >= thr
            stack[d * N_OMEGA + t] = lit.astype(np.uint8)
            kk, jj = np.nonzero(lit)
            o_rows.append(kk.astype(np.uint16))
            o_cols.append(jj.astype(np.uint16))
            o_vals.append(img[kk, jj].astype(np.float32))
            o_offs.append(o_offs[-1] + len(kk))
            contrib[t] = []
    stack.flush()
    del stack
    np.savez(store_dir / "intensity_clean.npz", frame=np.asarray(o_offs, np.int64),
             row=np.concatenate(o_rows), col=np.concatenate(o_cols), val=np.concatenate(o_vals))
    info = dict(tag=a.tag, alpha_deg=a.alpha, frame_width_deg=width_deg, threshold=thr,
                lit_clean=n_lit_clean, lit_rock=int(o_offs[-1]),
                share_pixels_window_over_1_frame=float(np.mean(n_spread)), stack=str(stack_path))
    (CM.OUT / "go7" / f"rock_{a.out}.json").write_text(json.dumps(info, indent=1) + "\n")
    print(json.dumps(info))


if __name__ == "__main__":
    main()
