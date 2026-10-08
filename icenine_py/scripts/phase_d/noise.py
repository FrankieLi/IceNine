"""Detector-noise model for the Phase D "realistic" full-sample images, and the ASCII image I/O.

The noise parameters are the ones of the earlier realistic windows (icenine.orientation_eval.
RealismConfig, variant "noise": p_miss 0.1, p_flip 0.05, p_hot 0.05, p_blob 0.1; overlap is NOT
added here because neighbouring grains now overlap physically). Those were applied per spot
window; here they are applied per connected component ("spot") of a whole detector frame:

  * miss : the whole spot is not recorded (all its pixels, grown pixels, hot pixel and blob)
  * flip : every lit pixel is dropped, and every unlit 4-neighbour of a lit pixel is lit, each
           with probability p_flip (threshold jitter at the spot edge); a lit neighbour is copied
  * hot  : one isolated hot pixel near the spot, uniform in a window x window box centred on the
           spot centroid (unlit pixels only)
  * blob : a spurious 2-4 px (per axis) block, same placement rule
Overlapping spots form one component, so they are missed together; that is the difference to the
per-window version. The random stream of a frame depends only on (seed, detector, frame).
"""

from dataclasses import dataclass
from typing import Tuple

import numpy as np
from scipy import ndimage


@dataclass(frozen=True)
class NoiseParams:
    p_miss: float = 0.1
    p_flip: float = 0.05
    p_hot: float = 0.05
    p_blob: float = 0.1
    window: int = 32  # placement box of hot pixels / blobs around the spot centroid


def add_detector_noise(
    img: np.ndarray, rng: np.random.Generator, p: NoiseParams = NoiseParams()
) -> Tuple[np.ndarray, dict]:
    """Return (noisy copy of img (rows, cols) float32, counts dict). Input is not modified."""
    img = np.asarray(img, dtype=np.float32)
    H, W = img.shape
    lit = img > 0
    lab, n = ndimage.label(lit)
    counts = {
        "n_spots": int(n),
        "n_missed": 0,
        "n_dropped_px": 0,
        "n_grown_px": 0,
        "n_hot": 0,
        "n_blob": 0,
    }
    counts["spot_sizes"] = np.bincount(lab.ravel())[1:].tolist() if n else []
    out = img.copy()
    if n == 0:
        return out, counts

    # flip: drop lit pixels, grow into unlit 4-neighbours
    pad = np.pad(img, 1)
    padl = np.pad(lab, 1)
    nb_val = np.maximum.reduce([pad[:-2, 1:-1], pad[2:, 1:-1], pad[1:-1, :-2], pad[1:-1, 2:]])
    nb_lab = np.maximum.reduce([padl[:-2, 1:-1], padl[2:, 1:-1], padl[1:-1, :-2], padl[1:-1, 2:]])
    cand = (~lit) & (nb_val > 0)
    cy, cx = np.nonzero(cand)
    grow_sel = rng.random(len(cy)) < p.p_flip
    ly, lx = np.nonzero(lit)
    drop_sel = rng.random(len(ly)) < p.p_flip
    out[ly[drop_sel], lx[drop_sel]] = 0
    gy, gx = cy[grow_sel], cx[grow_sel]
    out[gy, gx] = nb_val[gy, gx]
    lab = lab.copy()
    lab[gy, gx] = nb_lab[gy, gx]
    counts["n_dropped_px"] = int(drop_sel.sum())
    counts["n_grown_px"] = int(grow_sel.sum())

    # per-spot draws (fixed order so the stream is reproducible)
    miss = rng.random(n) < p.p_miss
    has_hot = rng.random(n) < p.p_hot
    hot_off = rng.random((n, 2))
    has_blob = rng.random(n) < p.p_blob
    blob_size = np.minimum((rng.random((n, 2)) * 3).astype(int), 2) + 2
    blob_off = rng.random((n, 2))

    idx = np.arange(1, n + 1)
    centroids = np.array(ndimage.center_of_mass(lit, lab, idx)).reshape(n, 2)
    mean_int = np.asarray(ndimage.mean(img, lab, idx), dtype=np.float32).reshape(n)

    miss_full = np.concatenate([[False], miss])
    out[miss_full[lab]] = 0
    counts["n_missed"] = int(miss.sum())

    half = p.window / 2.0
    for s in np.nonzero(~miss)[0]:
        cyy, cxx = centroids[s]
        if has_hot[s]:
            y = int(cyy - half + hot_off[s, 0] * p.window)
            x = int(cxx - half + hot_off[s, 1] * p.window)
            if 0 <= y < H and 0 <= x < W and out[y, x] == 0:
                out[y, x] = mean_int[s]
                counts["n_hot"] += 1
        if has_blob[s]:
            sy, sx = blob_size[s]
            y = int(cyy - half + blob_off[s, 0] * (p.window - 4))
            x = int(cxx - half + blob_off[s, 1] * (p.window - 4))
            y0, y1, x0, x1 = max(y, 0), min(y + sy, H), max(x, 0), min(x + sx, W)
            if y0 < y1 and x0 < x1:
                blk = out[y0:y1, x0:x1]
                blk[blk == 0] = mean_int[s]
                counts["n_blob"] += 1
    return out, counts


def write_ascii_image(path: str, img: np.ndarray) -> int:
    """Write non-zero pixels as 'j, k, intensity' lines (j = column, k = row), the C++
    CImageData ASCII format (no header). Returns the number of lines."""
    k, j = np.nonzero(img)
    vals = img[k, j]
    with open(path, "w") as f:
        f.write(
            "".join(f"{a}, {b}, {c:g}\n" for a, b, c in zip(j.tolist(), k.tolist(), vals.tolist()))
        )
    return len(k)


def read_ascii_image(path: str, rows: int, cols: int) -> np.ndarray:
    """Inverse of write_ascii_image (also reads the C++ files). Returns float32 (rows, cols)."""
    out = np.zeros((rows, cols), dtype=np.float32)
    data = np.loadtxt(path, delimiter=",", comments="#", ndmin=2) if _nonempty(path) else None
    if data is not None and len(data):
        out[data[:, 1].astype(int), data[:, 0].astype(int)] = data[:, 2]
    return out


def _nonempty(path: str) -> bool:
    with open(path) as f:
        return any(ln.strip() and not ln.startswith("#") for ln in f)
