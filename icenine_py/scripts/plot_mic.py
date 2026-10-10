"""Plot .mic orientation maps and misorientation maps against a reference.

Python 3 port of the plotting in MicAnalyzer/Python/xdm_toolkit/xdm_plot.py (triangle patches,
colour from the Rodrigues vector), without the hexrd dependency.

  uv run python scripts/plot_mic.py --ref truth.mic [--ref-index idx.npy] [--z Z] \
      --mic a.mic b.mic --labels "A" "B" --out maps.png [--title "..."] [--stats stats.json]

Top row: orientation maps (the reference, then each --mic). Orientation colour: the cubic
fundamental-zone Rodrigues vector, each component mapped linearly from [-(sqrt2-1), sqrt2-1] to
[0, 1] as (R, G, B). The scale is fixed, so the same orientation has the same colour in every panel
(xdm_plot renormalised each map separately).
Bottom row: symmetry-reduced misorientation (degrees) of each --mic against the reference, voxel by
voxel, on a blue-to-red scale from 0 to the cubic maximum, 62.8 degrees. Reference voxels with no
matching --mic voxel are drawn grey.

Matching: if a --mic has the same voxels in the same order as the reference (positions within
1e-4), they are paired by index. Otherwise each --mic voxel is paired with the reference voxel whose
triangle centroid is nearest, within 0.3 of the smaller triangle side; this needs the same grid.
--ref-index picks reference voxels (in that order), e.g. a region of a full-sample truth.
--z picks one layer; a file with several z values and no --z is an error.
"""

import argparse
import json
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.collections import PolyCollection  # noqa: E402
from scipy.spatial import cKDTree  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from icenine.mic_file import MicFile, Voxel  # noqa: E402
from icenine.orientation_search import get_symmetry_quaternions  # noqa: E402
from icenine.sampling import matrix_to_quaternion  # noqa: E402
from icenine.symmetry import create_cubic_symmetry  # noqa: E402

RF_MAX = np.sqrt(2.0) - 1.0  # cubic fundamental zone bound on each Rodrigues component
CUBIC_MAX_MISORIENTATION_DEG = 62.7994  # largest disorientation for cubic-cubic
WRONG_DEG = 1.0
UNMATCHED_GREY = (0.82, 0.82, 0.82)


def qmul(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    """Hamilton product of broadcastable (..., 4) quaternions [w, x, y, z]."""
    aw, ax, ay, az = np.moveaxis(a, -1, 0)
    bw, bx, by, bz = np.moveaxis(b, -1, 0)
    return np.stack(
        [
            aw * bw - ax * bx - ay * by - az * bz,
            aw * bx + ax * bw + ay * bz - az * by,
            aw * by - ax * bz + ay * bw + az * bx,
            aw * bz + ax * by - ay * bx + az * bw,
        ],
        axis=-1,
    )


def fz_rodrigues(q: np.ndarray, sym: np.ndarray) -> np.ndarray:
    """Rodrigues vector (N, 3) of the symmetry-equivalent with the smallest rotation angle."""
    eq = qmul(q[:, None, :], sym[None, :, :])  # (N, S, 4): R * S, crystal symmetry on the right
    best = np.abs(eq[..., 0]).argmax(axis=1)
    qf = eq[np.arange(len(q)), best]
    qf = qf * np.sign(qf[:, :1] + (qf[:, :1] == 0))
    return qf[:, 1:] / qf[:, :1]


def orientation_colors(q: np.ndarray, sym: np.ndarray) -> np.ndarray:
    rf = fz_rodrigues(q, sym)
    return np.clip((rf + RF_MAX) / (2.0 * RF_MAX), 0.0, 1.0)


def misorientation_deg(qa: np.ndarray, qb: np.ndarray, sym: np.ndarray) -> np.ndarray:
    """Symmetry-reduced misorientation (degrees) of row i of qa with row i of qb."""
    rel = qmul(qa * np.array([1.0, -1.0, -1.0, -1.0]), qb)
    w = np.abs(qmul(rel[:, None, :], sym[None, :, :])[..., 0]).max(axis=1)
    return np.degrees(2.0 * np.arccos(np.clip(w, 0.0, 1.0)))


def select_voxels(
    mic: MicFile, idx: Optional[np.ndarray] = None, z: Optional[float] = None, name: str = ""
) -> List[Voxel]:
    """Voxels of one layer, optionally the subset `idx` (in that order) first."""
    vox = mic.voxels if idx is None else [mic.voxels[i] for i in idx]
    zs = np.array([float(v.position[2]) for v in vox])
    if z is None:
        levels = np.unique(np.round(zs, 6))
        if len(levels) > 1:
            raise ValueError(f"{name}: {len(levels)} z layers ({levels[:6]} ...); pass --z")
        return list(vox)
    return [v for v, zz in zip(vox, zs) if abs(zz - z) < 1e-6]


def side_of(v: Voxel, mic: MicFile) -> float:
    return float(v.side_length or mic.initial_side_length / 2**v.generation)


def triangles(vox: Sequence[Voxel], mic: MicFile) -> np.ndarray:
    """(N, 3, 2) triangle vertices; position is the left vertex of the base."""
    out = np.empty((len(vox), 3, 2))
    for k, v in enumerate(vox):
        s = side_of(v, mic)
        x, y = float(v.position[0]), float(v.position[1])
        h = (1.0 if v.points_up else -1.0) * s * np.sqrt(3.0) / 2.0
        out[k] = [[x, y], [x + s, y], [x + s / 2.0, y + h]]
    return out


def quats(vox: Sequence[Voxel]) -> np.ndarray:
    return np.array(
        [matrix_to_quaternion(np.asarray(v.orientation, dtype=np.float64)) for v in vox]
    ).reshape(-1, 4)


def match_to_reference(tri_ref: np.ndarray, tri: np.ndarray, tol: float) -> np.ndarray:
    """For each reference voxel, the index of its matching voxel in `tri`, or -1."""
    if len(tri) == len(tri_ref) and np.abs(tri[:, 0] - tri_ref[:, 0]).max() <= 1e-4:
        return np.arange(len(tri_ref))
    c_ref, c = tri_ref.mean(axis=1), tri.mean(axis=1)
    d, j = cKDTree(c).query(c_ref, distance_upper_bound=tol)
    j = np.where(np.isfinite(d), j, -1)
    # keep one reference voxel per --mic voxel (the nearest)
    order = np.argsort(np.where(j >= 0, d, np.inf), kind="stable")
    seen = set()
    for i in order:
        if j[i] < 0:
            break
        if j[i] in seen:
            j[i] = -1
        else:
            seen.add(int(j[i]))
    return j


def draw(ax: plt.Axes, tri: np.ndarray, lim: np.ndarray, **kw: Any) -> PolyCollection:
    pc = PolyCollection(tri, edgecolors="face", linewidths=0.1, **kw)
    ax.add_collection(pc)
    ax.set_xlim(lim[0, 0], lim[1, 0])
    ax.set_ylim(lim[0, 1], lim[1, 1])
    ax.set_aspect("equal")
    ax.set_xticks([])
    ax.set_yticks([])
    return pc


def plot_maps(
    ref_path: str,
    mic_paths: List[str],
    labels: List[str],
    out: str,
    ref_index: Optional[np.ndarray] = None,
    ref_label: str = "Reference",
    title: str = "",
    z: Optional[float] = None,
) -> List[Dict[str, Any]]:
    """Draw the figure and return one stats dict per --mic."""
    sym = np.asarray(get_symmetry_quaternions(create_cubic_symmetry(1.0)), dtype=np.float64)
    ref = MicFile.read(ref_path)
    rv = select_voxels(ref, ref_index, z, ref_path)
    if not rv:
        raise ValueError(f"{ref_path}: no voxels selected")
    tri_ref = triangles(rv, ref)
    q_ref = quats(rv)
    s_ref = min(side_of(v, ref) for v in rv)
    lim = np.array([tri_ref.reshape(-1, 2).min(axis=0), tri_ref.reshape(-1, 2).max(axis=0)])

    n = len(mic_paths) + 1
    fig, axes = plt.subplots(
        2,
        n + 1,
        figsize=(3.6 * n + 0.8, 7.6),
        squeeze=False,
        gridspec_kw={"width_ratios": [1.0] * n + [0.05], "wspace": 0.08, "hspace": 0.25},
    )
    axes[0, n].axis("off")
    draw(axes[0, 0], tri_ref, lim, facecolors=orientation_colors(q_ref, sym))
    axes[0, 0].set_title(ref_label, fontsize=10)
    axes[1, 0].axis("off")
    axes[1, 0].text(
        0.5,
        0.5,
        "Top row: orientation\n(cubic Rodrigues vector\nas RGB, fixed scale)\n\n"
        "Bottom row: misorientation\nagainst the reference (deg)\ngrey = no matching voxel",
        ha="center",
        va="center",
        fontsize=9,
        transform=axes[1, 0].transAxes,
    )

    stats: List[Dict[str, Any]] = []
    mappable = None
    for j, (path, label) in enumerate(zip(mic_paths, labels), start=1):
        mic = MicFile.read(path)
        mv = select_voxels(mic, None, z, path)
        tri = triangles(mv, mic)
        q = quats(mv)
        tol = 0.3 * min([s_ref] + [side_of(v, mic) for v in mv])
        m = match_to_reference(tri_ref, tri, tol)
        ok = m >= 0
        err = np.full(len(rv), np.nan)
        if ok.any():
            err[ok] = misorientation_deg(q[m[ok]], q_ref[ok], sym)
        if len(mv):
            draw(axes[0, j], tri, lim, facecolors=orientation_colors(q, sym))
        else:
            draw(axes[0, j], tri_ref[:0], lim)
        axes[0, j].set_title(label, fontsize=10)
        if (~ok).any():
            draw(axes[1, j], tri_ref[~ok], lim, facecolors=[UNMATCHED_GREY])
        mappable = draw(
            axes[1, j],
            tri_ref[ok],
            lim,
            array=err[ok],
            cmap="jet",
            clim=(0.0, CUBIC_MAX_MISORIENTATION_DEG),
        )
        e = err[ok]
        st: Dict[str, Any] = {
            "label": label,
            "path": str(path),
            "n_reference": int(len(rv)),
            "n_voxels": int(len(mv)),
            "n_matched": int(ok.sum()),
            "n_over_1deg": int((e > WRONG_DEG).sum()),
            "median_deg": float(np.median(e)) if len(e) else None,
            "median_right_deg": (
                float(np.median(e[e <= WRONG_DEG])) if (e <= WRONG_DEG).any() else None
            ),
            "p90_deg": float(np.percentile(e, 90)) if len(e) else None,
            "max_deg": float(e.max()) if len(e) else None,
        }
        stats.append(st)
        sub = f"over 1 deg: {st['n_over_1deg']} of {st['n_matched']} voxels"
        if st["n_matched"] < len(rv):
            sub += f" ({len(rv) - st['n_matched']} unmatched)"
        if st["median_deg"] is not None:
            sub += f"\nmedian {st['median_deg']:.3f} deg"
        axes[1, j].set_title(sub, fontsize=8)
    if mappable is not None:
        cb = fig.colorbar(mappable, cax=axes[1, n])
        cb.set_label("misorientation (deg), max 62.8", fontsize=9)
    else:
        axes[1, n].axis("off")
    if title:
        fig.suptitle(title, fontsize=12)
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return stats


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    ap.add_argument("--ref", required=True, help="reference .mic (e.g. the truth)")
    ap.add_argument("--ref-index", help=".npy of reference voxel indices to use, in order")
    ap.add_argument("--ref-label", default="Reference")
    ap.add_argument("--z", type=float, help="z of the layer to plot (needed for multi-layer files)")
    ap.add_argument("--mic", nargs="+", required=True, help=".mic files to compare")
    ap.add_argument("--labels", nargs="+", help="panel labels, one per --mic")
    ap.add_argument("--title", default="")
    ap.add_argument("--out", required=True, help="output image (.png, .pdf)")
    ap.add_argument("--stats", help="write per-panel statistics to this JSON file")
    a = ap.parse_args()
    labels = a.labels or [Path(p).parent.name for p in a.mic]
    if len(labels) != len(a.mic):
        ap.error("--labels needs one label per --mic")
    idx = np.load(a.ref_index) if a.ref_index else None
    stats = plot_maps(a.ref, a.mic, labels, a.out, idx, a.ref_label, a.title, a.z)
    if a.stats:
        Path(a.stats).write_text(json.dumps(stats, indent=2))
    for st in stats:
        print(f"{st['label']}: over 1 deg {st['n_over_1deg']} of {st['n_matched']} matched")


if __name__ == "__main__":
    main()
