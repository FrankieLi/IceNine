"""Task B: three diagrams of the geometry behind the dominant edge types, drawn from real traces
(no schematic numbers): a spot crossing a frame boundary, a vertex crossing a pixel boundary, and a
spot at the detector edge (the three real peaks rejected at the truth).

  uv run python scripts/eval_reduction/cost_edges_diagrams.py --out
  benchmarks/eval_reduction/B/figures
"""

import argparse
import json
import sys
from pathlib import Path
from types import SimpleNamespace
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.patches import Polygon, Rectangle  # noqa: E402

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import cost_edges as CE  # noqa: E402
import cost_edges_run as RUN  # noqa: E402
import parity_common as PC  # noqa: E402


def frame_diagram(ctx: Any, base: SimpleNamespace, out: Path) -> dict:
    v = 15901
    verts, ph, Rt = ctx.verts(v), int(ctx.mic.voxels[v].phase), ctx.truth(v)
    d = CE.unit_directions(RUN.N_DIR, 1000 + v)[0]
    radii = np.arange(0.0, 0.601, 0.005)
    cf = CE.make_cf(base, "full", "cpp")
    trs = [
        CE.trace_eval(cf, R.astype(np.float32), verts, ph) for R in CE.rotate_about(Rt, radii, d)
    ]
    om = np.degrees(np.stack([t["omega"] for t in trs]))  # (T, P)
    fr = np.stack([t["frame"] for t in trs])
    pres = np.stack([t["status"] == CE.PRESENT for t in trs])
    rng = base.range_map
    low, width = np.degrees(rng.low), np.degrees(rng.width)
    moved = [k for k in range(om.shape[1]) if pres[:, k].all() and len(set(fr[:, k])) > 1]
    pick = moved[:5]
    fig, axs = plt.subplots(1, len(pick), figsize=(3.3 * len(pick), 3.9), sharex=True)
    for ax, k in zip(np.atleast_1d(axs), pick):
        j = int(np.nonzero(np.diff(fr[:, k]))[0][0])
        edge = low + width * np.round((om[j, k] + om[j + 1, k]) / 2 / width - low / width)
        ax.plot(radii, om[:, k] - edge, lw=1.6)
        ax.axhline(0.0, color="0.4", lw=1.0, ls="--")
        ax.plot([radii[j + 1]], [om[j + 1, k] - edge], "ko", ms=5)
        ax.set_title(f"peak {k}: image {fr[j, k]} then {fr[j + 1, k]}", fontsize=9)
        ax.set_xlabel("rotation from the truth (deg)")
    np.atleast_1d(axs)[0].set_ylabel("omega minus the frame edge (degrees)")
    fig.suptitle(
        "A frame crossing: the predicted spot's omega drifts across a frame edge (dashed) as the "
        "orientation rotates.\n"
        f"At the dot the peak is compared with the neighbouring image (frames are {width:.0f} "
        f"degree wide). Voxel {v}, one direction, five of {len(moved)} peaks that cross within "
        f"0.6 degree.",
        fontsize=8.5,
    )
    fig.tight_layout()
    fig.savefig(out / "diagram_frame_crossing.png", dpi=150)
    plt.close(fig)
    return dict(
        voxel=v,
        n_peaks_that_change_frame=len(moved),
        n_present_all_steps=int(pres.all(0).sum()),
        shown=len(pick),
    )


def draw_tri(ax: Any, tri: np.ndarray, color: str, label: str) -> None:
    ax.add_patch(Polygon(tri, closed=True, fill=False, ec=color, lw=1.6, label=label))
    ax.plot(tri[:, 0], tri[:, 1], "o", color=color, ms=3)


def pixel_diagram(ctx: Any, base: SimpleNamespace, out: Path) -> dict:
    v = 15901
    verts, ph, Rt = ctx.verts(v), int(ctx.mic.voxels[v].phase), ctx.truth(v)
    d = CE.unit_directions(RUN.N_DIR, 1000 + v)[0]
    cf = CE.make_cf(base, "full", "cpp")
    im0 = base.exp_data.get_image(0, 0)
    dims = (int(im0.num_rows), int(im0.num_cols))
    radii = np.arange(0.0, 0.2001, 0.0005)
    Rs = CE.rotate_about(Rt, radii, d)
    prev = CE.trace_eval(cf, Rs[0].astype(np.float32), verts, ph)
    best = None
    for t in range(1, len(radii)):
        cur = CE.trace_eval(cf, Rs[t].astype(np.float32), verts, ph)
        ev = CE.classify(prev, cur, dims, 0)
        for key in np.nonzero((ev == CE.EVENT_TYPES.index("pixel")) & (cur["q"] != prev["q"]))[0]:
            tri = PC.python_pixels(ctx, Rs[t], verts, int(key) // 2, int(key) % 2)[0]
            area = 0.5 * abs(np.cross(tri[1] - tri[0], tri[2] - tri[0]))
            if best is None or area > best[0]:
                best = (area, t, int(key), prev, cur)
        prev = cur
    assert best is not None
    _, t, key, a, b = best
    refl, branch = key // 2, key % 2
    frame = int(b["frame"][key])
    pa = PC.python_pixels(ctx, Rs[t - 1], verts, refl, branch)[0]  # detector 0, (3, 2)
    pb = PC.python_pixels(ctx, Rs[t], verts, refl, branch)[0]
    img = base.exp_data.get_image(frame, 0).get_binary_numpy()
    c0 = pa.mean(0)
    x0, x1 = int(c0[0]) - 7, int(c0[0]) + 8
    y0, y1 = int(c0[1]) - 7, int(c0[1]) + 8
    fig, axs = plt.subplots(1, 2, figsize=(9, 4.4))
    for ax, tri, tag, tr in (
        (axs[0], pa, f"{radii[t - 1]:.4f} deg", a),
        (axs[1], pb, f"{radii[t]:.4f} deg", b),
    ):
        ax.imshow(
            img[y0:y1, x0:x1],
            extent=(x0, x1, y1, y0),
            cmap="Greys",
            vmin=0,
            vmax=1.4,
            interpolation="nearest",
        )
        ax.set_xticks(np.arange(x0, x1 + 1, 2))
        ax.set_yticks(np.arange(y0, y1 + 1, 2))
        ax.set_xlim(x0, x1)
        ax.set_ylim(y1, y0)
        ax.grid(True, color="0.75", lw=0.4)
        ax.add_patch(Polygon(tri, closed=True, fill=False, ec="tab:red", lw=1.8))
        ax.plot(tri[:, 0], tri[:, 1], "o", color="tab:red", ms=3)
        ax.set_title(
            f"rotation {tag}: pixels lit by the triangle {int(tr['cnt'][key, 0])} of "
            f"{int(tr['cnt'][key, 1])}",
            fontsize=9,
        )
        ax.set_xlabel("detector column (pixel)")
    axs[0].set_ylabel("detector row (pixel)")
    fig.suptitle(
        f"A pixel crossing: voxel {v}, detector 0, frame {frame}, rotation about {radii[t]:.3f} "
        f"deg. Grey = lit pixels of the image, red = predicted triangle.\n"
        f"The cost of this peak changes from {a['q'][key]:.3f} to {b['q'][key]:.3f} when a "
        f"vertex moves into another pixel.",
        fontsize=9,
    )
    fig.tight_layout()
    fig.savefig(out / "diagram_pixel_crossing.png", dpi=150)
    plt.close(fig)
    return dict(
        voxel=v,
        radius_before=float(radii[t - 1]),
        radius_after=float(radii[t]),
        q_before=float(a["q"][key]),
        q_after=float(b["q"][key]),
        frame=frame,
    )


def detector_diagram(ctx: Any, base: SimpleNamespace, out: Path) -> dict:
    summ = json.loads((out.parent / "step3" / "summary.json").read_text())
    peaks = summ["peaks"]
    cfs = CE.make_cf(base, "screen", "cpp")
    cff = CE.make_cf(base, "full", "cpp")
    ph0 = int(ctx.mic.voxels[peaks[0]["voxel"]].phase)
    gs = cfs._phase_recip_vecs[ph0][0].numpy()
    gf = cff._phase_recip_vecs[ph0][0].numpy()
    fig, axs = plt.subplots(1, len(peaks), figsize=(4.3 * len(peaks), 4.4))
    for ax, p in zip(np.atleast_1d(axs), peaks):
        v = p["voxel"]
        fk = int(np.abs(gf - gs[p["refl"]]).sum(1).argmin())
        verts = ctx.verts(v)
        tri = PC.python_pixels(ctx, ctx.truth(v), verts, fk, p["branch"])[1]  # detector 1
        img = base.exp_data.get_image(p["frame"], 1).get_binary_numpy()
        cx, cy = p["screen_detectors"][1]["center"]
        x0, x1 = cx - 14, cx + 15
        y0, y1 = 2048 - 16, 2048 + 6
        sub = np.zeros((y1 - y0, x1 - x0))
        sub[: 2048 - y0, :] = img[y0:2048, max(x0, 0) : x1]
        ax.imshow(
            sub, extent=(x0, x1, y1, y0), cmap="Greys", vmin=0, vmax=1.4, interpolation="nearest"
        )
        ax.add_patch(Rectangle((x0, 2048), x1 - x0, y1 - 2048, color="tab:blue", alpha=0.18, lw=0))
        ax.axhline(2048, color="tab:blue", lw=1.5)
        ax.add_patch(Polygon(tri, closed=True, fill=False, ec="tab:red", lw=1.8))
        ax.add_patch(Rectangle((cx - 3, cy - 3), 7, 7, fill=False, ec="tab:green", lw=1.4, ls="--"))
        ax.set_xlim(x0, x1)
        ax.set_ylim(y1, y0)
        ax.set_title(f"voxel {v}, frame {p['frame']}: centre row {cy}", fontsize=9)
        ax.set_xlabel("detector 1 column (pixel)")
        ax.set_yticks([2034, 2038, 2042, 2046, 2050])
    np.atleast_1d(axs)[0].set_ylabel("detector 1 row (pixel)")
    fig.suptitle(
        "Detector edge: the three real peaks the C++ rule rejects at the true orientation. Blue "
        "= beyond the last row (2047),\n"
        "red = predicted triangle (entirely beyond the edge), green dashed = the screen's 7 x 7 "
        "search box, grey = lit pixels (none)",
        fontsize=9,
    )
    fig.tight_layout()
    fig.savefig(out / "diagram_detector_edge.png", dpi=150)
    plt.close(fig)
    return dict(n=len(peaks))


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    out = Path(a.out).resolve()
    out.mkdir(parents=True, exist_ok=True)
    RUN._init("full")
    ctx, base = RUN._ENV["ctx"], RUN._ENV["base"]
    info = dict(
        frame=frame_diagram(ctx, base, out),
        pixel=pixel_diagram(ctx, base, out),
        detector=detector_diagram(ctx, base, out),
    )
    (out / "diagrams.json").write_text(json.dumps(info, indent=1) + "\n")
    print(info)


if __name__ == "__main__":
    main()
