"""Task D figures (matplotlib PNG under benchmarks/eval_reduction/D/figures/).

pairing_geometry.png   diagram: two detectors, the ray, the origin on the sample plane, the
                       tolerance terms (side view and top view, true distances; tolerances and
                       the voxel are drawn enlarged where marked)
origin_error_map.png   origin error of the correct pairs over the sample plane, against the ray
                       z-slope, and its distribution
attribution.png        attributed true peaks and false attachments against the radius
sharpness.png          score of wrong orientations against misorientation, pixel vs attributed

uv run python scripts/eval_reduction/pairing_figures.py [geometry|error|attribution|sharpness]
"""

import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.patches import Circle, Polygon  # noqa: E402

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import pairing_common as pc  # noqa: E402

FIG = pc.BENCH / "figures"


def geometry() -> None:
    geo = pc.load_geometry()
    gs = pc.detector_geoms(geo)
    X0, X1 = gs[0].distance_x(), gs[1].distance_x()
    ratio = X1 / X0
    pw = gs[0].pw
    vox = float(geo["voxel_side"])
    om_step = float(np.mean(geo["omega_hi"] - geo["omega_lo"]))
    reg = float(geo["registration_mm"])
    R_s = 0.5  # sample radius drawn (mm)
    fig, ax = plt.subplots(1, 2, figsize=(16, 6.2), gridspec_kw={"width_ratios": [1.15, 1]})

    # ---- side view: x (beam) horizontal, z (rotation axis) vertical
    a = ax[0]
    tth, eta = np.radians(10.0), np.radians(60.0)  # a typical ray (2 theta 10 deg, eta 60 deg)
    bz = np.tan(tth) * np.sin(eta)  # dz/dx
    xo = 0.1
    xs = np.array([-0.7, X1 + 0.7])
    a.plot(xs, [0, 0], color="0.75", lw=1)
    a.add_patch(
        Polygon(
            [[-R_s, -0.012], [R_s, -0.012], [R_s, 0.012], [-R_s, 0.012]],
            color="tab:green",
            alpha=0.7,
        )
    )
    a.text(-0.55, -0.12, "sample layer, z = 0\n(1 mm wide)", fontsize=9, color="tab:green")
    for X, name in ((X0, "detector 0"), (X1, "detector 1")):
        a.plot([X, X], [-0.15, 1.15], color="k", lw=2)
        a.text(X + 0.04, 1.08, f"{name}\nx = {X:.2f} mm", fontsize=9)
    z0, z1 = bz * (X0 - xo), bz * (X1 - xo)
    a.plot([xo, X1], [0, z1], color="tab:red", lw=1.5)
    a.plot([X0, X1], [z0, z1], "o", color="tab:red")
    a.text(X0 - 0.6, z0 + 0.04, "spot centroid P0", color="tab:red", fontsize=9)
    a.text(X1 - 0.65, z1 + 0.04, "spot centroid P1", color="tab:red", fontsize=9)
    a.plot([xo], [0], "s", color="tab:red")
    a.annotate(
        "origin O: the line P0-P1 meets z = 0",
        (xo, 0),
        (0.6, 0.35),
        arrowprops=dict(arrowstyle="->"),
        fontsize=9,
    )
    # exaggerated tolerance wedge (x100): P0, P1 moved by one pixel in z
    k = 100
    for s0, s1 in ((+1, -1), (-1, +1)):
        zz0, zz1 = z0 + s0 * pw * k, z1 + s1 * pw * k
        slope = (zz1 - zz0) / (X1 - X0)
        a.plot(
            [X0, xo - 0.55], [zz0, zz0 + slope * (xo - 0.55 - X0)], color="tab:red", lw=0.8, ls="--"
        )
    a.text(
        -0.68,
        0.30,
        "centroid tolerance of 1 pixel\n(drawn x100): the intercept\n"
        "moves along the ray by\n(1 pixel) / (z-slope)",
        fontsize=8,
        color="tab:red",
    )
    a.set_xlim(-0.75, X1 + 1.5)
    a.set_ylim(-0.2, 1.25)
    a.set_xlabel("x along the beam (mm)")
    a.set_ylabel("z along the rotation axis (mm)")
    a.set_title("Side view (true distances): two spots, one line, one origin on the layer z = 0")

    # ---- top view: x horizontal, y vertical, the sample plane
    b = ax[1]
    b.add_patch(Circle((0, 0), R_s, fill=False, ec="tab:green", lw=1.5))
    b.text(-0.45, 0.46, "sample hull", color="tab:green", fontsize=9)
    b.add_patch(Circle((0, 0), 0.5, fill=False, ec="0.8", lw=0.5))
    for X in (X0, X1):
        b.plot([X, X], [-1.0, 1.0], color="k", lw=2)
    O = np.array([0.18, 0.22])
    ay = 0.12
    p0 = np.array([X0, O[1] + ay * (X0 - O[0])])
    p1 = np.array([X1, O[1] + ay * (X1 - O[0])])
    b.plot([O[0], p1[0]], [O[1], p1[1]], color="tab:red", lw=1.5)
    b.plot([O[0]], [O[1]], "s", color="tab:red")
    b.plot([p0[0], p1[0]], [p0[1], p1[1]], "o", color="tab:red")
    # arc motion within one frame: the origin could be anywhere on this arc
    r = np.linalg.norm(O)
    th = np.arctan2(O[1], O[0])
    tt = th + np.linspace(-om_step / 2, om_step / 2, 20)
    b.plot(r * np.cos(tt), r * np.sin(tt), color="tab:blue", lw=3)
    b.annotate(
        f"arc motion r x d_omega / 2 = {r * om_step / 2 * 1000:.1f} um\n"
        f"(this origin; up to {R_s * om_step / 2 * 1000:.1f} um\n"
        "at the rim): the exact omega inside\nthe 1 degree frame is unknown",
        (r * np.cos(th), r * np.sin(th)),
        (0.35, -0.2),
        arrowprops=dict(arrowstyle="->", color="tab:blue"),
        fontsize=8,
        color="tab:blue",
    )
    b.add_patch(
        Polygon(
            [[-0.6, -0.7], [-0.6 + vox * 15, -0.7], [-0.6 + vox * 7.5, -0.7 + vox * 15 * 0.866]],
            color="tab:orange",
        )
    )
    b.text(
        -0.6,
        -0.78,
        f"one voxel (drawn x15): side {vox * 1000:.2f} um",
        fontsize=8,
        color="tab:orange",
    )
    b.text(
        0.95,
        -0.3,
        "tolerance terms of the pairing\n"
        f"  centroid: 1 pixel on each detector = {pw * 1000:.2f} um\n"
        "  detector 1, plan presets only:\n"
        f"    voxel size x (L1/L0 - 1) = {vox * (ratio - 1) * 1000:.2f} um\n"
        f"    registration = {reg * 1000:.0f} um\n"
        "  hull margin: voxel size\n    + arc motion r x d_omega / 2\n"
        f"L1/L0 = {ratio:.3f} (detector distances)",
        fontsize=8,
        va="top",
    )
    b.set_xlim(-0.75, X1 + 0.6)
    b.set_ylim(-1.0, 1.0)
    b.set_xlabel("x along the beam (mm)")
    b.set_ylabel("y (mm)")
    b.set_title("Top view (y axis stretched): the origin must lie in the hull")
    fig.tight_layout()
    FIG.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIG / "pairing_geometry.png", dpi=130)
    plt.close(fig)


def error_map() -> None:
    z = np.load(pc.CACHE / "origin_error_arrays.npz")
    err, cond, xy = z["err_smp"], z["cond"], z["truth_xy"]
    fig, ax = plt.subplots(1, 3, figsize=(16, 5))
    hb = ax[0].hexbin(
        xy[:, 0],
        xy[:, 1],
        C=err,
        gridsize=24,
        reduce_C_function=np.median,
        cmap="viridis",
        vmin=0,
        vmax=10,
    )
    plt.colorbar(hb, ax=ax[0], label="median origin error (um)")
    ax[0].set_aspect("equal")
    ax[0].set_xlabel("x (mm), sample frame")
    ax[0].set_ylabel("y (mm)")
    ax[0].set_title(
        f"Origin error over the sample plane\n(correct pairs of one-grain spots, n={len(err)})"
    )
    bins = [0.0, 0.02, 0.05, 0.1, 0.2, 1.0]
    data, labels = [], []
    for lo, hi in zip(bins[:-1], bins[1:]):
        m = (cond >= lo) & (cond < hi)
        data.append(err[m])
        labels.append(f"{lo:g} to {hi:g}\n(n={int(m.sum())})")
    ax[1].boxplot(data, labels=labels, whis=(5, 95), showfliers=False)
    ax[1].set_yscale("log")
    ax[1].set_xlabel("ray z-slope |u_z| (unit direction)")
    ax[1].set_ylabel("sample-frame origin error (um)")
    ax[1].set_title("Error against the ray slope (box: quartiles, whiskers 5 to 95%)")
    ax[1].axhline(9.375, color="tab:orange", ls="--", lw=1)
    ax[1].text(0.6, 10.5, "voxel side 9.4 um", color="tab:orange", fontsize=8)
    ax[2].hist(np.clip(err, 0.1, 1000), bins=np.logspace(-1, 3, 50), color="0.4")
    ax[2].set_xscale("log")
    for q in (50, 90, 99):
        v = np.percentile(err, q)
        ax[2].axvline(v, color="tab:red", lw=1)
        ax[2].text(
            v,
            ax[2].get_ylim()[1] * 0.95,
            f"p{q}={v:.1f}",
            rotation=90,
            va="top",
            ha="right",
            fontsize=8,
        )
    ax[2].set_xlabel("sample-frame origin error (um)")
    ax[2].set_ylabel("pairs")
    ax[2].set_title("Distribution of the origin error")
    fig.tight_layout()
    fig.savefig(FIG / "origin_error_map.png", dpi=130)
    plt.close(fig)


def attribution() -> None:
    s = json.loads((pc.BENCH / "step3" / "attach_summary.json").read_text())
    fig, ax = plt.subplots(1, 2, figsize=(12, 4.6))
    for sname, col in (("voxels30", "tab:blue"), ("spread300", "tab:green")):
        d = s["sets"][sname]
        rs = s["radii_um"]
        rec = [d[f"A r={r} um"]["recall"]["rate"] for r in rs]
        fa = [d[f"A r={r} um"]["false_per_voxel"] for r in rs]
        ax[0].plot(rs, rec, "o-", color=col, label=f"{sname} radius A")
        ax[0].axhline(
            d["B footprint"]["recall"]["rate"], color=col, ls=":", label=f"{sname} footprint B"
        )
        ax[1].plot(rs, fa, "o-", color=col, label=f"{sname} radius A")
        ax[1].axhline(
            d["B footprint"]["false_per_voxel"], color=col, ls=":", label=f"{sname} footprint B"
        )
    ax[0].set_xscale("log")
    ax[1].set_xscale("log")
    ax[1].set_yscale("log")
    ax[0].set_xlabel("attribution radius r (um)")
    ax[1].set_xlabel("attribution radius r (um)")
    ax[0].set_ylabel("share of the voxel's true peaks attached")
    ax[1].set_ylabel("false pairs attached per voxel")
    ax[0].set_title("True peaks found")
    ax[1].set_title("False pairs attached")
    ax[0].legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(FIG / "attribution.png", dpi=130)
    plt.close(fig)


def sharpness() -> None:
    s = json.loads((pc.BENCH / "step3" / "separate_summary.json").read_text())
    fig, ax = plt.subplots(1, 1, figsize=(7, 4.6))
    names = {
        "pixel": "pixel overlap (C++ rule)",
        "A_c_0.05_0": "attributed, radius A 15.7 um",
        "A30_c_0.05_0": "attributed, radius A 30 um",
        "B_c_0.05_0": "attributed, footprint B",
        "oracle_c_0.05_0": "attributed, oracle (true pairs)",
    }
    for k, n in names.items():
        d = s["normalised_median_by_angle"][k]
        x = [float(a) for a in d]
        ax.plot(x, list(d.values()), "o-", label=n)
    ax.set_xscale("log")
    ax.axhline(0.5, color="0.7", lw=0.8)
    ax.set_xlabel("misorientation from the truth (degrees, random axis)")
    ax.set_ylabel("score / score at the truth (median of 30 voxels)")
    ax.set_title("Sharpness: direction tolerance 0.05 degrees, same frame")
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(FIG / "sharpness.png", dpi=130)
    plt.close(fig)


if __name__ == "__main__":
    which = sys.argv[1:] or ["geometry", "error", "attribution", "sharpness"]
    for w in which:
        {
            "geometry": geometry,
            "error": error_map,
            "attribution": attribution,
            "sharpness": sharpness,
        }[w]()
        print("wrote", w)
