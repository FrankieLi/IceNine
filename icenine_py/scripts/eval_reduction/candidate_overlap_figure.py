"""Task C diagram: the candidates of one level in orientation space.

  uv run python scripts/eval_reduction/candidate_overlap_figure.py --voxel 15901 --rule cpp

Two panels from a stored step-1 run (candidate_overlap.py):
  left   all candidates handed to the quick MC at the chosen level, as symmetry-reduced rotation
         vectors (the member of the symmetry class with the smallest rotation angle, projected on
         its x and y components, degrees)
  right  the candidates near the truth, as the rotation vector of the crystal-frame
         misorientation truth^T R (reduced the same way), so the truth sits at the origin and the
         axes are degrees of misorientation; circles show the level's diameter D and 0.5 D around
         the truth-nearest candidate
Colour = clique (the candidates of one FZ point, cycling through 20 colours); a black ring marks the
kept quarter (the first n_keep of the post-quick-MC ranking, drawn at the discrete-stage
position of the candidate); the red star is the truth.
"""

import argparse
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ICE = HERE.parents[1]
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(ICE / "scripts" / "findoptimal_robustness"))
RAW = ICE / "benchmarks" / "eval_reduction" / "C" / "cache" / "step1_raw"
FIG = ICE / "benchmarks" / "eval_reduction" / "C" / "figures"


def reduced_rotvec_deg(R: np.ndarray) -> np.ndarray:
    """Rotation vector (deg) of the class member R S with the smallest angle."""
    from csl import OPS
    from scipy.spatial.transform import Rotation

    M = R[:, None] @ OPS[None]  # (N, 24, 3, 3)
    rv = Rotation.from_matrix(M.reshape(-1, 3, 3)).as_rotvec().reshape(len(R), 24, 3)
    best = np.linalg.norm(rv, axis=2).argmin(axis=1)
    return np.degrees(rv[np.arange(len(R)), best])


def main() -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    from candidate_overlap import load

    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    ap.add_argument("--voxel", type=int, default=15901)
    ap.add_argument("--rule", default="cpp")
    ap.add_argument("--level", type=int, default=1)
    ap.add_argument("--set", dest="which", default="full")
    ap.add_argument("--zoom-deg", type=float, default=0.0, help="half width of the right panel")
    a = ap.parse_args()
    run = load(RAW / f"{a.which}_{a.rule}_v{a.voxel}.npz")
    L = run["levels"][a.level]
    Rt = np.array(run["meta"]["R_true"])
    D = float(np.degrees(L["diameter"]))
    R = L["R_disc"]
    clique = L["clique"].astype(int)
    nk = int(L["n_keep"])
    kept = np.zeros(len(R), bool)
    kept[L["mc_idx"][:nk]] = True
    _, cid = np.unique(clique, return_inverse=True)
    cmap = plt.get_cmap("tab20")
    col = cmap(cid % 20)

    fig, ax = plt.subplots(1, 2, figsize=(13, 6.2))
    g = reduced_rotvec_deg(R)
    gt = reduced_rotvec_deg(Rt[None])[0]
    for k, (m, ec, lw, s) in enumerate([(~kept, "none", 0.0, 9), (kept, "black", 0.9, 16)]):
        ax[0].scatter(g[m, 0], g[m, 1], c=col[m], s=s, edgecolors=ec, linewidths=lw, zorder=2 + k)
    ax[0].scatter(*gt[:2], marker="*", s=260, c="red", edgecolors="black", zorder=5, label="truth")
    ax[0].set_title(
        f"Level {a.level}: {len(R)} candidates handed to the quick MC, {int(len(np.unique(cid)))} "
        f"cliques,\nkept quarter ringed ({nk}); symmetry-reduced rotation vector, x-y projection"
    )
    ax[0].set_xlabel("rotation vector x (deg)")
    ax[0].set_ylabel("rotation vector y (deg)")
    ax[0].set_aspect("equal")
    ax[0].legend(loc="upper right")

    # local view: misorientation to the truth
    M = np.swapaxes(Rt, -1, -2)[None] @ R
    loc = reduced_rotvec_deg(M)
    dist = np.linalg.norm(loc, axis=1)
    half = a.zoom_deg or 2.5 * D
    near = np.abs(loc[:, 0]) < half
    near &= np.abs(loc[:, 1]) < half
    for k, (m, ec, lw, s) in enumerate([(~kept, "none", 0.0, 26), (kept, "black", 1.2, 44)]):
        mm = m & near
        ax[1].scatter(
            loc[mm, 0], loc[mm, 1], c=col[mm], s=s, edgecolors=ec, linewidths=lw, zorder=2 + k
        )
    ax[1].scatter(0, 0, marker="*", s=300, c="red", edgecolors="black", zorder=6)
    i0 = int(dist.argmin())
    for rad, ls, lab in (
        (D, "-", f"D = {D:.2f} deg"),
        (0.5 * D, "--", f"0.5 D = {0.5 * D:.2f} deg"),
    ):
        ax[1].add_patch(
            plt.Circle(loc[i0, :2], rad, fill=False, ls=ls, color="grey", lw=1.1, label=lab)
        )
    ax[1].set_xlim(-half, half)
    ax[1].set_ylim(-half, half)
    ax[1].set_aspect("equal")
    ax[1].set_title(
        f"Near the truth: {int(near.sum())} candidates; nearest is {dist[i0]:.2f} deg away"
        f"{' (kept)' if kept[i0] else ' (not kept)'}\n"
        "misorientation to the truth (deg); circles around the nearest candidate"
    )
    ax[1].set_xlabel("misorientation x (deg)")
    ax[1].set_ylabel("misorientation y (deg)")
    ax[1].legend(loc="upper right")
    fig.suptitle(
        f"Voxel {a.voxel}, {'C++' if a.rule == 'cpp' else 'lenient'} overlap rule, "
        f"clean image: colours = cliques (candidates of one FZ point), black ring = kept quarter",
        fontsize=11,
    )
    fig.tight_layout()
    FIG.mkdir(parents=True, exist_ok=True)
    fn = FIG / f"candidates_level{a.level}_v{a.voxel}_{a.rule}.png"
    fig.savefig(fn, dpi=130)
    stats = dict(
        voxel=a.voxel, rule=a.rule, level=a.level, n_candidates=int(len(R)),
        n_cliques=int(len(np.unique(cid))), n_keep=nk, diameter_deg=D,
        nearest_to_truth_deg=float(dist[i0]), nearest_is_kept=bool(kept[i0]),
        n_near_panel=int(near.sum()),
    )  # fmt: skip
    (FIG / (fn.stem + ".json")).write_text(json.dumps(stats, indent=1) + "\n")
    print(fn, stats)


if __name__ == "__main__":
    main()
