"""AR1: why the first candidate is not the truth, on the pilot's full-search voxels.

Uses the render provenance of the clean full-sample images (the forward model's ``peak_hook``
records, one entry per drawn spot: voxel, reflection, frame, detector) to name the true grain(s) of
every observed ray, so that each attached ray is one of: own pure (only this voxel's grain drew
it), own mixed (drawn also by another grain), other (not drawn by this grain).

Per voxel: counts of those kinds among the attached rays; the rays the true orientation explains;
the rank-1 candidate: which true grain it matches (misorientation to the grains' true
orientations under 1 degree) and how many attached rays it explains. A miss (rank-1 further than
2 degrees from the truth) gets a reason:
  too_few_rays   the truth explains fewer than 8 attached rays
  neighbour      rank 1 equals the true orientation of another grain, which explains more rays
  no_grain       rank 1 matches no grain of the sample (a coincidence of rays)

  uv run python scripts/autoresearch/ar1/misses.py [--tag base]
Writes benchmarks/autoresearch/ar1/misses_pilot_<tag>.json.
"""

import argparse
import json
import sys
from pathlib import Path
from typing import Any, Dict, List

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import ar1_common as AC  # noqa: E402

from icenine import ray_indexing as RI  # noqa: E402

PROV = AC.ICE.parent / "Examples" / "Example2.ManyGrains" / "ScatteringData_PhaseD" / "full"
GRAINMAP = (
    AC.ICE.parent
    / "Examples"
    / "Example2.ManyGrains"
    / "SimInput"
    / "rand_500grains_1mm_neworient_s0_grainmap.npy"
)
P = RI.IndexParams()


def ray_grains(sc: Any, z: Any, grain_of_voxel: np.ndarray, att: np.ndarray) -> List[set]:
    """Set of true grains that drew the spots of each attached ray (both detectors)."""
    rays = sc.observed_rays()
    key = z["comp_frame"].astype(np.int64) * 2 + z["comp_det"].astype(np.int64)
    order = np.argsort(key, kind="stable")
    ks = key[order]
    out: List[set] = []
    for i in att:
        gs: set = set()
        for det, box in ((0, rays.box0[i]), (1, rays.box1[i])):
            k = int(rays.frame[i]) * 2 + det
            lo, hi = np.searchsorted(ks, k), np.searchsorted(ks, k, side="right")
            comp = order[lo:hi]
            cx, cy = z["comp_c_col"][comp], z["comp_c_row"][comp]
            ok = (
                (cx >= box[0] - 0.5)
                & (cx <= box[1] + 0.5)
                & (cy >= box[2] - 0.5)
                & (cy <= box[3] + 0.5)
            )
            for c in comp[ok]:
                vox = z["_tri_voxel"][z["_tri_start"][c] : z["_tri_start"][c + 1]]
                gs.update(int(x) for x in np.unique(grain_of_voxel[vox]))
        out.append(gs)
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--tag", default="base")
    a = ap.parse_args()
    raw = np.load(PROV / "provenance" / "provenance.npz")
    z: Dict[str, Any] = {k: raw[k] for k in ("comp_frame", "comp_det", "comp_c_col", "comp_c_row")}
    tri_comp, tri_vox = raw["tri_comp"], raw["tri_voxel"]
    o = np.argsort(tri_comp, kind="stable")
    z["_tri_voxel"] = tri_vox[o]
    z["_tri_start"] = np.concatenate(
        [[0], np.cumsum(np.bincount(tri_comp, minlength=len(z["comp_frame"])))]
    )
    gmap = np.load(GRAINMAP)
    stored = AC.pilot_stored()
    cand = np.load(AC.OUT / f"cand_pilot_{a.tag}.npz")
    AC.init("full")
    ctx = AC.bench().RUN._ENV["ctx"]
    first_vox: Dict[int, int] = {}
    for v, gid in enumerate(gmap):
        first_vox.setdefault(int(gid), v)
    rows: List[Dict[str, Any]] = []
    for v in sorted(stored):
        e = AC.env("full", v)
        sc = e["sc"]
        crystal = RI.crystal_from_smooth_cost(sc, AC.symmetry_of(e), e["phase"], P.q_max)
        sym = crystal.sym
        att = sc.attached_rays(e["verts"])
        g = RI.scattering_vectors(sc, att)
        gs = ray_grains(sc, z, gmap, att)
        G = int(gmap[v])
        own = np.array([G in s for s in gs])
        pure = np.array([s == {G} for s in gs])
        none = np.array([len(s) == 0 for s in gs])
        Rt = e["R_true"]
        te = RI.explained_mask(Rt, g, crystal, P.ang_tol_deg, P.len_tol)
        R = cand[f"v{v}"]
        d = RI.misorientation_deg(R, Rt[None], sym)
        row: Dict[str, Any] = dict(
            voxel=v,
            grain=G,
            n_att=int(len(att)),
            n_own_pure=int(pure.sum()),
            n_own_mixed=int((own & ~pure).sum()),
            n_other=int((~own).sum()),
            n_no_provenance=int(none.sum()),
            n_true_explained=int(te.sum()),
            true_explained_pure=int((te & pure).sum()),
            true_explained_mixed=int((te & own & ~pure).sum()),
            true_explained_other=int((te & ~own).sum()),
            rank1_dist=float(d[0]),
            rank_of_truth=int(np.argmax(d < 2.0)) + 1 if (d < 2.0).any() else None,
        )
        if d[0] >= 2.0:
            e1 = RI.explained_mask(R[0], g, crystal, P.ang_tol_deg, P.len_tol)
            count: Dict[int, int] = {}
            for i in np.nonzero(e1)[0]:
                if len(gs[i]) == 1:
                    gid = next(iter(gs[i]))
                    count[gid] = count.get(gid, 0) + 1
            maj = max(count, key=count.get) if count else None
            match = None
            if maj is not None:
                Rg = ctx.truth(first_vox[maj])
                dm = float(RI.misorientation_deg(R[0][None], Rg[None], sym)[0])
                match = dict(grain=maj, dist=dm, pure_rays=count[maj])
            # any grain of the sample within 1 degree of rank 1 (not only the majority one)
            row["rank1_explained"] = int(e1.sum())
            row["rank1_majority_grain"] = match
            if row["n_true_explained"] < 8:
                row["reason"] = "too_few_rays"
            elif match is not None and match["dist"] < 1.0 and maj != G:
                row["reason"] = "neighbour"
            else:
                row["reason"] = "no_grain"
        rows.append(row)
    AC.OUT.mkdir(parents=True, exist_ok=True)
    (AC.OUT / f"misses_pilot_{a.tag}.json").write_text(json.dumps(rows) + "\n")
    miss = [r for r in rows if "reason" in r]
    print("n", len(rows), "rank1 misses", len(miss))
    for k in ("too_few_rays", "neighbour", "no_grain"):
        print(k, sum(r["reason"] == k for r in miss))
    print(
        "median attached",
        np.median([r["n_att"] for r in rows]),
        "own pure",
        np.median([r["n_own_pure"] for r in rows]),
        "own mixed",
        np.median([r["n_own_mixed"] for r in rows]),
        "other",
        np.median([r["n_other"] for r in rows]),
    )


if __name__ == "__main__":
    main()
