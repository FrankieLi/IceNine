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
For a miss whose rank 1 is another grain, whether that grain touches the voxel: a voxel of that
grain shares a vertex with it (``touches``), and the smallest vertex-to-vertex distance between
the voxel and that grain in voxel sides (``gap_sides``, 0 when touching). Grain map and vertices
from the simulation's input .mic.

  uv run python scripts/autoresearch/ar1/misses.py [--tag base]
Writes benchmarks/autoresearch/ar1/misses_pilot_<tag>.json and context_pilot.json (voxels in the
.mic, grains in the pilot region, the frame width, the full-search threshold).
"""

import argparse
import json
import sys
from pathlib import Path
from typing import Any, Dict, List, Tuple

import numpy as np
from scipy.spatial import cKDTree

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
MIC = GRAINMAP.parent / "rand_500grains_1mm_neworient_s0.mic"
P = RI.IndexParams()


def voxel_vertices_xy(mic_path: Path) -> Tuple[np.ndarray, np.ndarray]:
    """In-plane vertices (N, 3, 2) of the .mic's triangles (as ``_get_voxel_vertices``) and the
    side length of each voxel (N,)."""
    from icenine.mic_file import MicFile

    m = MicFile.read(str(mic_path))
    xy = np.array([v.position[:2] for v in m.voxels], dtype=np.float64)
    s = np.array([v.side_length for v in m.voxels], dtype=np.float64)
    up = np.array([v.points_up for v in m.voxels], dtype=bool)
    h = s * np.sqrt(3.0) / 2.0
    out = np.empty((len(xy), 3, 2))
    out[:, 0] = xy
    out[:, 1, 0] = np.where(up, xy[:, 0] + s, xy[:, 0] + s / 2.0)
    out[:, 1, 1] = np.where(up, xy[:, 1], xy[:, 1] - h)
    out[:, 2, 0] = np.where(up, xy[:, 0] + s / 2.0, xy[:, 0] + s)
    out[:, 2, 1] = np.where(up, xy[:, 1] + h, xy[:, 1])
    return out, s


def touch(voxel: int, grain: int, gmap: np.ndarray, vxy: np.ndarray, side: np.ndarray) -> Dict:
    """Whether ``grain`` touches ``voxel`` (a voxel of the grain shares a vertex with it, within
    1e-3 of a side), whether one shares an edge (two vertices), and the smallest vertex distance
    between them in voxel sides."""
    gv = vxy[gmap == grain]  # (m, 3, 2)
    d, j = cKDTree(gv.reshape(-1, 2)).query(vxy[voxel], k=1)
    gap = float(d.min() / side[voxel])
    tol = 1e-3 * side[voxel]
    near = np.linalg.norm(gv[:, :, None, :] - vxy[voxel][None, None, :, :], axis=-1) < tol
    shared = near.any(axis=2).sum(axis=1)  # vertices of each grain voxel on the voxel's vertices
    return dict(touches=bool(gap < 1e-3), shares_edge=bool((shared >= 2).any()), gap_sides=gap)


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
    vxy, side = voxel_vertices_xy(MIC)
    assert len(vxy) == len(gmap), "grain map and .mic differ in length"
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
            if match is not None and maj != G:
                row["rank1_grain_touch"] = touch(v, maj, gmap, vxy, side)
        rows.append(row)
    AC.OUT.mkdir(parents=True, exist_ok=True)
    (AC.OUT / f"misses_pilot_{a.tag}.json").write_text(json.dumps(rows) + "\n")
    region = np.load(AC.ICE / "benchmarks" / "phase_d_pilot" / "region_index.npy")
    sc0 = AC.env("full", sorted(stored)[0])["sc"]
    context = dict(
        mic_voxels=int(len(gmap)),
        region_voxels=int(len(region)),
        region_grains=int(len(np.unique(gmap[region]))),
        full_search_min_evals=int(AC.FULL_SEARCH_MIN_EVALS),
        frame_width_deg=float(np.degrees(np.median(np.diff(RI.frame_centre_omegas(sc0))))),
        voxel_side_mm=float(np.median(side)),
    )
    (AC.OUT / "context_pilot.json").write_text(json.dumps(context, indent=1) + "\n")
    miss = [r for r in rows if "reason" in r]
    print("n", len(rows), "rank1 misses", len(miss))
    for k in ("too_few_rays", "neighbour", "no_grain"):
        print(k, sum(r["reason"] == k for r in miss))
    print(
        "rank-1 grain touches the voxel",
        sum(r.get("rank1_grain_touch", {}).get("touches", False) for r in miss),
    )
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
