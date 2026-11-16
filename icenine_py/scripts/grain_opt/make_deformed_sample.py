"""GO3: write a deformed-grain truth map (random voxel-to-voxel misorientation), its grid-only copy
and the record JSON, in the layout of the AR3 gradient sample.

  uv run python scripts/grain_opt/make_deformed_sample.py --tag rand --corr 0
  uv run python scripts/grain_opt/make_deformed_sample.py --tag corr3 --corr 3

Reads the layer's one-orientation-per-grain truth and its voxel-to-grain map (as
make_gradient_sample.py). Levels: one draw per grain, uniformly from LEVELS_DEG (0.5, 1, 2, 5 deg
RMS rotation about the grain orientation), generator seed --seed; the same generator then draws
the fields (``icenine.deformed_sample``). ``--corr`` is the correlation length in voxel sides
(0: white noise). Writes ``benchmarks/grain_opt/samples/<tag>/``: ``truth.mic``, ``grid.mic``,
``sample.json`` (per grain: level, measured RMS and max angle, n voxels, R_centre; neighbour
misorientation quantiles inside grains; sha256 of the truth map). Then render with
  uv run python scripts/phase_d/render_full.py --tag <tag> --workers 4 --variants clean \\
      --mic <icenine_py>/benchmarks/grain_opt/samples/<tag>/truth.mic
and build the stores with build_store.py --tag <tag>.
"""

import argparse
import hashlib
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ICE = HERE.parents[1]
sys.path.insert(0, str(ICE / "scripts" / "autoresearch" / "ar3"))
sys.path.insert(0, str(ICE / "scripts" / "autoresearch" / "ar1"))

import ar3_common as C  # noqa: E402
import numpy as np  # noqa: E402


def main() -> None:
    import misses as M
    from icenine.deformed_sample import (
        LEVELS_DEG,
        neighbour_misorientation_deg,
        random_field_orientations,
    )
    from icenine.mic_file import MicFile
    from icenine.ray_indexing import misorientation_deg, proper_rotations
    from icenine.symmetry import create_cubic_symmetry

    ap = argparse.ArgumentParser()
    ap.add_argument("--tag", required=True)
    ap.add_argument("--corr", type=float, default=0.0, help="correlation length, voxel sides")
    ap.add_argument("--seed", type=int, default=0)
    a = ap.parse_args()
    out_dir = ICE / "benchmarks" / "grain_opt" / "samples" / a.tag
    out_dir.mkdir(parents=True, exist_ok=True)
    truth = MicFile.read(str(M.MIC))
    gmap = np.load(M.GRAINMAP)
    n = len(truth.voxels)
    vxy, _sides = M.voxel_vertices_xy(M.MIC)
    side = float(truth.voxels[0].side_length)
    centres = vxy.mean(axis=1)
    gids = sorted(int(g) for g in np.unique(gmap))
    R_all = np.array([v.orientation for v in truth.voxels], dtype=np.float64)
    R_old = {g: R_all[int(np.nonzero(gmap == g)[0][0])] for g in gids}
    for g in gids:
        assert np.array_equal(R_all[gmap == g], np.broadcast_to(R_old[g], R_all[gmap == g].shape))
    rng = np.random.default_rng(a.seed)
    lv_idx = rng.integers(0, len(LEVELS_DEG), size=len(gids))
    levels = {g: float(LEVELS_DEG[i]) for g, i in zip(gids, lv_idx)}
    R_new = random_field_orientations(centres, gmap, R_old, levels, a.corr * float(side), rng)
    C.write_mic(truth, {v: R_new[v] for v in range(n)}, out_dir / "truth.mic")

    back = MicFile.read(str(out_dir / "truth.mic"))
    Rb = np.array([v.orientation for v in back.voxels], dtype=np.float64)
    sym = proper_rotations(create_cubic_symmetry(3.6))
    ang = np.asarray(misorientation_deg(Rb, np.stack([R_old[int(g)] for g in gmap]), sym), float)
    nb = [[int(u) for u in back.get_neighbors(v, radius=2.0 * back.voxels[v].side_length)]
          for v in range(n)]
    pairs, mis = neighbour_misorientation_deg(Rb, nb, sym)
    same = gmap[pairs[:, 0]] == gmap[pairs[:, 1]]
    grid = MicFile.read(str(M.MIC.parent / "rand_500grains_1mm_neworient_s0_grid.mic"))
    gxy = np.array([v.position[:2] for v in grid.voxels])
    assert np.abs(gxy - np.array([v.position[:2] for v in back.voxels])).max() < 1e-6
    C.write_mic(back, {v: np.zeros((3, 3)) for v in range(n)}, out_dir / "grid.mic")
    sha = hashlib.sha256((out_dir / "truth.mic").read_bytes()).hexdigest()
    nbq = {}
    for L in LEVELS_DEG:
        sel = same & np.array([levels[int(gmap[p])] == L for p in pairs[:, 0]])
        q = np.quantile(mis[sel], [0.5, 0.9, 0.99]) if sel.any() else [np.nan] * 3
        nbq[str(L)] = dict(n_pairs=int(sel.sum()), median=float(q[0]), p90=float(q[1]),
                           p99=float(q[2]), share_over_2deg=float((mis[sel] > 2.0).mean()))
    info = dict(
        tag=a.tag, seed=a.seed, corr_len_voxels=a.corr, corr_len=a.corr * float(side),
        voxel_side=float(side), levels_deg=list(LEVELS_DEG), level_meaning="RMS rotation angle",
        n_voxels=n, n_grains=len(gids), truth_sha256=sha, neighbour_misorientation=nbq,
        grains=[
            dict(id=g, level_deg=levels[g],
                 measured_rms_deg=float(np.sqrt((ang[gmap == g] ** 2).mean())),
                 measured_max_deg=float(ang[gmap == g].max()),
                 n_voxels=int((gmap == g).sum()), R_centre=R_old[g].tolist())
            for g in gids
        ],
        n_by_level={str(L): int(sum(1 for g in gids if levels[g] == L)) for L in LEVELS_DEG},
    )
    (out_dir / "sample.json").write_text(json.dumps(info) + "\n")
    np.save(out_dir / "grainmap.npy", gmap)
    print(json.dumps(dict(tag=a.tag, n_by_level=info["n_by_level"], neighbour=nbq)))


if __name__ == "__main__":
    main()
