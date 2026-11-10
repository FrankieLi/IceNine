"""AR1 step 1: sanity of the scattering vectors, before any voting.

For a voxel with known true orientation R, take the observed rays attached to it, compute each
ray's scattering vector in the sample frame and measure how many lie near some R h (h: reciprocal
vectors up to Q-max 8), and the angular and length errors of those that do. The errors set every
tolerance of the indexer. A control rotation (the truth turned by 37 degrees about (1, 2, 3)) is
counted the same way, so that matches by chance are visible.

  uv run python scripts/autoresearch/ar1/sanity.py e0 [--part i/n]
  uv run python scripts/autoresearch/ar1/sanity.py pilot
Writes benchmarks/autoresearch/ar1/sanity_<kind>[_p<i>].json (per voxel counts) and .npz (the
angular and relative length errors of the rays within 1 degree and 1% of the truth).
"""

import argparse
import json
import sys
from pathlib import Path
from typing import Any, Dict, List

import numpy as np
from scipy.spatial.transform import Rotation

sys.path.insert(0, str(Path(__file__).resolve().parent))
import ar1_common as AC  # noqa: E402

from icenine import ray_indexing as RI  # noqa: E402

LEN_TOL = 0.01
THRESH = (0.25, 0.5, 1.0, 2.0)
CONTROL = Rotation.from_rotvec(np.radians(37.0) * np.array([1.0, 2.0, 3.0]) / np.sqrt(14.0))


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("kind", choices=["e0", "pilot"])
    ap.add_argument("--n", type=int, default=0)
    ap.add_argument("--part", default="0/1")
    a = ap.parse_args()
    stored = AC.e0_stored() if a.kind == "e0" else AC.pilot_stored()
    kind = "e0" if a.kind == "e0" else "full"
    voxels = sorted(stored)
    if a.n:
        voxels = voxels[: a.n]
    pi, pn = (int(x) for x in a.part.split("/"))
    voxels = voxels[pi::pn]
    AC.init(kind)
    rows: List[Dict[str, Any]] = []
    all_ang: List[np.ndarray] = []
    all_len: List[np.ndarray] = []
    for v in voxels:
        e = AC.env(kind, v)
        sc = e["sc"]
        crystal = RI.crystal_from_smooth_cost(sc, AC.symmetry_of(e), e["phase"], 8.0)
        att = sc.attached_rays(e["verts"])
        g = RI.scattering_vectors(sc, att)
        if len(g) == 0:
            rows.append(dict(voxel=v, n_att=0))
            continue
        Rt = e["R_true"]
        _, ang, rel = RI.match_rays(Rt, g, crystal, 1.0, LEN_TOL)
        _, angc, relc = RI.match_rays(CONTROL.as_matrix() @ Rt, g, crystal, 1.0, LEN_TOL)
        ok = (ang <= 1.0) & (rel <= LEN_TOL)
        all_ang.append(ang[ok])
        all_len.append(rel[ok])
        rows.append(
            dict(
                voxel=v,
                n_att=int(len(g)),
                truth={f"{t:g}": int(((ang <= t) & (rel <= LEN_TOL)).sum()) for t in THRESH},
                control={f"{t:g}": int(((angc <= t) & (relc <= LEN_TOL)).sum()) for t in THRESH},
            )
        )
        print(v, rows[-1]["n_att"], rows[-1]["truth"], rows[-1]["control"], flush=True)
    AC.OUT.mkdir(parents=True, exist_ok=True)
    suffix = f"_p{pi}" if pn > 1 else ""
    (AC.OUT / f"sanity_{a.kind}{suffix}.json").write_text(json.dumps(rows) + "\n")
    np.savez(
        AC.OUT / f"sanity_{a.kind}{suffix}.npz",
        ang=np.concatenate(all_ang),
        rel=np.concatenate(all_len),
    )


if __name__ == "__main__":
    main()
