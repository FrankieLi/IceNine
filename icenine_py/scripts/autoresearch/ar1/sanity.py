"""AR1 step 1: sanity of the scattering vectors, before any voting.

For a voxel with known true orientation R, take the observed rays attached to it, compute each
ray's scattering vector in the sample frame and measure how many lie near some R h (h: reciprocal
vectors up to Q-max 8), and the angular and length errors of those that do. The errors set every
tolerance of the indexer.

  uv run python scripts/autoresearch/ar1/sanity.py e0 --n 200
  uv run python scripts/autoresearch/ar1/sanity.py pilot
Writes benchmarks/autoresearch/ar1/sanity_<kind>.json (per voxel and pooled).
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

WIDE_DEG, WIDE_LEN = 5.0, 0.05


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("kind", choices=["e0", "pilot"])
    ap.add_argument("--n", type=int, default=0)
    a = ap.parse_args()
    if a.kind == "e0":
        stored = AC.e0_stored()
        voxels = sorted(stored)
        AC.init("e0")
        kind = "e0"
    else:
        stored = AC.pilot_stored()
        voxels = sorted(stored)
        AC.init("full")
        kind = "full"
    if a.n:
        voxels = voxels[: a.n]
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
        j, ang, rel = RI.match_rays(e["R_true"], g, crystal, WIDE_DEG, WIDE_LEN)
        near = (ang <= WIDE_DEG) & (rel <= WIDE_LEN)
        all_ang.append(ang[near])
        all_len.append(rel[near])
        rows.append(
            dict(
                voxel=v,
                n_att=int(len(g)),
                n_within={
                    f"{t:g}": int(((ang <= t) & (rel <= WIDE_LEN)).sum())
                    for t in (0.25, 0.5, 1, 2, 3, 5)
                },
                med_ang=float(np.median(ang[near])) if near.any() else None,
            )
        )
        print(v, rows[-1]["n_att"], rows[-1]["n_within"], flush=True)
    ang = np.concatenate(all_ang)
    rel = np.concatenate(all_len)
    pooled = dict(
        n_voxels=len(rows),
        n_attached=int(sum(r["n_att"] for r in rows)),
        n_matched_wide=int(len(ang)),
        ang_pct={str(q): float(np.percentile(ang, q)) for q in (50, 75, 90, 95, 99)},
        len_pct={str(q): float(np.percentile(rel, q)) for q in (50, 75, 90, 95, 99)},
        share_within={f"{t:g}": float((ang <= t).mean()) for t in (0.25, 0.5, 1, 2, 3, 5)},
    )
    print(json.dumps(pooled, indent=1))
    AC.OUT.mkdir(parents=True, exist_ok=True)
    (AC.OUT / f"sanity_{a.kind}.json").write_text(json.dumps(dict(pooled=pooled, rows=rows)) + "\n")


if __name__ == "__main__":
    main()
