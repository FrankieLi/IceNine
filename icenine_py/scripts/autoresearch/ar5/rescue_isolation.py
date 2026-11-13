"""AR5 step 9, review fixes 5 and 6: which voxels differ between arm A and arm G (draw isolation).

  RESEED_SAMPLE=<held|mid|grad|big> uv run python scripts/autoresearch/ar5/rescue_isolation.py

Reads the stored maps of arm ``A`` (step 8, linked) and arm ``G`` (before and after the boundary
pass) and G's ``rescue_log.json.gz``. Writes ``rescue/isolation_<sample>.json``:

* ``differ_before`` / ``differ_after``: voxels whose orientation differs between A and G (any bit of
  the stored map), the rescued voxels among them and the others; the others in a grain that holds a
  rescued voxel (truth grain map, scoring only) and the largest A-to-G angle among the others;
* ``rescued``: at each rescued voxel, the angle between G's voxel and A's voxel (before the pass):
  median and largest; the largest among the voxels where A was right; the voxels where A was wrong
  (more than 1 degree) with A's and G's errors.

This is the check of draw isolation the review asks for: with aligned draws only the rescued
voxels, and voxels grown from them, may change. ``big`` has no failed seed step, so G = A there
says only that the idle wrapper adds no draws. Angles in degrees (symmetry-reduced). No p-values.
"""

import gzip
import json
from typing import Any, Dict, List

import rescue_common as RCM  # noqa: E402  (first)

import numpy as np  # noqa: E402

import gradient_common as G  # noqa: E402
import reseed_common as RC  # noqa: E402

SAMPLE = RC.SAMPLE
N = 24570


def differ_part(
    d: np.ndarray, same: np.ndarray, rescued: List[int], gmap: np.ndarray
) -> Dict[str, Any]:
    """Counts of the voxels that differ (``~same``), split by rescued and other."""
    diff = set(np.where(~same)[0].tolist())
    rv = set(rescued)
    others = sorted(diff - rv)
    grains = {int(gmap[v]) for v in rv}
    return dict(
        n=len(diff),
        rescued_in=len(diff & rv),
        rescued_not_in=len(rv - diff),
        others=len(others),
        others_in_a_rescued_grain=int(sum(1 for v in others if int(gmap[v]) in grains)),
        others_max_deg=float(max((d[v] for v in others), default=0.0)),
        equals_rescued_set=diff == rv,
        n_over_1e6_deg=int((d > 1e-6).sum()),
    )


def main() -> None:
    C = RCM.RCACHE
    T = G.Truth()
    Ab, _ = G.read_map(C / f"A_{SAMPLE}" / "recon.mic.gz", N)
    Gb, _ = G.read_map(C / f"G_{SAMPLE}" / "recon.mic.gz", N)
    Aa = np.load(C / f"passes_{SAMPLE}" / "after_A.npy")
    Ga = np.load(C / f"passes_{SAMPLE}" / "after_G.npy")
    with gzip.open(C / f"G_{SAMPLE}" / "rescue_log.json.gz", "rt") as f:
        log = json.load(f)
    rv = sorted(int(r["voxel"]) for r in log if r["rescued"])
    db, da = T.mis(Ab, Gb), T.mis(Aa, Ga)
    same_b = np.all(Ab.reshape(N, -1) == Gb.reshape(N, -1), axis=1)
    same_a = np.all(Aa.reshape(N, -1) == Ga.reshape(N, -1), axis=1)
    eA, eG = T.err(Ab), T.err(Gb)
    resc: Dict[str, Any] = dict(n=len(rv))
    if rv:
        dv = np.array([db[v] for v in rv])
        right = [v for v in rv if eA[v] <= G.WRONG]
        resc.update(
            median_deg=float(np.median(dv)),
            max_deg=float(dv.max()),
            max_voxel=int(rv[int(dv.argmax())]),
            a_right=len(right),
            max_deg_where_a_right=float(max(db[v] for v in right)) if right else None,
            a_wrong=[
                dict(voxel=v, err_a=float(eA[v]), err_g=float(eG[v])) for v in rv if eA[v] > G.WRONG
            ],
        )
    out = dict(
        sample=SAMPLE,
        n_voxels=N,
        rescued_voxels=rv,
        differ_before=differ_part(db, same_b, rv, T.gmap),
        differ_after=differ_part(da, same_a, rv, T.gmap),
        rescued=resc,
        note=(
            "voxels whose stored orientation differs between arm A and arm G; with aligned draws "
            "only the rescued voxels and voxels grown from them may change"
        ),
    )
    (RCM.RESCUE / f"isolation_{SAMPLE}.json").write_text(
        json.dumps(out, indent=1, default=float) + "\n"
    )
    print(json.dumps({k: out[k] for k in ("differ_before", "differ_after", "rescued")}, indent=1))


if __name__ == "__main__":
    main()
