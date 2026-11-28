"""AR5 step 20, phase B (review suggestions): what the repair reached on ``held4`` and ``mid``.

  uv run python scripts/autoresearch/ar5/cmarepair_reach.py

Per sample, from the stored J5 map after the pass, the J5C map and the repair's fits (truth is used
for scoring only; ``held4`` is read after its marker): wrong voxels after the pass that the flag
caught, how many of those were fixed, wrong voxels of J5C that were never flagged; for every fixed
voxel which start gave the adopted end (own end or a candidate) and how far the voxel moved; the
voxels that moved less than 1 degree (the sub-degree refinements: how many, median change in
error); the moves over 1 degree that are not fixes. Writes ``reach.json`` in
``benchmarks/autoresearch/ar5/cmarepair``.
"""

import json
import sys
from pathlib import Path
from typing import Any, Dict

import cmarepair_data as RD
import cmarepair_guard as CG
import cmarepair_rules as CR
import numpy as np

C = CG.cache_dir()


def one(s: str) -> Dict[str, Any]:
    y = RD.sym()
    d = RD.load(s, "J5")
    Rt = d["R_t"]
    Rc = np.load(CG.guarded(C / f"J5C_{s}" / "after.npy"))
    fits = np.load(CG.guarded(C / f"J5C_{s}" / "fits.npz"))
    flagged = np.zeros(len(Rt), bool)
    flagged[d["flagged"]] = True
    e_j5 = np.asarray(CR.misorientation_deg(d["R"], Rt, y), float)
    e_c = np.asarray(CR.misorientation_deg(Rc, Rt, y), float)
    w5, wc = e_j5 > CR.WRONG_DEG, e_c > CR.WRONG_DEG
    moved = np.asarray(CR.misorientation_deg(Rc, d["R"], y), float)
    changed = np.where(moved > 0)[0]
    fixed = np.where(w5 & ~wc)[0]
    # which fit gave each adopted end (the first fit of the voxel whose end equals the adopted one)
    via: Dict[int, Any] = {}
    for v in changed:
        rows = np.where(fits["voxel"] == v)[0]
        hit = [r for r in rows if np.array_equal(fits["end"][r], Rc[v])]
        via[int(v)] = int(fits["start"][hit[0]]) if hit else None
    fx = [dict(voxel=int(v), start=via[int(v)], move_deg=float(moved[v])) for v in fixed]
    sub = [int(v) for v in changed if moved[v] < 1.0]
    over = [int(v) for v in changed if moved[v] >= 1.0 and v not in set(map(int, fixed))]
    return dict(
        sample=s,
        wrong_J5=int(w5.sum()),
        wrong_J5C=int(wc.sum()),
        flagged=int(flagged.sum()),
        wrong_J5_flagged=int((w5 & flagged).sum()),
        fixed=int(len(fixed)),
        fixed_by_candidate=int(sum(1 for f in fx if f["start"] is not None and f["start"] >= 0)),
        fixed_by_own_end=int(sum(1 for f in fx if f["start"] == -1)),
        fixed_move_deg_min=float(min(f["move_deg"] for f in fx)),
        fixed_move_deg_max=float(max(f["move_deg"] for f in fx)),
        fixed_J5_error_deg_min=float(e_j5[fixed].min()),
        fixed_J5_error_deg_max=float(e_j5[fixed].max()),
        fixed_J5C_error_deg_min=float(e_c[fixed].min()),
        fixed_J5C_error_deg_max=float(e_c[fixed].max()),
        wrong_J5C_never_flagged=int((wc & ~flagged).sum()),
        n_changed=int(len(changed)),
        n_sub_degree_moves=len(sub),
        sub_degree_lower_error=int(sum(1 for v in sub if e_c[v] < e_j5[v])),
        sub_degree_median_error_change_deg=float(np.median([e_c[v] - e_j5[v] for v in sub])),
        moves_over_1deg_not_fixes=[
            dict(voxel=v, move_deg=float(moved[v]), err_J5=float(e_j5[v]), err_J5C=float(e_c[v]))
            for v in over
        ],
    )


def main() -> int:
    out = {s: one(s) for s in ("held4", "mid")}
    Path(CG.RD / "reach.json").write_text(json.dumps(out, indent=1) + "\n")
    print(json.dumps(out, indent=1))
    return 0


if __name__ == "__main__":
    sys.exit(main())
