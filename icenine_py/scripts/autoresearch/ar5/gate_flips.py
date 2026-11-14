"""AR5 step 7, Phase C: where do the live pass and the full pass differ?

  uv run python scripts/autoresearch/ar5/gate_flips.py

For each live skipped pass (``gate/cache/live_<sample>_<run>_h0.95/after_live.npy``, not committed)
and the stored full pass: the voxels whose end differs (by more than 1e-6 in the matrix entries),
split into skipped and kept voxels, and the own hit ratio of every voxel that is wrong in one
map and right in the other. Also the key-off parity run on ``mid``: the largest difference
between its map and the stored map after the pass. Truth is used for scoring only. Writes
``gate/flips.json``. Describes; it tests no cause.
"""

import json
from typing import Any, Dict

import gate_common as GC
import gate_rules as GR
import numpy as np

LIVE = [("mid", "S11n"), ("grad", "S11n"), ("big", "S11"), ("held", "S11v2")]


def main() -> None:
    out: Dict[str, Any] = {}
    for s, run in LIVE:
        P = GC.load_pass(s, run)
        live = np.load(GC.GATE / "cache" / f"live_{s}_{run}_h0.95" / "after_live.npy")
        T = GC.Truth(s)
        own = P["old_score"][:, 1]
        old_cost = P["old_score"][:, 0]
        rep = GR.replay(own, old_cost, P["voxel"], P["hit"], P["cost"], P["evals"], GR.H_S)
        sk = rep["skipped"]
        diff = np.abs(live - P["after"]).reshape(GC.N, -1).max(1) > 1e-6
        w_full = T.err(P["after"]) > GC.WRONG
        w_live = T.err(live) > GC.WRONG
        fixed = np.flatnonzero(w_full & ~w_live)
        broken = np.flatnonzero(~w_full & w_live)
        boundary = np.zeros(GC.N, dtype=bool)
        boundary[P["voxel"]] = True
        kept_b = boundary & ~sk
        out[f"{s} {run}"] = dict(
            boundary_kept=int(kept_b.sum()),
            skipped=int(sk.sum()),
            end_differs_skipped=int((diff & sk).sum()),
            end_differs_kept=int((diff & kept_b).sum()),
            end_differs_not_boundary=int((diff & ~boundary).sum()),
            fixed_vs_full=int(len(fixed)),
            broken_vs_full=int(len(broken)),
            fixed_own_hit=[round(float(x), 3) for x in sorted(own[fixed])],
            broken_own_hit=[round(float(x), 3) for x in sorted(own[broken])],
            flips_skipped=int(sk[fixed].sum() + sk[broken].sum()),
            wrong_full_and_live_both=int((w_full & w_live).sum()),
        )
        print(s, run, {k: v for k, v in out[f"{s} {run}"].items() if "own_hit" not in k})
    a = np.load(GC.GATE / "cache" / "live_mid_S11n_h0" / "after_live.npy")
    b = GC.load_pass("mid", "S11n")
    d = np.abs(a - b["after"]).reshape(GC.N, -1).max(1)
    ch = np.flatnonzero((b["after"] != b["before"]).reshape(GC.N, -1).any(1))
    out["parity_map"] = dict(
        max_abs_diff=float(d.max()),
        voxels_differing=int((d > 0).sum()),
        stored_changed=int(len(ch)),
        differing_among_unchanged=int((np.delete(d, ch) > 0).sum()),
    )
    print(out["parity_map"])
    (GC.GATE / "flips.json").write_text(json.dumps(out, indent=1) + "\n")


if __name__ == "__main__":
    main()
