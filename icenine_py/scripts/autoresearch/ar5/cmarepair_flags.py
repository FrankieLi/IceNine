"""AR5 step 20, phase A step 4 (gate g0): step 19's flag by arithmetic on the stored maps and pass
tries of the five disclosed samples and of step 18's J5L run of mid (no image, no repair, no fit).

  uv run python scripts/autoresearch/ar5/cmarepair_flags.py [--out FILE]

Applies step 19's frozen flag (imported through ``cmarepair_rules``; ``cmarepair_freeze.py
--check`` must pass first: this script refuses otherwise). Checks the flagged voxels and fits per
sample (the plan's table), the stored J5 wrong counts and totals (``stats_<s>.json``), and that
the after-pass map rebuilt from the stored tries is the stored one. Writes
``benchmarks/autoresearch/ar5/cmarepair/flags.json`` and exits 1 when anything does not reproduce
(the stop rule; nothing is adjusted).
"""

import argparse
import json
import sys
from typing import Any, Dict

import cmarepair_common as RC  # noqa: E402  (first: guards; RESEED_SAMPLE is not needed here)
import cmarepair_data as RD  # noqa: E402
import cmarepair_rules as CR  # noqa: E402
import numpy as np  # noqa: E402


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=str(RC.RD / "flags.json"))
    a = ap.parse_args()
    RC.CF.check_or_exit()
    y = RD.sym()
    rows: Dict[str, Any] = {}
    ok = True
    for s, arm in [(s, "J5") for s in CR.SAMPLES] + [("mid", "J5L")]:
        d = RD.load(s, arm)
        st = RD.stored_stats(s, arm)
        w = CR.wrong_mask(d["R"], d["R_t"], y)
        if arm == "J5L":
            want = (CR.J5L["flagged"], CR.J5L["fits"], CR.J5L["wrong"], CR.J5L["total"])
        else:
            want = (CR.FLAGGED[s], CR.FITS[s], CR.J5_WRONG[s], CR.J5_TOTAL[s])
        got = (len(d["flagged"]), d["fits"], int(w.sum()), st["total"])
        row = dict(
            flagged=got[0],
            fits=got[1],
            wrong=got[2],
            total_with_pass=got[3],
            expected=dict(flagged=want[0], fits=want[1], wrong=want[2], total=want[3]),
            stored_stats_wrong=st["wrong"],
            wrong_caught=int(w[d["flagged"]].sum()),
            after_map_max_abs_diff=d["after_equals_stored"],
            ok=bool(got == want and st["wrong"] == got[2] and d["after_equals_stored"] < 1e-6),
        )
        # the repair run reads the grid, not the truth map: its neighbour lists must be the same
        g = RD.read_grid(s)
        row["grid_nbrs_equal_truth_nbrs"] = bool(g["nbrs"] == d["nbrs"])
        row["ok"] = bool(row["ok"] and row["grid_nbrs_equal_truth_nbrs"])
        rows[f"{s}_{arm}" if arm != "J5" else s] = row
        ok = ok and row["ok"]
        print(s, arm, row, flush=True)
    out = dict(g0_ok=ok, rows=rows)
    with open(a.out, "w") as f:
        json.dump(out, f, indent=1)
        f.write("\n")
    print("g0", "OK" if ok else "FAILED")
    return 0 if ok else 1


if __name__ == "__main__":
    np.seterr(all="warn")
    sys.exit(main())
