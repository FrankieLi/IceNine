"""AR5 step 19, phase A step 4 (gate g0): the flag by arithmetic on the stored maps and pass tries
of the five disclosed samples and of step 18's J5L run of mid (no image, no repair, no fit).

  uv run python scripts/autoresearch/ar5/repair_flags.py [--out FILE]

Applies the frozen ``repair_rules`` (``repair_freeze.py --check`` must pass first: this script
refuses otherwise). Checks the amended prediction 1 (D-20261115-N7): flagged and fits per sample,
the stored J5 wrong counts and totals (``stats_<s>.json``), and that the after-pass map rebuilt from
the stored tries is the stored one. Writes ``benchmarks/autoresearch/ar5/repair/flags.json`` and
exits 1 when anything does not reproduce (the stop rule; nothing is adjusted).
"""

import argparse
import json
import sys
from typing import Any, Dict

import repair_common as RC  # noqa: E402  (first: guards; RESEED_SAMPLE is not needed here)
import repair_data as RD
import repair_rules as RR

import numpy as np


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=str(RC.RD / "flags.json"))
    a = ap.parse_args()
    RC.RF.check_or_exit()
    y = RD.sym()
    rows: Dict[str, Any] = {}
    ok = True
    for s, arm in [(s, "J5") for s in RR.SAMPLES] + [("mid", "J5L")]:
        d = RD.load(s, arm)
        st = RD.stored_stats(s, arm)
        w = RR.wrong_mask(d["R"], d["R_t"], y)
        if arm == "J5L":
            want = (RR.J5L["flagged"], RR.J5L["fits"], RR.J5L["wrong"], RR.J5L["total"])
        else:
            want = (RR.FLAGGED[s], RR.FITS[s], RR.J5_WRONG[s], RR.J5_TOTAL[s])
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
        rows[f"{s}_{arm}" if arm != "J5" else s] = row
        ok = ok and row["ok"]
        print(s, arm, row, flush=True)
    out = dict(g0_ok=ok, rows=rows, amendment="D-20261115-N7")
    with open(a.out, "w") as f:
        json.dump(out, f, indent=1)
        f.write("\n")
    print("g0", "OK" if ok else "FAILED")
    return 0 if ok else 1


if __name__ == "__main__":
    np.seterr(all="warn")
    sys.exit(main())
