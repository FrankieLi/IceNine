"""AR5 step 1, part 3a: the hit-ratio flag on the two pick samples, and the pick.

  uv run python scripts/autoresearch/ar5/stop_part3_pick.py

For thresholds 0.50 to 0.95 in steps of 0.05: voxels under the threshold, the wrong ones among them,
the wrong ones above it; also by source and for boundary voxels alone. The pick (plan): the lowest
threshold that flags at least 80% of the wrong voxels pooled over ``big`` and ``grad``. Writes
``part3_pick.json``; committed before any script reads ``mid`` for Part 3.
"""

import json
from typing import Any, Dict

import stop_part3_lib as L
import stop_rules as SR

PICK_SAMPLES = ("big", "grad")


def main() -> None:
    M = {s: L.load(s) for s in PICK_SAMPLES}
    tab: Dict[str, Any] = {s: L.table(M[s]) for s in PICK_SAMPLES}
    thr = SR.pick_threshold([M[s]["hit"] for s in M], [M[s]["wrong"] for s in M])
    pooled = {}
    for t in SR.THRESHOLDS:
        k = f"{t:.2f}"
        fl = sum(tab[s][k]["all"]["flagged"] for s in M)
        wf = sum(tab[s][k]["all"]["wrong_flagged"] for s in M)
        n = sum(tab[s][k]["all"]["n"] for s in M)
        nw = sum(tab[s][k]["all"]["n_wrong"] for s in M)
        pooled[k] = dict(flagged=fl, wrong_flagged=wf, n=n, n_wrong=nw)
    out = dict(pick_samples=list(PICK_SAMPLES), threshold=thr, pooled=pooled, tables=tab)
    L.P.AR5.joinpath("part3_pick.json").write_text(json.dumps(out, indent=1) + "\n")
    print("picked threshold:", thr)
    for k, v in pooled.items():
        print(k, v)


if __name__ == "__main__":
    main()
