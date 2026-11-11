"""AR5 step 1, part 3b: the hit-ratio flag on ``mid`` at the picked threshold; prediction 4.

  uv run python scripts/autoresearch/ar5/stop_part3_test.py

Reads ``part3_pick.json`` (committed before ``mid`` was read). Writes ``part3_flags.json``: the
tables of all thresholds for the three samples, the flag at the picked threshold (voxels flagged,
wrong flagged, wrong left unflagged with a Wilson 95% interval: what a user who stops at the flag
leaves behind), and prediction 4 (the flag catches at least 80% of the wrong voxels and flags at
most 5% of all voxels, on the pick samples and on mid). Intervals treat voxels as independent
although they cluster in grains.
"""

import json
from typing import Any, Dict

import stop_part3_lib as L
import stop_rules as SR
from stop_part3_pick import PICK_SAMPLES


def main() -> None:
    pick = json.loads((L.P.AR5 / "part3_pick.json").read_text())
    thr = pick["threshold"]
    M = {s: L.load(s) for s in ("big", "grad", "mid")}
    tables: Dict[str, Any] = {s: L.table(M[s]) for s in M}
    key = f"{thr:.2f}"
    at: Dict[str, Any] = {}
    for s in M:
        a = tables[s][key]["all"]
        un = a["n"] - a["flagged"]
        at[s] = dict(
            **a,
            share_flagged=a["flagged"] / a["n"],
            recall=a["wrong_flagged"] / a["n_wrong"],
            unflagged=un,
            wrong_unflagged_interval=SR.wilson_row(a["wrong_unflagged"], un),
            recall_interval=SR.wilson_row(a["wrong_flagged"], a["n_wrong"]),
        )
    pooled_pick = dict(
        flagged=sum(at[s]["flagged"] for s in PICK_SAMPLES),
        n=sum(at[s]["n"] for s in PICK_SAMPLES),
        wrong_flagged=sum(at[s]["wrong_flagged"] for s in PICK_SAMPLES),
        n_wrong=sum(at[s]["n_wrong"] for s in PICK_SAMPLES),
    )
    pooled_pick["share_flagged"] = pooled_pick["flagged"] / pooled_pick["n"]
    pooled_pick["recall"] = pooled_pick["wrong_flagged"] / pooled_pick["n_wrong"]
    ok = lambda r, sh: bool(r >= 0.8 - SR.TOL and sh <= 0.05 + SR.TOL)  # noqa: E731
    p4 = dict(
        pick_samples=ok(pooled_pick["recall"], pooled_pick["share_flagged"]),
        mid=ok(at["mid"]["recall"], at["mid"]["share_flagged"]),
    )
    p4["holds"] = bool(p4["pick_samples"] and p4["mid"])
    out = dict(
        threshold=thr,
        at_threshold=at,
        pooled_pick=pooled_pick,
        prediction4=p4,
        tables=tables,
    )
    L.P.AR5.joinpath("part3_flags.json").write_text(json.dumps(out, indent=1) + "\n")
    for s in M:
        print(s, {k: v for k, v in at[s].items() if "interval" not in k})
    print("pooled pick", pooled_pick)
    print("prediction 4", p4)


if __name__ == "__main__":
    main()
