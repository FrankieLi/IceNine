"""AR5 step 1, part 1b: every setting of the try test on the two pick samples, and the pick.

  uv run python scripts/autoresearch/ar5/stop_part1_pick.py

Pick samples: ``big`` and ``grad`` (the first gradient sample). The pick (plan): the setting with
the fewest evaluations among those that lose at most 5% of the net gain of the pass, pooled over
the two samples. P3 takes the ``d`` of the P1 pick (if no P1 setting qualifies, the P1 setting with
the highest net gain). Writes ``part1_pick.json``; it is committed before ``stop_part1_test.py``
(``mid``) is run, so that git shows the order. Nothing here reads ``mid``.
"""

import json
from typing import Any, Dict, List

import stop_part1_lib as P
import stop_rules as SR

PICK_SAMPLES = ("big", "grad")


def main() -> None:
    data = {s: P.load(s) for s in PICK_SAMPLES}
    full = {s: P.evaluate(data[s], "none", None, None) for s in PICK_SAMPLES}
    assert (full["big"]["wrong_after"], full["grad"]["wrong_after"]) == (24, 63)
    pooled_full_gain = sum(f["net_gain"] for f in full.values())

    def pooled(kind: str, d: Any, m: Any, name: str) -> Dict[str, Any]:
        per = {s: P.evaluate(data[s], kind, d, m) for s in PICK_SAMPLES}
        return dict(
            name=name,
            kind=kind,
            d=d,
            m=m,
            per_sample=per,
            evals=sum(p["evals"] for p in per.values()),
            net_gain=sum(p["net_gain"] for p in per.values()),
            tries_kept=sum(p["tries_kept"] for p in per.values()),
        )

    p1 = [pooled("P1", d, None, f"P1 d={d:g}") for d in SR.D_VALUES]
    p2 = [pooled("P2", None, m, f"P2 m={m}") for m in SR.M_VALUES]
    d_row = SR.pick_d(p1, pooled_full_gain)
    d_pick = d_row["d"]
    p1_qualifies = SR.pick_setting(p1, pooled_full_gain) is not None
    p3 = [pooled("P3", d_pick, m, f"P3 d={d_pick:g} m={m}") for m in SR.M_VALUES]
    rows: List[Dict[str, Any]] = p1 + p2 + p3
    pick = SR.pick_setting(rows, pooled_full_gain)
    out = dict(
        pick_samples=list(PICK_SAMPLES),
        full_pass=dict(
            per_sample=full, evals=sum(f["evals"] for f in full.values()), net_gain=pooled_full_gain
        ),
        min_net_gain=(1 - SR.MAX_LOSS) * pooled_full_gain,
        p1_pick_d=d_pick,
        p1_has_qualifying_setting=p1_qualifies,
        rows=rows,
        pick=(
            None
            if pick is None
            else dict(name=pick["name"], kind=pick["kind"], d=pick["d"], m=pick["m"])
        ),
        note="chosen on big and grad only; mid is read by stop_part1_test.py after this commit",
    )
    P.AR5.joinpath("part1_pick.json").write_text(json.dumps(out, indent=1) + "\n")
    print("pooled full pass: evals", out["full_pass"]["evals"], "net gain", pooled_full_gain)
    for r in rows:
        print(
            f"{r['name']:<16} evals {r['evals']:>9} net gain {r['net_gain']:>4} "
            f"kept {r['tries_kept']:>6}"
        )
    print("pick:", out["pick"])


if __name__ == "__main__":
    main()
