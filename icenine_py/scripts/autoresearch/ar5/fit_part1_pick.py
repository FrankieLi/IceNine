"""AR5 step 2, Part 1: the pick, from the replays of the two pick samples (``big`` and ``grad``).

  uv run python scripts/autoresearch/ar5/fit_part1_pick.py

Pick rule (plan): on the pick samples pooled, the cheapest arm (total seed-step evaluations) for
which (a) the steps that pass the seed test are at most 1% of all steps fewer than with ``M`` and
(b) the steps that end within 1 degree of the truth are at most 1% of all steps fewer than with
``M``. Counts are net per arm; the paired split per flag is written beside them. This script reads
the pick samples only; it never reads ``mid`` and has no way to (the sample list is fixed here).
Writes ``benchmarks/autoresearch/ar5/fit/pick.json`` (tracked) and is committed before ``mid`` runs.
"""

import json
import os

os.environ.setdefault("RESEED_SAMPLE", "big")

import fit_common as FC  # noqa: E402
import fit_part1_lib as L  # noqa: E402
import fit_rules as FR  # noqa: E402

PICK_SAMPLES = ("big", "grad")


def main() -> None:
    per = {s: L.load(s) for s in PICK_SAMPLES}
    rows = [r for s in PICK_SAMPLES for r in per[s]]
    n = len(rows)
    for s in PICK_SAMPLES:
        assert per[s], f"no replay rows for {s}"
    pooled = L.counts(rows)
    pick = FR.pick_arm(pooled, n)
    out = dict(
        rule="cheapest arm with (a) passes and (b) right-within-1-degree each at most 1% of all "
        "steps fewer than M (net counts, pooled over the pick samples)",
        pick_samples=list(PICK_SAMPLES),
        n_steps=n,
        n_steps_by_sample={s: len(per[s]) for s in PICK_SAMPLES},
        log_mismatches_M={s: sum(1 for r in per[s] if not r["log_ok"]) for s in PICK_SAMPLES},
        counts_pooled=pooled,
        counts_by_sample={s: L.counts(per[s]) for s in PICK_SAMPLES},
        paired_passed={a: L.paired(rows, a, "passed") for a in FC.ARMS if a != "M"},
        paired_right={a: L.paired(rows, a, "right") for a in FC.ARMS if a != "M"},
        one_percent_of_steps=n / 100.0,
        rule_a={a: FR.mid_pass(pooled[a], pooled["M"], n)["a"] for a in FC.ARMS if a != "M"},
        rule_b={a: FR.mid_pass(pooled[a], pooled["M"], n)["b"] for a in FC.ARMS if a != "M"},
        evals_fall_vs_M={a: FR.fall(pooled[a]["evals"], pooled["M"]["evals"]) for a in FC.ARMS},
        pick=pick,
    )
    FC.PICK_FILE.write_text(json.dumps(out, indent=1) + "\n")
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
