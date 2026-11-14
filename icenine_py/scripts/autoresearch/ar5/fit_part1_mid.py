"""AR5 step 2, Part 1: the pass rule on ``mid`` (run after ``pick.json`` is committed).

  scripts/autoresearch/ar5/fit_bg.sh mid part1_mid fit_part1_mid.py

Reads the replay shards of ``mid`` (``fit_part1_replay.py``; ``fit_common`` refuses to import with
``RESEED_SAMPLE=mid`` while ``pick.json`` is not committed) and writes ``part1_mid.json`` (arm
counts and paired splits, as for the pick samples) and ``part1_mid_pass.json`` (for every arm:
rules (a),
(b), the fall of the seed-step evaluations and rule (c), by ``fit_rules.mid_pass``; the plan's
verdict is the one of the picked arm, ``pick`` in ``pick.json``). Records the hash of the commit
that holds ``pick.json`` and the time at which the first ``mid`` file was written.
"""

import json

import fit_common as FC
import fit_part1_lib as L
import fit_rules as FR


def main() -> None:
    assert FC.SAMPLE == "mid"
    pick = json.loads(FC.PICK_FILE.read_text())["pick"]
    rows = L.load("mid")
    n = len(rows)
    counts = L.counts(rows)
    part1 = dict(
        n_steps=n,
        counts=counts,
        paired_passed={a: L.paired(rows, a, "passed") for a in FC.ARMS if a != "M"},
        paired_right={a: L.paired(rows, a, "right") for a in FC.ARMS if a != "M"},
        log_mismatches_M=sum(1 for r in rows if not r["log_ok"]),
        err_q=L.err_quantiles(rows),
    )
    (FC.FIT / "part1_mid.json").write_text(json.dumps(part1, indent=1) + "\n")
    arms = {a: FR.mid_pass(counts[a], counts["M"], n) for a in FC.ARMS if a != "M"}
    out = dict(
        pick=pick,
        pick_commit=FC.pick_commit(),
        verdict_of_pick=arms.get(pick),
        arms=arms,
        n_steps=n,
    )
    (FC.FIT / "part1_mid_pass.json").write_text(json.dumps(out, indent=1) + "\n")
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
