"""AR5 step 7, Phase A step 3: write ``benchmarks/autoresearch/ar5/gate/frozen.json``.

  uv run python scripts/autoresearch/ar5/gate_freeze.py

Frozen before any script of this step opens a pass file or a map of ``held``: the skip cut ``h_s``,
the decision rule (the numbers in ``gate_rules.py``), the predictions (quoted word for word from
the plan) and the SHA-256 of the disclosed passes' files (tries, maps, pass.json, run.json). It
refuses to overwrite an existing file. It does not touch any file of ``held``.
"""

import datetime
import json
import subprocess

import gate_common as GC
import gate_rules as GR

PREDICTIONS = [
    "1. On `held` the replay at 0.95 cuts 55% to 85% of the pass and fixed / broken against the "
    "full pass is at most 0 / 2.",
    "2. On `held`, `S11v2` with the skipped pass needs 5% to 12% fewer evaluations than with the "
    "full pass (its pass is a smaller share: an estimate of 3.5 of about 36 million, from step 4's "
    "totals; this draft did not open `held`'s pass files).",
    "3. The live pass spends exactly the replay's evaluations (61 per kept try) on every sample.",
    "4. Live wrong voxels are within 2 of the replay's on every sample (one generator for the "
    "pass, so skipped tries shift the random draws of later tries).",
    "5. Median error rises by at most 0.0005 degree on every sample.",
]
RULE_TEXT = (
    "The skip passes on a sample if: the pass is cut by at least 50%; broken minus fixed against "
    "the full pass is at most 5% of the pass's net gain, rounded up, at least 1 voxel; median "
    "error is at most 0.001 degree higher; and the run with the skipped pass meets the house rules "
    "against `S8` (20% fewer, broken minus fixed at most 5, median at most 0.005 degree higher). "
    "It passes the step if it passes on `held` live and on every disclosed sample. Then the owner "
    "decides whether the documentation recommends the key. If it fails on `held`: the key stays, "
    "off, the page says so; `h_s` is not tuned on `held`."
)


def main() -> None:
    out = GC.GATE / "frozen.json"
    if out.exists():
        raise SystemExit(f"refusing to overwrite {out}")
    files = {}
    for sample, run in GC.DISCLOSED:
        for k, p in GC.pass_files(sample, run).items():
            files[f"{sample}/{run}/{k}"] = GC.sha256(p)
    head = subprocess.run(
        ["git", "rev-parse", "HEAD"], capture_output=True, text=True, check=True
    ).stdout.strip()
    GC.GATE.mkdir(parents=True, exist_ok=True)
    out.write_text(
        json.dumps(
            dict(
                frozen_at=datetime.datetime.now().astimezone().isoformat(timespec="seconds"),
                code_head=head,
                h_s=GR.H_S,
                h_s_note="the library's SeedIndexStopHitRatio and step 1's flag pick on big and "
                "grad; no other value is a decision; the 0.90-0.99 curve is a what-if",
                rule=dict(
                    text=RULE_TEXT,
                    min_cut=GR.MIN_CUT,
                    loss_share_pct=GR.LOSS_SHARE_PCT,
                    loss_min_voxels=1,
                    median_tol_full_deg=GR.MEDIAN_TOL_FULL,
                    house_fewer_pct=GR.HOUSE_FEWER_PCT,
                    house_broken_minus_fixed=GR.HOUSE_BROKEN_MINUS_FIXED,
                    house_median_tol_deg=GR.HOUSE_MEDIAN_TOL,
                ),
                predictions=PREDICTIONS,
                held_read=False,
                disclosed_sha256=files,
            ),
            indent=1,
        )
        + "\n"
    )
    print("wrote", out, len(files), "file hashes")


if __name__ == "__main__":
    main()
