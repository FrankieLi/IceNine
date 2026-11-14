"""AR5 step 11, phase A step 4: the gate, from the replay files (disclosed samples only).

  uv run python scripts/autoresearch/ar5/patience_gate.py

Reads ``patience/replay_<sample>[_sIofN].json`` (``patience_replay.py``; every file of a sample
that is not tagged; the shards must cover every seed step that ran an MC refinement) and writes
``patience/gate.json``: the plan's gate (g1 to g4) and predictions 1 and 2. Exit 1 when the gate
fails (stop, tune nothing). Wording: plan "AR5 step 11", "Decision rules".

* g1 every uncut step, and every refinement cut after its last improvement, equals G's stored one
  (orientation, cost): no fit is ``bad_structure``; a fit is ``moved`` only where
  ``patience_rules.cut_point`` on the stored trace says the stop comes before a later gain; a CMA
  fit never changes.
* g2 taken ends changed at most 30 (mid), 3 (grad), 6 (big).
* g3 no changed taken end is wrong (over 1 degree) where G's was right, and at most 1 more failed
  seed test than G on a sample.
* g4 the replay saves at least 8% of G's run (the total with the pass) on mid, every extra
  evaluation counted (a rescue or a grid search a changed step now makes is in ``new_total``).
"""

import glob
import json
import sys
from pathlib import Path
from typing import Any, Dict, List

import numpy as np

HERE = Path(__file__).resolve().parent
PAT = HERE.parents[2] / "benchmarks" / "autoresearch" / "ar5" / "patience"
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent.parent / "common"))

import stats as ST  # noqa: E402

TOTAL = {"mid": 13004457, "grad": 10196692, "big": 3807747}
STEPS_WITH_MC = {"mid": 770, "grad": 478, "big": 99}
MAX_CHANGED = {"mid": 30, "grad": 3, "big": 6}
PRED_CHANGED = {"mid": (20, 28), "grad": (0, 2), "big": (3, 5)}
PRED_SAVE = {"mid": (0.09, 0.114), "grad": (0.0, 0.003), "big": (0.04, 0.061)}
G4_MID = 0.08


def load(s: str) -> List[Dict[str, Any]]:
    files = sorted(
        f
        for f in glob.glob(str(PAT / f"replay_{s}*.json"))
        if "_pilot" not in f and "_n" not in Path(f).name.replace(f"replay_{s}", "")
    )
    rows: List[Dict[str, Any]] = []
    for f in files:
        rows += json.loads(Path(f).read_text())["rows"]
    return rows


def sample_gate(s: str, rows: List[Dict[str, Any]]) -> Dict[str, Any]:
    n = len(rows)
    complete = n == STEPS_WITH_MC[s] and len({r["step"] for r in rows}) == n
    bad_struct = [r["step"] for r in rows if r["bad"] or "bad_structure" in r["classes"]]
    moved_unpredicted = []
    moved = after_last = same = 0
    for r in rows:
        moved += r["classes"].count("moved")
        after_last += r["classes"].count("after_last")
        same += r["classes"].count("same")
    # a "moved" fit is allowed only where the stored trace predicts it: the replay counted
    # prediction mismatches per step (live kind differs from cut_point's kind)
    pred_mismatch_steps = [r["step"] for r in rows if r["pred_mismatch"]]
    for r in rows:
        if r["pred_mismatch"] and "moved" in r["classes"]:
            moved_unpredicted.append(r["step"])
    g1_ok = bool(complete and not bad_struct and not moved_unpredicted)
    changed = [r for r in rows if not r["taken_same"]]
    wrong_new = [
        r["step"]
        for r in changed
        if r["err_new"] is not None
        and r["err_new"] > 1.0
        and (r["err_g"] is None or r["err_g"] <= 1.0)
    ]
    err_g_unknown = [r["step"] for r in changed if r["err_g"] is None]
    failed_g = int(sum(1 for r in rows if r["g_failed"]))
    failed_new = int(sum(1 for r in rows if r["new_failed"]))
    g_total = int(sum(r["g_total"] for r in rows))
    new_total = int(sum(r["new_total"] for r in rows))
    saved = g_total - new_total
    share = saved / TOTAL[s]
    extra_rows = [r["step"] for r in rows if r["new_total"] > r["g_total"]]
    mis = [r["mis_to_g_deg"] for r in changed]
    stop_type = [r for r in changed if r["g_stopped"]]
    still = [r for r in stop_type if r["new_hit"] >= 0.95]
    lo, hi = ST.wilson(len(changed), n) if n else (0.0, 0.0)
    g2_ok = len(changed) <= MAX_CHANGED[s]
    g3_ok = (not wrong_new) and failed_new - failed_g <= 1
    pc = PRED_CHANGED[s]
    ps = PRED_SAVE[s]
    return dict(
        sample=s,
        steps_replayed=n,
        steps_expected=STEPS_WITH_MC[s],
        complete=complete,
        fits=dict(same=same, after_last=after_last, moved=moved),
        g1_ok=g1_ok,
        g1_bad_structure_steps=bad_struct,
        g1_moved_where_not_predicted=moved_unpredicted,
        prediction_mismatch_steps=pred_mismatch_steps,
        taken_changed=len(changed),
        taken_changed_wilson95_of_replayed_steps=[lo, hi],
        taken_changed_limit=MAX_CHANGED[s],
        g2_ok=bool(g2_ok),
        changed_wrong_where_g_right=wrong_new,
        changed_with_g_error_unknown=err_g_unknown,
        failed_seed_tests=dict(g=failed_g, s=failed_new, extra=failed_new - failed_g),
        g3_ok=bool(g3_ok),
        g_total_of_replayed_steps=g_total,
        new_total_of_replayed_steps=new_total,
        saved_evals=saved,
        saved_share_of_g_run=share,
        steps_costing_more_than_g=extra_rows,
        changed_paths={
            p: sum(1 for r in changed if r["path"] == p) for p in ("seed", "rescued", "grid")
        },
        grid_searches_new=int(sum(1 for r in rows if r["path"] == "grid" and not r["path_same"])),
        g4_ok=bool(share >= G4_MID) if s == "mid" else None,
        changed_stop_type=len(stop_type),
        changed_stop_type_still_095=len(still),
        changed_mis_to_g_deg=dict(
            min=float(min(mis)) if mis else None,
            median=float(np.median(mis)) if mis else None,
            max=float(max(mis)) if mis else None,
        ),
        changed_err_new_max=float(max(r["err_new"] for r in changed)) if changed else None,
        prediction1=dict(
            g1=g1_ok,
            changed_in_range=bool(pc[0] <= len(changed) <= pc[1]),
            changed_range=list(pc),
            stop_type_still_095_12_to_20=(bool(12 <= len(still) <= 20) if s == "mid" else None),
            mis_in_0_01_to_0_3_median_under_0_1=bool(
                mis and min(mis) >= 0.01 and max(mis) <= 0.3 and float(np.median(mis)) < 0.1
            ),
            none_wrong=bool(
                not any(r["err_new"] is not None and r["err_new"] > 1.0 for r in changed)
            ),
            seed_test_failures_within_1=bool(abs(failed_new - failed_g) <= 1),
        ),
        prediction2=dict(
            saving_in_range=bool(ps[0] <= share <= ps[1]), range=list(ps), saving=share
        ),
    )


def main() -> None:
    per: Dict[str, Any] = {}
    for s in ("mid", "grad", "big"):
        per[s] = sample_gate(s, load(s))
    ok = all(per[s]["g1_ok"] and per[s]["g2_ok"] and per[s]["g3_ok"] for s in per) and bool(
        per["mid"]["g4_ok"]
    )
    out = dict(
        gate_passed=bool(ok),
        verdict="READY TO RUN" if ok else "STOP: the gate failed (no tuning, no other k)",
        items=dict(
            g1={s: per[s]["g1_ok"] for s in per},
            g2={s: per[s]["g2_ok"] for s in per},
            g3={s: per[s]["g3_ok"] for s in per},
            g4_mid=per["mid"]["g4_ok"],
        ),
        per_sample=per,
    )
    (PAT / "gate.json").write_text(json.dumps(out, indent=1) + "\n")
    print(json.dumps(out["items"], indent=1))
    print(out["verdict"])
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
