"""AR5 step 22: the summary's pure decision functions (frozen rules of ``defer_rules.JUDGE``) and
the stored result's consistency. Nothing here opens a sample."""

import json
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
AR5 = ROOT / "scripts" / "autoresearch" / "ar5"
sys.path.insert(0, str(AR5))
sys.path.insert(0, str(ROOT / "scripts" / "common"))

import defer_rules as DR  # noqa: E402
import defer_summary as DS  # noqa: E402


def test_rules_against_j5s_use_the_frozen_thresholds() -> None:
    assert DS.judge_vs_j5s("held", 0.07, 3, 0.001)["ok"]
    assert not DS.judge_vs_j5s("held", 0.0699, 0, 0.0)["r1"]
    assert not DS.judge_vs_j5s("held", 0.2, 4, 0.0)["r2"]
    assert not DS.judge_vs_j5s("held", 0.2, 0, 0.0011)["r3"]
    assert DS.judge_vs_j5s("grad", 0.015, 0, 0.0)["r1"]
    assert not DS.judge_vs_j5s("grad", 0.0149, 0, 0.0)["r1"]
    assert (
        DS.judge_vs_j5s("big", 0.03, 0, 0.0)["r1"] and not DS.judge_vs_j5s("big", 0.029, 0, 0)["r1"]
    )
    # held6 is a held type: 7%
    assert (
        not DS.judge_vs_j5s("held6", 0.069, 0, 0.0)["r1"]
        and DS.judge_vs_j5s("held6", 0.07, 0, 0)["r1"]
    )


def test_house_rule_against_s8() -> None:
    assert DS.judge_vs_s8(0.20, 5, 0.005)["ok"]
    assert not DS.judge_vs_s8(0.19, 0, 0.0)["fewer"]
    assert not DS.judge_vs_s8(0.5, 6, 0.0)["broken"]
    assert not DS.judge_vs_s8(0.5, 0, 0.0051)["median"]


def test_step_passes_only_if_held6_and_every_disclosed_sample_pass() -> None:
    ok = {"vs_j5s": {"ok": True}, "vs_s8": {"ok": True}}
    bad = {"vs_j5s": {"ok": False}, "vs_s8": {"ok": True}}
    rows = {s: ok for s in DR.DISCLOSED + ("held6",)}
    assert DS.step_verdict(rows)["step_passes"]
    assert not DS.step_verdict({**rows, "mid": bad})["step_passes"]
    assert not DS.step_verdict({**rows, "held6": bad})["step_passes"]
    part = {s: ok for s in DR.DISCLOSED}  # no held6: not a pass
    assert not DS.step_verdict(part)["step_passes"]


def test_in_range_is_inclusive() -> None:
    assert (
        DS.in_range(0.3, 0.3, 0.7)
        and DS.in_range(0.7, 0.3, 0.7)
        and not DS.in_range(0.72, 0.3, 0.7)
    )


def test_stored_summary_is_consistent_when_present() -> None:
    p = ROOT / "benchmarks" / "autoresearch" / "ar5" / "defer" / "summary.json"
    if not p.exists():
        pytest.skip("summary.json not written")
    d = json.loads(p.read_text())
    for s, r in d["rows"].items():
        t = r["totals"]
        assert r["fewer"]["j5s"] == pytest.approx(1 - t["j5sd"] / t["j5s"])
        s_ = r["stages"]
        # the stages add up to the J5SD total (index and grid fall-back evaluations are the search)
        assert sum(s_.values()) == t["j5sd"], s
    assert d["verdict"]["step_passes"] == all(d["verdict"]["per_sample"].values())
