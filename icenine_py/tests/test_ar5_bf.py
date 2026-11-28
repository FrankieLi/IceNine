"""AR5 step 24 (task ar5x): rule ``BF`` (ball first, the stock second stage only when no ball end
reaches 0.95). Tests of the pure arithmetic (the fallback's extra evaluations, the step's total,
the judges against the frozen numbers), the two-pass seed step on a scripted pick (run in a
subprocess: importing ``bf_pick`` sets the study's caches and replaces older freeze checks in its
process): the ball pass stops at the first end at 0.95, the fallback runs only when the ball pass
has none, no CMA fit is run twice, the first stage is checked equal to the kept one, every
evaluation is counted (the ball stage of a step that falls back included), rule off gives
``J5D``'s step; the guards on ``held2``, ``held7`` and ``fine7`` and the link script's refusals
(plan "AR5 step 24", phase A step 2). Written before the scripts were frozen; nothing here opens a
sample."""

import os
import subprocess
import sys
from pathlib import Path
from typing import Any

import pytest

ROOT = Path(__file__).resolve().parent.parent
AR5 = ROOT / "scripts" / "autoresearch" / "ar5"
sys.path.insert(0, str(AR5))

import bf_guard as BG  # noqa: E402
import bf_rules as BRU  # noqa: E402


# ---- the arithmetic ------------------------------------------------------------------------------
def fits_row(pos: int, j: int, kind: str, evals: int) -> list:
    return [pos, j, kind, 0.9, 1.0, evals, [0.0] * 9]


def test_fallback_extra_is_refinements_minus_first_stages() -> None:
    fits = [
        fits_row(0, 3, "C", 201),
        fits_row(0, 3, "M", 1000),
        fits_row(1, 7, "C", 201),
        fits_row(1, 7, "M", 700),
    ]
    out = BRU.fallback_extra(fits, {3: 208, 7: 190})
    assert out == dict(
        extra=(1000 - 208) + (700 - 190), refinements=2
    )  # the CMA fits are not in it


def test_fallback_extra_refuses_a_refinement_the_ball_pass_never_ran() -> None:
    with pytest.raises(AssertionError):
        BRU.fallback_extra([fits_row(0, 5, "M", 900)], {3: 208})


def test_step_evals_ball_found_has_no_extra() -> None:
    assert BRU.step_evals(5000, 2, None) == 5000
    assert BRU.step_evals(5000, 2, 0) == 5000
    with pytest.raises(AssertionError):
        BRU.step_evals(5000, 2, 10)
    assert BRU.step_evals(5000, -1, 1234) == 6234


def test_lost_end_and_feature_row() -> None:
    assert BRU.lost_end(1, -1) and not BRU.lost_end(-1, -1) and not BRU.lost_end(1, 0)
    r = BRU.feature_row(4, [dict(hit=0.8, cost=0.5)], [dict(hit=0.7, cost=0.6)], True, [0.1])
    assert r["ball_hit"] == [0.8] and r["cma_cost"] == [0.6] and r["fallback_found"] is True


def test_judges_follow_the_frozen_numbers() -> None:
    j = BRU.JUDGE
    assert j["r1_fewer_min"]["held"] == 0.04 and j["r1_fewer_min"]["mid"] == 0.04
    assert j["r1_fewer_min"]["grad"] == 0.015 and j["r1_fewer_min"]["big"] == 0.01
    assert j["r2_broken_minus_fixed_max"] == 3
    assert BRU.judge_vs_j5sd("held", 0.0536, 3, 0.001)["ok"]
    assert not BRU.judge_vs_j5sd("held", 0.039, 0, 0.0)["ok"]
    assert not BRU.judge_vs_j5sd("big", 0.02, 4, 0.0)["ok"]
    assert not BRU.judge_vs_j5sd("big", 0.02, 0, 0.0011)["ok"]
    assert BRU.judge_vs_s8(0.2, 5, 0.005)["ok"] and not BRU.judge_vs_s8(0.19, 0, 0.0)["ok"]
    g = BRU.GATE
    assert g["g2_fewer_min"] == 0.04 and g["g2_broken_minus_fixed_max"] == 3
    assert g["g2_grid_more_max"] == 2 and g["g2_seed_test_more_max"] == 10
    assert g["g1_fewer_min"] == 0.10 and g["g1_lost_max"] == 0.03
    assert BRU.CAPS_MIN == dict(S8=480, rest=300)
    assert BRU.FAIL_HALVE == 10 and BRU.CAP_PROPOSALS == 300 and BRU.FLOOR_DEG == 0.01


# ---- the guards ---------------------------------------------------------------------------------
def test_guards_refuse_held2_fine7_and_held7_early(tmp_path: Path) -> None:
    for bad in ("held2", "fine7", "myheld2x"):
        with pytest.raises(SystemExit) as e:
            BG.refuse_sample(bad, tmp_path)
        assert e.value.code == 2
    with pytest.raises(SystemExit):
        BG.refuse_sample("nonsense", tmp_path)
    with pytest.raises(SystemExit):
        BG.refuse_sample("held7", tmp_path)  # no gate_passed yet
    for ok in ("mid", "held", "held6", "grad", "big", "held5", ""):
        BG.refuse_sample(ok, tmp_path)
    (tmp_path / "logs").mkdir()
    (tmp_path / "logs" / "gate_passed").write_text("ok\n")
    BG.refuse_sample("held7", tmp_path)


def test_guarded_paths(tmp_path: Path) -> None:
    for p in (
        tmp_path / "held2" / "x.json",
        tmp_path / "a" / "fine7" / "y",
        ROOT / "benchmarks" / "autoresearch" / "fine1" / "z.json",
    ):
        with pytest.raises(SystemExit):
            BG.guarded(p, tmp_path)
    p = tmp_path / "J5DB_held7" / "run.json"
    with pytest.raises(SystemExit):
        BG.guarded(p, tmp_path)  # held7's runs have not ended
    assert BG.guarded(tmp_path / "J5DB_held6" / "run.json", tmp_path).name == "run.json"
    assert BG.guarded(tmp_path / "J5D_mid" / "run.json", tmp_path).name == "run.json"
    (tmp_path / "logs").mkdir()
    (tmp_path / "logs" / "bf_runs.done").write_text("exit=0\n")
    assert BG.guarded(p, tmp_path) == p


def test_producer_may_open_only_its_inputs(tmp_path: Path, monkeypatch: Any) -> None:
    monkeypatch.setenv("SCREEN_PRODUCER", "1")
    ok = tmp_path / "passes_held7" / "before_J5DB.npy"
    assert BG.guarded(ok, tmp_path) == ok
    with pytest.raises(SystemExit):
        BG.guarded(tmp_path / "J5DB_held7" / "run.json", tmp_path)
    with pytest.raises(SystemExit):
        BG.guarded(tmp_path / "held2" / "tries.npz", tmp_path)  # held2 is never excepted


def test_link_script_refuses_other_samples() -> None:
    for bad in ("held2", "held7", "fine7"):
        r = subprocess.run(
            ["bash", str(AR5 / "bf_link_cache.sh"), "--dry-run", bad],
            capture_output=True,
            text=True,
        )
        assert r.returncode == 2 and "refused" in r.stderr


def test_freeze_lists_this_steps_scripts() -> None:
    import bf_freeze as BFZ

    cur = BFZ.current()
    for s in BFZ.OWN:
        assert f"scripts/autoresearch/ar5/{s}" in cur["files"]
    assert any(k.endswith("ball_pick.py") for k in cur["files"])  # step 23's pinned too
    assert any(k.endswith("defer_pick.py") for k in cur["files"])  # step 22's too
    assert BFZ.CONSTANTS["fail_halve"] == 10 and BFZ.CONSTANTS["cap_proposals"] == 300


# ---- the two-pass step on a scripted pick (subprocess: the study's set-up must not leak) ---------
CODE = r"""
import sys
sys.path.insert(0, {ar5!r})
import bf_freeze as BFZ
BFZ.problems = lambda: ""
import numpy as np
import ball_freeze as _BF
_BF.problems = lambda: ""
import bf_pick as BFP
import ball_pick as BP
import defer_pick as DP
import high500_pick as HP
import order_run as ORUN
from types import SimpleNamespace

R = [0.0] * 9
S = dict(ball_stop=-1, fb_stop=-1, bad_first_stage=False, c_calls=[], fb_calls=0)


def fit_c(self, kind, rec, j, Rj, cj, vertices, phase):
    S["c_calls"].append(int(j))
    return dict(j=int(j), kind="C", hit=0.7, cost=1.0, evals=201, R=np.zeros((3, 3)))


HP.High500Pick._fit = fit_c


def trace(bad):
    return [[1, 3.0], [50, 2.0 if not bad else 2.5]]


def fake_e1(self, vertices, phase_index, n_cost=5, n_refine=2, tau=-1):
    ORUN.COUNTS["seed_steps"] += 1
    DP.STATS["steps"] += 1
    n = ORUN.COUNTS["seed_steps"]
    ball = BP.BALL["on"]
    fits, tr_rows, fin = [], [], []
    stop = S["ball_stop"] if ball else S["fb_stop"]
    if not ball:
        S["fb_calls"] += 1
    for p, j in enumerate((4, 9)):
        c = self._fit("C", None, j, None, 0.0, None, 0)
        fits.append([p, j, "C", c["hit"], c["cost"], c["evals"], R])
        if ball:
            ev = 208 + 31 + 1  # first stage 208, ball stage 30 proposals + start, final 1
            BP.LOG.append(dict(voxel=self.cur, step=n, j=j, evals=ev, stage1=208, proposals=30,
                               gains=2, halvings=1, end="floor", r_final_deg=0.005,
                               start_cost=1.0, end_cost=0.5, box_deg=0.329, trace=trace(False)))
            tr_rows.append(dict(voxel=self.cur, step=n, j=j, evals=ev, trace=trace(False)))
        else:
            ev = 900
            tr_rows.append(dict(voxel=self.cur, step=n, j=j, evals=ev,
                                trace=trace(S["bad_first_stage"] and j == 9)))
        fits.append([p, j, "M", 0.8, 0.5, ev, R])
        fin.append(201 + ev)
    ORUN.TRACES.extend(tr_rows)
    costed = 20
    total = costed + sum(fin)
    entry = dict(
        voxel=self.cur, reseed=self.reseed, stop=stop, fits=fits, evals=total, costed=costed
    )
    self.e1_log.append(entry)
    DP.LOG.append(dict(voxel=int(self.cur), evals=total))
    return SimpleNamespace(evals=total, finisher_evals=list(fin), refined=[4, 9])


DP.DeferPick.e1_step = fake_e1


def mk():
    p = object.__new__(BFP.BFPick)
    p.e1_log, p.cur, p.reseed = [], 11, False
    return p


def reset():
    for k in list(BFP.STATS):
        BFP.STATS[k] = 0
    del BFP.LOG[:], BFP.FEATURES[:], BP.LOG[:], ORUN.TRACES[:], DP.LOG[:]
    ORUN.COUNTS["seed_steps"] = 0
    DP.STATS["steps"] = 0
    S.update(c_calls=[], fb_calls=0)


# 1. a step that ends in the ball pass: no fallback, the step is the ball pass's
reset(); S.update(ball_stop=1, fb_stop=-1, bad_first_stage=False)
p = mk(); out = p.e1_step(None, 0)
assert S["fb_calls"] == 0 and len(p.e1_log) == 1 and out.evals == p.e1_log[0]["evals"]
assert BFP.STATS["ball_found"] == 1 and BFP.STATS["fallback_steps"] == 0
assert ORUN.COUNTS["seed_steps"] == 1 and len(BFP.FEATURES) == 0
assert p.e1_log[0]["bf"]["fallback"] is False

# 2. no ball end at 0.95: the fallback runs, no CMA fit twice, every evaluation counted once
reset(); S.update(ball_stop=-1, fb_stop=-1, bad_first_stage=False)
p = mk(); out = p.e1_step(None, 0)
assert S["fb_calls"] == 1 and sorted(S["c_calls"]) == [4, 9], S["c_calls"]  # CMA fits run once
assert len(p.e1_log) == 1 and ORUN.COUNTS["seed_steps"] == 1 and DP.STATS["steps"] == 1
assert len(DP.LOG) == 1 and len(ORUN.TRACES) == 2  # one row per step; the fallback's traces only
ball_evals = 20 + 2 * (201 + 240)
extra = 2 * (900 - 208)
assert out.evals == ball_evals + extra == p.e1_log[0]["evals"], (out.evals, ball_evals, extra)
assert BFP.STATS["fallback_none"] == 1 and BFP.STATS["fallback_extra_evals"] == extra
assert BFP.STATS["ball_waste_evals"] == 2 * 31 and BFP.STATS["stage1_checked"] == 2
assert len(BFP.FEATURES) == 1 and BFP.FEATURES[0]["fallback_found"] is False
assert p.e1_log[0]["bf"]["fallback"] is True
assert BP.BALL["on"] is True and p._cma_use is None  # the switch is put back

# 3. a fallback that finds 0.95
reset(); S.update(ball_stop=-1, fb_stop=0, bad_first_stage=False)
p = mk(); out = p.e1_step(None, 0)
assert BFP.STATS["fallback_found"] == 1 and p.e1_log[0]["stop"] == 0

# 4. a first stage that differs from the kept one stops the run
reset(); S.update(ball_stop=-1, fb_stop=-1, bad_first_stage=True)
p = mk()
try:
    p.e1_step(None, 0)
    raise SystemExit("no mismatch raised")
except BFP.BFStageMismatch:
    pass
assert BFP.STATS["stage1_mismatch"] == 1 and BP.BALL["on"] is True

# 5. rule off: the seed step is J5D's (the stock stage only, CMA fits run, no ball)
reset(); BFP.BF["on"] = False; S.update(ball_stop=-1, fb_stop=-1, bad_first_stage=False)
p = mk(); out = p.e1_step(None, 0)
assert BP.BALL["on"] is False and S["fb_calls"] == 1 and sorted(S["c_calls"]) == [4, 9]
assert out.evals == 20 + 2 * (201 + 900) and len(BFP.LOG) == 0 and len(BP.LOG) == 0
BFP.BF["on"] = True
print("OK")
"""


def test_two_pass_step_on_a_scripted_pick(tmp_path: Path) -> None:
    r = subprocess.run(
        [sys.executable, "-c", CODE.format(ar5=str(AR5))],
        capture_output=True,
        text=True,
        cwd=str(ROOT),
        env={**os.environ, "RESEED_SAMPLE": "mid", "BF_CACHE": str(tmp_path)},
    )
    assert r.returncode == 0, r.stderr[-3000:]
    assert "OK" in r.stdout


def test_pick_inherits_the_library_stages_unchanged() -> None:
    """The ball stage is step 23's (``BallPick._fit_m``) and the fallback's refinement is step 22's
    (``DeferPick._fit_m`` through the ball switch): ``BFPick`` overrides only ``_fit`` (the CMA
    cache) and ``e1_step`` (the two passes)."""
    src = (AR5 / "bf_pick.py").read_text()
    assert "def _fit_m" not in src and "def _fit(" in src and "def e1_step(" in src
    assert "default_rng" not in src and "random" not in src.replace("Draws no random number", "")
