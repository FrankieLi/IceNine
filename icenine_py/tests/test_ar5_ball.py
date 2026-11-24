"""AR5 step 23 (task ar5w): rule ``CB`` (converging ball) for the seed refinement's second MC stage.
Tests of the pure rule (the ball's draws, the halving, the floor and the cap, strict-less
acceptance, the generator and the evaluation count), the stage's wrapper on a real ``MCOptimizer``
with a stub cost function (run in a subprocess: importing ``ball_pick`` sets the study's caches and
replaces older freeze checks in its process), the guards on ``held2``, ``held7`` and ``fine7`` and
the link script's refusals (plan "AR5 step 23", phase A step 2). Written before the scripts were
frozen; nothing here opens a sample."""

import math
import os
import subprocess
import sys
from pathlib import Path
from typing import Any, List, Tuple

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parent.parent
AR5 = ROOT / "scripts" / "autoresearch" / "ar5"
sys.path.insert(0, str(AR5))

import ball_guard as BG  # noqa: E402
import ball_rules as BR  # noqa: E402

R0 = math.radians(0.329)


# ---- the ball's draws ----------------------------------------------------------------------------
def test_ball_angle_never_over_r_and_median() -> None:
    rng = np.random.default_rng(1)
    r = 0.01
    ang = np.array([BR.draw_ball(rng, r)[1] for _ in range(20000)])
    assert ang.max() <= r
    assert abs(np.median(ang) / r - 0.5 ** (1 / 3)) < 0.01  # 0.794 r
    assert abs(np.mean(ang < 0.5 * r) - 0.125) < 0.01  # the cube-root law: a ball, not a shell


def test_ball_axes_even_over_the_sphere() -> None:
    rng = np.random.default_rng(2)
    ax = []
    for _ in range(20000):
        m, a = BR.draw_ball(rng, 0.3)
        k = np.array([m[2, 1] - m[1, 2], m[0, 2] - m[2, 0], m[1, 0] - m[0, 1]]) / (2 * math.sin(a))
        ax.append(k)
    ax = np.array(ax)
    assert np.allclose(np.linalg.norm(ax, axis=1), 1.0, atol=1e-9)
    assert np.all(np.abs(ax.mean(axis=0)) < 0.02)
    assert np.allclose((ax**2).mean(axis=0), 1 / 3, atol=0.01)  # even over the sphere
    assert abs(np.mean(ax[:, 2] > 0.5) - 0.25) < 0.01  # a cap of height 0.5 holds a quarter


def test_ball_rotation_is_proper_and_four_numbers() -> None:
    rng = np.random.default_rng(3)
    m, _ = BR.draw_ball(rng, 0.2)
    assert np.allclose(m @ m.T, np.eye(3), atol=1e-12)
    assert abs(np.linalg.det(m) - 1.0) < 1e-12
    a, b = np.random.default_rng(3), np.random.default_rng(3)
    BR.draw_ball(a, 0.2)
    b.normal(size=3)
    b.uniform()
    assert a.uniform() == b.uniform()  # exactly three normals and one uniform


# ---- the stage ----------------------------------------------------------------------------------
class Rec:
    """An evaluate callback that records the calls and the radius of every proposal."""

    def __init__(self, costs: Any) -> None:
        self.costs = costs  # callable(call index) -> cost
        self.calls: List[np.ndarray] = []

    def __call__(self, R: np.ndarray) -> Tuple[float, Any]:
        i = len(self.calls)
        self.calls.append(np.array(R))
        return float(self.costs(i)), ("info", i)


def test_all_fail_halves_every_ten_and_ends_at_the_floor(monkeypatch: Any) -> None:
    radii: List[float] = []
    real = BR.draw_ball

    def spy(rng: Any, r: float) -> Any:
        radii.append(r)
        return real(rng, r)

    monkeypatch.setattr(BR, "draw_ball", spy)
    ev = Rec(lambda i: 1.0)
    out = BR.ball_stage(ev, np.eye(3), R0, np.random.default_rng(0))
    assert out["end"] == "floor" and out["gains"] == 0 and out["halvings"] == 6
    assert out["proposals"] == 60 == BR.floor_proposals_min(R0)
    # r halves after exactly 10 failures in a row and at no other time
    for k, r in enumerate(radii):
        assert r == pytest.approx(R0 / 2 ** (k // 10))
    assert out["r_final"] < math.radians(0.01) <= out["r_final"] * 2
    assert out["evals"] == len(ev.calls) == 61  # the start evaluation and 60 proposals
    assert np.array_equal(out["R"], np.eye(3))  # no gain: the start


def test_gain_resets_the_count_and_never_changes_r(monkeypatch: Any) -> None:
    radii: List[float] = []
    real = BR.draw_ball
    monkeypatch.setattr(BR, "draw_ball", lambda rng, r: (radii.append(r), real(rng, r))[1])
    # a gain at proposal 10, 20, ...: nine failures in a row each time, never ten
    ev = Rec(lambda i: 1000.0 - i if i % 10 == 0 else 5000.0)
    out = BR.ball_stage(ev, np.eye(3), R0, np.random.default_rng(0))
    assert out["end"] == "cap" and out["proposals"] == 300
    assert out["halvings"] == 0 and out["gains"] == 30
    assert set(radii) == {R0}
    # a gain after nine failures, then ten failures: the halving comes at the tenth after the gain
    radii.clear()
    ev = Rec(lambda i: 1.0 - 0.5 if i == 10 else 1.0)  # start 1.0, one gain at proposal 10
    out = BR.ball_stage(ev, np.eye(3), R0, np.random.default_rng(0), cap=40)
    assert out["gains"] == 1 and out["halvings"] == 3
    assert radii[:10] == [R0] * 10 and radii[10:20] == [R0] * 10  # gain at 10; next ten fail
    assert radii[20] == R0 / 2 and radii[30] == R0 / 4


def test_floor_is_tested_before_the_cap() -> None:
    out = BR.ball_stage(Rec(lambda i: 1.0), np.eye(3), R0, np.random.default_rng(0), cap=60)
    assert out["end"] == "floor" and out["proposals"] == 60
    out = BR.ball_stage(Rec(lambda i: 1.0), np.eye(3), R0, np.random.default_rng(0), cap=59)
    assert out["end"] == "cap" and out["proposals"] == 59


def test_strict_less_moves_the_centre_equal_does_not() -> None:
    ev = Rec(lambda i: 1.0)  # equal cost every time: nothing is accepted
    out = BR.ball_stage(ev, np.eye(3), R0, np.random.default_rng(5), cap=5)
    assert out["gains"] == 0 and np.array_equal(out["R"], np.eye(3))
    # every proposal improves: each is the previous best turned by a fresh draw (left product)
    ev = Rec(lambda i: 10.0 - i)
    out = BR.ball_stage(ev, np.eye(3), R0, np.random.default_rng(5), cap=7)
    assert out["gains"] == 7 and out["cost"] == 3.0 and out["info"] == ("info", 7)
    g = np.random.default_rng(5)
    best = np.eye(3)
    for k in range(1, 8):
        best = BR.draw_ball(g, R0)[0] @ best
        assert np.allclose(ev.calls[k], best, atol=1e-15)
    assert np.allclose(out["R"], best, atol=1e-15)


def test_generator_is_continued_and_every_evaluation_counted() -> None:
    a = np.random.default_rng(11)
    a.uniform(size=17)  # stage 1 drew something
    ev = Rec(lambda i: 1.0)
    out = BR.ball_stage(ev, np.eye(3), R0, a)
    # the draws used the same generator: it is where a replay with the same radii leaves it
    g = np.random.default_rng(11)
    g.uniform(size=17)
    r = R0
    for k in range(out["proposals"]):
        BR.draw_ball(g, r)
        if (k + 1) % 10 == 0:
            r *= 0.5
    assert a.uniform() == g.uniform()
    assert out["evals"] == len(ev.calls) == 1 + out["proposals"]


def test_start_below_the_floor_ends_at_once() -> None:
    out = BR.ball_stage(Rec(lambda i: 1.0), np.eye(3), 1e-5, np.random.default_rng(0))
    assert out["end"] == "floor" and out["proposals"] == 0 and out["evals"] == 1


def test_split_and_first_stage_equal_and_lost_end() -> None:
    assert BR.split_refinement(300, 100, 197) == dict(stage1=100, stage2=198, final=2)
    assert BR.split_refinement(300, 100, 198)["final"] == 1
    a = [[3, 5.0], [40, 4.0], [250, 3.0]]
    b = [[3, 5.0], [40, 4.0], [260, 2.5]]
    assert BR.first_stage_equal(208, a, b)
    assert not BR.first_stage_equal(255, a, b)
    assert BR.lost_end(0, -1) and not BR.lost_end(-1, -1) and not BR.lost_end(2, 1)


def test_judges_follow_the_frozen_numbers() -> None:
    assert BR.judge_vs_j5sd("held", 0.10, 4, 0.001)["ok"]
    assert not BR.judge_vs_j5sd("held", 0.0999, 4, 0.0)["r1"]
    assert not BR.judge_vs_j5sd("mid", 0.2, 5, 0.0)["r2"]
    assert not BR.judge_vs_j5sd("grad", 0.2, 0, 0.0011)["r3"]
    assert (
        BR.judge_vs_j5sd("big", 0.01, 0, 0.0)["r1"]
        and not BR.judge_vs_j5sd("grad", 0.019, 0, 0)["r1"]
    )
    assert BR.judge_vs_s8(0.2, 5, 0.005)["ok"] and not BR.judge_vs_s8(0.19, 0, 0)["ok"]
    assert BR.FAIL_HALVE == 10 and BR.FLOOR_DEG == 0.01 and BR.CAP_PROPOSALS == 300


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
        BG.refuse_sample(ok, tmp_path)  # held6 is disclosed here
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
    p = tmp_path / "J5DC_held7" / "run.json"
    with pytest.raises(SystemExit):
        BG.guarded(p, tmp_path)  # held7's runs have not ended
    assert BG.guarded(tmp_path / "J5DC_held6" / "run.json", tmp_path).name == "run.json"
    assert BG.guarded(tmp_path / "J5D_mid" / "run.json", tmp_path).name == "run.json"
    (tmp_path / "logs").mkdir()
    (tmp_path / "logs" / "ball_runs.done").write_text("exit=0\n")
    assert BG.guarded(p, tmp_path) == p


def test_producer_may_open_only_its_inputs(tmp_path: Path, monkeypatch: Any) -> None:
    monkeypatch.setenv("SCREEN_PRODUCER", "1")
    ok = tmp_path / "passes_held7" / "before_J5DC.npy"
    assert BG.guarded(ok, tmp_path) == ok
    with pytest.raises(SystemExit):
        BG.guarded(tmp_path / "J5DC_held7" / "run.json", tmp_path)
    with pytest.raises(SystemExit):
        BG.guarded(tmp_path / "held2" / "tries.npz", tmp_path)  # held2 is never excepted


def test_link_script_refuses_other_samples() -> None:
    for bad in ("held2", "held7", "fine7"):
        r = subprocess.run(
            ["bash", str(AR5 / "ball_link_cache.sh"), "--dry-run", bad],
            capture_output=True,
            text=True,
        )
        assert r.returncode == 2 and "refused" in r.stderr


def test_freeze_lists_this_steps_scripts() -> None:
    import ball_freeze as BF

    cur = BF.current()
    for s in BF.OWN:
        assert f"scripts/autoresearch/ar5/{s}" in cur["files"]
    assert any(k.endswith("defer_pick.py") for k in cur["files"])  # step 22's pinned too
    assert BF.CONSTANTS["fail_halve"] == 10 and BF.CONSTANTS["cap_proposals"] == 300


# ---- the wrapper on a real MCOptimizer (subprocess: the study's set-up must not leak) ----------
CODE = r'''
import sys
sys.path.insert(0, {ar5!r})
import ball_freeze as BF
BF.problems = lambda: ""
import numpy as np
import ball_pick as BP
import high500_pick as HP
import ball_rules as BR
from icenine.orientation_search import MCOptimizer

class Info:
    def __init__(self, c): self.cost = c

class Cost:
    """cost = the angle (radians) from a fixed target rotation; counts calls"""
    def __init__(self):
        self.n = 0
        a = 0.004
        self.T = np.array([[np.cos(a), -np.sin(a), 0], [np.sin(a), np.cos(a), 0], [0, 0, 1.0]])
    def evaluate(self, R, v, ph):
        self.n += 1
        return Info(float(np.arccos(np.clip((np.trace(self.T.T @ R) - 1) / 2, -1, 1))))

class Tr:
    on = True
    n = 37

HP.High500Pick.ptracer = Tr()
BP.install()
box = np.radians(0.329)

def run(active):
    BP.BALL["active"] = active
    c = Cost()
    mc = MCOptimizer(c, None, 0, np.random.default_rng(4))
    mc._rng.uniform(size=5)
    out = mc.variance_minimizing_optimize(np.eye(3), box, 200, 3, 0.0, 0.02 ** 2)
    BP.BALL["active"] = False
    return c.n, out

n_on, out_on = run(True)
rec = BP.RECORD[-1]
assert rec["stage1"] == 37 and rec["end"] in ("floor", "cap")
assert n_on == 1 + rec["proposals"], (n_on, rec)
assert out_on.cost <= Cost().evaluate(np.eye(3), None, 0).cost  # never worse than the start
assert abs(rec["box_deg"] - 0.329) < 1e-9
# inactive: the library's own stage, untouched (the same generator, the same calls)
n_off, out_off = run(False)
c2 = Cost(); mc2 = MCOptimizer(c2, None, 0, np.random.default_rng(4)); mc2._rng.uniform(size=5)
ref = BP._ORIG(mc2, np.eye(3), box, 200, 3, 0.0, 0.02 ** 2)
assert n_off == c2.n and out_off.cost == ref.cost
assert np.array_equal(out_off.orientation, ref.orientation)
assert n_on != n_off
print("OK", n_on, n_off, rec["end"], rec["proposals"])
'''


def test_wrapper_replaces_only_the_active_stage(tmp_path: Path) -> None:
    r = subprocess.run(
        [sys.executable, "-c", CODE.format(ar5=str(AR5))],
        capture_output=True,
        text=True,
        cwd=str(ROOT),
        env={**os.environ, "RESEED_SAMPLE": "mid", "BALL_CACHE": str(tmp_path)},
    )
    assert r.returncode == 0, r.stderr[-2500:]
    assert "OK" in r.stdout
