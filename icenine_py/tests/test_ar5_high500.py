"""AR5 step 16 (task ar5p): rule PH(500), the high stop alone at a 500-evaluation stall in the seed
step's MC refinement, on synthetic sequences. Plan "AR5 step 16", phase A step 2: stop at 0.95 only
(never low), go on under 0.95, re-check after a gain, the first lowest cost kept, one evaluation per
check, ``PH(3000)`` absent, other draws untouched, guards on ``held2`` and on ``held3``'s logs.
"""

import os
import subprocess
import sys
from pathlib import Path
from typing import Any, Callable, Dict, List

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "scripts" / "autoresearch" / "ar5"))

import band_rules as BR  # noqa: E402
import high500_guard as HG  # noqa: E402
import high500_rules as HR  # noqa: E402
import stallhit_rules as SR  # noqa: E402


class _Out:
    def __init__(self, cost: float) -> None:
        self.cost = cost


class _FakeCost:
    """A cost function whose cost is the first coordinate of the orientation."""

    def evaluate(self, orientation: Any, *args: Any, **kw: Any) -> _Out:
        return _Out(float(np.asarray(orientation).ravel()[0]))


def _drive(costs: List[float], k: int, hit: Callable[[float], float]) -> Dict[str, Any]:
    """Feed ``costs`` through the wrapped cost function exactly as ``high500_pick`` arms the meter
    (``begin(k, None, ())``: no low side); ``hit(best cost)`` plays the hit ratio."""
    meter = BR.BandMeter()
    calls: List[float] = []

    def hit_fn(R: Any) -> float:
        calls.append(float(np.asarray(R).ravel()[0]))
        return hit(float(np.asarray(R).ravel()[0]))

    meter.hit_fn = hit_fn
    patched = BR.patch_evaluate(_FakeCost.evaluate, meter)
    fake = _FakeCost()
    meter.begin(k, None, ())
    stopped = False
    try:
        for c in costs:
            patched(fake, np.full((3, 3), c))
    except SR.StallHitStop:
        stopped = True
    return dict(stopped=stopped, meter=meter, calls=calls)


def test_constants_are_the_frozen_values() -> None:
    assert HR.TAU == 0.95 and HR.K == 500 and HR.PH_K == 3000 and HR.SEED_TEST == 0.80
    assert "PH(500)" in HR.RULE and "No low side" in HR.RULE and "0.95" in HR.RULE
    assert HR.EXPECT["held"] == dict(total=24_457_773, mc=11_381_609)
    assert HR.EXPECT["mid"] == dict(total=11_295_158, mc=3_795_130)
    assert HR.EXPECT["grad"] == dict(total=9_175_867, mc=1_078_188)
    assert HR.EXPECT["big"] == dict(total=3_635_256, mc=424_450)


def test_ph_3000_is_absent_from_this_rule() -> None:
    """The rule has one form: k = 500; no 3000, 1000, 250 or low side in the rule's own code."""
    src = Path(HR.__file__).read_text()
    assert "RECORD_KS" not in src and "FORMS" not in src
    assert HR.RECORD_K == HR.K == 500
    # arm J5's tracer is armed with k=500 and low=None, never with 3000
    pick_src = (ROOT / "scripts" / "autoresearch" / "ar5" / "high500_pick.py").read_text()
    assert 'tr.begin(int(H5["k"]), None, tuple(H5["aux"]))' in pick_src
    assert "BAND" not in pick_src and "LOW" not in pick_src


def test_decide_high_only_with_closed_end() -> None:
    assert HR.decide(0.95) == "high" and HR.decide(0.9499999) == "go"
    assert HR.decide(0.0) == "go" and HR.decide(0.5) == "go" and HR.decide(0.7999) == "go"


def test_stops_high_at_0_95_and_the_check_is_the_final_evaluation() -> None:
    costs = [0.9, 0.5] + [0.7] * 10 + [0.4]
    r = _drive(costs, 5, lambda c: 0.96)
    m = r["meter"]
    assert r["stopped"] and not m.on and m.n == 7 and m.stop_kind == "high"
    assert m.checks == [[7, 0.96, True]] and m.continued == 0
    assert float(m.best_R[0, 0]) == 0.5  # the stopped end is the running best, not the later 0.4


def test_never_stops_low_goes_on_under_0_95_counting_the_check() -> None:
    costs = [0.9, 0.5] + [0.7] * 10
    for h in (0.0, 0.5, 0.79, 0.80, 0.88, 0.9499999):
        r = _drive(costs, 5, lambda c, h=h: h)
        m = r["meter"]
        assert not r["stopped"] and m.on and m.n == len(costs) and m.stop_kind == ""
        assert m.checks == [[7, h, False]] and m.continued == 1


def test_stops_at_exactly_0_95_not_just_under() -> None:
    costs = [0.9, 0.5] + [0.7] * 10
    assert _drive(costs, 5, lambda c: 0.95)["stopped"]
    assert not _drive(costs, 5, lambda c: 0.9499999)["stopped"]


def test_rechecks_after_a_gain_and_stops_when_the_new_best_reaches_0_95() -> None:
    costs = [0.9, 0.6] + [0.7] * 6 + [0.4] + [0.8] * 10  # stall at 7, gain at 9, stall at 14

    def hit(c: float) -> float:
        return 0.88 if c == 0.6 else 0.97

    r = _drive(costs, 5, hit)
    m = r["meter"]
    assert r["stopped"] and m.n == 14 and m.stop_kind == "high"
    assert m.checks == [[7, 0.88, False], [14, 0.97, True]] and m.continued == 1
    assert float(m.best_R[0, 0]) == 0.4 and r["calls"] == [0.6, 0.4]


def test_no_second_check_without_a_gain() -> None:
    costs = [0.9, 0.5] + [0.7] * 40
    r = _drive(costs, 5, lambda c: 0.9)
    assert len(r["calls"]) == len(r["meter"].checks) == 1  # one hit evaluation per stall check


def test_first_lowest_cost_is_kept_on_a_tie() -> None:
    costs = [0.9, 0.3, 0.3, 0.3, 0.5, 0.5, 0.5, 0.3, 0.3]
    r = _drive(costs, 3, lambda c: 1.0)
    m = r["meter"]
    assert r["stopped"] and m.best == 0.3 and m.trace == [[1, 0.9], [2, 0.3]]
    assert m.counter.best_n == 2


def test_a_gap_of_exactly_k_does_not_check() -> None:
    costs = [0.9, 0.5, 0.7, 0.7, 0.7, 0.7, 0.4]
    r = _drive(costs, 5, lambda c: 1.0)
    assert not r["stopped"] and r["calls"] == []


def test_stall_of_500_checks_at_the_501st_evaluation_after_the_last_gain() -> None:
    costs = [0.9, 0.5] + [0.7] * 600
    r = _drive(costs, 500, lambda c: 0.97)
    assert r["stopped"]
    assert r["meter"].n == 502  # gain at 2, then 500 evaluations with no gain: check at 502
    ok = _drive([0.9, 0.5] + [0.7] * 499, 500, lambda c: 0.97)
    assert not ok["stopped"] and ok["calls"] == []


def test_walk_events_stops_high_or_never_with_no_low_side() -> None:
    ev = [[500, 0.5, 0.40], [1200, 0.3, 0.97]]
    c = HR.refine_cuts([dict(pos=0, evals=4000, cont_j=1, final_cost=0.2, events={"500": ev})])
    assert c[0]["kind"] == "moved" and c[0]["evals"] == 1200 + 1 + 1  # one check went on
    assert c[0]["w"]["kind"] == "high" and c[0]["cost"] == 0.3
    # a hit ratio far under 0.80 is not a stop here (step 15's two-sided rule would stop)
    c = HR.refine_cuts(
        [dict(pos=0, evals=4000, cont_j=0, final_cost=0.2, events={"500": [[500, 0.5, 0.40]]})]
    )
    assert c[0]["kind"] == "none" and c[0]["evals"] == 4001
    # no 500 event: no check, J's evaluations
    c = HR.refine_cuts([dict(pos=0, evals=4000, cont_j=0, final_cost=0.2, events={"500": []})])
    assert c[0]["kind"] == "none" and c[0]["evals"] == 4000
    # events at other sizes are never read
    c = HR.refine_cuts(
        [dict(pos=0, evals=4000, cont_j=0, final_cost=0.2, events={"1000": [[1000, 0.5, 0.99]]})]
    )
    assert c[0]["kind"] == "none"


def test_a_stop_at_js_own_end_keeps_the_path_and_a_moved_stop_is_marked() -> None:
    refs = [
        dict(
            pos=0,
            evals=5000,
            cont_j=1,
            final_cost=0.30,
            end_hit=0.97,
            events={"500": [[1500, 0.38, 0.96]]},
        ),
        dict(
            pos=1,
            evals=4000,
            cont_j=0,
            final_cost=0.35,
            end_hit=0.90,
            events={"500": [[1500, 0.35, 0.96]]},
        ),
    ]
    c = HR.refine_cuts(refs)
    assert c[0]["kind"] == "moved" and c[1]["kind"] == "after_last"
    assert HR.path_change_refs(refs) == []  # moved at 0.95+, J's own end is also at 0.95+
    refs[0]["end_hit"] = 0.90  # J's end under 0.95: the stop at 0.96 turns a miss into a stop
    assert HR.path_change_refs(refs) == [0]


def test_whole_step_walk_for_a_moved_end_at_0_95_keeps_the_taken_candidate() -> None:
    fits = [
        dict(pos=0, kind="C", hit=0.60, cost=0.40, evals=100),
        dict(pos=0, kind="M", hit=0.97, cost=0.30, evals=5000),
        dict(pos=1, kind="C", hit=0.50, cost=0.45, evals=100),
        dict(pos=1, kind="M", hit=0.70, cost=0.35, evals=4000),
    ]
    refs = [
        dict(
            pos=0,
            evals=5000,
            cont_j=0,
            final_cost=0.30,
            end_hit=0.97,
            events={"500": [[1500, 0.38, 0.96]]},
        )
    ]
    r = SR.walk_step(fits, HR.refine_cuts(refs))
    assert r["stopped"] and r["taken"] == 0 and not r["lost"]


def test_guard_refuses_held2_as_a_sample() -> None:
    for name in ("held2", "held2x"):
        with pytest.raises(SystemExit):
            HG.refuse_sample(name)
    with pytest.raises(SystemExit):
        HG.refuse_sample("other")
    for ok in ("big", "grad", "mid", "held", "held3", ""):
        HG.refuse_sample(ok)


def test_guard_refuses_any_path_that_mentions_held2(tmp_path: Path) -> None:
    (tmp_path / HG.H3_MARKER).write_text("{}")  # held3 open here: held2 stays refused
    for p in (
        tmp_path / "J_held2" / "run.json",
        tmp_path / "passes_held2" / "D4_J" / "pass.json",
        Path("a/b/pass_J_held2.log"),
        tmp_path / "stats_held2.json",
        tmp_path / "x" / "held2_clean.npy",
    ):
        with pytest.raises(SystemExit):
            HG.guarded(p, tmp_path)
    real = tmp_path / "passes_held2"
    real.mkdir()
    link = tmp_path / "link"
    link.symlink_to(real)
    with pytest.raises(SystemExit):
        HG.guarded(link / "x.json", tmp_path)  # a link into held2 is refused as well
    assert HG.guarded(tmp_path / "J_held" / "run.json", tmp_path).name == "run.json"


def test_guard_refuses_held3_logs_until_the_marker_exists(tmp_path: Path) -> None:
    for p in (
        tmp_path / "J_held3" / "run.json",
        tmp_path / "passes_held3" / "D4_J5" / "pass.json",
        tmp_path / "logs" / "pass_J5_held3.log",
        tmp_path / "x" / "held3_clean.npy",
    ):
        with pytest.raises(SystemExit):
            HG.guarded(p, tmp_path)
    assert HG.guarded(tmp_path / "J_held" / "run.json", tmp_path).name == "run.json"
    (tmp_path / HG.H3_MARKER).write_text("{}")
    assert HG.guarded(tmp_path / "J_held3" / "run.json", tmp_path).name == "run.json"


def test_freeze_lists_every_own_script_and_reports_a_problem_without_a_file(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    import high500_freeze as HF

    cur = HF.current()
    for s in HF.OWN + HF.BAND_USED:
        assert f"scripts/autoresearch/ar5/{s}" in cur["files"]
    assert "icenine/reconstructor.py" in cur["library"] and "uv.lock" in cur["library"]
    assert HF.CONSTANTS["k"] == 500 and HF.CONSTANTS["low_side"] == "none"
    monkeypatch.setattr(HF, "FROZEN", tmp_path / "frozen.json")
    assert "does not exist" in HF.problems()


def test_link_script_links_only_the_four_disclosed_samples() -> None:
    sh = str(ROOT / "scripts" / "autoresearch" / "ar5" / "high500_link_cache.sh")
    for name in ("held2", "held3", "other"):
        r = subprocess.run([sh, "--dry-run", name], capture_output=True, text=True)
        assert r.returncode == 2 and "refused" in r.stderr
    src = Path(sh).read_text()
    assert '*held2*) echo "refusing held2' in src  # a source or destination with held2
    assert "realpath" in src and "refusing held2 (resolved)" in src


@pytest.mark.parametrize("name", ["held2", "weird"])
def test_importing_the_common_module_with_a_bad_sample_exits_before_anything_else(
    name: str,
) -> None:
    env = dict(os.environ, RESEED_SAMPLE=name)
    code = "import sys; sys.path.insert(0, %r); import high500_common" % str(
        ROOT / "scripts" / "autoresearch" / "ar5"
    )
    r = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True, env=env)
    assert r.returncode == 2 and (name != "held2" or "held2" in r.stderr)


def test_gate_judges_by_the_fixed_rule_in_a_subprocess() -> None:
    """``high500_gate.judge`` on synthetic results (a subprocess: importing the gate imports the
    study's common module, which sets environment variables)."""
    import json

    code = """
import json, sys
sys.path.insert(0, %r)
import high500_gate as G
import high500_rules as HR
def mk(s, moved=None, pc=None, walk=None, share=0.07, broken=(), fj=0, fn=0, gj=0, gn=0):
    e = HR.PRED1[s]
    return dict(moved_steps=list(range(e["moved"] if moved is None else moved)),
                path_change_steps=list(range(e["path_change"] if pc is None else pc)),
                saving_walk=e["saving"] if walk is None else walk, net_saving_share=share,
                live_broken=list(broken), live_failed_j=fj, live_failed_new=fn,
                live_grid_j=gj, live_grid_new=gn)
def allres(**over):
    r = {s: mk(s) for s in HR.GATE["samples"]}
    r.update(over)
    return r
p1 = {s: dict(ok=True) for s in HR.GATE["samples"]}
out = {}
out["good"] = G.judge(allres(), p1)
out["g0_moved"] = G.judge(allres(held=mk("held", moved=9)), p1)["g0"]
out["g0_saving"] = G.judge(allres(big=mk("big", walk=88145)), p1)["g0"]
out["g0_p1"] = G.judge(allres(), dict(p1, mid=dict(ok=False)))["g0"]
held11 = mk("held", moved=11)
out["g1_held11"] = G.judge(allres(held=held11), p1)["g1"]
out["g1_pc3"] = G.judge(allres(mid=mk("mid", pc=3)), p1)["g1"]
out["g1_grad3"] = G.judge(allres(grad=mk("grad", moved=3)), p1)["g1"]
out["g2_broken"] = G.judge(allres(mid=mk("mid", broken=[466])), p1)["g2"]
out["g2_failed"] = G.judge(allres(mid=mk("mid", fj=0, fn=1)), p1)["g2"]
out["g2_grid"] = G.judge(allres(mid=mk("mid", gj=0, gn=1)), p1)["g2"]
out["g3_mid"] = G.judge(allres(mid=mk("mid", share=0.0549)), p1)["g3"]
out["g3_held_ok"] = G.judge(allres(held=mk("held", share=0.055)), p1)["g3"]
out["g3_grad_free"] = G.judge(allres(grad=mk("grad", share=0.0147)), p1)["g3"]
print(json.dumps(out))
""" % str(
        ROOT / "scripts" / "autoresearch" / "ar5"
    )
    env = dict(os.environ, RESEED_SAMPLE="big")
    r = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True, env=env)
    assert r.returncode == 0, r.stderr[-800:]
    o = json.loads(r.stdout.strip().splitlines()[-1])
    assert o["good"]["passes"] is True
    assert o["g0_moved"] is False and o["g0_saving"] is False and o["g0_p1"] is False
    assert o["g1_held11"] is False and o["g1_pc3"] is False and o["g1_grad3"] is False
    assert o["g2_broken"] is False and o["g2_failed"] is False and o["g2_grid"] is False
    assert o["g3_mid"] is False and o["g3_held_ok"] is True and o["g3_grad_free"] is True


def test_held3_support_changes_no_older_script_and_patches_exact_lines() -> None:
    import band_samples as BS

    for which, name in BS.SCRIPTS.items():
        old = (BS.AR3 / name).read_text().splitlines()
        new = BS.patched_source(which).splitlines()
        assert len(old) == len(new) and sum(1 for a, b in zip(old, new) if a != b) == 1


def test_other_random_draws_are_untouched_by_a_stop() -> None:
    """The meter draws no random number: a global generator is in the same state after a run
    with a stop as after a run without one."""
    state = np.random.get_state()
    _drive([0.9, 0.5] + [0.7] * 10, 5, lambda c: 0.96)
    a = np.random.random()
    np.random.set_state(state)
    _drive([0.9, 0.5] + [0.7] * 10, 0, lambda c: 0.96)
    b = np.random.random()
    assert a == b
    for mod in (HR, BR):
        src = open(mod.__file__).read()
        assert "np.random" not in src and "default_rng" not in src


def test_run_all_script_refuses_bad_arguments_and_names_its_markers() -> None:
    sh = ROOT / "scripts" / "autoresearch" / "ar5" / "high500_run_all.sh"
    for args in (["--bogus"], ["--jobs", "5"], ["--jobs", "0"]):
        r = subprocess.run([str(sh)] + args, capture_output=True, text=True)
        assert r.returncode == 2, args
    src = sh.read_text()
    # no setsid (not on macOS), the marker, the held3 marker, a resume that moves old markers aside
    assert "setsid scripts" not in src.replace("#   (macOS has no setsid", "")
    assert "high500_runs.done" in src and "held3_runs_done.json" in src
    assert 'mv "$f" "$idir/' in src and "held2" in src
    # at most four heavy processes; the stop rule on a cap; the pilot gates the held3 runs
    assert '"$jobs" -gt 4' in src and "STOP RULE" in src and "pilot_held3.done" in src


def test_run_all_dry_run_resume_lists_an_old_failed_marker_to_move(tmp_path: Path) -> None:
    logs = tmp_path / "logs"
    logs.mkdir()
    (logs / "high500_runs.done").write_text("exit=1\nreason=a run failed\n")
    sh = ROOT / "scripts" / "autoresearch" / "ar5" / "high500_run_all.sh"
    env = dict(os.environ)
    env.pop("VIRTUAL_ENV", None)
    chk = subprocess.run(
        ["uv", "run", "python", "scripts/autoresearch/ar5/high500_freeze.py", "--check"],
        capture_output=True,
        text=True,
        cwd=str(ROOT),
        env=env,
    )
    if chk.returncode != 0:  # a later step's library edit (step 16 phase B, AR5 step 20), by design
        tail = chk.stderr.strip().splitlines()[-1] if chk.stderr.strip() else ""
        ok = tail.startswith("freeze check failed: icenine/") and tail.endswith(
            "differs from the frozen hash"
        )
        if not ok:  # any other mismatch or a crash is a real failure
            pytest.fail(f"unexpected freeze result: {chk.stderr[-300:]}")
        pytest.skip("the step 16 freeze disagrees only on a later step's library edit")
    r = subprocess.run(
        [str(sh), "--smoke", "--resume", "--dry-run", "--cache", str(tmp_path)],
        capture_output=True,
        text=True,
        cwd=str(ROOT),
        env=env,
    )
    assert r.returncode == 0, r.stderr[-600:]
    assert "would move" in r.stdout and "old marker: exit=1" in r.stdout
    # without --resume the same old marker is a refusal
    r = subprocess.run(
        [str(sh), "--smoke", "--dry-run", "--cache", str(tmp_path)],
        capture_output=True,
        text=True,
        cwd=str(ROOT),
        env=env,
    )
    assert r.returncode == 2 and "marker" in r.stderr


def test_make_held3_refuses_without_a_committed_freeze_and_lists_three_stages() -> None:
    import high500_make_held3 as MH

    assert [s["name"] for s in MH.stages()] == ["build", "render", "stack"]
    src = Path(MH.__file__).read_text()
    assert "check_or_exit" in src and "shutil.move" in src and "build_state.json" in src
    assert "--workers" in src and '"held3"' in src


def test_judge_rules_match_the_frozen_pass_rules() -> None:
    import high500_judge as HJ

    j = HR.JUDGE
    assert {k: v / 100 for k, v in HJ.FEWER_MIN.items()} == j["r1_total_fewer"]
    assert HJ.MORE_MAX["grad"] / 100 == j["r1_grad_more_max"]
    assert (
        HJ.NET == j["r2_broken_minus_fixed_max"] and HJ.MED == j["r3_median_error_higher_max_deg"]
    )
    assert HJ.HOUSE_MED == j["r5_vs_s8"]["median_error_higher_max_deg"]
    h = HJ.house(100, 200, 5, 0, 0.01, 0.01)
    assert h["passes"] and HJ.house(81, 100, 5, 0, 0.01, 0.01)["a"] is False
    ok = dict(a=True, b=True, c=True, passes=True)
    bad = dict(a=False, b=True, c=True, passes=False)
    r = HJ.rules("mid", 95, 100, 5, 0.02, 0.02, 1, 3, 4, ok, ok)
    assert r["passes"] and r["r1"]
    assert not HJ.rules("mid", 96, 100, 0, 0.02, 0.02, 0, 3, 3, ok, ok)["r1"]  # 4% < 5% on mid
    assert not HJ.rules("mid", 90, 100, 6, 0.02, 0.02, 0, 3, 3, ok, ok)["r2"]
    assert not HJ.rules("mid", 90, 100, 0, 0.0211, 0.02, 0, 3, 3, ok, ok)["r3"]
    assert not HJ.rules("mid", 90, 100, 0, 0.02, 0.02, 2, 3, 3, ok, ok)["r4"]
    assert not HJ.rules("mid", 90, 100, 0, 0.02, 0.02, 0, 3, 5, ok, ok)["r4"]
    # J failing rule 5 means J5 is judged by 1 to 4 only
    assert HJ.rules("mid", 90, 100, 0, 0.02, 0.02, 0, 3, 3, bad, bad)["passes"]
    assert not HJ.rules("mid", 90, 100, 0, 0.02, 0.02, 0, 3, 3, bad, ok)["passes"]
    assert HJ.rules("grad", 100, 100, 0, 0.02, 0.02, 0, 1, 1, ok, ok)["r1"]
    assert not HJ.rules("grad", 101, 100, 0, 0.02, 0.02, 0, 1, 1, ok, ok)["r1"]


def test_the_only_frozen_file_that_differs_is_the_phase_b_docstring_note() -> None:
    """Phase B edited one comment in ``icenine/orientation_search.py`` (the frozen hash is the one
    of the freeze commit). Every other frozen script and library file must still agree."""
    import json

    import high500_freeze as HF

    saved = json.loads(HF.FROZEN.read_text())
    cur = HF.current()
    differ = [
        k for part in ("files", "library") for k, v in cur[part].items() if saved[part].get(k) != v
    ]
    # AR5 step 20 (BoundaryPassRepair) added library files and edited these four: allowed here
    allowed = {
        "icenine/orientation_search.py",
        "icenine/boundary_pass.py",
        "icenine/boundary_repair.py",
        "icenine/config_file.py",
        "icenine/reconstructor.py",
    }
    assert set(differ) <= allowed, differ
