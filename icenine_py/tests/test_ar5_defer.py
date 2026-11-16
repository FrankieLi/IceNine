"""AR5 step 22 (task ar5v): rule ``DEF``, defer a stalled seed refinement while untried candidates
remain. Tests of the pure rule (the walk of a stored seed step, the seed step's state machine), the
guards on ``held2`` and ``held6`` and the freeze and run-all layers (plan "AR5 step 22", phase A
step 2). Written before the scripts were frozen; nothing here opens a sample."""

import itertools
import subprocess
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional

import pytest

ROOT = Path(__file__).resolve().parent.parent
AR5 = ROOT / "scripts" / "autoresearch" / "ar5"
sys.path.insert(0, str(AR5))

import defer_guard as DG  # noqa: E402
import defer_rules as DR  # noqa: E402

TAU = DR.TAU


# ---- a fake seed step: scripted fits -------------------------------------------------------------
class Fake:
    """Scripted candidates. ``cma[p]`` = (hit, cost); ``mc[p]`` = dict(pause_at=(n, hit_at_check)
    or None, final=(hit, cost), evals=J5's total evaluations of the refinement). A refinement with
    ``pause_at`` goes on at its check (hit under 0.95) and ends at ``final``; one without ends at
    ``final`` with no check that goes on."""

    def __init__(self, cma: List[Any], mc: List[Optional[Dict[str, Any]]]) -> None:
        self.cma, self.mc = cma, mc
        self.log: List[Any] = []

    def run_c(self, p: int) -> Dict[str, Any]:
        self.log.append(("c", p))
        h, c = self.cma[p]
        return dict(hit=h, cost=c, evals=201, R=("cma", p))

    def run_m(self, p: int, allow: bool) -> Dict[str, Any]:
        spec = self.mc[p]
        assert spec is not None, f"no MC scripted for {p}"
        self.log.append(("m", p, allow))
        if allow and spec.get("pause_at"):
            n, h = spec["pause_at"]
            assert h < TAU
            return dict(
                hit=h,
                cost=spec["cost_at_pause"],
                evals=n + 1,
                R=("part", p),
                paused=True,
                state=dict(n=n),
            )
        h, c = spec["final"]
        return dict(hit=h, cost=c, evals=spec["evals"], R=("mc", p))

    def resume_m(self, p: int, prior: Dict[str, Any]) -> Dict[str, Any]:
        spec = self.mc[p]
        self.log.append(("r", p))
        assert prior["paused"] and prior["state"]["n"] == spec["pause_at"][0]
        h, c = spec["final"]
        return dict(hit=h, cost=c, evals=spec["evals"], R=("mc", p))


def mc(
    final_hit: float, final_cost: float, evals: int, pause_at: Any = None, cost_at: float = 1.0
) -> Dict[str, Any]:
    return dict(
        final=(final_hit, final_cost), evals=evals, pause_at=pause_at, cost_at_pause=cost_at
    )


def step(f: Fake, n: int = 5, defer: bool = True) -> Dict[str, Any]:
    return DR.run_seed_step(n, f.run_c, f.run_m, f.resume_m, TAU, defer=defer)


def counted(res: Dict[str, Any]) -> int:
    return sum(k["evals"] for k in res["kept"])


def test_no_check_that_goes_on_means_no_pause_and_the_step_is_j5s() -> None:
    f = Fake([(0.5, 3.0), (0.97, 1.0)], [mc(0.6, 2.0, 900), None])
    r = step(f, 2)
    assert r["paused"] == [] and r["stop"] == 1 and r["found"] == "tried"
    assert [c[0] for c in r["calls"]] == ["C", "M", "C"]
    assert counted(r) == 201 + 900 + 201


def test_a_pause_only_at_a_check_that_goes_on_and_never_at_the_last_candidate() -> None:
    cma = [(0.5, 3.0)] * 5
    pa = [mc(0.6, 2.0, 2000, pause_at=(500, 0.6)) for _ in range(5)]
    r = step(Fake(cma, pa))
    assert r["paused"] == [0, 1, 2, 3]  # candidate 4 is last: never paused
    m4 = [k for k in r["kept"] if k["pos"] == 4][0]["m"]
    assert not m4.get("paused")
    # all paused ones resumed in candidate order, then every end taken as J5's
    assert r["resumed"] == [0, 1, 2, 3] and r["found"] == "none" and r["stop"] is None
    assert [c for c in r["calls"] if c[0] == "R"] == [("R", 0), ("R", 1), ("R", 2), ("R", 3)]


def test_a_refinement_that_stops_at_its_first_check_is_not_paused() -> None:
    # J5's check stops it (hit at 0.95 or more): the fit comes back complete, found at once
    f = Fake([(0.5, 3.0)], [mc(0.97, 1.0, 600)])
    r = step(f, 5)
    assert r["paused"] == [] and r["stop"] == 0 and r["found"] == "tried"
    assert counted(r) == 201 + 600


def test_defer_off_pauses_nothing() -> None:
    cma = [(0.5, 3.0)] * 5
    pa = [mc(0.6, 2.0, 2000, pause_at=(500, 0.6)) for _ in range(5)]
    f = Fake(cma, pa)
    r = step(f, 5, defer=False)
    assert r["paused"] == [] and r["resumed"] == []
    assert all(c[2] is False for c in f.log if c[0] == "m")
    assert counted(r) == 5 * (201 + 2000)


def test_type_a_found_at_a_later_candidate_than_every_paused_refinement() -> None:
    cma = [(0.5, 3.0), (0.5, 3.0), (0.5, 3.0)]
    pa = [mc(0.6, 2.0, 2000, pause_at=(500, 0.6), cost_at=2.5), mc(0.97, 1.0, 700), None]
    r = step(Fake(cma, pa), 3)
    assert r["stop"] == 1 and r["found"] == "tried" and r["paused"] == [0] and r["resumed"] == []
    k0 = r["kept"][0]
    assert k0["paused"] and k0["evals"] == 201 + 501  # the part before the pause, check counted
    assert k0["e"]["hit"] < TAU
    assert counted(r) == 201 + 501 + 201 + 700


def test_resume_goes_in_candidate_order_and_stops_at_the_first_end_at_the_threshold() -> None:
    cma = [(0.5, 3.0)] * 4
    # candidates 0 and 1 are paused; 2 and 3 end under 0.95; resumed 0 ends under, resumed 1 reaches
    pa = [
        mc(0.7, 2.0, 1500, pause_at=(500, 0.6)),
        mc(0.96, 1.0, 1800, pause_at=(500, 0.6)),
        mc(0.8, 2.0, 400),
        mc(0.8, 2.0, 400),
    ]
    r = step(Fake(cma, pa), 4)
    assert r["paused"] == [0, 1] and r["resumed"] == [0, 1]
    assert r["stop"] == 1 and r["found"] == "resume"
    assert [c[0] for c in r["calls"]][-2:] == ["R", "R"]
    assert counted(r) == 4 * 201 + 1500 + 1800 + 400 + 400  # each fit counted once, as J5


def test_a_later_paused_refinement_is_not_resumed_when_an_earlier_one_reaches() -> None:
    cma = [(0.5, 3.0)] * 3
    pa = [
        mc(0.97, 1.0, 1500, pause_at=(500, 0.6)),  # a type C step: J5 stopped here
        mc(0.7, 2.0, 1500, pause_at=(500, 0.6)),
        mc(0.7, 2.0, 400),
    ]
    f = Fake(cma, pa)
    r = step(f, 3)
    assert r["stop"] == 0 and r["found"] == "resume" and r["resumed"] == [0]
    assert [c for c in f.log if c[0] == "r"] == [("r", 0)]
    k1 = [k for k in r["kept"] if k["pos"] == 1][0]
    assert k1["paused"] and k1["e"]["hit"] < TAU  # never resumed: never the taken end


def test_a_c_step_whose_later_candidate_reaches_moves_the_end() -> None:
    cma = [(0.5, 3.0)] * 3
    pa = [mc(0.97, 1.0, 1500, pause_at=(500, 0.6)), mc(0.96, 1.5, 900), None]
    r = step(Fake(cma, pa), 3)
    assert r["stop"] == 1 and r["found"] == "tried" and r["resumed"] == []


def test_every_evaluation_is_counted_once_in_any_order() -> None:
    cma = [(0.5, 3.0)] * 5
    for flags in itertools.product([0, 1], repeat=4):
        pa: List[Optional[Dict[str, Any]]] = []
        for p in range(5):
            if flags[p % 4] and p < 4:
                pa.append(mc(0.8, 2.0, 1200 + p, pause_at=(500, 0.6)))
            else:
                pa.append(mc(0.8, 2.0, 1200 + p))
        f = Fake(cma, pa)
        r = step(f)
        j5 = step(Fake(cma, pa), defer=False)
        # a none step (no end at 0.95 anywhere): exactly J5's count and ends
        assert counted(r) == counted(j5)
        assert [k["e"]["R"] for k in r["kept"]] == [k["e"]["R"] for k in j5["kept"]]
        assert r["stop"] is None and j5["stop"] is None


# ---- the walk of a stored step (the replay's arithmetic) ----------------------------------------
def tr(voxel: int, evals: int, checks: List[List[Any]]) -> Dict[str, Any]:
    return dict(voxel=voxel, evals=evals, checks=checks)


def entry(stop: int, evals: int, voxel: int = 7) -> Dict[str, Any]:
    return dict(voxel=voxel, stop=stop, evals=evals)


def test_walk_type_a_exact_saving() -> None:
    # candidate 0 stalls at 632 and goes on (hit 0.69), 1100 evaluations in all; found at 1
    w = DR.walk_step(entry(1, 3000), {0: tr(7, 1100, [[632, 0.69, False]]), 1: tr(7, 500, [])})
    assert w["type"] == "A" and w["paused"] == {0: 1100 - 633}
    assert w["saving"] == 467 and w["expected_evals"] == 3000 - 467


def test_walk_check_that_stops_is_no_pause() -> None:
    w = DR.walk_step(entry(0, 900), {0: tr(7, 700, [[500, 0.97, True]])})
    assert w["type"] == "-" and w["saving"] == 0 and w["expected_evals"] == 900


def test_walk_last_candidate_never_paused_and_type_c_and_n() -> None:
    n_pos = 5
    c = DR.walk_step(
        entry(1, 3000),
        {0: tr(7, 1000, [[600, 0.6, False]]), 1: tr(7, 1500, [[500, 0.7, False]])},
        n_pos,
    )
    assert c["type"] == "C" and c["saving"] == 1000 - 601 and c["expected_evals"] is None
    last = DR.walk_step(
        entry(4, 4000), {4: tr(7, 1500, [[500, 0.7, False], [900, 0.97, True]])}, n_pos
    )
    assert last["type"] == "-" and last["paused"] == {}
    n = DR.walk_step(entry(-1, 9000), {0: tr(7, 1000, [[600, 0.6, False]])}, n_pos)
    assert n["type"] == "N" and n["expected_evals"] == 9000 and n["saving"] == 0
    assert n["cand_cost"][0] == 201 + 601 and n["cand_cost"][1] == 201


def test_size_sample_bounds_and_plan_match() -> None:
    steps = [
        dict(type="A", saving=400, stop=1, cand_cost=[]),
        dict(type="C", saving=100, stop=2, cand_cost=[]),
        dict(type="N", saving=0, stop=-1, cand_cost=[300, 500]),
    ]
    r = DR.size_sample(steps, 10_000)
    assert (r["a_steps"], r["c_steps"], r["n_steps"]) == (1, 1, 1)
    assert r["a_saving"] == 400 and r["a_pct"] == pytest.approx(4.0)
    assert r["mean_cand"] == 400 and r["extra_hi"] == (4 - 2) * 400
    assert r["extra_lo"] == 201
    assert r["lo"] == pytest.approx(100 * (500 - 800) / 10_000)
    assert r["hi"] == pytest.approx(100 * (500 - 201) / 10_000)


def test_the_plan_table_is_in_the_rules() -> None:
    assert set(DR.PLAN_TABLE) == set(DR.DISCLOSED)
    t = DR.PLAN_TABLE["held"]
    assert (t["a_steps"], t["c_steps"], t["n_steps"]) == (1010, 238, 512)
    assert (
        DR.matches_plan(
            "held",
            dict(
                total=19_457_817,
                a_steps=1010,
                a_pct=12.28,
                c_steps=238,
                n_steps=512,
                lo=10.2,
                hi=14.1,
            ),
        )
        == []
    )
    assert DR.matches_plan(
        "held",
        dict(
            total=19_457_817, a_steps=1011, a_pct=12.28, c_steps=238, n_steps=512, lo=10.2, hi=14.1
        ),
    )


def test_decision_rules_are_the_plans() -> None:
    g, j = DR.GATE, DR.JUDGE
    assert g["g3_fewer_min"] == 0.07 and g["g3_broken_minus_fixed_max"] == 3
    assert g["g3_median_error_higher_max_deg"] == 0.001 and g["g3_grid_more_max"] == 2
    assert j["r1_fewer_min"]["grad"] == 0.015 and j["r1_fewer_min"]["big"] == 0.03
    assert j["r1_fewer_min"]["held"] == 0.07 and j["r1_fewer_min"]["mid"] == 0.07
    assert j["r4_vs_s8"] == dict(
        fewer=0.20, broken_minus_fixed_max=5, median_error_higher_max_deg=0.005
    )
    assert DR.CAPS_MIN == dict(S8=480, rest=300)
    assert (DR.TAU, DR.K, DR.N_CAND) == (0.95, 500, 5)


# ---- guards --------------------------------------------------------------------------------------
def test_held2_is_refused_everywhere(tmp_path: Path) -> None:
    with pytest.raises(SystemExit) as e:
        DG.refuse_sample("held2")
    assert e.value.code == 2
    for p in ("a/held2/b", "J_held2", "stats_held2.json", "held2_clean.npy"):
        with pytest.raises(SystemExit):
            DG.guarded(p, tmp_path)
    with pytest.raises(SystemExit):
        DG.refuse_sample("held9")
    for s in DR.DISCLOSED:
        DG.refuse_sample(s, tmp_path)


def test_held6_is_closed_until_the_gate_and_then_until_its_runs_end(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    (tmp_path / "logs").mkdir()
    with pytest.raises(SystemExit):
        DG.refuse_sample("held6", tmp_path)
    (tmp_path / "logs" / "gate_passed").write_text("ok\n")
    DG.refuse_sample("held6", tmp_path)
    p = tmp_path / "J5D_held6" / "run.json"
    with pytest.raises(SystemExit):
        DG.guarded(p, tmp_path)
    with pytest.raises(SystemExit):
        DG.guarded(tmp_path / "logs" / "J5D_held6.log", tmp_path)
    # nothing of a disclosed sample is closed
    assert DG.guarded(tmp_path / "J5D_held5" / "run.json", tmp_path)
    # a producer job opens its own inputs by name only
    monkeypatch.setenv("SCREEN_PRODUCER", "1")
    assert DG.guarded(tmp_path / "passes_held6" / "before_J5.npy", tmp_path)
    with pytest.raises(SystemExit):
        DG.guarded(tmp_path / "J5D_held6" / "run.json", tmp_path)
    with pytest.raises(SystemExit):
        DG.guarded(tmp_path / "x" / "gradient_truth.mic" / "held6", tmp_path)
    monkeypatch.delenv("SCREEN_PRODUCER")
    (tmp_path / "logs" / "defer_runs.done").write_text("exit=0\n")
    assert DG.guarded(p, tmp_path) == p
    with pytest.raises(SystemExit):  # held2 stays closed whatever the marker says
        DG.guarded(tmp_path / "held2" / "x", tmp_path)


# ---- freeze and scripts --------------------------------------------------------------------------
def test_freeze_lists_exist_pin_files_and_hash_no_data() -> None:
    import defer_freeze as DF

    for name in DF.OWN:
        assert (AR5 / name).exists(), name
    cur = DF.current()
    assert "scripts/autoresearch/ar5/defer_rules.py" in cur["files"]
    assert "scripts/autoresearch/ar5/defer_pick.py" in cur["files"]
    assert "scripts/autoresearch/ar5/screen_rules.py" in cur["files"]
    assert "icenine/boundary_pass.py" in cur["library"] and "uv.lock" in cur["library"]
    assert all(len(v) == 64 for part in cur.values() for v in part.values())  # sha256, not names
    assert not any("held6" in k for k in cur["files"])
    assert all(
        k.startswith(("scripts/", "icenine/")) or k in ("pyproject.toml", "uv.lock")
        for part in cur.values()
        for k in part
    )  # scripts and library only: no data file is hashed (that would read it)
    assert DF.CONSTANTS["plan_table"] == DR.PLAN_TABLE and DF.CONSTANTS["caps_min"] == DR.CAPS_MIN


def test_the_appliers_import_the_common_layer_first() -> None:
    for src in ("defer_pick.py", "defer_run.py", "defer_launch.py"):
        text = (AR5 / src).read_text()
        assert "import defer_common" in text, src
        assert (
            text.index("import defer_common") < text.index("import high500_")
            if "import high500_" in text
            else True
        )
    common = (AR5 / "defer_common.py").read_text()
    # our freeze first, then the older checks are switched off, then step 21's set-up
    assert common.index("DF.problems()") < common.index("_SCF.problems = lambda")
    assert common.index("_SCF.problems = lambda") < common.index("import screen_common")
    assert common.index("DG.refuse_sample(SAMPLE)") < common.index("DF.problems()")


def test_the_pick_never_pauses_the_last_candidate_and_counts_the_rerun_apart() -> None:
    text = (AR5 / "defer_pick.py").read_text()
    assert "allow_pause" in text and "rerun_evals" in text and "DeferRerunMismatch" in text
    rules = (AR5 / "defer_rules.py").read_text()
    assert "p < n_pos - 1" in rules


def test_the_run_all_script_parses_and_refuses_the_rest_stage_before_the_gate(
    tmp_path: Path,
) -> None:
    p = AR5 / "defer_run_all.sh"
    if not p.exists():
        pytest.skip("defer_run_all.sh not written yet")
    assert subprocess.run(["bash", "-n", str(p)], capture_output=True, text=True).returncode == 0
    r = subprocess.run(
        ["bash", str(p), "--dry-run", "--stage", "rest", "--cache", str(tmp_path)],
        capture_output=True,
        text=True,
    )
    assert r.returncode == 2 and "gate_passed" in r.stderr
    r = subprocess.run(["bash", str(p), "--stage", "nonsense"], capture_output=True, text=True)
    assert r.returncode == 2
    text = p.read_text()
    assert "defer_gate.done" in text and "defer_runs.done" in text and "gate_passed" in text
    assert "held2" in text  # the guard
