"""AR5 step 17 (task ar5q): the rule family for an early skip of the BFS local optimizations of
arm J5, on synthetic fits and a fake BFS object. Plan "AR5 step 17", phase A step 2: each form
skips or stops only below its line; checks are counted; a skipped fit is rejected and keeps the
stated end; discarded draws leave other fits unchanged; recording evaluations are kept apart;
guards on ``held2`` and on ``held4``.
"""

import json
import os
import subprocess
import sys
from collections import defaultdict
from fractions import Fraction
from pathlib import Path
from typing import Any, Dict, List, Optional

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "scripts" / "autoresearch" / "ar5"))

import fitskip_guard as FG  # noqa: E402
import fitskip_hook as FH  # noqa: E402
import fitskip_rules as FR  # noqa: E402
import order_draws as OD  # noqa: E402

I9 = np.eye(3).ravel().tolist()


def rot_z(deg: float) -> np.ndarray:
    a = np.radians(deg)
    c, s = np.cos(a), np.sin(a)
    return np.array([[c, -s, 0.0], [s, c, 0.0], [0.0, 0.0, 1.0]])


SYM = np.stack([np.eye(3), rot_z(90), rot_z(180), rot_z(270)])


def fit(
    voxel: int = 0,
    kind: str = "neighbour",
    best: float = 0.8,
    start_hit: float = 0.5,
    hit20: float = 0.5,
    evals: int = 61,
    accepted: bool = False,
    n_revisits: int = 0,
    start_R: Optional[Any] = None,
    final_before: bool = False,
) -> Dict[str, Any]:
    return dict(
        voxel=voxel,
        kind=kind,
        best=best,
        start_hit=start_hit,
        hit20=hit20,
        evals=evals,
        accepted=accepted,
        n_revisits=n_revisits,
        start_R=I9 if start_R is None else np.asarray(start_R).ravel().tolist(),
        final_before=final_before,
    )


# ---- the family
def test_family_has_27_rules_and_unique_names() -> None:
    fam = FR.family()
    assert len(fam) == 27
    assert len({r.name for r in fam}) == 27
    assert sum(r.form == "D" for r in fam) == 3
    assert sum(r.form == "S" for r in fam) == 12 and sum(r.form == "C" for r in fam) == 12
    assert FR.rule_by_name("S(0.3,both)").t == 0.3
    with pytest.raises(ValueError):
        FR.rule_by_name("S(0.7,both)")


def test_below_line_is_strict() -> None:
    assert FR.below_line(0.239, 0.8, 0.3)
    assert not FR.below_line(0.24, 0.8, 0.3)  # 0.3 * 0.8 = 0.24: on the line, not under
    assert not FR.below_line(0.9, 0.8, 0.3)
    assert not FR.below_line(0.0, 0.0, 0.3)  # no best yet: nothing is under the line


def test_accepted_test_is_the_callers() -> None:
    assert FR.accepted_test(0.73, 0.8)  # 0.73 / 0.8 = 0.9125
    assert not FR.accepted_test(0.72, 0.8)  # 0.9 exactly is not over 0.9
    assert FR.accepted_test(0.9, 0.5)  # a new best is accepted
    assert not FR.accepted_test(0.0, 0.0)


# ---- replay: S and C
def test_s_skips_only_under_the_line_and_counts_every_check() -> None:
    fits = [
        fit(start_hit=0.23),  # under 0.3 * 0.8 = 0.24: skipped, J5 rejected -> saved
        fit(start_hit=0.24),  # on the line: not skipped
        fit(start_hit=0.5),
        fit(kind="revisit", start_hit=0.1),  # out of scope for a neighbour rule
    ]
    r = FR.replay(FR.rule_by_name("S(0.3,neighbour)"), fits)
    assert (r["n_fits"], r["skipped"], r["wrong"]) == (4, 1, 0)
    assert r["checks"] == 3 and r["saved"] == 61 and r["net"] == 61 - 3
    both = FR.replay(FR.rule_by_name("S(0.3,both)"), fits)
    assert both["checks"] == 4 and both["skipped"] == 2 and both["net"] == 122 - 4
    rev = FR.replay(FR.rule_by_name("S(0.3,revisit)"), fits)
    assert rev["checks"] == 1 and rev["skipped"] == 1


def test_wrong_skip_earns_nothing_and_is_counted() -> None:
    fits = [fit(start_hit=0.1, accepted=True), fit(start_hit=0.1, accepted=False)]
    r = FR.replay(FR.rule_by_name("S(0.3,both)"), fits)
    assert (r["skipped"], r["wrong"], r["saved"], r["checks"], r["net"]) == (2, 1, 61, 2, 59)


def test_final_answer_changes_are_counted() -> None:
    fits = [fit(start_hit=0.1, final_before=True), fit(start_hit=0.1)]
    assert FR.replay(FR.rule_by_name("S(0.3,both)"), fits)["final_changed"] == 1


def test_c_uses_the_hit_after_20_and_saves_evals_beyond_20() -> None:
    fits = [fit(hit20=0.2, start_hit=0.9), fit(hit20=0.5, start_hit=0.0)]
    r = FR.replay(FR.rule_by_name("C(0.3,both)"), fits)
    assert (r["skipped"], r["checks"], r["saved"], r["net"]) == (1, 2, 61 - 20, 41 - 2)


# ---- replay: D
def test_d_skips_a_revisit_near_a_rejected_start_only() -> None:
    near = rot_z(0.05)  # 0.05 degrees from the identity
    far = rot_z(0.4)
    fits = [
        fit(voxel=3, kind="neighbour", start_R=np.eye(3)),  # rejected: its start is remembered
        fit(voxel=3, kind="revisit", n_revisits=1, start_R=near),  # within 0.1: skipped
        fit(voxel=3, kind="revisit", n_revisits=2, start_R=far),  # 0.4 away from both: fitted
        fit(voxel=4, kind="revisit", n_revisits=1, start_R=near),  # other voxel: fitted
    ]
    r = FR.replay(FR.rule_by_name("D(0.1)"), fits, SYM)
    assert (r["skipped"], r["checks"], r["saved"], r["net"]) == (1, 0, 61, 61)
    r25 = FR.replay(FR.rule_by_name("D(0.5)"), fits, SYM)
    assert r25["skipped"] == 2  # the 0.4 degree start is now near one that was rejected


def test_d_never_skips_a_neighbour_fit_and_resets_with_the_counter() -> None:
    fits = [
        fit(voxel=1, kind="neighbour", start_R=np.eye(3)),
        fit(voxel=1, kind="neighbour", start_R=np.eye(3)),  # still a neighbour fit
        fit(voxel=1, kind="revisit", n_revisits=2, start_R=np.eye(3)),  # skipped
        fit(voxel=1, kind="revisit", n_revisits=1, start_R=np.eye(3)),  # counter went down: reset
    ]
    r = FR.replay(FR.rule_by_name("D(0.1)"), fits, SYM)
    assert r["skipped"] == 1


def test_mis_deg_uses_the_symmetry() -> None:
    assert FR.mis_deg(np.eye(3), rot_z(90), SYM) == pytest.approx(0.0, abs=1e-6)
    assert FR.mis_deg(np.eye(3), rot_z(10), SYM) == pytest.approx(10.0, abs=1e-6)
    assert FR.mis_deg(np.eye(3), rot_z(100), SYM) == pytest.approx(10.0, abs=1e-6)


# ---- selection
def _rep(wrong: int, n: int, net: int, fin: int = 0) -> Dict[str, Any]:
    return dict(wrong=wrong, n_fits=n, net=net, final_changed=fin)


def _all(**kw: Any) -> Dict[str, Dict[str, Any]]:
    out = {}
    for s in FR.SAMPLES:
        out[s] = _rep(
            kw.get("wrong", 0), kw.get("n", 50_000), kw.get("net", 1000), kw.get("fin", 0)
        )
    return out


def test_admissible_limits() -> None:
    assert FR.admissible(_all(wrong=50, n=50_000))  # exactly 0.1%
    assert not FR.admissible(_all(wrong=51, n=50_000))
    assert not FR.admissible(_all(fin=3))
    assert FR.admissible(_all(fin=2))
    per = _all()
    per["big"] = _rep(0, 50_000, 1, 3)  # one sample fails (b): the rule is out
    assert not FR.admissible(per)


def test_select_largest_net_then_fewer_wrong_then_form_then_t() -> None:
    res = {
        "S(0.3,both)": _all(net=2000, wrong=10),
        "D(0.1)": _all(net=2000, wrong=10),  # tie on net and wrong: D before S
        "S(0.4,both)": _all(net=1000),
        "C(0.3,both)": _all(net=2000, wrong=5),  # fewer wrong skips beats the form order
        "S(0.5,both)": _all(net=9999, wrong=999),  # not admissible
    }
    out = FR.select(res)
    assert out["pick"] == "C(0.3,both)"
    assert [x["rule"] for x in out["admissible"]][:3] == ["C(0.3,both)", "D(0.1)", "S(0.3,both)"]
    res.pop("C(0.3,both)")
    assert FR.select(res)["pick"] == "D(0.1)"
    res.pop("D(0.1)")
    assert FR.select(res)["pick"] == "S(0.3,both)"
    assert FR.select({"S(0.5,both)": _all(net=9999, wrong=999)})["pick"] is None


def test_select_smaller_t_breaks_the_last_tie() -> None:
    res = {"S(0.4,both)": _all(net=7), "S(0.3,both)": _all(net=7)}
    assert FR.select(res)["pick"] == "S(0.3,both)"


def test_net_share_is_exact_and_gate_g2() -> None:
    per = {s: _rep(0, 1, -(-FR.EXPECT[s]["total"] * 4 // 100)) for s in FR.SAMPLES}
    assert FR.net_share_sum(per) == sum(
        (Fraction(-(-FR.EXPECT[s]["total"] * 4 // 100), FR.EXPECT[s]["total"]) for s in FR.SAMPLES),
        Fraction(0),
    )
    g = FR.gate_g2(per)
    assert g["held"] and g["mid"] and g["held3"] and g["grad"] and g["big"]
    per["held"] = _rep(0, 1, int(FR.EXPECT["held"]["total"] * 0.039))
    assert not FR.gate_g2(per)["held"]


def test_counts_ok_against_the_plan_table() -> None:
    stats = dict(
        n_revisit_attempted=8010,
        n_revisit_rejected=2873,
        n_neighbor_fits=24533,
        n_neighbor_accepted_first=19109,
    )
    assert FR.counts_ok("big", stats, 3_567_098)["ok"]
    assert not FR.counts_ok("big", stats, 3_567_099)["ok"]
    stats["n_revisit_rejected"] = 2874
    assert not FR.counts_ok("big", stats, 3_567_098)["ok"]


def test_plan_numbers_are_in_the_stored_table() -> None:
    assert [FR.EXPECT[s]["total"] for s in FR.SAMPLES] == [
        23_503_385,
        22_939_812,
        10_589_493,
        9_041_686,
        3_567_098,
    ]
    for e in FR.EXPECT.values():
        assert e["rej_revisits"] < e["revisits"] and e["nbr_rej"] < e["nbr"]


# ---- the hook, on a fake BFS object
class _Info:
    def __init__(self, hit: float, conf: float, cost: float) -> None:
        self.hit, self.conf, self.cost = hit, conf, cost


class _Voxel:
    def __init__(self, R: np.ndarray, conf: float = 0.4, cost: float = 5.0, phase: int = 0) -> None:
        self.orientation = R
        self.confidence = conf
        self.cost = cost
        self.phase = phase


class _Rec:
    def __init__(self, idx: int, n_revisits: int = 0) -> None:
        self.voxel_idx = idx
        self.n_revisits = n_revisits
        self.n_evals = 0


class _Recon:
    def __init__(self) -> None:
        self.last_local_optimization_evals = 61
        self.overlap_calls: List[np.ndarray] = []
        self.hit_of = lambda R: 0.5  # the hit ratio of an orientation (a test device)

    def evaluate_overlap(self, R: Any, vertices: Any, phase: int) -> _Info:
        self.overlap_calls.append(np.array(R, np.float64))
        h = float(self.hit_of(np.asarray(R, np.float64)))
        return _Info(h, h, 3.0 - h)


class _Cr:
    def __init__(self) -> None:
        self.calls: List[Any] = []


class _Wrap:
    def __init__(self) -> None:
        self.cr = _Cr()


class _FakeBFS:
    def __init__(self) -> None:
        self.stats: Dict[str, int] = defaultdict(int)
        self.reconstructor = _Recon()
        self.fit_calls: List[Any] = []
        self.fit_hit = 0.8
        self.fit_R = rot_z(1.0)
        self.wrap: Optional[_Wrap] = None

    @staticmethod
    def _ratios(info: _Info) -> Any:
        return info.hit, info.conf

    def _local_fit(
        self, voxel: Any, vertices: Any, rng: Any, rec: Any, phase: str, sigma0: Any = None
    ) -> Any:
        self.fit_calls.append((rec.voxel_idx, phase))
        rec.n_evals += 61
        self.stats[f"n_evals_{phase}"] += 61
        if self.wrap is not None:  # 61 evaluations: costs fall, orientation = z rotation by i
            self.wrap.cr.calls = [(10.0 - 0.1 * i, rot_z(0.01 * i)) for i in range(61)]
        res = type("R", (), dict(orientation=self.fit_R, cost=1.5, overlap_info=None))()
        return res, self.fit_hit, self.fit_hit

    def _expand(
        self, mic: Any, start_idx: int, best_conf: float, rng: Any, solution: Any, **kw: Any
    ) -> None:
        pass

    def _insert_seed(self, mic: Any, voxel_idx: int, queue: Any, *a: Any, **k: Any) -> None:
        queue.append(voxel_idx + 100)


def _hook(rule: Optional[str], tmp_path: Optional[Path] = None, **kw: Any) -> Any:
    bfs = _FakeBFS()
    book = OD.DrawBook()
    wrap = _Wrap() if tmp_path is not None else None
    bfs.wrap = wrap
    r = FR.rule_by_name(rule) if rule else None
    h = FH.FitHook(
        bfs,
        wrap=wrap,
        rule=r,
        sym=SYM,
        out_path=(tmp_path / "rows.jsonl") if tmp_path is not None else None,
        book=book,
        voxel_key=OD.voxel_key,
        **kw,
    )
    return bfs, h, book


def _call(
    bfs: Any, kind: str, idx: int, start: np.ndarray, best: float, n_rev: int = 0, conf: float = 0.4
) -> Any:
    bfs._expand(None, 0, best, None, [])
    rec = _Rec(idx, n_rev)
    vox = _Voxel(start, conf=conf)
    out = bfs._local_fit(vox, np.zeros((3, 3)), None, rec, kind)
    return out, rec


def test_s_skip_holds_the_start_of_a_neighbour_and_counts_the_check() -> None:
    bfs, h, book = _hook("S(0.3,both)")
    bfs.reconstructor.hit_of = lambda R: 0.1  # the start's hit ratio: under 0.3 * 0.8
    (res, hit, conf), rec = _call(bfs, "neighbor", 5, rot_z(3.0), 0.8)
    assert bfs.fit_calls == []  # the fit did not run
    assert np.allclose(res.orientation, rot_z(3.0)) and hit == 0.1 and conf == 0.1
    assert (
        rec.n_evals == 1 and bfs.stats["n_evals_neighbor"] == 1
    )  # one counted evaluation: the check
    assert not FR.accepted_test(hit, 0.8)  # a rejected fit
    assert h.counts["skipped"] == 1 and h.counts["checks"] == 1


def test_s_skip_of_a_revisit_keeps_the_old_end_through_a_negative_confidence() -> None:
    bfs, h, _ = _hook("S(0.3,both)")
    bfs.reconstructor.hit_of = lambda R: 0.1
    (res, hit, conf), rec = _call(bfs, "revisit", 5, rot_z(3.0), 0.8, n_rev=1, conf=0.4)
    assert conf == -1.0 < 0.4  # ``_expand``'s push rule (conf >= old) keeps the old fit
    assert bfs.fit_calls == [] and rec.n_evals == 1


def test_s_no_skip_runs_the_fit_and_still_counts_the_check() -> None:
    bfs, h, _ = _hook("S(0.3,both)")
    bfs.reconstructor.hit_of = lambda R: 0.6
    (res, hit, conf), rec = _call(bfs, "neighbor", 5, rot_z(3.0), 0.8)
    assert bfs.fit_calls == [(5, "neighbor")] and hit == 0.8
    assert rec.n_evals == 61 + 1 and bfs.stats["n_evals_neighbor"] == 62  # the check is on top
    assert h.counts["skipped"] == 0 and h.counts["checks"] == 1


def test_s_on_the_line_is_not_skipped_and_scope_is_respected() -> None:
    bfs, h, _ = _hook("S(0.3,revisit)")
    bfs.reconstructor.hit_of = lambda R: 0.0
    _call(bfs, "neighbor", 1, rot_z(0.0), 0.8)  # a neighbour fit: out of scope, no check
    assert bfs.fit_calls == [(1, "neighbor")] and h.counts["checks"] == 0
    bfs2, h2, _ = _hook("S(0.3,both)")
    bfs2.reconstructor.hit_of = lambda R: 0.24
    _call(bfs2, "neighbor", 1, rot_z(0.0), 0.8)  # 0.24 = 0.3 * 0.8: on the line
    assert bfs2.fit_calls == [(1, "neighbor")]


def test_a_skipped_fit_takes_and_discards_its_aligned_draw() -> None:
    """The kth fit at a voxel gets the same generator whether or not an earlier fit was skipped."""
    key = OD.voxel_key(np.zeros((3, 3)))
    plain = OD.DrawBook()
    first, second = plain.rng(key, "_local_fit"), plain.rng(key, "_local_fit")
    bfs, h, book = _hook("S(0.3,both)")
    bfs.reconstructor.hit_of = lambda R: 0.1
    _call(bfs, "neighbor", 5, rot_z(3.0), 0.8)  # skipped: takes draw 0 and drops it
    assert book.counts[(key, "_local_fit")] == 1
    got = book.rng(key, "_local_fit")  # the next fit at the voxel is draw number 1
    assert got.random() == second.random() and first is not None


def test_d_skip_needs_no_evaluation_and_a_nearby_rejected_start() -> None:
    bfs, h, book = _hook("D(0.25)")
    bfs.fit_hit = 0.1  # the fits are rejected (hit 0.1 against a best of 0.8)
    _call(bfs, "neighbor", 7, rot_z(0.0), 0.8)  # rejected, its start remembered
    assert len(bfs.fit_calls) == 1
    (res, hit, conf), rec = _call(bfs, "revisit", 7, rot_z(0.1), 0.8, n_rev=1)
    assert len(bfs.fit_calls) == 1  # skipped
    assert conf == -1.0 and hit == 0.0 and rec.n_evals == 0  # no evaluation, no cost
    assert h.counts["checks"] == 0 and h.counts["skipped"] == 1
    _call(bfs, "revisit", 7, rot_z(2.0), 0.8, n_rev=2)  # far: fitted
    assert len(bfs.fit_calls) == 2


def test_rule_c_has_no_live_form() -> None:
    with pytest.raises(ValueError):
        _hook("C(0.3,both)")


def test_the_hook_changes_nothing_when_no_rule_acts(tmp_path: Path) -> None:
    bfs, h, book = _hook(None, tmp_path)
    out, rec = _call(bfs, "neighbor", 2, rot_z(2.0), 0.8)
    assert out[1] == 0.8 and rec.n_evals == 61 and bfs.stats["n_evals_neighbor"] == 61
    assert book.n_issued == 0  # the hook takes no draw


def test_recording_row_and_recording_evaluations_apart(tmp_path: Path) -> None:
    bfs, h, _ = _hook(None, tmp_path)
    bfs.reconstructor.hit_of = lambda R: 0.55
    bfs._insert_seed(None, 2, [])  # voxel 102 gets expander 2 (the fake queue is a list)
    out, rec = _call(bfs, "revisit", 102, rot_z(2.0), 0.8, n_rev=3, conf=0.9)
    h.close()
    rows = FH.read_rows(tmp_path / "rows.jsonl")
    assert len(rows) == 1
    r = rows[0]
    assert set(FH.ROW_KEYS) <= set(r)
    assert r["kind"] == "revisit" and r["voxel"] == 102 and r["n_revisits"] == 3
    assert r["best"] == 0.8 and r["evals"] == 61 and r["n_calls"] == 61 and len(r["runbest"]) == 61
    assert r["runbest"] == sorted(r["runbest"], reverse=True)  # a running best never rises
    assert r["start_hit"] == 0.55 and r["hit20"] == 0.55 and r["rec_evals"] == 2
    # the running best after 20 evaluations is the 20th call (costs fall): it was scored
    assert any(np.allclose(c, rot_z(0.19)) for c in bfs.reconstructor.overlap_calls)
    # accepted: 0.8 / 0.8 = 1 > 0.9; kept: yes (accepted)
    assert r["accepted"] is True and r["kept"] is True
    # the counted run is unchanged: the fit's 61 evaluations only, the 2 recording ones are apart
    assert rec.n_evals == 61 and bfs.stats["n_evals_revisit"] == 61
    assert h.counts["rec_evals"] == 2 and h.counts["checks"] == 0


def test_recording_a_rejected_revisit_is_not_kept_when_its_confidence_is_lower(
    tmp_path: Path,
) -> None:
    bfs, h, _ = _hook(None, tmp_path)
    bfs.fit_hit = 0.3
    _call(bfs, "revisit", 9, rot_z(2.0), 0.8, n_rev=1, conf=0.5)  # 0.3 < old 0.5: old end stays
    _call(bfs, "neighbor", 10, rot_z(2.0), 0.8)  # a first fit is always kept
    h.close()
    rows = FH.read_rows(tmp_path / "rows.jsonl")
    assert [r["accepted"] for r in rows] == [False, False]
    assert [r["kept"] for r in rows] == [False, True]


def test_expansion_best_follows_the_fits_of_one_expansion(tmp_path: Path) -> None:
    bfs, h, _ = _hook(None, tmp_path)
    bfs._expand(None, 0, 0.6, None, [])
    rec = _Rec(1)
    bfs._local_fit(_Voxel(rot_z(0.0)), np.zeros((3, 3)), None, rec, "neighbor")  # hit 0.8
    bfs._local_fit(_Voxel(rot_z(0.0)), np.zeros((3, 3)), None, _Rec(2), "neighbor")
    h.close()
    rows = FH.read_rows(tmp_path / "rows.jsonl")
    assert [r["best"] for r in rows] == [0.6, 0.8]  # the best rose with the first fit
    assert [r["pos"] for r in rows] == [0, 1]


# ---- guards
def test_refuse_sample_held2_and_foreign_names(tmp_path: Path) -> None:
    for bad in ("held2", "xheld2", "other"):
        with pytest.raises(SystemExit) as e:
            FG.refuse_sample(bad, gate=tmp_path / "gate.json")
        assert e.value.code == 2
    FG.refuse_sample("held3", gate=tmp_path / "gate.json")
    FG.refuse_sample("", gate=tmp_path / "gate.json")


def test_held4_is_refused_until_the_gate_file_exists(tmp_path: Path) -> None:
    gate = tmp_path / "gate.json"
    with pytest.raises(SystemExit):
        FG.refuse_sample("held4", gate=gate)
    gate.write_text("{}")
    FG.refuse_sample("held4", gate=gate)


def test_guarded_refuses_the_second_held_sample_paths(tmp_path: Path) -> None:
    for name in ("J_held2/run.json", "passes_held2", "pass_J_held2.log", "stats_held2.json"):
        with pytest.raises(SystemExit):
            FG.guarded(tmp_path / name, cache=tmp_path)
    assert FG.guarded(tmp_path / "J5_held3" / "run.json", cache=tmp_path)


def test_guarded_refuses_held4_until_its_marker_exists(tmp_path: Path) -> None:
    p = tmp_path / "J5_held4" / "recon.mic.gz"
    with pytest.raises(SystemExit):
        FG.guarded(p, cache=tmp_path)
    with pytest.raises(SystemExit):
        FG.guarded(tmp_path / "logs" / "S8_held4.log", cache=tmp_path)
    (tmp_path / FG.H4_MARKER).write_text("{}")
    assert FG.guarded(p, cache=tmp_path) == p


def test_common_refuses_held2_before_importing_anything_else() -> None:
    env = dict(os.environ, RESEED_SAMPLE="held2")
    out = subprocess.run(
        [sys.executable, "-c", "import fitskip_common"],
        cwd=str(ROOT / "scripts" / "autoresearch" / "ar5"),
        env=env,
        capture_output=True,
        text=True,
    )
    assert out.returncode == 2
    assert "held2" in out.stderr


def test_common_refuses_held4_before_the_gate() -> None:
    env = dict(os.environ, RESEED_SAMPLE="held4")
    out = subprocess.run(
        [sys.executable, "-c", "import fitskip_common"],
        cwd=str(ROOT / "scripts" / "autoresearch" / "ar5"),
        env=env,
        capture_output=True,
        text=True,
    )
    if (FG.FS_DIR / "gate.json").exists():
        pytest.skip("the gate file exists")
    assert out.returncode == 2 and "gate" in out.stderr


def test_freeze_lists_every_script_and_the_rules() -> None:
    import fitskip_freeze as FF

    cur = FF.current()
    for s in FF.OWN:
        assert f"scripts/autoresearch/ar5/{s}" in cur["files"]
    assert "icenine/reconstructor.py" in cur["library"]
    assert FF.CONSTANTS["rules"] == [r.name for r in FR.family()]
    assert json.dumps(FF.CONSTANTS)  # serialisable
