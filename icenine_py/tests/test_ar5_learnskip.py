"""AR5 step 18 (task ar5r): the learned skip of BFS local optimizations, on synthetic fits and a
fake BFS object. Plan "AR5 step 18", phase A step 2: each feature on hand-made rows; the JSON
evaluator equals scikit-learn; a fit is skipped only under the threshold; one counted check a fit;
a skipped fit is rejected and keeps the stated end; discarded draws; the ``held2`` and ``held4``
guards.
"""

import json
import sys
from collections import defaultdict
from fractions import Fraction
from pathlib import Path
from typing import Any, Dict, List, Optional

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "scripts" / "autoresearch" / "ar5"))

import fitskip_rules as FR  # noqa: E402
import learnskip_features as LFE  # noqa: E402
import learnskip_guard as LG  # noqa: E402
import learnskip_hook as LH  # noqa: E402
import learnskip_model as LM  # noqa: E402
import learnskip_rules as LR  # noqa: E402
import order_draws as OD  # noqa: E402


def rot_z(deg: float) -> np.ndarray:
    a = np.radians(deg)
    c, s = np.cos(a), np.sin(a)
    return np.array([[c, -s, 0.0], [s, c, 0.0], [0.0, 0.0, 1.0]])


SYM4 = np.stack([np.eye(3), rot_z(90), rot_z(180), rot_z(270)])
NAMES = LFE.FEATURES


def feat(st: LFE.FeatureState, **kw: Any) -> Dict[str, float]:
    a: Dict[str, Any] = dict(
        voxel=1,
        kind="neighbour",
        n_revisits=0,
        best=0.8,
        start_R=np.eye(3),
        start_hit=0.4,
        start_cost=5.0,
        start_conf=0.3,
        parent=-1,
        pos=0,
    )
    a.update(kw)
    return dict(zip(NAMES, st.features(**a)))


# ---- the 16 features on hand-made rows
def test_there_are_16_features_in_the_planned_order() -> None:
    assert len(NAMES) == 16 == LFE.N_FEATURES == LR.FEATURE_COUNT
    assert NAMES[0] == "kind" and NAMES[2] == "best" and NAMES[15] == "place"


def test_features_of_a_first_neighbour_fit() -> None:
    f = feat(LFE.FeatureState(SYM4), pos=3)
    assert f["kind"] == 0.0 and f["n_revisits"] == 0.0 and f["best"] == 0.8
    assert f["start_hit"] == 0.4 and f["start_cost"] == 5.0 and f["start_conf"] == 0.3
    assert f["start_over_best"] == pytest.approx(0.5)
    assert f["dist_rejected"] == 10.0  # none rejected: 10
    assert f["parent_unknown"] == 1.0 and f["parent_hit"] == 0.8  # best if unknown
    assert f["parent_cost"] == -1.0 and f["start_over_parent"] == 1.0
    assert f["own_hit"] == -1.0 and f["fits_so_far"] == 0.0 and f["rejected_so_far"] == 0.0
    assert f["place"] == 3.0


def test_best_zero_gives_ratio_zero() -> None:
    assert feat(LFE.FeatureState(SYM4), best=0.0)["start_over_best"] == 0.0


def test_parent_features_come_from_the_parents_kept_fit_end() -> None:
    st = LFE.FeatureState(SYM4)
    st.update(9, np.eye(3), True, True, 0.7, 2.5)  # parent 9: kept end hit 0.7, cost 2.5
    f = feat(st, voxel=1, parent=9, start_hit=0.35)
    assert f["parent_unknown"] == 0.0 and f["parent_hit"] == 0.7 and f["parent_cost"] == 2.5
    assert f["start_over_parent"] == pytest.approx(0.5)
    st.update(8, np.eye(3), True, True, 0.0, 1.0)  # a parent whose kept hit is 0: ratio 1
    assert feat(st, parent=8)["start_over_parent"] == 1.0


def test_a_skipped_fit_or_an_unkept_fit_is_not_a_kept_end() -> None:
    st = LFE.FeatureState(SYM4)
    st.update(4, np.eye(3), False, False, None, None)  # skipped: no fit end
    assert feat(st, voxel=4, kind="revisit", n_revisits=1)["own_hit"] == -1.0
    st.update(4, np.eye(3), False, False, 0.9, 1.0)  # a fit end that the push rule dropped
    assert feat(st, voxel=4, kind="revisit", n_revisits=2)["own_hit"] == -1.0
    st.update(4, np.eye(3), False, True, 0.31, 1.0)
    assert feat(st, voxel=4, kind="revisit", n_revisits=3)["own_hit"] == 0.31
    f = feat(st, voxel=5, parent=4)  # an unrelated voxel; parent 4's kept end
    assert f["parent_hit"] == 0.31


def test_counts_of_fits_and_rejected_fits_at_a_voxel() -> None:
    st = LFE.FeatureState(SYM4)
    st.update(2, np.eye(3), False, True, 0.1, 1.0)
    st.update(2, np.eye(3), True, True, 0.8, 1.0)
    st.update(2, np.eye(3), False, False, None, None)
    f = feat(st, voxel=2, kind="revisit", n_revisits=1)
    assert f["fits_so_far"] == 3.0 and f["rejected_so_far"] == 2.0
    assert feat(st, voxel=3)["fits_so_far"] == 0.0


def test_distance_to_rejected_starts_uses_symmetry_cap_and_reset() -> None:
    st = LFE.FeatureState(SYM4)
    st.update(7, rot_z(0.0), False, True, 0.1, 1.0)  # a rejected start at the identity
    st.update(7, rot_z(5.0), False, True, 0.1, 1.0)  # another at 5 degrees
    f = feat(st, voxel=7, kind="revisit", n_revisits=1, start_R=rot_z(92.0))
    assert f["dist_rejected"] == pytest.approx(2.0, abs=1e-6)  # 92 is 2 degrees from 90 = identity
    f = feat(st, voxel=7, kind="revisit", n_revisits=2, start_R=rot_z(40.0))
    assert f["dist_rejected"] == pytest.approx(10.0)  # 35 degrees away: capped at 10
    f = feat(st, voxel=7, kind="revisit", n_revisits=2, start_R=rot_z(1.0))
    assert f["dist_rejected"] == pytest.approx(10.0)  # count did not go up: reset, none left
    # an accepted fit's start is not a rejected start
    st2 = LFE.FeatureState(SYM4)
    st2.update(7, rot_z(0.0), True, True, 0.9, 1.0)
    assert feat(st2, voxel=7, kind="revisit", n_revisits=1)["dist_rejected"] == 10.0


def test_min_misorientation_agrees_with_step_17s_function() -> None:
    rng = np.random.default_rng(1)
    for _ in range(5):
        a, b = rot_z(rng.uniform(0, 360)), rot_z(rng.uniform(0, 360))
        assert LFE.min_misorientation_deg(a, [b], SYM4) == pytest.approx(
            FR.mis_deg(a, b, SYM4), abs=1e-9
        )
    assert LFE.min_misorientation_deg(np.eye(3), [], SYM4) == 99.0


def test_cubic_rotations_are_the_librarys_for_a_cubic_crystal() -> None:
    from icenine import ray_indexing as RI
    from icenine.symmetry import create_cubic_symmetry

    lib = RI.proper_rotations(create_cubic_symmetry(3.61))
    mine = LFE.cubic_rotations()
    assert len(lib) == len(mine) == 24
    for m in mine:
        assert any(np.abs(m - k).max() < 1e-9 for k in lib)


def test_rows_to_arrays_gives_the_same_features_as_a_walk_with_the_state() -> None:
    rows = []
    for i in range(6):
        rows.append(
            dict(
                voxel=i % 3,
                kind="neighbour" if i < 3 else "revisit",
                n_revisits=0 if i < 3 else i - 2,
                best=0.8,
                start_R=rot_z(7.0 * i).ravel().tolist(),
                start_hit=0.1 * i,
                start_cost=3.0 - 0.1 * i,
                start_conf=0.2,
                parent=-1 if i == 0 else 0,
                pos=i,
                accepted=i % 2 == 0,
                kept=True,
                end_hit=0.5 + 0.05 * i,
                end_cost=2.0,
                evals=61,
            )
        )
    out = LFE.rows_to_arrays(rows, SYM4)
    st = LFE.FeatureState(SYM4)
    for j, r in enumerate(rows):
        x = st.features(
            r["voxel"],
            r["kind"],
            r["n_revisits"],
            r["best"],
            r["start_R"],
            r["start_hit"],
            r["start_cost"],
            r["start_conf"],
            r["parent"],
            r["pos"],
        )
        assert np.array_equal(out["X"][j], x)
        st.update(r["voxel"], r["start_R"], r["accepted"], r["kept"], r["end_hit"], r["end_cost"])
    assert out["accepted"].tolist() == [True, False] * 3 and out["revisit"].sum() == 3


# ---- the model: JSON evaluator equals scikit-learn
def _toy(n: int = 4000) -> Any:
    rng = np.random.default_rng(0)
    X = rng.normal(size=(n, LFE.N_FEATURES))
    y = (X[:, 0] + 0.5 * X[:, 3] * X[:, 5] + 0.3 * rng.normal(size=n)) > 0.2
    return X, y


def test_json_evaluator_equals_scikit_learn(tmp_path: Path) -> None:
    X, y = _toy()
    m = LM.train(X, y)
    model = LM.export(m, threshold=0.25, extra=dict(features=list(NAMES)))
    p = tmp_path / "m.json"
    LM.dump(model, p)
    pm = LM.PlainModel.load(p)
    assert pm.threshold == 0.25 and len(pm.trees) == m.n_iter_
    got = pm.score_many(X[:500])
    want = m.predict_proba(X[:500])[:, 1]
    assert np.abs(got - want).max() < 1e-9


def test_model_spec_is_the_frozen_one() -> None:
    from sklearn.ensemble import HistGradientBoostingClassifier

    X, y = _toy(300)
    m = LM.train(X, y)
    assert isinstance(m, HistGradientBoostingClassifier)
    assert (m.max_iter, m.learning_rate, m.random_state) == (300, 0.1, 0)
    assert m.max_leaf_nodes == 31 and m.min_samples_leaf == 20 and m.l2_regularization == 0.0
    assert m.max_bins == 255 and m.early_stopping == "auto" and m.validation_fraction == 0.1


def test_skip_only_strictly_under_the_threshold() -> None:
    pm = LM.PlainModel(dict(baseline=0.0, threshold=0.5, trees=[]))
    assert pm.score([0.0] * 16) == 0.5
    assert not pm.skips(0.5) and pm.skips(0.4999)
    assert not LM.PlainModel(dict(baseline=0.0, threshold=None, trees=[])).skips(0.0)
    assert not LR.skip_vector(np.array([0.0, 1.0]), None).any()
    assert LR.skip_vector(np.array([0.1, 0.5, 0.9]), 0.5).tolist() == [True, False, False]


# ---- accounting and the threshold procedure
def sample(name: str, n: int = 1000, seed: int = 0, total: Optional[int] = None) -> LR.SampleData:
    rng = np.random.default_rng(seed)
    acc = rng.random(n) < 0.4
    d = LR.SampleData(name, np.zeros((n, 16)), acc, np.full(n, 61), np.zeros(n, bool))
    if total:
        d.total = total
    return d


def test_account_counts_wrong_skips_saved_and_every_check() -> None:
    d = sample("big", 10)
    d.accepted = np.array([True, False, False, True, False] * 2)
    skip = np.array([True, True, False, False, True] + [False] * 5)
    a = LR.account(skip, d)
    assert a["wrong"] == 1 and a["skipped"] == 3 and a["saved"] == 2 * 61
    assert a["checks"] == 10 and a["net"] == 2 * 61 - 10  # a wrong skip earns nothing
    assert a["share"] == Fraction(2 * 61 - 10, d.total)
    d.final_before = np.array([False, True] + [False] * 8)
    assert LR.account(skip, d)["changed"] == 1  # a skipped J5-rejected fit that was final
    d.final_before = np.array([True] + [False] * 9)  # the wrong skip's end is not "changed"
    assert LR.account(skip, d)["changed"] == 0


def test_limits_are_a_tenth_of_a_percent_and_two() -> None:
    assert LR.within_limits(dict(n=1000, wrong=1, changed=2))
    assert not LR.within_limits(dict(n=1000, wrong=2, changed=0))
    assert not LR.within_limits(dict(n=1000, wrong=0, changed=3))
    assert LR.within_limits(dict(n=72458, wrong=72, changed=0))
    assert not LR.within_limits(dict(n=72458, wrong=73, changed=0))


def test_threshold_procedure_picks_the_largest_net_within_the_limits() -> None:
    names = ("big", "mid")
    data = {s: sample(s, 2000, i) for i, s in enumerate(names)}
    oof: Dict[str, np.ndarray] = {}
    rng = np.random.default_rng(5)
    for s in names:  # a score that ranks accepted fits above rejected ones, with some mixing
        oof[s] = np.where(data[s].accepted, 0.6, 0.2) + 0.3 * rng.random(2000)
    out = LR.choose_threshold(oof, data)
    t = out["threshold"]
    assert t is not None and out["n_candidates"] > 10
    for s in names:  # the pick is inside the limits
        assert LR.within_limits(LR.account(oof[s] < t, data[s]))
    # no other candidate inside the limits has a larger sum of shares
    best = sum(LR.account(oof[s] < t, data[s])["share"] for s in names)
    grid = np.unique(np.quantile(np.concatenate(list(oof.values())), LR.QUANTILES))
    for g in grid:
        accs = [LR.account(oof[s] < g, data[s]) for s in names]
        if all(LR.within_limits(a) for a in accs):
            assert sum(a["share"] for a in accs) <= best
    assert out["share_sum"] == float(best)


def test_no_candidate_inside_the_limits_means_no_skip() -> None:
    d = sample("big", 1000)
    out = LR.choose_threshold(
        {"big": np.where(d.accepted, 0.0, 1.0)}, {"big": d}
    )  # backwards score
    assert out["threshold"] is not None  # the lowest quantile skips nothing and is inside
    d.final_before[:] = True
    out = LR.choose_threshold({"big": np.where(d.accepted, 0.0, 1.0)}, {"big": d})
    assert out["threshold"] is not None
    # a score that is the same for all fits cannot skip anything either
    assert LR.skip_vector(np.full(5, 0.3), None).sum() == 0


def test_first_wrong_skip_and_s_rule_helpers() -> None:
    skip = np.array([False, True, True, False])
    acc = np.array([False, False, True, True])
    assert LR.first_wrong_skip(skip, acc) == 2
    assert LR.first_wrong_skip(skip, np.zeros(4, bool)) is None
    d = sample("big", 4)
    d.X[:, 2] = 0.8
    d.X[:, 3] = [0.1, 0.24, 0.3, 0.7]
    assert LR.s_rule_skip(d).tolist() == [
        True,
        False,
        False,
        False,
    ]  # 0.24 = 0.3 * 0.8: on the line


def test_plan_numbers_in_the_stored_table() -> None:
    assert [LR.EXPECT[s]["total"] for s in LR.SAMPLES] == [
        23_503_385,
        22_939_812,
        10_589_493,
        9_041_686,
        3_567_098,
    ]
    assert [LR.FITS[s] for s in LR.SAMPLES] == [72_458, 76_911, 52_027, 53_398, 32_543]
    for s in LR.SAMPLES:
        e = LR.EXPECT[s]
        assert e["revisits"] + e["nbr"] == LR.FITS[s]
        assert LR.PLAN_TABLE[s]["allowed"] == LR.FITS[s] // 1000


def test_gate_g1_needs_limits_and_net_share() -> None:
    a = dict(n=1000, wrong=1, changed=0, share=Fraction(5, 100))
    assert LR.gate_g1(dict(mid=a))["mid"]
    assert not LR.gate_g1(dict(grad=a))["grad"]  # grad needs 9%
    assert not LR.gate_g1(dict(mid=dict(a, wrong=2)))["mid"]
    assert not LR.gate_g1(dict(mid=dict(a, share=Fraction(44, 1000))))["mid"]


def test_judge_rules() -> None:
    base = dict(total=1000, median=0.010, wrong_seeds=3, grid=5, reseeds=100)
    j5 = dict(base, vs_s8=dict(broken=13, fixed=9))
    s8 = dict(total=5000, median=0.009)
    ok = dict(base, total=950, vs_j5=dict(broken=2, fixed=1), vs_s8=dict(broken=13, fixed=9))
    assert LR.judge(ok, j5, s8)["passes"]
    assert not LR.judge(dict(ok, total=980), j5, s8)["passes"]  # 2% fewer: rule 1
    assert not LR.judge(dict(ok, vs_j5=dict(broken=7, fixed=1)), j5, s8)["passes"]  # rule 2
    assert not LR.judge(dict(ok, median=0.0111), j5, s8)["passes"]  # rule 3
    assert not LR.judge(dict(ok, wrong_seeds=5), j5, s8)["passes"]  # rule 4
    assert not LR.judge(dict(ok, grid=7), j5, s8)["passes"]  # grid searches within 1, both ways
    assert not LR.judge(dict(ok, grid=3), j5, s8)["passes"]
    assert LR.judge(dict(ok, grid=4), j5, s8)["passes"]
    # J5 itself failing rule 5: judged by 1 to 4 only
    j5_bad = dict(j5, vs_s8=dict(broken=20, fixed=0))
    bad5 = dict(ok, vs_s8=dict(broken=21, fixed=0))
    assert LR.judge(bad5, j5_bad, s8)["passes"]
    assert not LR.judge(bad5, j5, s8)["passes"]


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
        self.hit_of = lambda R: 0.5
        self.n_calls = 0

    def evaluate_overlap(self, R: Any, vertices: Any, phase: int) -> _Info:
        self.n_calls += 1
        h = float(self.hit_of(np.asarray(R, np.float64)))
        return _Info(h, h / 2, 3.0 - h)


class _FakeBFS:
    def __init__(self) -> None:
        self.stats: Dict[str, int] = defaultdict(int)
        self.reconstructor = _Recon()
        self.fit_calls: List[Any] = []
        self.fit_hit = 0.8
        self.fit_R = rot_z(1.0)

    @staticmethod
    def _ratios(info: _Info) -> Any:
        return info.hit, info.conf

    def _local_fit(
        self, voxel: Any, vertices: Any, rng: Any, rec: Any, phase: str, sigma0: Any = None
    ) -> Any:
        self.fit_calls.append((rec.voxel_idx, phase, sigma0))
        rec.n_evals += 61
        self.stats[f"n_evals_{phase}"] += 61
        res = type("R", (), dict(orientation=self.fit_R, cost=1.5, overlap_info=None))()
        return res, self.fit_hit, self.fit_hit / 2

    def _expand(
        self, mic: Any, start_idx: int, best_conf: float, rng: Any, solution: Any, **kw: Any
    ) -> None:
        pass

    def _insert_seed(self, mic: Any, voxel_idx: int, queue: Any, *a: Any, **k: Any) -> None:
        queue.append(voxel_idx + 100)


class _Const:
    """A model that returns a fixed score (a test device with the PlainModel interface)."""

    def __init__(self, score: float, threshold: Optional[float]) -> None:
        self._s, self.threshold = score, threshold
        self.d: Dict[str, Any] = {}
        self.xs: List[List[float]] = []

    def score(self, x: Any) -> float:
        self.xs.append(list(x))
        return self._s

    def skips(self, s: float) -> bool:
        return self.threshold is not None and s < self.threshold


def _hook(score: float, threshold: Optional[float], tmp_path: Optional[Path] = None) -> Any:
    bfs = _FakeBFS()
    book = OD.DrawBook()
    model = _Const(score, threshold)
    h = LH.LearnHook(
        bfs,
        model,  # type: ignore[arg-type]
        SYM4,
        wrap=None,
        out_path=(tmp_path / "rows.jsonl") if tmp_path is not None else None,
        book=book,
        voxel_key=OD.voxel_key,
    )
    return bfs, h, book, model


def _call(
    bfs: Any, kind: str, idx: int, start: np.ndarray, best: float, n_rev: int = 0, conf: float = 0.4
) -> Any:
    bfs._expand(None, 0, best, None, [])
    rec = _Rec(idx, n_rev)
    vox = _Voxel(start, conf=conf)
    out = bfs._local_fit(vox, np.zeros((3, 3)), None, rec, kind)
    return out, rec


def test_every_fit_is_checked_once_and_a_skip_costs_nothing_more() -> None:
    bfs, h, book, model = _hook(0.1, 0.5)
    bfs.reconstructor.hit_of = lambda R: 0.1
    (res, hit, conf), rec = _call(bfs, "neighbor", 5, rot_z(3.0), 0.8)
    assert bfs.fit_calls == []  # skipped: the fit did not run
    assert rec.n_evals == 1 and bfs.stats["n_evals_neighbor"] == 1  # one counted check
    assert bfs.reconstructor.n_calls == 1
    assert h.counts["skipped"] == 1 and h.counts["checks"] == 1
    # a fit that is not skipped runs on top of its check
    bfs2, h2, _, _ = _hook(0.9, 0.5)
    (res, hit, conf), rec = _call(bfs2, "neighbor", 5, rot_z(3.0), 0.8)
    assert bfs2.fit_calls == [(5, "neighbor", None)] and hit == 0.8
    assert rec.n_evals == 62 and bfs2.stats["n_evals_neighbor"] == 62
    assert h2.counts["skipped"] == 0 and h2.counts["checks"] == 1


def test_skipped_neighbour_fit_holds_its_start_and_is_rejected() -> None:
    bfs, h, _, _ = _hook(0.1, 0.5)
    bfs.reconstructor.hit_of = lambda R: 0.1
    (res, hit, conf), _ = _call(bfs, "neighbor", 5, rot_z(3.0), 0.8)
    assert np.allclose(res.orientation, rot_z(3.0)) and res.cost == pytest.approx(2.9)
    assert hit == 0.1 and conf == 0.05
    assert not FR.accepted_test(hit, 0.8)  # a rejected fit
    assert h.counts["skipped_start_passes_test"] == 0


def test_skipped_revisit_keeps_the_old_end_through_a_negative_confidence() -> None:
    bfs, h, _, _ = _hook(0.1, 0.5)
    bfs.reconstructor.hit_of = lambda R: 0.1
    (res, hit, conf), rec = _call(bfs, "revisit", 5, rot_z(3.0), 0.8, n_rev=1, conf=0.4)
    assert conf == -1.0 < 0.4  # ``_expand``'s push rule (conf >= old) keeps the old fit
    assert bfs.fit_calls == [] and rec.n_evals == 1


def test_a_skipped_fit_has_no_fit_end_and_counts_as_a_rejected_fit() -> None:
    bfs, h, _, model = _hook(0.1, 0.5)
    bfs.reconstructor.hit_of = lambda R: 0.1
    _call(bfs, "neighbor", 5, rot_z(3.0), 0.8)
    _call(bfs, "revisit", 5, rot_z(9.0), 0.8, n_rev=1)
    x = dict(zip(NAMES, model.xs[1]))
    assert x["own_hit"] == -1.0  # no kept fit end
    assert x["fits_so_far"] == 1.0 and x["rejected_so_far"] == 1.0
    assert x["dist_rejected"] == pytest.approx(6.0, abs=1e-6)  # the skipped start was remembered


def test_the_hook_features_are_the_replays_features() -> None:
    """On a sequence with no skip, the live features equal ``rows_to_arrays`` fit by fit."""
    bfs, h, _, model = _hook(0.9, None)  # no threshold: nothing is skipped
    rows = []
    starts = [rot_z(3.0), rot_z(20.0), rot_z(23.0)]
    bfs.reconstructor.hit_of = lambda R: 0.4
    plan = [("neighbor", 5, 0), ("neighbor", 6, 0), ("revisit", 5, 1)]
    for k, ((phase, v, nrv), R) in enumerate(zip(plan, starts)):
        bfs.fit_hit = 0.5 + 0.1 * k
        bfs._expand(None, 0, 0.8, None, [])
        rec = _Rec(v, nrv)
        vox = _Voxel(R, conf=0.2)
        res, hit, conf = bfs._local_fit(vox, np.zeros((3, 3)), None, rec, phase)
        best = max(hit, 0.8)
        acc = bool(hit / best > 0.9)
        rows.append(
            dict(
                voxel=v,
                kind="revisit" if phase == "revisit" else "neighbour",
                n_revisits=nrv,
                best=0.8,
                start_R=R.ravel().tolist(),
                start_hit=0.4,
                start_cost=2.6,
                start_conf=0.2,
                parent=-1,
                pos=0,
                accepted=acc,
                kept=bool(phase != "revisit" or acc or conf >= 0.2),
                end_hit=hit,
                end_cost=1.5,
                evals=61,
            )
        )
    out = LFE.rows_to_arrays(rows, SYM4)
    assert np.allclose(np.array(model.xs), out["X"], atol=1e-12)


def test_a_model_call_counts_zero_evaluations_and_timers_run() -> None:
    bfs, h, _, _ = _hook(0.9, 0.5)
    _call(bfs, "neighbor", 5, rot_z(3.0), 0.8)
    assert bfs.stats["n_evals_neighbor"] == 62  # 61 + one check; the model call adds none
    s = h.summary()
    assert s["model_wall_s"] >= 0.0 and s["feature_wall_s"] >= 0.0 and s["counts"]["scored"] == 1


def test_restart_retries_pass_through_unchanged() -> None:
    bfs, h, _, model = _hook(0.1, 0.5)
    bfs._expand(None, 0, 0.8, None, [])
    rec = _Rec(3, 0)
    bfs._local_fit(_Voxel(rot_z(1.0)), np.zeros((3, 3)), None, rec, "neighbor", sigma0_deg=5.0)
    assert bfs.fit_calls == [(3, "neighbor", 5.0)] and model.xs == [] and h.counts["checks"] == 0


def test_a_skipped_fit_takes_and_discards_its_aligned_draw() -> None:
    key = OD.voxel_key(np.zeros((3, 3)))
    plain = OD.DrawBook()
    plain.rng(key, "_local_fit")
    second = plain.rng(key, "_local_fit")
    bfs, h, book, _ = _hook(0.1, 0.5)
    bfs.reconstructor.hit_of = lambda R: 0.1
    _call(bfs, "neighbor", 5, rot_z(3.0), 0.8)  # skipped: takes draw 0 and drops it
    assert book.counts[(key, "_local_fit")] == 1
    assert book.rng(key, "_local_fit").random() == second.random()


def test_rows_are_written_with_features_score_and_skip(tmp_path: Path) -> None:
    bfs, h, _, _ = _hook(0.1, 0.5, tmp_path)
    bfs.reconstructor.hit_of = lambda R: 0.1
    _call(bfs, "neighbor", 5, rot_z(3.0), 0.8)
    bfs.reconstructor.hit_of = lambda R: 0.6
    h.model.threshold = 0.05  # the next fit is not skipped
    _call(bfs, "revisit", 5, rot_z(9.0), 0.8, n_rev=1)
    h.close()
    rows = [json.loads(line) for line in (tmp_path / "rows.jsonl").read_text().splitlines()]
    assert [r["skipped"] for r in rows] == [True, False] and len(rows[0]["x"]) == 16
    assert rows[0]["evals"] == 0 and rows[1]["evals"] == 61 and rows[0]["end_hit"] is None
    assert rows[1]["accepted"] and rows[1]["end_hit"] == 0.8


def test_a_skip_whose_start_passes_the_acceptance_test_is_visible() -> None:
    bfs, h, _, _ = _hook(0.1, 0.5)
    bfs.reconstructor.hit_of = lambda R: 0.79  # 0.79 / 0.8 over 0.9: the library would accept it
    _call(bfs, "neighbor", 5, rot_z(3.0), 0.8)
    assert h.counts["skipped"] == 1 and h.counts["skipped_start_passes_test"] == 1


# ---- guards
def test_held2_is_refused_everywhere(tmp_path: Path) -> None:
    for name in ("held2", "xheld2"):
        with pytest.raises(SystemExit) as e:
            LG.refuse_sample(name)
        assert e.value.code == 2
    for p in (tmp_path / "J_held2" / "run.json", tmp_path / "passes_held2", "x/held2_clean.npy"):
        with pytest.raises(SystemExit):
            LG.guarded(p, cache=tmp_path)
    with pytest.raises(SystemExit):
        LG.refuse_sample("nonsense")


def test_held4_is_refused_before_the_gate_and_its_files_before_the_marker(tmp_path: Path) -> None:
    with pytest.raises(SystemExit) as e:
        LG.refuse_sample("held4", gate=tmp_path / "gate.json")
    assert e.value.code == 2
    (tmp_path / "gate.json").write_text("{}")
    LG.refuse_sample("held4", gate=tmp_path / "gate.json")  # now allowed (the build, the runs)
    for s in LG.DISCLOSED:
        LG.refuse_sample(s)
    p = tmp_path / "J5L_held4" / "run.json"
    with pytest.raises(SystemExit):
        LG.guarded(p, cache=tmp_path)
    assert LG.guarded(tmp_path / "J5L_mid" / "run.json", cache=tmp_path)
    (tmp_path / LG.H4_MARKER).write_text("{}")
    assert LG.guarded(p, cache=tmp_path) == p


def test_the_guard_module_names_step_18s_gate_and_marker() -> None:
    assert LG.gate_file().parts[-2:] == ("learnskip", "gate.json")
    assert LG.H4_MARKER == "held4_runs_done.json"
    assert LG.cache_dir().parts[-2:] == ("learnskip", "cache")


def test_samples_support_has_held3_and_held4_seeds() -> None:
    import learnskip_samples as LS

    assert LS.SEEDS == {"held3": 3, "held4": 4}
    src = LS.patched_source("make", "held4")
    assert "SEED = 4" in src and "SEED = {" not in src
    assert 'K = "held4"' in LS.patched_source("stack", "held4")
