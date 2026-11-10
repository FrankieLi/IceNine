"""Unit tests of icenine.grain_paint on a small grid with two grains and a known boundary (a fake
cost; no physics)."""

from typing import Dict, List, Optional, Set, Tuple

import numpy as np

from icenine.grain_paint import PaintResult, paint_grains

W = 6  # grid width and height; voxel v = row * W + col


def rz(deg: float) -> np.ndarray:
    a = np.deg2rad(deg)
    return np.array([[np.cos(a), -np.sin(a), 0], [np.sin(a), np.cos(a), 0], [0, 0, 1.0]])


class FakeHooks:
    """Grain A (0 degrees) holds columns 0-2, grain B (30 degrees) columns 3-5. An orientation
    within 0.5 degrees of the voxel's truth scores cost 0.05 and hit 0.95; otherwise 0.8 and 0.2."""

    min_seed_hit = 0.8

    def __init__(
        self,
        truth_deg: Optional[Dict[int, float]] = None,
        weak: Optional[Dict[int, float]] = None,  # voxel -> hit its own grain gets (a "bad" voxel)
        both: Optional[Dict[int, Tuple[float, float]]] = None,  # voxel -> costs of A, B (both pass)
        no_seed: Optional[Set[int]] = None,  # voxels whose seed paths fail (hit 0.5)
        index_fail: Optional[Set[int]] = None,  # voxels whose index path fails, grid passes
        fit_hit: Optional[Dict[int, float]] = None,  # voxel -> hit of its local fit (else 0.95)
        seed_hit: Optional[Dict[int, float]] = None,  # voxel -> hit of its index path (else 0.95)
        override: Optional[Dict[Tuple[int, int], Tuple[float, float]]] = None,  # (v, deg) -> score
    ) -> None:
        self.truth = {v: (0.0 if v % W < 3 else 30.0) for v in range(W * W)}
        self.truth.update(truth_deg or {})
        self.weak = weak or {}
        self.both = both or {}
        self.no_seed = no_seed or set()
        self.index_fail = index_fail or set()
        self.fit_hit = fit_hit or {}
        self.seed_hit = seed_hit or {}
        self.override = override or {}
        self.n_scores = 0
        self.fit_calls: List[Tuple[int, float]] = []

    def neighbours(self, v: int) -> List[int]:
        r, c = divmod(v, W)
        out = []
        for dr, dc in ((-1, 0), (1, 0), (0, -1), (0, 1)):
            if 0 <= r + dr < W and 0 <= c + dc < W:
                out.append((r + dr) * W + c + dc)
        return out

    def _ang(self, R: np.ndarray) -> float:
        return float(np.rad2deg(np.arctan2(R[1, 0], R[0, 0])))

    def score(self, R: np.ndarray, v: int) -> Tuple[float, float]:
        self.n_scores += 1
        a = self._ang(R)
        if (v, int(round(a))) in self.override:
            return self.override[(v, int(round(a)))]
        if v in self.both and abs(a - 0.0) < 0.5:
            return self.both[v][0], 0.95
        if v in self.both and abs(a - 30.0) < 0.5:
            return self.both[v][1], 0.95
        if abs(a - self.truth[v]) < 0.5:
            return 0.05, self.weak.get(v, 0.95)
        return 0.8, 0.2

    def local_fit(self, v: int, R0: np.ndarray) -> Tuple[np.ndarray, float, float, int]:
        self.fit_calls.append((v, self._ang(R0)))
        return rz(self.truth[v]), 0.04, self.fit_hit.get(v, 0.95), 100

    def index_seed(self, v: int) -> Tuple[Optional[np.ndarray], float, float, int]:
        if v in self.no_seed or v in self.index_fail:
            return rz(self.truth[v]), 0.5, 0.5, 10
        return rz(self.truth[v]), 0.05, self.seed_hit.get(v, 0.95), 10

    def grid_seed(self, v: int) -> Tuple[np.ndarray, float, float, int]:
        if v in self.no_seed:
            return rz(self.truth[v]), 0.5, 0.5, 1000
        return rz(self.truth[v]), 0.05, 0.95, 1000

    def misorientation_deg(self, Ra: np.ndarray, Rb: np.ndarray) -> float:
        return abs(self._ang(Ra) - self._ang(Rb))


def painted(res: PaintResult) -> np.ndarray:
    out = np.full(W * W, np.nan)
    for v in range(W * W):
        R, _ = res.final_orientation(v)
        if R is not None:
            out[v] = float(np.rad2deg(np.arctan2(R[1, 0], R[0, 0])))
    return out


def truth(h: FakeHooks) -> np.ndarray:
    return np.array([h.truth[v] for v in range(W * W)])


def test_two_grains_painted_without_fits() -> None:
    h = FakeHooks()
    order = [0, 5] + list(range(W * W))  # a seed in each grain first
    res = paint_grains(order, W * W, h)
    assert len(res.entries) == 2
    assert np.allclose(painted(res), truth(h))
    assert res.counts["fits_accepted"] == 0 and res.counts["failed"] == 0
    assert res.evals["seed_index"] == 20
    # each grain scored on its own voxels (17 besides the seed) and on the other grain's boundary
    # column once (6 voxels): 2 * 17 + 2 * 6
    assert res.evals["score"] == h.n_scores == 2 * 17 + 2 * 6
    # no voxel is passed by both grains
    assert all(len(res.passing(v)) == 1 for v in range(W * W))
    assert res.evals["fit_accepted"] == 0 and res.evals["seed_grid"] == 0


def test_weak_voxel_is_scored_but_does_not_pass() -> None:
    # voxel 8 (row 1, col 2) is a weak voxel of grain A: hit 0.5 < 0.9 * 0.95: A does not pass there
    h = FakeHooks(weak={8: 0.5})
    res = paint_grains([0] + list(range(W * W)), W * W, h)
    assert 8 in res.scores[0] and res.scores[0][8][2] is False
    assert res.scores[0][8][1] == 0.5


def test_flood_best_is_raised_before_the_test() -> None:
    # the seed (voxel 0) has hit 0.8. Its first neighbour in queue order, voxel 6, scores 0.95 and
    # raises the flood's best to 0.95. Voxel 1, scored next with hit 0.8, then fails
    # (0.8 / 0.95 < 0.9), although it would pass against the seed's own hit (0.8 / 0.8).
    h = FakeHooks(seed_hit={0: 0.8}, weak={1: 0.8})
    res = paint_grains([0], W * W, h)
    assert res.entries[0].hit == 0.8
    assert res.scores[0][6] == (0.05, 0.95, True)
    assert res.scores[0][1] == (0.05, 0.8, False)
    assert res.entries[0].best_hit == 0.95


def test_settled_by_a_fit_when_scored_but_not_passed() -> None:
    # voxel 8 is weak for A (scored, not passed): it gets one fit from A, which passes
    h = FakeHooks(weak={8: 0.5})
    res = paint_grains([0, 5, 8] + list(range(W * W)), W * W, h)
    assert res.stage[8] == "fit"
    assert res.fits[8].from_grain == 0
    assert (8, 0.0) in h.fit_calls
    assert res.evals["fit_accepted"] == 100 * res.counts["fits_accepted"]
    assert painted(res)[8] == 0.0


def test_rejected_fit_then_the_voxel_becomes_a_seed() -> None:
    # voxel 8 is weak for A (scored, not passed); its fit from A also gets hit 0.5, so the fit is
    # rejected (0.5 / 0.95 < 0.9) and the voxel goes on to the seed path, which accepts it
    h = FakeHooks(weak={8: 0.5}, fit_hit={8: 0.5})
    res = paint_grains([0, 5, 8] + list(range(W * W)), W * W, h)
    assert (8, 0.0) in h.fit_calls
    assert res.counts["fits_rejected"] == 1 and res.counts["fits_accepted"] == 0
    assert res.evals["fit_rejected"] == 100 and res.evals["fit_accepted"] == 0
    assert res.stage[8] == "seed" and 8 not in res.fits
    row = next(r for r in res.log if r["voxel"] == 8)
    assert row["fit_from"] == 0 and row["fit_ok"] is False and row["seed"] == "index"
    assert res.entries[2].voxel == 8 and res.entries[2].near_earlier == 0
    assert res.evals["seed_index"] == 3 * 10


def test_fit_settled_voxel_keeps_its_fit_over_a_later_passing_grain() -> None:
    # voxel 8 (A by truth) is weak for A, so it is settled by a fit from A (cost 0.04). Grain B,
    # seeded later at voxel 5, then passes on voxel 8 with a lower cost (0.01). The fit is kept.
    h = FakeHooks(override={(8, 0): (0.05, 0.5), (8, 30): (0.01, 0.95)})
    res = paint_grains([0, 8, 5] + list(range(W * W)), W * W, h)
    assert res.stage[8] == "fit" and res.fits[8].cost == 0.04
    assert res.scores[1][8] == (0.01, 0.95, True)
    assert res.painted_grain(8) == 1
    R, src = res.final_orientation(8)
    assert src == "fit" and R is not None and np.allclose(R, rz(0.0))


def test_grid_fallback_accepted_when_the_index_path_fails() -> None:
    h = FakeHooks(index_fail={2})
    res = paint_grains([2] + list(range(W * W)), W * W, h)
    assert res.entries[0].voxel == 2 and res.entries[0].kind == "grid"
    assert res.stage[2] == "seed" and res.log[0]["seed"] == "grid"
    assert res.evals["index_rejected"] == 10
    assert res.evals["seed_grid"] == 1000
    assert res.evals["grid_rejected"] == 0
    assert np.allclose(painted(res), truth(h))


def test_flood_does_not_depend_on_the_seed_order() -> None:
    def hooks() -> FakeHooks:  # a contested voxel and a weak voxel, so the floods are not trivial
        return FakeHooks(both={3: (0.30, 0.10)}, weak={8: 0.5})

    a = paint_grains([0, 5] + list(range(W * W)), W * W, hooks())
    b = paint_grains([5, 0] + list(range(W * W)), W * W, hooks())
    assert [e.voxel for e in a.entries][:2] == [0, 5]
    assert [e.voxel for e in b.entries][:2] == [5, 0]
    assert a.scores[0] == b.scores[1]  # grain A
    assert a.scores[1] == b.scores[0]  # grain B
    assert a.entries[0].best_hit == b.entries[1].best_hit


def test_contested_voxel_takes_the_lowest_cost() -> None:
    # voxel 3 (row 0, col 3) is B's by truth, but A passes there too with a higher cost
    h = FakeHooks(both={3: (0.30, 0.10)})
    res = paint_grains([0, 5] + list(range(W * W)), W * W, h)
    p = res.ranked(3)
    assert [g for g, _, _ in p] == [1, 0]  # B (entry 1) first
    assert res.painted_grain(3) == 1
    # swap the costs: A wins
    h2 = FakeHooks(both={3: (0.05, 0.30)})
    res2 = paint_grains([0, 5] + list(range(W * W)), W * W, h2)
    assert res2.painted_grain(3) == 0
    assert painted(res2)[3] == 0.0


def test_failed_voxel_is_final_and_unpainted_unless_a_grain_passes() -> None:
    # voxel 2 fails its seed path and the grid search; it is then passed by the flood of the grain
    # seeded next (at voxel 0)
    h = FakeHooks(no_seed={2})
    res = paint_grains([2] + list(range(W * W)), W * W, h)
    assert res.counts["failed"] >= 1 and res.stage[2] == "failed"
    # the failed voxel 2 is passed by the flood of the grain seeded at voxel 0 (its truth is A)
    assert res.painted_grain(2) is not None
    R, src = res.final_orientation(2)
    assert src == "grain"
    assert res.evals["index_rejected"] == 10
    assert res.evals["grid_rejected"] == 1000
    assert res.evals["seed_grid"] == 0


def test_near_earlier_entry_is_counted() -> None:
    # every voxel of column 1 is weak for grain A: the flood from column 0 scores them and stops, so
    # column 2 is reached by nobody and seeds a second entry with the same orientation
    h = FakeHooks(weak={r * W + 1: 0.1 for r in range(W)})
    res = paint_grains([0] + list(range(W * W)), W * W, h)
    a = [e for e in res.entries if e.R[0, 0] > 0.99]  # orientation of grain A
    assert len(a) >= 2
    assert a[0].near_earlier is None
    assert a[1].near_earlier is not None and a[1].near_deg < 1.0
    # the first entry's flood stops at the weak column: it never scores column 2
    assert not any(r * W + 2 in res.scores[0] for r in range(W))


def test_rejected_seed_voxel_with_unrelated_orientation_counts_both_paths() -> None:
    h = FakeHooks(no_seed={35})
    res = paint_grains([35, 0], W * W, h)
    assert res.stage[35] == "failed"
    assert res.log[0]["seed"] == "failed"
    assert res.evals["index_rejected"] == 10
    assert res.evals["grid_rejected"] == 1000
    # no grain passes on voxel 35 (grain A's flood stays in columns 0-2): it keeps the rejected
    # search's result
    assert res.painted_grain(35) is None
    R, src = res.final_orientation(35)
    assert src == "failed" and R is not None
