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
    ) -> None:
        self.truth = {v: (0.0 if v % W < 3 else 30.0) for v in range(W * W)}
        self.truth.update(truth_deg or {})
        self.weak = weak or {}
        self.both = both or {}
        self.no_seed = no_seed or set()
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
        if v in self.both and abs(a - 0.0) < 0.5:
            return self.both[v][0], 0.95
        if v in self.both and abs(a - 30.0) < 0.5:
            return self.both[v][1], 0.95
        if abs(a - self.truth[v]) < 0.5:
            return 0.05, self.weak.get(v, 0.95)
        return 0.8, 0.2

    def local_fit(self, v: int, R0: np.ndarray) -> Tuple[np.ndarray, float, float, int]:
        self.fit_calls.append((v, self._ang(R0)))
        return rz(self.truth[v]), 0.04, 0.95, 100

    def index_seed(self, v: int) -> Tuple[Optional[np.ndarray], float, float, int]:
        if v in self.no_seed:
            return rz(self.truth[v]), 0.5, 0.5, 10
        return rz(self.truth[v]), 0.05, 0.95, 10

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


def test_flood_best_hit_and_stop_at_failed_neighbour() -> None:
    # voxel 8 (row 1, col 2) is a weak voxel of grain A: hit 0.5 < 0.9 * 0.95: A does not pass there
    h = FakeHooks(weak={8: 0.5})
    res = paint_grains([0] + list(range(W * W)), W * W, h)
    assert 8 in res.scores[0] and res.scores[0][8][2] is False
    assert res.scores[0][8][1] == 0.5


def test_settled_by_a_fit_when_scored_but_not_passed() -> None:
    # voxel 8 is weak for A (scored, not passed): it gets one fit from A, which passes
    h = FakeHooks(weak={8: 0.5})
    res = paint_grains([0, 5, 8] + list(range(W * W)), W * W, h)
    assert res.stage[8] == "fit"
    assert res.fits[8].from_grain == 0
    assert (8, 0.0) in h.fit_calls
    assert res.evals["fit_accepted"] == 100 * res.counts["fits_accepted"]
    assert painted(res)[8] == 0.0


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
    # a column of failed voxels (col 2) hides A's interior from the flood of the seed in col 0;
    # voxel 2 fails its seed paths and the grid search
    h = FakeHooks(no_seed={2})
    res = paint_grains([2] + list(range(W * W)), W * W, h)
    assert res.counts["failed"] >= 1 and res.stage[2] == "failed"
    # the failed voxel 2 is passed by the flood of the grain seeded at voxel 0 (its truth is A)
    assert res.painted_grain(2) is not None
    R, src = res.final_orientation(2)
    assert src == "grain"
    assert res.evals["seed_rejected"] == 10 + 1000


def test_near_earlier_entry_is_counted() -> None:
    # every voxel of column 1 is weak for grain A: the flood from column 0 scores them and stops, so
    # column 2 is reached by nobody and seeds a second entry with the same orientation
    h = FakeHooks(weak={r * W + 1: 0.1 for r in range(W)})
    res = paint_grains([0] + list(range(W * W)), W * W, h)
    a = [e for e in res.entries if e.R[0, 0] > 0.99]  # orientation of grain A
    assert len(a) >= 2
    assert a[0].near_earlier is None
    assert a[1].near_earlier is not None and a[1].near_deg < 1.0


def test_rejected_seed_voxel_with_unrelated_orientation_counts_both_paths() -> None:
    h = FakeHooks(no_seed={35})
    res = paint_grains([35, 0], W * W, h)
    assert res.stage[35] == "failed"
    assert res.log[0]["seed"] == "failed"
    assert res.evals["seed_rejected"] == 1010
