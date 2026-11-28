"""AR5 step 29 (rule RU): before a seed step, try the orientations of the 8 nearest accepted voxels
with a longer CMA fit (first step 0.5 degree, 200 evaluations; a second fit from the lowest-cost
start over 0.5 degree from the first), and index only if both fail.

Checks: no try when nothing is accepted; the 8 nearest by distance, ties by index; rejected voxels
left out; the second start more than 0.5 degree from the first; a success skips the seed step and
returns the fit's end; a failure runs the PXF seed step unchanged (equal to a direct step); the
evaluations are counted; a try never accepts an end under 0.95; the hexagonal fixture runs; the
judges and the frozen constants.
"""

import sys
from pathlib import Path
from types import SimpleNamespace
from typing import Any, Callable, Dict, List, Tuple

import numpy as np
import pytest
from scipy.spatial import cKDTree
from scipy.spatial.transform import Rotation

from icenine import ray_indexing as RI

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts" / "autoresearch" / "ar5"))
import ru_rules as RR  # noqa: E402


def grid(n: int = 7) -> np.ndarray:
    """A square grid of voxel centres (unit spacing), index = row * n + column."""
    return np.array([(float(i % n), float(i // n)) for i in range(n * n)])


def near(pos: np.ndarray, v: int, ok: Callable[[int], bool], k: int = 8) -> List[int]:
    return RR.nearest_accepted(cKDTree(pos), pos, v, ok, k)


def test_no_try_when_nothing_is_accepted() -> None:
    pos = grid()
    assert near(pos, 24, lambda j: False) == []
    calls: List[str] = []
    end, row = RR.run_try(
        [],
        lambda j: np.eye(3),
        lambda R: calls.append("cost") or 0.0,  # type: ignore[func-returns-value]
        lambda R: calls.append("fit") or (R, 1, 0.0),  # type: ignore[func-returns-value]
        lambda R: 1.0,
        lambda a, b: 0.0,
    )
    assert end is None and not row["tried"] and row["evals"] == 0 and calls == []


def test_the_8_nearest_by_distance_ties_by_index() -> None:
    pos = grid()
    v = 24  # the centre of the 7 x 7 grid
    got = near(pos, v, lambda j: True)
    # distance 1: 17, 23, 25, 31 (ties by index); distance sqrt2: 16, 18, 30, 32
    assert got == [17, 23, 25, 31, 16, 18, 30, 32]
    # the query itself is never a start; fewer if fewer exist
    assert v not in got
    assert near(pos, v, lambda j: j in (0, 48)) == [0, 48]
    # all accepted voxels at one distance beyond the first query size: still ties by index
    big = grid(31)  # 961 voxels, the centre is 480
    got = near(big, 480, lambda j: True)
    assert got[:4] == [449, 479, 481, 511] and got[4:] == [448, 450, 510, 512]


def test_rejected_voxels_are_left_out_and_any_distance_counts() -> None:
    pos = grid()
    v = 24
    rejected = {17, 23, 25, 31}
    got = near(pos, v, lambda j: j not in rejected)
    assert not (set(got) & rejected) and got[:4] == [16, 18, 30, 32] and len(got) == 8
    far = near(pos, v, lambda j: j == 0)  # only a distant voxel is accepted: it is the start
    assert far == [0]


def test_second_start_is_more_than_half_a_degree_from_the_first() -> None:
    def rz(deg: float) -> np.ndarray:
        return np.asarray(Rotation.from_euler("z", deg, degrees=True).as_matrix())

    sym = np.eye(3)[None]
    mis = lambda a, b: float(RI.misorientation_deg(a[None], b[None], sym)[0])  # noqa: E731
    R = [rz(0.0), rz(0.3), rz(0.5), rz(0.7), rz(5.0)]
    # cost order 0, 1, 2, 3, 4: 0.3 and 0.5 are not over 0.5 degree from the first; 0.7 is
    assert RR.second_start([0, 1, 2, 3, 4], R, mis) == 3
    assert RR.second_start([0, 1, 2], R, mis) is None
    assert RR.second_start([3, 2, 0], R, mis) == 0  # measured from the lowest-cost one (3)
    assert RR.second_start([4, 0], R, mis) == 0


def fixture(
    hits: Dict[int, float],
    angles: Dict[int, float],
    costs: Dict[int, float],
) -> Tuple[Any, ...]:
    """Starts are z-rotations by ``angles[j]``; a fit from start ``j`` ends at the start itself
    with hit ``hits[j]``; the strict cost of start ``j`` is ``costs[j]``."""
    rz = lambda d: np.asarray(Rotation.from_euler("z", d, degrees=True).as_matrix())  # noqa: E731
    sym = np.eye(3)[None]
    R_of = {j: rz(a) for j, a in angles.items()}

    def key(R: np.ndarray) -> int:
        return min(R_of, key=lambda j: float(np.abs(R_of[j] - R).max()))

    fits: List[int] = []

    def fit(R: np.ndarray) -> Tuple[np.ndarray, int, float]:
        fits.append(key(R))
        return R, 201, 0.25

    return (
        lambda j: R_of[j],
        lambda R: costs[key(R)],
        fit,
        lambda R: hits[key(R)],
        lambda a, b: float(RI.misorientation_deg(a[None], b[None], sym)[0]),
        fits,
    )


def test_a_success_returns_the_fit_end_and_counts_evaluations() -> None:
    orient, cost, fit, hit, mis, fits = fixture(
        hits={1: 0.99, 2: 0.99}, angles={1: 0.0, 2: 3.0}, costs={1: 0.4, 2: 0.2}
    )
    end, row = RR.run_try([1, 2], orient, cost, fit, hit, mis)
    assert fits == [2]  # the lowest-cost start (voxel 2) is fitted first and succeeds
    assert row["success"] and not row["second"] and end is not None
    assert row["evals"] == 2 + 201  # two strict costs and one fit
    assert row["fits"][0]["hit"] >= RR.ACCEPT and row["fits"][0]["start"] == 2


def test_the_second_start_runs_only_after_a_first_fit_under_095() -> None:
    orient, cost, fit, hit, mis, fits = fixture(
        hits={1: 0.94999, 2: 0.96, 3: 0.99},
        angles={1: 0.0, 2: 0.2, 3: 4.0},
        costs={1: 0.1, 2: 0.2, 3: 0.3},
    )
    end, row = RR.run_try([1, 2, 3], orient, cost, fit, hit, mis)
    # first fit from voxel 1 ends at 0.94999 (under 0.95); voxel 2 is 0.2 degree away (not a
    # second start), voxel 3 is 4 degrees away: the second fit
    assert fits == [1, 3] and row["success"] and row["second"]
    assert row["evals"] == 3 + 2 * 201
    assert [f["second"] for f in row["fits"]] == [False, True]
    # both fail: no end, both fits counted, the step falls back
    orient, cost, fit, hit, mis, fits = fixture(
        hits={1: 0.5, 2: 0.5}, angles={1: 0.0, 2: 3.0}, costs={1: 0.1, 2: 0.2}
    )
    end, row = RR.run_try([1, 2], orient, cost, fit, hit, mis)
    assert end is None and not row["success"] and row["evals"] == 2 + 2 * 201
    # no second start (all starts the same orientation): one fit
    orient, cost, fit, hit, mis, fits = fixture(
        hits={1: 0.5, 2: 0.5}, angles={1: 0.0, 2: 0.1}, costs={1: 0.1, 2: 0.2}
    )
    end, row = RR.run_try([1, 2], orient, cost, fit, hit, mis)
    assert end is None and fits == [1] and row["evals"] == 2 + 201


def test_exactly_095_is_accepted_and_nothing_under_it() -> None:
    orient, cost, fit, hit, mis, fits = fixture(hits={1: 0.95}, angles={1: 0.0}, costs={1: 0.1})
    end, row = RR.run_try([1], orient, cost, fit, hit, mis)
    assert end is not None and row["success"]


class FakeResult:
    def __init__(self, evals: int, total_s: float, orientation: Any = None) -> None:
        self.evals, self.total_s, self.orientation = evals, total_s, orientation
        self.n_cand = 3
        self.taken = 1


def run_wrapper(near_ids: List[int], hits: Dict[int, float]) -> Tuple[Any, List[Any], Any, Any]:
    orient, cost, fit, hit, mis, fits = fixture(
        hits=hits,
        angles={1: 0.0, 2: 3.0},
        costs={1: 0.1, 2: 0.2},
    )
    direct = FakeResult(500, 1.5, np.eye(3))
    calls: List[Any] = []

    def orig(self: Any, vertices: Any, phase: int, *a: Any, **k: Any) -> Any:
        calls.append((vertices, phase, a, k))
        return direct

    log: List[Dict[str, Any]] = []

    def ctx(self: Any, vertices: Any, phase: int) -> Dict[str, Any]:
        return dict(
            voxel=7,
            reseed=False,
            near=near_ids,
            orient_of=orient,
            strict_cost=cost,
            fit=fit,
            hit_of=hit,
            mis_deg=mis,
        )

    step = RR.make_index_seed(orig, log, ctx, RI.SeedResult)
    return step, calls, log, direct


def test_a_success_skips_the_seed_step() -> None:
    step, calls, log, _ = run_wrapper([1, 2], {1: 0.99, 2: 0.99})
    out = step(SimpleNamespace(), "V", 0, 20, 5)
    assert calls == []  # no indexing, no seed step
    assert isinstance(out, RI.SeedResult) and out.orientation is not None
    assert out.evals == 2 + 201 and out.n_cand == 0 and out.index_s == 0.0
    assert out.cost == pytest.approx(0.25)  # the fit's own end cost; no extra evaluation
    assert log[-1]["success"] and log[-1]["indexing_calls"] == 0 and log[-1]["voxel"] == 7


def test_a_failure_runs_the_pxf_step_unchanged_and_adds_the_try() -> None:
    step, calls, log, direct = run_wrapper([1, 2], {1: 0.5, 2: 0.5})
    ev0, t0 = direct.evals, direct.total_s
    out = step(SimpleNamespace(), "V", 0, 20, 5)
    assert out is direct and calls == [("V", 0, (20, 5), {})]  # the same call, the same result
    assert out.evals == ev0 + 2 + 2 * 201 and out.total_s >= t0  # the try counts as the seed's
    assert out.n_cand == 3 and out.taken == 1  # nothing else of the step changed
    assert not log[-1]["success"] and log[-1]["indexing_calls"] == 1 and log[-1]["tried"]


def test_no_start_runs_the_pxf_step_with_no_added_evaluations() -> None:
    step, calls, log, direct = run_wrapper([], {})
    out = step(SimpleNamespace(), "V", 0)
    assert out is direct and out.evals == 500 and len(calls) == 1
    assert not log[-1]["tried"] and log[-1]["indexing_calls"] == 1


def test_indexing_calls_equal_steps_minus_successes() -> None:
    rows: List[Dict[str, Any]] = []
    for ids, hits in (([1, 2], {1: 0.99, 2: 0.99}), ([1, 2], {1: 0.1, 2: 0.1}), ([], {})):
        step, calls, log, _ = run_wrapper(ids, hits)
        step(SimpleNamespace(), "V", 0)
        rows += log
    assert len(rows) == 3
    assert sum(r["indexing_calls"] for r in rows) == len(rows) - sum(r["success"] for r in rows)


def test_hexagonal_fixture_runs() -> None:
    from pymatgen.core.lattice import Lattice
    from pymatgen.core.structure import Structure

    from icenine.symmetry import CrystalSymmetry

    st = Structure(
        Lattice.hexagonal(2.95, 4.68), ["Ti", "Ti"], [[1 / 3, 2 / 3, 0.25], [2 / 3, 1 / 3, 0.75]]
    )
    sym = RI.proper_rotations(CrystalSymmetry(st))
    assert len(sym) == 12
    mis = lambda a, b: float(RI.misorientation_deg(a[None], b[None], sym)[0])  # noqa: E731
    R0 = np.asarray(Rotation.random(random_state=3).as_matrix())
    # a start 60 degrees about the c axis is a symmetry copy (0 degree); 0.4 and 1.0 degree are not
    cz = lambda d: np.asarray(Rotation.from_euler("z", d, degrees=True).as_matrix())  # noqa: E731
    Rs = [R0, R0 @ cz(60.0), R0 @ cz(0.4), R0 @ cz(1.0)]
    assert mis(Rs[0], Rs[1]) < 1e-6
    assert RR.second_start([0, 1, 2, 3], Rs, mis) == 3  # 1 is a copy, 2 is within 0.5 degree
    hex_pos = np.array(
        [(i + 0.5 * (j % 2), j * np.sqrt(3) / 2) for j in range(6) for i in range(6)]
    )
    got = near(hex_pos, 14, lambda j: True)  # a hex lattice: 6 equidistant neighbours first
    assert len(got) == 8 and len(set(got[:6])) == 6
    d = np.linalg.norm(hex_pos[got] - hex_pos[14], axis=1)
    assert np.allclose(d[:6], 1.0) and d[6] > 1.0


def test_judges_and_constants() -> None:
    assert RR.judge_r1(19_892_329)["ok"] and not RR.judge_r1(19_892_330)["ok"]  # 15% of 23,402,740
    assert RR.OLD_TOTAL == 6_822_192 + 16_580_548 == 23_402_740
    assert RR.judge_r2(3, 0.001, 2)["ok"] and not RR.judge_r2(4, 0.0, 0)["ok"]
    assert not RR.judge_r2(0, 0.0011, 0)["ok"] and not RR.judge_r2(0, 0.0, 3)["ok"]
    assert RR.judge_r3(3_863)["ok"] and not RR.judge_r3(3_864)["ok"]  # 15% lower than 4,545 s
    assert RR.judge_h1(0.2, 5, 0.005)["ok"] and not RR.judge_h1(0.19, 0, 0.0)["ok"]
    assert RR.judge_h2(0.10, 5, 0.001)["ok"] and not RR.judge_h2(0.099, 0, 0.0)["ok"]
    assert not RR.judge_h2(0.2, 6, 0.0)["ok"] and not RR.judge_h2(0.2, 0, 0.0011)["ok"]
    assert (RR.K_NEAREST, RR.SIGMA0_DEG, RR.MAX_EVALS, RR.SEP_DEG, RR.ACCEPT) == (
        8,
        0.5,
        200,
        0.5,
        0.95,
    )
