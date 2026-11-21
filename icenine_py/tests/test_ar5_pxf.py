"""AR5 step 28 (rule PXF): a seed step that fails the seed test after a subsampled indexing
(50,000 pairs) is run again with the full sample (600,000 pairs).

Checks: no fallback at or under 50,000 pairs; no fallback when the first step passes; the fallback
step equals a direct 600,000-pair step at the voxel; no candidate also triggers it; both steps'
evaluations are counted; the rescue (``G1``) sees the second step's kept ends; the pick's rows are
tagged; the hexagonal fixture runs; the judges.
"""

import sys
from pathlib import Path
from types import SimpleNamespace
from typing import Any, Dict, List, Tuple

import numpy as np
import pytest
from scipy.spatial.transform import Rotation

from icenine import ray_indexing as RI
from icenine.symmetry import create_cubic_symmetry

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts" / "autoresearch" / "ar5"))
import pxf_rules as PR  # noqa: E402

A_LAT = 3.61
Q0 = 2.0 * np.pi / A_LAT


def fcc_vectors(max_n2: int = 20) -> np.ndarray:
    out = []
    r = range(-4, 5)
    for h in r:
        for k in r:
            for ll in r:
                n2 = h * h + k * k + ll * ll
                if n2 == 0 or n2 > max_n2:
                    continue
                if len({h % 2, k % 2, ll % 2}) == 1:
                    out.append((h, k, ll))
    return Q0 * np.array(out, dtype=np.float64)


@pytest.fixture(scope="module")
def crystal() -> RI.CrystalModel:
    sym = RI.proper_rotations(create_cubic_symmetry(A_LAT))
    return RI.build_crystal(fcc_vectors(), sym, q_max=8.0)


def rand_rot(seed: int) -> np.ndarray:
    return np.asarray(Rotation.random(random_state=seed).as_matrix())


def make_rays(R: np.ndarray, crystal: RI.CrystalModel, n: int, seed: int) -> np.ndarray:
    rng = np.random.default_rng(seed)
    idx = rng.choice(len(crystal.h), n, replace=False)
    return crystal.h[idx] @ R.T


def clutter(crystal: RI.CrystalModel, n: int, seed: int) -> np.ndarray:
    rng = np.random.default_rng(seed)
    d = rng.normal(size=(n, 3))
    d /= np.linalg.norm(d, axis=1, keepdims=True)
    return d * rng.choice(crystal.q, n)[:, None]


def noisy(g: np.ndarray, deg: float, seed: int) -> np.ndarray:
    rng = np.random.default_rng(seed)
    ax = rng.normal(size=g.shape)
    ax /= np.linalg.norm(ax, axis=1, keepdims=True)
    rot = Rotation.from_rotvec(ax * np.radians(rng.uniform(0, deg, (len(g), 1))))
    return np.asarray(rot.apply(g) * (1 + rng.uniform(-0.002, 0.002, (len(g), 1))))


def small_g(crystal: RI.CrystalModel) -> np.ndarray:
    return np.concatenate(
        [noisy(make_rays(rand_rot(2), crystal, 25, 3), 0.3, 4), clutter(crystal, 20, 5)]
    )


def big_g(crystal: RI.CrystalModel) -> np.ndarray:
    return np.concatenate(
        [noisy(make_rays(rand_rot(7), crystal, 45, 6), 0.3, 7), clutter(crystal, 400, 8)]
    )


def fake_step(
    g: np.ndarray, crystal: RI.CrystalModel, ctl: PR.Controller, none_first: bool = False
) -> Tuple[Any, List[Any]]:
    """A seed step stand-in: index ``g`` through the wrapped indexer, 7 evaluations per step."""
    seen: List[Any] = []

    def index_voxel(
        sc: Any, verts: Any, cr: Any, p: Any, top_k: int = 20, keep: bool = False
    ) -> Any:
        res = RI.index_rays(g, cr, p, top_k)
        return res, np.arange(len(g))

    wrapped = PR.make_index_voxel(index_voxel, ctl, "pxf")

    def step() -> Any:
        res, _ = wrapped(None, None, crystal, RI.IndexParams(), 20)
        seen.append(res)
        no = none_first and len(seen) == 1
        return SimpleNamespace(
            orientation=None if no else res.R[0],
            evals=7,
            total_s=1.0,
            candidates=res.R,
            step_no=len(seen),
        )

    return step, seen


def hit_by_stage(ctl: PR.Controller, first: float, second: float) -> Any:
    return lambda out: first if out.step_no == 1 else second


def test_no_fallback_at_or_under_50000_pairs(crystal: RI.CrystalModel) -> None:
    ctl = PR.Controller()
    step, seen = fake_step(small_g(crystal), crystal, ctl)
    out, row = PR.run_pxf_step(step, hit_by_stage(ctl, 0.1, 0.1), ctl)
    assert seen[0].n_pairs < 50_000
    assert not row["fallback"] and row["failed_first"] and not row["subsampled"]
    assert len(seen) == 1 and len(ctl.log) == 1 and out.evals == 7


def test_no_fallback_when_first_step_passes(crystal: RI.CrystalModel) -> None:
    ctl = PR.Controller()
    step, seen = fake_step(big_g(crystal), crystal, ctl)
    out, row = PR.run_pxf_step(step, hit_by_stage(ctl, 0.95, 0.1), ctl)
    assert row["subsampled"] and not row["failed_first"] and not row["fallback"]
    assert len(seen) == 1 and ctl.log[0]["n_pairs"] == 50_000 and out.evals == 7


def test_fallback_equals_direct_full_step(crystal: RI.CrystalModel) -> None:
    g = big_g(crystal)
    ctl = PR.Controller()
    step, seen = fake_step(g, crystal, ctl)
    out, row = PR.run_pxf_step(step, hit_by_stage(ctl, 0.5, 0.9), ctl)
    direct = RI.index_rays(g, crystal, RI.IndexParams())  # 600,000: the old setting
    assert row["fallback"] and row["passed_second"]
    assert len(seen) == 2 and seen[0].n_pairs == 50_000 and seen[1].n_pairs == direct.n_pairs
    assert direct.n_pairs > 50_000
    assert np.array_equal(out.candidates, direct.R)  # the second step is the old step
    assert (seen[1].n_pairs, seen[1].n_votes) == (direct.n_pairs, direct.n_votes)
    assert [r["stage"] for r in ctl.log] == ["first", "fallback"]
    assert [r["max_pairs"] for r in ctl.log] == [50_000, 600_000]
    assert ctl.max_pairs == 50_000 and ctl.stage == "first"  # reset for the next voxel


def test_both_steps_evaluations_counted(crystal: RI.CrystalModel) -> None:
    ctl = PR.Controller()
    step, _ = fake_step(big_g(crystal), crystal, ctl)
    out, row = PR.run_pxf_step(step, hit_by_stage(ctl, 0.5, 0.9), ctl)
    assert out.evals == 14 and out.total_s == 2.0
    assert row["first"]["evals"] == 7 and row["second"]["evals"] == 7
    assert row["first"]["n_pairs"] == 50_000 and row["second"]["n_pairs"] > 50_000
    assert row["second"]["n_votes"] > row["first"]["n_votes"]
    assert PR.summed(ctl.log)["pairs"] == row["first"]["n_pairs"] + row["second"]["n_pairs"]


def test_second_fails_too_hands_on_second(crystal: RI.CrystalModel) -> None:
    ctl = PR.Controller()
    step, seen = fake_step(big_g(crystal), crystal, ctl)
    out, row = PR.run_pxf_step(step, hit_by_stage(ctl, 0.5, 0.6), ctl)
    assert row["fallback"] and not row["passed_second"]
    assert out.candidates is seen[1].R or np.array_equal(out.candidates, seen[1].R)
    assert out.evals == 14


def test_no_candidate_after_subsample_falls_back(crystal: RI.CrystalModel) -> None:
    ctl = PR.Controller()
    step, seen = fake_step(big_g(crystal), crystal, ctl, none_first=True)
    out, row = PR.run_pxf_step(step, hit_by_stage(ctl, 0.0, 0.9), ctl)
    assert row["first"]["hit"] is None and row["failed_first"] and row["fallback"]
    assert len(seen) == 2 and out.orientation is not None


class FakePick:
    """A stand-in for ``DeferPick``: ``e1_step`` leaves ``last_step`` and one ``e1_log`` row."""

    def __init__(self, crystal: RI.CrystalModel, g: np.ndarray, ctl: PR.Controller) -> None:
        self.e1_log: List[Dict[str, Any]] = []
        self.last_step: Dict[str, Any] = {}
        self.cur, self.reseed = 12, False
        self.calls = 0
        self.step_fn, _ = fake_step(g, crystal, ctl)

    def e1_step(self, vertices: Any, phase_index: int, n_cost: int = 5, n_refine: int = 2,
                tau: float = -1) -> Any:  # fmt: skip
        out = self.step_fn()
        self.calls += 1
        self.e1_log.append(dict(voxel=self.cur, call=self.calls))
        self.last_step = dict(voxel=self.cur, call=self.calls, kept=[out.orientation])
        return out


def test_rescue_sees_second_steps_ends_and_rows_are_tagged(crystal: RI.CrystalModel) -> None:
    ctl = PR.Controller()
    pick = FakePick(crystal, big_g(crystal), ctl)
    log: List[Dict[str, Any]] = []
    hits = {1: 0.4, 2: 0.9}
    wrapped = PR.make_e1_step(
        FakePick.e1_step,
        ctl,
        log,
        hit_of=lambda self, out, v, ph: hits[self.calls],
        threshold_of=lambda self: 0.80,
    )
    out = wrapped(pick, None, 0, n_cost=20, n_refine=5)
    assert pick.calls == 2
    assert pick.last_step["call"] == 2  # G1 starts from the second step's kept ends
    assert [x["pxf_stage"] for x in pick.e1_log] == ["first", "fallback"]
    assert log[0]["fallback"] and log[0]["voxel"] == 12 and out.evals == 14


def test_wrapper_without_fallback_leaves_one_row(crystal: RI.CrystalModel) -> None:
    ctl = PR.Controller()
    pick = FakePick(crystal, small_g(crystal), ctl)
    log: List[Dict[str, Any]] = []
    wrapped = PR.make_e1_step(FakePick.e1_step, ctl, log, lambda s, o, v, p: 0.1, lambda s: 0.80)
    out = wrapped(pick, None, 0)
    assert pick.calls == 1 and out.evals == 7 and not log[0]["fallback"]
    assert [x["pxf_stage"] for x in pick.e1_log] == ["first"]


def test_wrapper_passes_setting_and_logs_in_both_modes() -> None:
    seen: List[int] = []

    def fake(sc: Any, verts: Any, crystal: Any, p: Any, top_k: int = 20, keep: bool = False) -> Any:
        seen.append(p.max_pairs)
        return SimpleNamespace(n_rays=5, n_pairs=7, n_votes=11), np.arange(5)

    for mode, want in (("old", 600_000), ("pxf", 50_000)):
        ctl = PR.Controller()
        w = PR.make_index_voxel(fake, ctl, mode)
        w(None, None, None, RI.IndexParams(), 20)
        assert seen[-1] == want
        assert ctl.log[0]["n_pairs"] == 7 and ctl.log[0]["max_pairs"] == want
    with pytest.raises(ValueError):
        PR.make_index_voxel(fake, PR.Controller(), "fast")


def test_hexagonal_runs_with_the_fallback() -> None:
    from pymatgen.core.lattice import Lattice
    from pymatgen.core.structure import Structure

    from icenine.symmetry import CrystalSymmetry

    st = Structure(
        Lattice.hexagonal(2.95, 4.68), ["Ti", "Ti"], [[1 / 3, 2 / 3, 0.25], [2 / 3, 1 / 3, 0.75]]
    )
    symmetry = CrystalSymmetry(st)
    rl = st.lattice.reciprocal_lattice
    r = range(-6, 7)
    hkl = np.array([(h, k, ll) for h in r for k in r for ll in r if (h, k, ll) != (0, 0, 0)])
    vec = rl.get_cartesian_coords(hkl.astype(np.float64))
    vec = vec[np.linalg.norm(vec, axis=1) <= 5.0]
    cr = RI.build_crystal(vec, RI.proper_rotations(symmetry), q_max=5.0)
    R = rand_rot(41)
    big = np.concatenate([noisy(make_rays(R, cr, 40, 6), 0.3, 7), clutter(cr, 600, 8)])
    ctl = PR.Controller()
    step, seen = fake_step(big, cr, ctl)
    out, row = PR.run_pxf_step(step, hit_by_stage(ctl, 0.2, 0.9), ctl)
    assert row["fallback"] and len(seen) == 2 and len(out.candidates) >= 1
    assert seen[0].n_pairs == 50_000 and seen[1].n_pairs > 50_000


def test_judges() -> None:
    assert PR.judge_c1(50, 50, 100, 100)["ok"] and not PR.judge_c1(51, 50, 100, 100)["ok"]
    assert PR.judge_c2(24_497_957)["ok"] and not PR.judge_c2(25_000_000)["ok"]
    assert PR.judge_c3(3, 0.001, 2)["ok"] and not PR.judge_c3(4, 0.0, 0)["ok"]
    assert PR.judge_t2(6_164)["ok"] and not PR.judge_t2(6_165)["ok"]
    assert PR.judge_h1(0.2, 5, 0.005)["ok"] and not PR.judge_h1(0.19, 0, 0.0)["ok"]
    assert PR.fewer(40, 100) == pytest.approx(0.6)
    assert PR.STORED_TOTALS == dict(mid=7_383_368, held=17_114_589)
