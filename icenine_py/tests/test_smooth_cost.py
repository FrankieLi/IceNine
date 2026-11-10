"""Task E: smooth ranking costs (icenine/smooth_cost.py) and the SmoothCostRank key.

The key parses and defaults to off; off leaves the level loop unchanged (no ranker is built, the
recorded search is identical to a run without the attribute, and equal bit for bit to a fixture
recorded on the commit before Task E); on, the coarse levels are ordered by the smooth key, and the
level's recorded best stays the best by the strict cost. The arms are checked on the ThreeVoxels toy
problem (real two-detector images): the control arm S reproduces the strict cost near the truth, the
smoothed arms change in small increments along a fine ray where the strict cost steps, the shared
multi-spec evaluation equals separate calls, and the E4 autograd gradient matches finite
differences."""

import json
import math
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
import pytest
import torch
from scipy.spatial.transform import Rotation

from icenine import cost_functions as CF
from icenine.orientation_search import SearchParameters
from icenine.reconstructor import _get_voxel_vertices
from icenine.smooth_cost import SmoothCost, SmoothSpec, smoothstep

needs_c = pytest.mark.skipif(not CF._HAS_C_RASTERIZE, reason="C extension missing")
ROOT = Path(__file__).resolve().parents[1]
KEYOFF_FIXTURE = ROOT / "tests" / "data" / "smooth_keyoff_ref.npz"
sys.path.insert(0, str(ROOT / "scripts" / "eval_reduction"))


# ------------------------------------------------------------------------------ spec and config
def test_spec_parse() -> None:
    assert not SmoothSpec.parse("off").active
    assert not SmoothSpec.parse("").active
    s = SmoothSpec.parse("E1+E2,sigma=8,omega=auto,levels=2")
    assert s.arms == frozenset({"E1", "E2"}) and s.sigma_px == 8.0 and s.omega_half == 0.0
    assert s.levels == 2 and s.label() == "E1+E2"
    assert SmoothSpec.parse("e4,tau=0.25").tau_deg == 0.25
    for bad in ("E9", "E1,sigma=-1", "E1,banana=1", "E1,sigma", "E1, sigma=2", "E1,levels=0"):
        with pytest.raises(ValueError, match="SmoothCostRank"):
            SmoothSpec.parse(bad)


def test_config_key_default_off_and_parse(tmp_path: Any) -> None:
    from icenine.config_file import ConfigFile

    base = Path(__file__).resolve().parents[1] / "scripts" / "phase_d" / "configs"
    text = (base / "ReconstructPhaseD_mc_clean.config").read_text()
    sp0 = SearchParameters.from_config(
        ConfigFile.from_file(str(base / "ReconstructPhaseD_mc_clean.config"))
    )
    assert sp0.smooth_rank == "off"
    f = tmp_path / "on.config"
    f.write_text(text + "\nSmoothCostRank E1+E2+E3,sigma=auto\n")
    sp = SearchParameters.from_config(ConfigFile.from_file(str(f)))
    assert SmoothSpec.parse(sp.smooth_rank).arms == frozenset({"E1", "E2", "E3"})
    f.write_text(text + "\nSmoothCostRank E7\n")
    with pytest.raises(ValueError, match="SmoothCostRank"):
        ConfigFile.from_file(str(f))
    with pytest.raises(ValueError, match="SmoothCostRank"):
        SearchParameters(smooth_rank="E1,sigma=0")


def test_smoothstep() -> None:
    x = torch.tensor([-1.0, 0.0, 0.5, 1.0, 2.0])
    assert torch.allclose(smoothstep(x), torch.tensor([0.0, 0.0, 0.5, 1.0, 1.0]))


# ------------------------------------------------------------------------------ toy problem
def _toy() -> Tuple[Any, Any, np.ndarray]:
    from test_findoptimal_refactor import _build

    return _build()


def _sc(rec: Any) -> SmoothCost:
    s = rec.setup
    return SmoothCost(
        s.simulator,
        s.detector_list,
        s.range_map,
        s.exp_data,
        s.sample,
        s.structure_list,
        eta_limit=s.exp_setup.get_eta_limit(),
    )


def _strict(rec: Any) -> Any:
    s = rec.setup
    return CF.VoxelCostFunction(
        simulator=s.simulator,
        detector_list=s.detector_list,
        range_map=s.range_map,
        exp_data=s.exp_data,
        sample=s.sample,
        structure_list=s.structure_list,
        eta_limit=s.exp_setup.get_eta_limit(),
        pixel_radius=0,
    )


def _ray(R_true: np.ndarray, radii_deg: np.ndarray, seed: int = 4) -> np.ndarray:
    d = np.random.default_rng(seed).normal(size=3)
    d /= np.linalg.norm(d)
    return Rotation.from_rotvec(np.radians(radii_deg)[:, None] * d).as_matrix() @ R_true


@needs_c
def test_control_arm_tracks_strict_cost_near_truth() -> None:
    rec, voxel, R_true = _toy()
    sc, cf = _sc(rec), _strict(rec)
    verts = _get_voxel_vertices(voxel)
    Rs = _ray(R_true, np.array([0.0, 0.02, 0.05, 0.1, 0.2]))
    strict = np.array([cf.evaluate(R.astype(np.float32), verts, voxel.phase).cost for R in Rs])
    s = sc.evaluate_many(Rs, verts, voxel.phase, SmoothSpec.parse("S"))
    # identical at the truth; four sample points against the full rasterisation drift apart as
    # the spots move (0.04 at 0.05 deg, about 0.1 at 0.1 deg on this problem)
    assert abs(s[0] - strict[0]) < 1e-9
    assert np.all(np.abs(s[:3] - strict[:3]) < 0.05), (s, strict)
    # both rise away from the truth
    assert s[-1] > s[0] + 0.2 and strict[-1] > strict[0] + 0.2


@needs_c
def test_smoothed_arms_change_gradually_where_the_strict_cost_steps() -> None:
    rec, voxel, R_true = _toy()
    sc = _sc(rec)
    verts = _get_voxel_vertices(voxel)
    radii = np.linspace(0.0, 0.4, 801)
    Rs = _ray(R_true, radii)
    names = ["S", "E1,sigma=3", "E1+E2,sigma=3,omega=1", "E1+E2+E3,sigma=3,omega=1"]
    out = sc.evaluate_specs(Rs, verts, voxel.phase, [SmoothSpec.parse(n) for n in names])
    steps = {n: np.abs(np.diff(c)) for n, c in zip(names, out)}
    big = steps["S"].max()
    assert big > 0.01  # the control steps
    # E1 alone smooths the pixel edge only (frame crossings still step): more, smaller steps
    assert (steps["E1,sigma=3"] > 0).sum() > (steps["S"] > 0).sum()
    # with omega blending as well, no step above a quarter of the control's largest
    for n in names[2:]:
        assert steps[n].max() < 0.25 * big, (n, steps[n].max(), big)


@needs_c
def test_evaluate_specs_equals_separate_calls() -> None:
    rec, voxel, R_true = _toy()
    sc = _sc(rec)
    verts = _get_voxel_vertices(voxel)
    Rs = _ray(R_true, np.linspace(0.0, 2.0, 37), seed=9)
    names = ["S", "E1,sigma=3", "E1,sigma=12", "E2", "E3", "E1+E2+E3", "E4,tau=0.25", "E3+E4"]
    specs = [SmoothSpec.parse(n) for n in names]
    together = sc.evaluate_specs(Rs, verts, voxel.phase, specs, diameter_rad=math.radians(5.0))
    for sp, t in zip(specs, together):
        alone = sc.evaluate_many(Rs, verts, voxel.phase, sp, diameter_rad=math.radians(5.0))
        np.testing.assert_array_equal(t, alone)
        assert np.all((t >= 0.0) & (t <= 1.0))


@needs_c
def test_auto_widths_follow_the_level() -> None:
    rec, _, _ = _toy()
    sc = _sc(rec)
    sp = SmoothSpec.parse("E1+E2+E4")
    w0 = sc.level_widths(sp, math.radians(5.0))
    w2 = sc.level_widths(sp, math.radians(5.0) / 1.5**2)
    assert w0[0] > w2[0] and w0[2] > w2[2] and w0[1] >= w2[1] >= 1.0
    assert sc.level_widths(SmoothSpec.parse("E1,sigma=4"), math.radians(5.0))[0] == 4.0


@needs_c
def test_e4_rays_attached_and_gradient_matches_finite_difference() -> None:
    rec, voxel, R_true = _toy()
    sc = _sc(rec)
    verts = _get_voxel_vertices(voxel)
    rays = sc.observed_rays()
    assert len(rays.frame) > 10
    assert len(sc.attached_rays(verts)) > 10
    c = sc.evaluate_many(
        _ray(R_true, np.array([0.0, 1.0])), verts, voxel.phase, SmoothSpec.parse("E4,tau=0.25")
    )
    assert c[0] < c[1]
    R0 = torch.as_tensor(_ray(R_true, np.array([0.1]), seed=2)[0])

    def f(w: torch.Tensor) -> torch.Tensor:
        th = torch.linalg.norm(w)
        k = w / th
        K = torch.zeros(3, 3, dtype=torch.float64)
        K[0, 1], K[0, 2], K[1, 2] = -k[2], k[1], -k[0]
        K = K - K.T
        Rw = torch.eye(3, dtype=torch.float64) + torch.sin(th) * K + (1 - torch.cos(th)) * K @ K
        return sc.e4_cost_torch(Rw @ R0, verts, voxel.phase, tau_deg=0.25)

    w = torch.tensor([1e-4, -2e-4, 1.5e-4], dtype=torch.float64, requires_grad=True)
    val = f(w)
    (g,) = torch.autograd.grad(val, w)
    h = 1e-7
    fd = []
    for i in range(3):
        e = torch.zeros(3, dtype=torch.float64)
        e[i] = h
        fd.append(float((f(w.detach() + e) - f(w.detach() - e)) / (2 * h)))
    np.testing.assert_allclose(g.numpy(), np.array(fd), rtol=1e-3, atol=1e-3)
    assert float(torch.linalg.norm(g)) > 0


# ------------------------------------------------------------------------------ the level loop
def _run(smooth_rank: Optional[str] = None) -> Tuple[Any, List[Tuple[str, Dict[str, Any]]], Any]:
    """The toy search of the key-off fixture (about 20 candidates at level 0, two levels)."""
    from smooth_keyoff_fixture import toy_run

    return toy_run(smooth_rank)


def _same_events(a: List[Tuple[str, Dict[str, Any]]], b: List[Tuple[str, Dict[str, Any]]]) -> None:
    assert [n for n, _ in a] == [n for n, _ in b]
    for (_, da), (_, db) in zip(a, b):
        assert da.keys() == db.keys()
        for k in da:
            np.testing.assert_array_equal(np.asarray(da[k]), np.asarray(db[k]))


@needs_c
def test_key_off_is_bit_identical_to_the_default() -> None:
    res0, ev0, rec0 = _run()
    res1, ev1, rec1 = _run("off")
    np.testing.assert_array_equal(res0.orientation, res1.orientation)
    assert res0.cost == res1.cost
    _same_events(ev0, ev1)
    assert rec0.last_smooth_ranker is None and rec1.last_smooth_ranker is None


@needs_c
def test_key_off_equals_the_base_commit_fixture_exactly() -> None:
    """The key-off search equals the search recorded on the commit before Task E
    (``tests/data/smooth_keyoff_ref.npz``, ``scripts/eval_reduction/smooth_keyoff_fixture.py``):
    final orientation and cost, every recorder event and ``last_level_best``, bit for bit."""
    import platform

    from smooth_keyoff_fixture import arrays_of
    from test_s1_evaluate import platform_difference

    if not KEYOFF_FIXTURE.exists():
        pytest.skip("tests/data/smooth_keyoff_ref.npz missing")
    ref = np.load(KEYOFF_FIXTURE, allow_pickle=False)
    meta = json.loads(str(ref["meta"]))
    diff = platform_difference(meta, platform.system(), platform.machine(), torch.__version__)
    if diff is not None:
        pytest.skip(f"exact comparison only on the fixture's platform and torch build ({diff})")
    prev = torch.get_num_threads()
    torch.set_num_threads(int(meta["threads"]))
    try:
        res, ev, rec = _run()
        got = arrays_of(res, ev, rec.last_level_best)
        # negative control: the key on (E1+E2+E3 at level 0) changes the recorded level-0 order
        res_on, ev_on, rec_on = _run("E1+E2+E3,levels=1")
        on = arrays_of(res_on, ev_on, rec_on.last_level_best)
    finally:
        torch.set_num_threads(prev)
    assert sorted(got) == sorted(k for k in ref.files if k != "meta")
    perm0 = [k for k in got if k.endswith("_quick_mc/perm")][0]
    assert len(got[perm0]) >= 10  # a list long enough for the ranking to matter
    for k, v in got.items():
        want = ref[k]
        assert v.dtype == want.dtype and np.array_equal(v, want), k
    assert not np.array_equal(on[perm0], ref[perm0])


@needs_c
def test_key_on_orders_the_coarse_level_by_the_smooth_cost() -> None:
    spec = "E1+E2+E3,levels=1"
    res, ev, rec = _run(spec)
    ranker = rec.last_smooth_ranker
    assert ranker is not None and ranker.n_evals > 0
    q = [d for n, d in ev if n == "quick_mc" and d["level"] == 0][0]
    disc = [d for n, d in ev if n == "discrete" and d["level"] == 0][0]
    key = ranker.cost.evaluate_many(
        q["R"],
        ranker.voxel_vertices,
        ranker.phase,
        SmoothSpec.parse(spec),
        diameter_rad=float(disc["diameter"]),
    )
    assert len(key) >= 10 and np.all(np.diff(key) >= 0)  # the recorded order is the smooth key's
    # the strict cost is still what every candidate carries, and its order is not the recorded one
    assert np.all((q["cost"] >= 0) & (q["cost"] <= 1))
    assert not np.all(np.diff(q["cost"]) >= 0)
    assert np.isfinite(res.cost)


@needs_c
def test_level_best_is_the_strict_best_under_smooth_ranking() -> None:
    """At a smooth-ranked level, ``last_level_best`` holds the best candidate by the strict cost
    (the first of equal costs in ranked order), not the first of the smooth-ranked list."""
    res, ev, rec = _run("E1+E2+E3")
    lists = [d for n, d in ev if n == "quick_mc"]
    assert len(lists) == len(rec.last_level_best) == 2
    first_differs = 0
    for q, (lv, R_best, c_best) in zip(lists, rec.last_level_best):
        i = int(np.argmin(q["cost"]))
        assert lv == q["level"] and c_best == float(q["cost"][i])
        np.testing.assert_array_equal(R_best, q["R"][i])
        first_differs += int(q["cost"][0] > q["cost"][i])
    assert first_differs >= 1  # the toy has a level where the smooth key's first is not the best


# ------------------------------------------------------------------------------ SmoothCost cache
def _toy_search(rec: Any, voxel: Any) -> Tuple[Any, List[Tuple[str, Dict[str, Any]]]]:
    import contextlib
    import io

    events: List[Tuple[str, Dict[str, Any]]] = []
    rec.recorder = lambda name, data: events.append((name, data))
    with contextlib.redirect_stdout(io.StringIO()):
        res = rec.reconstruct_voxel(
            _get_voxel_vertices(voxel), voxel.phase, rng=np.random.default_rng(7)
        )
    return res, events


@needs_c
def test_smooth_cost_cache_builds_once_per_image_set() -> None:
    """Same objects and settings: the cached ``SmoothCost`` is returned. A new image-set object,
    another ``min_sin_eta`` or eta limit, or a sample-to-lab transform changed in place give a new
    build."""
    import copy

    from icenine.smooth_cost import SmoothCostCache

    rec, _, _ = _toy()
    s = rec.setup
    eta = s.exp_setup.get_eta_limit()
    cache = SmoothCostCache()
    a = cache.get(s, eta, 0.0)
    assert cache.get(s, eta, 0.0) is a
    assert (cache.n_builds, cache.n_hits) == (1, 1)
    b = cache.get(s, eta, 0.3)  # another min_sin_eta
    assert b is not a and b.min_sin_eta == 0.3
    c = cache.get(s, eta * 0.5, 0.3)  # another eta limit
    assert c is not b and c.eta_limit == pytest.approx(eta * 0.5)
    assert cache.n_builds == 3
    old = s.exp_data
    try:
        s.exp_data = copy.copy(old)  # a new image-set object (equal content)
        d = cache.get(s, eta * 0.5, 0.3)
        assert d is not c and d.exp_data is s.exp_data and cache.n_builds == 4
    finally:
        s.exp_data = old
    e = cache.get(s, eta, 0.0)  # back to the first image set: one more build, then hits
    assert cache.get(s, eta, 0.0) is e and cache.n_builds == 5
    m = s.sample.sample_to_lab_matrix
    keep = m.clone()
    try:
        m[0, 3] += 1e-3  # changed in place
        assert cache.get(s, eta, 0.0) is not e and cache.n_builds == 6
    finally:
        m.copy_(keep)
    cache.clear()
    assert cache.get(s, eta, 0.0) is not e and cache.n_builds == 7


@needs_c
def test_cached_smooth_cost_search_equals_a_fresh_build_bit_for_bit() -> None:
    """Key on (all four arms, so the lit-pixel trees, the observed rays and the attached rays are
    all cached): a second search on the same reconstructor, served by the warm cache after the cost
    has also been used for another voxel position, equals a search whose ``SmoothCost`` is built
    fresh (``smooth_cost_cache`` None, the behaviour before the cache), bit for bit."""
    from smooth_keyoff_fixture import arrays_of, build_toy

    spec = "E1+E2+E3+E4"
    rec, voxel = build_toy()
    rec.setup.search_params.smooth_rank = spec
    res1, ev1 = _toy_search(rec, voxel)
    first = arrays_of(res1, ev1, rec.last_level_best)
    sc = rec.last_smooth_ranker.cost
    sc.attached_rays(_get_voxel_vertices(voxel) + 0.01)  # warm state of another position
    res2, ev2 = _toy_search(rec, voxel)
    second = arrays_of(res2, ev2, rec.last_level_best)
    assert rec.last_smooth_ranker.cost is sc
    assert (rec.smooth_cost_cache.n_builds, rec.smooth_cost_cache.n_hits) == (1, 1)
    assert sc._rays is not None and len(sc._attached) == 2

    fresh_rec, voxel_f = build_toy()
    fresh_rec.setup.search_params.smooth_rank = spec
    fresh_rec.smooth_cost_cache = None
    res3, ev3 = _toy_search(fresh_rec, voxel_f)
    fresh = arrays_of(res3, ev3, fresh_rec.last_level_best)
    assert fresh_rec.last_smooth_ranker.n_evals == rec.last_smooth_ranker.n_evals > 0
    for got in (first, second):
        assert sorted(got) == sorted(fresh)
        for k, v in fresh.items():
            assert got[k].dtype == v.dtype and np.array_equal(got[k], v), k


@needs_c
def test_key_off_builds_no_smooth_cost() -> None:
    from smooth_keyoff_fixture import build_toy

    rec, voxel = build_toy()
    _toy_search(rec, voxel)
    assert rec.last_smooth_ranker is None
    assert rec.smooth_cost_cache is not None and rec.smooth_cost_cache.n_builds == 0
