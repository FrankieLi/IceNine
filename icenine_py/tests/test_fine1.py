"""Tests of the pure parts of fine1 (a whole layer one generation finer): the grid split on
hand-made triangles, the gradient cap, the neighbour-pair measure and the guards."""

import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts" / "autoresearch" / "fine1"))
import fine1_lib as F  # noqa: E402

from icenine.gradient_sample import gradient_orientations, rotation_matrices  # noqa: E402

S = 0.6 / 64  # the side of the generation-6 grid


def _centres(xy: np.ndarray, up: np.ndarray, side: float) -> np.ndarray:
    return F.triangle_vertices(xy, up, side).mean(axis=1)


# ------------------------------------------------------------------------------------ the split


def test_split_up_triangle_has_the_cpp_children() -> None:
    xy, up, c = F.split_once(np.array([[0.0, 0.0]]), np.array([True]), 1.0)
    h = np.sqrt(3.0) / 4.0
    assert c == 0.5 and up.tolist() == [True, True, True, False]
    assert np.allclose(xy, [[0, 0], [0.5, 0], [0.25, h], [0.25, h]])


def test_split_down_triangle_has_the_cpp_children() -> None:
    xy, up, c = F.split_once(np.array([[0.0, 0.0]]), np.array([False]), 1.0)
    h = np.sqrt(3.0) / 4.0
    assert up.tolist() == [False, False, False, True]
    assert np.allclose(xy, [[0, 0], [0.5, 0], [0.25, -h], [0.25, -h]])


@pytest.mark.parametrize("pointing_up", [True, False])
def test_children_tile_the_parent(pointing_up: bool) -> None:
    pxy, pup = np.array([[0.3, -0.2]]), np.array([pointing_up])
    cxy, cup, cs = F.split_once(pxy, pup, 0.8)
    chk = F.grid_checks(pxy, pup, 0.8, cxy, cup, cs)
    assert chk["ok"] and chk["n_child"] == 4 and chk["area_rel_diff"] < 1e-12
    # the four children are pairwise different triangles
    cen = _centres(cxy, cup, cs)
    assert len({tuple(np.round(p, 9)) for p in cen}) == 4


def test_split_equals_the_literal_cpp_port_for_one_generation() -> None:
    rng = np.random.default_rng(3)
    xy = rng.uniform(-0.3, 0.3, size=(7, 2))
    up = rng.random(7) < 0.5
    cxy, cup, _ = F.split_once(xy, up, 0.6 / 64)
    k = 0
    for i in range(7):
        kids = F.cpp_subdivide(xy[i, 0], xy[i, 1], bool(up[i]), 0.6 / 64, 6, 7)
        assert len(kids) == 4
        for kx, ky, ku in kids:
            assert abs(kx - cxy[k, 0]) < 1e-12 and abs(ky - cxy[k, 1]) < 1e-12
            assert ku == bool(cup[k])
            k += 1


def test_split_twice_equals_the_cpp_port_for_two_generations() -> None:
    xy, up = np.array([[0.1, 0.2], [0.0, 0.0]]), np.array([True, False])
    a_xy, a_up, a_s = F.split_once(xy, up, 0.4)
    b_xy, b_up, _ = F.split_once(a_xy, a_up, a_s)
    got = sorted((round(x, 9), round(y, 9), bool(u)) for (x, y), u in zip(b_xy, b_up))
    want = []
    for i in range(2):
        for kx, ky, ku in F.cpp_subdivide(xy[i, 0], xy[i, 1], bool(up[i]), 0.4, 3, 5):
            want.append((round(kx, 9), round(ky, 9), ku))
    assert len(got) == 32 and got == sorted(want)


def test_grid_checks_refuse_a_wrong_count_and_a_child_outside() -> None:
    pxy, pup = np.array([[0.0, 0.0], [0.5, 0.0]]), np.array([True, False])
    cxy, cup, cs = F.split_once(pxy, pup, 0.5)
    assert F.grid_checks(pxy, pup, 0.5, cxy, cup, cs)["ok"]
    assert not F.grid_checks(pxy, pup, 0.5, cxy[:-1], cup[:-1], cs)["ok"]  # count
    moved = cxy.copy()
    moved[5] += 0.01
    bad = F.grid_checks(pxy, pup, 0.5, moved, cup, cs)
    assert not bad["children_inside_parent"] and not bad["ok"]
    assert not F.grid_checks(pxy, pup, 0.5, cxy, cup, cs * 0.9)["ok"]  # area


def test_generation_7_side_is_found_by_the_cpp_formula() -> None:
    """C++ finds the generation as ceil(log(0.6 / side) / log 2); the grid is accepted as is
    (nothing is split again) if that is 7 for MinSideLength 0.0046875."""
    gen = np.ceil(np.log(0.6 / 0.0046875) / np.log(2.0))
    assert int(gen) == 7
    assert int(np.ceil(np.log(0.6 / 0.009375) / np.log(2.0))) == 6


# ------------------------------------------------------------------------- the neighbour pairs


def test_edge_pairs_of_one_split_triangle_are_the_three_corner_middle_pairs() -> None:
    cxy, cup, cs = F.split_once(np.array([[0.0, 0.0]]), np.array([True]), 1.0)
    pairs = F.edge_pairs(_centres(cxy, cup, cs), cs)
    assert sorted(map(tuple, np.sort(pairs, axis=1).tolist())) == [(0, 3), (1, 3), (2, 3)]


def test_edge_pairs_skip_vertex_and_same_direction_neighbours() -> None:
    # two up triangles side by side (same direction, one side apart) share only a vertex
    xy, up = np.array([[0.0, 0.0], [1.0, 0.0]]), np.array([True, True])
    assert len(F.edge_pairs(_centres(xy, up, 1.0), 1.0)) == 0
    # an up and a down triangle across an edge share it
    xy, up = np.array([[0.0, 0.0], [0.5, np.sqrt(3) / 2]]), np.array([True, False])
    assert len(F.edge_pairs(_centres(xy, up, 1.0), 1.0)) == 1


def test_edge_pairs_of_a_larger_patch_count_interior_edges() -> None:
    # two splits of one triangle: 16 triangles, (3 * 16 - 12 boundary edges) / 2 = 18 shared edges
    xy, up, s = F.split_once(np.array([[0.0, 0.0]]), np.array([True]), 1.0)
    xy, up, s = F.split_once(xy, up, s)
    pairs = F.edge_pairs(_centres(xy, up, s), s)
    # the brute-force count: two triangles that share two vertices share an edge
    verts = F.triangle_vertices(xy, up, s)
    key = lambda v: tuple(np.round(v, 9))  # noqa: E731
    sets = [{key(p) for p in tri} for tri in verts]
    want = {(i, j) for i in range(16) for j in range(i + 1, 16) if len(sets[i] & sets[j]) == 2}
    assert {tuple(sorted(p)) for p in pairs.tolist()} == want and len(want) == 18


def test_relative_angle_of_two_rotations_about_z() -> None:
    a = rotation_matrices(np.array([[0.0, 0.0, np.radians(0.2)]]))
    b = rotation_matrices(np.array([[0.0, 0.0, np.radians(-0.15)]]))
    assert F.relative_angle_deg(a, b)[0] == pytest.approx(0.35, abs=1e-9)
    assert F.relative_angle_deg(a, a)[0] == pytest.approx(0.0, abs=1e-9)
    tiny = rotation_matrices(np.array([[0.0, 1e-9, 0.0]]))
    assert F.relative_angle_deg(tiny, np.eye(3)[None])[0] == pytest.approx(
        np.degrees(1e-9), rel=1e-6
    )


# ------------------------------------------------------------------------------------ the cap


def _line_grain(n: int = 30, step: float = 0.01) -> tuple:
    """One grain of ``n`` voxels on a line, plus edge pairs between neighbours on the line."""
    xy = np.stack([np.arange(n) * step, np.zeros(n)], axis=1)
    pairs = np.stack([np.arange(n - 1), np.arange(1, n)], axis=1)
    return xy, np.zeros(n, dtype=int), pairs


def test_cap_leaves_a_small_gradient_equal_to_gradient_orientations() -> None:
    xy, ids, pairs = _line_grain()
    R0 = {0: np.eye(3)}
    a = gradient_orientations(xy, ids, R0, {0: 0.5}, np.random.default_rng(5))
    rng = np.random.default_rng(5)
    b, capped = F.capped_gradient_orientations(xy, ids, R0, {0: 0.5}, rng, pairs)
    assert capped == [] and np.array_equal(a, b)


def test_cap_scales_a_steep_grain_down_to_the_cap() -> None:
    xy, ids, pairs = _line_grain(n=10)
    R0 = {0: np.eye(3)}
    out, capped = F.capped_gradient_orientations(
        xy, ids, R0, {0: 4.0}, np.random.default_rng(1), pairs, cap_deg=0.45
    )
    assert len(capped) == 1
    rec = capped[0]
    assert rec["level_drawn_deg"] == 4.0 and rec["max_step_drawn_deg"] >= 0.45
    assert rec["max_step_kept_deg"] == pytest.approx(0.45, abs=1e-9)
    assert rec["max_step_kept_deg"] <= 0.45 and rec["scale"] < 1.0
    steps = F.relative_angle_deg(out[pairs[:, 0]], out[pairs[:, 1]])
    assert steps.max() == pytest.approx(0.45, abs=1e-9) and steps.max() <= 0.45
    # the kept level is the angle to the centre orientation
    spread = F.relative_angle_deg(out, np.broadcast_to(np.eye(3), out.shape)).max()
    assert rec["level_kept_deg"] == pytest.approx(spread, abs=1e-9)


def test_cap_keeps_the_random_draws_and_ignores_pairs_between_grains() -> None:
    xy = np.concatenate([_line_grain(10)[0], _line_grain(10)[0] + [0.0, 1.0]])
    ids = np.repeat([0, 1], 10)
    pairs = np.array([[0, 1], [10, 11], [9, 10]])  # the last joins two grains
    R0 = {0: np.eye(3), 1: np.eye(3)}
    rng_a, rng_b = np.random.default_rng(2), np.random.default_rng(2)
    gradient_orientations(xy, ids, R0, {0: 4.0, 1: 0.2}, rng_a)
    _, capped = F.capped_gradient_orientations(xy, ids, R0, {0: 4.0, 1: 0.2}, rng_b, pairs)
    assert rng_a.standard_normal() == rng_b.standard_normal()  # same number of draws
    assert [c["id"] for c in capped] == [0]


def test_cap_triggers_at_exactly_the_cap_and_not_below() -> None:
    xy, ids, pairs = _line_grain(n=2, step=1.0)  # one pair: its step is the whole spread range
    R0 = {0: np.eye(3)}
    out, capped = F.capped_gradient_orientations(
        xy, ids, R0, {0: 0.9}, np.random.default_rng(0), pairs, cap_deg=0.45
    )
    # two voxels: each is level 0.9 from the centre at most; the pair step is large -> capped
    assert len(capped) == 1
    out2, capped2 = F.capped_gradient_orientations(
        xy, ids, R0, {0: 0.2}, np.random.default_rng(0), pairs, cap_deg=0.45
    )
    assert capped2 == []


# ------------------------------------------------------------------------------------ the gate


def test_gate_needs_both_conditions_and_some_pairs() -> None:
    assert F.gate(1000, 0.45, 4.0)["passed"]
    assert not F.gate(1000, 0.5, 4.0)["passed"]  # the limit is strict
    assert not F.gate(1000, 0.45, 2.99)["passed"]
    assert not F.gate(0, 0.0, 4.0)["passed"]  # no pairs: nothing was measured


def test_percentiles_and_histogram_account_for_every_pair() -> None:
    v = np.linspace(0.0, 0.499, 1000)
    h = F.histogram(v)
    assert sum(b["count"] for b in h) == 1000 and h[-1]["count"] == 0
    p = F.percentiles(v)
    assert p["p50"] == pytest.approx(np.median(v)) and p["p99.9"] > p["p99"] > p["p90"]


def test_nearest_centre_rows_are_sorted_and_the_right_count() -> None:
    rng = np.random.default_rng(0)
    xy = rng.uniform(-1, 1, size=(500, 2))
    rows = F.nearest_centre_rows(xy, 50)
    assert len(rows) == 50 and np.all(np.diff(rows) > 0)
    d = np.linalg.norm(xy, axis=1)
    assert d[rows].max() <= np.sort(d)[50]


# ------------------------------------------------------------------------- configs, runs, launcher

import subprocess  # noqa: E402

import fine1_config as CFG  # noqa: E402
import fine1_freeze as FZ  # noqa: E402
import fine1_paths as PP  # noqa: E402
import fine1_run_py as RP  # noqa: E402
import fine1_stack as ST  # noqa: E402

ICE = Path(__file__).resolve().parents[1]
LAUNCH = ICE / "scripts" / "autoresearch" / "fine1" / "fine1_launch.sh"
FRAMES = "/x/fine7/clean/500Grains.sim"


def test_ps_config_changes_only_the_five_keys() -> None:
    text = CFG.python_config("PS", Path("/x/grid.mic"), FRAMES, "out")
    assert CFG.changed_keys(CFG.BASE_PS.read_text(), text) == sorted(CFG.CHANGED_KEYS)
    k = CFG.keys_of(text)
    assert k["MinSideLength"] == k["MaxInitSideLength"] == "0.0046875"
    assert {"BFSRevisitRefit", "BFSReseed", "FusedCost"} <= set(k)


def test_j5s_config_adds_the_twelve_library_keys_of_lib_check() -> None:
    text = CFG.python_config("J5S", Path("/x/grid.mic"), FRAMES, "out")
    added = {s.split()[0] for s in CFG.J5S_KEYS}
    assert CFG.changed_keys(CFG.BASE_J5S.read_text(), text) == sorted(set(CFG.CHANGED_KEYS) | added)
    assert CFG.J5S_KEYS == CFG.lib_check_keys()  # step 21's key list, same order
    k = CFG.keys_of(text)
    assert k["SeedByIndexing"] == "1" and "BoundaryPassSkipHitRatio" not in k


def test_j5s_config_parses_and_satisfies_the_screen_preconditions(tmp_path: Path) -> None:
    from icenine.config_file import ConfigFile
    from icenine.reconstructor import check_bfs_keys

    p = tmp_path / "j5s.config"
    p.write_text(CFG.python_config("J5S", Path("/x/grid.mic"), FRAMES, "out"))
    cfg = ConfigFile.from_file(str(p))
    assert cfg.boundary_pass and cfg.boundary_pass_repair and cfg.boundary_pass_screen
    assert cfg.boundary_pass_skip_hit == 0.0 and cfg.seed_by_indexing and cfg.fused_cost
    check_bfs_keys(cfg)  # raises on a violated precondition


def test_cpp_config_drops_the_python_lines_and_sets_the_fine_side() -> None:
    text = CFG.cpp_config(Path("/x/grid.mic"), FRAMES, "outc")
    k = CFG.keys_of(text)
    assert not ({"BFSRevisitRefit", "BFSReseed", "FusedCost"} & set(k))
    assert k["MinSideLength"] == k["MaxInitSideLength"] == "0.0046875"
    assert k["SampleFilename"] == "/x/grid.mic" and k["InfileBasename"] == FRAMES
    assert k["OutStructureBasename"] == "outc"
    # the search lines are those of the Python arm PS
    ps = CFG.keys_of(CFG.python_config("PS", Path("/x/grid.mic"), FRAMES, "outc"))
    for key in ("MaxMCSteps", "MinLocalResolution", "MaxQ", "SuccessiveRestarts", "EtaLimit"):
        assert k[key] == ps[key]


def test_counts_summary_splits_global_searches_from_the_rest() -> None:
    st = dict(
        n_seeds=3,
        n_reseeds=2,
        n_evals_seed=300,
        n_evals_reseed=200,
        n_evals_neighbor=50,
        n_evals_revisit=10,
        n_evals_boundary_pass=5,
        n_evals_boundary_repair=7,
        wall_seed_s=30.0,
        wall_reseed_s=20.0,
    )
    s = RP.counts_summary(st, 100.0)
    assert s["evals_all"] == 572 and s["evals_global"] == 500 and s["evals_other"] == 72
    assert s["global_searches"] == 5 and s["wall_per_global_search_s"] == 10.0
    assert s["wall_other_s"] == 50.0


def test_freeze_check_refuses_a_missing_or_uncommitted_file(tmp_path: Path, monkeypatch) -> None:
    monkeypatch.setattr(PP, "FROZEN", tmp_path / "frozen.json")
    assert "does not exist" in FZ.problems(False, False)
    (tmp_path / "frozen.json").write_text("{}")
    assert "not committed" in FZ.problems(False, False)


def test_stack_check_reads_the_ascii_frame_as_row_column(tmp_path: Path) -> None:
    f = tmp_path / "frame"
    f.write_text("886, 1344, 158593\n887, 1344, 0\n5, 7, 3\n")
    lit = ST.ascii_lit(f)
    assert sorted(map(tuple, lit.tolist())) == [(7, 5), (1344, 886)]  # (k = row, j = column)


def test_launcher_syntax_and_argument_guards() -> None:
    assert subprocess.run(["bash", "-n", str(LAUNCH)]).returncode == 0
    bad = subprocess.run(["bash", str(LAUNCH), "--stage", "mid"], capture_output=True)
    assert bad.returncode == 64
    nocap = subprocess.run(["bash", str(LAUNCH), "--stage", "whole"], capture_output=True)
    assert nocap.returncode == 64 and b"cap-h" in nocap.stderr
    unk = subprocess.run(["bash", str(LAUNCH), "--stage", "slice", "--bogus"], capture_output=True)
    assert unk.returncode == 64
