"""Tests of the tex1 scripts (textured finer layer ``tex24``): pure parts and config text, on small
hand-made inputs; no data file of the sample is read."""

import json
import subprocess
import sys
from pathlib import Path
from typing import Any, Dict

import numpy as np
import pytest

ICE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ICE / "scripts" / "autoresearch" / "tex1"))
sys.path.insert(0, str(ICE / "scripts"))
import plot_mic as PM  # noqa: E402
import tex1_config as CFG  # noqa: E402
import tex1_freeze as FZ  # noqa: E402
import tex1_lib as F  # noqa: E402
import tex1_paths as PP  # noqa: E402
import tex1_plot as PL  # noqa: E402
import tex1_run_py as RP  # noqa: E402
import tex1_sample_check as SCH  # noqa: E402
import tex1_score as SC  # noqa: E402

LAUNCH = ICE / "scripts" / "autoresearch" / "tex1" / "tex1_launch.sh"
FRAMES = "/x/tex24/clean/500Grains.sim"
H = np.sqrt(3.0) / 2.0


def patch(n: int = 12, side: float = 1.0) -> np.ndarray:
    """Centres of a triangle patch: ``n`` by ``n`` up triangles with the down triangle on the right
    edge of each (up centre (x + 0.5, y + h/3), down centre (x + 1, y + 2h/3))."""
    out = []
    for r in range(n):
        for i in range(n):
            x, y = (i + 0.5 * r) * side, r * H * side
            out.append((x + 0.5 * side, y + H * side / 3))
            out.append((x + side, y + 2 * H * side / 3))
    return np.array(out)


# ------------------------------------------------------------------------------- neighbour pairs


def test_interior_triangle_has_three_edge_and_twelve_touching_neighbours() -> None:
    c = patch(12)
    mid = int(np.argmin(np.linalg.norm(c - c.mean(axis=0), axis=1)))
    e, t = F.edge_pairs(c, 1.0), F.touching_pairs(c, 1.0)
    assert int(((e[:, 0] == mid) | (e[:, 1] == mid)).sum()) == 3
    assert int(((t[:, 0] == mid) | (t[:, 1] == mid)).sum()) == 12


def test_edge_pairs_are_a_subset_of_touching_pairs() -> None:
    c = patch(6)
    e = {tuple(p) for p in F.edge_pairs(c, 1.0).tolist()}
    t = {tuple(p) for p in F.touching_pairs(c, 1.0).tolist()}
    assert e < t


# ------------------------------------------------------------------------------- classes and draws


def test_classes_are_20_and_50_percent_drawn_independently_and_reproducibly() -> None:
    ids = list(range(250))
    noise, grad = F.draw_classes(ids, 24, 0.2, 0.5)
    assert (len(noise), len(grad)) == (50, 125)
    assert F.draw_classes(ids, 24, 0.2, 0.5) == (noise, grad)
    assert 0 < len(set(noise) & set(grad)) < 50  # independent: some grains carry both
    n2, g2 = F.draw_classes(ids[:100], 24, 0.2, 0.5)  # another count: 20% and 50% of it
    assert (len(n2), len(g2)) == (20, 50)


def test_class_of_numbers_the_four_classes() -> None:
    assert [F.class_of(a, b) for a, b in ((0, 0), (1, 0), (0, 1), (1, 1))] == [0, 1, 2, 3]


def test_ball_vectors_stay_in_the_ball_and_fill_it() -> None:
    v = F.ball_vectors(20000, np.random.default_rng(0), 1.0)
    r = np.degrees(np.linalg.norm(v, axis=1))
    assert r.max() <= 1.0 and r.max() > 0.99
    assert np.median(r) == pytest.approx(0.5 ** (1 / 3), abs=0.01)  # uniform in volume


def test_scale_to_edge_max_hits_the_target_from_either_side() -> None:
    c = patch(8)
    pairs = F.edge_pairs(c, 1.0)
    rng = np.random.default_rng(1)
    rv = F.field_vectors(c, rng.standard_normal((3, 2)), rng.standard_normal((3, 3)))
    for target in (0.1, 0.5, 3.0):
        for pre in (1e-4, 1e-1):  # a vector set far too small or too large to start
            s, cur = F.scale_to_edge_max(rv * pre, pairs[:, 0], pairs[:, 1], target)
            assert cur == pytest.approx(target, rel=1e-9) and cur <= target
            M = F.rotation_matrices(rv * pre * s)
            assert F.relative_angle_deg(M[pairs[:, 0]], M[pairs[:, 1]]).max() <= target


def two_grain_patch() -> Any:
    c = patch(12)
    ids = (c[:, 0] > np.median(c[:, 0])).astype(int)
    rng = np.random.default_rng(0)
    Rg = {}
    for g in (0, 1):
        a = rng.standard_normal(3)
        Rg[g] = F.rotation_matrices(a / np.linalg.norm(a) * 1.0)[0]
    return c, ids, Rg, F.edge_pairs(c, 1.0)


def test_textured_orientations_respect_the_three_bounds() -> None:
    c, ids, Rg, pairs = two_grain_patch()
    R, Rs, rec = F.textured_orientations(c, ids, Rg, [0], [0, 1], pairs, 24)
    ref = np.stack([Rg[int(g)] for g in ids])
    assert len(rec) == 2 and all(0.1 <= r["t_g_deg"] < 0.5 for r in rec)
    same = ids[pairs[:, 0]] == ids[pairs[:, 1]]
    smooth_step = F.relative_angle_deg(Rs[pairs[same, 0]], Rs[pairs[same, 1]])
    for r in rec:  # the largest edge step of a gradient grain equals its strength
        m = ids[pairs[same, 0]] == r["id"]
        assert smooth_step[m].max() == pytest.approx(r["t_g_deg"], rel=1e-9)
        assert r["edge_max_deg"] == pytest.approx(r["t_g_deg"], rel=1e-9)
    # noise grain 0: within 1 degree of the smooth field; grain 1 equals it
    assert F.relative_angle_deg(R[ids == 0], Rs[ids == 0]).max() <= 1.0
    assert F.relative_angle_deg(R[ids == 0], Rs[ids == 0]).max() > 0.5
    assert np.array_equal(R[ids == 1], Rs[ids == 1])
    # any two voxels of the noise grain: noise parts differ by at most 2 degrees
    in0 = pairs[same][ids[pairs[same, 0]] == 0]
    ang0 = F.relative_angle_deg(R[in0[:, 0]], R[in0[:, 1]])
    assert ang0.max() <= 2.0 + smooth_step.max() + 1e-9 and ang0.max() > 1.0
    assert F.relative_angle_deg(Rs, ref).max() > 0.5  # the smooth field leaves R_g


def test_no_class_means_every_voxel_equals_its_grain_orientation() -> None:
    c, ids, Rg, pairs = two_grain_patch()
    R, Rs, rec = F.textured_orientations(c, ids, Rg, [], [], pairs, 24)
    assert rec == [] and np.array_equal(R, Rs)
    assert np.array_equal(R, np.stack([Rg[int(g)] for g in ids]))


def test_rebuild_from_recorded_field_parameters_equals_the_draw() -> None:
    c, ids, Rg, pairs = two_grain_patch()
    R, _, rec = F.textured_orientations(c, ids, Rg, [0], [0, 1], pairs, 24)
    params = {
        r["id"]: dict(M1=r["M1"], M2=r["M2"], t_g_deg=r["t_g_deg"], scale=r["scale"]) for r in rec
    }
    R2, _, _ = F.textured_orientations(c, ids, Rg, [0], [0, 1], pairs, 24, field_params=params)
    assert F.relative_angle_deg(R, R2).max() < 1e-9


def test_a_gradient_grain_without_edge_pairs_keeps_one_orientation() -> None:
    c = np.array([[0.0, 0.0], [5.0, 0.0]])  # two voxels far apart: no edge pair
    Rg = {0: np.eye(3)}
    _, Rs, rec = F.textured_orientations(c, np.array([0, 0]), Rg, [], [0], np.zeros((0, 2)), 24)
    assert rec[0]["scale"] == 0.0 and np.array_equal(Rs[0], np.eye(3))


def test_nearest_centre_rows_are_sorted_and_the_right_count() -> None:
    xy = np.array([[3.0, 0.0], [0.1, 0.0], [0.0, 0.2], [9.0, 9.0]])
    assert F.nearest_centre_rows(xy, 2).tolist() == [1, 2]


# ------------------------------------------------------------------------------- the sample gate


def gate_inputs() -> Any:
    c, ids, Rg, pairs = two_grain_patch()
    touch = F.touching_pairs(c, 1.0)
    cls = {0: F.NOISE_ONLY, 1: F.GRADIENT_ONLY}
    R, Rs, _ = F.textured_orientations(c, ids, Rg, [0], [1], pairs, 24)
    from icenine.ray_indexing import proper_rotations
    from icenine.symmetry import create_cubic_symmetry

    sym = proper_rotations(create_cubic_symmetry(3.6))
    return R, Rs, ids, cls, Rg, pairs, touch, sym


def test_gate_passes_on_a_textured_patch_and_reports_the_classes() -> None:
    R, Rs, ids, cls, Rg, e, t, sym = gate_inputs()
    res = SCH.check(R, Rs, ids, cls, Rg, e, t, sym)
    assert res["gate"]["passed"], res["gate"]
    assert res["classes"]["noise"]["n_grains"] == 1 and res["classes"]["gradient"]["n_grains"] == 1
    assert res["pairs"]["edge"]["gradient"]["max"] <= 0.5 + F.FILE_ROUNDING_DEG
    assert res["pairs"]["edge"]["noise"]["max"] > 0.5


def test_gate_fails_on_a_bad_sharp_voxel_a_wide_noise_voxel_or_a_steep_gradient() -> None:
    R, Rs, ids, cls, Rg, e, t, sym = gate_inputs()
    cls_sharp = {0: F.SHARP, 1: F.GRADIENT_ONLY}
    bad = R.copy()
    bad[0] = F.rotation_matrices(np.array([[0.0, 0.0, np.radians(0.01)]]))[0] @ bad[0]
    assert not SCH.check(bad, Rs, ids, cls_sharp, Rg, e, t, sym)["gate"]["checks"][
        "sharp_voxels_equal_R_g"
    ]
    wide = R.copy()
    k = int(np.nonzero(ids == 0)[0][0])
    wide[k] = F.rotation_matrices(np.array([[0.0, 0.0, np.radians(1.5)]]))[0] @ Rs[k]
    chk = SCH.check(wide, Rs, ids, cls, Rg, e, t, sym)["gate"]["checks"]
    assert not chk["noise_within_1_degree_of_smooth_field"]
    steep = R.copy()
    j = int(np.nonzero(ids == 1)[0][0])
    steep[j] = F.rotation_matrices(np.array([[0.0, 0.0, np.radians(2.0)]]))[0] @ steep[j]
    assert not SCH.check(steep, Rs, ids, cls, Rg, e, t, sym)["gate"]["checks"][
        "gradient_only_edge_pairs_at_most_0.5"
    ]


def test_boundary_minima_find_the_least_misorientation_across_a_boundary() -> None:
    c, ids, Rg, _ = two_grain_patch()
    t = F.touching_pairs(c, 1.0)
    from icenine.ray_indexing import proper_rotations
    from icenine.symmetry import create_cubic_symmetry

    sym = proper_rotations(create_cubic_symmetry(3.6))
    R = np.stack([Rg[int(g)] for g in ids])
    best = SCH.boundary_minima(ids, R, t, sym)
    assert list(best) == [(0, 1)]
    from icenine.ray_indexing import misorientation_deg

    assert best[(0, 1)] == pytest.approx(float(misorientation_deg(Rg[0], Rg[1], sym)))


# ------------------------------------------------------------------------------- configs


def test_cpp_config_sets_sample_radius_2_and_the_fine_side_and_drops_python_lines() -> None:
    text = CFG.cpp_config(Path("/x/grid.mic"), FRAMES, "outc")
    k = CFG.keys_of(text)
    assert not ({"BFSRevisitRefit", "BFSReseed", "FusedCost"} & set(k))
    assert k["SampleRadius"] == "2.0"
    assert k["MinSideLength"] == k["MaxInitSideLength"] == "0.0046875"
    assert k["SampleFilename"] == "/x/grid.mic" and k["InfileBasename"] == FRAMES
    ps = CFG.keys_of(CFG.BASE_PS.read_text())  # every search line is the stock base's
    for key in ("MaxMCSteps", "MinLocalResolution", "MaxQ", "SuccessiveRestarts", "EtaLimit"):
        assert k[key] == ps[key]
    assert ps["SampleRadius"] == "0.50"  # what fine1 ran with


def test_rcp_config_adds_the_sixteen_keys_and_sets_sample_radius() -> None:
    text = CFG.python_config("RCP", Path("/x/grid.mic"), FRAMES, "out")
    added = {s.split()[0] for s in CFG.RCP_KEYS}
    assert len(CFG.RCP_KEYS) == 16 and len(added) == 16
    assert CFG.changed_keys(CFG.BASE_J5S.read_text(), text) == sorted(set(CFG.CHANGED_KEYS) | added)
    assert CFG.J5S_KEYS == CFG.lib_check_keys()  # step 21's key list, same order
    k = CFG.keys_of(text)
    assert k["SeedIndexMCDefer"] == "1" and k["RayIndexMaxPairs"] == "50000"
    assert k["RayIndexFallbackFull"] == "1" and k["SeedReuseKnown"] == "1"
    assert k["SampleRadius"] == "2.0"


def test_rcp_config_parses_and_check_bfs_keys_accepts_all_sixteen(tmp_path: Path) -> None:
    from icenine.config_file import ConfigFile
    from icenine.reconstructor import check_bfs_keys

    p = tmp_path / "rcp.config"
    p.write_text(CFG.python_config("RCP", Path("/x/grid.mic"), FRAMES, "out"))
    cfg = ConfigFile.from_file(str(p))
    assert cfg.seed_index_mc_defer == 1 and cfg.seed_reuse_known == 1
    assert cfg.ray_index_max_pairs == 50000 and cfg.ray_index_fallback_full == 1
    assert cfg.boundary_pass and cfg.boundary_pass_repair and cfg.boundary_pass_screen
    check_bfs_keys(cfg)  # raises on a violated precondition


def test_counts_summary_splits_global_searches_from_the_rest() -> None:
    st = dict(n_seeds=3, n_reseeds=2, n_evals_seed=300, n_evals_reseed=200, n_evals_neighbor=50)
    st.update(n_evals_boundary_repair=7, wall_seed_s=30.0, wall_reseed_s=20.0)
    s = RP.counts_summary(st, 100.0)
    assert s["evals_all"] == 557 and s["evals_global"] == 500 and s["global_searches"] == 5
    assert s["wall_other_s"] == 50.0


def test_run_script_dry_run_writes_the_config_and_starts_nothing(tmp_path: Path) -> None:
    r = subprocess.run(
        ["uv", "run", "python", str(ICE / "scripts/autoresearch/tex1/tex1_run_py.py")]
        + ["--arm", "RCP", "--grid", "slice", "--out", str(tmp_path), "--dry-run"],
        capture_output=True,
        text=True,
        cwd=str(ICE),
    )
    assert r.returncode == 0, r.stderr[-500:]
    assert (tmp_path / "arm.config").exists() and "dry run" in r.stdout
    assert not (tmp_path / "run.json").exists()


# ----------------------------------------------------------------------- freeze and launcher


def test_freeze_check_refuses_a_missing_or_uncommitted_file(tmp_path: Path, monkeypatch) -> None:
    monkeypatch.setattr(PP, "FROZEN", tmp_path / "frozen.json")
    assert "does not exist" in FZ.problems(False, False)
    (tmp_path / "frozen.json").write_text("{}")
    assert "not committed" in FZ.problems(False, False)


def test_launcher_syntax_defaults_and_argument_guards() -> None:
    assert subprocess.run(["bash", "-n", str(LAUNCH)]).returncode == 0
    for name in ("tex1_launch_slice.sh", "tex1_launch_whole.sh"):
        assert subprocess.run(["bash", "-n", str(LAUNCH.parent / name)]).returncode == 0
    text = LAUNCH.read_text()
    assert "--continue-contended" in text and "exit 5" in text
    assert "ARMS=RCP || ARMS=CSP,RCP" in text
    bad = subprocess.run(["bash", str(LAUNCH), "--stage", "mid"], capture_output=True)
    assert bad.returncode == 64
    nocap = subprocess.run(["bash", str(LAUNCH), "--stage", "whole"], capture_output=True)
    assert nocap.returncode == 64 and b"cap-h" in nocap.stderr
    unk = subprocess.run(["bash", str(LAUNCH), "--stage", "slice", "--bogus"], capture_output=True)
    assert unk.returncode == 64


# ------------------------------------------------------------------------------- scoring


def test_score_matches_output_voxels_to_grid_rows_by_position_and_direction() -> None:
    g_xy = np.array([[0.0, 0.0], [0.0, 0.0001], [0.5, 0.0]])
    g_dir = np.array([1, 2, 1])
    idx = SC.match(np.array([[0.5, 0.0], [0.0, 0.0001]]), np.array([1, 2]), g_xy, g_dir)
    assert idx.tolist() == [2, 1]
    with pytest.raises(AssertionError):
        SC.match(np.array([[0.2, 0.0]]), np.array([1]), g_xy, g_dir)
    with pytest.raises(AssertionError):
        SC.match(np.array([[0.0, 0.0], [0.0, 0.0]]), np.array([1, 1]), g_xy, g_dir)


def test_class_row_counts_wrong_no_output_and_closer_than_the_trivial_answer() -> None:
    err = np.array([0.1, 0.5, 2.0, np.inf, 0.2])
    has = np.array([True, True, True, False, True])
    triv = np.array([0.4, 0.3, 0.9, 0.5, 0.2])
    r = SC.class_row(err, has, triv, np.ones(5, dtype=bool))
    assert (r["n"], r["n_wrong"], r["n_no_output"]) == (5, 2, 1)  # 2.0 wrong, no output wrong
    assert r["n_closer_than_trivial"] == 1  # only 0.1 < 0.4 (0.5 > 0.3, 0.2 is not < 0.2)
    assert r["closer_fraction"] == pytest.approx(0.25)
    assert r["trivial_wrong"] == 0 and r["median_err_deg"] == pytest.approx(0.35)
    sub = SC.class_row(err, has, triv, np.array([1, 1, 0, 0, 0], dtype=bool))
    assert sub["n"] == 2 and sub["n_wrong"] == 0


def test_class_rows_split_by_class() -> None:
    err = np.array([0.1, 3.0, 0.2, 0.3])
    rows = SC.class_rows(
        err, np.ones(4, dtype=bool), np.zeros(4), np.array([0, 1, 2, 3]), np.ones(4, dtype=bool)
    )
    assert rows["all"]["n_wrong"] == 1 and rows["noise"]["n_wrong"] == 1
    assert rows["sharp"]["n"] == rows["gradient"]["n"] == rows["both"]["n"] == 1


def test_pair_counts_by_voxel() -> None:
    a = np.array([True, True, False, False])
    b = np.array([True, False, True, False])
    assert SC.pair_counts(a, b) == dict(both_wrong=1, only_a_wrong=1, only_b_wrong=1, neither=1)


def test_spread_rows_compare_recovered_with_true_spread() -> None:
    from icenine.ray_indexing import proper_rotations
    from icenine.symmetry import create_cubic_symmetry

    sym = proper_rotations(create_cubic_symmetry(3.6))
    rv = np.zeros((4, 3))
    rv[:, 2] = np.radians([0.0, 1.0, 2.0, 4.0])
    truth = F.rotation_matrices(rv)
    rec = F.rotation_matrices(rv * 0.5)  # recovered spread is half the true one
    rows = SC.spread_rows(
        rec, np.ones(4, dtype=bool), np.zeros(4, dtype=int), {0: np.eye(3)}, truth, [0], sym
    )
    assert rows[0]["true_max_deg"] == pytest.approx(4.0) and rows[0][
        "rec_max_deg"
    ] == pytest.approx(2.0)
    assert SC.spread_summary(rows)["ratio_p90_median"] == pytest.approx(0.5)


def test_cost_evaluation_lines_sum_over_the_clients() -> None:
    text = "noise\n[cost evaluations] client 1: 5\n[cost evaluations] client 2: 7\nmore\n"
    assert SC.eval_lines_from_log(text) == dict(n_lines=2, total=12)
    assert SC.eval_lines_from_log("a stock build prints nothing\n") == dict(n_lines=0, total=0)


def test_score_refuses_to_run_before_both_end_markers_exist(tmp_path: Path, monkeypatch) -> None:
    monkeypatch.setattr(SC, "RUNS", tmp_path)
    monkeypatch.setattr(sys, "argv", ["tex1_score.py"])
    with pytest.raises(SystemExit) as e:
        SC.main()
    assert "missing" in str(e.value)


# ------------------------------------------------------------------------------- plots


def tiny_mic_files(tmp_path: Path) -> Any:
    from icenine.mic_file import MicFile

    src = MicFile.read(
        str(
            ICE.parent
            / "Examples/Example2.ManyGrains/SimInput/rand_500grains_1mm_neworient_s0_grid.mic"
        )
    )
    src.voxels = src.voxels[:60]
    truth = tmp_path / "truth.mic"
    src.write(str(truth))
    rng = np.random.default_rng(0)
    for k, v in enumerate(src.voxels):
        v.orientation = F.rotation_matrices(rng.normal(size=(1, 3)) * 0.01)[0] @ np.asarray(
            v.orientation
        )
    arm = tmp_path / "arm.mic"
    src.write(str(arm))
    gm = tmp_path / "grain.npz"
    np.savez(gm, grain=(np.arange(60) // 20).astype(np.int32))
    return truth, arm, gm


def test_plot_mic_without_the_new_flags_is_byte_equal(tmp_path: Path) -> None:
    truth, arm, gm = tiny_mic_files(tmp_path)
    a, b = tmp_path / "a.png", tmp_path / "b.png"
    sa = PM.plot_maps(str(truth), [str(arm)], ["x"], str(a))
    sb = PM.plot_maps(str(truth), [str(arm)], ["x"], str(b), boundaries=None, mis_max=None)
    assert a.read_bytes() == b.read_bytes() and sa == sb
    c = tmp_path / "c.png"
    PM.plot_maps(str(truth), [str(arm)], ["x"], str(c), boundaries=str(gm), mis_max=2.0)
    assert c.read_bytes() != a.read_bytes()


def test_boundary_segments_are_the_shared_edges_between_different_grains() -> None:
    c = patch(3)
    # four triangles: 0, 1 share an edge (up and its right down); give them two grains
    tri = np.array(
        [
            [[0, 0], [1, 0], [0.5, H]],  # up
            [[0.5, H], [1, 0], [1.5, H]],  # down on its right edge
            [[1, 0], [2, 0], [1.5, H]],  # up, shares the edge (1,0)-(1.5,H) with the down one
        ],
        dtype=float,
    )
    assert len(PM.boundary_segments(tri, np.array([0, 0, 0]))) == 0
    s = PM.boundary_segments(tri, np.array([0, 1, 1]))
    assert len(s) == 1 and sorted(map(tuple, s[0].tolist())) == [(0.5, H), (1.0, 0.0)]
    assert len(PM.boundary_segments(tri, np.array([0, 1, 2]))) == 2
    assert c.shape[1] == 2


def test_histogram_fractions_keep_every_voxel() -> None:
    e = np.array([0.0, 1e-5, 0.05, 1.0, 62.0, 90.0])
    f = PL.histogram_fractions(e)
    assert (
        f.sum() == pytest.approx(1.0)
        and f[0] == pytest.approx(2 / 6)
        and f[-1] == pytest.approx(2 / 6)
    )
    assert PL.histogram_fractions(np.array([])).sum() == 0.0


def synthetic_score() -> Dict[str, Any]:
    def cls(n: int, w: int) -> Dict[str, Any]:
        return dict(n=n, n_wrong=w, wilson95=[max(w - 3, 0) / n, (w + 3) / n])

    arms = {}
    for arm, w, ev, wall, wk in (("CSP", 400, 2.0e9, 20000.0, 12), ("RCP", 30, 5.0e7, 9000.0, 1)):
        arms[arm] = dict(
            workers=wk,
            evals_all=ev,
            wall_s=wall,
            common_n=98000,
            common_classes={"all": cls(98000, w)},
            classes={"all": cls(98280, w + 5)},
        )
    return dict(arms=arms)


def test_chart_rows_and_the_two_charts_on_synthetic_scores(tmp_path: Path) -> None:
    panel = PL.chart_rows(synthetic_score())
    assert panel["id"] == "tex24_98280" and [r["arm"] for r in panel["rows"]] == ["CSP", "RCP"]
    assert panel["rows"][0]["workers"] == 12 and "13 ranks" in panel["rows"][0]["label"]
    grey = [
        dict(arm="CSP", evaluations=1.4e9, wrong_on_common_85676=4400, wall_s=10000.0, workers=12)
    ]
    PL.cost_accuracy_figure(panel, grey, tmp_path / "a")
    PL.wall_time_figure(panel, grey, tmp_path / "b")
    for stem in ("a", "b"):
        assert (tmp_path / f"{stem}.png").stat().st_size > 5000
        assert (tmp_path / f"{stem}.svg").stat().st_size > 1000
    json.dumps(panel)  # serialisable


def test_distribution_figure_on_synthetic_errors(tmp_path: Path) -> None:
    rng = np.random.default_rng(0)
    n = 400
    err = 10 ** rng.uniform(-3, 1, n)
    has = rng.random(n) > 0.05
    err[~has] = np.inf
    PL.distribution_figure(
        err, has, np.abs(rng.normal(0, 0.5, n)), rng.integers(0, 4, n), "RCP", tmp_path / "d"
    )
    assert (tmp_path / "d.png").exists() and (tmp_path / "d.svg").exists()
