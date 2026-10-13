"""scripts/phase_d scoring helpers (score_bfs.py) and the pilot region selection."""

import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(ROOT / "scripts" / "phase_d"))

import grains as G  # noqa: E402
import make_region as MR  # noqa: E402
import score_bfs as S  # noqa: E402
from icenine.orientation_search import get_symmetry_quaternions  # noqa: E402
from icenine.symmetry import create_cubic_symmetry  # noqa: E402

SYM = get_symmetry_quaternions(create_cubic_symmetry(1.0))


def _q(*euler: float) -> np.ndarray:
    return G.euler_to_quat(np.array([euler], dtype=float))[0]


def test_pair_misorientation_matches_matrix_version_and_symmetry() -> None:
    qa = np.stack([_q(30, 40, 50), _q(10, 20, 30), _q(5, 5, 5)])
    qb = np.stack([_q(35, 40, 50), G._qmul(_q(10, 20, 30), SYM[7]), _q(5, 5, 5)])
    d = S.pair_misorientation_deg(qa, qb, SYM)
    ref = np.array(
        [G.misorientation_matrix_deg(qa[i : i + 1], qb[i : i + 1], SYM)[0, 0] for i in range(3)]
    )
    assert np.allclose(d, ref, atol=1e-6)
    assert d[0] == pytest.approx(5.0, abs=1e-3)
    assert d[1] < 1e-4 and d[2] < 1e-4


def test_wrong_summary_counts_and_median_of_right() -> None:
    err = np.array([0.01, 0.02, 0.03, 5.0, 30.0, 1.0])  # 1.0 is not wrong (> 1 is)
    s = S.wrong_summary(err)
    assert (s["n"], s["wrong"], s["n_right"]) == (6, 2, 4)
    assert s["median_right_err_deg"] == pytest.approx(0.025)
    assert 0 < s["wilson_lo"] < 2 / 6 < s["wilson_hi"] < 1
    empty = S.wrong_summary(np.array([]))
    assert empty["n"] == 0 and empty["wrong_rate"] is None


def test_clustered_paired_p_exact_and_degenerate() -> None:
    grain = np.repeat(np.arange(4), 5)
    a = np.zeros(20, dtype=bool)
    b = np.zeros(20, dtype=bool)
    assert S.clustered_paired_p(a, b, grain)["p"] == 1.0
    # all differences in one grain: only 2 sign patterns, p = 1.0 however many voxels differ
    b[:5] = True
    r = S.clustered_paired_p(a, b, grain)
    assert r["n_clusters_differing"] == 1 and r["p"] == 1.0 and r["wrong_b"] == 5
    # the same differences spread over four grains: 2 of 16 patterns reach the statistic
    b = np.zeros(20, dtype=bool)
    b[[0, 5, 10, 15]] = True
    assert S.clustered_paired_p(a, b, grain)["p"] == pytest.approx(2 / 16)


def test_cluster_bootstrap_wider_than_independent_for_clustered_errors() -> None:
    grain = np.repeat(np.arange(20), 10)
    wrong = np.zeros(200, dtype=bool)
    wrong[:50] = True  # 5 whole grains wrong
    ci = S.cluster_bootstrap_rate(wrong, grain, n_boot=500)
    lo, hi = S.wilson(50, 200)
    assert ci["cluster_boot_hi"] - ci["cluster_boot_lo"] > hi - lo


def test_grain_status_found_partial_lost_fragmented() -> None:
    qa, qb = _q(10, 20, 30), _q(60, 70, 80)
    grain = np.array([0] * 4 + [1] * 4 + [2] * 4 + [3] * 4)
    # grain 0 all right; grain 1: 1 of 4 right (partial) with two other orientations;
    # grain 2 all wrong with one common wrong orientation (lost, not fragmented);
    # grain 3: 3 right + 1 wrong (found, fragmented)
    q = np.stack([qa] * 4 + [qa, qb, qb, _q(100, 100, 100)] + [qb] * 4 + [qa, qa, qa, qb])
    err = np.array([0, 0, 0, 0, 0, 9, 9, 9, 9, 9, 9, 9, 0, 0, 0, 9], dtype=float)
    st = S.grain_status(err, grain, q, SYM)
    assert (st["found"], st["partial"], st["lost"]) == (2, 1, 1)
    assert st["lost_ids"] == [2]
    assert st["fragmented_ids"] == [1, 3]


def test_cross_boundary_wrong_flags_neighbour_orientation() -> None:
    pos = np.array([[0.0, 0.0], [0.01, 0.0], [0.02, 0.0]])
    grain = np.array([0, 0, 1])
    qa, qb, qc = _q(10, 20, 30), _q(60, 70, 80), _q(100, 20, 130)
    truth = np.stack([qa, qa, qb])
    rec = np.stack([qa, qb, qc])  # voxel 1 carries grain 1's orientation, voxel 2 a random one
    err = S.pair_misorientation_deg(rec, truth, SYM)
    flag = S.cross_boundary_wrong(err, grain, pos, rec, truth, SYM, radius=0.011)
    assert flag.tolist() == [False, True, False]


def test_select_disc_and_boundary_flags() -> None:
    cent = np.array([[x * 0.01, y * 0.01] for x in range(10) for y in range(10)], dtype=float)
    idx = MR.select_disc(cent, (0.045, 0.045), 4)
    assert idx.tolist() == [44, 45, 54, 55]  # file order
    assert len(MR.select_disc(cent, (0.045, 0.045), 9)) == 9
    grain = (cent[:, 0] >= 0.05).astype(int)  # two halves
    flag = MR.boundary_flags(cent, grain, side=0.01)
    assert flag.reshape(10, 10)[4:6].all() and not flag.reshape(10, 10)[:3].any()


def test_clustered_paired_p_exact_for_up_to_20_clusters_and_random_beyond() -> None:
    # 19 differing grains all in the same direction: exact p = 2 / 2^19
    grain = np.repeat(np.arange(19), 2)
    a = np.zeros(38, dtype=bool)
    b = np.zeros(38, dtype=bool)
    b[::2] = True
    r = S.clustered_paired_p(a, b, grain)
    assert r["exact"] and r["p"] == pytest.approx(2 / 2**19)
    # 24 differing grains: the Monte-Carlo branch, floor (1 + hits) / (1 + n_perm)
    grain = np.repeat(np.arange(24), 2)
    a = np.zeros(48, dtype=bool)
    b = np.zeros(48, dtype=bool)
    b[::2] = True
    r = S.clustered_paired_p(a, b, grain, n_perm=2000)
    assert "exact" not in r and r["p"] == pytest.approx(1 / 2001)
    # balanced differences give a large p in the random branch too
    b[::4] = False
    a[::4] = True
    assert S.clustered_paired_p(a, b, grain, n_perm=2000)["p"] > 0.5


def test_grain_status_reports_fragmented_among_found() -> None:
    qa, qb = _q(10, 20, 30), _q(60, 70, 80)
    grain = np.array([0] * 4 + [1] * 4)
    q = np.stack([qa, qa, qa, qb] + [qb] * 4)  # grain 0 found + fragmented; grain 1 lost
    err = np.array([0, 0, 0, 9, 9, 9, 9, 9], dtype=float)
    st = S.grain_status(err, grain, q, SYM)
    assert st["fragmented"] == 1 and st["fragmented_found"] == 1 and st["lost"] == 1
