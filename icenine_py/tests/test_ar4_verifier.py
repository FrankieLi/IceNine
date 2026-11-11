"""Pure functions of AR4 step 1 (``scripts/autoresearch/ar4/verify_lib.py``): the peak-family rank,
the summary numbers of a pair, the threshold pick, the flag counts, the rotation about an axis,
the padding of per-peak lists and the logistic fit. Written before the module (plan "AR4 step 1").
"""

import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "scripts" / "autoresearch" / "ar4"))

import verify_lib as VL  # noqa: E402


def test_family_rank_is_dense_rank_by_length() -> None:
    g = np.array([2.0, 1.0, 2.0, 3.5, 1.0])
    assert VL.family_rank(g).tolist() == [1, 0, 1, 2, 0]


def test_family_rank_merges_lengths_equal_to_rounding() -> None:
    g = np.array([1.73205081, 1.73205083, 2.0])
    assert VL.family_rank(g, decimals=4).tolist() == [0, 0, 1]


def test_family_rank_orders_by_global_ladder() -> None:
    # the rank is looked up in the ladder of the phase, not of the pair
    ladder = np.array([1.0, 2.0, 3.0, 4.0])
    assert VL.family_rank_in(np.array([3.0, 1.0]), ladder).tolist() == [2, 0]


def test_lit_share_counts_peaks_with_any_lit_pixel() -> None:
    share = np.array([0.0, 0.5, 1.0, 0.0])
    assert VL.lit_share(share) == 0.5
    assert VL.lit_share(np.zeros(0)) == 0.0


def test_low_quarter_is_mean_of_lowest_quarter() -> None:
    share = np.array([0.9, 0.0, 0.8, 0.1, 0.7, 0.6, 0.5, 1.0])  # 8 peaks: lowest 2
    assert np.isclose(VL.low_quarter(share), 0.05)
    # 5 peaks: ceil(5 / 4) = 2 lowest
    assert np.isclose(VL.low_quarter(np.array([0.4, 0.2, 0.6, 0.8, 1.0])), 0.3)
    assert VL.low_quarter(np.zeros(0)) == 0.0


def test_summary_vector_order_and_values() -> None:
    share = np.array([0.0, 1.0, 1.0, 0.5])
    v = VL.summary_vector(hit=0.7, cost=0.3, conf=0.6, share=share, n_explained=12, best_votes=4)
    assert v.shape == (len(VL.SUMMARY_NAMES),)
    d = dict(zip(VL.SUMMARY_NAMES, v.tolist()))
    assert d["hit"] == 0.7 and d["cost"] == 0.3 and d["conf"] == 0.6
    assert d["lit_share"] == 0.75
    assert np.isclose(d["low_quarter"], 0.0)
    assert d["n_explained"] == 12 and d["best_votes"] == 4


def test_pick_threshold_flags_fewest_for_95_percent_of_wrong() -> None:
    # 20 wrong voxels with scores 0.01 .. 0.20; 95% of 20 is 19 -> cut at the 19th smallest
    wrong_scores = np.arange(1, 21) / 100.0
    right_scores = np.array([0.195, 0.5, 0.9])  # one right voxel sits under the cut
    score = np.concatenate([wrong_scores, right_scores])
    wrong = np.concatenate([np.ones(20, bool), np.zeros(3, bool)])
    cut = VL.pick_threshold(score, wrong, 0.95)
    assert np.isclose(cut, 0.19)
    flagged, wf = VL.flag_counts(score, wrong, cut)
    assert wf == 19 and flagged == 19  # the right voxel at 0.195 is above the cut
    assert VL.pick_threshold(score, wrong, 1.0) == 0.20


def test_pick_threshold_ties_at_the_cut_are_all_flagged() -> None:
    score = np.array([0.5, 0.5, 0.5, 0.9])
    wrong = np.array([True, True, False, False])
    cut = VL.pick_threshold(score, wrong, 0.5)  # need 1 of 2 wrong: cut 0.5, tie flags all three
    assert cut == 0.5
    assert VL.flag_counts(score, wrong, cut) == (3, 2)


def test_pick_threshold_pooled_over_samples() -> None:
    s1, w1 = np.array([0.1, 0.9]), np.array([True, False])
    s2, w2 = np.array([0.3, 0.2, 0.8]), np.array([True, True, False])
    score, wrong = VL.pool([s1, s2]), VL.pool([w1, w2])
    assert score.shape == (5,) and wrong.dtype == bool
    assert VL.pick_threshold(score, wrong, 0.95) == 0.3  # need ceil(.95 * 3) = 3 wrong


def test_pick_threshold_without_wrong_voxels_raises() -> None:
    try:
        VL.pick_threshold(np.array([0.5]), np.array([False]), 0.95)
    except ValueError:
        return
    raise AssertionError("expected ValueError")


def test_flag_counts_is_score_at_or_under_cut() -> None:
    score = np.array([0.2, 0.4, 0.4, 0.6])
    wrong = np.array([True, True, False, False])
    assert VL.flag_counts(score, wrong, 0.4) == (3, 2)
    assert VL.flag_counts(score, wrong, 0.1) == (0, 0)


def test_rotate_about_axis_turns_by_the_angle() -> None:
    rng = np.random.default_rng(1)
    R = np.linalg.qr(rng.normal(size=(3, 3)))[0]
    if np.linalg.det(R) < 0:
        R[:, 0] *= -1
    axis = np.array([0.3, -0.5, 0.8])
    for deg in (0.1, 0.3, 0.5):
        Q = VL.rotate_about_axis(R, deg, axis)
        rel = Q @ R.T
        ang = np.degrees(np.arccos(np.clip((np.trace(rel) - 1) / 2, -1, 1)))
        assert np.isclose(ang, deg, atol=1e-6)
        assert np.allclose(Q @ Q.T, np.eye(3), atol=1e-12)


def test_random_axes_are_unit_and_seeded() -> None:
    a = VL.random_axes(5, seed=0)
    b = VL.random_axes(5, seed=0)
    assert np.allclose(np.linalg.norm(a, axis=1), 1.0) and np.array_equal(a, b)


def test_pad_ragged_pads_with_zeros_and_mask() -> None:
    flat = np.arange(10, dtype=np.float32).reshape(5, 2)  # 5 peaks, 2 features
    offsets = np.array([0, 2, 2, 5])  # pairs of 2, 0 and 3 peaks
    x, m = VL.pad_ragged(flat, offsets, max_len=4)
    assert x.shape == (3, 4, 2) and m.shape == (3, 4)
    assert m.sum(axis=1).tolist() == [2, 0, 3]
    assert np.array_equal(x[0, :2], flat[:2]) and np.all(x[0, 2:] == 0)
    assert np.array_equal(x[2, :3], flat[2:5])


def test_pad_ragged_truncates_to_max_len() -> None:
    flat = np.ones((6, 1), np.float32)
    x, m = VL.pad_ragged(flat, np.array([0, 6]), max_len=4)
    assert x.shape == (1, 4, 1) and m.sum() == 4


def test_logistic_fit_separates_a_clear_signal() -> None:
    rng = np.random.default_rng(0)
    x = rng.normal(size=(400, 3))
    y = (x[:, 0] + 0.1 * rng.normal(size=400) > 0).astype(np.float64)
    w, b, mu, sd = VL.fit_logistic(x, y, epochs=300, lr=0.1)
    p = VL.predict_logistic(x, w, b, mu, sd)
    assert ((p > 0.5) == (y > 0.5)).mean() > 0.95
    assert p.min() >= 0.0 and p.max() <= 1.0


def test_wrong_label_is_error_above_one_degree() -> None:
    err = np.array([0.0, 0.99, 1.0, 1.01, 30.0])
    assert VL.right_label(err).tolist() == [True, True, True, False, False]
