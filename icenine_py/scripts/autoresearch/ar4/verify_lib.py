"""AR4 step 1: pure functions of the learned verifier (no cost function, no sample).

A *pair* is one voxel with one orientation; it is *right* when the orientation is within 1 degree
of the voxel's truth. A score is in [0, 1]; a pair is *flagged* when its score is at or under a
cut. Nothing here is specific to a crystal system: the peak family is a rank by the length of the
reciprocal vector, not Miller indices. Tests: ``tests/test_ar4_verifier.py``.
"""

from typing import Dict, List, Sequence, Tuple

import numpy as np

RIGHT_DEG = 1.0  # a pair is right when its error is at most this (degrees)
# the summary numbers of a pair, in this order (the first three come from one evaluation of the
# library's cost function; the others from the per-peak shares and the stored ray lists)
SUMMARY_NAMES = (
    "hit",
    "cost",
    "conf",
    "lit_share",
    "low_quarter",
    "n_explained",
    "best_votes",
)


def right_label(err_deg: np.ndarray) -> np.ndarray:
    """True where the error against the truth is at most 1 degree."""
    return np.asarray(err_deg) <= RIGHT_DEG


def family_rank(g_mag: np.ndarray, decimals: int = 4) -> np.ndarray:
    """Dense 0-based rank of each reciprocal-vector length among the lengths in ``g_mag``
    (lengths equal to ``decimals`` places share a rank)."""
    r = np.round(np.asarray(g_mag, dtype=np.float64), decimals)
    _, inv = np.unique(r, return_inverse=True)
    return inv.reshape(r.shape).astype(np.int64)


def family_rank_in(g_mag: np.ndarray, ladder: np.ndarray, decimals: int = 4) -> np.ndarray:
    """Rank of each length in a fixed ascending ``ladder`` of lengths (the phase's families)."""
    lad = np.unique(np.round(np.asarray(ladder, dtype=np.float64), decimals))
    return np.searchsorted(lad, np.round(np.asarray(g_mag, dtype=np.float64), decimals)).astype(
        np.int64
    )


def lit_share(share: np.ndarray) -> float:
    """Share of the predicted peaks with any lit pixel (0.0 for no peaks)."""
    share = np.asarray(share)
    return float((share > 0).mean()) if share.size else 0.0


def low_quarter(share: np.ndarray) -> float:
    """Mean of the lowest quarter (rounded up) of the per-peak shares (0.0 for no peaks)."""
    share = np.sort(np.asarray(share, dtype=np.float64))
    if share.size == 0:
        return 0.0
    k = int(np.ceil(share.size / 4))
    return float(share[:k].mean())


def summary_vector(
    hit: float,
    cost: float,
    conf: float,
    share: np.ndarray,
    n_explained: float,
    best_votes: float,
) -> np.ndarray:
    """The summary numbers of a pair in the order of ``SUMMARY_NAMES``."""
    return np.array(
        [hit, cost, conf, lit_share(share), low_quarter(share), n_explained, best_votes],
        dtype=np.float64,
    )


def pool(parts: Sequence[np.ndarray]) -> np.ndarray:
    """Concatenate the arrays of several samples."""
    return np.concatenate([np.asarray(p) for p in parts])


def pick_threshold(score: np.ndarray, wrong: np.ndarray, frac: float = 0.95) -> float:
    """The cut that flags the fewest voxels while flagging at least ``frac`` of the wrong ones.

    A voxel is flagged when ``score <= cut``. The cut is the ``ceil(frac * n_wrong)``-th smallest
    score among the wrong voxels; voxels tied with the cut are all flagged (tie rule)."""
    score = np.asarray(score, dtype=np.float64)
    wrong = np.asarray(wrong, dtype=bool)
    n_wrong = int(wrong.sum())
    if n_wrong == 0:
        raise ValueError("no wrong voxel to flag")
    need = int(np.ceil(frac * n_wrong - 1e-9))
    return float(np.sort(score[wrong])[need - 1])


def flag_counts(score: np.ndarray, wrong: np.ndarray, cut: float) -> Tuple[int, int]:
    """(voxels flagged, wrong voxels flagged) for ``score <= cut``."""
    f = np.asarray(score, dtype=np.float64) <= cut
    return int(f.sum()), int((f & np.asarray(wrong, dtype=bool)).sum())


def auc(score: np.ndarray, right: np.ndarray) -> float:
    """Share of (right, wrong) pairs in which the right one scores higher (a tie counts one half);
    NaN when there is no right or no wrong voxel. Counts ordered pairs; not a test."""
    score = np.asarray(score, dtype=np.float64)
    right = np.asarray(right, dtype=bool)
    a, b = score[right], score[~right]
    if len(a) == 0 or len(b) == 0:
        return float("nan")
    b_sorted = np.sort(b)
    lo = np.searchsorted(b_sorted, a, side="left")  # wrong scores strictly under each right one
    hi = np.searchsorted(b_sorted, a, side="right")  # wrong scores at or under
    return float(((lo + hi) / 2.0).sum() / (len(a) * len(b)))


def curve_at_wrong_scores(
    score: np.ndarray, wrong: np.ndarray
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Measured points of the flag curve: for each distinct score of a wrong voxel as the cut,
    (cuts, voxels flagged, wrong voxels flagged), in rising cut order."""
    score = np.asarray(score, dtype=np.float64)
    wrong = np.asarray(wrong, dtype=bool)
    cuts = np.unique(score[wrong])
    flagged = np.array([int((score <= c).sum()) for c in cuts], dtype=np.int64)
    wf = np.array([int(((score <= c) & wrong).sum()) for c in cuts], dtype=np.int64)
    return cuts, flagged, wf


def random_axes(n: int, seed: int = 0) -> np.ndarray:
    """``n`` random unit axes (normal vectors, seeded)."""
    a = np.random.default_rng(seed).normal(size=(n, 3))
    return a / np.linalg.norm(a, axis=1, keepdims=True)


def rotate_about_axis(R: np.ndarray, deg: float, axis: np.ndarray) -> np.ndarray:
    """``R`` turned by ``deg`` degrees about ``axis`` in the sample frame (Rodrigues, on the left).
    The turn is a rotation of the orientation by exactly ``deg`` (misorientation ``deg``)."""
    k = np.asarray(axis, dtype=np.float64)
    k = k / np.linalg.norm(k)
    t = np.radians(deg)
    K = np.array([[0.0, -k[2], k[1]], [k[2], 0.0, -k[0]], [-k[1], k[0], 0.0]])
    Q = np.eye(3) + np.sin(t) * K + (1.0 - np.cos(t)) * (K @ K)
    return Q @ np.asarray(R, dtype=np.float64)


def pad_ragged(
    flat: np.ndarray, offsets: np.ndarray, max_len: int
) -> Tuple[np.ndarray, np.ndarray]:
    """Pad per-pair lists stored flat (``flat[offsets[i]:offsets[i+1]]`` is pair ``i``) to
    ``(n_pairs, max_len, n_features)`` plus a boolean mask; longer lists are cut to ``max_len``."""
    n = len(offsets) - 1
    nf = flat.shape[1]
    x = np.zeros((n, max_len, nf), dtype=flat.dtype)
    m = np.zeros((n, max_len), dtype=bool)
    for i in range(n):
        seg = flat[offsets[i] : offsets[i + 1]][:max_len]
        x[i, : len(seg)] = seg
        m[i, : len(seg)] = True
    return x, m


def fit_logistic(
    x: np.ndarray, y: np.ndarray, epochs: int = 500, lr: float = 0.1, l2: float = 1e-4
) -> Tuple[np.ndarray, float, np.ndarray, np.ndarray]:
    """Logistic regression by full-batch gradient descent on standardized inputs, class-balanced
    (each class carries half the weight). Returns (weights, bias, mean, std) of the inputs."""
    x = np.asarray(x, dtype=np.float64)
    y = np.asarray(y, dtype=np.float64)
    mu = x.mean(axis=0)
    sd = x.std(axis=0)
    sd[sd < 1e-12] = 1.0
    z = (x - mu) / sd
    pos = max(float(y.sum()), 1.0)
    neg = max(float(len(y) - y.sum()), 1.0)
    wt = np.where(y > 0.5, 0.5 / pos, 0.5 / neg)
    w = np.zeros(z.shape[1])
    b = 0.0
    for _ in range(epochs):
        p = 1.0 / (1.0 + np.exp(-(z @ w + b)))
        r = (p - y) * wt
        w -= lr * (z.T @ r + l2 * w)
        b -= lr * r.sum()
    return w, float(b), mu, sd


def predict_logistic(
    x: np.ndarray, w: np.ndarray, b: float, mu: np.ndarray, sd: np.ndarray
) -> np.ndarray:
    """Probability of right from ``fit_logistic``'s parameters."""
    z = (np.asarray(x, dtype=np.float64) - mu) / sd
    return 1.0 / (1.0 + np.exp(-(z @ w + b)))


def wilson_interval(k: int, n: int) -> List[float]:
    """Wilson 95% interval as a list (the shared ``stats.wilson``)."""
    import sys
    from pathlib import Path

    sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "common"))
    import stats

    lo, hi = stats.wilson(k, n)
    return [float(lo), float(hi)]


# ---- use 2 and the sum-back check (added in the review fixes; used by ``verify_summary.py`` only)


def grid_only_evals(grid_evals: int, reseed: bool, seed_step_evals: int) -> int:
    """Grid-search evaluations of one grid-fallback case.

    AR5 step 1's ``part2_cases.json`` field ``grid_evals`` is, for a reseed case, the stored
    ``full_search_evaluations``: the whole fit from the seed, so the seed step plus the grid
    search. The seed step runs before any accept decision and is never skipped, so it is
    subtracted here. For a non-reseed case the field is already the grid search alone."""
    return int(grid_evals) - int(seed_step_evals) if reseed else int(grid_evals)


def use2_counts(
    fb_score: np.ndarray,
    fb_right: np.ndarray,
    fb_grid: np.ndarray,
    ps_score: np.ndarray,
    ps_right: np.ndarray,
    cut: float,
) -> Dict[str, int]:
    """Accept and refuse counts of one score at one cut.

    ``fb_*``: the ends of the seed steps that failed the seed test (the grid-fallback cases) and
    their grid-search evaluations; ``ps_*``: the ends of the seed steps that passed. A fallback
    end is accepted when its score is at or above the cut; a passing end is refused when its
    score is under the cut."""
    fs, fr = np.asarray(fb_score, np.float64), np.asarray(fb_right, bool)
    fg = np.asarray(fb_grid, np.int64)
    ps, pr = np.asarray(ps_score, np.float64), np.asarray(ps_right, bool)
    acc = fs >= cut
    ref = ps < cut
    return dict(
        right_accepted=int((acc & fr).sum()),
        wrong_accepted=int((acc & ~fr).sum()),
        grid_evals_of_right_accepted=int(fg[acc & fr].sum()),
        passed_refused=int(ref.sum()),
        passed_refused_right=int((ref & pr).sum()),
        passed_refused_wrong=int((ref & ~pr).sum()),
    )


def use2_rule_costs(counts: Dict[str, int], mean_grid_evals: float) -> Dict[str, float]:
    """The two readings of use 2, priced by arithmetic from stored numbers (not a run).

    OR rule: a seed is accepted when the seed test passes or the score is at or above the cut;
    nothing is refused, so the skipped grid searches are the saving. Replace rule: the score
    replaces the seed test; each refused passing seed then pays one grid search, priced at
    ``mean_grid_evals`` (the mean of the stored grid-fallback cases)."""
    saved = float(counts["grid_evals_of_right_accepted"])
    extra_right = counts["passed_refused_right"] * mean_grid_evals
    extra_wrong = counts["passed_refused_wrong"] * mean_grid_evals
    return dict(
        or_saved=saved,
        or_refused=0.0,
        replace_saved=saved,
        replace_extra_right=float(extra_right),
        replace_extra_wrong=float(extra_wrong),
        replace_net_right_only=float(saved - extra_right),
        replace_net_all=float(saved - extra_right - extra_wrong),
    )


def sumback_differs(peaks: np.ndarray, offsets: np.ndarray, counts: np.ndarray) -> np.ndarray:
    """Per pair: do the per-peak pixel counts sum to other totals than ``counts``?

    ``peaks`` holds rows ``(rank, dist, omega, share, pixels)`` per peak, pair ``i`` owning rows
    ``offsets[i]:offsets[i + 1]``; ``counts`` holds ``(pixel_overlap, pixel_on_detector, ...)``
    of the fused evaluation. The overlap of a peak is ``share * pixels`` rounded."""
    p = np.asarray(peaks, np.float64)
    off = np.asarray(offsets, np.int64)
    c = np.asarray(counts)
    cs_pod = np.concatenate([[0.0], np.cumsum(p[:, 4])])
    cs_po = np.concatenate([[0.0], np.cumsum(np.rint(p[:, 3] * p[:, 4]))])
    pod = np.rint(cs_pod[off[1:]] - cs_pod[off[:-1]])
    po = np.rint(cs_po[off[1:]] - cs_po[off[:-1]])
    out: np.ndarray = (pod != c[:, 1]) | (po != c[:, 0])
    return out
