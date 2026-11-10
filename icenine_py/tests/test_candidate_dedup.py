"""Task C (evaluation reduction): fast symmetry-reduced clustering against the O(N^2 x 24)
reference ``csl.dedup_orientations``, and the CandidateDedup specification and stage functions."""

import sys
from pathlib import Path
from typing import List, Set

import numpy as np
import pytest
from scipy.spatial.transform import Rotation

sys.path.insert(0, str(Path(__file__).parent.parent / "scripts" / "findoptimal_robustness"))

import csl  # noqa: E402

from icenine import candidate_dedup as CD  # noqa: E402
from icenine.sampling import matrix_to_quaternion  # noqa: E402

CUBIC_Q = np.array([matrix_to_quaternion(m) for m in csl.OPS])
RNG = np.random.default_rng(7)


def _random_set(n: int, n_clusters: int, spread_deg: float, rng: np.random.Generator) -> np.ndarray:
    """Clusters of orientations, each member re-expressed by a random cubic operator."""
    centres = Rotation.random(n_clusters, random_state=int(rng.integers(1 << 30))).as_matrix()
    which = rng.integers(0, n_clusters, n)
    jitter = Rotation.from_rotvec(
        rng.normal(size=(n, 3)) * np.radians(spread_deg) / np.sqrt(3)
    ).as_matrix()
    ops = csl.OPS[rng.integers(0, 24, n)]
    out: np.ndarray = centres[which] @ jitter @ ops
    return out


def _brute_pairs(R: np.ndarray, tol_deg: float) -> List[Set[int]]:
    n = len(R)
    out: List[Set[int]] = [set() for _ in range(n)]
    for i in range(n):
        d = csl.reduced_misorientation_deg(R[i], R)
        for j in np.nonzero(d <= tol_deg)[0]:
            if j != i:
                out[i].add(int(j))
    return out


@pytest.mark.parametrize("tol_deg", [0.25, 1.0, 3.0, 8.0])
def test_greedy_nms_equals_reference_random(tol_deg: float) -> None:
    R = _random_set(260, 12, 4.0, np.random.default_rng(1))
    kept, leader = CD.greedy_nms(R, np.radians(tol_deg), CUBIC_Q)
    ref = csl.dedup_orientations(R, tol_deg)
    assert kept.tolist() == ref.tolist()
    assert set(leader.tolist()) == set(kept.tolist())
    assert all(leader[k] == k for k in kept)


def test_neighbour_lists_equal_brute_force() -> None:
    R = _random_set(150, 8, 3.0, np.random.default_rng(2))
    for tol_deg in (0.5, 2.0, 6.0):
        got = CD.neighbour_lists(R, np.radians(tol_deg), CUBIC_Q)
        want = _brute_pairs(R, tol_deg)
        assert [set(g.tolist()) for g in got] == want


def test_fz_boundary_sets() -> None:
    """Orientations on both sides of the cubic fundamental-zone boundary (rotation angle 45 deg
    about <001>, 60 deg about <111>, ... are class-ambiguous), each given as two equivalent forms,
    and small jitters across the boundary."""
    rng = np.random.default_rng(3)
    axes = np.array([[0, 0, 1], [1, 1, 1], [1, 1, 0], [0, 1, 0]], dtype=float)
    axes /= np.linalg.norm(axes, axis=1)[:, None]
    angles = np.array([45.0, 60.0, 62.8, 90.0, 44.0, 46.0])
    base = []
    for ax in axes:
        for a in angles:
            base.append(Rotation.from_rotvec(ax * np.radians(a)).as_matrix())
    base = np.stack(base)
    rows = []
    for b in base:
        for _ in range(4):
            jit = Rotation.from_rotvec(rng.normal(size=3) * np.radians(0.3)).as_matrix()
            rows.append(b @ jit @ csl.OPS[rng.integers(0, 24)])
    R = np.stack(rows)
    for tol_deg in (0.2, 1.0, 2.5):
        got = CD.neighbour_lists(R, np.radians(tol_deg), CUBIC_Q)
        assert [set(g.tolist()) for g in got] == _brute_pairs(R, tol_deg)
        kept, _ = CD.greedy_nms(R, np.radians(tol_deg), CUBIC_Q)
        assert kept.tolist() == csl.dedup_orientations(R, tol_deg).tolist()


def test_priority_order_is_respected() -> None:
    R = _random_set(120, 6, 2.0, np.random.default_rng(4))
    order = np.random.default_rng(5).permutation(len(R))
    kept, _ = CD.greedy_nms(R, np.radians(2.0), CUBIC_Q, order=order)
    ref = csl.dedup_orientations(R[order], 2.0)
    assert kept.tolist() == order[ref].tolist()


def test_other_symmetries_and_no_symmetry() -> None:
    """A 4-fold group about z and no symmetry: neighbours equal a brute-force quaternion min."""
    rng = np.random.default_rng(6)
    R = Rotation.random(120, random_state=11).as_matrix()
    R = np.concatenate([R, R[:30] @ Rotation.from_rotvec([0, 0, np.pi / 2]).as_matrix()])
    c4 = np.array([[np.cos(k * np.pi / 4), 0, 0, np.sin(k * np.pi / 4)] for k in range(4)])
    for sym in (c4, None):
        ops = (
            np.array([np.eye(3)])
            if sym is None
            else np.stack([Rotation.from_quat(np.r_[q[1:], q[0]]).as_matrix() for q in sym])
        )
        for tol_deg in (1.0, 20.0, 60.0):
            got = CD.neighbour_lists(R, np.radians(tol_deg), sym)
            for i in range(len(R)):
                M = (R[i][None] @ ops)[None] @ np.swapaxes(R, -1, -2)[:, None]
                ang = np.degrees(Rotation.from_matrix(M.reshape(-1, 3, 3)).magnitude())
                want = set(np.nonzero(ang.reshape(len(R), -1).min(axis=1) <= tol_deg)[0]) - {i}
                assert set(got[i].tolist()) == want


def test_speed_is_not_quadratic() -> None:
    """A guard against a Python-level pairwise loop, not a timing measurement: the bound is loose
    so that a loaded machine does not fail it."""
    import time

    R = _random_set(4000, 300, 1.5, np.random.default_rng(8))
    t0 = time.perf_counter()
    CD.greedy_nms(R, np.radians(1.0), CUBIC_Q)
    assert time.perf_counter() - t0 < 30.0


def test_components_and_merged_cliques() -> None:
    # three tight pairs on a line of small rotations: A-B close, B-C close, D far
    ax = np.array([0.0, 0.0, 1.0])
    angs = np.array([0.0, 1.0, 2.0, 40.0])
    R = np.stack([Rotation.from_rotvec(ax * np.radians(a)).as_matrix() for a in angs])
    clique = np.array([0, 1, 2, 3])
    lab = CD.merged_clique_labels(R, clique, np.radians(1.2), CUBIC_Q)
    assert lab[0] == lab[1] == lab[2] != lab[3]
    lab = CD.merged_clique_labels(R, clique, np.radians(0.5), CUBIC_Q)
    assert len(set(lab.tolist())) == 4
    nb = CD.neighbour_lists(R, np.radians(1.2), CUBIC_Q)
    assert CD.component_labels(4, nb).tolist().count(CD.component_labels(4, nb)[0]) == 3


def test_spec_parse_roundtrip_and_errors() -> None:
    assert not CD.DedupSpec.parse("off").active
    assert not CD.DedupSpec.parse("").active
    s = CD.DedupSpec.parse("clique=merged:2@0.5,pre=0.25,post=0.5,quarter=filtered")
    assert (s.clique, s.clique_k, s.clique_link, s.pre, s.post, s.quarter) == (
        "merged",
        2,
        0.5,
        0.25,
        0.5,
        "filtered",
    )
    assert CD.DedupSpec.parse(s.to_text()) == s
    assert CD.DedupSpec.parse("clique=cpp").clique_k == 1
    for bad in (
        "clique=foo:1", "pre=-1", "quarter=x", "bogus=1", "pre", "pre=nan", "post=inf",
        "clique=merged:1@nan", "clique=merged:1@-1",
    ):  # fmt: skip
        with pytest.raises(ValueError):
            CD.DedupSpec.parse(bad)


def test_best_k_per_group() -> None:
    group = np.array([0, 0, 0, 1, 1, 2])
    score = np.array([0.5, 0.1, 0.3, 0.2, 0.2, 0.9])
    assert CD.best_k_per_group(group, score, 1).tolist() == [False, True, False, True, False, True]
    assert CD.best_k_per_group(group, score, 2).tolist() == [False, True, True, True, True, True]


def test_keep_before_quick_mc_clique_and_pre() -> None:
    ax = np.array([0.0, 0.0, 1.0])
    angs = np.array([0.0, 0.2, 5.0, 5.1, 30.0])
    R = np.stack([Rotation.from_rotvec(ax * np.radians(a)).as_matrix() for a in angs])
    score = np.array([0.3, 0.1, 0.5, 0.4, 0.2])
    clique = np.array([0, 0, 1, 1, 2])
    diameter = np.radians(1.0)
    # off: everything kept
    assert CD.keep_before_quick_mc(CD.DedupSpec(), R, score, clique, diameter, CUBIC_Q).all()
    s = CD.DedupSpec.parse("clique=cpp:1")
    assert CD.keep_before_quick_mc(s, R, score, clique, diameter, CUBIC_Q).tolist() == [
        False,
        True,
        False,
        True,
        True,
    ]
    s = CD.DedupSpec.parse("pre=0.5")  # tolerance 0.5 deg: the pairs (0,1) and (2,3) collapse
    assert CD.keep_before_quick_mc(s, R, score, clique, diameter, CUBIC_Q).tolist() == [
        False,
        True,
        False,
        True,
        True,
    ]
    s = CD.DedupSpec.parse("clique=merged:1@6")  # link radius 6 deg merges the first four
    assert CD.keep_before_quick_mc(s, R, score, clique, diameter, CUBIC_Q).tolist() == [
        False,
        True,
        False,
        False,
        True,
    ]
    # an unknown clique (-1) is its own group
    s = CD.DedupSpec.parse("clique=cpp:1")
    unk = np.array([-1, -1, -1, -1, -1])
    assert CD.keep_before_quick_mc(s, R, score, unk, diameter, CUBIC_Q).all()


def test_distinct_order() -> None:
    ax = np.array([0.0, 0.0, 1.0])
    angs = np.array([0.0, 0.2, 5.0, 0.3, 30.0])  # ranked best first
    R = np.stack([Rotation.from_rotvec(ax * np.radians(a)).as_matrix() for a in angs])
    s = CD.DedupSpec.parse("post=0.5")
    order, n_distinct = CD.distinct_order(s, R, np.radians(1.0), CUBIC_Q)
    assert n_distinct == 3
    assert order.tolist() == [0, 2, 4, 1, 3]
