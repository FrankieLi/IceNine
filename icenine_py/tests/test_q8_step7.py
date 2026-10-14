"""scripts/phase_d/q8_step7.py (Q8 step 7): the width rule, the shared random axes, the rotated
truths, the pairing groups and the bootstrap, on synthetic data and the tracked per-voxel results.
"""

import sys
from pathlib import Path

import numpy as np
import pytest
from scipy.spatial.transform import Rotation

ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(ROOT / "scripts" / "phase_d"))
sys.path.insert(0, str(ROOT / "scripts" / "findoptimal_robustness"))
sys.path.insert(0, str(ROOT / "scripts" / "common"))
sys.path.insert(0, str(ROOT / "scripts"))
sys.path.insert(0, str(ROOT / "benchmarks"))

import q8_step7 as Q  # noqa: E402


def test_width_of_grid_and_interpolated() -> None:
    radii = np.asarray(Q.RADII)  # 0.03 0.1 0.3 0.5 1 2 3
    curve = np.array([0.0, 0.05, 0.2, 0.4, 0.6, 0.9, 1.0, 1.0])  # index 0 = truth
    grid, interp = Q.width_of(curve, radii)
    assert grid == 0.5  # first radius whose value reaches halfway (0.5) is 0.5 deg (0.6)
    # linear in log r between 0.3 deg (0.4) and 0.5 deg (0.6): halfway is the midpoint
    assert interp == pytest.approx(np.exp(0.5 * (np.log(0.3) + np.log(0.5))))
    # a curve that does not rise above its truth value has no width
    assert np.isnan(Q.width_of(np.zeros(8), radii)[0])
    # halfway is relative to the truth value: a higher floor needs a higher level
    g2, _ = Q.width_of(np.array([0.5, 0.5, 0.55, 0.7, 0.8, 0.9, 1.0, 1.0]), radii)
    assert g2 == 0.5  # halfway between 0.5 and 1.0 is 0.75


def test_axes_depend_only_on_the_voxel() -> None:
    a1, a2 = Q.axes_of(123), Q.axes_of(123)
    assert a1.shape == (Q.N_AXES, 3) and np.array_equal(a1, a2)
    assert np.allclose(np.linalg.norm(a1, axis=1), 1.0)
    assert not np.array_equal(a1, Q.axes_of(124))


def test_rotated_truth_has_the_requested_angles() -> None:
    R = Rotation.random(random_state=5).as_matrix()
    Rs = Q.rotated_truth(R, 77)
    assert Rs.shape == (1 + len(Q.RADII) * Q.N_AXES, 3, 3)
    assert np.array_equal(Rs[0], R)
    ang = np.degrees(Rotation.from_matrix(Rs[1:] @ R.T).magnitude()).reshape(len(Q.RADII), -1)
    for i, r in enumerate(Q.RADII):
        assert np.allclose(ang[i], r, atol=1e-9)
    # the same axes on any image set: the rotated truths depend on the voxel only
    assert np.array_equal(Q.rotated_truth(R, 77), Rs)


def test_pairing_groups_partition_the_tracked_step4_voxels() -> None:
    vox = Q.voxels()
    grp = Q.pairing(vox)
    assert len(vox) == 200
    for arm in Q.ARMS:
        total = sum(int(m.sum()) for m in grp[arm].values())
        assert total == 200
        assert (
            sum(int((grp[arm][a] & grp[arm][b]).sum()) for a in Q.GROUPS for b in Q.GROUPS if a < b)
            == 0
        )

    # the Step 4 splits of no drop against clean (lenient 31 / 12 / 10, ray 35 / 9 / 5, C++ 32 / 14
    # / 24)
    def count(arm: str) -> tuple:
        return tuple(int(grp[arm][g].sum()) for g in ("fixed", "broken", "both wrong"))

    assert count("off") == (31, 12, 10)
    assert count("rayR2") == (35, 9, 5)
    assert count("on") == (32, 14, 24)


def test_bootstrap_of_identical_image_sets_is_one() -> None:
    rng = np.random.default_rng(0)
    base = np.sort(rng.random((12, 1 + len(Q.RADII), 5)), axis=1)  # rising curves
    r = Q.boot_width_diff(base, base.copy(), n_boot=20, seed=1)
    assert r["ratio"] == pytest.approx(1.0) and r["diff"] == pytest.approx(0.0)
    assert r["ratio_lo"] == pytest.approx(1.0) and r["ratio_hi"] == pytest.approx(1.0)


def test_pass_curve_counts_positive_overlaps() -> None:
    a = np.zeros((4, 1 + len(Q.RADII), 3))
    a[:, 0, :] = 5
    a[0, 1, :] = 1
    pc = Q.pass_curve(a)
    assert pc[0] == 1.0 and pc[1] == pytest.approx(0.25) and pc[2] == 0.0
    assert Q.pass_curve(a, np.array([True, False, False, False]))[1] == 1.0


def test_group_top_quarter_counts_carry_wilson_intervals() -> None:
    """The group-level top-quarter counts are given with a Wilson interval (small groups)."""
    import json

    assert Q.split4(dict(both=20, only_nodrop=3, only_clean=2, neither=6)) == "20 / 3 / 2 / 6"
    assert Q.count_ci(20, 31) == "20 of 31 (47 to 79%)"
    summ = json.loads((ROOT / "benchmarks" / "q8_nodrop" / "step7" / "summary.json").read_text())
    for arm, groups in summ["groups_paired"].items():
        for g, rec in groups.items():
            k, n = rec["keep"], rec["n_vox"]
            assert k["both"] + k["only_nodrop"] + k["only_clean"] + k["neither"] == n
            assert rec["keep_both_ci"] == Q.count_ci(k["both"], n)
            assert rec["keep_clean_ci"] == Q.count_ci(k["both"] + k["only_clean"], n)
            assert rec["keep_nodrop_ci"] == Q.count_ci(k["both"] + k["only_nodrop"], n)
    assert summ["groups_paired"]["off"]["fixed"]["keep_clean_ci"] == "22 of 31 (53 to 84%)"
