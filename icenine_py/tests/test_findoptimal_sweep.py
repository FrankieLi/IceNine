"""scripts/findoptimal_sweep.py: the cubic-symmetry-reduced misorientation used for the global
reconstruction's error."""

import sys
from pathlib import Path

import numpy as np
from scipy.spatial.transform import Rotation

sys.path.insert(0, str(Path(__file__).parent.parent / "scripts"))
import findoptimal_sweep as fso  # noqa: E402

from icenine.sampling import get_misorientation, matrix_to_quaternion  # noqa: E402


def test_symmetry_equivalent_orientation_has_zero_reduced_error():
    ops = fso.cubic_rotations()
    R_true = Rotation.random(random_state=5).as_matrix()
    for S in ops:
        assert fso.misorientation_deg(R_true @ S, R_true) < 1e-5
    # an equivalent orientation that is also slightly rotated: reduced = the small rotation
    small = Rotation.from_rotvec(np.radians([0.3, -0.2, 0.1])).as_matrix()
    R_est = small @ R_true @ ops[7]
    expect = np.degrees(np.linalg.norm(np.radians([0.3, -0.2, 0.1])))
    assert abs(float(fso.misorientation_deg(R_est, R_true)) - expect) < 1e-5
    # the plain angle of a non-identity operator is >= 90 deg
    plain = fso.misorientation_deg(R_true @ ops[7], R_true, reduce=False)
    assert plain >= 90.0 - 1e-6 or np.allclose(ops[7], np.eye(3))


def test_reduced_equals_plain_for_small_errors_and_matches_project_misorientation():
    rng = np.random.default_rng(0)
    R_true = Rotation.random(20, random_state=1).as_matrix()
    e = Rotation.from_rotvec(np.radians(rng.normal(size=(20, 3)) * 2.0)).as_matrix() @ R_true
    sym = fso.misorientation_deg(e, R_true)
    plain = fso.misorientation_deg(e, R_true, reduce=False)
    np.testing.assert_allclose(sym, plain, atol=1e-9)
    # a larger random pair: agree with the project's quaternion routine (cubic, 24 proper ops)

    quats = np.array([matrix_to_quaternion(S) for S in fso.cubic_rotations()])
    A = Rotation.random(10, random_state=3).as_matrix()
    B = Rotation.random(10, random_state=4).as_matrix()
    ours = fso.misorientation_deg(A, B)
    ref = [
        np.degrees(get_misorientation(matrix_to_quaternion(b), matrix_to_quaternion(a), quats))
        for a, b in zip(A, B)
    ]
    np.testing.assert_allclose(ours, ref, atol=1e-6)


def test_csl_classification_of_known_boundaries():
    import findoptimal_sweep_summary as fss

    ops = fso.cubic_rotations()
    R_true = Rotation.random(random_state=11).as_matrix()
    cases = {
        3: (60.0, [1, 1, 1]),
        5: (36.8699, [1, 0, 0]),
        7: (38.2132, [1, 1, 1]),
        11: (50.4788, [1, 1, 0]),
    }
    for sigma, (angle, axis) in cases.items():
        axis = np.array(axis, dtype=float) / np.linalg.norm(axis)
        M = Rotation.from_rotvec(np.radians(angle) * axis).as_matrix()
        R_final = R_true @ M @ ops[13]  # a symmetry-equivalent copy of the boundary orientation
        c = fss.classify_csl(R_true, R_final)
        assert c["sigma"] == sigma, (sigma, c)
        assert abs(c["angle"] - angle) < 1e-3 and c["deviation_deg"] < 1e-2
    # a generic rotation (37.9 deg about an axis 25 deg from every low-index direction) matches none
    ax = np.array([0.873, 0.358, 0.332])
    M = Rotation.from_rotvec(np.radians(37.91) * ax / np.linalg.norm(ax)).as_matrix()
    assert fss.classify_csl(R_true, R_true @ M)["sigma"] == 0
