"""S3 Step 2 side checks (not part of the library): how the lockstep quick MC maps its uniforms and
builds its trial rotations.

1. ``uniform`` mapping. Recomputing ``low + (high - low) * u`` with numpy ufuncs (two roundings)
   against the scalar ``Generator.uniform(-r, r)`` of the sequential path, against the same
   expression rounded once (exact fused multiply-add, via fractions), and against
   ``Generator.uniform(-r, r, size=n)`` (what ``draw_trial_offsets`` uses).
2. Vectorised trial-rotation build against the scalar one that ``MCOptimizer._mc_block`` runs:
   the quaternion product and the batch quaternion-to-matrix conversion, and a numpy-vectorised
   ``QuaternionGrid.get_near_identity_point``.

  uv run python scripts/speed_parity/s3_draw_check.py --out benchmarks/s3_steps1_5
"""

import argparse
import json
import math
import sys
from fractions import Fraction
from pathlib import Path
from typing import Any, Dict

import numpy as np
from scipy.spatial.transform import Rotation

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1]))

from icenine.sampling import QuaternionGrid, _quat_multiply, quaternion_to_matrix  # noqa: E402


def mapping_check(n: int = 20000, step_deg: float = 0.2) -> Dict[str, Any]:
    r = math.tan(math.radians(step_deg)) / math.sqrt(12.0)
    a, b, c = (np.random.default_rng(123) for _ in range(3))
    u = a.random(n)
    scalar = np.array([b.uniform(-r, r) for _ in range(n)])
    sized = c.uniform(-r, r, size=n)
    low, high = -r, r
    separate = low + (high - low) * u
    fused = np.array([float(Fraction(high - low) * Fraction(x) + Fraction(low)) for x in u])
    return dict(
        n=n,
        radius=r,
        separate_roundings_differ_from_scalar=int((separate != scalar).sum()),
        fused_multiply_add_differs_from_scalar=int((fused != scalar).sum()),
        uniform_size_differs_from_scalar=int((sized != scalar).sum()),
        generator_states_equal=bool(a.bit_generator.state == b.bit_generator.state),
    )


def _vslerp(q1: np.ndarray, q2: np.ndarray, t: np.ndarray) -> np.ndarray:
    q1 = np.broadcast_to(q1, (len(t), 4))
    q2 = np.broadcast_to(q2, (len(t), 4))
    cos_theta = np.einsum("ij,ij->i", q1, q2)
    flip = cos_theta < 0.0
    cos_theta = np.where(flip, -cos_theta, cos_theta)
    near = 1.0 - cos_theta < 0.01
    theta = np.arccos(np.where(near, 0.0, cos_theta))
    sin_theta = np.sin(theta)
    with np.errstate(all="ignore"):
        beta = np.where(near, 1.0 - t, np.sin(theta - t * theta) / sin_theta)
        alpha = np.where(near, t, np.sin(t * theta) / sin_theta)
    alpha = np.where(flip, -alpha, alpha)
    return beta[:, None] * q1 + alpha[:, None] * q2


def build_check(n: int = 4000, step_deg: float = 0.2) -> Dict[str, Any]:
    rng = np.random.default_rng(0)
    r = math.tan(math.radians(step_deg)) / math.sqrt(12.0)
    xyz = rng.uniform(-r, r, size=(n, 3))
    q = Rotation.random(n, random_state=1).as_quat()[:, [3, 0, 1, 2]]
    grid = QuaternionGrid()
    delta = np.stack([grid.get_near_identity_point(*xyz[i]) for i in range(n)])
    prod = np.stack([_quat_multiply(delta[i], q[i]) for i in range(n)])
    mat = np.stack([quaternion_to_matrix(prod[i]) for i in range(n)])
    w1, x1, y1, z1 = delta.T
    w2, x2, y2, z2 = q.T
    prod_v = np.stack(
        [
            w1 * w2 - x1 * x2 - y1 * y2 - z1 * z2,
            w1 * x2 + x1 * w2 + y1 * z2 - z1 * y2,
            w1 * y2 + y1 * w2 + z1 * x2 - x1 * z2,
            w1 * z2 + z1 * w2 + x1 * y2 - y1 * x2,
        ],
        axis=1,
    )
    mat_v = Rotation.from_quat(prod[:, [1, 2, 3, 0]]).as_matrix()
    v = grid._hyperface_vertices[0]
    pos = xyz + 0.5
    a, b, c = pos[:, 0], pos[:, 1], pos[:, 2]
    y1_ = _vslerp(_vslerp(v[0], v[1], a), _vslerp(v[4], v[5], a), b)
    y2_ = _vslerp(_vslerp(v[2], v[3], a), _vslerp(v[6], v[7], a), b)
    res = _vslerp(y1_, y2_, c)
    delta_v = np.where((res[:, 0] < 0)[:, None], -res, res)
    return dict(
        n=n,
        quaternion_product_vectorised_differs=int((prod_v != prod).any(axis=1).sum()),
        quaternion_to_matrix_batch_differs=int((mat_v != mat).any(axis=(1, 2)).sum()),
        near_identity_vectorised_differs=int((delta_v != delta).any(axis=1).sum()),
        near_identity_max_abs_difference=float(np.abs(delta_v - delta).max()),
    )


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    res = dict(mapping=mapping_check(), build=build_check())
    out = Path(a.out) / "draw_check.json"
    out.write_text(json.dumps(res, indent=1))
    print(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()
