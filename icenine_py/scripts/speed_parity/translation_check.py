"""Check: for fixed orientation, moving the voxel by dx translates each peak's projected
vertices rigidly, by an amount linear in dx, identical on parallel detectors.

Derivation: docs/spot_translation_derivation.md.
Usage (from icenine_py/): uv run python scripts/speed_parity/translation_check.py
"""

import sys
from pathlib import Path
import numpy as np
import torch

sys.path.insert(0, str(Path(__file__).resolve().parent))
import prof_eval as P
from icenine import cost_functions as cf

fn = P.make(0, 8.0, True)
cap = {}
oc = cf._c_stage_d_overlap


def spy(*a):
    cap["verts"] = np.array(a[4], dtype=np.float64)  # (M, n_det, 3, 2) [col,row]
    cap["wedge"] = np.array(a[1])
    return oc(*a)


cf._c_stage_d_overlap = spy


def project(o, shift_um):
    v = P.verts + torch.tensor(shift_um, dtype=torch.float32) * 1e-3  # mm
    fn.evaluate(o, v, P.vox.phase)
    return cap["verts"].copy(), cap["wedge"].copy()


rng = np.random.default_rng(0)
worst = dict(rigid=0.0, lin=0.0, add=0.0, det=0.0)
for o in P.orients[:20]:
    base, w0 = project(o, [0, 0, 0])
    a = rng.uniform(-150, 150, 3)
    a[2] = 0.0  # in-layer move, um
    b = rng.uniform(-150, 150, 3)
    b[2] = 0.0
    va, wa = project(o, a)
    v2a, _ = project(o, 2 * a)
    vb, _ = project(o, b)
    vab, _ = project(o, a + b)
    assert (wa == w0).all(), "frame indices changed"
    da, d2a, db, dab = va - base, v2a - base, vb - base, vab - base  # (M, det, 3, 2)
    worst["rigid"] = max(
        worst["rigid"], np.abs(da - da[:, :, :1, :]).max()
    )  # same shift for all 3 vertices
    worst["lin"] = max(worst["lin"], np.abs(d2a - 2 * da).max())
    worst["add"] = max(worst["add"], np.abs(dab - da - db).max())
    worst["det"] = max(worst["det"], np.abs(da[:, 0] - da[:, 1]).max())  # det 0 vs det 1 shift
    typ = np.abs(da).max()
print(f"frame indices unchanged under in-layer moves: yes (20 orientations)")
print(f"typical max shift for a ~150 um move: {typ:.1f} px")
for k, v in worst.items():
    print(f"max deviation [{k}]: {v:.2e} px")
