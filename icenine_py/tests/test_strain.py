"""Strain in the forward model: ``SpotModel.predict(..., strain=)`` (GO5).

Convention: crystal-frame symmetric strain eps [e11 e22 e33 e23 e13 e12], M = (I + eps)^-1,
g = O M g0, Bragg angle and omega from the strained |g|."""

import sys
from pathlib import Path
from typing import Any, Tuple

import numpy as np
import torch

from icenine import grain_fit as GF
from icenine.diffraction_core import KEV_OVER_HBAR_C_IN_ANG

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))
DT = torch.float64


def _setup() -> Tuple[Any, Any, torch.Tensor, torch.Tensor]:
    from test_grain_fit import _sc, _toy

    from icenine.reconstructor import _get_voxel_vertices

    rec, voxel, _ = _toy()
    sc = _sc(rec)
    model = GF.SpotModel(sc, phase=int(voxel.phase))
    mic = rec.setup.sample.get_mic()
    Rs = np.stack([np.asarray(v.orientation, np.float64) for v in mic.voxels])
    Vs = np.stack([_get_voxel_vertices(v).numpy().astype(np.float64) for v in mic.voxels])
    return sc, model, torch.as_tensor(Rs), torch.as_tensor(Vs)


def _same(a: Any, b: Any) -> bool:
    return all(
        torch.equal(getattr(a, f), getattr(b, f))
        for f in ("vox", "refl", "branch", "fpos", "frame", "col", "row", "hit")
    )


def _eps(n: int, **comp: float) -> torch.Tensor:
    idx = {"e11": 0, "e22": 1, "e33": 2, "e23": 3, "e13": 4, "e12": 5}
    e = torch.zeros(n, 6, dtype=DT)
    for k, v in comp.items():
        e[:, idx[k]] = v
    return e


def test_strain_to_m_roundtrip_and_component_layout() -> None:
    e = torch.tensor([[1e-3, 2e-3, 3e-3, 4e-4, 5e-4, 6e-4]], dtype=DT)
    M = GF.strain_to_m(e)[0]
    E = torch.linalg.inv(M) - torch.eye(3, dtype=DT)
    assert torch.allclose(E, E.T)
    assert abs(float(E[1, 2]) - 4e-4) < 1e-15
    assert abs(float(E[0, 2]) - 5e-4) < 1e-15 and abs(float(E[0, 1]) - 6e-4) < 1e-15
    assert np.abs(GF.m_to_strain(M.numpy()[None]) - e.numpy()).max() < 1e-15


def test_none_bit_identical_and_zero_strain_equal() -> None:
    _, model, Rs, Vs = _setup()
    a = model.predict(Rs, Vs)
    assert _same(a, model.predict(Rs, Vs, None))
    z = model.predict(Rs, Vs, torch.zeros(len(Rs), 6, dtype=DT))
    assert torch.equal(a.vox, z.vox) and torch.equal(a.refl, z.refl)
    assert torch.allclose(a.col, z.col, atol=1e-9) and torch.allclose(a.row, z.row, atol=1e-9)


def test_pure_dilation_scales_g_and_shrinks_two_theta() -> None:
    sc, model, Rs, Vs = _setup()
    d = 1e-3
    M = GF.strain_to_m(_eps(1, e11=d, e22=d, e33=d))
    g_s = torch.linalg.norm(torch.einsum("ij,kj->ki", Rs[0] @ M[0], model.g_hkl), dim=1)
    g0 = torch.linalg.norm(torch.einsum("ij,kj->ki", Rs[0], model.g_hkl), dim=1)
    assert torch.allclose(g_s, g0 / (1 + d), rtol=1e-12)
    assert torch.allclose(g0, model.g_mag, rtol=1e-5)
    k = KEV_OVER_HBAR_C_IN_ANG * sc.energy
    assert bool((torch.asin(g_s / (2 * k)) < torch.asin(g0 / (2 * k))).all())
    # the spot model with eps = d I equals the unstrained model with g0 scaled by 1/(1+d)
    strained = model.predict(Rs, Vs, _eps(len(Rs), e11=d, e22=d, e33=d))
    base = model.predict(Rs, Vs)
    saved = (model.g_hkl, model.g_mag, model._consts)
    model.g_hkl, model.g_mag, model._consts = saved[0] / (1 + d), saved[1] / (1 + d), {}
    ref = model.predict(Rs, Vs)
    model.g_hkl, model.g_mag, model._consts = saved
    assert torch.equal(strained.vox, ref.vox) and torch.equal(strained.refl, ref.refl)
    assert torch.allclose(strained.col, ref.col, atol=1e-6)
    assert torch.allclose(strained.row, ref.row, atol=1e-6)
    # the spots really moved relative to the unstrained prediction
    if torch.equal(base.refl, strained.refl):
        assert float((base.col - strained.col).abs().max()) > 1e-2


def test_shear_matches_finite_differences() -> None:
    """Spot columns/rows move linearly with a small shear: the central difference equals the
    autograd derivative (differentiable in eps) and predicts a 1e-4 step."""
    _, model, Rs, Vs = _setup()
    n = len(Rs)
    h = 1e-6

    def pos(eps: float, comp: str) -> Tuple[Any, torch.Tensor]:
        r = model.predict(Rs, Vs, _eps(n, **{comp: eps}))
        return r, torch.cat([r.col.reshape(-1), r.row.reshape(-1)])

    def pos_t(x: torch.Tensor, comp: str) -> torch.Tensor:
        sel = torch.zeros(6, dtype=DT)
        sel[{"e23": 3, "e12": 5}[comp]] = 1.0
        r = model.predict(Rs, Vs, (x * sel)[None].expand(n, 6))
        return torch.cat([r.col.reshape(-1), r.row.reshape(-1)])

    for comp in ("e23", "e12"):
        r0, p0 = pos(0.0, comp)
        rp, pp = pos(h, comp)
        rm, pm = pos(-h, comp)
        assert torch.equal(rp.refl, r0.refl) and torch.equal(rm.refl, r0.refl)
        fd = (pp - pm) / (2 * h)
        _, jv = torch.autograd.functional.jvp(
            lambda x: pos_t(x, comp), torch.zeros((), dtype=DT), torch.ones((), dtype=DT)
        )
        assert float(fd.abs().max()) > 1.0
        assert float((jv - fd).abs().max()) < 1e-6 * float(fd.abs().max()) + 1e-4
        rb, pb = pos(1e-4, comp)
        if torch.equal(rb.refl, r0.refl):
            lin = 1e-4 * fd
            assert float((pb - p0 - lin).abs().max()) < 0.05 * float(lin.abs().max()) + 1e-3
