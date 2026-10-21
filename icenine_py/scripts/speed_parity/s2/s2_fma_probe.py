"""S2 Step 0: which float32 operation order reproduces each small torch product bit for bit.

``s2fused.c`` computes the products of ``evaluate`` in float32 with ``fmaf`` where torch's kernel
fuses a multiply and an add. This probe is the evidence for that choice. For each product, on the
inputs ``evaluate`` gives it (ThreeVoxels voxel 0, the triangle cost at Q-max 8, ``--n``
orientations of S1's set), it runs the torch call the library makes and compares the result,
element by element and bit for bit, with candidate orders computed with libm's ``fmaf`` (through
ctypes, so the fused multiply-add is exact) and plain float32 products and sums:

* ``plain``: ``(a0*b0 + a1*b1) + a2*b2`` (each step rounded to float32);
* ``fma_chain``: ``fmaf(a2, b2, fmaf(a1, b1, a0*b0))`` (the order ``s2fused.c`` uses);
* ``fma_rev``: ``fmaf(a0, b0, fmaf(a1, b1, a2*b2))``.

For the 4 x 4 vertex einsum the chain has a fourth term ``fmaf(t, 1, ...)`` (the translation), for
``torch.norm`` the terms are the squares. Products: g_lab ``(R @ g_hkl.T).T``, ``Rz @ base``,
``bmm(full_rot, normals)``, the vertex einsum, ``torch.norm`` of the g vectors and of the reflected
directions. Not a timing run (no preflight).

  uv run python scripts/speed_parity/s2/s2_fma_probe.py --out benchmarks/s2_step0/fmaf_probe.json
"""

import argparse
import contextlib
import ctypes
import ctypes.util
import io
import json
import platform
import sys
from pathlib import Path
from typing import Any, Callable, Dict, List

import numpy as np
import torch

HERE = Path(__file__).resolve().parent
ICE = HERE.parents[2]
sys.path.insert(0, str(ICE))
sys.path.insert(0, str(HERE.parent))
sys.path.insert(0, str(HERE))

import s1_parity as S1  # noqa: E402
import s2_parity as SP  # noqa: E402
from icenine.diffraction_core import BraggConstants, scattering_omegas_fixed_beam  # noqa: E402

_libm = ctypes.CDLL(ctypes.util.find_library("m"))
_libm.fmaf.restype = ctypes.c_float
_libm.fmaf.argtypes = [ctypes.c_float, ctypes.c_float, ctypes.c_float]

F32 = np.float32


def fmaf(a: np.ndarray, b: np.ndarray, c: np.ndarray) -> np.ndarray:
    """Element-wise exact float32 fused multiply-add (libm ``fmaf``)."""
    a, b, c = np.broadcast_arrays(F32(a), F32(b), F32(c))
    out = np.empty(a.shape, F32)
    fa, fb, fc, fo = a.reshape(-1), b.reshape(-1), c.reshape(-1), out.reshape(-1)
    for i in range(fo.size):
        fo[i] = _libm.fmaf(float(fa[i]), float(fb[i]), float(fc[i]))
    return out


def dot3(a: List[np.ndarray], b: List[np.ndarray], order: str) -> np.ndarray:
    """a0*b0 + a1*b1 + a2*b2 (+ a3*b3 if given) in float32 in the named order."""
    a = [F32(x) for x in a]
    b = [F32(x) for x in b]
    if order == "plain":
        s = a[0] * b[0]
        for i in range(1, len(a)):
            s = F32(s + a[i] * b[i])
        return s
    if order == "fma_chain":
        s = a[0] * b[0]
        for i in range(1, len(a)):
            s = fmaf(a[i], b[i], s)
        return s
    if order == "fma_rev":
        s = a[-1] * b[-1]
        for i in range(len(a) - 2, -1, -1):
            s = fmaf(a[i], b[i], s)
        return s
    raise ValueError(order)


ORDERS = ("plain", "fma_chain", "fma_rev")


def compare(torch_out: np.ndarray, cand: Callable[[str], np.ndarray]) -> Dict[str, Dict[str, Any]]:
    t = np.ascontiguousarray(torch_out, F32)
    res = {}
    for o in ORDERS:
        c = np.ascontiguousarray(cand(o), F32)
        eq = int((t.view(np.int32) == c.view(np.int32)).sum())
        res[o] = dict(n=int(t.size), equal=eq, max_ulp=int(SP.ulp_diff(t, c).max()))
    return res


def main() -> int:
    torch.set_num_threads(1)
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--out", required=True)
    ap.add_argument("--n", type=int, default=20, help="orientations (half near, half FZ)")
    a = ap.parse_args()
    out_path = Path(a.out).resolve()  # load_dataset changes the working directory
    with contextlib.redirect_stdout(io.StringIO()):
        ds = S1.load_dataset("three0")
    R_all = S1.orientation_set(ds, a.n // 2, a.n - a.n // 2)
    cf = ds.cost_fn(SP.MODES["local"])
    g_hkl, g_mag = cf._phase_recip_vecs[ds.phase]
    consts = BraggConstants(g_mag, cf.simulator.beam_energy, cf.simulator.beam_deflection_chi)
    s2l = cf.sample.sample_to_lab_matrix
    base = s2l[:3, :3]
    trans = s2l[:3, 3]
    beam = cf.simulator.beam_direction
    verts = ds.verts
    acc: Dict[str, Dict[str, Dict[str, int]]] = {}

    def add(name: str, r: Dict[str, Dict[str, Any]]) -> None:
        d = acc.setdefault(name, {o: dict(n=0, equal=0, max_ulp=0) for o in ORDERS})
        for o in ORDERS:
            d[o]["n"] += r[o]["n"]
            d[o]["equal"] += r[o]["equal"]
            d[o]["max_ulp"] = max(d[o]["max_ulp"], r[o]["max_ulp"])

    G = g_hkl.numpy()
    for R in R_all:
        Rt = torch.from_numpy(R).float()
        Rn = Rt.numpy()
        g_lab = (Rt @ g_hkl.T).T  # cost_functions.evaluate
        gl = g_lab.numpy()
        add("g_lab (R @ g_hkl.T).T", compare(
            gl, lambda o: np.stack([dot3([Rn[i, j] for j in range(3)],
                                         [G[:, j] for j in range(3)], o)
                                    for i in range(3)], axis=1)))  # fmt: skip
        br = scattering_omegas_fixed_beam(g_lab, consts)
        obs = br.observable
        if not obs.any():
            continue
        og = g_lab[obs]
        mag = torch.norm(og, dim=1, keepdim=True)
        ogn = og.numpy()
        add("torch.norm(g_lab)", compare(
            mag.numpy()[:, 0],
            lambda o: np.sqrt(dot3([ogn[:, j] for j in range(3)],
                                   [ogn[:, j] for j in range(3)], o)).astype(F32)))  # fmt: skip
        nrm = torch.cat([og / mag, og / mag])
        om = torch.cat([br.omega1[obs], br.omega2[obs]])
        N = om.shape[0]
        Rz = torch.zeros(N, 3, 3)
        c, s = torch.cos(om), torch.sin(om)
        Rz[:, 0, 0], Rz[:, 0, 1], Rz[:, 1, 0], Rz[:, 1, 1], Rz[:, 2, 2] = c, -s, s, c, 1.0
        full_rot = Rz @ base
        Rzn, Bn, Fn = Rz.numpy(), base.numpy(), full_rot.numpy()
        add("Rz @ base", compare(
            Fn, lambda o: np.stack([np.stack([dot3([Rzn[:, i, k] for k in range(3)],
                                                   [Bn[k, j] for k in range(3)], o)
                                              for j in range(3)], axis=1)
                                    for i in range(3)], axis=1)))  # fmt: skip
        ln = torch.bmm(full_rot, nrm.unsqueeze(-1)).squeeze(-1)
        nn = nrm.numpy()
        add("bmm(full_rot, normals)", compare(
            ln.numpy(), lambda o: np.stack([dot3([Fn[:, i, k] for k in range(3)],
                                                 [nn[:, k] for k in range(3)], o)
                                            for i in range(3)], axis=1)))  # fmt: skip
        dot = (beam.unsqueeze(0) * ln).sum(dim=1, keepdim=True)
        refl = beam.unsqueeze(0) - 2.0 * dot * ln
        rf = refl.numpy()
        add("torch.norm(reflected)", compare(
            torch.norm(refl, dim=1).numpy(),
            lambda o: np.sqrt(dot3([rf[:, j] for j in range(3)],
                                   [rf[:, j] for j in range(3)], o)).astype(F32)))  # fmt: skip
        full_4x4 = torch.zeros(N, 4, 4)
        full_4x4[:, :3, :3] = full_rot
        full_4x4[:, :3, 3] = trans
        full_4x4[:, 3, 3] = 1.0
        v4 = torch.cat([verts, torch.ones(3, 1)], dim=1)
        lv = torch.einsum("mij,vj->mvi", full_4x4, v4)[:, :, :3]
        M4, V4 = full_4x4.numpy(), v4.numpy()
        add("einsum vertices (4 x 4)", compare(
            lv.numpy(), lambda o: np.stack([np.stack([dot3([M4[:, i, k] for k in range(4)],
                                                           [V4[v, k] for k in range(4)], o)
                                                      for i in range(3)], axis=1)
                                            for v in range(3)], axis=1)))  # fmt: skip
    rows = [
        dict(product=name, order=o, n=d[o]["n"], equal=d[o]["equal"],
             equal_frac=d[o]["equal"] / max(d[o]["n"], 1), max_ulp=d[o]["max_ulp"])
        for name, d in acc.items() for o in ORDERS
    ]  # fmt: skip
    out = dict(
        meta=dict(
            head=S1.git_head(),
            dirty=S1.git_dirty(),
            dirty_files=SP.dirty_files(),
            torch=torch.__version__,
            numpy=np.__version__,
            machine=platform.machine(),
            platform=platform.platform(),
            dataset="three0",
            mode="local",
            n_orientations=len(R_all),
        ),
        rows=rows,
    )
    out_path.write_text(json.dumps(out, indent=1) + "\n")
    for r in rows:
        print(f"  {r['product']:28s} {r['order']:9s} {r['equal']:7d} / {r['n']:7d}  "
              f"max {r['max_ulp']} ulp")  # fmt: skip
    return 0


if __name__ == "__main__":
    sys.exit(main())
