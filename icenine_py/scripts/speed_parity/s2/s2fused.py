"""S2 Step 0: build and wrap the scratch fused-evaluate extension (``s2fused.c``).

``build()`` compiles ``s2fused.c`` with the system ``cc`` and the flags setuptools uses for
``icenine._rasterize`` (``-O3 -fwrapv -DNDEBUG``, native architecture only) into
``scripts/speed_parity/cache/s2/`` (gitignored). The C file includes the library's ``_rasterize.c``,
so the stage-D helpers (triangle and pixel-radius overlap, qualified-peak count) are the library's;
the per-peak stage-D loop around them is copied (without the ray-rule block). A change to either
file triggers a rebuild.

``FusedEval(cost_fn)`` is the prototype of the planned key: ``evaluate(orientation, verts, phase)``
returns the ``OverlapInfo`` of ``cost_fn.evaluate`` from one C call. Python keeps the key
dispatch, ``eval_count``, building the ``OverlapInfo`` and the lazy ``FrameTable`` fill (the C call
returns the frames it needs; ``ensure`` fills them and the call is repeated). Covered: hard mode
without a ray rule (pixel radius 0 and > 0, centroid and vertex 0, default rule, ``CppOverlapRule``,
``CppEdgeRules``, ``min_sin_eta``). Anything else raises ``Unsupported``.
"""

import importlib
import subprocess
import sys
import sysconfig
from pathlib import Path
from typing import Any, Dict, List, Optional

import numpy as np
import torch

HERE = Path(__file__).resolve().parent
ICE = HERE.parents[2]
sys.path.insert(0, str(ICE))

from icenine.cost_functions import OverlapInfo, VoxelCostFunction  # noqa: E402
from icenine.diffraction_core import BraggConstants  # noqa: E402

SRC = HERE / "s2fused.c"
RAST = ICE / "icenine" / "_rasterize.c"
CACHE = ICE / "scripts" / "speed_parity" / "cache" / "s2"
CFLAGS = ["-O3", "-fwrapv", "-DNDEBUG", "-Wall", "-Wno-unused-function"]


class Unsupported(ValueError):
    """The cost function uses a mode the fused path does not cover."""


#: build variants: "portable" computes g_lab = R g_hkl with fmaf in C (not bit-identical to torch's
#: sgemm); "accelerate" (macOS only) calls Accelerate's cblas_sgemm as torch does on this machine
VARIANTS = {
    "portable": ("_s2fused", [], []),
    "accelerate": ("_s2fused_acc", ["-DS2_ACCELERATE"], ["-framework", "Accelerate"]),
}


def lib_path(variant: str = "portable") -> Path:
    return CACHE / f"{VARIANTS[variant][0]}{sysconfig.get_config_var('EXT_SUFFIX')}"


def build(variant: str = "portable", force: bool = False) -> Path:
    lib = lib_path(variant)
    lib.parent.mkdir(parents=True, exist_ok=True)
    newest = max(SRC.stat().st_mtime, RAST.stat().st_mtime)
    if force or not lib.exists() or lib.stat().st_mtime < newest:
        _, defs, libs = VARIANTS[variant]
        cmd = [
            "cc",
            *CFLAGS,
            *defs,
            "-shared",
            "-undefined",
            "dynamic_lookup",
            "-I",
            sysconfig.get_paths()["include"],
            "-I",
            np.get_include(),
            "-o",
            str(lib),
            str(SRC),
            *libs,
        ]
        subprocess.run(cmd, check=True)
    return lib


_ext: Dict[str, Any] = {}


def ext(variant: str = "portable") -> Any:
    if variant not in _ext:
        build(variant)
        if str(CACHE) not in sys.path:
            sys.path.insert(0, str(CACHE))
        _ext[variant] = importlib.import_module(VARIANTS[variant][0])
    return _ext[variant]


def unsupported_reason(cf: VoxelCostFunction) -> Optional[str]:
    if cf.mode != "hard":
        return f"mode {cf.mode!r}"
    if cf.ray_rule is not None:
        return "RayRule"
    if len(cf.detector_list) > 8:
        return "more than 8 detectors"
    return None


def frame_mode(low: Any, width: Any) -> int:
    """How numpy evaluates ``(omegas - low) / width`` for float32 omegas (NEP 50: a numpy float64
    scalar promotes, a Python float does not): 0 both steps in float64, 1 the subtraction in
    float32 and the division in float64, 2 both in float32."""
    a = np.zeros(1, np.float32) - low
    if a.dtype == np.float64:
        return 0
    return 1 if (a / width).dtype == np.float64 else 2


def det_geometry(cf: VoxelCostFunction) -> np.ndarray:
    rows = []
    for det in cf.detector_list:
        pl = det._detector_plane
        rows.append(
            np.concatenate(
                [
                    pl.normal.numpy().astype(np.float32),
                    np.array([float(pl.d)], np.float32),
                    det._position.numpy().astype(np.float32),
                    det._lab_frame_coord_origin.numpy().astype(np.float32),
                    det._lab_frame_basis_j.numpy().astype(np.float32),
                    det._lab_frame_basis_k.numpy().astype(np.float32),
                    np.array(
                        [
                            det.pixel_half_width,
                            det.pixel_width,
                            det.pixel_half_height,
                            det.pixel_height,
                        ],
                        np.float32,
                    ),
                ]
            )
        )
    return np.ascontiguousarray(np.stack(rows), dtype=np.float32)


class _Ctx:
    """The per-(cost function, phase) context and the arrays it points to (kept alive here)."""

    def __init__(self, cf: VoxelCostFunction, phase: int, variant: str) -> None:
        g_hkl_t, g_mag_t = cf._phase_recip_vecs[phase]
        consts = BraggConstants(g_mag_t, cf.simulator.beam_energy, cf.simulator.beam_deflection_chi)
        self.g_hkl = np.ascontiguousarray(g_hkl_t.numpy(), dtype=np.float32)
        self.sin_theta = np.ascontiguousarray(consts.sin_theta.numpy(), dtype=np.float32)
        self.g_mag_eps = np.ascontiguousarray(consts.g_mag_eps.numpy(), dtype=np.float32)
        scal = np.array(
            [
                float(consts.sin_chi_laue.reshape(-1)[0]),
                float(consts.cos_chi_laue.reshape(-1)[0]),
                np.float32(consts.epsilon),
                np.float32(cf.eta_limit),
                np.float32(cf.min_sin_eta),
            ],
            np.float32,
        )
        # all entries of the expanded chi tensors are the same value
        assert bool((consts.sin_chi_laue == consts.sin_chi_laue.reshape(-1)[0]).all())
        s2l = cf.sample.sample_to_lab_matrix
        beam = cf.simulator.beam_direction.numpy().astype(np.float32)
        base = np.ascontiguousarray(s2l[:3, :3].numpy(), dtype=np.float32)
        trans = np.ascontiguousarray(s2l[:3, 3].numpy(), dtype=np.float32)
        self.table = cf.frame_table()
        self.frame_index = np.ascontiguousarray(cf._frame_index, dtype=np.int64)
        rmap = cf.range_map
        rule_bits = (1 if cf.cpp_overlap_rule else 0) | (2 if cf.cpp_edge_rules else 0)
        ints = np.array(
            [
                cf.pixel_radius,
                int(bool(cf.center_screen)),
                rule_bits,
                int(cf.min_sin_eta > 0.0),
                frame_mode(rmap.low, rmap.width),
                len(cf.detector_list),
            ],
            np.int32,
        )
        self.capsule = ext(variant).make_context(
            self.g_hkl,
            self.sin_theta,
            self.g_mag_eps,
            scal,
            beam,
            base,
            trans,
            ints,
            np.array([float(rmap.low), float(rmap.width)], np.float64),
            self.frame_index,
            det_geometry(cf),
            self.table,
            self.table.images,
            self.table._filled,
        )


class FusedEval:
    """``cost_fn.evaluate`` through one C call (prototype of the ``FusedCost`` key)."""

    def __init__(self, cf: VoxelCostFunction, variant: str = "portable") -> None:
        reason = unsupported_reason(cf)
        if reason is not None:
            raise Unsupported(reason)
        self.cf = cf
        self.variant = variant
        self._ctx: Dict[int, _Ctx] = {}
        self._vt: Optional[torch.Tensor] = None
        self._vn: Optional[np.ndarray] = None
        self._evaluate = ext(variant).evaluate

    def context(self, phase: int) -> Optional[_Ctx]:
        ctx = self._ctx.get(phase)
        if ctx is None or ctx.table is not self.cf.frame_table():
            if phase not in self.cf._phase_recip_vecs:
                return None
            ctx = self._ctx[phase] = _Ctx(self.cf, phase, self.variant)
        return ctx

    def _verts(self, voxel_vertices: torch.Tensor) -> np.ndarray:
        if voxel_vertices is not self._vt:
            self._vt = voxel_vertices
            self._vn = np.ascontiguousarray(voxel_vertices.detach().numpy(), dtype=np.float32)
        assert self._vn is not None
        return self._vn

    def evaluate(
        self, orientation: np.ndarray, voxel_vertices: torch.Tensor, phase_index: int = 0
    ) -> OverlapInfo:
        self.cf.eval_count += 1
        ctx = self.context(phase_index)
        if ctx is None:
            return OverlapInfo()
        v = self._verts(voxel_vertices)
        r = self._evaluate(ctx.capsule, orientation, v)
        while type(r) is list:  # frames to fill first (lazy FrameTable), then call again
            ctx.table.ensure(np.asarray(r, dtype=np.int64))
            r = self._evaluate(ctx.capsule, orientation, v)
        oi = OverlapInfo()
        (
            oi.pixel_overlap,
            oi.pixel_on_detector,
            oi.peak_overlap,
            oi.peak_on_detector,
            oi.quality,
            oi.n_quality_points,
            oi.detectors_overlap,
        ) = r
        return oi

    def debug(
        self, orientation: np.ndarray, voxel_vertices: torch.Tensor, phase_index: int = 0
    ) -> Dict[str, Any]:
        """Intermediates of every (solution, reflection) row (row s * K + k) and the result; the
        frames are filled first, so ``result`` is always set."""
        ctx = self.context(phase_index)
        assert ctx is not None
        v = self._verts(voxel_vertices)
        d: Dict[str, Any] = self._evaluate(ctx.capsule, orientation, v, 1)
        while "missing" in d:
            ctx.table.ensure(np.asarray(d["missing"], dtype=np.int64))
            d = self._evaluate(ctx.capsule, orientation, v, 1)
        return d


def result_tuple(oi: OverlapInfo) -> tuple:
    return (
        oi.pixel_overlap,
        oi.pixel_on_detector,
        oi.peak_overlap,
        oi.peak_on_detector,
        oi.quality,
        oi.n_quality_points,
        oi.detectors_overlap,
    )


FIELDS: List[str] = [
    "pixel_overlap",
    "pixel_on_detector",
    "peak_overlap",
    "peak_on_detector",
    "quality",
    "n_quality_points",
    "detectors_overlap",
]


if __name__ == "__main__":
    for v in VARIANTS:
        print(build(v, force="--force" in sys.argv))
