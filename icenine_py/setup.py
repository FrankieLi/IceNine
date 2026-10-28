"""Setup script for C extensions (numpy include dirs)."""

from setuptools import setup, Extension
import numpy as np

# No extension-wide -ffp-contract=off (D-20261027-I9, superseding part 4 of D-20261027-I8): the
# flag also changes the key-off _rasterize.c (on clang arm64 it removes the multiply-adds of the
# Sutherland-Hodgman intersections, so triangle_overlap moves on some border-crossing triangles).
# Contraction is off in _fused_eval.c only, by its per-file pragmas. Validated on macOS clang only.
ext_modules = [
    Extension(
        "icenine._rasterize",
        # _fused_eval.c: the fused evaluate of the FusedCost key (S2); ICENINE_FUSED adds its
        # method table entries and the HAS_FUSED constant in _rasterize.c
        sources=["icenine/_rasterize.c", "icenine/_fused_eval.c"],
        # _stage_d.h: the declarations both files share (an edit rebuilds both)
        depends=["icenine/_stage_d.h"],
        include_dirs=[np.get_include()],
        define_macros=[("ICENINE_FUSED", "1")],
    ),
]

setup(ext_modules=ext_modules)
