"""Setup script for C extensions (numpy include dirs)."""

import sys

from setuptools import setup, Extension
import numpy as np

# no multiply-add contraction anywhere in the extension (D-20261027-I8): GCC contracts by default
# (-ffp-contract=fast) and ignores the STDC FP_CONTRACT pragma, so the float order of
# _fused_eval.c is pinned by the flag (clang already honours the file's pragma). In _rasterize.c
# the flag also removes clang's default in-expression contraction (on arm64: the multiply-adds of
# the Sutherland-Hodgman intersections), so the clip uses the plain double order of the Python
# reference. Validated on macOS clang only.
_NO_CONTRACT = [] if sys.platform == "win32" else ["-ffp-contract=off"]

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
        extra_compile_args=_NO_CONTRACT,
    ),
]

setup(ext_modules=ext_modules)
