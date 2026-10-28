"""Setup script for C extensions (numpy include dirs)."""
from setuptools import setup, Extension
import numpy as np

ext_modules = [
    Extension(
        "icenine._rasterize",
        # _fused_eval.c: the fused evaluate of the FusedCost key (S2); ICENINE_FUSED adds its
        # method table entries and the HAS_FUSED constant in _rasterize.c
        sources=["icenine/_rasterize.c", "icenine/_fused_eval.c"],
        include_dirs=[np.get_include()],
        define_macros=[("ICENINE_FUSED", "1")],
    ),
]

setup(ext_modules=ext_modules)
