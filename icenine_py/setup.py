"""Setup script for C extensions (numpy include dirs)."""
from setuptools import setup, Extension
import numpy as np

ext_modules = [
    Extension(
        "icenine._rasterize",
        sources=["icenine/_rasterize.c"],
        include_dirs=[np.get_include()],
    ),
]

setup(ext_modules=ext_modules)
