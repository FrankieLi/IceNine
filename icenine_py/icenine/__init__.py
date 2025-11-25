"""
IceNine Python/PyTorch Implementation

A Python port of IceNine crystallography and diffraction primitives,
focused on forward model calculations for X-ray diffraction.

Modules:
    constants: Physical constants (KEV_OVER_HBAR_C_IN_ANG, etc.)
    symmetry: Crystal symmetry operations using pymatgen
    crystal_structure: Crystal structure definitions and reciprocal lattice
    diffraction_core: Core diffraction calculations (scattering vectors, omega angles)
    geometry: Geometric primitives (Euler angles, Plane, Ray) for transformations
    detector: Detector geometry and coordinate transformations
    image_data: Detector image container with dual-mode storage (dense/sparse)
    mic_file: MIC file I/O for microstructure voxel data
    simulation_range: Omega range system for discontinuous data collection
    sample: Sample with voxel grid and coordinate transformations
    config_file: Configuration file parser for IceNine experiments
    file_io: File I/O utilities for detector and structure files
    experiment_setup: Experiment setup and parameter management
"""

__version__ = "0.1.0"
__author__ = "S. F. Li"

from . import constants
from . import symmetry
from . import crystal_structure
from . import diffraction_core
from . import geometry
from . import detector
from . import image_data
from . import mic_file
from . import simulation_range
from . import sample
from . import config_file
from . import file_io
from . import experiment_setup

__all__ = [
    "constants",
    "symmetry",
    "crystal_structure",
    "diffraction_core",
    "geometry",
    "detector",
    "image_data",
    "mic_file",
    "simulation_range",
    "sample",
    "config_file",
    "file_io",
    "experiment_setup",
]
