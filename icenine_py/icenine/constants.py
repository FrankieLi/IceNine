"""
Physical constants for X-ray diffraction calculations.

This module contains fundamental physical constants used throughout
the IceNine diffraction calculations, ported from PhysicalConstants.h

All constants use the same units as the C++ implementation:
- Lengths in Angstroms (Å)
- Energies in keV
- Angles in radians
"""

# keV/(ℏc) in 1/Angstrom
# This is chosen for convenience since lattice constants are measured in Angstroms
# Fundamental constant: ℏc ≈ 1.9732705 keV·Å
# Therefore: 1 keV / (ℏc) = 1 / 1.9732705 ≈ 0.506773182 Å⁻¹
KEV_OVER_HBAR_C_IN_ANG = 0.506773182

# Conversion factor: Angstroms to millimeters
ANGSTROMS_TO_MM = 1.0e-7

# For backwards compatibility with C++ naming conventions
keV_over_hbar_c_in_ang = KEV_OVER_HBAR_C_IN_ANG
angstroms_to_mm = ANGSTROMS_TO_MM


def angstrom_to_mm(ang: float) -> float:
    """
    Convert length from Angstroms to millimeters.

    Args:
        ang: Length in Angstroms

    Returns:
        Length in millimeters

    Example:
        >>> angstrom_to_mm(1e7)  # 1 mm in Angstroms
        1.0
    """
    return ang * ANGSTROMS_TO_MM


def wavelength_from_energy(energy_kev: float) -> float:
    """
    Calculate X-ray wavelength from beam energy using E = hc/λ.

    Args:
        energy_kev: X-ray energy in keV

    Returns:
        Wavelength in Angstroms

    Example:
        >>> wavelength_from_energy(50.02099)  # Beam energy from ReconstructTest.config
        0.24780...
    """
    # E = hc/λ  =>  λ = hc/E
    # λ (Å) = (ℏc in keV·Å) / E(keV)
    return 1.0 / (KEV_OVER_HBAR_C_IN_ANG * energy_kev)


def energy_from_wavelength(wavelength_ang: float) -> float:
    """
    Calculate X-ray energy from wavelength using E = hc/λ.

    Args:
        wavelength_ang: X-ray wavelength in Angstroms

    Returns:
        Energy in keV

    Example:
        >>> energy_from_wavelength(0.2478)  # Wavelength from ReconstructTest.config
        50.020...
    """
    # E = hc/λ
    # E (keV) = (ℏc in keV·Å) * (1/λ in Å)
    return 1.0 / (KEV_OVER_HBAR_C_IN_ANG * wavelength_ang)


# Test parameters from ConfigFiles/ReconstructTest.config
TEST_PARAMS = {
    "element": "Au",
    "lattice_type": "FCC",
    "a": 4.0782,  # Angstroms
    "beam_energy": 50.02099,  # keV
    "wavelength": 0.2478,  # Angstroms
    "max_q": 8.0,  # Å⁻¹
    "space_group": 225,  # Fm-3m (face-centered cubic)
}
