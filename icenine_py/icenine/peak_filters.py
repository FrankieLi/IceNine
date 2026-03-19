"""
Peak acceptance filters for diffraction simulation.

Filters determine which diffraction peaks are observable and calculates their
intensities based on geometric and physical constraints.

Python port of Src/PeakFilters.h

Author: S. F. Li
"""

import torch
from typing import Tuple


class XDMEtaAcceptFn:
    """
    Peak acceptance filter based on eta angle with Lorentz-polarization correction.

    The eta angle (η) describes the angular position of a diffraction peak on
    the Debye-Scherrer ring. This filter rejects peaks outside the acceptable
    eta range and calculates intensity with appropriate geometric corrections.

    C++ Reference:
        PeakFilters.h:60-87 struct XDMEtaAcceptFn

    Attributes:
        min_eta: Minimum acceptable eta angle (radians)
        max_eta: Maximum acceptable eta angle (radians)
        form_intensity: Form factor intensity (structure factor contribution)
        sin_2theta: sin(2θ) where θ is the Bragg angle

    Physics:
        - Eta angle: η = atan2(|Gy|, |Gz|) where G is scattering vector
        - Lorentz-polarization correction: I = I_form / (|sin(η)| × sin(2θ))
        - This accounts for the geometric probability of observing the reflection
    """

    def __init__(
        self,
        min_eta: float,
        max_eta: float,
        form_intensity: float,
        sin_2theta: float
    ):
        """
        Initialize peak acceptance filter.

        Args:
            min_eta: Minimum eta angle (radians)
            max_eta: Maximum eta angle (radians)
            form_intensity: Form factor intensity (arbitrary units)
            sin_2theta: sin(2θ) for this reflection

        Example:
            >>> import numpy as np
            >>> # Accept peaks with eta < 60°, for (111) reflection at θ=15°
            >>> filter_fn = XDMEtaAcceptFn(
            ...     min_eta=0.0,
            ...     max_eta=np.deg2rad(60),
            ...     form_intensity=100.0,
            ...     sin_2theta=np.sin(2 * np.deg2rad(15))
            ... )
        """
        self.min_eta = min_eta
        self.max_eta = max_eta
        self.form_intensity = form_intensity
        self.sin_2theta = sin_2theta

    def __call__(
        self,
        scattering_dir: torch.Tensor
    ) -> Tuple[bool, float]:
        """
        Evaluate filter for given scattering direction.

        Args:
            scattering_dir: Scattering vector direction in sample frame, shape (3,)
                           G = [Gx, Gy, Gz] (should be normalized)

        Returns:
            Tuple of (accept, intensity):
                - accept: True if peak passes eta angle filter
                - intensity: Corrected intensity (valid only if accept=True)

        Algorithm:
            1. Calculate eta: η = atan2(|Gy|, |Gz|)
            2. Check if min_eta < η < max_eta
            3. Calculate intensity: I = I_form / (|sin(η)| × sin(2θ))

        C++ Reference:
            PeakFilters.h:69-83 operator()

        Example:
            >>> import torch
            >>> filter_fn = XDMEtaAcceptFn(0.0, 1.0, 100.0, 0.5)
            >>>
            >>> # Scattering along +Z (eta = 0)
            >>> accept, intensity = filter_fn(torch.tensor([0., 0., 1.]))
            >>> accept
            True
            >>>
            >>> # Scattering along +Y (eta = π/2)
            >>> accept, intensity = filter_fn(torch.tensor([0., 1., 0.]))
            >>> accept
            False  # eta exceeds max_eta
        """
        # Extract components (C++: oScatteringDir.m_fY, oScatteringDir.m_fZ)
        gy = scattering_dir[1]
        gz = scattering_dir[2]

        # Calculate eta angle: atan2(|Gy|, |Gz|)
        # C++: Float fEta = atan2( fabs( oScatteringDir.m_fY ),
        #                          fabs( oScatteringDir.m_fZ ) );
        eta = torch.atan2(torch.abs(gy), torch.abs(gz))

        # Check acceptance
        # C++: Bool bAccept = ( fEta < fMaxEta );
        # Note: C++ only checks max, assumes min=0. We check both for generality.
        accept = (eta >= self.min_eta) and (eta < self.max_eta)

        # Calculate intensity with Lorentz-polarization correction
        # C++: Float fIntensity = fFormIntensity / ( fabs(sin(fEta)) * fSin2Theta );
        sin_eta = torch.sin(eta)
        intensity = self.form_intensity / (torch.abs(sin_eta) * self.sin_2theta + 1e-10)

        return accept.item() if isinstance(accept, torch.Tensor) else accept, intensity.item()

    def __repr__(self):
        """String representation for debugging."""
        return (
            f"XDMEtaAcceptFn("
            f"min_eta={self.min_eta:.4f}, "
            f"max_eta={self.max_eta:.4f}, "
            f"form_intensity={self.form_intensity:.2f}, "
            f"sin_2theta={self.sin_2theta:.4f})"
        )


class TrivialAcceptFn:
    """
    Trivial acceptance filter that accepts all peaks with constant intensity.

    Useful for testing and debugging when geometric effects should be ignored.

    Example:
        >>> filter_fn = TrivialAcceptFn(intensity=1.0)
        >>> accept, intensity = filter_fn(torch.tensor([0., 0., 1.]))
        >>> accept
        True
        >>> intensity
        1.0
    """

    def __init__(self, intensity: float = 1.0):
        """
        Initialize trivial filter.

        Args:
            intensity: Constant intensity for all peaks
        """
        self.intensity = intensity

    def __call__(self, scattering_dir: torch.Tensor) -> Tuple[bool, float]:
        """
        Accept all peaks with constant intensity.

        Args:
            scattering_dir: Scattering direction (ignored)

        Returns:
            (True, intensity)
        """
        return True, self.intensity

    def __repr__(self):
        """String representation for debugging."""
        return f"TrivialAcceptFn(intensity={self.intensity})"
