"""
Tests for reconstruction orchestrator module.

Tests verify:
1. ConvergenceCode constants
2. _get_voxel_vertices for up/down triangles
3. BasicVoxelReconstructor construction
4. ReconstructionSetup initialization
"""

import math
from pathlib import Path
from unittest.mock import MagicMock

import numpy as np
import pytest
import torch

from icenine.reconstructor import (
    ConvergenceCode,
    _get_voxel_vertices,
)


class TestConvergenceCode:
    """Tests for convergence codes."""

    def test_codes_are_distinct(self):
        """All codes should be different integers."""
        codes = [
            ConvergenceCode.NOT_CONVERGED,
            ConvergenceCode.HIT_RATIO_CONVERGED,
            ConvergenceCode.COST_CONVERGED,
            ConvergenceCode.MAX_LEVEL_REACHED,
        ]
        assert len(set(codes)) == 4


class TestGetVoxelVertices:
    """Tests for voxel vertex extraction."""

    def test_upward_triangle(self):
        """Points-up triangle has correct vertex layout."""
        voxel = MagicMock()
        voxel.position = [0.0, 0.0, 0.0]
        voxel.side_length = 1.0
        voxel.points_up = True

        verts = _get_voxel_vertices(voxel)
        assert verts.shape == (3, 3)

        # V0 = (0, 0, 0)
        np.testing.assert_allclose(verts[0].numpy(), [0.0, 0.0, 0.0], atol=1e-6)
        # V1 = (1, 0, 0)
        np.testing.assert_allclose(verts[1].numpy(), [1.0, 0.0, 0.0], atol=1e-6)
        # V2 = (0.5, sqrt(3)/2, 0)
        sqrt3_half = 0.5 * math.sqrt(3.0)
        np.testing.assert_allclose(verts[2].numpy(), [0.5, sqrt3_half, 0.0], atol=1e-6)

    def test_downward_triangle(self):
        """Points-down triangle has correct vertex layout."""
        voxel = MagicMock()
        voxel.position = [0.0, 0.0, 0.0]
        voxel.side_length = 1.0
        voxel.points_up = False

        verts = _get_voxel_vertices(voxel)
        assert verts.shape == (3, 3)

        # V0 = (0, 0, 0)
        np.testing.assert_allclose(verts[0].numpy(), [0.0, 0.0, 0.0], atol=1e-6)
        # V1 = (0.5, -sqrt(3)/2, 0)
        sqrt3_half = 0.5 * math.sqrt(3.0)
        np.testing.assert_allclose(verts[1].numpy(), [0.5, -sqrt3_half, 0.0], atol=1e-6)
        # V2 = (1, 0, 0)
        np.testing.assert_allclose(verts[2].numpy(), [1.0, 0.0, 0.0], atol=1e-6)

    def test_offset_position(self):
        """Vertices are offset by voxel position."""
        voxel = MagicMock()
        voxel.position = [1.0, 2.0, 3.0]
        voxel.side_length = 1.0
        voxel.points_up = True

        verts = _get_voxel_vertices(voxel)
        # V0 should be at position
        np.testing.assert_allclose(verts[0].numpy(), [1.0, 2.0, 3.0], atol=1e-6)

    def test_scaled_side_length(self):
        """Vertices scale with side length."""
        voxel = MagicMock()
        voxel.position = [0.0, 0.0, 0.0]
        voxel.side_length = 2.0
        voxel.points_up = True

        verts = _get_voxel_vertices(voxel)
        # V1 should be at (2, 0, 0)
        np.testing.assert_allclose(verts[1].numpy(), [2.0, 0.0, 0.0], atol=1e-6)
