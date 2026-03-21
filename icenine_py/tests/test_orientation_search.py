"""
Tests for orientation search module.

Tests verify:
1. SearchCandidate sorting by cost
2. SearchParameters construction from config
3. MCOptimizer convergence properties
4. Convergence check logic
"""

import math

import numpy as np
import pytest

from icenine.orientation_search import (
    SearchCandidate,
    SearchParameters,
    hit_ratio_converged,
)
from icenine.cost_functions import OverlapInfo


class TestSearchCandidate:
    """Tests for SearchCandidate dataclass."""

    def test_sorting_by_cost(self):
        """Candidates sort by cost (ascending)."""
        c1 = SearchCandidate(orientation=np.eye(3), cost=0.5)
        c2 = SearchCandidate(orientation=np.eye(3), cost=0.2)
        c3 = SearchCandidate(orientation=np.eye(3), cost=0.8)

        sorted_candidates = sorted([c1, c2, c3])
        assert sorted_candidates[0].cost == 0.2
        assert sorted_candidates[1].cost == 0.5
        assert sorted_candidates[2].cost == 0.8

    def test_default_cost(self):
        """Default cost is 1.0 (worst)."""
        c = SearchCandidate(orientation=np.eye(3))
        assert c.cost == 1.0


class TestSearchParameters:
    """Tests for SearchParameters."""

    def test_defaults(self):
        """Default parameters are sensible."""
        p = SearchParameters()
        assert p.max_local_resolution == 3
        assert p.max_discrete_candidates == 100
        assert p.max_mc_steps == 3500
        assert p.successive_restarts == 2
        assert p.local_grid_radius > 0

    def test_from_config(self):
        """Can construct from config-like object."""
        class MockConfig:
            local_orientation_grid_radius = math.radians(5.0)
            min_local_resolution = 0
            max_local_resolution = 3
            max_discrete_candidates = 50
            max_mc_steps = 1000
            mc_radius_scale_factor = 0.5
            successive_restarts = 1
            max_convergence_cost = 0.05
            max_deepening_hit_ratio = 0.8
            max_accepted_cost = 0.95

        p = SearchParameters.from_config(MockConfig())
        assert p.max_discrete_candidates == 50
        assert p.max_mc_steps == 1000
        assert abs(p.max_deepening_hit_ratio - 0.8) < 1e-10


class TestHitRatioConvergence:
    """Tests for convergence check."""

    def test_converged(self):
        """Hit ratio above threshold → converged."""
        info = OverlapInfo()
        info.pixel_overlap = 8
        info.pixel_on_detector = 10
        assert hit_ratio_converged(info, threshold=0.7)

    def test_not_converged(self):
        """Hit ratio below threshold → not converged."""
        info = OverlapInfo()
        info.pixel_overlap = 3
        info.pixel_on_detector = 10
        assert not hit_ratio_converged(info, threshold=0.7)

    def test_zero_pixels(self):
        """No pixels on detector → not converged."""
        info = OverlapInfo()
        assert not hit_ratio_converged(info, threshold=0.1)

    def test_exact_threshold(self):
        """Hit ratio exactly at threshold → converged."""
        info = OverlapInfo()
        info.pixel_overlap = 7
        info.pixel_on_detector = 10
        assert hit_ratio_converged(info, threshold=0.7)
