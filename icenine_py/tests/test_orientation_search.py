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


class TestRiemannianAdamOptimizer:
    """Tests for RiemannianAdamOptimizer."""

    def test_requires_geoopt_error(self, monkeypatch):
        """RiemannianAdamOptimizer raises RuntimeError when geoopt is unavailable."""
        import icenine.orientation_search as os_mod
        monkeypatch.setattr(os_mod, "_GEOOPT_AVAILABLE", False)

        from icenine.orientation_search import RiemannianAdamOptimizer

        with pytest.raises(RuntimeError, match="geoopt is required"):
            RiemannianAdamOptimizer(
                hard_cost_fn=None,
                diff_cost_fn=None,
                voxel_vertices=None,
            )

    @pytest.mark.skipif(
        not __import__("importlib").util.find_spec("geoopt"),
        reason="geoopt not installed",
    )
    def test_returns_search_candidate(self):
        """RiemannianAdamOptimizer.optimize() returns SearchCandidate with cost < 1.0."""
        import torch

        from icenine.orientation_search import RiemannianAdamOptimizer, SearchCandidate
        from icenine.cost_functions import OverlapInfo

        # Hard cost mock: always returns cost=0.3
        class MockHardInfo:
            cost = 0.3
            hit_ratio = 0.8
            peak_overlap = 5
            peak_on_detector = 6
            pixel_overlap = 8
            pixel_on_detector = 10

        class MockHardCost:
            def evaluate(self, orientation, voxel_vertices, phase_index):
                return MockHardInfo()

        # Differentiable cost mock: returns constant tensor (no grad)
        class MockDiffInfo:
            cost = torch.tensor(0.3)

        class MockDiffCost:
            def evaluate(self, R, voxel_vertices, phase_index=0, scale=2):
                return MockDiffInfo()

        opt = RiemannianAdamOptimizer(
            hard_cost_fn=MockHardCost(),
            diff_cost_fn=MockDiffCost(),
            voxel_vertices=None,
            phase_index=0,
            rng=np.random.default_rng(42),
        )
        result = opt.optimize(
            initial_orientation=np.eye(3),
            angular_box_side=math.radians(2.0),
            n_steps=5,
            lr=1e-4,
            scale=2,
            max_restarts=1,
        )
        assert isinstance(result, SearchCandidate)
        assert result.cost < 1.0

    @pytest.mark.skipif(
        not __import__("importlib").util.find_spec("geoopt"),
        reason="geoopt not installed",
    )
    def test_zero_restarts_single_hard_eval_after_adam(self):
        """With max_restarts=0, hard_cost_fn.evaluate called exactly twice (initial + post-Adam)."""
        import torch
        from icenine.orientation_search import RiemannianAdamOptimizer

        eval_count = [0]

        class MockHardInfo:
            cost = 0.5
            hit_ratio = 0.0
            peak_overlap = 1
            peak_on_detector = 2
            pixel_overlap = 5
            pixel_on_detector = 10

        class MockHardCost:
            def evaluate(self, orientation, voxel_vertices, phase_index):
                eval_count[0] += 1
                return MockHardInfo()

        class MockDiffInfo:
            cost = torch.tensor(0.5)

        class MockDiffCost:
            def evaluate(self, R, voxel_vertices, phase_index=0, scale=2):
                return MockDiffInfo()

        opt = RiemannianAdamOptimizer(
            hard_cost_fn=MockHardCost(),
            diff_cost_fn=MockDiffCost(),
            voxel_vertices=None,
            rng=np.random.default_rng(0),
        )
        opt.optimize(
            initial_orientation=np.eye(3),
            angular_box_side=math.radians(2.0),
            n_steps=3,
            max_restarts=0,
        )
        # 1 initial eval + 1 post-Adam eval = 2 total
        assert eval_count[0] == 2

    def test_search_params_hybrid_defaults(self):
        """SearchParameters defaults have use_hybrid_optimizer=False."""
        from icenine.orientation_search import SearchParameters

        p = SearchParameters()
        assert p.use_hybrid_optimizer is False
        assert p.adam_n_steps == 100
        assert abs(p.adam_lr - 1e-4) < 1e-12
        assert p.adam_scale == 2
