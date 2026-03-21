"""
Tests for cost functions module.

Tests verify:
1. OverlapInfo incremental quality averaging
2. Qualified peak counting (contiguous detector logic)
3. VoxelCostFunction integration with forward sim output
"""

from pathlib import Path

import numpy as np
import pytest
import torch

from icenine.cost_functions import (
    OverlapInfo,
    count_qualified_peaks,
)


# ============================================================================
# Test: OverlapInfo
# ============================================================================

class TestOverlapInfo:
    """Tests for OverlapInfo dataclass."""

    def test_initial_state(self):
        """All counters start at zero."""
        info = OverlapInfo()
        assert info.pixel_overlap == 0
        assert info.quality == 0.0
        assert info.cost == 1.0
        assert info.hit_ratio == 0.0

    def test_update_counts(self):
        """Counts accumulate correctly."""
        info = OverlapInfo()
        info.update_counts(10, 20, 1, 1)
        assert info.pixel_overlap == 10
        assert info.pixel_on_detector == 20
        assert info.peak_overlap == 1
        assert info.peak_on_detector == 1

        info.update_counts(5, 10, 1, 1)
        assert info.pixel_overlap == 15
        assert info.pixel_on_detector == 30
        assert info.peak_overlap == 2

    def test_hit_ratio(self):
        """Hit ratio = pixel_overlap / pixel_on_detector."""
        info = OverlapInfo()
        info.update_counts(10, 20, 1, 1)
        assert abs(info.hit_ratio - 0.5) < 1e-10

    def test_update_quality_single(self):
        """Quality for single peak with full overlap on 2 detectors."""
        info = OverlapInfo()
        info.update_quality(
            peak_pixel_overlap=10, peak_pixel_on_detector=10,
            peak_detectors_overlap=2, n_detectors_total=2,
        )

        # (10/10) * (2/2) = 1.0
        assert abs(info.quality - 1.0) < 1e-10
        assert abs(info.cost - 0.0) < 1e-10

    def test_update_quality_incremental(self):
        """Quality uses incremental (Welford) running mean with per-peak values."""
        info = OverlapInfo()

        # First peak: 100% overlap, 2/2 detectors
        info.update_quality(
            peak_pixel_overlap=10, peak_pixel_on_detector=10,
            peak_detectors_overlap=2, n_detectors_total=2,
        )
        assert abs(info.quality - 1.0) < 1e-10

        # Second peak: 5/10 pixel overlap, 1/2 detectors
        info.update_quality(
            peak_pixel_overlap=5, peak_pixel_on_detector=10,
            peak_detectors_overlap=1, n_detectors_total=2,
        )

        # cur_quality = (5/10) * (1/2) = 0.25
        # running_mean = 1.0 + (0.25 - 1.0) / 2 = 0.625
        assert abs(info.quality - 0.625) < 1e-10

    def test_update_quality_zero_pixels(self):
        """Quality unchanged when no pixels on detector."""
        info = OverlapInfo()
        info.update_quality(
            peak_pixel_overlap=0, peak_pixel_on_detector=0,
            peak_detectors_overlap=0, n_detectors_total=2,
        )
        assert info.quality == 0.0


# ============================================================================
# Test: Qualified Peak Counting
# ============================================================================

class TestCountQualifiedPeaks:
    """Tests for contiguous detector validation."""

    def test_single_detector_lit_and_overlap(self):
        """Single detector: lit + overlap → valid."""
        peak_on, peak_ovlp, n_det = count_qualified_peaks(
            [True], [True], 1
        )
        assert peak_on == 1
        assert peak_ovlp == 1
        assert n_det == 1

    def test_single_detector_lit_no_overlap(self):
        """Single detector: lit, no overlap → peak valid but overlap invalid."""
        peak_on, peak_ovlp, n_det = count_qualified_peaks(
            [True], [False], 1
        )
        assert peak_on == 1
        assert peak_ovlp == 0
        assert n_det == 0

    def test_two_detectors_contiguous(self):
        """Two contiguous detectors with overlap."""
        peak_on, peak_ovlp, n_det = count_qualified_peaks(
            [True, True], [True, True], 2
        )
        assert peak_on == 1
        assert peak_ovlp == 1
        assert n_det == 2

    def test_two_detectors_first_only(self):
        """Only first detector lit → valid (trailing unlit OK)."""
        peak_on, peak_ovlp, n_det = count_qualified_peaks(
            [True, False], [True, False], 2
        )
        assert peak_on == 1
        assert peak_ovlp == 1
        assert n_det == 1

    def test_gap_in_detectors_invalid(self):
        """Gap in lit detectors → invalid peak."""
        peak_on, peak_ovlp, n_det = count_qualified_peaks(
            [True, False, True], [True, False, True], 3
        )
        assert peak_on == 0
        assert peak_ovlp == 0
        assert n_det == 0

    def test_not_starting_at_zero_invalid(self):
        """Lit pattern not starting at detector 0 → invalid."""
        peak_on, peak_ovlp, n_det = count_qualified_peaks(
            [False, True, True], [False, True, True], 3
        )
        assert peak_on == 0
        assert peak_ovlp == 0

    def test_empty_detectors(self):
        """No detectors → all zeros."""
        peak_on, peak_ovlp, n_det = count_qualified_peaks(
            [], [], 0
        )
        assert peak_on == 0
        assert peak_ovlp == 0
        assert n_det == 0

    def test_all_dark(self):
        """All detectors dark → invalid."""
        peak_on, peak_ovlp, n_det = count_qualified_peaks(
            [False, False], [False, False], 2
        )
        assert peak_on == 0
        assert peak_ovlp == 0
