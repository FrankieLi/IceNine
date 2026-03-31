"""
Differentiable cost function infrastructure for gradient-based optimization.

Provides:
- ExperimentalImageStack: Pre-stacked experimental images as a single contiguous
  tensor for GPU-friendly batch access via grid_sample.
- MultiScaleImageStack: Max-pool downsampled image pyramid for coarse-to-fine
  optimization. For binary images, max_pool2d is morphological dilation;
  bilinear grid_sample on the downsampled images provides smooth gradients.
- DifferentiableCostFunction: Replaces Stage D (sequential overlap counting)
  with differentiable bilinear sampling, enabling PyTorch autograd through
  the full orientation → cost pipeline.
"""

import math
from dataclasses import dataclass
from typing import List, Optional, Tuple

import numpy as np
import torch
import torch.nn.functional as F

from .crystal_structure import CrystalStructure
from .detector import Detector
from .diffraction_core import get_scattering_omegas_torch
from .experimental_data import ExperimentalData
from .image_data import ImageData
from .sample import Sample
from .simulation import Simulation
from .simulation_range import SimulationRange


# ---------------------------------------------------------------------------
# ExperimentalImageStack — pre-stacked image tensor
# ---------------------------------------------------------------------------


class ExperimentalImageStack:
    """
    Pre-stacked experimental images as a single contiguous tensor.

    Stores all (omega_interval, detector) images in a flat 4D tensor of shape
    (n_omega * n_det, 1, H, W), suitable for batch indexing and grid_sample.

    Flat index: flat_idx = omega_idx * n_det + det_idx

    Args:
        exp_data: ExperimentalData with 2D list of ImageData objects
        binary: If True (default), binarize images to 0.0/1.0 matching
                the existing get_binary_numpy() reconstruction pipeline.
                If False, preserve original intensity values.
        dtype: Storage dtype. Must be float for grid_sample. Default float32.
    """

    def __init__(
        self,
        exp_data: ExperimentalData,
        binary: bool = True,
        dtype: torch.dtype = torch.float32,
    ):
        self.n_omega = exp_data.n_omega_intervals
        self.n_det = exp_data.n_detectors
        self.binary = binary
        self.dtype = dtype

        # Determine image dimensions from first available image
        first_img = exp_data.get_image(0, 0)
        self.H = first_img.num_rows
        self.W = first_img.num_cols

        # Build flat tensor: (n_omega * n_det, 1, H, W)
        n_total = self.n_omega * self.n_det
        self.images = torch.zeros(n_total, 1, self.H, self.W, dtype=dtype)

        for omega_idx in range(self.n_omega):
            for det_idx in range(self.n_det):
                flat_idx = omega_idx * self.n_det + det_idx
                img = exp_data.get_image(omega_idx, det_idx)
                self._fill_image(flat_idx, img, binary)

    def _fill_image(self, flat_idx: int, img: ImageData, binary: bool) -> None:
        """Extract dense tensor from ImageData and store in the stack."""
        if img.mode == "dense":
            tensor = img._pixels_dense
        else:
            tensor = img._pixels_sparse.to_dense()

        if binary:
            tensor = (tensor > 0).to(self.dtype)
        else:
            tensor = tensor.to(self.dtype)

        self.images[flat_idx, 0] = tensor

    @property
    def device(self) -> torch.device:
        return self.images.device

    def to(self, device: torch.device) -> "ExperimentalImageStack":
        """Move image stack to device (GPU/CPU). Returns self for chaining."""
        self.images = self.images.to(device)
        return self

    def flat_index(self, omega_idx: int, det_idx: int) -> int:
        """Compute flat index from (omega, detector) pair."""
        return omega_idx * self.n_det + det_idx

    def get_image(self, omega_idx: int, det_idx: int) -> torch.Tensor:
        """Get single image as (1, H, W) tensor."""
        flat_idx = self.flat_index(omega_idx, det_idx)
        return self.images[flat_idx]

    def get_images_batch(self, flat_indices: torch.Tensor) -> torch.Tensor:
        """
        Batch-gather images by flat index.

        Args:
            flat_indices: (M,) int tensor of flat image indices

        Returns:
            (M, 1, H, W) tensor of images
        """
        return self.images[flat_indices]

    @property
    def memory_bytes(self) -> int:
        """Total memory in bytes."""
        return self.images.nelement() * self.images.element_size()

    def __repr__(self) -> str:
        mb = self.memory_bytes / (1024 * 1024)
        return (
            f"ExperimentalImageStack("
            f"n_omega={self.n_omega}, n_det={self.n_det}, "
            f"H={self.H}, W={self.W}, "
            f"binary={self.binary}, dtype={self.dtype}, "
            f"device={self.device}, "
            f"memory={mb:.1f} MB)"
        )


# ---------------------------------------------------------------------------
# SparseImageStack — memory-efficient sparse coordinate storage
# ---------------------------------------------------------------------------


class SparseImageStack:
    """
    Memory-efficient image stack using sparse coordinate storage.

    For binary diffraction images with <0.001% fill rate (typical: ~9 bright
    pixels per 2048×2048 image), stores only the (row, col) coordinates of
    bright pixels instead of dense float32 tensors. This reduces storage from
    ~5.6 GB to ~26 KB for a typical 180×2 detector dataset.

    Dense float32 images are materialized on-demand for grid_sample via
    get_images_batch(), then freed after use.

    Args:
        exp_data: ExperimentalData with 2D list of ImageData objects
        binary: If True (default), store only pixel coordinates (values=1.0).
                If False, store coordinates and intensity values.
        dtype: Output dtype for densified images. Default float32.
    """

    def __init__(
        self,
        exp_data: ExperimentalData,
        binary: bool = True,
        dtype: torch.dtype = torch.float32,
    ):
        self.n_omega = exp_data.n_omega_intervals
        self.n_det = exp_data.n_detectors
        self.binary = binary
        self.dtype = dtype

        first_img = exp_data.get_image(0, 0)
        self.H = first_img.num_rows
        self.W = first_img.num_cols

        n_total = self.n_omega * self.n_det
        self._pixel_coords: List[torch.Tensor] = []  # (nnz, 2) int16 per image
        self._pixel_values: List[Optional[torch.Tensor]] = []  # (nnz,) or None if binary

        for omega_idx in range(self.n_omega):
            for det_idx in range(self.n_det):
                img = exp_data.get_image(omega_idx, det_idx)
                self._extract_sparse(img, binary)

    def _extract_sparse(self, img: ImageData, binary: bool) -> None:
        """Extract nonzero pixel coordinates (and values) from an ImageData."""
        if img.mode == "dense":
            tensor = img._pixels_dense
        else:
            tensor = img._pixels_sparse.to_dense()

        nz = torch.nonzero(tensor, as_tuple=False)  # (nnz, 2)

        if nz.shape[0] == 0:
            self._pixel_coords.append(torch.empty(0, 2, dtype=torch.int16))
            self._pixel_values.append(None if binary else torch.empty(0, dtype=self.dtype))
        else:
            self._pixel_coords.append(nz.to(torch.int16))
            if binary:
                self._pixel_values.append(None)
            else:
                rows, cols = nz[:, 0].long(), nz[:, 1].long()
                self._pixel_values.append(tensor[rows, cols].to(self.dtype))

    @property
    def device(self) -> torch.device:
        if self._pixel_coords and self._pixel_coords[0].numel() > 0:
            return self._pixel_coords[0].device
        return torch.device("cpu")

    def to(self, device: torch.device) -> "SparseImageStack":
        """Move sparse data to device. Returns self for chaining."""
        self._pixel_coords = [c.to(device) for c in self._pixel_coords]
        self._pixel_values = [
            v.to(device) if v is not None else None for v in self._pixel_values
        ]
        return self

    def flat_index(self, omega_idx: int, det_idx: int) -> int:
        """Compute flat index from (omega, detector) pair."""
        return omega_idx * self.n_det + det_idx

    def get_image(self, omega_idx: int, det_idx: int) -> torch.Tensor:
        """Get single image as (1, H, W) tensor (densified on demand)."""
        flat_idx = self.flat_index(omega_idx, det_idx)
        return self._densify(torch.tensor([flat_idx]))[0]

    def get_images_batch(self, flat_indices: torch.Tensor) -> torch.Tensor:
        """
        Materialize dense float32 images for the requested indices.

        Deduplicates internally: if many peaks share the same omega wedge,
        each unique image is densified only once.

        Args:
            flat_indices: (M,) int tensor of flat image indices

        Returns:
            (M, 1, H, W) float32 tensor
        """
        unique_idx, inverse = torch.unique(flat_indices, return_inverse=True)
        dense_unique = self._densify(unique_idx)  # (U, 1, H, W)
        return dense_unique[inverse]  # (M, 1, H, W)

    def _densify(self, flat_indices: torch.Tensor) -> torch.Tensor:
        """Scatter sparse pixels into dense float32 tensors."""
        U = flat_indices.shape[0]
        result = torch.zeros(U, 1, self.H, self.W, dtype=self.dtype)
        for i, idx in enumerate(flat_indices.tolist()):
            coords = self._pixel_coords[idx]
            if coords.numel() == 0:
                continue
            rows = coords[:, 0].long()
            cols = coords[:, 1].long()
            if self._pixel_values[idx] is not None:
                result[i, 0, rows, cols] = self._pixel_values[idx]
            else:
                result[i, 0, rows, cols] = 1.0
        return result

    @property
    def memory_bytes(self) -> int:
        """Total memory of sparse coordinate storage."""
        total = 0
        for coords in self._pixel_coords:
            total += coords.nelement() * coords.element_size()
        for values in self._pixel_values:
            if values is not None:
                total += values.nelement() * values.element_size()
        return total

    @classmethod
    def from_image_directory(
        cls,
        directory: str,
        basename: str,
        ext: str,
        serial_length: int,
        n_omega: int,
        n_detectors: int,
        num_rows: int,
        num_cols: int,
        binary: bool = True,
        file_start: int = 0,
        det_offset: int = 0,
        dtype: torch.dtype = torch.float32,
    ) -> "SparseImageStack":
        """
        Load directly from ASCII image files into sparse coordinate storage.

        Reads one file at a time, extracts nonzero pixel coordinates, and
        discards the raw data immediately. Peak memory is O(one image) instead
        of O(all images), making this suitable for large datasets that would
        OOM with ExperimentalData.from_image_directory().

        Args:
            directory: Directory containing the image files
            basename: File basename (e.g., "3Grains.sim")
            ext: File extension (e.g., "d")
            serial_length: Digits in serial number (e.g., 5 for "00000")
            n_omega: Number of omega intervals
            n_detectors: Number of detectors
            num_rows: Detector image height in pixels
            num_cols: Detector image width in pixels
            binary: If True (default), store only pixel coordinates (values=1.0).
            file_start: Starting file number (default 0)
            det_offset: Detector numbering offset (default 0)
            dtype: Output dtype for densified images. Default float32.

        Returns:
            SparseImageStack loaded directly from files
        """
        from pathlib import Path
        import time as _time

        directory_path = Path(directory)
        obj = cls.__new__(cls)
        obj.n_omega = n_omega
        obj.n_det = n_detectors
        obj.H = num_rows
        obj.W = num_cols
        obj.binary = binary
        obj.dtype = dtype
        obj._pixel_coords = []
        obj._pixel_values = []

        total_files = n_omega * n_detectors
        loaded = 0
        t_start = _time.time()

        # Iterate omega-first to match flat_index = omega_idx * n_det + det_idx
        for omega_idx in range(n_omega):
            for det_idx in range(n_detectors):
                file_num = file_start + omega_idx
                file_num_str = str(file_num).zfill(serial_length)
                filename = f"{basename}{file_num_str}.{ext}{det_offset + det_idx}"
                filepath = directory_path / filename

                if not filepath.exists():
                    raise FileNotFoundError(f"Image file not found: {filepath}")

                # Parse ASCII file directly into sparse coords (no dense allocation)
                rows_list, cols_list, vals_list = [], [], []
                with open(filepath, "r") as f:
                    for line in f:
                        line = line.strip()
                        if not line or line.startswith("#") or line.startswith(","):
                            continue
                        parts = line.split(",")
                        if len(parts) >= 3:
                            j = int(parts[0].strip())
                            k = int(parts[1].strip())
                            intensity = float(parts[2].strip())
                            if intensity > 0:
                                rows_list.append(k)
                                cols_list.append(j)
                                vals_list.append(intensity)

                if rows_list:
                    coords = torch.tensor(
                        list(zip(rows_list, cols_list)), dtype=torch.int16
                    )
                    obj._pixel_coords.append(coords)
                    if binary:
                        obj._pixel_values.append(None)
                    else:
                        obj._pixel_values.append(
                            torch.tensor(vals_list, dtype=dtype)
                        )
                else:
                    obj._pixel_coords.append(torch.empty(0, 2, dtype=torch.int16))
                    obj._pixel_values.append(
                        None if binary else torch.empty(0, dtype=dtype)
                    )

                loaded += 1
                if loaded % 60 == 0 or loaded == total_files:
                    elapsed = _time.time() - t_start
                    print(
                        f"  Loading sparse images: {loaded}/{total_files} "
                        f"({elapsed:.1f}s)",
                        flush=True,
                    )

        return obj

    def __repr__(self) -> str:
        total_nnz = sum(c.shape[0] for c in self._pixel_coords)
        kb = self.memory_bytes / 1024
        return (
            f"SparseImageStack("
            f"n_omega={self.n_omega}, n_det={self.n_det}, "
            f"H={self.H}, W={self.W}, "
            f"binary={self.binary}, "
            f"total_bright_pixels={total_nnz}, "
            f"memory={kb:.1f} KB)"
        )


# ---------------------------------------------------------------------------
# MultiScaleImageStack — max_pool downsampling pyramid
# ---------------------------------------------------------------------------


class MultiScaleImageStack:
    """
    Max-pool downsampled image pyramid for coarse-to-fine optimization.

    **Why this broadens the cost function (same goal as Gaussian blur):**

    The IceNine cost function is extremely sharp: quality drops near zero within
    ~0.3° of the correct orientation because Bragg peaks are only ~3 pixels wide.
    A gradient optimizer starting >0.3° away sees zero gradient and cannot converge.

    The original plan used Gaussian blur to spread each ~3px spot to ~23px,
    widening the angular basin from ~0.3° to ~2°. However, Gaussian blur via
    F.conv2d causes a memory explosion (~13 GB per 2048×2048 image due to im2col
    unfolding) and is unnecessary for binary images.

    **Why max_pool downsampling achieves the same effect:**

    1. **Morphological dilation**: For binary images (0/1 pixels), max_pool2d with
       kernel size k is exactly morphological dilation by a k×k square structuring
       element, followed by downsampling. Any bright pixel in the k×k neighborhood
       survives. A factor-8 downsample on a 3px spot produces a ~1px spot in the
       256×256 output — but that 1px now represents an 8px-wide region in the
       original space, so any simulated peak landing within 4px of the real peak
       will "hit" it.

    2. **Smooth gradients via bilinear grid_sample**: The cost function evaluates
       soft overlap by sampling the experimental image at the simulated peak's
       centroid with bilinear interpolation. On a downsampled image, the bilinear
       kernel spans ~8× more angular space, giving a smooth gradient signal even
       when the simulated peak is several pixels away in the original image.

    3. **No F.conv2d needed**: max_pool2d is O(H*W) in memory (sliding window),
       not O(H*W*kernel²) like im2col convolution. Peak memory per image drops
       from ~13 GB (Gaussian, 31×31 kernel) to ~14 MB (max_pool, factor=8).

    Each scale is produced by max_pool2d at a power-of-2 downsample factor.

    **Omega-direction broadening (optional):**

    In addition to spatial downsampling, `omega_window > 0` applies morphological
    dilation along the omega axis: each downsampled frame is replaced by the
    element-wise max over ±`omega_window` neighboring frames.  This extends the
    effective per-frame coverage from 1° to (2*omega_window+1)°, complementing the
    spatial dilation from max_pool2d.

    Omega blending is only applied to downsampled scales (factor > 1).  At full
    resolution the images are stored sparsely — densifying all 360 frames to blend
    them would cost ~5.6 GB and provides little benefit since the spatial basin is
    already very sharp at that scale.

    Boundary semantics: zero-padding (no circular wrap-around).  Frame 0 does NOT
    pick up signal from frame n_omega-1, which is physically correct because the
    beginning and end of the omega range are not adjacent.

    Args:
        image_stack: Base SparseImageStack or ExperimentalImageStack (scale 0)
        downsample_factors: Power-of-2 factors. Default [1, 4, 8].
                            factor=1: original resolution (no broadening).
                            factor=4: 512×512, ~4px effective spot width.
                            factor=8: 256×256, ~8px effective spot width (~2° basin).
        omega_window: Number of neighboring omega frames to include on each side
                      via element-wise max (morphological dilation in omega direction).
                      Default 0 (disabled).  omega_window=1 gives ±1 frame = 3-frame
                      coverage, extending the omega basin from ~1° to ~3°.
    """

    def __init__(
        self,
        image_stack,
        downsample_factors: Optional[List[int]] = None,
        omega_window: int = 0,
        _prebuilt_downsampled: Optional[List] = None,
    ):
        """
        Args:
            image_stack: Base SparseImageStack (scale 0, shared by reference).
            downsample_factors: Downsampling factors; default [1, 4, 8].
            omega_window: Omega blending radius (0 = disabled).
            _prebuilt_downsampled: Optional list of already-densified
                ExperimentalImageStack objects for factor > 1 scales, in the
                same order as downsample_factors (excluding the factor=1 entry).
                When provided, _downsample_stack is skipped and the supplied
                stacks are used directly (omega blending still applied if
                omega_window > 0).  Use MultiScaleImageStack.build_shared_base()
                to build a single set of downsampled stacks and then construct
                multiple MultiScaleImageStack objects with different omega_windows
                without re-densifying.
        """
        if downsample_factors is None:
            downsample_factors = [1, 4, 8]

        self.downsample_factors = downsample_factors
        self.omega_window = omega_window
        self.n_scales = len(downsample_factors)
        self.scales: List[ExperimentalImageStack] = []

        prebuilt_iter = iter(_prebuilt_downsampled or [])
        for factor in downsample_factors:
            if factor == 1:
                # Scale 0: share original stack (no copy), never omega-blend
                self.scales.append(image_stack)
            else:
                if _prebuilt_downsampled is not None:
                    ds = next(prebuilt_iter)
                else:
                    ds = self._downsample_stack(image_stack, factor)
                if omega_window > 0:
                    ds = self._omega_blend(ds, omega_window)
                self.scales.append(ds)

    @classmethod
    def build_shared_base(
        cls,
        image_stack,
        downsample_factors: Optional[List[int]] = None,
    ) -> List:
        """Build the downsampled (but unblended) stacks once for sharing.

        Returns a list of ExperimentalImageStack objects (one per factor > 1
        entry in downsample_factors) that can be passed to multiple
        MultiScaleImageStack constructors via _prebuilt_downsampled=.

        Usage:
            shared = MultiScaleImageStack.build_shared_base(sparse, [1, 4, 8])
            ms0 = MultiScaleImageStack(sparse, [1,4,8], omega_window=0,
                                       _prebuilt_downsampled=shared)
            ms1 = MultiScaleImageStack(sparse, [1,4,8], omega_window=1,
                                       _prebuilt_downsampled=shared)
            ms2 = MultiScaleImageStack(sparse, [1,4,8], omega_window=2,
                                       _prebuilt_downsampled=shared)
        """
        if downsample_factors is None:
            downsample_factors = [1, 4, 8]
        tmp = cls.__new__(cls)
        tmp.downsample_factors = downsample_factors
        tmp.omega_window = 0
        tmp.n_scales = len(downsample_factors)
        tmp.scales = []
        result = []
        for factor in downsample_factors:
            if factor > 1:
                ds = tmp._downsample_stack(image_stack, factor)
                result.append(ds)
        return result

    def _downsample_stack(self, stack, factor: int) -> ExperimentalImageStack:
        """Densify and max_pool2d-downsample all images. No Gaussian blur.

        max_pool2d with kernel=factor on binary images is morphological dilation
        by a (factor×factor) structuring element, then downsampling.
        """
        n_total = stack.n_omega * stack.n_det
        out_H = stack.H // factor
        out_W = stack.W // factor
        is_sparse = isinstance(stack, SparseImageStack)

        downsampled = torch.zeros(n_total, 1, out_H, out_W, dtype=torch.float32)

        # Process one image at a time to minimize peak memory.
        # IMPORTANT: reuse a single full-resolution buffer rather than allocating
        # a new (H×W) tensor each iteration. Without this, 360 separate 16 MB
        # allocations accumulate in PyTorch's memory pool even after being freed,
        # resulting in ~6 GB of pooled memory for a 2048² image stack.
        with torch.no_grad():
            buf = torch.zeros(1, 1, stack.H, stack.W, dtype=torch.float32)
            for i in range(n_total):
                buf.zero_()
                if is_sparse:
                    coords = stack._pixel_coords[i]
                    if coords.numel() > 0:
                        rows = coords[:, 0].long()
                        cols = coords[:, 1].long()
                        vals = stack._pixel_values[i]
                        buf[0, 0, rows, cols] = 1.0 if vals is None else vals
                else:
                    buf[0:1].copy_(stack.images[i : i + 1])
                downsampled[i : i + 1] = F.max_pool2d(buf, factor)

        result = ExperimentalImageStack.__new__(ExperimentalImageStack)
        result.n_omega = stack.n_omega
        result.n_det = stack.n_det
        result.H = out_H
        result.W = out_W
        result.binary = True  # max_pool on binary input stays binary
        result.dtype = torch.float32
        result.images = downsampled
        return result

    def _omega_blend(self, stack: ExperimentalImageStack, window: int) -> ExperimentalImageStack:
        """Morphological dilation along the omega axis (max over ±window frames).

        For binary images this is equivalent to a 1D flat structuring element of
        width 2*window+1 in the omega direction, applied before each frame is
        evaluated.  Zero-padding semantics: out-of-range frames contribute nothing.

        Args:
            stack: Dense ExperimentalImageStack produced by _downsample_stack.
            window: Number of frames on each side to include.

        Returns:
            New ExperimentalImageStack with the same shape but omega-blended images.
        """
        n_omega = stack.n_omega
        n_det = stack.n_det
        H = stack.H
        W = stack.W

        # Reshape to (n_omega, n_det, 1, H, W) for clean omega-axis indexing.
        # blended_orig stays unmodified (source for all shifts).
        # blended accumulates the maximum — written in-place via out= to avoid
        # allocating a temporary tensor on each shift iteration.
        with torch.no_grad():
            imgs = stack.images  # (n_omega * n_det, 1, H, W)
            blended_orig = imgs.reshape(n_omega, n_det, 1, H, W).clone()
            blended = blended_orig.clone()

            for shift in range(1, window + 1):
                # Use out= to write result directly into blended, no temp alloc.
                torch.maximum(blended[shift:], blended_orig[:-shift],
                               out=blended[shift:])
                torch.maximum(blended[:-shift], blended_orig[shift:],
                               out=blended[:-shift])

        result = ExperimentalImageStack.__new__(ExperimentalImageStack)
        result.n_omega = n_omega
        result.n_det = n_det
        result.H = H
        result.W = W
        result.binary = True
        result.dtype = torch.float32
        result.images = blended.reshape(n_omega * n_det, 1, H, W)
        return result

    def get_at_scale(self, scale_idx: int) -> ExperimentalImageStack:
        """Get image stack at the given downsample scale."""
        return self.scales[scale_idx]

    def to(self, device: torch.device) -> "MultiScaleImageStack":
        """Move all scales to device. Returns self for chaining."""
        for stack in self.scales:
            stack.to(device)
        return self

    def __repr__(self) -> str:
        base = self.scales[0]
        mb_total = sum(s.memory_bytes for s in self.scales) / (1024 * 1024)
        return (
            f"MultiScaleImageStack("
            f"n_scales={self.n_scales}, downsample_factors={self.downsample_factors}, "
            f"omega_window={self.omega_window}, "
            f"shape=({base.n_omega}×{base.n_det}, {base.H}×{base.W}), "
            f"total_memory={mb_total:.1f} MB)"
        )


# ---------------------------------------------------------------------------
# DifferentiableOverlapInfo — gradient-carrying overlap result
# ---------------------------------------------------------------------------


@dataclass
class DifferentiableOverlapInfo:
    """
    Differentiable overlap metrics — quality and cost carry autograd grad_fn.

    Unlike OverlapInfo (which uses plain floats), these tensors preserve
    the computation graph for backpropagation through the cost function.
    """

    quality: torch.Tensor  # scalar, has grad_fn
    cost: torch.Tensor  # 1 - quality, has grad_fn
    n_peaks: int  # number of peaks evaluated (informational)


# ---------------------------------------------------------------------------
# DifferentiableCostFunction — full pipeline with grid_sample Stage D
# ---------------------------------------------------------------------------


class DifferentiableCostFunction:
    """
    Differentiable cost function for gradient-based orientation optimization.

    Reimplements the VoxelCostFunction pipeline with two key changes:
    1. Accepts torch.Tensor orientation (with requires_grad) instead of np.ndarray
    2. Replaces Stage D (binary overlap counting) with differentiable bilinear
       sampling via grid_sample against pre-stacked image tensors

    The omega-to-wedge routing (Stage A) and eta filtering are treated as
    non-differentiable discrete routing — gradients flow only through the
    pixel coordinate computation (Stages B-C) and the soft overlap (Stage D').

    Usage:
        image_stack = ExperimentalImageStack(exp_data, binary=True)
        multi_stack = MultiScaleImageStack(image_stack)
        diff_cost = DifferentiableCostFunction(
            simulator, detector_list, range_map, multi_stack,
            sample, structure_list, max_q=8.0,
        )

        orientation = torch.tensor(orient_matrix, requires_grad=True)
        info = diff_cost.evaluate(orientation, voxel_vertices, scale=0)
        info.cost.backward()
        # orientation.grad now contains dcost/dorientation
    """

    def __init__(
        self,
        simulator: Simulation,
        detector_list: List[Detector],
        range_map: SimulationRange,
        image_stack: "MultiScaleImageStack",
        sample: Sample,
        structure_list: List[CrystalStructure],
        eta_limit: float = math.pi / 2.0,
        max_q: float = 0.0,
    ):
        self.simulator = simulator
        self.detector_list = detector_list
        self.range_map = range_map
        self.image_stack = image_stack
        self.sample = sample
        self.structure_list = structure_list
        self.eta_limit = eta_limit
        self.eval_count = 0

        # Pre-compute reciprocal vectors per phase (same as VoxelCostFunction)
        self._phase_recip_vecs: dict = {}
        for phase_idx, structure in enumerate(structure_list):
            recp_vecs = structure.get_reflection_vectors()
            if recp_vecs:
                if max_q > 0:
                    recp_vecs = [rv for rv in recp_vecs if rv.q_mag <= max_q]
                if not recp_vecs:
                    continue
                g_hkl_batch = torch.stack(
                    [torch.from_numpy(rv.q_vec).float() for rv in recp_vecs]
                )
                g_mag_batch = torch.tensor(
                    [rv.q_mag for rv in recp_vecs], dtype=torch.float32
                )
                self._phase_recip_vecs[phase_idx] = (g_hkl_batch, g_mag_batch)

    def evaluate(
        self,
        orientation: torch.Tensor,
        voxel_vertices: torch.Tensor,
        phase_index: int = 0,
        scale: int = 0,
    ) -> DifferentiableOverlapInfo:
        """
        Evaluate differentiable cost for a candidate orientation.

        Args:
            orientation: 3x3 rotation matrix as torch.Tensor (can have requires_grad=True)
            voxel_vertices: Triangle vertices in sample frame, shape (3, 3)
            phase_index: Crystal phase index
            scale: Blur scale index for multi-scale image stack

        Returns:
            DifferentiableOverlapInfo with quality/cost tensors carrying grad_fn
        """
        self.eval_count += 1
        zero_result = DifferentiableOverlapInfo(
            quality=torch.tensor(0.0), cost=torch.tensor(1.0), n_peaks=0
        )

        if phase_index not in self._phase_recip_vecs:
            return zero_result

        g_hkl_batch, g_mag_batch = self._phase_recip_vecs[phase_index]
        stack = self.image_stack.get_at_scale(scale)

        # --- Orientation → observable peaks (differentiable through orientation) ---
        # g_lab = orientation @ g_hkl^T  (gradients flow through orientation)
        g_lab_batch = (orientation @ g_hkl_batch.T).T  # (K, 3)

        bragg_result = get_scattering_omegas_torch(
            g_lab_batch,
            g_mag_batch,
            self.simulator.beam_energy,
            self.simulator.beam_deflection_chi,
        )

        obs_mask = bragg_result.observable
        if not obs_mask.any():
            return zero_result

        beam_dir = self.simulator.beam_direction
        base_rot = self.sample.sample_to_lab_matrix[:3, :3]

        obs_g = g_lab_batch[obs_mask]
        obs_mag = torch.norm(obs_g, dim=1, keepdim=True)
        obs_normals = obs_g / obs_mag

        obs_omega1 = bragg_result.omega1[obs_mask]
        obs_omega2 = bragg_result.omega2[obs_mask]
        all_omegas = torch.cat([obs_omega1, obs_omega2])
        all_normals = torch.cat([obs_normals, obs_normals])
        N = all_omegas.shape[0]

        # Batch Rz(omega)
        cos_w = torch.cos(all_omegas)
        sin_w = torch.sin(all_omegas)
        Rz = torch.zeros(N, 3, 3)
        Rz[:, 0, 0] = cos_w
        Rz[:, 0, 1] = -sin_w
        Rz[:, 1, 0] = sin_w
        Rz[:, 1, 1] = cos_w
        Rz[:, 2, 2] = 1.0

        full_rot = Rz @ base_rot
        lab_normals = torch.bmm(full_rot, all_normals.unsqueeze(-1)).squeeze(-1)
        dot = (beam_dir.unsqueeze(0) * lab_normals).sum(dim=1, keepdim=True)
        reflected = beam_dir.unsqueeze(0) - 2.0 * dot * lab_normals

        # Eta filter (non-differentiable routing)
        rd_norms = torch.norm(reflected, dim=1)
        safe_norms = torch.where(rd_norms > 0, rd_norms, torch.ones_like(rd_norms))
        ry = torch.abs(reflected[:, 1]) / safe_norms
        rz_val = torch.abs(reflected[:, 2]) / safe_norms
        eta = torch.atan2(ry, rz_val)
        valid = (eta < self.eta_limit) & (rd_norms > 0)

        peak_omegas = all_omegas[valid]
        peak_normals = all_normals[valid]
        M_total = peak_omegas.shape[0]

        if M_total == 0:
            return zero_result

        # --- Stage A: omega-to-wedge mapping (non-differentiable routing) ---
        low = self.range_map.low
        width = self.range_map.width
        # Detach omegas for routing — no gradients through wedge selection
        omega_np = peak_omegas.detach().numpy()
        bin_indices = ((omega_np - low) / width).astype(int)
        index_list = self.range_map.index_list
        n_bins = len(index_list)

        wedge_indices = np.full(M_total, -1, dtype=np.int64)
        for i in range(M_total):
            n = bin_indices[i]
            if 0 <= n < n_bins and index_list[n] is not None:
                wedge_indices[i] = index_list[n]

        valid_mask = wedge_indices >= 0
        if not np.any(valid_mask):
            return zero_result

        vi = np.where(valid_mask)[0]
        vi_t = torch.from_numpy(vi).long()
        v_omegas = peak_omegas[vi_t]
        v_normals = peak_normals[vi_t]
        v_wedge = wedge_indices[vi]
        M = len(vi)

        # --- Stages B-C: batch geometry + ray-detector intersection ---
        orig_matrix = self.sample.sample_to_lab_matrix

        cos_w2 = torch.cos(v_omegas)
        sin_w2 = torch.sin(v_omegas)
        Rz2 = torch.zeros(M, 3, 3)
        Rz2[:, 0, 0] = cos_w2
        Rz2[:, 0, 1] = -sin_w2
        Rz2[:, 1, 0] = sin_w2
        Rz2[:, 1, 1] = cos_w2
        Rz2[:, 2, 2] = 1.0

        base_rot2 = orig_matrix[:3, :3]
        full_rot2 = Rz2 @ base_rot2
        lab_normals2 = torch.bmm(full_rot2, v_normals.unsqueeze(-1)).squeeze(-1)
        dot2 = (beam_dir.unsqueeze(0) * lab_normals2).sum(dim=1, keepdim=True)
        reflected2 = beam_dir.unsqueeze(0) - 2.0 * dot2 * lab_normals2

        # Build 4x4 transforms for vertex projection
        full_4x4 = torch.zeros(M, 4, 4)
        full_4x4[:, :3, :3] = full_rot2
        full_4x4[:, :3, 3] = orig_matrix[:3, 3]
        full_4x4[:, 3, 3] = 1.0

        verts_4d = torch.cat([voxel_vertices, torch.ones(3, 1)], dim=1)
        lab_verts = torch.einsum("mij,vj->mvi", full_4x4, verts_4d)[:, :, :3]

        # --- Stage D': Differentiable soft overlap via grid_sample ---
        n_det = len(self.detector_list)
        # Use detector pixel dimensions for bounds and normalization —
        # the stack may be downsampled (blurred), but grid_sample's [-1,1]
        # normalized coords handle the mapping automatically.
        det_H = self.detector_list[0].num_rows
        det_W = self.detector_list[0].num_cols

        # Accumulate soft quality across detectors
        quality_accum = torch.tensor(0.0)
        n_valid_samples = 0

        for det_idx in range(n_det):
            detector = self.detector_list[det_idx]
            plane = detector._detector_plane
            plane_n = plane.normal
            plane_d = plane.d

            # Ray-plane intersection for all M peaks × 3 vertices
            origins = lab_verts.reshape(M * 3, 3)
            dirs = reflected2.unsqueeze(1).expand(M, 3, 3).reshape(M * 3, 3)

            denom = (dirs * plane_n).sum(dim=1)
            numer = -((origins * plane_n).sum(dim=1) + plane_d)
            parallel = torch.abs(denom) < 1e-8
            safe_denom = torch.where(parallel, torch.ones_like(denom), denom)
            t = torch.where(parallel, torch.zeros_like(denom), numer / safe_denom)
            hits = (~parallel) & (t > 0)

            pts = origins + t.unsqueeze(1) * dirs

            # Detector coordinate transform (inlined for differentiability)
            relative = pts - detector._position
            pixel_loc = relative - detector._lab_frame_coord_origin
            j = (pixel_loc * detector._lab_frame_basis_j).sum(dim=1)
            k = (pixel_loc * detector._lab_frame_basis_k).sum(dim=1)

            col = (j + detector.pixel_half_width) / detector.pixel_width
            row = (k + detector.pixel_half_height) / detector.pixel_height

            hits_mv = hits.reshape(M, 3)
            all_hit = hits_mv.all(dim=1)  # (M,) bool
            cols_mv = col.reshape(M, 3)
            rows_mv = row.reshape(M, 3)

            # Centroid of triangle vertices (differentiable)
            centroids_col = cols_mv.mean(dim=1)  # (M,)
            centroids_row = rows_mv.mean(dim=1)  # (M,)

            # Bounds check: only count peaks whose centroids are within detector
            # (The hard cost function does this implicitly in Stage D)
            in_bounds = (
                (centroids_col >= 0) & (centroids_col < det_W)
                & (centroids_row >= 0) & (centroids_row < det_H)
            )
            valid_det = all_hit & in_bounds  # (M,) bool

            # Build flat indices for image lookup
            flat_idx = torch.from_numpy(
                (v_wedge * n_det + det_idx).astype(np.int64)
            )

            # Gather images for all M peaks (handles both dense and sparse stacks)
            batch_imgs = stack.get_images_batch(flat_idx)  # (M, 1, H, W)

            # Normalize to [-1, 1] for grid_sample (align_corners=True)
            # pixel 0 → -1, pixel det_W-1 → +1  (works for any stack resolution)
            grid_x = 2.0 * centroids_col / (det_W - 1) - 1.0
            grid_y = 2.0 * centroids_row / (det_H - 1) - 1.0
            grid = torch.stack([grid_x, grid_y], dim=1)  # (M, 2)
            grid = grid.unsqueeze(1).unsqueeze(1)  # (M, 1, 1, 2)

            # Bilinear sampling (differentiable w.r.t. grid coordinates)
            sampled = F.grid_sample(
                batch_imgs, grid, mode="bilinear",
                padding_mode="zeros", align_corners=True,
            )  # (M, 1, 1, 1)
            sampled = sampled.squeeze(3).squeeze(2).squeeze(1)  # (M,)

            # Mask by valid_det (non-differentiable mask, but multiplied differentiably)
            hit_float = valid_det.float()
            masked = sampled * hit_float

            quality_accum = quality_accum + masked.sum()
            n_valid_samples += int(hit_float.sum().item())

        # Aggregate quality
        if n_valid_samples == 0:
            return zero_result

        quality = quality_accum / n_valid_samples
        cost = 1.0 - quality

        return DifferentiableOverlapInfo(
            quality=quality,
            cost=cost,
            n_peaks=M,
        )
