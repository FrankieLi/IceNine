"""
Detector image data container with dual-mode storage support.

This module provides ImageData, a container for X-ray detector images with two storage modes:
1. Dense mode: Full differentiability for vision models and gradient-based optimization
2. Sparse mode: Memory-efficient storage for large-scale reconstruction

The API is unified - all operations work identically regardless of storage mode.

Author: S. F. Li
Date: 2025-01-16

C++ Reference: Src/ImageData.h, Src/ImageData.cpp
"""

from typing import Optional, Tuple, Union
from dataclasses import dataclass
import torch
import numpy as np

try:
    from ._rasterize import triangle_overlap as _c_triangle_overlap
    from ._rasterize import pixel_radius_overlap as _c_pixel_radius_overlap
    _HAS_C_RASTERIZE = True
except ImportError:
    _HAS_C_RASTERIZE = False


@dataclass
class ImageDataParameters:
    """
    Parameters for ImageData configuration.

    Attributes:
        num_rows: Number of pixel rows (K direction)
        num_cols: Number of pixel columns (J direction)
        mode: Storage mode ('dense' or 'sparse')
        dtype: PyTorch data type
        device: PyTorch device ('cpu' or 'cuda')
    """
    num_rows: int
    num_cols: int
    mode: str  # 'dense' or 'sparse'
    dtype: torch.dtype
    device: str


class ImageData:
    """
    Detector image container with dual-mode storage.

    Storage Modes:
    --------------
    - 'dense': torch.Tensor (num_rows, num_cols)
      * Full differentiability through all operations
      * Higher memory usage
      * Optimal for gradient-based optimization and vision models
      * Memory: num_rows × num_cols × sizeof(dtype)

    - 'sparse': Sparse coordinate format
      * Memory-efficient (only stores non-zero pixels)
      * Limited differentiability (PyTorch sparse tensors support some gradients)
      * Optimal for large-scale reconstruction with sparse detector data
      * Memory: ~num_nonzero × (2 × int64 + sizeof(dtype))

    Coordinate System:
    ------------------
    - J-axis: Horizontal (columns), 0 to num_cols-1
    - K-axis: Vertical (rows), 0 to num_rows-1
    - Indexing: pixels[k, j] (row-major, consistent with NumPy/PyTorch)

    Usage:
    ------
    # Dense mode (differentiable)
    image = ImageData(2048, 2048, mode='dense')
    image.add_triangle(v0, v1, v2, intensity=1.0, mode='soft')

    # Sparse mode (memory efficient)
    image = ImageData(2048, 2048, mode='sparse')
    image.set_pixel(100, 200, 1.5)

    # Convert between modes
    dense_image = sparse_image.to_dense()
    sparse_image = dense_image.to_sparse()

    C++ Reference:
    --------------
    Src/ImageData.h class CImageData
    Src/ImageData.cpp
    """

    def __init__(
        self,
        num_rows: int,
        num_cols: int,
        mode: str = 'dense',
        dtype: torch.dtype = torch.float32,
        device: str = 'cpu'
    ):
        """
        Initialize ImageData.

        Args:
            num_rows: Number of pixel rows (K direction)
            num_cols: Number of pixel columns (J direction)
            mode: Storage mode - 'dense' or 'sparse'
            dtype: PyTorch data type for pixel values
            device: PyTorch device ('cpu' or 'cuda')

        Raises:
            ValueError: If mode is not 'dense' or 'sparse'
        """
        if mode not in ['dense', 'sparse']:
            raise ValueError(f"mode must be 'dense' or 'sparse', got '{mode}'")

        self.num_rows = num_rows
        self.num_cols = num_cols
        self._mode = mode
        self.dtype = dtype
        self.device = torch.device(device)

        # Storage: only one is active depending on mode
        self._pixels_dense: Optional[torch.Tensor] = None
        self._pixels_sparse: Optional[torch.sparse.Tensor] = None

        # Binary cache for reconstruction (uint8, lazily populated)
        self._binary_cache: Optional[np.ndarray] = None

        # Initialize storage based on mode
        if mode == 'dense':
            self._pixels_dense = torch.zeros(
                (num_rows, num_cols),
                dtype=dtype,
                device=self.device
            )
        else:  # sparse
            # Create empty sparse tensor in COO format
            indices = torch.empty((2, 0), dtype=torch.long, device=self.device)
            values = torch.empty(0, dtype=dtype, device=self.device)
            self._pixels_sparse = torch.sparse_coo_tensor(
                indices,
                values,
                (num_rows, num_cols),
                dtype=dtype,
                device=self.device
            )

    # =========================================================================
    # Properties
    # =========================================================================

    @property
    def mode(self) -> str:
        """Get current storage mode ('dense' or 'sparse')."""
        return self._mode

    @property
    def shape(self) -> Tuple[int, int]:
        """Get image shape (num_rows, num_cols)."""
        return (self.num_rows, self.num_cols)

    @property
    def num_nonzero(self) -> int:
        """Count number of non-zero pixels."""
        if self._mode == 'dense':
            return int((self._pixels_dense != 0).sum().item())
        else:
            # Sparse tensor automatically stores only non-zeros
            return self._pixels_sparse._nnz()

    @property
    def density(self) -> float:
        """Fraction of non-zero pixels (0.0 to 1.0)."""
        total_pixels = self.num_rows * self.num_cols
        return self.num_nonzero / total_pixels if total_pixels > 0 else 0.0

    # =========================================================================
    # Binary Cache (for reconstruction cost functions)
    # =========================================================================

    def get_binary_numpy(self) -> np.ndarray:
        """
        Return a C-contiguous uint8 array where 1 = pixel > 0, 0 = otherwise.

        Lazily computed and cached. The cache is invalidated by any mutator
        method (set_pixel, add_to_pixel, clear, set_pixels, load_ascii,
        add_triangle, add_triangle_scanline).

        Used by reconstruction cost functions which only need binary overlap
        checks, not intensity values. uint8 is 4x smaller than float32.
        """
        if self._binary_cache is None:
            self.ensure_binary_cache()
        return self._binary_cache

    def ensure_binary_cache(self) -> None:
        """Pre-compute and cache the binary numpy array."""
        if self._mode == 'dense':
            self._binary_cache = (self._pixels_dense.numpy() > 0).astype(np.uint8)
        else:
            self._binary_cache = (self._pixels_sparse.to_dense().numpy() > 0).astype(np.uint8)
        # Guarantee C-contiguous for C extension
        if not self._binary_cache.flags['C_CONTIGUOUS']:
            self._binary_cache = np.ascontiguousarray(self._binary_cache)

    # =========================================================================
    # Basic Pixel Operations
    # =========================================================================

    def set_pixel(
        self,
        j: Union[int, torch.Tensor],
        k: Union[int, torch.Tensor],
        value: Union[float, torch.Tensor]
    ) -> None:
        """
        Set pixel value at position (j, k).

        Args:
            j: Column index (J-axis), scalar or tensor
            k: Row index (K-axis), scalar or tensor
            value: Pixel value to set

        C++ Reference:
            ImageData.cpp CImageData::SetPixel
        """
        self._binary_cache = None
        j = self._to_tensor(j, dtype=torch.long)
        k = self._to_tensor(k, dtype=torch.long)
        value = self._to_tensor(value, dtype=self.dtype)

        if self._mode == 'dense':
            self._pixels_dense[k, j] = value
        else:
            # For sparse: need to rebuild tensor with new value
            # This is less efficient for single-pixel updates but maintains sparsity
            self._set_pixel_sparse(j, k, value)

    def get_pixel(
        self,
        j: Union[int, torch.Tensor],
        k: Union[int, torch.Tensor]
    ) -> torch.Tensor:
        """
        Get pixel value at position (j, k).

        Args:
            j: Column index (J-axis), scalar or tensor
            k: Row index (K-axis), scalar or tensor

        Returns:
            Pixel value as tensor

        C++ Reference:
            ImageData.cpp CImageData::At
        """
        j = self._to_tensor(j, dtype=torch.long)
        k = self._to_tensor(k, dtype=torch.long)

        if self._mode == 'dense':
            return self._pixels_dense[k, j]
        else:
            # Convert to dense temporarily for indexing
            # TODO: Optimize for batched queries
            return self._pixels_sparse.to_dense()[k, j]

    def add_to_pixel(
        self,
        j: Union[int, torch.Tensor],
        k: Union[int, torch.Tensor],
        value: Union[float, torch.Tensor]
    ) -> None:
        """
        Add value to pixel at position (j, k) (accumulate).

        Args:
            j: Column index (J-axis), scalar or tensor
            k: Row index (K-axis), scalar or tensor
            value: Value to add to pixel

        C++ Reference:
            ImageData.cpp CImageData::AddToPixel
        """
        self._binary_cache = None
        j = self._to_tensor(j, dtype=torch.long)
        k = self._to_tensor(k, dtype=torch.long)
        value = self._to_tensor(value, dtype=self.dtype)

        if self._mode == 'dense':
            self._pixels_dense[k, j] += value
        else:
            current = self.get_pixel(j, k)
            self.set_pixel(j, k, current + value)

    def clear(self) -> None:
        """
        Reset all pixels to zero.

        C++ Reference:
            ImageData.cpp CImageData::ClearImage
        """
        self._binary_cache = None
        if self._mode == 'dense':
            self._pixels_dense.zero_()
        else:
            # Create new empty sparse tensor
            indices = torch.empty((2, 0), dtype=torch.long, device=self.device)
            values = torch.empty(0, dtype=self.dtype, device=self.device)
            self._pixels_sparse = torch.sparse_coo_tensor(
                indices,
                values,
                (self.num_rows, self.num_cols),
                dtype=self.dtype,
                device=self.device
            )

    # =========================================================================
    # Batched Operations
    # =========================================================================

    def set_pixels(
        self,
        j_coords: torch.Tensor,
        k_coords: torch.Tensor,
        values: torch.Tensor
    ) -> None:
        """
        Set multiple pixels at once (batched operation).

        Args:
            j_coords: Column indices, shape (N,)
            k_coords: Row indices, shape (N,)
            values: Pixel values, shape (N,)
        """
        self._binary_cache = None
        if self._mode == 'dense':
            self._pixels_dense[k_coords, j_coords] = values
        else:
            # Batch update for sparse mode
            for j, k, v in zip(j_coords, k_coords, values):
                self._set_pixel_sparse(j, k, v)

    def get_pixels(
        self,
        j_coords: torch.Tensor,
        k_coords: torch.Tensor
    ) -> torch.Tensor:
        """
        Get multiple pixel values at once (batched operation).

        Args:
            j_coords: Column indices, shape (N,)
            k_coords: Row indices, shape (N,)

        Returns:
            Pixel values, shape (N,)
        """
        if self._mode == 'dense':
            return self._pixels_dense[k_coords, j_coords]
        else:
            dense = self._pixels_sparse.to_dense()
            return dense[k_coords, j_coords]

    # =========================================================================
    # Bounds Checking
    # =========================================================================

    def is_in_bounds(
        self,
        j: Union[int, torch.Tensor],
        k: Union[int, torch.Tensor]
    ) -> torch.Tensor:
        """
        Check if pixel coordinates are within image bounds.

        Args:
            j: Column index (J-axis), scalar or tensor
            k: Row index (K-axis), scalar or tensor

        Returns:
            Boolean tensor: True if in bounds, False otherwise

        C++ Reference:
            ImageData.cpp CImageData::IsInBound
        """
        j = self._to_tensor(j)
        k = self._to_tensor(k)

        in_j = (j >= 0) & (j < self.num_cols)
        in_k = (k >= 0) & (k < self.num_rows)
        return in_j & in_k

    def is_dark(
        self,
        j: Union[int, torch.Tensor],
        k: Union[int, torch.Tensor]
    ) -> torch.Tensor:
        """
        Check if pixel is dark (value <= 0).

        Args:
            j: Column index (J-axis), scalar or tensor
            k: Row index (K-axis), scalar or tensor

        Returns:
            Boolean tensor: True if pixel <= 0

        C++ Reference:
            ImageData.cpp CImageData::IsDark
        """
        return self.get_pixel(j, k) <= 0

    def is_bright(
        self,
        j: Union[int, torch.Tensor],
        k: Union[int, torch.Tensor]
    ) -> torch.Tensor:
        """
        Check if pixel is bright (value > 0).

        Args:
            j: Column index (J-axis), scalar or tensor
            k: Row index (K-axis), scalar or tensor

        Returns:
            Boolean tensor: True if pixel > 0

        C++ Reference:
            ImageData.cpp CImageData::IsBright
        """
        return self.get_pixel(j, k) > 0

    # =========================================================================
    # Mode Conversion
    # =========================================================================

    def to_dense(self) -> 'ImageData':
        """
        Convert to dense mode (returns new ImageData instance).

        If already dense, returns a copy.

        Returns:
            New ImageData in dense mode with same pixel data
        """
        if self._mode == 'dense':
            # Return a copy
            new_image = ImageData(
                self.num_rows, self.num_cols,
                mode='dense', dtype=self.dtype, device=self.device.type
            )
            new_image._pixels_dense = self._pixels_dense.clone()
            return new_image
        else:
            # Convert sparse to dense
            new_image = ImageData(
                self.num_rows, self.num_cols,
                mode='dense', dtype=self.dtype, device=self.device.type
            )
            new_image._pixels_dense = self._pixels_sparse.to_dense()
            return new_image

    def to_sparse(self) -> 'ImageData':
        """
        Convert to sparse mode (returns new ImageData instance).

        If already sparse, returns a copy.

        Returns:
            New ImageData in sparse mode with same pixel data
        """
        if self._mode == 'sparse':
            # Return a copy
            new_image = ImageData(
                self.num_rows, self.num_cols,
                mode='sparse', dtype=self.dtype, device=self.device.type
            )
            new_image._pixels_sparse = self._pixels_sparse.clone()
            return new_image
        else:
            # Convert dense to sparse
            new_image = ImageData(
                self.num_rows, self.num_cols,
                mode='sparse', dtype=self.dtype, device=self.device.type
            )
            new_image._pixels_sparse = self._pixels_dense.to_sparse()
            return new_image

    # =========================================================================
    # I/O Operations
    # =========================================================================

    def save_ascii(self, filename: str) -> None:
        """
        Save image to ASCII file (j, k, intensity format).

        Only non-zero pixels are saved (efficient for sparse data).
        Format matches C++ CImageData::WritePixelsToASCII.

        Args:
            filename: Output file path

        File Format:
            # Comment lines start with #
            j1, k1, intensity1
            j2, k2, intensity2
            ...

        C++ Reference:
            ImageData.cpp CImageData::WritePixelsToASCII
        """
        with open(filename, 'w') as f:
            f.write(f"# ImageData ASCII format\n")
            f.write(f"# num_rows={self.num_rows}, num_cols={self.num_cols}\n")
            f.write(f"# Format: j, k, intensity\n")

            if self._mode == 'dense':
                # Write non-zero pixels
                nonzero = torch.nonzero(self._pixels_dense, as_tuple=False)
                for idx in nonzero:
                    k, j = idx[0].item(), idx[1].item()
                    intensity = self._pixels_dense[k, j].item()
                    f.write(f"{j}, {k}, {intensity}\n")
            else:
                # Sparse: write all stored values
                indices = self._pixels_sparse._indices()
                values = self._pixels_sparse._values()
                for i in range(indices.shape[1]):
                    k = indices[0, i].item()
                    j = indices[1, i].item()
                    intensity = values[i].item()
                    f.write(f"{j}, {k}, {intensity}\n")

    def load_ascii(self, filename: str) -> None:
        """
        Load image from ASCII file (j, k, intensity format).

        Args:
            filename: Input file path

        C++ Reference:
            ImageData.cpp CImageData::ReadCXDMSimulationDataFile
        """
        self._binary_cache = None
        j_list, k_list, intensity_list = [], [], []

        with open(filename, 'r') as f:
            for line in f:
                line = line.strip()
                # Skip comments and empty lines
                if not line or line.startswith('#') or line.startswith(','):
                    continue

                parts = line.split(',')
                if len(parts) >= 3:
                    j = int(parts[0].strip())
                    k = int(parts[1].strip())
                    intensity = float(parts[2].strip())
                    j_list.append(j)
                    k_list.append(k)
                    intensity_list.append(intensity)

        # Set pixels
        if j_list:
            j_tensor = torch.tensor(j_list, dtype=torch.long, device=self.device)
            k_tensor = torch.tensor(k_list, dtype=torch.long, device=self.device)
            values = torch.tensor(intensity_list, dtype=self.dtype, device=self.device)

            if self._mode == 'dense':
                self._pixels_dense[k_tensor, j_tensor] = values
            else:
                # Build sparse tensor
                indices = torch.stack([k_tensor, j_tensor], dim=0)
                self._pixels_sparse = torch.sparse_coo_tensor(
                    indices,
                    values,
                    (self.num_rows, self.num_cols),
                    dtype=self.dtype,
                    device=self.device
                ).coalesce()

    def save_binary(self, filename: str) -> None:
        """
        Save image to binary file using PyTorch serialization.

        Saves mode and tensor data. More efficient than ASCII.

        Args:
            filename: Output file path
        """
        data = {
            'num_rows': self.num_rows,
            'num_cols': self.num_cols,
            'mode': self._mode,
            'dtype': str(self.dtype),
            'pixels_dense': self._pixels_dense if self._mode == 'dense' else None,
            'pixels_sparse': self._pixels_sparse if self._mode == 'sparse' else None,
        }
        torch.save(data, filename)

    @staticmethod
    def load_binary(filename: str, device: str = 'cpu', weights_only: bool = False) -> 'ImageData':
        """
        Load image from binary file.

        Args:
            filename: Input file path
            device: Device to load tensors onto
            weights_only: If True, only load weights (PyTorch 2.6+ security feature)

        Returns:
            ImageData instance
        """
        data = torch.load(filename, map_location=device, weights_only=weights_only)

        image = ImageData(
            num_rows=data['num_rows'],
            num_cols=data['num_cols'],
            mode=data['mode'],
            dtype=eval(data['dtype']),  # Convert string back to torch.dtype
            device=device
        )

        if data['mode'] == 'dense':
            image._pixels_dense = data['pixels_dense'].to(device)
        else:
            image._pixels_sparse = data['pixels_sparse'].to(device)

        return image

    # =========================================================================
    # NumPy Conversions
    # =========================================================================

    def to_numpy(self) -> np.ndarray:
        """
        Convert to dense NumPy array.

        Returns:
            NumPy array, shape (num_rows, num_cols)
        """
        if self._mode == 'dense':
            return self._pixels_dense.cpu().numpy()
        else:
            return self._pixels_sparse.to_dense().cpu().numpy()

    @staticmethod
    def from_numpy(
        array: np.ndarray,
        mode: str = 'dense',
        dtype: torch.dtype = torch.float32,
        device: str = 'cpu'
    ) -> 'ImageData':
        """
        Create ImageData from NumPy array.

        Args:
            array: NumPy array, shape (num_rows, num_cols)
            mode: Storage mode ('dense' or 'sparse')
            dtype: PyTorch data type
            device: PyTorch device

        Returns:
            ImageData instance
        """
        num_rows, num_cols = array.shape
        image = ImageData(num_rows, num_cols, mode=mode, dtype=dtype, device=device)

        tensor = torch.from_numpy(array).to(dtype=dtype, device=device)

        if mode == 'dense':
            image._pixels_dense = tensor
        else:
            image._pixels_sparse = tensor.to_sparse()

        return image

    # =========================================================================
    # Helper Methods
    # =========================================================================

    def _to_tensor(
        self,
        value: Union[int, float, torch.Tensor],
        dtype: Optional[torch.dtype] = None
    ) -> torch.Tensor:
        """Convert value to tensor if not already."""
        if isinstance(value, torch.Tensor):
            if dtype is not None and value.dtype != dtype:
                return value.to(dtype=dtype, device=self.device)
            return value.to(device=self.device)
        else:
            if dtype is None:
                dtype = self.dtype
            return torch.tensor(value, dtype=dtype, device=self.device)

    def _set_pixel_sparse(
        self,
        j: torch.Tensor,
        k: torch.Tensor,
        value: torch.Tensor
    ) -> None:
        """
        Set pixel in sparse mode (internal helper).

        Rebuilds sparse tensor with new value.
        """
        # Get existing indices and values
        indices = self._pixels_sparse._indices()
        values = self._pixels_sparse._values()

        # Check if pixel already exists
        mask = (indices[0] == k) & (indices[1] == j)

        if mask.any():
            # Update existing value
            values = values.clone()
            values[mask] = value
        else:
            # Add new entry
            new_idx = torch.tensor([[k.item()], [j.item()]], dtype=torch.long, device=self.device)
            indices = torch.cat([indices, new_idx], dim=1)
            values = torch.cat([values, value.unsqueeze(0)])

        # Rebuild sparse tensor and coalesce
        self._pixels_sparse = torch.sparse_coo_tensor(
            indices,
            values,
            (self.num_rows, self.num_cols),
            dtype=self.dtype,
            device=self.device
        ).coalesce()

    def get_parameters(self) -> ImageDataParameters:
        """
        Get image parameters as dataclass.

        Returns:
            ImageDataParameters with current configuration
        """
        return ImageDataParameters(
            num_rows=self.num_rows,
            num_cols=self.num_cols,
            mode=self._mode,
            dtype=self.dtype,
            device=self.device.type
        )

    # =========================================================================
    # Geometric Rasterization
    # =========================================================================

    def add_triangle(
        self,
        v0: torch.Tensor,
        v1: torch.Tensor,
        v2: torch.Tensor,
        intensity: Union[float, torch.Tensor] = 1.0,
        mode: str = 'soft',
        temperature: float = 1.0
    ) -> None:
        """
        Rasterize triangle onto image.

        This is a key operation for forward simulation: projects voxel faces
        (triangles) onto the detector image.

        Modes:
        ------
        - 'soft': Differentiable soft assignment based on barycentric coordinates
          Enables gradient flow through triangle vertices for optimization.
          Uses sigmoid-based soft thresholding controlled by temperature.

        - 'hard': Discrete binary rasterization (inside/outside test)
          Faster but non-differentiable. Uses strict barycentric coordinate test.

        Args:
            v0, v1, v2: Triangle vertices in pixel coordinates, shape (2,)
                       Format: [j, k] where j=column, k=row
            intensity: Intensity value for triangle pixels (scalar or tensor)
            mode: 'soft' (differentiable) or 'hard' (discrete)
            temperature: Softness parameter for 'soft' mode (lower = sharper)
                        - 0.1: Sharp boundaries, close to hard rasterization
                        - 1.0: Smooth gradients (default)
                        - 10.0: Very soft, wide gradients

        Example:
            >>> image = ImageData(100, 100, mode='dense')
            >>> v0 = torch.tensor([10.0, 10.0])
            >>> v1 = torch.tensor([50.0, 10.0])
            >>> v2 = torch.tensor([30.0, 50.0])
            >>> image.add_triangle(v0, v1, v2, intensity=1.0, mode='soft')

        C++ Reference:
            ImageData.cpp CImageData::AddTriangle
        """
        self._binary_cache = None
        v0 = self._to_tensor(v0, dtype=torch.float32)
        v1 = self._to_tensor(v1, dtype=torch.float32)
        v2 = self._to_tensor(v2, dtype=torch.float32)
        intensity = self._to_tensor(intensity, dtype=self.dtype)

        # Compute bounding box
        j_min = torch.floor(torch.min(torch.stack([v0[0], v1[0], v2[0]]))).long()
        j_max = torch.ceil(torch.max(torch.stack([v0[0], v1[0], v2[0]]))).long()
        k_min = torch.floor(torch.min(torch.stack([v0[1], v1[1], v2[1]]))).long()
        k_max = torch.ceil(torch.max(torch.stack([v0[1], v1[1], v2[1]]))).long()

        # Clamp to image bounds
        j_min = torch.clamp(j_min, 0, self.num_cols - 1)
        j_max = torch.clamp(j_max, 0, self.num_cols - 1)
        k_min = torch.clamp(k_min, 0, self.num_rows - 1)
        k_max = torch.clamp(k_max, 0, self.num_rows - 1)

        # Create grid of pixel centers
        j_range = torch.arange(j_min, j_max + 1, dtype=torch.float32, device=self.device)
        k_range = torch.arange(k_min, k_max + 1, dtype=torch.float32, device=self.device)
        j_grid, k_grid = torch.meshgrid(j_range, k_range, indexing='xy')
        j_grid = j_grid.flatten()
        k_grid = k_grid.flatten()

        # Pixel centers (add 0.5 for center of pixel)
        pixel_points = torch.stack([j_grid + 0.5, k_grid + 0.5], dim=1)

        # Compute barycentric coordinates for all pixels
        bary = self._barycentric_coordinates(pixel_points, v0, v1, v2)

        if mode == 'hard':
            # Hard rasterization: binary inside/outside test
            # Pixel is inside if all barycentric coordinates >= 0
            inside = (bary >= 0).all(dim=1)
            weights = inside.to(dtype=self.dtype)
        else:  # soft
            # Soft rasterization: sigmoid-based soft assignment
            # Use minimum barycentric coordinate as "signed distance"
            # Positive = inside, negative = outside
            min_bary = bary.min(dim=1)[0]

            # Soft thresholding with temperature
            weights = torch.sigmoid(min_bary / temperature)

        # Add weighted intensity to pixels
        if self._mode == 'dense':
            j_idx = j_grid.long()
            k_idx = k_grid.long()
            self._pixels_dense[k_idx, j_idx] += weights * intensity
        else:
            # For sparse mode: add non-zero pixels
            nonzero = weights > 0
            if nonzero.any():
                j_idx = j_grid[nonzero].long()
                k_idx = k_grid[nonzero].long()
                values = weights[nonzero] * intensity
                for j, k, v in zip(j_idx, k_idx, values):
                    self.add_to_pixel(j, k, v)

    def add_triangle_scanline(
        self,
        v0: torch.Tensor,
        v1: torch.Tensor,
        v2: torch.Tensor,
        intensity: Union[float, torch.Tensor] = 1.0,
    ) -> None:
        """
        Rasterize triangle using scanline fill matching C++ exactly.

        Replicates the C++ rasterization pipeline:
        1. Sutherland-Hodgman polygon clipping against detector bounds (float)
        2. round() to convert clipped vertices to integer pixels
        3. Bresenham edge tracing to build edge table
        4. Scanline fill between left/right edges

        This produces identical pixel sets to the C++ implementation.

        Args:
            v0, v1, v2: Triangle vertices in pixel coordinates, shape (2,)
                        Format: [j, k] where j=column, k=row (float)
            intensity: Intensity value for all triangle pixels

        C++ Reference:
            Raster.tmpl.cpp GeneralRasterizePolygon + CalculateScanline
            SutherlandHodgman.h (Sutherland-Hodgman polygon clipping)
        """
        v0 = self._to_tensor(v0, dtype=torch.float32)
        v1 = self._to_tensor(v1, dtype=torch.float32)
        v2 = self._to_tensor(v2, dtype=torch.float32)
        self._binary_cache = None
        intensity_val = float(intensity) if isinstance(intensity, (int, float)) else intensity.item()

        # Step 0: Truncate pixel coordinates to integers, matching C++
        # C++ ToRowPixel/ToColPixel return Int via (Int) cast (truncation toward zero).
        # Negative values become -1 in C++. We replicate this behavior so that
        # the vertices entering the SH clipper have the same integer values as C++.
        def truncate_pixel(val):
            """Match C++ ToRowPixel/ToColPixel: truncate positive, -1 for negative."""
            if val < 0:
                return -1.0
            return float(int(val))

        # Step 1: Sutherland-Hodgman clipping against [0, W-1] x [0, H-1]
        # Vertices are integer-valued floats (matching C++ flow)
        polygon = [
            (truncate_pixel(v0[0].item()), truncate_pixel(v0[1].item())),
            (truncate_pixel(v1[0].item()), truncate_pixel(v1[1].item())),
            (truncate_pixel(v2[0].item()), truncate_pixel(v2[1].item())),
        ]
        clipped = self._sutherland_hodgman_clip(
            polygon, 0.0, float(self.num_cols - 1), 0.0, float(self.num_rows - 1)
        )

        if len(clipped) < 3:
            return  # Triangle fully clipped away

        # Step 2: Round to integer pixel coordinates (matching C++ OutputStage)
        # Since input vertices are already integer-valued, round() is ~no-op for
        # in-bounds vertices. For clipped vertices at boundaries, round() gives
        # the correct integer.
        int_vertices = [(round(x), round(y)) for x, y in clipped]

        # Step 3+4: Bresenham edge table + scanline fill
        pixels = self._scanline_fill(int_vertices)

        if not pixels:
            return

        # Step 5: Accumulate intensity
        for col, row in pixels:
            if 0 <= col < self.num_cols and 0 <= row < self.num_rows:
                if self._mode == 'dense':
                    self._pixels_dense[row, col] += intensity_val
                else:
                    self.add_to_pixel(col, row, intensity_val)

    @staticmethod
    def _sutherland_hodgman_clip(polygon, x_min, x_max, y_min, y_max):
        """
        Sutherland-Hodgman polygon clipping against rectangular boundary.

        Matches C++ SutherlandHodgman.h exactly, including boundary conventions:
        - Left (x >= x_min): inside
        - Right (x < x_max): inside  (C++ uses less, excludes right boundary)
        - Top (y >= y_min): inside
        - Bottom (y < y_max): inside  (C++ uses less, excludes bottom boundary)

        Args:
            polygon: List of (x, y) tuples
            x_min, x_max, y_min, y_max: Clipping boundaries

        Returns:
            List of (x, y) tuples for clipped polygon
        """
        def clip_edge(vertices, is_inside, intersect):
            if not vertices:
                return []
            output = []
            prev = vertices[-1]
            prev_inside = is_inside(prev)
            for curr in vertices:
                curr_inside = is_inside(curr)
                if curr_inside:
                    if not prev_inside:
                        output.append(intersect(prev, curr))
                    output.append(curr)
                elif prev_inside:
                    output.append(intersect(prev, curr))
                prev = curr
                prev_inside = curr_inside
            return output

        def intersect_left(p0, p1):
            dx = p1[0] - p0[0]
            if abs(dx) < 0.01:
                return (x_min, p0[1])
            slope = (p1[1] - p0[1]) / dx
            return (x_min, p0[1] + slope * (x_min - p0[0]))

        def intersect_right(p0, p1):
            dx = p1[0] - p0[0]
            if abs(dx) < 0.01:
                return (x_max, p0[1])
            slope = (p1[1] - p0[1]) / dx
            return (x_max, p0[1] + slope * (x_max - p0[0]))

        def intersect_top(p0, p1):
            dy = p1[1] - p0[1]
            if abs(dy) < 0.01:
                return (p0[0], y_min)
            slope = (p1[0] - p0[0]) / dy
            return (p0[0] + slope * (y_min - p0[1]), y_min)

        def intersect_bottom(p0, p1):
            dy = p1[1] - p0[1]
            if abs(dy) < 0.01:
                return (p0[0], y_max)
            slope = (p1[0] - p0[0]) / dy
            return (p0[0] + slope * (y_max - p0[1]), y_max)

        # C++ clips in order: Right, Top, Left, Bottom
        # Boundary conventions from SutherlandHodgman.h:
        #   Right: less<REAL> → x < x_max
        #   Top: greater_equal<REAL> → y >= y_min
        #   Left: greater_equal<REAL> → x >= x_min
        #   Bottom: less<REAL> → y < y_max
        result = polygon
        result = clip_edge(result, lambda p: p[0] < x_max, intersect_right)
        result = clip_edge(result, lambda p: p[1] >= y_min, intersect_top)
        result = clip_edge(result, lambda p: p[0] >= x_min, intersect_left)
        result = clip_edge(result, lambda p: p[1] < y_max, intersect_bottom)
        return result

    @staticmethod
    def _scanline_fill(int_vertices):
        """
        Scanline fill using Bresenham edge tracing, matching C++ exactly.

        Args:
            int_vertices: List of (x, y) integer tuples (polygon vertices)

        Returns:
            List of (x, y) integer tuples for all filled pixels

        C++ Reference:
            Raster.tmpl.cpp CalculateScanline + GeneralRasterizePolygon fill loop
        """
        if len(int_vertices) < 3:
            return []

        # Find Y range
        y_min = min(v[1] for v in int_vertices)
        y_max = max(v[1] for v in int_vertices)

        if y_min == y_max:
            # Degenerate: horizontal line
            x_min = min(v[0] for v in int_vertices)
            x_max = max(v[0] for v in int_vertices)
            return [(x, y_min) for x in range(x_min, x_max + 1)]

        # Edge table: for each row, [left_x, right_x], initialized to -1
        edge_table = {}
        for y in range(y_min, y_max + 1):
            edge_table[y] = [-1, -1]

        # Trace edges using Bresenham (matching C++ CalculateScanline)
        n = len(int_vertices)
        for i in range(n):
            v0 = int_vertices[(i - 1) % n] if i > 0 else int_vertices[n - 1]
            v1 = int_vertices[i]
            ImageData._bresenham_edge(edge_table, v0[0], v0[1], v1[0], v1[1])

        # Fill scanlines (matching C++ GeneralRasterizePolygon fill loop)
        pixels = []
        for y in range(y_min, y_max + 1):
            left, right = edge_table[y]
            if left >= 0 or right >= 0:
                if left < 0:
                    # Single pixel
                    pixels.append((right, y))
                elif right < 0:
                    # Single pixel
                    pixels.append((left, y))
                else:
                    if left > right:
                        left, right = right, left
                    for x in range(left, right + 1):
                        pixels.append((x, y))
        return pixels

    @staticmethod
    def _bresenham_edge(edge_table, v0x, v0y, v1x, v1y):
        """
        Bresenham line algorithm for edge table, matching C++ CalculateScanline exactly.

        Records left/right extremes for each scanline row.

        C++ Reference:
            Raster.tmpl.cpp CalculateScanline (Bresenham from Wikipedia)
        """
        steep = abs(v1y - v0y) > abs(v1x - v0x)

        if steep:
            v0x, v0y = v0y, v0x
            v1x, v1y = v1y, v1x

        if v0x > v1x:
            v0x, v1x = v1x, v0x
            v0y, v1y = v1y, v0y

        delta_x = v1x - v0x
        delta_y = abs(v1y - v0y)
        error = delta_x
        y_step = 1 if v0y < v1y else -1
        y = v0y

        for x in range(v0x, v1x + 1):
            if steep:
                # plot(y, x) → edge_table[x] records y
                row = x
                col = y
            else:
                # plot(x, y) → edge_table[y] records x
                row = y
                col = x

            if row in edge_table:
                if edge_table[row][0] < 0:
                    edge_table[row][0] = col
                elif edge_table[row][0] > col:
                    edge_table[row][0] = col

                if edge_table[row][1] < col:
                    edge_table[row][1] = col

            error -= 2 * delta_y
            if error < 0:
                y += y_step
                error += 2 * delta_x

    def add_polygon(
        self,
        vertices: torch.Tensor,
        intensity: Union[float, torch.Tensor] = 1.0,
        mode: str = 'soft',
        temperature: float = 1.0
    ) -> None:
        """
        Rasterize polygon onto image by triangulation.

        Uses fan triangulation from first vertex. Assumes simple polygon
        (no self-intersections, vertices in order).

        Args:
            vertices: Polygon vertices in pixel coordinates, shape (N, 2)
                     Format: [..., [j, k], ...] where j=column, k=row
                     Vertices should be in winding order
            intensity: Intensity value for polygon pixels
            mode: 'soft' (differentiable) or 'hard' (discrete)
            temperature: Softness parameter for 'soft' mode

        Example:
            >>> # Square polygon
            >>> vertices = torch.tensor([[10, 10], [50, 10], [50, 50], [10, 50]])
            >>> image.add_polygon(vertices, intensity=1.0)

        C++ Reference:
            ImageData.cpp CImageData::AddPolygon
        """
        vertices = self._to_tensor(vertices, dtype=torch.float32)

        if vertices.shape[0] < 3:
            raise ValueError("Polygon must have at least 3 vertices")

        # Fan triangulation: split into triangles from first vertex
        # Triangle 1: v0, v1, v2
        # Triangle 2: v0, v2, v3
        # Triangle 3: v0, v3, v4
        # ...
        v0 = vertices[0]
        for i in range(1, vertices.shape[0] - 1):
            v1 = vertices[i]
            v2 = vertices[i + 1]
            self.add_triangle(v0, v1, v2, intensity, mode, temperature)

    def get_num_pixels_lit(
        self,
        v0: torch.Tensor,
        v1: torch.Tensor,
        v2: torch.Tensor,
        mode: str = 'soft',
        temperature: float = 1.0
    ) -> torch.Tensor:
        """
        Count number of pixels lit by triangle.

        Creates temporary image, rasterizes triangle, counts bright pixels.

        Args:
            v0, v1, v2: Triangle vertices in pixel coordinates, shape (2,)
            mode: 'soft' or 'hard' rasterization
            temperature: Softness parameter for 'soft' mode

        Returns:
            Number of pixels with intensity > 0 (tensor scalar)

        C++ Reference:
            ImageData.cpp CImageData::GetNumPixelsLit
        """
        # Create temporary image
        temp = ImageData(
            self.num_rows, self.num_cols,
            mode='dense',  # Always use dense for temporary
            dtype=self.dtype,
            device=self.device.type
        )

        # Rasterize triangle
        temp.add_triangle(v0, v1, v2, intensity=1.0, mode=mode, temperature=temperature)

        # Count bright pixels
        if mode == 'hard':
            return (temp._pixels_dense > 0).sum()
        else:
            # For soft mode, sum weights (fractional count)
            return temp._pixels_dense.sum()

    def get_triangle_overlap(
        self,
        v0: torch.Tensor,
        v1: torch.Tensor,
        v2: torch.Tensor,
        mode: str = 'soft',
        temperature: float = 1.0
    ) -> torch.Tensor:
        """
        Count pixels overlapping between triangle and existing image data.

        Args:
            v0, v1, v2: Triangle vertices in pixel coordinates, shape (2,)
            mode: 'soft' or 'hard' rasterization
            temperature: Softness parameter for 'soft' mode

        Returns:
            Number of overlapping pixels (tensor scalar)

        C++ Reference:
            ImageData.cpp CImageData::GetTriangleOverlap
        """
        num_overlap, _ = self.get_triangle_overlap_property(v0, v1, v2, mode, temperature)
        return num_overlap

    def get_triangle_overlap_property(
        self,
        v0: torch.Tensor,
        v1: torch.Tensor,
        v2: torch.Tensor,
        mode: str = 'soft',
        temperature: float = 1.0
    ) -> Tuple[torch.Tensor, torch.Tensor]:
        """
        Calculate overlap between triangle and existing image data.

        This is the CRITICAL operation for cost functions in reconstruction.
        Compares simulated triangle projection against experimental detector image.

        Returns:
            num_overlap_pixels: Pixels where BOTH triangle AND existing data are bright
                               (intersection of simulated and experimental)
            num_lit_pixels: Total pixels lit by triangle
                           (simulated projection size)

        The ratio num_overlap_pixels / num_lit_pixels measures how well the
        simulated projection matches the experimental data.

        Args:
            v0, v1, v2: Triangle vertices in pixel coordinates, shape (2,)
            mode: 'soft' (differentiable) or 'hard' (discrete)
            temperature: Softness parameter for 'soft' mode

        Returns:
            (num_overlap_pixels, num_lit_pixels): Tuple of tensor scalars

        Example:
            >>> # Experimental image with peaks
            >>> experimental = ImageData(100, 100, mode='dense')
            >>> experimental.add_triangle(...)  # Add experimental peaks
            >>>
            >>> # Simulated triangle projection
            >>> v0 = torch.tensor([10.0, 10.0], requires_grad=True)
            >>> v1 = torch.tensor([50.0, 10.0], requires_grad=True)
            >>> v2 = torch.tensor([30.0, 50.0], requires_grad=True)
            >>>
            >>> overlap, total = experimental.get_triangle_overlap_property(v0, v1, v2)
            >>> match_ratio = overlap / total  # How well does it match?
            >>> loss = -overlap  # Maximize overlap
            >>> loss.backward()  # Gradients flow through v0, v1, v2

        C++ Reference:
            ImageData.cpp CImageData::GetTriangleOverlapProperty
            CostFunctions.h (usage in overlap calculations)
        """
        v0 = self._to_tensor(v0, dtype=torch.float32)
        v1 = self._to_tensor(v1, dtype=torch.float32)
        v2 = self._to_tensor(v2, dtype=torch.float32)

        if mode == 'hard':
            # Fast path: use C extension for rasterization + overlap counting
            if _HAS_C_RASTERIZE and self._mode == 'dense':
                image_np = self._pixels_dense.numpy()
                if not image_np.flags['C_CONTIGUOUS']:
                    image_np = np.ascontiguousarray(image_np)
                n_overlap_int, n_lit_int = _c_triangle_overlap(
                    image_np,
                    v0[0].item(), v0[1].item(),
                    v1[0].item(), v1[1].item(),
                    v2[0].item(), v2[1].item(),
                )
                num_overlap = torch.tensor(
                    float(n_overlap_int), dtype=self.dtype, device=self.device
                )
                num_lit = torch.tensor(
                    float(n_lit_int), dtype=self.dtype, device=self.device
                )
            else:
                # Fallback: Python scanline rasterization matching C++
                # GeneralRasterizePolygon. Handles degenerate triangles correctly.
                def truncate_pixel(val):
                    if val < 0:
                        return -1.0
                    return float(int(val))

                polygon = [
                    (truncate_pixel(v0[0].item()), truncate_pixel(v0[1].item())),
                    (truncate_pixel(v1[0].item()), truncate_pixel(v1[1].item())),
                    (truncate_pixel(v2[0].item()), truncate_pixel(v2[1].item())),
                ]
                clipped = self._sutherland_hodgman_clip(
                    polygon, 0.0, float(self.num_cols - 1), 0.0, float(self.num_rows - 1)
                )

                if len(clipped) < 3:
                    zero = torch.tensor(0.0, dtype=self.dtype, device=self.device)
                    return zero, zero

                int_vertices = [(round(x), round(y)) for x, y in clipped]
                pixels = self._scanline_fill(int_vertices)

                if not pixels:
                    zero = torch.tensor(0.0, dtype=self.dtype, device=self.device)
                    return zero, zero

                num_lit_int = 0
                num_overlap_int = 0
                for col, row in pixels:
                    if 0 <= col < self.num_cols and 0 <= row < self.num_rows:
                        num_lit_int += 1
                        if self._mode == 'dense':
                            if self._pixels_dense[row, col].item() > 0:
                                num_overlap_int += 1
                        else:
                            exp_dense = self._pixels_sparse.to_dense()
                            if exp_dense[row, col].item() > 0:
                                num_overlap_int += 1

                num_overlap = torch.tensor(
                    float(num_overlap_int), dtype=self.dtype, device=self.device
                )
                num_lit = torch.tensor(
                    float(num_lit_int), dtype=self.dtype, device=self.device
                )
        else:
            # Soft mode: barycentric coordinates with sigmoid weighting
            all_j = torch.stack([v0[0], v1[0], v2[0]])
            all_k = torch.stack([v0[1], v1[1], v2[1]])
            j_min = max(int(torch.floor(all_j.min()).item()), 0)
            j_max = min(int(torch.ceil(all_j.max()).item()), self.num_cols - 1)
            k_min = max(int(torch.floor(all_k.min()).item()), 0)
            k_max = min(int(torch.ceil(all_k.max()).item()), self.num_rows - 1)

            if j_min > j_max or k_min > k_max:
                zero = torch.tensor(0.0, dtype=self.dtype, device=self.device)
                return zero, zero

            j_range = torch.arange(j_min, j_max + 1, dtype=torch.float32, device=self.device)
            k_range = torch.arange(k_min, k_max + 1, dtype=torch.float32, device=self.device)
            j_grid, k_grid = torch.meshgrid(j_range, k_range, indexing='xy')
            pixel_points = torch.stack([j_grid.flatten() + 0.5, k_grid.flatten() + 0.5], dim=1)

            bary = self._barycentric_coordinates(pixel_points, v0, v1, v2)
            min_bary = bary.min(dim=1)[0]
            weights = torch.sigmoid(min_bary / temperature)

            num_lit = weights.sum()

            j_idx = j_grid.flatten().long()
            k_idx = k_grid.flatten().long()

            if self._mode == 'dense':
                exp_vals = self._pixels_dense[k_idx, j_idx]
            else:
                exp_dense = self._pixels_sparse.to_dense()
                exp_vals = exp_dense[k_idx, j_idx]

            exp_weights = torch.clamp(exp_vals, 0, 1)
            num_overlap = (weights * exp_weights).sum()

        return num_overlap, num_lit

    # =========================================================================
    # Geometric Helper Methods
    # =========================================================================

    def _barycentric_coordinates(
        self,
        points: torch.Tensor,
        v0: torch.Tensor,
        v1: torch.Tensor,
        v2: torch.Tensor
    ) -> torch.Tensor:
        """
        Compute barycentric coordinates for points with respect to triangle.

        Barycentric coordinates (u, v, w) satisfy:
        - point = u*v0 + v*v1 + w*v2
        - u + v + w = 1
        - Point is inside triangle if u, v, w >= 0

        Args:
            points: Query points, shape (N, 2)
            v0, v1, v2: Triangle vertices, shape (2,)

        Returns:
            Barycentric coordinates, shape (N, 3)
            Each row is [u, v, w] for corresponding point
        """
        # Edge vectors
        v0v1 = v1 - v0
        v0v2 = v2 - v0

        # Compute barycentric coordinates using cross products
        # Area of triangle ABC = 0.5 * |AB × AC|
        # For 2D cross product: [x1, y1] × [x2, y2] = x1*y2 - y1*x2

        # Total triangle area (2x)
        denom = v0v1[0] * v0v2[1] - v0v1[1] * v0v2[0]

        # Handle degenerate triangles
        if torch.abs(denom) < 1e-8:
            # Degenerate triangle: return coordinates that indicate "outside"
            return torch.full((points.shape[0], 3), -1.0, dtype=torch.float32, device=self.device)

        # For each point
        v0p = points - v0.unsqueeze(0)  # shape (N, 2)

        # Compute v and w using cross products
        v = (v0p[:, 0] * v0v2[1] - v0p[:, 1] * v0v2[0]) / denom
        w = (v0v1[0] * v0p[:, 1] - v0v1[1] * v0p[:, 0]) / denom
        u = 1.0 - v - w

        # Stack into (N, 3)
        bary = torch.stack([u, v, w], dim=1)

        return bary

    def __repr__(self) -> str:
        """String representation."""
        return (f"ImageData(num_rows={self.num_rows}, num_cols={self.num_cols}, "
                f"mode='{self._mode}', dtype={self.dtype}, "
                f"nonzero={self.num_nonzero}, density={self.density:.4f})")
