"""
Detector geometry and coordinate transformations for X-ray diffraction.

This module implements the Detector class, which represents a 2D area detector
used in High-Energy Diffraction Microscopy (HEDM) experiments. It handles:

- Coordinate transformations between lab frame, detector frame, and pixel coordinates
- Detector positioning and orientation in 3D space
- Ray-detector intersection calculations
- Differentiable operations for integration with PyTorch-based simulators

Coordinate Systems:
    Lab Frame: 3D global coordinate system
        - X-axis: beam direction (beam travels along +X)
        - Y-axis: horizontal (perpendicular to beam)
        - Z-axis: vertical (perpendicular to beam)
        - Origin: typically at sample center

    Detector Frame: 2D coordinate system on detector surface (in mm)
        - J-axis: horizontal direction on detector (mm)
        - K-axis: vertical direction on detector (mm)
        - Origin: defined by beam center and coordinate origin offset

    Pixel Coordinates: 2D discrete grid indices
        - Column (col): horizontal pixel index (0 to num_cols-1)
        - Row (row): vertical pixel index (0 to num_rows-1)
        - Origin: (0, 0) at top-left corner

Convention:
    - All angles in degrees (user-facing API)
    - Euler angles: Bunge convention (ZXZ intrinsic rotations)
    - Distances in mm
    - Detector plane normal points towards the source (+X direction when unrotated)

C++ Reference:
    Src/Detector.h
    Src/Detector.cpp
    XDM++/libXDM/3dMath.h (Euler angle conventions)
"""

import torch
from typing import Tuple, Optional, Union
from dataclasses import dataclass

from .geometry import euler_to_matrix_torch, Plane, Ray


@dataclass
class DetectorParameters:
    """
    Parameters for detector geometry.

    This dataclass holds all geometric parameters needed to define a detector.
    It's useful for saving/loading detector configurations.

    Attributes:
        num_rows: Number of pixel rows (K direction)
        num_cols: Number of pixel columns (J direction)
        pixel_height: Physical height of each pixel in mm (K direction)
        pixel_width: Physical width of each pixel in mm (J direction)
        beam_center_j: Beam center position in PIXELS along J-axis
        beam_center_k: Beam center position in PIXELS along K-axis
        position: 3D position of detector rotation center in lab frame (mm)
        orientation: 3x3 rotation matrix (lab to detector frame)
    """
    num_rows: int
    num_cols: int
    pixel_height: float
    pixel_width: float
    beam_center_j: float  # in pixels
    beam_center_k: float  # in pixels
    position: torch.Tensor  # shape (3,)
    orientation: torch.Tensor  # shape (3, 3)


class Detector:
    """
    X-ray area detector with geometric transformations.

    The Detector class represents a 2D area detector (e.g., CCD camera) used in
    synchrotron X-ray diffraction experiments. It provides:

    - Coordinate transformations between different frames of reference
    - Detector positioning and orientation in 3D lab frame
    - Ray-plane intersection for diffraction calculations
    - Differentiable operations for gradient-based optimization

    The detector can be positioned and oriented arbitrarily in the lab frame,
    allowing for tilted and translated detector geometries.

    Example:
        >>> # Create a detector at origin, perpendicular to beam
        >>> detector = Detector(
        ...     num_rows=2048, num_cols=2048,
        ...     pixel_height=0.2, pixel_width=0.2,
        ...     beam_center_j=204.8, beam_center_k=204.8,
        ...     position=torch.tensor([100.0, 0.0, 0.0]),
        ...     dtype=torch.float32, device='cpu'
        ... )
        >>>
        >>> # Transform a point from lab frame to pixel coordinates
        >>> lab_point = torch.tensor([100.0, 10.0, 5.0])
        >>> row, col = detector.lab_to_pixel(lab_point)

    C++ Reference:
        Src/Detector.h (class CDetector, lines 103-332)
        Src/Detector.cpp (implementation)
    """

    def __init__(
        self,
        num_rows: int,
        num_cols: int,
        pixel_height: float,
        pixel_width: float,
        beam_center_j: float,
        beam_center_k: float,
        position: Optional[torch.Tensor] = None,
        orientation: Optional[torch.Tensor] = None,
        j_unit_vector: Optional[torch.Tensor] = None,
        k_unit_vector: Optional[torch.Tensor] = None,
        dtype: torch.dtype = torch.float32,
        device: Union[str, torch.device] = 'cpu'
    ):
        """
        Initialize a Detector.

        Args:
            num_rows: Number of pixel rows (K direction)
            num_cols: Number of pixel columns (J direction)
            pixel_height: Physical height of each pixel in mm (K direction)
            pixel_width: Physical width of each pixel in mm (J direction)
            beam_center_j: Beam center in PIXELS along J-axis (e.g., 1024.0 for center of 2048px detector)
            beam_center_k: Beam center in PIXELS along K-axis (e.g., 1024.0 for center of 2048px detector)
            position: 3D position of detector center in lab frame (mm), shape (3,).
                     Default: origin (0, 0, 0)
            orientation: 3x3 rotation matrix (lab to detector), shape (3, 3).
                        Default: identity (detector perpendicular to beam)
            j_unit_vector: J-axis direction in detector frame, shape (3,).
                          Default: (1, 0, 0) matching C++ convention from detector file.
            k_unit_vector: K-axis direction in detector frame, shape (3,).
                          Default: (0, -1, 0) matching C++ convention from detector file.
            dtype: PyTorch data type for tensors
            device: PyTorch device ('cpu' or 'cuda')

        Note:
            beam_center_j and beam_center_k are in PIXELS, not mm. This matches the C++ API.
            They are multiplied by pixel size internally to get the physical offset in mm.

        C++ Reference:
            Detector.cpp CDetector::CDetector (lines 66-79)
            Detector.cpp CDetector::SetImageParameter (lines 84-132)
            Detector.cpp line 102: "fBeamCenterX and Y are measured in pixels"
        """
        self.dtype = dtype
        self.device = torch.device(device)

        # Pixel dimensions
        self.num_rows = num_rows
        self.num_cols = num_cols
        self.pixel_height = pixel_height
        self.pixel_width = pixel_width
        self.pixel_half_height = pixel_height / 2.0
        self.pixel_half_width = pixel_width / 2.0

        # Beam center in detector coordinates (mm)
        self.beam_center_j = beam_center_j
        self.beam_center_k = beam_center_k

        # Position: detector rotation center in lab frame
        if position is None:
            self._position = torch.zeros(3, dtype=dtype, device=self.device)
        else:
            self._position = position.to(dtype=dtype, device=self.device)

        # Orientation: rotation matrix from lab frame to detector frame
        if orientation is None:
            self._orientation = torch.eye(3, dtype=dtype, device=self.device)
        else:
            self._orientation = orientation.to(dtype=dtype, device=self.device)

        # Detector coordinate system basis vectors (in detector frame)
        # These come from the detector file (JUnitVector, KUnitVector)
        # C++ Reference: Detector.cpp lines 105-106, 109-110
        if j_unit_vector is not None:
            self._det_frame_basis_j = j_unit_vector.to(dtype=dtype, device=self.device)
        else:
            self._det_frame_basis_j = torch.tensor([1.0, 0.0, 0.0], dtype=dtype, device=self.device)
        if k_unit_vector is not None:
            self._det_frame_basis_k = k_unit_vector.to(dtype=dtype, device=self.device)
        else:
            self._det_frame_basis_k = torch.tensor([0.0, -1.0, 0.0], dtype=dtype, device=self.device)

        # Calculate coordinate origin in detector frame
        # This is the offset from beam center to (0,0) pixel
        # C++ Reference: Detector.cpp lines 603-607 (GetCoordOrigin)
        self._det_frame_coord_origin = (
            -beam_center_j * pixel_width * self._det_frame_basis_j
            - beam_center_k * pixel_height * self._det_frame_basis_k
        )

        # Lab frame versions (will be updated by _calculate_image_plane)
        self._lab_frame_basis_j = self._det_frame_basis_j.clone()
        self._lab_frame_basis_k = self._det_frame_basis_k.clone()
        self._lab_frame_coord_origin = self._det_frame_coord_origin.clone()

        # Detector plane (in lab frame)
        # C++ Reference: Detector.h lines 131, 172-182
        self._detector_plane: Optional[Plane] = None

        # Calculate initial detector plane in lab frame
        self._calculate_image_plane()

    def _calculate_image_plane(self) -> None:
        """
        Calculate detector plane equation in lab frame.

        This method updates the detector plane based on current position and
        orientation. It's called internally whenever position or orientation changes.

        The detector plane is defined by three points on the detector face
        (in the x-y plane of detector frame), rotated by the orientation matrix
        and translated to the detector position.

        C++ Reference:
            Detector.cpp CDetector::CalculateImagePlane (lines 211-244)
            Detector.cpp lines 122-131: detector plane is the y-z plane (by convention),
            defined by points (1,0,0), (0,1,0), (0,0,0)
        """
        # Rotate the basis vectors to lab frame
        # C++ Reference: Detector.cpp lines 155-158
        self._lab_frame_basis_j = self._orientation @ self._det_frame_basis_j
        self._lab_frame_basis_k = self._orientation @ self._det_frame_basis_k
        self._lab_frame_coord_origin = self._orientation @ self._det_frame_coord_origin

        # Define detector plane using 3 corner points, matching C++ exactly
        # C++ Detector.cpp lines 125-129:
        #   v1 = (1, 0, 0), v2 = (0, 1, 0), v3 = (0, 0, 0)
        det_pt1 = torch.tensor([1.0, 0.0, 0.0], dtype=self.dtype, device=self.device)
        det_pt2 = torch.tensor([0.0, 1.0, 0.0], dtype=self.dtype, device=self.device)
        det_pt3 = torch.tensor([0.0, 0.0, 0.0], dtype=self.dtype, device=self.device)

        # Rotate then translate
        # C++ lines 215-223
        p1 = self._orientation @ det_pt1 + self._position
        p2 = self._orientation @ det_pt2 + self._position
        p3 = self._orientation @ det_pt3 + self._position

        # Compute plane normal from edges
        # C++ lines 225-229: edge1 = p3 - p1, edge2 = p2 - p1, normal = cross(edge2, edge1)
        edge1 = p3 - p1
        edge2 = p2 - p1
        normal = torch.linalg.cross(edge2, edge1)
        normal = normal / torch.norm(normal)

        # Plane equation: A*x + B*y + C*z + D = 0
        # C++ line 237: D = -dot(normal, p1)
        d = -torch.dot(normal, p1)

        # Create plane (coefficients: A, B, C, D)
        plane_coeffs = torch.cat([normal, d.unsqueeze(0)])
        self._detector_plane = Plane(coeffs=plane_coeffs)

    # =============================================================================
    # Detector Transformations
    # =============================================================================

    def set_position(self, position: torch.Tensor) -> None:
        """
        Set detector position in lab frame.

        Args:
            position: 3D position vector in lab frame (mm), shape (3,)

        C++ Reference:
            Detector.cpp CDetector::SetLocation (lines 139-143)
        """
        self._position = position.to(dtype=self.dtype, device=self.device)
        self._calculate_image_plane()

    def set_orientation(self, orientation: torch.Tensor) -> None:
        """
        Set detector orientation using a rotation matrix.

        Args:
            orientation: 3x3 rotation matrix (lab to detector frame), shape (3, 3)

        C++ Reference:
            Detector.cpp CDetector::SetOrientation (lines 150-162)
        """
        self._orientation = orientation.to(dtype=self.dtype, device=self.device)
        self._calculate_image_plane()

    def set_orientation_euler(self, phi: float, theta: float, psi: float) -> None:
        """
        Set detector orientation using Euler angles (Bunge convention).

        The Euler angles define a rotation using the ZXZ intrinsic convention:
        R = Rz(phi) @ Rx(theta) @ Rz(psi)

        Args:
            phi: First Euler angle in degrees (0-360)
            theta: Second Euler angle in degrees (0-180)
            psi: Third Euler angle in degrees (0-360)

        C++ Reference:
            Detector.cpp CDetector::SetOrientation (lines 169-173)
            XDM++/libXDM/3dMath.cpp BuildActiveEulerMatrix (lines 152-174)
        """
        # Convert to tensors
        phi_t = torch.tensor(phi, dtype=self.dtype, device=self.device)
        theta_t = torch.tensor(theta, dtype=self.dtype, device=self.device)
        psi_t = torch.tensor(psi, dtype=self.dtype, device=self.device)

        # Build rotation matrix using Euler angles
        orientation = euler_to_matrix_torch(phi_t, theta_t, psi_t)
        self.set_orientation(orientation)

    def translate(self, translation: torch.Tensor) -> None:
        """
        Translate detector by a displacement vector.

        Args:
            translation: 3D displacement vector in lab frame (mm), shape (3,)

        C++ Reference:
            Detector.cpp CDetector::Translate (lines 180-184)
        """
        self._position = self._position + translation.to(dtype=self.dtype, device=self.device)
        self._calculate_image_plane()

    def rotate(self, phi: float, theta: float, psi: float) -> None:
        """
        Rotate detector by additional Euler angles.

        This applies an additional rotation on top of the current orientation.
        The rotation is: new_orientation = R(phi, theta, psi) @ old_orientation

        Args:
            phi: First Euler angle in degrees
            theta: Second Euler angle in degrees
            psi: Third Euler angle in degrees

        C++ Reference:
            Detector.cpp CDetector::Rotate (lines 191-204)
        """
        # Convert to tensors
        phi_t = torch.tensor(phi, dtype=self.dtype, device=self.device)
        theta_t = torch.tensor(theta, dtype=self.dtype, device=self.device)
        psi_t = torch.tensor(psi, dtype=self.dtype, device=self.device)

        # Build rotation matrix
        rotation = euler_to_matrix_torch(phi_t, theta_t, psi_t)

        # Apply rotation: new = rotation @ old
        self._orientation = rotation @ self._orientation
        self._calculate_image_plane()

    # =============================================================================
    # Coordinate Transformations
    # =============================================================================

    def lab_to_detector_coordinate(self, lab_pos: torch.Tensor) -> Tuple[torch.Tensor, torch.Tensor]:
        """
        Convert lab frame position to detector frame coordinates (J, K) in mm.

        Args:
            lab_pos: 3D position in lab frame (mm), shape (3,) or (N, 3)

        Returns:
            Tuple of (j, k) detector coordinates in mm, each shape () or (N,)

        C++ Reference:
            Detector.cpp CDetector::LabToDetectorCoordinate (lines 266-274)
        """
        # Convert to detector frame relative to rotation center
        relative_pos = lab_pos - self._position

        # Project onto image coordinate system
        pixel_loc = relative_pos - self._lab_frame_coord_origin

        # Get J and K coordinates by dot product with basis vectors
        j = torch.sum(pixel_loc * self._lab_frame_basis_j, dim=-1) if pixel_loc.dim() > 1 else torch.dot(pixel_loc, self._lab_frame_basis_j)
        k = torch.sum(pixel_loc * self._lab_frame_basis_k, dim=-1) if pixel_loc.dim() > 1 else torch.dot(pixel_loc, self._lab_frame_basis_k)

        return j, k

    def detector_to_lab_coordinate(self, j: torch.Tensor, k: torch.Tensor) -> torch.Tensor:
        """
        Convert detector frame coordinates (J, K) in mm to lab frame position.

        Args:
            j: J coordinate in mm, shape () or (N,)
            k: K coordinate in mm, shape () or (N,)

        Returns:
            3D position in lab frame (mm), shape (3,) or (N, 3)

        C++ Reference:
            Detector.cpp CDetector::DetectorToLabCoordinate (lines 281-289)
        """
        # Handle scalar vs. batched inputs
        is_batched = j.dim() > 0 and j.shape[0] > 1

        if is_batched:
            # Batched operation
            j_component = j.unsqueeze(-1) * self._lab_frame_basis_j.unsqueeze(0)
            k_component = k.unsqueeze(-1) * self._lab_frame_basis_k.unsqueeze(0)
            result = j_component + k_component + self._position + self._lab_frame_coord_origin
        else:
            # Scalar operation
            result = j * self._lab_frame_basis_j + k * self._lab_frame_basis_k
            result = result + self._position + self._lab_frame_coord_origin

        return result

    def lab_to_pixel(self, lab_pos: torch.Tensor) -> Tuple[torch.Tensor, torch.Tensor]:
        """
        Convert lab frame position to pixel coordinates (row, col).

        Args:
            lab_pos: 3D position in lab frame (mm), shape (3,) or (N, 3)

        Returns:
            Tuple of (row, col) pixel coordinates (continuous), each shape () or (N,)
            Note: Pixel coordinates are continuous (float) - use floor/round for discrete indices

        C++ Reference:
            Detector.cpp CDetector::LabToPixel (lines 296-311)
        """
        # First convert to detector coordinates
        j, k = self.lab_to_detector_coordinate(lab_pos)

        # Convert to pixel coordinates
        # C++ Reference: Detector.cpp ToRowPixel, ToColPixel (lines 528-548)
        row = self._to_row_pixel(k)
        col = self._to_col_pixel(j)

        return row, col

    def pixel_to_lab_coordinate(self, col: torch.Tensor, row: torch.Tensor) -> torch.Tensor:
        """
        Convert pixel coordinates to lab frame position.

        Args:
            col: Column pixel coordinate (X), shape () or (N,)
            row: Row pixel coordinate (Y), shape () or (N,)

        Returns:
            3D position in lab frame (mm), shape (3,) or (N, 3)

        C++ Reference:
            Detector.cpp CDetector::PixelToLabCoordinate (lines 317-326)
        """
        # Convert pixel to detector coordinates (mm)
        j = self._col_pixel_to_image_j(col)
        k = self._row_pixel_to_image_k(row)

        # Convert detector coordinates to lab frame
        return self.detector_to_lab_coordinate(j, k)

    def _col_pixel_to_image_j(self, col: torch.Tensor) -> torch.Tensor:
        """Convert column pixel to J coordinate (mm)."""
        # C++ Reference: Detector.cpp ColPixelToImageJ (lines 515-518)
        return col * self.pixel_width

    def _row_pixel_to_image_k(self, row: torch.Tensor) -> torch.Tensor:
        """Convert row pixel to K coordinate (mm)."""
        # C++ Reference: Detector.cpp RowPixelToImageK (lines 507-510)
        return row * self.pixel_height

    def _to_row_pixel(self, k: torch.Tensor) -> torch.Tensor:
        """
        Convert K coordinate (mm) to row pixel coordinate.

        C++ Reference: Detector.cpp ToRowPixel (lines 528-535)
        """
        # Note: C++ version truncates to integer and handles negative values
        # Python version returns continuous coordinates for differentiability
        return (k + self.pixel_half_height) / self.pixel_height

    def _to_col_pixel(self, j: torch.Tensor) -> torch.Tensor:
        """
        Convert J coordinate (mm) to column pixel coordinate.

        C++ Reference: Detector.cpp ToColPixel (lines 542-548)
        """
        # Note: C++ version truncates to integer and handles negative values
        # Python version returns continuous coordinates for differentiability
        return (j + self.pixel_half_width) / self.pixel_width

    # =============================================================================
    # Ray Intersection
    # =============================================================================

    def intersect_ray(self, ray: Ray) -> Tuple[torch.Tensor, torch.Tensor]:
        """
        Compute ray-detector plane intersection.

        Args:
            ray: Ray to intersect with detector plane

        Returns:
            Tuple of (intersects, t):
                intersects: Boolean tensor indicating if ray intersects (shape ())
                t: Parameter value where ray intersects plane (shape ())
                   Ray point at intersection: ray.origin + t * ray.direction

        C++ Reference:
            Detector.cpp CDetector::Intersects (lines 254-259)
        """
        if self._detector_plane is None:
            raise RuntimeError("Detector plane not initialized")

        return ray.intersect_plane(self._detector_plane)

    # =============================================================================
    # Properties and Accessors
    # =============================================================================

    @property
    def position(self) -> torch.Tensor:
        """
        Detector position (rotation center) in lab frame (mm).

        Returns:
            3D position vector, shape (3,)

        C++ Reference:
            Detector.cpp CDetector::GetLocation (lines 373-377)
            Detector.cpp CDetector::GetRotationCenter (lines 390-394)
        """
        return self._position.clone()

    @property
    def orientation(self) -> torch.Tensor:
        """
        Detector orientation matrix (lab to detector frame).

        Returns:
            3x3 rotation matrix, shape (3, 3)

        C++ Reference:
            Detector.cpp CDetector::GetOrientationMatrix (lines 382-386)
        """
        return self._orientation.clone()

    @property
    def detector_width(self) -> float:
        """
        Total detector width in mm.

        C++ Reference:
            Detector.cpp CDetector::GetDetectorWidth (lines 467-470)
        """
        return float(self.num_cols * self.pixel_width)

    @property
    def detector_height(self) -> float:
        """
        Total detector height in mm.

        C++ Reference:
            Detector.cpp CDetector::GetDetectorHeight (lines 475-478)
        """
        return float(self.num_rows * self.pixel_height)

    @property
    def coordinate_origin(self) -> torch.Tensor:
        """
        Detector coordinate origin in lab frame.

        This is the position of the (0,0) pixel corner in the lab frame.

        Returns:
            3D position vector, shape (3,)

        C++ Reference:
            Detector.cpp CDetector::GetDetectorCoordinateOrigin (lines 432-436)
        """
        return self._lab_frame_coord_origin + self._position

    @property
    def basis_vectors(self) -> Tuple[torch.Tensor, torch.Tensor]:
        """
        Detector coordinate basis vectors in lab frame.

        Returns:
            Tuple of (j_basis, k_basis), each shape (3,)

        C++ Reference:
            Detector.cpp CDetector::GetCoordinateBasis (lines 458-462)
        """
        return self._lab_frame_basis_j.clone(), self._lab_frame_basis_k.clone()

    @property
    def detector_plane(self) -> Plane:
        """
        Detector plane in lab frame.

        Returns:
            Plane object representing detector surface

        C++ Reference:
            Detector.cpp CDetector::GetLabDetectorPlane (lines 499-502)
        """
        if self._detector_plane is None:
            raise RuntimeError("Detector plane not initialized")
        return self._detector_plane

    def in_range(self, col: torch.Tensor, row: torch.Tensor) -> torch.Tensor:
        """
        Check if pixel coordinates are within detector bounds.

        Args:
            col: Column pixel coordinates, any shape
            row: Row pixel coordinates, any shape

        Returns:
            Boolean tensor indicating if coordinates are in range, same shape as input

        C++ Reference:
            Detector.h CDetector::InRange (lines 308-311)
        """
        in_bounds = (
            (row >= 0) & (row < self.num_rows) &
            (col >= 0) & (col < self.num_cols)
        )
        return in_bounds

    def get_parameters(self) -> DetectorParameters:
        """
        Get detector parameters as a dataclass.

        Returns:
            DetectorParameters containing all geometric parameters

        C++ Reference:
            Detector.cpp CXDMDetectorFactory::GetImageParameters (lines 672-687)
        """
        return DetectorParameters(
            num_rows=self.num_rows,
            num_cols=self.num_cols,
            pixel_height=self.pixel_height,
            pixel_width=self.pixel_width,
            beam_center_j=self.beam_center_j,
            beam_center_k=self.beam_center_k,
            position=self.position,
            orientation=self.orientation
        )

    # =============================================================================
    # Factory Methods
    # =============================================================================

    @classmethod
    def from_beam_center(
        cls,
        num_rows: int,
        num_cols: int,
        pixel_height: float,
        pixel_width: float,
        beam_center_j: float,
        beam_center_k: float,
        position: torch.Tensor,
        orientation: Optional[torch.Tensor] = None,
        dtype: torch.dtype = torch.float32,
        device: Union[str, torch.device] = 'cpu'
    ) -> 'Detector':
        """
        Create a Detector from beam center position.

        This is a convenience factory method that matches the C++ factory pattern.
        It's equivalent to using the regular constructor but provides a more
        explicit interface for HEDM experiments.

        Args:
            num_rows: Number of pixel rows
            num_cols: Number of pixel columns
            pixel_height: Pixel height in mm
            pixel_width: Pixel width in mm
            beam_center_j: Beam center in PIXELS along J-axis
            beam_center_k: Beam center in PIXELS along K-axis
            position: Detector position in lab frame (mm)
            orientation: Detector orientation matrix (optional)
            dtype: PyTorch data type
            device: PyTorch device

        Returns:
            Detector instance

        C++ Reference:
            Detector.cpp CXDMDetectorFactory::MakeDetector (lines 613-638)
        """
        return cls(
            num_rows=num_rows,
            num_cols=num_cols,
            pixel_height=pixel_height,
            pixel_width=pixel_width,
            beam_center_j=beam_center_j,
            beam_center_k=beam_center_k,
            position=position,
            orientation=orientation,
            dtype=dtype,
            device=device
        )

    # =============================================================================
    # Image Operations
    # =============================================================================

    def add_direct_beam(
        self,
        image: 'ImageData',
        beam_height: float,
        beam_width: float,
        intensity: float = 1.0,
        mode: str = 'soft'
    ) -> None:
        """
        Add beam aperture polygon to detector image.

        Rasterizes a rectangular beam aperture region onto the detector image.
        The beam is centered at the detector's beam center position.

        This is used in forward simulation to mark the direct beam footprint
        on the detector.

        Args:
            image: ImageData instance to rasterize onto
                  Must match detector dimensions (num_rows, num_cols)
            beam_height: Beam aperture height in mm
            beam_width: Beam aperture width in mm
            intensity: Intensity value for beam region (default 1.0)
            mode: Rasterization mode ('soft' or 'hard'), passed to add_polygon

        Raises:
            ValueError: If image dimensions don't match detector

        Example:
            >>> from icenine.detector import Detector
            >>> from icenine.image_data import ImageData
            >>>
            >>> detector = Detector(
            ...     num_rows=2048, num_cols=2048,
            ...     pixel_height=0.2, pixel_width=0.2,
            ...     beam_center_j=1024.0, beam_center_k=1024.0,
            ...     position=torch.tensor([100.0, 0.0, 0.0])
            ... )
            >>>
            >>> image = ImageData(2048, 2048, mode='dense')
            >>> detector.add_direct_beam(image, beam_height=2.0, beam_width=2.0)

        C++ Reference:
            Detector.cpp CDetector::AddDirectBeam (lines 563-596)
        """
        # Import here to avoid circular dependency
        from .image_data import ImageData

        # Validate image dimensions match detector
        if image.num_rows != self.num_rows or image.num_cols != self.num_cols:
            raise ValueError(
                f"Image dimensions ({image.num_rows}, {image.num_cols}) "
                f"do not match detector ({self.num_rows}, {self.num_cols})"
            )

        # Calculate beam center in detector coordinates (mm)
        # beam_center_j and beam_center_k are in PIXELS
        beam_center_j_mm = self.beam_center_j * self.pixel_width
        beam_center_k_mm = self.beam_center_k * self.pixel_height

        # Calculate polygon corners in detector frame (mm)
        # Beam aperture is centered at beam center
        half_width = beam_width / 2.0
        half_height = beam_height / 2.0

        corners_mm = torch.tensor([
            [beam_center_j_mm - half_width, beam_center_k_mm - half_height],
            [beam_center_j_mm + half_width, beam_center_k_mm - half_height],
            [beam_center_j_mm + half_width, beam_center_k_mm + half_height],
            [beam_center_j_mm - half_width, beam_center_k_mm + half_height],
        ], dtype=self.dtype, device=self.device)

        # Convert mm to pixel coordinates
        # j = column, k = row
        corners_pixels = corners_mm.clone()
        corners_pixels[:, 0] /= self.pixel_width   # j (cols) - mm to pixels
        corners_pixels[:, 1] /= self.pixel_height  # k (rows) - mm to pixels

        # Rasterize polygon onto image
        image.add_polygon(corners_pixels, intensity=intensity, mode=mode)

    def __repr__(self) -> str:
        """String representation of Detector."""
        return (
            f"Detector(\n"
            f"  num_rows={self.num_rows}, num_cols={self.num_cols},\n"
            f"  pixel_size=({self.pixel_width:.4f}, {self.pixel_height:.4f}) mm,\n"
            f"  beam_center=({self.beam_center_j:.2f}, {self.beam_center_k:.2f}) mm,\n"
            f"  position={self.position.tolist()},\n"
            f"  device={self.device}, dtype={self.dtype}\n"
            f")"
        )
