"""
Experimental data loader for reconstruction.

Loads detector images from disk (ASCII .d files) or from forward simulation
output into a 2D array indexed by [omega_interval][detector_index].

C++ Reference:
    Src/SimulationData.h  — CSimulationData / ImageMapT
    Src/ReconstructionSetup.cpp — ReadExperimentalData loop
"""

from pathlib import Path
from typing import List, Optional, TYPE_CHECKING, Union

import numpy as np

from icenine.config_file import ConfigFile
from icenine.experiment_setup import XDMExperimentSetup
from icenine.image_data import ImageData

if TYPE_CHECKING:
    from icenine.differentiable_cost import ExperimentalImageStack, SparseImageStack


class ExperimentalData:
    """
    Container for experimental detector images used during reconstruction.

    Wraps a 2D array of ImageData objects:
        images[omega_interval_index][detector_index]

    C++ Reference:
        SimulationData.h — boost::multi_array<CSearchableImageData, 2> mImageMap
    """

    def __init__(
        self,
        images: List[List[ImageData]],
        n_omega_intervals: int,
        n_detectors: int,
    ):
        """
        Initialize with pre-loaded images.

        Args:
            images: 2D list [omega_interval][detector] of ImageData
            n_omega_intervals: Number of omega intervals
            n_detectors: Number of detectors
        """
        self.images = images
        self.n_omega_intervals = n_omega_intervals
        self.n_detectors = n_detectors

    def get_image(self, omega_index: int, detector_index: int) -> ImageData:
        """
        Get the experimental image for a given omega interval and detector.

        Args:
            omega_index: Omega interval index
            detector_index: Detector index

        Returns:
            ImageData at the specified position
        """
        return self.images[omega_index][detector_index]

    @classmethod
    def from_forward_simulation(cls, forward_sim) -> "ExperimentalData":
        """
        Create ExperimentalData from ForwardSimulation output.

        This is the primary method for synthetic test data: run a forward
        simulation, then use its output directly as "experimental" data
        for reconstruction validation.

        Args:
            forward_sim: ForwardSimulation instance (must have .images populated)

        Returns:
            ExperimentalData wrapping the simulation images
        """
        images = forward_sim.images
        if not images:
            raise ValueError("ForwardSimulation has no images. Run simulate_detector_images() first.")

        n_omega = len(images)
        n_det = len(images[0])
        result = cls(images=images, n_omega_intervals=n_omega, n_detectors=n_det)
        result.prepare_for_reconstruction()
        return result

    @classmethod
    def from_ascii_files(
        cls,
        config: ConfigFile,
        exp_setup: XDMExperimentSetup,
    ) -> "ExperimentalData":
        """
        Load experimental data from ASCII detector image files.

        Reads files matching the naming convention:
            {basename}{file_number:0{serial_length}d}.{ext}{detector_offset + det_idx}

        Args:
            config: ConfigFile with InfileBasename, InfileExtension, InfileSerialLength
            exp_setup: Initialized XDMExperimentSetup (provides detectors, omega ranges)

        Returns:
            ExperimentalData loaded from disk

        C++ Reference:
            ReconstructionSetup.cpp:55-108 (loading loop)
        """
        detector_list = exp_setup.get_detector_list()
        omega_ranges = exp_setup.get_omega_range_list()
        file_ranges = exp_setup.get_file_range_list()

        n_omega = len(omega_ranges)
        n_det = len(detector_list)

        basename = config.in_file_basename
        ext = config.in_file_ext
        serial_length = config.in_file_serial_length
        det_offset = config.bc_peak_detector_offset

        import time as _time

        # Initialize 2D image array
        images: List[List[Optional[ImageData]]] = [
            [None for _ in range(n_det)] for _ in range(n_omega)
        ]

        total_files = n_omega * n_det
        loaded = 0
        t_start = _time.time()

        for det_idx, detector in enumerate(detector_list):
            file_range = file_ranges[det_idx]
            for omega_idx in range(n_omega):
                file_num = file_range.low + omega_idx
                file_num_str = str(file_num).zfill(serial_length)
                filename = f"{basename}{file_num_str}.{ext}{det_offset + det_idx}"

                filepath = Path(filename)
                if not filepath.exists():
                    raise FileNotFoundError(
                        f"Experimental data file not found: {filepath}"
                    )

                image = ImageData(detector.num_rows, detector.num_cols)
                image.load_ascii(str(filepath))
                images[omega_idx][det_idx] = image
                loaded += 1

                if loaded % 60 == 0 or loaded == total_files:
                    elapsed = _time.time() - t_start
                    print(f"  Loading images: {loaded}/{total_files} "
                          f"({elapsed:.1f}s)", flush=True)

        result = cls(
            images=images,  # type: ignore[arg-type]
            n_omega_intervals=n_omega,
            n_detectors=n_det,
        )
        print(f"  Preparing binary caches...", flush=True)
        result.prepare_for_reconstruction()
        return result

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
        file_start: int = 0,
        det_offset: int = 0,
        mode: str = "dense",
    ) -> "ExperimentalData":
        """
        Load experimental data from a directory of ASCII image files.

        Simpler interface than from_ascii_files when you don't have a full
        config/experiment setup (e.g., loading forward sim output for testing).

        Args:
            directory: Directory containing the image files
            basename: File basename (e.g., "3Grains.sim")
            ext: File extension (e.g., "d")
            serial_length: Digits in serial number (e.g., 5 for "00000")
            n_omega: Number of omega intervals
            n_detectors: Number of detectors
            num_rows: Detector image height in pixels
            num_cols: Detector image width in pixels
            file_start: Starting file number (default 0)
            det_offset: Detector numbering offset (default 0)
            mode: ImageData storage mode - 'dense' (default) or 'sparse'.
                  Use 'sparse' for large images to reduce memory from
                  O(n_images × H × W) to O(total_nonzero_pixels).
                  Note: prepare_for_reconstruction() (which eagerly builds the
                  uint8 binary caches VoxelCostFunction's hard cost path reads)
                  only runs for mode='dense', below. Combining mode='sparse'
                  with the hard VoxelCostFunction is not supported.

        Returns:
            ExperimentalData loaded from directory
        """
        import time as _time

        directory = Path(directory)
        total_files = n_omega * n_detectors
        loaded = 0

        images: List[List[Optional[ImageData]]] = [
            [None for _ in range(n_detectors)] for _ in range(n_omega)
        ]

        t_start = _time.time()
        for det_idx in range(n_detectors):
            for omega_idx in range(n_omega):
                file_num = file_start + omega_idx
                file_num_str = str(file_num).zfill(serial_length)
                filename = f"{basename}{file_num_str}.{ext}{det_offset + det_idx}"
                filepath = directory / filename

                if not filepath.exists():
                    raise FileNotFoundError(
                        f"Image file not found: {filepath}"
                    )

                image = ImageData(num_rows, num_cols, mode=mode)
                image.load_ascii(str(filepath))
                images[omega_idx][det_idx] = image
                loaded += 1

                if loaded % 60 == 0 or loaded == total_files:
                    elapsed = _time.time() - t_start
                    print(f"  Loading images: {loaded}/{total_files} "
                          f"({elapsed:.1f}s)", flush=True)

        result = cls(
            images=images,  # type: ignore[arg-type]
            n_omega_intervals=n_omega,
            n_detectors=n_detectors,
        )
        if mode == "dense":
            print(f"  Preparing binary caches...", flush=True)
            result.prepare_for_reconstruction()
        return result

    @classmethod
    def from_binary_memmap(
        cls,
        path: Union[str, Path],
        n_omega: int,
        n_detectors: int,
    ) -> "ExperimentalData":
        """Opt-in low-memory loader: binary (pixel > 0) frames from a uint8 ``.npy`` stack.

        The stack (see ``write_binary_stack``) has shape (n_detectors * n_omega, H, W), frame
        ``det * n_omega + omega``. It is opened with ``np.load(mmap_mode="r")`` and every
        ``ImageData._binary_cache`` is set to its read-only C-contiguous slice, so the float32
        stack is never built (per process: only the pages the cost function touches, shared
        through the OS page cache between concurrent runs). The images are *sparse-mode, empty*
        ImageData: only the hard-cost path (``get_binary_numpy``) is supported. Intensities,
        ``count_bright_pixels`` and the soft cost see no pixels. Default loaders are unchanged.
        """
        stack = np.load(str(path), mmap_mode="r")
        if stack.dtype != np.uint8 or stack.ndim != 3:
            raise ValueError(
                f"expected a uint8 (frames, H, W) stack, got {stack.dtype} {stack.shape}"
            )
        if stack.shape[0] != n_omega * n_detectors:
            raise ValueError(f"stack has {stack.shape[0]} frames, expected {n_omega * n_detectors}")
        rows, cols = stack.shape[1], stack.shape[2]
        images: List[List[ImageData]] = [
            [None for _ in range(n_detectors)] for _ in range(n_omega)  # type: ignore[misc]
        ]
        for d in range(n_detectors):
            for o in range(n_omega):
                im = ImageData(rows, cols, mode="sparse")
                frame = stack[d * n_omega + o]
                assert frame.flags["C_CONTIGUOUS"]
                im._binary_cache = frame
                images[o][d] = im
        return cls(images=images, n_omega_intervals=n_omega, n_detectors=n_detectors)

    def prepare_for_reconstruction(self) -> None:
        """
        Pre-compute binary caches for all images.

        Call this before reconstruction to eagerly populate the uint8 binary
        arrays, avoiding lazy computation during the hot cost function loop.
        """
        for omega_idx in range(self.n_omega_intervals):
            for det_idx in range(self.n_detectors):
                self.images[omega_idx][det_idx].ensure_binary_cache()

    def to_image_stack(
        self, binary: bool = True
    ) -> "ExperimentalImageStack":
        """
        Convert to a pre-stacked contiguous tensor for batch access.

        This creates an ExperimentalImageStack with all images in a single
        (n_omega * n_det, 1, H, W) tensor, enabling batch grid_sample
        and GPU acceleration.

        Args:
            binary: If True (default), binarize to 0.0/1.0.
                    If False, preserve original intensities.

        Returns:
            ExperimentalImageStack ready for differentiable cost function.
        """
        from .differentiable_cost import ExperimentalImageStack
        return ExperimentalImageStack(self, binary=binary)

    def to_sparse_image_stack(
        self, binary: bool = True
    ) -> "SparseImageStack":
        """
        Convert to a memory-efficient sparse image stack.

        Stores only the coordinates of bright pixels (~26 KB for typical
        diffraction data vs ~5.6 GB for dense float32). Dense images are
        materialized on-demand for grid_sample.

        Args:
            binary: If True (default), store only pixel coordinates (values=1.0).
                    If False, store coordinates and intensity values.

        Returns:
            SparseImageStack ready for differentiable cost function.
        """
        from .differentiable_cost import SparseImageStack
        return SparseImageStack(self, binary=binary)

    def count_bright_pixels(self) -> int:
        """Count total bright pixels across all images."""
        total = 0
        for omega_idx in range(self.n_omega_intervals):
            for det_idx in range(self.n_detectors):
                image = self.images[omega_idx][det_idx]
                if image._mode == 'dense':
                    total += (image._pixels_dense > 0).sum().item()
                else:
                    total += image._pixels_sparse._values().numel()
        return int(total)

    def __repr__(self) -> str:
        return (
            f"ExperimentalData("
            f"omega_intervals={self.n_omega_intervals}, "
            f"detectors={self.n_detectors})"
        )


def write_binary_stack(
    directory: Union[str, Path],
    basename: str,
    ext: str,
    serial_length: int,
    n_omega: int,
    n_detectors: int,
    num_rows: int,
    num_cols: int,
    out_path: Union[str, Path],
    file_start: int = 0,
    det_offset: int = 0,
) -> Path:
    """Convert a directory of ASCII frames (the ``from_image_directory`` layout) into one
    uint8 ``.npy`` stack (frame ``det * n_omega + omega``) for ``from_binary_memmap``.

    Written through ``open_memmap`` frame by frame, so peak memory is one frame.
    """
    directory, out_path = Path(directory), Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    stack = np.lib.format.open_memmap(
        str(out_path), mode="w+", dtype=np.uint8, shape=(n_omega * n_detectors, num_rows, num_cols)
    )
    for d in range(n_detectors):
        for o in range(n_omega):
            name = f"{basename}{str(file_start + o).zfill(serial_length)}.{ext}{det_offset + d}"
            fn = directory / name
            if not fn.exists():
                raise FileNotFoundError(f"Image file not found: {fn}")
            ks, js = [], []
            with open(fn) as f:
                for line in f:
                    line = line.strip()
                    if not line or line.startswith("#") or line.startswith(","):
                        continue
                    parts = line.split(",")
                    if len(parts) >= 3 and float(parts[2]) > 0:
                        js.append(int(parts[0]))
                        ks.append(int(parts[1]))
            frame = np.zeros((num_rows, num_cols), dtype=np.uint8)
            frame[ks, js] = 1
            stack[d * n_omega + o] = frame
    stack.flush()
    del stack
    return out_path
