"""
Forward diffraction simulation engine.

Implements the main simulation loop for generating synthetic detector images
from polycrystalline samples. Uses the core Simulation class to project
diffraction peaks onto detectors across omega rotation ranges.

Python port of Src/ForwardSimulation.h/cpp

Author: S. F. Li
"""

from typing import List, Optional
from pathlib import Path
import numpy as np
import torch

from .config_file import ConfigFile
from .experiment_setup import XDMExperimentSetup
from .simulation import Simulation, PeakInfo
from .sample import Sample
from .detector import Detector
from .image_data import ImageData
from .peak_filters import XDMEtaAcceptFn
from .constants import KEV_OVER_HBAR_C_IN_ANG


class ForwardSimulation:
    """
    Forward diffraction simulation engine.

    Generates synthetic detector images from polycrystalline samples by
    simulating X-ray diffraction peaks across omega rotation ranges.

    C++ Reference:
        ForwardSimulation.h:64-123 class CXDMForwardSimulation

    Attributes:
        config_file: Configuration file with experiment parameters
        exp_setup: Initialized experiment setup
        simulator: Core simulation engine
        images: 2D list of detector images [omega_index][detector_index]

    Algorithm:
        Triple nested loop over:
        1. Voxels in sample
        2. Reflections from crystal structure
        3. Omega angles where Bragg condition is satisfied

        For each (voxel, reflection, omega):
            - Rotate sample to omega angle
            - Project voxel onto detector(s)
            - Accumulate intensity in image

    Example:
        >>> from icenine.config_file import ConfigFile
        >>> config = ConfigFile.from_file("experiment.config")
        >>> simulator = ForwardSimulation(config)
        >>> simulator.simulate_detector_images()
    """

    def __init__(self, config_file: ConfigFile):
        """
        Initialize forward simulation.

        Args:
            config_file: Configuration file with experiment parameters

        C++ Reference:
            ForwardSimulation.cpp:54-58 CXDMForwardSimulation constructor
        """
        self.config_file = config_file
        self.exp_setup = XDMExperimentSetup(config_file)
        self.simulator = Simulation()
        self.images: List[List[ImageData]] = []

    def simulate_detector_images(
        self,
        sample: Optional[Sample] = None,
        output_dir: Optional[Path] = None
    ) -> List[List[ImageData]]:
        """
        Main entry point for forward simulation.

        Generates detector images for all omega ranges and detectors,
        optionally saving to disk.

        Args:
            sample: Sample to simulate (if None, loaded from config)
            output_dir: Output directory for images (if None, uses config)

        Returns:
            2D list of detector images [omega_index][detector_index]

        C++ Reference:
            ForwardSimulation.cpp:65-128 SimulateDetectorImagesOptimized

        Algorithm:
            1. Initialize experiment setup (read files)
            2. Initialize sample (or use provided)
            3. Create image storage
            4. Run simulation loop
            5. Save images to disk

        Example:
            >>> simulator = ForwardSimulation(config)
            >>> images = simulator.simulate_detector_images()
            >>> print(f"Generated {len(images)} omega steps")
        """
        print("Initializing experiment...")

        # Initialize experiment (read detector files, omega ranges, etc.)
        # C++: oExpSetup.InitializeExperiment()
        self.exp_setup.initialize_experiment()

        # Get experimental parameters
        # C++: const vector<SRange> & vOmegaRangeList = oExpSetup.GetOmegaRangeList()
        omega_ranges = self.exp_setup.get_omega_range_list()
        file_ranges = self.exp_setup.get_file_range_list()
        detector_list = self.exp_setup.get_detector_list()
        range_map = self.exp_setup.get_range_to_index_map()

        # Initialize simulator
        # C++: oSimulator.Initialize(oExpSetup)
        self.simulator = Simulation(self.exp_setup)

        # Initialize or use provided sample
        if sample is None:
            # C++: oExpSetup.InitializeSample(oCurrentLayer, oDetectorList[0])
            sample = Sample()
            self.exp_setup.initialize_sample(sample, detector_list[0])

        # Create image storage
        # C++: ImageMap oSimData;
        # C++: oSimData.resize(boost::extents[vOmegaRangeList.size()][oDetectorList.size()])
        print(f"Creating {len(omega_ranges)} x {len(detector_list)} image array...")

        self.images = []
        for i in range(len(omega_ranges)):
            detector_images = []
            for detector in detector_list:
                # C++: oSimData[i][j].Resize(oDetectorList[j].GetNumCols(),
                #                            oDetectorList[j].GetNumRows())
                # C++: oSimData[i][j].Fill(0)
                image = ImageData(detector.num_rows, detector.num_cols)
                detector_images.append(image)
            self.images.append(detector_images)

        # Run simulation
        print("Begin Simulation")
        self._simulate_peaks(
            self.images,
            detector_list,
            sample,
            range_map
        )
        print("Finished Simulation")

        # Output images
        if output_dir is None:
            output_dir = Path(".")

        self._save_images(
            self.images,
            omega_ranges,
            file_ranges,
            detector_list,
            output_dir
        )

        return self.images

    def _simulate_peaks(
        self,
        images: List[List[ImageData]],
        detector_list: List[Detector],
        sample: Sample,
        range_map
    ):
        """
        Core simulation loop: iterate over voxels, reflections, and omegas.

        This is the computational kernel of the forward simulation.

        Args:
            images: 2D list of detector images to accumulate into
            detector_list: List of detector geometries
            sample: Sample with voxel grid
            range_map: Omega range to index mapping

        C++ Reference:
            ForwardSimulation.cpp:200-284 SimulatePeaks

        Algorithm:
            FOR each voxel v in sample:
                Get crystal structure for voxel phase
                Get reciprocal vectors G_hkl from structure

                FOR each reflection G_hkl:
                    Transform to lab frame: G' = O * G_hkl
                    Solve Bragg condition for omega angles

                    IF peak is observable:
                        FOR each omega solution (ω₁, ω₂):
                            Find omega index in range map
                            IF omega in valid range:
                                Rotate sample to omega
                                Project voxel onto detector(s)
                                Add intensity to image[omega_index][detector_index]
                                Reset sample orientation

        Performance:
            - C++ cannot parallelize due to sparse matrix memory access
            - Python version could potentially use thread-safe image accumulation
            - Progress reported every 10000 voxels
        """
        # Create peak acceptance filter
        # C++: HEDM::XDMEtaAcceptFn FAcceptFn(-oExpSetup.GetEtaLimit(), oExpSetup.GetEtaLimit())
        eta_limit = self.exp_setup.get_eta_limit()

        # Calculate wavenumber for Bragg angle calculations
        # C++: Float fWavenumber = PhysicalConstants::keV_over_hbar_c_in_ang * oExpSetup.GetBeamEnergy()
        wavenumber = KEV_OVER_HBAR_C_IN_ANG * self.exp_setup.beam_energy

        # Get crystal structures
        # C++: const vector<CUnitCell> & oCryStructList = oCurrentLayer.GetStructureList()
        structure_list = sample.get_structure_list()

        # Get voxel grid
        # C++: std::shared_ptr<CMic> pMic = std::dynamic_pointer_cast<CMic>(oCurrentLayer.GetMic())
        mic = sample.get_mic()
        voxel_list = mic.voxels

        print(f"Simulating {len(voxel_list)} voxels...")

        # Main simulation loop
        # C++: for(vector<SVoxel>::const_iterator pCurVoxel = pMic->VoxelListBegin(); ...)
        voxel_count = 0

        for voxel in voxel_list:
            voxel_count += 1

            # Progress reporting
            # C++: if(nVoxelCount % 10000 == 0) std::cout << nVoxelCount << std::endl
            if voxel_count % 10000 == 0:
                print(f"  Voxel {voxel_count}/{len(voxel_list)}")

            # Get crystal structure for this voxel's phase
            # C++: Int nCryStructIndex = pCurVoxel->nPhase
            # C++: const vector<CRecpVector> & oRecipVectors = oCryStructList[nCryStructIndex].GetReflectionVectorList()
            phase_index = voxel.phase
            if phase_index >= len(structure_list):
                continue  # Skip invalid phase

            crystal_structure = structure_list[phase_index]
            reciprocal_vectors = crystal_structure.get_reflection_vectors()

            # Get voxel orientation as PyTorch tensor
            # C++: pCurVoxel->oOrientMatrix
            voxel_orientation = torch.from_numpy(voxel.orientation).float()

            # Process each reflection
            # C++: for(Size_Type nRecipIndex = 0; nRecipIndex < oRecipVectors.size(); nRecipIndex++)
            for recp_idx, recp_vector in enumerate(reciprocal_vectors):
                # Transform scattering vector to lab frame
                # C++: SVector3 oScatteringVec = oRecipVectors[nRecipIndex].v
                # C++: oScatteringVec.Transform(pCurVoxel->oOrientMatrix)  // g_hkl' = O * g_hkl
                g_hkl = torch.from_numpy(recp_vector.q_vec).float()
                g_lab = voxel_orientation @ g_hkl
                g_magnitude = recp_vector.q_mag

                # Solve Bragg condition for omega angles
                # C++: bool bPeakObservable = oSimulator.GetScatteringOmegas(fOmegaRes[0], fOmegaRes[1], ...)
                from .diffraction_core import get_scattering_omegas_torch

                result = get_scattering_omegas_torch(
                    g_lab.unsqueeze(0),  # Add batch dimension
                    torch.tensor([g_magnitude]),
                    self.exp_setup.beam_energy,
                    self.exp_setup.get_beam_deflection_chi_laue()
                )

                if not result.observable[0]:
                    continue  # Peak not observable

                # Calculate sin(2θ) for Lorentz-polarization correction
                # C++: Float fSinTheta = oRecipVectors[nRecipIndex].fMag / (Float(2.0) * fWavenumber)
                # C++: FAcceptFn.fSin2Theta = sin(Float(2) * asin(fSinTheta))
                sin_theta = g_magnitude / (2.0 * wavenumber)
                sin_2theta = np.sin(2.0 * np.arcsin(sin_theta))

                # Normalize scattering direction
                # C++: oScatteringVec.Normalize()
                # C++: const SVector3 & oScatteringDir = oScatteringVec
                scattering_dir = g_lab / torch.norm(g_lab)

                # Process both omega solutions
                # C++: for(int i = 0; i < 2; i++)
                omega_solutions = [result.omega1[0].item(), result.omega2[0].item()]

                for omega in omega_solutions:
                    # Find omega index in range map
                    # C++: Size_Type nOmegaIndex = oRangeToIndexMap(fOmegaRes[i])
                    # C++: if(nOmegaIndex != XDMSimulation::NoMatch)
                    omega_index = range_map.angle_to_wedge_index(omega)

                    if omega_index is None:
                        continue  # Omega outside valid ranges

                    # Save current sample orientation
                    # C++: const SVector3 oCurOrientation = oCurrentLayer.GetOrientation()
                    current_orientation = sample.get_orientation()

                    # Rotate sample to omega angle
                    # C++: oCurrentLayer.RotateZ(fOmegaRes[i])
                    sample.rotate_z(omega)

                    # Create peak filter with this reflection's intensity
                    # C++: FAcceptFn.fFormIntensity = oRecipVectors[nRecipIndex].fIntensity
                    peak_filter = XDMEtaAcceptFn(
                        min_eta=-eta_limit,
                        max_eta=eta_limit,
                        form_intensity=recp_vector.intensity,
                        sin_2theta=sin_2theta
                    )

                    # Project voxel onto all detectors
                    # C++: oSimulator.ProjectVoxel(oCurImageList, vDetectorList, oCurrentLayer,
                    #                              *pCurVoxel, oScatteringDir, FAcceptFn)
                    for det_idx, detector in enumerate(detector_list):
                        image = images[omega_index][det_idx]

                        # Get voxel vertices in sample frame
                        vertices = self._get_voxel_vertices(voxel)

                        # Project voxel
                        self.simulator.project_voxel(
                            image,
                            detector,
                            sample,
                            vertices,
                            scattering_dir,
                            peak_filter
                        )

                    # Restore sample orientation
                    # C++: oCurrentLayer.SetOrientation(oCurOrientation.m_fX, ...)
                    sample.set_orientation(*current_orientation)

    def _get_voxel_vertices(self, voxel) -> torch.Tensor:
        """
        Get triangular vertices for voxel projection.

        Computes the 3 vertices of the equilateral triangle voxel,
        matching the C++ implementation in MicIO.h lines 300-315.

        Args:
            voxel: Voxel with position, side_length, and points_up fields

        Returns:
            Vertices tensor, shape (3, 3) - counter-clockwise winding

        C++ Reference:
            XDM++/libXDM/MicIO.h lines 300-315
        """
        import math

        x = float(voxel.position[0])
        y = float(voxel.position[1])
        z = float(voxel.position[2])
        s = float(voxel.side_length)

        if voxel.points_up:
            # UP triangle (direction=1) - counter-clockwise winding
            # C++ MicIO.h lines 302-305
            vertices = torch.tensor([
                [x,           y,                          z],
                [x + s,       y,                          z],
                [x + s / 2.0, y + s / 2.0 * math.sqrt(3.0), z],
            ], dtype=torch.float32)
        else:
            # DOWN triangle (direction=2) - counter-clockwise winding
            # C++ MicIO.h lines 310-313
            vertices = torch.tensor([
                [x,           y,                            z],
                [x + s / 2.0, y - s / 2.0 * math.sqrt(3.0), z],
                [x + s,       y,                            z],
            ], dtype=torch.float32)

        return vertices

    def _save_images(
        self,
        images: List[List[ImageData]],
        omega_ranges: List,
        file_ranges: List,
        detector_list: List[Detector],
        output_dir: Path
    ):
        """
        Save detector images to disk.

        Args:
            images: 2D list of detector images
            omega_ranges: List of omega angle ranges
            file_ranges: List of file numbering ranges
            detector_list: List of detectors
            output_dir: Output directory

        C++ Reference:
            ForwardSimulation.cpp:112-126 Output section

        File naming convention:
            {basename}{file_number:0{length}d}.{ext}{detector_number}

        Example:
            Output_0000.tiff0, Output_0001.tiff0, ...
        """
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        print(f"Saving images to {output_dir}...")

        # C++: for(Size_Type nDetNum = 0; nDetNum < oDetectorList.size(); nDetNum++)
        for det_num in range(len(detector_list)):
            # C++: for(Size_Type i = 0; i < vOmegaRangeList.size(); i++)
            for omega_idx in range(len(omega_ranges)):
                # Calculate file index
                # C++: Int nCurrentFileIndex = vFileRangeList[nDetNum].nLow + i
                current_file_index = file_ranges[det_num].low + omega_idx

                # Construct filename
                # C++: tmpSS << oSetupFile.OutFileBasename << InitFileIO::NumToSuffix(...)
                #           << "." << oSetupFile.OutFileExt << nDetNum
                basename = self.config_file.out_file_basename or "Output"
                serial_length = self.config_file.out_file_serial_length or 4
                ext = self.config_file.out_file_ext or "tiff"

                file_number_str = str(current_file_index).zfill(serial_length)
                filename = f"{basename}{file_number_str}.{ext}{det_num}"
                filepath = output_dir / filename

                # Create parent directories if needed
                filepath.parent.mkdir(parents=True, exist_ok=True)

                # Save image
                # C++: oSimData[i][nDetNum].PrintRaster(tmpSS.str())
                image = images[omega_idx][det_num]
                image.save_ascii(str(filepath))

                print(f"  Saved {filename}")

        print(f"Saved {len(omega_ranges) * len(detector_list)} images")

    def __repr__(self):
        """String representation for debugging."""
        num_omegas = len(self.images) if self.images else 0
        num_detectors = len(self.images[0]) if self.images and self.images[0] else 0

        return (
            f"ForwardSimulation("
            f"images={num_omegas}x{num_detectors}, "
            f"beam_energy={self.exp_setup.beam_energy if hasattr(self.exp_setup, 'beam_energy') else 0:.2f} keV)"
        )
