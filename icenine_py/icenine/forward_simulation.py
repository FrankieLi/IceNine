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
import math
import time
import numpy as np
import torch

from .config_file import ConfigFile
from .experiment_setup import XDMExperimentSetup
from .simulation import Simulation, PeakInfo
from .sample import Sample
from .detector import Detector
from .image_data import ImageData
from .peak_filters import XDMEtaAcceptFn, batch_eta_filter
from .constants import KEV_OVER_HBAR_C_IN_ANG
from .diffraction_core import (
    get_scattering_omegas_torch,
    get_reflection_vector,
)
from .geometry import Ray


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
        output_dir: Optional[Path] = None,
        batched: bool = False,
        batch_size: Optional[int] = None,
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
        if batched:
            self._simulate_peaks_batched(
                self.images,
                detector_list,
                sample,
                range_map,
                batch_size=batch_size,
            )
        else:
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

        Optimized version:
        - Batches Bragg solving per-voxel (all reflections at once)
        - Uses functional rotations (no sample mutation/restore)
        - Pure Python floats in inner loop (avoids torch scalar overhead)
        - Pre-computes all detector geometry and reciprocal vectors

        C++ Reference:
            ForwardSimulation.cpp:200-284 SimulatePeaks
        """
        eta_limit = self.exp_setup.get_eta_limit()
        wavenumber = KEV_OVER_HBAR_C_IN_ANG * self.exp_setup.beam_energy
        beam_energy = self.exp_setup.beam_energy
        beam_deflection = self.exp_setup.get_beam_deflection_chi_laue()
        beam_dir_t = self.simulator.beam_direction
        bd0, bd1, bd2 = beam_dir_t[0].item(), beam_dir_t[1].item(), beam_dir_t[2].item()
        structure_list = sample.get_structure_list()
        mic = sample.get_mic()
        voxel_list = mic.voxels
        num_detectors = len(detector_list)

        # Pre-compute detector geometry as plain Python floats
        det_geom = []
        for det in detector_list:
            plane = det.detector_plane
            pn = plane.normal
            dp = det._position
            dco = det._lab_frame_coord_origin
            bj = det._lab_frame_basis_j
            bk = det._lab_frame_basis_k
            det_geom.append((
                pn[0].item(), pn[1].item(), pn[2].item(), plane.d.item(),
                dp[0].item(), dp[1].item(), dp[2].item(),
                dco[0].item(), dco[1].item(), dco[2].item(),
                bj[0].item(), bj[1].item(), bj[2].item(),
                bk[0].item(), bk[1].item(), bk[2].item(),
                det.pixel_half_width, det.pixel_half_height,
                det.pixel_width, det.pixel_height,
            ))

        # Pre-compute reciprocal vector tensors per phase (torch for batched Bragg)
        # and plain float arrays for inner loop
        phase_data_torch = {}
        phase_data_float = {}
        for phase_idx, structure in enumerate(structure_list):
            recp_vecs = structure.get_reflection_vectors()
            if not recp_vecs:
                continue
            g_hkl_batch = torch.stack(
                [torch.from_numpy(rv.q_vec).float() for rv in recp_vecs]
            )
            g_mag_batch = torch.tensor(
                [rv.q_mag for rv in recp_vecs], dtype=torch.float32
            )
            intensities_list = [rv.intensity for rv in recp_vecs]
            sin_theta = g_mag_batch / (2.0 * wavenumber)
            sin_2theta = torch.sin(2.0 * torch.asin(sin_theta))
            sin_2theta_list = sin_2theta.tolist()

            phase_data_torch[phase_idx] = (g_hkl_batch, g_mag_batch)
            phase_data_float[phase_idx] = (intensities_list, sin_2theta_list)

        # Get base sample-to-lab matrix as plain floats
        bm = sample.sample_to_lab_matrix
        br00, br01, br02 = bm[0, 0].item(), bm[0, 1].item(), bm[0, 2].item()
        br10, br11, br12 = bm[1, 0].item(), bm[1, 1].item(), bm[1, 2].item()
        br20, br21, br22 = bm[2, 0].item(), bm[2, 1].item(), bm[2, 2].item()
        bt0, bt1, bt2 = bm[0, 3].item(), bm[1, 3].item(), bm[2, 3].item()

        _cos = math.cos
        _sin = math.sin
        _sqrt = math.sqrt
        _atan2 = math.atan2
        _fabs = math.fabs

        print(f"Simulating {len(voxel_list)} voxels...")

        for voxel_count, voxel in enumerate(voxel_list):
            if (voxel_count + 1) % 10000 == 0:
                print(f"  Voxel {voxel_count + 1}/{len(voxel_list)}")

            phase_index = voxel.phase
            if phase_index not in phase_data_torch:
                continue

            g_hkl_batch, g_mag_batch = phase_data_torch[phase_index]
            intensities_list, sin_2theta_list = phase_data_float[phase_index]

            # Transform all reciprocal vectors to lab frame at once (torch, batched)
            voxel_orientation = torch.from_numpy(voxel.orientation).float()
            g_lab_batch = (voxel_orientation @ g_hkl_batch.T).T  # (N_refl, 3)

            # Batch Bragg solving for ALL reflections at once
            bragg_result = get_scattering_omegas_torch(
                g_lab_batch, g_mag_batch, beam_energy, beam_deflection
            )

            # Extract observable results to plain Python lists
            obs_mask = bragg_result.observable
            if not obs_mask.any():
                continue

            obs_indices = torch.where(obs_mask)[0].tolist()
            g_lab_np = g_lab_batch.numpy()
            omega1_np = bragg_result.omega1.numpy()
            omega2_np = bragg_result.omega2.numpy()

            # Pre-compute scattering directions as plain floats
            obs_data = []
            for idx in obs_indices:
                gx, gy, gz = float(g_lab_np[idx, 0]), float(g_lab_np[idx, 1]), float(g_lab_np[idx, 2])
                gnorm = _sqrt(gx * gx + gy * gy + gz * gz)
                inv_gnorm = 1.0 / gnorm if gnorm > 0 else 0.0
                sdx, sdy, sdz = gx * inv_gnorm, gy * inv_gnorm, gz * inv_gnorm
                obs_data.append((
                    sdx, sdy, sdz,
                    float(omega1_np[idx]), float(omega2_np[idx]),
                    intensities_list[idx], sin_2theta_list[idx]
                ))

            # Pre-compute voxel vertices as plain floats
            x = float(voxel.position[0])
            y = float(voxel.position[1])
            z = float(voxel.position[2])
            s = float(voxel.side_length)
            sqrt3 = 1.7320508075688772
            if voxel.points_up:
                verts = (
                    (x, y, z),
                    (x + s, y, z),
                    (x + s * 0.5, y + s * 0.5 * sqrt3, z),
                )
            else:
                verts = (
                    (x, y, z),
                    (x + s * 0.5, y - s * 0.5 * sqrt3, z),
                    (x + s, y, z),
                )

            # Process all observable peaks
            for sdx, sdy, sdz, omega1, omega2, form_intensity, sin_2theta in obs_data:
                for omega in (omega1, omega2):
                    omega_index = range_map.angle_to_wedge_index(omega)
                    if omega_index is None:
                        continue

                    # Compute Rz(omega) @ base_rotation as plain floats
                    cos_w = _cos(omega)
                    sin_w = _sin(omega)
                    r00 = cos_w * br00 - sin_w * br10
                    r01 = cos_w * br01 - sin_w * br11
                    r02 = cos_w * br02 - sin_w * br12
                    r10 = sin_w * br00 + cos_w * br10
                    r11 = sin_w * br01 + cos_w * br11
                    r12 = sin_w * br02 + cos_w * br12
                    r20 = br20
                    r21 = br21
                    r22 = br22

                    # Transform scattering direction to lab frame
                    lnx = r00 * sdx + r01 * sdy + r02 * sdz
                    lny = r10 * sdx + r11 * sdy + r12 * sdz
                    lnz = r20 * sdx + r21 * sdy + r22 * sdz

                    # Reflection: r_out = beam - 2*(beam·n)*n
                    dot_bn = bd0 * lnx + bd1 * lny + bd2 * lnz
                    two_dot = 2.0 * dot_bn
                    rdx = bd0 - two_dot * lnx
                    rdy = bd1 - two_dot * lny
                    rdz = bd2 - two_dot * lnz

                    # Peak filter: eta acceptance + intensity
                    rd_norm = _sqrt(rdx * rdx + rdy * rdy + rdz * rdz)
                    inv_rd_norm = 1.0 / rd_norm if rd_norm > 0 else 0.0
                    rdy_n = rdy * inv_rd_norm
                    rdz_n = rdz * inv_rd_norm
                    eta = _atan2(_fabs(rdy_n), _fabs(rdz_n))

                    if eta >= eta_limit:
                        continue

                    sin_eta = _sin(eta)
                    intensity = form_intensity / (_fabs(sin_eta) * sin_2theta + 1e-10)

                    # Project 3 vertices onto all detectors
                    all_hit = True
                    projected_pixels = [None] * (num_detectors * 3)

                    for vi in range(3):
                        vx, vy, vz = verts[vi]
                        lab_vx = r00 * vx + r01 * vy + r02 * vz + bt0
                        lab_vy = r10 * vx + r11 * vy + r12 * vz + bt1
                        lab_vz = r20 * vx + r21 * vy + r22 * vz + bt2

                        for di in range(num_detectors):
                            dg = det_geom[di]
                            # Ray-plane intersection
                            denom = dg[0] * rdx + dg[1] * rdy + dg[2] * rdz
                            if _fabs(denom) < 1e-8:
                                all_hit = False
                                break
                            numer = -(dg[0] * lab_vx + dg[1] * lab_vy + dg[2] * lab_vz + dg[3])
                            t = numer / denom
                            if t <= 0:
                                all_hit = False
                                break

                            # Intersection point
                            ix = lab_vx + t * rdx
                            iy = lab_vy + t * rdy
                            iz = lab_vz + t * rdz

                            # Lab to detector pixel coordinates
                            px = ix - dg[4] - dg[7]
                            py = iy - dg[5] - dg[8]
                            pz = iz - dg[6] - dg[9]
                            j_coord = px * dg[10] + py * dg[11] + pz * dg[12]
                            k_coord = px * dg[13] + py * dg[14] + pz * dg[15]
                            col = (j_coord + dg[16]) / dg[18]
                            row = (k_coord + dg[17]) / dg[19]

                            projected_pixels[di * 3 + vi] = (col, row)

                        if not all_hit:
                            break

                    if not all_hit:
                        continue

                    # Rasterize onto each detector
                    for di in range(num_detectors):
                        base = di * 3
                        v0 = torch.tensor(projected_pixels[base])
                        v1 = torch.tensor(projected_pixels[base + 1])
                        v2 = torch.tensor(projected_pixels[base + 2])
                        images[omega_index][di].add_triangle_scanline(
                            v0, v1, v2, intensity
                        )

    # ------------------------------------------------------------------
    # Batched (differentiable) forward simulation pipeline
    # ------------------------------------------------------------------

    def _simulate_peaks_batched(
        self,
        images: List[List[ImageData]],
        detector_list: List[Detector],
        sample: Sample,
        range_map,
        batch_size: Optional[int] = None,
    ):
        """Fully batched forward simulation using torch tensor operations.

        Processes all voxels via large tensor ops (stages 1-5 are differentiable).
        Only the final rasterization (stage 6) is sequential / non-differentiable.

        Args:
            images: 2D list [omega_idx][det_idx] of ImageData to accumulate into
            detector_list: list of Detector objects
            sample: Sample with voxels and crystal structures
            range_map: SimulationRange for omega-to-wedge lookup
            batch_size: Max voxels per chunk (None = all at once)
        """
        # ---- Stage 0: data preparation ----
        t0 = time.time()
        (
            orientations, vertices, phase_indices,
            phase_data, base_rot, translation,
            det_tensors, beam_dir, beam_energy, beam_deflection,
            eta_limit, omega_lookup,
        ) = self._prepare_batched_data(sample, detector_list, range_map)

        V = orientations.shape[0]
        if batch_size is None:
            batch_size = V

        num_detectors = len(detector_list)
        t_prep = time.time() - t0
        print(f"  Batched data prep: {t_prep:.2f}s  ({V} voxels)")

        # ---- Chunk loop ----
        total_rasterized = 0
        for chunk_start in range(0, V, batch_size):
            chunk_end = min(chunk_start + batch_size, V)
            chunk_orient = orientations[chunk_start:chunk_end]
            chunk_verts = vertices[chunk_start:chunk_end]
            chunk_phases = phase_indices[chunk_start:chunk_end]

            for phase_idx, pd in phase_data.items():
                phase_mask = chunk_phases == phase_idx
                if not phase_mask.any():
                    continue

                orient_p = chunk_orient[phase_mask]       # (Vp, 3, 3)
                verts_p = chunk_verts[phase_mask]          # (Vp, 3, 3)
                g_hkl = pd["g_hkl"]                        # (R, 3)
                g_mag = pd["g_mag"]                         # (R,)
                intensities = pd["intensities"]             # (R,)
                sin_2theta = pd["sin_2theta"]               # (R,)
                Vp = orient_p.shape[0]
                R = g_hkl.shape[0]

                # ---- Stage 1: batched Bragg solving ----
                g_lab = torch.bmm(
                    orient_p,
                    g_hkl.unsqueeze(0).expand(Vp, -1, -1).transpose(1, 2),
                ).transpose(1, 2)  # (Vp, R, 3)

                g_lab_flat = g_lab.reshape(Vp * R, 3)
                g_mag_flat = g_mag.unsqueeze(0).expand(Vp, R).reshape(Vp * R)

                bragg = get_scattering_omegas_torch(
                    g_lab_flat, g_mag_flat, beam_energy, beam_deflection
                )
                observable = bragg.observable.reshape(Vp, R)
                omega1 = bragg.omega1.reshape(Vp, R)
                omega2 = bragg.omega2.reshape(Vp, R)

                # ---- Stage 2: omega filter + compaction ----
                omegas = torch.stack([omega1, omega2], dim=2)  # (Vp, R, 2)
                obs_exp = observable.unsqueeze(2).expand_as(omegas)

                # Vectorized omega-to-wedge lookup (float64 to match serial precision)
                idx_tensor, ol_low, ol_width, ol_n = omega_lookup
                omega_flat = omegas.reshape(-1)
                f_vals = (omega_flat.double() - ol_low) / ol_width
                bin_idx = f_vals.long()
                # Must check f_vals >= 0 (not just bin_idx >= 0) because
                # .long() truncates -0.001 to 0, which would falsely pass
                bin_valid = (f_vals >= 0) & (bin_idx < ol_n)
                bin_idx_clamped = bin_idx.clamp(0, ol_n - 1)
                wedge_idx = idx_tensor[bin_idx_clamped]
                omega_valid = bin_valid & (wedge_idx >= 0)

                combined = obs_exp.reshape(-1) & omega_valid
                valid_flat = torch.where(combined)[0]
                N_valid = valid_flat.shape[0]
                if N_valid == 0:
                    continue

                # Index arrays
                voxel_idx = valid_flat // (R * 2)
                remainder = valid_flat % (R * 2)
                refl_idx = remainder // 2

                valid_omegas = omega_flat[valid_flat]
                valid_wedge_idx = wedge_idx[valid_flat]

                # ---- Stage 3: batched rotation + reflection ----
                cos_w = torch.cos(valid_omegas)
                sin_w = torch.sin(valid_omegas)
                zeros = torch.zeros_like(cos_w)
                ones = torch.ones_like(cos_w)

                Rz = torch.stack([
                    cos_w, -sin_w, zeros,
                    sin_w,  cos_w, zeros,
                    zeros,  zeros, ones,
                ], dim=1).reshape(N_valid, 3, 3)

                full_rot = torch.bmm(Rz, base_rot.unsqueeze(0).expand(N_valid, -1, -1))

                # Scattering directions: normalized g_lab vectors
                g_lab_valid = g_lab[voxel_idx, refl_idx]  # (N_valid, 3)
                g_norms = torch.norm(g_lab_valid, dim=1, keepdim=True)
                scat_dir = g_lab_valid / (g_norms + 1e-10)

                # Transform to lab frame
                lab_normal = torch.bmm(
                    full_rot, scat_dir.unsqueeze(2)
                ).squeeze(2)  # (N_valid, 3)

                # Reflection: r_out = beam - 2*(beam·n)*n
                beam_exp = beam_dir.unsqueeze(0).expand(N_valid, -1)
                dot_bn = torch.sum(beam_exp * lab_normal, dim=1, keepdim=True)
                ray_dir = beam_exp - 2.0 * dot_bn * lab_normal  # (N_valid, 3)

                # ---- Stage 4: batched eta filter ----
                valid_intensities = intensities[refl_idx]
                valid_sin2theta = sin_2theta[refl_idx]
                eta_accept, intensity = batch_eta_filter(
                    ray_dir, eta_limit, valid_intensities, valid_sin2theta
                )

                # ---- Stage 5: batched vertex projection ----
                valid_verts = verts_p[voxel_idx]  # (N_valid, 3, 3)

                # Transform vertices to lab frame
                lab_verts = torch.bmm(
                    full_rot,
                    valid_verts.transpose(1, 2),
                ).transpose(1, 2) + translation.unsqueeze(0).unsqueeze(1)
                # lab_verts: (N_valid, 3_verts, 3_xyz)

                # Per-detector ray-plane intersection + pixel projection
                all_det_hit = eta_accept.clone()
                pixel_coords_list = []

                for di in range(num_detectors):
                    dn = det_tensors["normals"][di]       # (3,)
                    dd = det_tensors["d"][di]              # scalar
                    dp = det_tensors["pos"][di]            # (3,)
                    do = det_tensors["origin"][di]         # (3,)
                    dbj = det_tensors["basis_j"][di]       # (3,)
                    dbk = det_tensors["basis_k"][di]       # (3,)
                    phw = det_tensors["pixel_params"][di, 0].item()
                    phh = det_tensors["pixel_params"][di, 1].item()
                    pw = det_tensors["pixel_params"][di, 2].item()
                    ph = det_tensors["pixel_params"][di, 3].item()

                    # denom = normal · ray_dir  → (N_valid,)
                    denom = torch.sum(dn.unsqueeze(0) * ray_dir, dim=1)
                    denom_ok = torch.abs(denom) > 1e-8

                    # numer = -(normal · vert + d)  → (N_valid, 3)
                    numer = -(
                        torch.sum(
                            dn.unsqueeze(0).unsqueeze(0) * lab_verts, dim=2
                        ) + dd
                    )

                    # t = numer / denom  → (N_valid, 3)
                    safe_denom = denom.clone()
                    safe_denom[~denom_ok] = 1.0
                    t = numer / safe_denom.unsqueeze(1)

                    t_ok = t > 0  # (N_valid, 3)
                    hit_d = denom_ok.unsqueeze(1) & t_ok
                    all_verts_hit = hit_d.all(dim=1)  # (N_valid,)
                    all_det_hit = all_det_hit & all_verts_hit

                    # Intersection points  → (N_valid, 3, 3)
                    intersect = lab_verts + t.unsqueeze(2) * ray_dir.unsqueeze(1)

                    # Lab-to-pixel: project onto detector basis
                    offset = dp + do  # position + coord_origin
                    pixel_loc = intersect - offset.unsqueeze(0).unsqueeze(0)
                    j_coord = torch.sum(pixel_loc * dbj.unsqueeze(0).unsqueeze(0), dim=2)
                    k_coord = torch.sum(pixel_loc * dbk.unsqueeze(0).unsqueeze(0), dim=2)

                    col = (j_coord + phw) / pw  # (N_valid, 3)
                    row = (k_coord + phh) / ph  # (N_valid, 3)

                    pixel_coords_list.append(
                        torch.stack([col, row], dim=2)  # (N_valid, 3, 2)
                    )

                # ---- Stage 6: rasterization (sequential) ----
                final_valid = torch.where(all_det_hit)[0]
                n_raster = final_valid.shape[0]
                total_rasterized += n_raster

                # Detach and convert to numpy for fast sequential access
                wedge_np = valid_wedge_idx[final_valid].numpy()
                intensity_np = intensity[final_valid].detach().numpy()
                px_np = [pc[final_valid].detach().numpy() for pc in pixel_coords_list]

                for i in range(n_raster):
                    oi = int(wedge_np[i])
                    inten = float(intensity_np[i])
                    for di in range(num_detectors):
                        v0c, v0r = float(px_np[di][i, 0, 0]), float(px_np[di][i, 0, 1])
                        v1c, v1r = float(px_np[di][i, 1, 0]), float(px_np[di][i, 1, 1])
                        v2c, v2r = float(px_np[di][i, 2, 0]), float(px_np[di][i, 2, 1])
                        images[oi][di].add_triangle_scanline(
                            (v0c, v0r), (v1c, v1r), (v2c, v2r), inten
                        )

            if chunk_end < V:
                print(
                    f"  Chunk {chunk_start}-{chunk_end}/{V},"
                    f" rasterized so far: {total_rasterized}"
                )

        t_total = time.time() - t0
        print(f"  Batched simulation: {t_total:.2f}s, {total_rasterized} triangles rasterized")

    def _prepare_batched_data(self, sample, detector_list, range_map):
        """Stage 0: collect all data into contiguous tensors."""
        eta_limit = self.exp_setup.get_eta_limit()
        wavenumber = KEV_OVER_HBAR_C_IN_ANG * self.exp_setup.beam_energy
        beam_energy = self.exp_setup.beam_energy
        beam_deflection = self.exp_setup.get_beam_deflection_chi_laue()
        beam_dir = self.simulator.beam_direction.float()

        structure_list = sample.get_structure_list()
        mic = sample.get_mic()
        voxel_list = mic.voxels
        V = len(voxel_list)
        sqrt3 = 1.7320508075688772

        # Voxel tensors
        orient_list = []
        vert_list = []
        phase_list = []
        for voxel in voxel_list:
            orient_list.append(torch.from_numpy(voxel.orientation).float())
            x, y, z = float(voxel.position[0]), float(voxel.position[1]), float(voxel.position[2])
            s = float(voxel.side_length)
            if voxel.points_up:
                v = [[x, y, z], [x + s, y, z], [x + s * 0.5, y + s * 0.5 * sqrt3, z]]
            else:
                v = [[x, y, z], [x + s * 0.5, y - s * 0.5 * sqrt3, z], [x + s, y, z]]
            vert_list.append(v)
            phase_list.append(voxel.phase)

        orientations = torch.stack(orient_list)                         # (V, 3, 3)
        vertices = torch.tensor(vert_list, dtype=torch.float32)         # (V, 3, 3)
        phase_indices = torch.tensor(phase_list, dtype=torch.long)      # (V,)

        # Phase data
        phase_data = {}
        for phase_idx, structure in enumerate(structure_list):
            recp_vecs = structure.get_reflection_vectors()
            if not recp_vecs:
                continue
            g_hkl = torch.stack([torch.from_numpy(rv.q_vec).float() for rv in recp_vecs])
            g_mag = torch.tensor([rv.q_mag for rv in recp_vecs], dtype=torch.float32)
            intensities_t = torch.tensor([rv.intensity for rv in recp_vecs], dtype=torch.float32)
            sin_theta = g_mag / (2.0 * wavenumber)
            sin_2theta_t = torch.sin(2.0 * torch.asin(sin_theta))
            phase_data[phase_idx] = {
                "g_hkl": g_hkl,
                "g_mag": g_mag,
                "intensities": intensities_t,
                "sin_2theta": sin_2theta_t,
            }

        # Base rotation (3x3 part of sample_to_lab)
        bm = sample.sample_to_lab_matrix
        base_rot = bm[:3, :3].float()
        translation = bm[:3, 3].float()

        # Detector tensors
        det_tensors = {
            "normals": [], "d": [], "pos": [], "origin": [],
            "basis_j": [], "basis_k": [], "pixel_params": [],
        }
        for det in detector_list:
            plane = det.detector_plane
            det_tensors["normals"].append(plane.normal.float())
            det_tensors["d"].append(plane.d.float())
            det_tensors["pos"].append(det._position.float())
            det_tensors["origin"].append(det._lab_frame_coord_origin.float())
            det_tensors["basis_j"].append(det._lab_frame_basis_j.float())
            det_tensors["basis_k"].append(det._lab_frame_basis_k.float())
            det_tensors["pixel_params"].append(torch.tensor([
                det.pixel_half_width, det.pixel_half_height,
                det.pixel_width, det.pixel_height,
            ], dtype=torch.float32))
        det_tensors["pixel_params"] = torch.stack(det_tensors["pixel_params"])

        # Omega lookup
        omega_lookup = range_map.to_lookup_tensor()

        return (
            orientations, vertices, phase_indices,
            phase_data, base_rot, translation,
            det_tensors, beam_dir, beam_energy, beam_deflection,
            eta_limit, omega_lookup,
        )

    def _get_voxel_vertices(self, voxel) -> torch.Tensor:
        """
        Get triangular vertices for voxel projection.

        Computes the 3 vertices of the equilateral triangle voxel,
        matching the C++ implementation in MicIO.h lines 300-315.

        C++ Reference:
            XDM++/libXDM/MicIO.h lines 300-315
        """
        x = float(voxel.position[0])
        y = float(voxel.position[1])
        z = float(voxel.position[2])
        s = float(voxel.side_length)

        if voxel.points_up:
            vertices = torch.tensor([
                [x,           y,                          z],
                [x + s,       y,                          z],
                [x + s / 2.0, y + s / 2.0 * math.sqrt(3.0), z],
            ], dtype=torch.float32)
        else:
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
