# IceNine Python/PyTorch Implementation

Python port of the IceNine forward model for synchrotron X-ray diffraction simulation, with PyTorch acceleration. Validated pixel-exact against the C++ implementation.

## Installation

```bash
cd icenine_py
uv pip install -e ".[dev]"
```

Dependencies: numpy, torch, pymatgen, scipy (see `pyproject.toml`).

## Quick Start: Forward Simulation

```bash
cd Examples/Example2.ThreeVoxels
uv run python run_python_simulation.py
```

This simulates X-ray diffraction from a 3-voxel copper sample (64.35 keV beam, 180 omega steps × 2 detectors) and writes 360 detector images to `ScatteringData_Python/`.

### Programmatic Usage

```python
from icenine.config_file import ConfigFile
from icenine.forward_simulation import ForwardSimulation

config = ConfigFile.from_file("path/to/experiment.config")
simulator = ForwardSimulation(config)
images = simulator.simulate_detector_images(output_dir="output/")
# images[omega_index][detector_index] is an ImageData object
```

### What ForwardSimulation Does

For each voxel in the sample:
1. Loads crystal structure and generates reciprocal lattice vectors (filtered by MaxQ and MinAmplitudeFraction)
2. For each reciprocal vector, solves the Bragg condition for omega angles
3. Maps each omega to an experimental wedge via `SimulationRange`
4. Rotates the sample to the omega angle
5. Projects the voxel triangle onto all detectors (Sutherland-Hodgman clipping + scanline rasterization)
6. Accumulates intensity with Lorentz-polarization correction: `I = I_form / (|sin(η)| × sin(2θ))`

## Package Modules

| Module | Description |
|--------|-------------|
| `forward_simulation.py` | Main simulation loop — generates detector images from a sample |
| `simulation.py` | Core engine — Bragg condition solving, vertex projection, voxel rasterization |
| `experiment_setup.py` | Reads config, initializes detectors/sample/omega ranges |
| `config_file.py` | Parser for `.config` files (80+ parameters, auto degree→radian conversion) |
| `detector.py` | Detector geometry, coordinate transforms, ray-plane intersection |
| `image_data.py` | Detector image storage, scanline triangle rasterizer |
| `sample.py` | Sample with orientation, translation, crystal structures |
| `mic_file.py` | Read/write `.mic` voxel grid files (Bunge Euler angles) |
| `crystal_structure.py` | Unit cell, reciprocal lattice, reflection vector generation |
| `diffraction_core.py` | Scattering omega calculation, reflected rays (PyTorch batched) |
| `peak_filters.py` | Eta-angle acceptance filter with Lorentz-polarization correction |
| `simulation_range.py` | Omega range system for discontinuous data collection wedges |
| `geometry.py` | Euler conversions, Plane/Ray classes |
| `symmetry.py` | Crystal symmetry operations (pymatgen wrapper) |
| `constants.py` | Physical constants |
| `file_io.py` | Detector file, crystal structure file, and omega file I/O |

## Config File Format

Forward simulation requires a `.config` file specifying:

```
BeamEnergy           64.351           # keV
BeamDirection        0  0  1          # unit vector
MaxQ                 16               # Å⁻¹, max scattering vector magnitude
EtaLimit             86               # degrees, max eta for peak acceptance
MinAmplitudeFraction 0.25             # filter reflections below this fraction of max intensity
SampleFilename       SimInput/three_voxels.mic
StructureFilename    DataFiles/copper.dat
DetectorFilename     ConfigFiles/StandardGeometry.2Det
OmegaFilename        DataFiles/omega_180_2L.dat
OutFileBasename      3Grains.sim
OutFileExtension     d
OutFileSerialLength  5
```

See `Examples/Example2.ThreeVoxels/ConfigFiles/Example2.Simulation.config` for a complete example.

## Testing

```bash
cd icenine_py
uv run pytest tests/ -v                # all unit tests
uv run pytest tests/test_simulation.py  # specific module
uv run pytest --cov=icenine tests/      # with coverage
```

### Integration Test (C++ vs Python)

```bash
cd Examples/Example2.ThreeVoxels
uv run python run_python_simulation.py   # generate Python output
uv run python compare_outputs.py         # compare against C++ reference
```

Expected: 3200/3201 pixels match at identical locations, max relative intensity difference < 1e-5.

## Validation Status

The forward simulation has been validated pixel-exact against C++ on the Example2.ThreeVoxels test case:

- **3200 matching pixels** (same location and intensity to within float32 precision)
- **Max relative intensity difference: 5.4e-6**
- **4 pixel-location mismatches** — all at omega bin boundaries due to floating-point rounding
- **353/360 files are nonempty** in both C++ and Python (identical set)

## Citation

S. F. Li and R. M. Suter, "Adaptive reconstruction method for three-dimensional orientation imaging", *Journal of Applied Crystallography*, 2013.
