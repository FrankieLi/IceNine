"""Write the four Phase D reconstruction configs (mc/cma x clean/realistic) into configs/.

Derived from Examples/Example2.ThreeVoxels/ConfigFiles/ReconstructQ8.config (the BFS settings:
MaxQ 8, MaxMCSteps 200, MaxLocalResolution 3, SuccessiveRestarts 2, ...). Paths are relative to
Examples/Example2.ManyGrains (run with that as the working directory).
Usage (from icenine_py/): uv run python scripts/phase_d/make_configs.py
"""

from pathlib import Path

HERE = Path(__file__).resolve().parent
EXAMPLE_CFG = HERE.parents[2] / "Examples" / "Example2.ThreeVoxels" / "ConfigFiles"

# Replace the three data lines and the header of ReconstructQ8.config; keep everything else.
OPTIMIZER = {
    "mc": "# classic optimizer (C++ parity): no LocalOptimizer key\n",
    "cma": "LocalOptimizer cma\nCMASigma0 0.2\nCMAMaxEvals 1000\n",
}
PLACEHOLDER = (
    "# Keys being added to the library (NOT used yet; uncomment when merged):\n"
    "#CMANeighborMaxEvals 250\n#CMARetrySigma0 1.5\n#BFSRefit 1\n"
)


def make(opt: str, variant: str) -> str:
    text = (EXAMPLE_CFG / "ReconstructQ8.config").read_text()
    text = text.replace(
        "#  Reconstruction config — MaxQ=8, for benchmarking single voxel",
        f"#  Phase D full-sample BFS reconstruction ({opt}, {variant} images), MaxQ=8.\n"
        "#  Run with cwd = Examples/Example2.ManyGrains. Python-only (LocalOptimizer keys).",
    )
    text = text.replace(
        "InfileBasename     \t      ScatteringData/3Grains.sim",
        f"InfileBasename     \t      ScatteringData_PhaseD/full/{variant}/500Grains.sim",
    )
    text = text.replace(
        "OutfileBasename\t\t      ScatteringData/3Grains.sim", "OutfileBasename\t\t      None"
    )
    text = text.replace(
        "OutStructureBasename \t  ReconstructedQ8",
        f"OutStructureBasename \t  PhaseD_{opt}_{variant}",
    )
    text = text.replace("SimInput/three_voxels.mic", "SimInput/rand_500grains_1mm_neworient_s0.mic")
    text = text.replace("LazyBFS\n", "LazyBFS\n\n" + OPTIMIZER[opt] + PLACEHOLDER, 1)
    assert "ScatteringData_PhaseD" in text and "neworient" in text and "LazyBFS" in text
    return text


if __name__ == "__main__":
    out = HERE / "configs"
    out.mkdir(exist_ok=True)
    for opt in ("mc", "cma"):
        for variant in ("clean", "realistic"):
            (out / f"ReconstructPhaseD_{opt}_{variant}.config").write_text(make(opt, variant))
            print("wrote", f"ReconstructPhaseD_{opt}_{variant}.config")
