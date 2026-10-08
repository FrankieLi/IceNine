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
    "cma": (
        "LocalOptimizer cma\nCMASigma0 0.2\nCMAMaxEvals 1000\n"
        "CMANeighborMaxEvals 250\nCMARetrySigma0 1.5\n"
    ),
    # CMA without the retry: separates the effect of the retry from that of the optimizer
    "cma_noretry": (
        "LocalOptimizer cma\nCMASigma0 0.2\nCMAMaxEvals 1000\n"
        "CMANeighborMaxEvals 250\nCMARetrySigma0 0\n"
    ),
}
IMAGE_DIR = {
    "clean": "full/clean",
    "realistic": "full/realistic",
    "realistic_q16": "full_q16/realistic",
}
PLACEHOLDER = (
    "# Planned library key (NOT used yet; enable in the mc and cma arms once the library branch\n"
    "# lands): REFIT voxels revisited within the BFS instead of in a post-pass.\n"
    "#BFSRevisitRefit 1\n"
    "# NOTE: CMANeighborMaxEvals, CMARetrySigma0 are library-branch keys (rejected here).\n"
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
        f"InfileBasename     \t      ScatteringData_PhaseD/{IMAGE_DIR[variant]}/500Grains.sim",
    )
    text = text.replace(
        "OutfileBasename\t\t      ScatteringData/3Grains.sim", "OutfileBasename\t\t      None"
    )
    text = text.replace(
        "OutStructureBasename \t  ReconstructedQ8",
        f"OutStructureBasename \t  PhaseD_{opt}_{variant}",
    )
    text = text.replace(
        "SimInput/three_voxels.mic", "SimInput/rand_500grains_1mm_neworient_s0_grid.mic"
    )
    text = text.replace("MaxInitSideLength      0.004000", "MaxInitSideLength      0.009375")
    text = text.replace("MinSideLength          0.004000", "MinSideLength          0.009375")
    text = text.replace("LazyBFS\n", "LazyBFS\n\n" + OPTIMIZER[opt] + PLACEHOLDER, 1)
    assert "ScatteringData_PhaseD" in text and "_grid.mic" in text and "LazyBFS" in text
    assert "0.004" not in text
    return text


if __name__ == "__main__":
    out = HERE / "configs"
    out.mkdir(exist_ok=True)
    for opt in OPTIMIZER:
        for variant in tuple(IMAGE_DIR):
            (out / f"ReconstructPhaseD_{opt}_{variant}.config").write_text(make(opt, variant))
            print("wrote", f"ReconstructPhaseD_{opt}_{variant}.config")
