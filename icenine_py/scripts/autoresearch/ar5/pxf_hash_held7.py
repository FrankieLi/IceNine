"""AR5 step 28: freeze ``held7`` after its build, before any method runs on it.

  uv run python scripts/autoresearch/ar5/pxf_hash_held7.py

Writes ``ar5/pxf/held7_frozen.json``: the sha256 of the built files (grid, truth, grain map, sample
record, stack, stack record; the bytes are hashed, nothing is parsed or analysed) and the counts
printed by ``pxf_build_held7.py`` (seed, voxels, grains). Refuses to overwrite an existing file.
"""

import hashlib
import json
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
ICE = HERE.parents[2]
H7 = ICE / "benchmarks" / "autoresearch" / "ar3" / "held7" / "gradient"
STACK = ICE.parent / "Examples" / "Example2.ManyGrains" / "ScatteringData_PhaseD" / "stacks"
OUT = ICE / "benchmarks" / "autoresearch" / "ar5" / "pxf" / "held7_frozen.json"
FILES = {
    "sample/gradient_grid.mic": H7 / "sample" / "gradient_grid.mic",
    "sample/gradient_truth.mic": H7 / "sample" / "gradient_truth.mic",
    "sample/grainmap.npy": H7 / "sample" / "grainmap.npy",
    "sample/sample.json": H7 / "sample" / "sample.json",
    "held7_stack.json": H7 / "held7_stack.json",
    "stacks/held7_clean.npy": STACK / "held7_clean.npy",
}


def sha(p: Path) -> str:
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for blk in iter(lambda: f.read(1 << 20), b""):
            h.update(blk)
    return h.hexdigest()


def main() -> int:
    if OUT.exists():
        print("refusing: held7_frozen.json exists", file=sys.stderr)
        return 2
    OUT.write_text(
        json.dumps(
            dict(
                sample="held7",
                generator_seed=7,
                n_voxels=24570,
                n_grains=250,
                stack_matches_render=True,
                counts_from="the end of pxf_build_held7.py (exit 0)",
                sha256={k: sha(p) for k, p in FILES.items()},
                frozen_at=time.strftime("%Y-%m-%dT%H:%M:%S%z"),
                note="written before S8, the new J5SD or the plain J5SD ran on held7; no truth, "
                "orientation or map of held7 was analysed",
            ),
            indent=1,
        )
        + "\n"
    )
    print("wrote", OUT)
    return 0


if __name__ == "__main__":
    sys.exit(main())
