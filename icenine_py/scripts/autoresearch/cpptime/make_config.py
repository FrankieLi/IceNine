"""Write the C++ pilot config from the Python arm A config (cpptime).

  uv run python scripts/autoresearch/cpptime/make_config.py --out DIR/cpp_pilot.config \
      [--sample region_grid.mic] [--basename cpp_pilot]

Drops the three Python-only lines (BFSRevisitRefit, BFSReseed, FusedCost), makes every path
absolute, points SampleFilename at the pilot region grid, sets OutStructureBasename and
SecondsBetweenSave 100000 (no partial file). The search lines are untouched.
"""

import argparse
import re
from pathlib import Path
from typing import List

ICE = Path(__file__).resolve().parents[3]
EX = ICE.parent / "Examples" / "Example2.ManyGrains"
SRC = ICE / "scripts" / "autoresearch" / "ar2" / "key_off_fused.config"
REGION_GRID = ICE / "benchmarks" / "phase_d_pilot" / "region_grid.mic"
PYTHON_ONLY = ("BFSRevisitRefit", "BFSReseed", "FusedCost")
PATH_KEYS = (
    "InfileBasename",
    "DetectorFilename",
    "StructureFilename",
    "RotationRangeFilename",
    "FundamentalZoneFilename",
    "OptimizationFilename",
    "DetectionLimitFilename",
    "OptimizationConstrainFilename",
)


def convert(text: str, sample: Path, out_basename: str, ex: Path = EX) -> str:
    """The C++ config text made from the Python config ``text``."""
    lines: List[str] = []
    seen = set()
    for ln in text.splitlines():
        tok = ln.split()
        key = tok[0] if tok else ""
        if key in PYTHON_ONLY:
            continue
        if key in PATH_KEYS:
            ln = f"{key} {ex / tok[1]}"
        elif key == "SampleFilename":
            ln = f"SampleFilename {sample}"
        elif key == "OutStructureBasename":
            ln = f"OutStructureBasename {out_basename}"
        elif key == "SecondsBetweenSave":
            ln = "SecondsBetweenSave 100000"
        if key:
            seen.add(key)
        lines.append(ln)
    for k in PATH_KEYS + ("SampleFilename", "OutStructureBasename", "SecondsBetweenSave"):
        assert k in seen, f"{k} missing from the source config"
    out = "\n".join(lines) + "\n"
    assert not re.search(r"(?m)^(BFSRevisitRefit|BFSReseed|FusedCost)\b", out)
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--sample", type=Path, default=REGION_GRID)
    ap.add_argument("--basename", default="cpp_pilot")
    ap.add_argument("--src", type=Path, default=SRC)
    a = ap.parse_args()
    a.out.parent.mkdir(parents=True, exist_ok=True)
    a.out.write_text(convert(a.src.read_text(), a.sample.resolve(), a.basename))
    print(f"wrote {a.out}")


if __name__ == "__main__":
    main()
