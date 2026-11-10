"""Task A step 0: where do the triangles that produce no pixel lie? Reads the C++ driver output
(DIR/cpp_<voxel>.txt) and counts, over the complete-cost arm, the peaks that have a lit detector
with no overlap and no produced pixel on the whole peak, and how many of them have a vertex
coordinate at the last image row or column (2047).

  uv run python scripts/eval_reduction/pixelless_check.py --dir DIR \
      --out benchmarks/eval_reduction/A/step0
"""

import argparse
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import parity_common as PC  # noqa: E402


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--dir", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    total = at_edge = 0
    for v in PC.PARITY_VOXELS:
        cpp = PC.parse_cpp(Path(a.dir) / f"cpp_{v}.txt")
        for rec in cpp["recs"].values():
            for p in rec["arms"]["EXF"]["peaks"]:
                if p["nomatch"] or p["pixon"] != 0:
                    continue
                for d in p["dets"]:
                    if d["hit"] == 3 and d["lit"] and not d["ov"]:
                        total += 1
                        at_edge += int((d["xy"] == 2047).any())
    res = dict(peaks_lit_without_pixels=total, of_which_with_vertex_at_2047=at_edge)
    (Path(a.out) / "pixelless_triangles.json").write_text(json.dumps(res) + "\n")
    print(res)


if __name__ == "__main__":
    main()
