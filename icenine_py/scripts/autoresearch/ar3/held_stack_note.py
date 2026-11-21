"""AR3 step 11: how the image stack of the sample ``held`` was built, and a check of the stack.

  AR3_SAMPLE=held uv run python scripts/autoresearch/ar3/start_stack_note.py [--build]

The images of tag ``held`` were rendered as those of ``mid`` were (``start_stack_note.py``;
1 worker, clean frames only; the record is ``scripts/phase_d/results/render_held_nodrop.json``):

  uv run python scripts/phase_d/render_full.py --tag held --workers 1 --variants clean \\
      --mic <absolute path of benchmarks/autoresearch/ar3/held/gradient/sample/gradient_truth.mic>

The clean frames in ``ScatteringData_PhaseD/held/clean`` are packed into the uint8 stack
``ScatteringData_PhaseD/stacks/held_clean.npy`` with ``write_binary_stack``, with the arguments of
``big_stack_note.py``. The stack and the images are gitignored and not committed. The check: the
stack's non-zero pixel count equals the render record's ``lit_clean_total``. Writes
``benchmarks/autoresearch/ar3/held/gradient/held_stack.json``.
"""

import argparse
import hashlib
import json
from pathlib import Path
from typing import Any, Dict

import ar3_common as C
import gradient_common as G
import numpy as np

assert C.SAMPLE_KIND in ("held", "held2"), "run with AR3_SAMPLE=held or held2"
K = C.SAMPLE_KIND  # AR5 step 12: the same note for held2
STACK = C.EX / "ScatteringData_PhaseD" / "stacks" / G.STACK
FRAMES = C.EX / "ScatteringData_PhaseD" / K / "clean"
RES = C.ICE / "scripts" / "phase_d" / "results"


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 24), b""):
            h.update(chunk)
    return h.hexdigest()


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--build", action="store_true", help="build the stack if it is missing")
    a = ap.parse_args()
    built = False
    if not STACK.exists():
        if not a.build:
            raise SystemExit(f"{STACK} is missing; run with --build")
        from icenine.experimental_data import write_binary_stack

        write_binary_stack(FRAMES, "500Grains.sim", "d", 5, 180, 2, 2048, 2048, STACK)
        built = True
    st = np.load(STACK, mmap_mode="r")
    nz = 0
    for i in range(st.shape[0]):
        nz += int(np.count_nonzero(st[i]))
    rec = json.loads((RES / f"render_{K}_nodrop.json").read_text())
    info = json.loads((G.SAMPLE / "sample.json").read_text())
    truth_sha = sha256(G.TRUTH_MIC)
    assert truth_sha == info["truth_sha256"]
    out: Dict[str, Any] = dict(
        render_command=(
            f"uv run python scripts/phase_d/render_full.py --tag {K} --workers 1 --variants clean"
            f" --mic <absolute path of benchmarks/autoresearch/ar3/{K}/gradient/sample/"
            "gradient_truth.mic>"
        ),
        render_record=f"scripts/phase_d/results/render_{K}_nodrop.json",
        mic=f"benchmarks/autoresearch/ar3/{K}/gradient/sample/gradient_truth.mic",
        mic_sha256=truth_sha,
        stack="Examples/Example2.ManyGrains/ScatteringData_PhaseD/stacks/" + G.STACK,
        stack_build=(
            f"write_binary_stack(ScatteringData_PhaseD/{K}/clean, '500Grains.sim', 'd', 5, 180, 2,"
            f" 2048, 2048, stacks/{K}_clean.npy)"
        ),
        stack_built_now=built,
        stack_shape=list(st.shape),
        stack_nonzero=nz,
        render_lit_clean_total=int(rec["lit_clean_total"]),
        render_spots_total=int(rec["spots_total"]),
        render_wall_s=float(rec["wall_s"]),
        render_timing_label=rec["timing_label"],
        stack_matches_render=bool(nz == int(rec["lit_clean_total"])),
    )
    (G.GRAD / f"{K}_stack.json").write_text(json.dumps(out, indent=1) + "\n")
    print(json.dumps(out, indent=1))
    assert out["stack_matches_render"]


if __name__ == "__main__":
    main()
