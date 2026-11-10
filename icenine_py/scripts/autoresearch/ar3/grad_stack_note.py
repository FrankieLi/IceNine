"""AR3 step 3: how the image stack of the gradient sample was built, and a check of the stack.

  uv run python scripts/autoresearch/ar3/grad_stack_note.py [--build]

The images of tag ``grad`` were rendered with a command of this form (rebuilt from the record
``scripts/phase_d/results/render_grad_nodrop.json``: tag ``grad``, 6 workers, clean frames only,
the default noise variant ``nodrop``):

  uv run python scripts/phase_d/render_full.py --tag grad --workers 6 --variants clean \\
      --mic <icenine_py>/benchmarks/autoresearch/ar3/gradient/sample/gradient_truth.mic

The record does not hold the ``--mic`` path, so the path and the sha256 of the truth map are
recorded here; the path is the one the run scripts use, and the stack check below is what shows
that the frames come from the gradient map rather than the stored one. The clean frames
in ``ScatteringData_PhaseD/grad/clean`` were then packed into the uint8 stack
``ScatteringData_PhaseD/stacks/grad_clean.npy`` with ``write_binary_stack``, with the arguments
that ``scripts/phase_d/memmap_check.py`` uses for the stored images (``--build`` does this when
the stack is missing). The stack and the images are gitignored and not committed.

The check: the stack's non-zero pixel count must equal the render record's ``lit_clean_total``
(and differ from the stored images' count in ``render_full.json``). Writes
``benchmarks/autoresearch/ar3/gradient/grad_stack.json``. Reads the 1.5 GB stack once.
"""

import argparse
import hashlib
import json
from pathlib import Path
from typing import Any, Dict

import ar3_common as C
import gradient_common as G
import numpy as np

STACK = C.EX / "ScatteringData_PhaseD" / "stacks" / G.STACK
FRAMES = C.EX / "ScatteringData_PhaseD" / "grad" / "clean"
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
    rec = json.loads((RES / "render_grad_nodrop.json").read_text())
    old = json.loads((RES / "render_full.json").read_text())
    info = json.loads((G.SAMPLE / "sample.json").read_text())
    truth_sha = sha256(G.TRUTH_MIC)
    assert truth_sha == info["truth_sha256"]
    out: Dict[str, Any] = dict(
        render_command=(
            "uv run python scripts/phase_d/render_full.py --tag grad --workers 6 --variants clean"
            " --mic <icenine_py>/benchmarks/autoresearch/ar3/gradient/sample/gradient_truth.mic"
        ),
        render_command_note="rebuilt from the render record; the --mic path was not recorded",
        render_record="scripts/phase_d/results/render_grad_nodrop.json",
        mic="benchmarks/autoresearch/ar3/gradient/sample/gradient_truth.mic",
        mic_sha256=truth_sha,
        mic_path_recorded_at_render=False,
        stack="Examples/Example2.ManyGrains/ScatteringData_PhaseD/stacks/" + G.STACK,
        stack_build=(
            "write_binary_stack(ScatteringData_PhaseD/grad/clean, '500Grains.sim', 'd', 5, 180, 2,"
            " 2048, 2048, stacks/grad_clean.npy)"
        ),
        stack_built_now=built,
        stack_shape=list(st.shape),
        stack_nonzero=nz,
        render_lit_clean_total=int(rec["lit_clean_total"]),
        render_spots_total=int(rec["spots_total"]),
        stored_lit_clean_total=int(old["lit_clean_total"]),
        stored_spots_total=int(old["spots_total"]),
        stack_matches_render=bool(nz == int(rec["lit_clean_total"])),
        stack_differs_from_stored=bool(nz != int(old["lit_clean_total"])),
    )
    (G.GRAD / "grad_stack.json").write_text(json.dumps(out, indent=1) + "\n")
    print(json.dumps(out, indent=1))
    assert out["stack_matches_render"] and out["stack_differs_from_stored"]


if __name__ == "__main__":
    main()
