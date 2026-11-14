"""AR5 step 12: build and render the new test sample ``held2`` (``held`` with generator seed 2).

  uv run python scripts/autoresearch/ar5/stallhit_make_held2.py [--dry-run] [--resume]

``held2`` is ``held`` (AR3 step 11: 250 Voronoi grains on the 24,570-voxel gradient grid, gradients
0 to 4 degrees, clean images, Q-max 8) built by the same script with the generator seed changed
from 1 to 2: a new sample for a carry-over test, not an extra random seed for a comparison. Three
stages, each the command of ``held``'s own build, with ``AR3_SAMPLE=held2``:

  build   ``scripts/autoresearch/ar3/held_make_sample.py``    -> ar3/held2/gradient/sample/
  render  ``scripts/phase_d/render_full.py --tag held2 --workers 1 --variants clean --mic <truth>``
          (about 10 minutes, one worker)         -> ScatteringData_PhaseD/held2/clean
  stack   ``scripts/autoresearch/ar3/held_stack_note.py --build`` -> stacks/held2_clean.npy and
          ar3/held2/gradient/held2_stack.json (the stack's lit pixels equal the render record's)

It refuses until ``ar5/stallhit/frozen.json`` is committed and agrees (``stallhit_freeze``). A stage
whose output exists is refused unless ``--resume`` (then it is skipped, and the stages after it
run). Each stage writes ``<cache>/logs/held2_<stage>.log`` and ``.rc``; the end prints the sample's
seed, voxel and grain counts and the stack check and nothing else: no grain, orientation or map of
``held2`` is analysed by this script. ``--dry-run`` prints the commands only.
"""

import argparse
import json
import os
import subprocess
import sys
from pathlib import Path
from typing import List

HERE = Path(__file__).resolve().parent
ICE = HERE.parents[2]
sys.path.insert(0, str(HERE))
import stallhit_freeze as SF  # noqa: E402

CACHE = Path(os.environ.get("STALLHIT_CACHE", str(SF.FROZEN.parent / "cache")))
EX = ICE.parent / "Examples" / "Example2.ManyGrains"
H2 = ICE / "benchmarks" / "autoresearch" / "ar3" / "held2" / "gradient"
TRUTH = H2 / "sample" / "gradient_truth.mic"
FRAMES = EX / "ScatteringData_PhaseD" / "held2" / "clean"
STACK = EX / "ScatteringData_PhaseD" / "stacks" / "held2_clean.npy"


def stages() -> List[dict]:
    return [
        dict(
            name="build",
            out=TRUTH,
            cmd=["uv", "run", "python", "scripts/autoresearch/ar3/held_make_sample.py"],
        ),
        dict(
            name="render",
            out=FRAMES,
            cmd=[
                "uv",
                "run",
                "python",
                "scripts/phase_d/render_full.py",
                "--tag",
                "held2",
                "--workers",
                "1",
                "--variants",
                "clean",
                "--mic",
                str(TRUTH),
            ],
        ),
        dict(
            name="stack",
            out=STACK,
            cmd=[
                "uv",
                "run",
                "python",
                "scripts/autoresearch/ar3/held_stack_note.py",
                "--build",
            ],
        ),
    ]


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--resume", action="store_true")
    a = ap.parse_args()
    if not a.dry_run:
        SF.check_or_exit()
    logs = CACHE / "logs"
    env = dict(os.environ, AR3_SAMPLE="held2")
    env.pop("VIRTUAL_ENV", None)
    for st in stages():
        out: Path = st["out"]
        if out.exists():
            if not a.resume:
                print(f"refusing: {out} already exists (use --resume to skip finished stages)")
                return 2
            print(f"skip {st['name']}: {out} exists")
            continue
        print(f"{st['name']}: {' '.join(st['cmd'])}", flush=True)
        if a.dry_run:
            continue
        logs.mkdir(parents=True, exist_ok=True)
        with open(logs / f"held2_{st['name']}.log", "w") as f:
            rc = subprocess.run(
                st["cmd"], cwd=str(ICE), env=env, stdout=f, stderr=subprocess.STDOUT
            ).returncode
        (logs / f"held2_{st['name']}.rc").write_text(f"{rc}\n")
        if rc != 0:
            print(f"stage {st['name']} failed with exit {rc}")
            return rc
    if not a.dry_run:
        info = json.loads((H2 / "sample" / "sample.json").read_text())
        stk = json.loads((H2 / "held2_stack.json").read_text())
        print(
            "held2 built:",
            dict(
                seed=info["seed"],
                n_voxels=info["n_voxels"],
                n_grains=info["n_grains"],
                stack_matches_render=stk["stack_matches_render"],
            ),
        )
        if info["seed"] != 2 or not stk["stack_matches_render"]:
            return 3
    return 0


if __name__ == "__main__":
    sys.exit(main())
