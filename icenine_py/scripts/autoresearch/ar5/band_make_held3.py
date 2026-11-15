"""AR5 step 15, phase A step 6: build and render the new clean test sample ``held3`` (``held`` with
generator seed 3).

  uv run python scripts/autoresearch/ar5/band_make_held3.py [--dry-run] [--resume]

``held3`` is ``held`` (AR3 step 11: 250 Voronoi grains on the 24,570-voxel gradient grid, gradients
0 to 4 degrees, clean images, Q-max 8) built by the same script with the generator seed changed to
3: a new sample for a carry-over test, not an extra random seed for a comparison. Three stages, each
the command of ``held``'s own build, with ``AR3_SAMPLE=held3`` (``band_samples.py`` makes the AR3
modules accept it at run time and replaces exactly one line of ``held_make_sample.py`` and one of
``held_stack_note.py``; no older script is edited):

  build   ``band_samples.py make``    (``held_make_sample.py``)     -> ar3/held3/gradient/sample/
  render  ``scripts/phase_d/render_full.py --tag held3 --workers 1 --variants clean --mic <truth>``
          (about 10 minutes, one worker)         -> ScatteringData_PhaseD/held3/clean
  stack   ``band_samples.py stack --build`` (``held_stack_note.py``) -> stacks/held3_clean.npy and
          ar3/held3/gradient/held3_stack.json (the stack's lit pixels equal the render record's)

It refuses until ``ar5/band/frozen.json`` is committed and agrees (``band_freeze``). A stage whose
output exists is refused unless ``--resume`` (then it is skipped, and the stages after it run). Each
stage writes ``<cache>/logs/held3_<stage>.log`` and ``.rc``; the end prints the sample's seed, voxel
and grain counts and the stack check and nothing else: no grain, orientation or map of ``held3`` is
analysed by this script. ``--dry-run`` prints the commands only.
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
import band_freeze as BF  # noqa: E402

CACHE = Path(os.environ.get("BAND_CACHE", str(BF.BD / "cache")))
EX = ICE.parent / "Examples" / "Example2.ManyGrains"
H3 = ICE / "benchmarks" / "autoresearch" / "ar3" / "held3" / "gradient"
TRUTH = H3 / "sample" / "gradient_truth.mic"
FRAMES = EX / "ScatteringData_PhaseD" / "held3" / "clean"
STACK = EX / "ScatteringData_PhaseD" / "stacks" / "held3_clean.npy"
SAMPLES_PY = "scripts/autoresearch/ar5/band_samples.py"


def stages() -> List[dict]:
    return [
        dict(name="build", out=TRUTH, cmd=["uv", "run", "python", SAMPLES_PY, "make"]),
        dict(
            name="render",
            out=FRAMES,
            cmd=[
                "uv",
                "run",
                "python",
                "scripts/phase_d/render_full.py",
                "--tag",
                "held3",
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
            cmd=["uv", "run", "python", SAMPLES_PY, "stack", "--build"],
        ),
    ]


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--resume", action="store_true")
    a = ap.parse_args()
    if not a.dry_run:
        BF.check_or_exit()
    logs = CACHE / "logs"
    env = dict(os.environ, AR3_SAMPLE="held3")
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
        with open(logs / f"held3_{st['name']}.log", "w") as f:
            rc = subprocess.run(
                st["cmd"], cwd=str(ICE), env=env, stdout=f, stderr=subprocess.STDOUT
            ).returncode
        (logs / f"held3_{st['name']}.rc").write_text(f"{rc}\n")
        if rc != 0:
            print(f"stage {st['name']} failed with exit {rc}")
            return rc
    if not a.dry_run:
        info = json.loads((H3 / "sample" / "sample.json").read_text())
        stk = json.loads((H3 / "held3_stack.json").read_text())
        print(
            "held3 built:",
            dict(
                seed=info["seed"],
                n_voxels=info["n_voxels"],
                n_grains=info["n_grains"],
                stack_matches_render=stk["stack_matches_render"],
            ),
        )
        if info["seed"] != 3 or not stk["stack_matches_render"]:
            return 3
    return 0


if __name__ == "__main__":
    sys.exit(main())
