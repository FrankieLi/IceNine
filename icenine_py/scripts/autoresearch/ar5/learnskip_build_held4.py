"""AR5 step 18, phase A step 7: build and render the clean test sample ``held4`` (``held`` with
generator seed 4). Step 16's ``high500_make_held3.py`` with this step's freeze, gate, cache and
sample support (``learnskip_samples.py``).

  uv run python scripts/autoresearch/ar5/learnskip_build_held4.py [--dry-run] [--resume]

``held4`` is ``held`` (AR3 step 11: 250 Voronoi grains on the 24,570-voxel gradient grid, gradients
0 to 4 degrees, clean images, Q-max 8) built by the same script with the generator seed changed to
4: a new sample for a carry-over test, not an extra random seed for a comparison. Three stages, each
the command of ``held``'s own build, with ``AR3_SAMPLE=held4``:

  build   ``learnskip_samples.py make``    (``held_make_sample.py``)  -> ar3/held4/gradient/sample/
  render  ``scripts/phase_d/render_full.py --tag held4 --workers 1 --variants clean --mic <truth>``
          (about 10 minutes, one worker)  -> ScatteringData_PhaseD/held4/clean
  stack   ``learnskip_samples.py stack --build`` (``held_stack_note.py``) -> stacks/held4_clean.npy
          and ar3/held4/gradient/held4_stack.json (the stack's lit pixels equal the render record's)

It refuses until ``ar5/learnskip/gate.json`` exists (the gate passed on the disclosed samples) and
``ar5/learnskip/frozen.json`` is committed and agrees. A stage counts as finished only when
``ar3/held4/build_state.json`` (written here after each stage with exit 0) says so and its output
exists. Without ``--resume`` any existing output is refused. With ``--resume`` a finished stage is
skipped, and the output of an unfinished stage is MOVED to
``<cache>/interrupted_<time>/held4_<stage>/`` and the stage is run again. Each stage writes
``<cache>/logs/held4_<stage>.log`` and ``.rc``; the end prints the sample's seed, voxel and grain
counts and the stack check and nothing else: no grain, orientation or map of ``held4`` is analysed
by this script. ``--dry-run`` prints the commands only.
"""

import argparse
import json
import os
import shutil
import subprocess
import sys
import time
from pathlib import Path
from typing import List

HERE = Path(__file__).resolve().parent
ICE = HERE.parents[2]
sys.path.insert(0, str(HERE))
import learnskip_freeze as LF  # noqa: E402
import learnskip_guard as LG  # noqa: E402

CACHE = LG.cache_dir()
EX = ICE.parent / "Examples" / "Example2.ManyGrains"
H4 = ICE / "benchmarks" / "autoresearch" / "ar3" / "held4" / "gradient"
STATE = ICE / "benchmarks" / "autoresearch" / "ar3" / "held4" / "build_state.json"
TRUTH = H4 / "sample" / "gradient_truth.mic"
FRAMES = EX / "ScatteringData_PhaseD" / "held4" / "clean"
STACK = EX / "ScatteringData_PhaseD" / "stacks" / "held4_clean.npy"
SAMPLES_PY = "scripts/autoresearch/ar5/learnskip_samples.py"


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
                "held4",
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


def read_state() -> dict:
    return json.loads(STATE.read_text()) if STATE.exists() else {}


def write_state(st: dict) -> None:
    STATE.parent.mkdir(parents=True, exist_ok=True)
    STATE.write_text(json.dumps(st, indent=1) + "\n")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--resume", action="store_true")
    a = ap.parse_args()
    if not a.dry_run:
        LG.refuse_sample("held4")  # the gate file must exist
        LF.check_or_exit()
    logs = CACHE / "logs"
    env = dict(os.environ, AR3_SAMPLE="held4")
    env.pop("VIRTUAL_ENV", None)
    state = read_state()
    stamp = time.strftime("%Y%m%d-%H%M%S")
    for st in stages():
        out: Path = st["out"]
        finished = state.get(st["name"]) == 0 and out.exists()
        if out.exists():
            if not a.resume:
                print(f"refusing: {out} already exists (use --resume to skip finished stages)")
                return 2
            if finished:
                print(f"skip {st['name']}: finished earlier, {out} exists")
                continue
            dest = CACHE / f"interrupted_{stamp}" / f"held4_{st['name']}"
            print(f"{st['name']}: unfinished output moved to {dest}")
            if not a.dry_run:
                dest.parent.mkdir(parents=True, exist_ok=True)
                shutil.move(str(out), str(dest))
        print(f"{st['name']}: {' '.join(st['cmd'])}", flush=True)
        if a.dry_run:
            continue
        logs.mkdir(parents=True, exist_ok=True)
        with open(logs / f"held4_{st['name']}.log", "w") as f:
            rc = subprocess.run(
                st["cmd"], cwd=str(ICE), env=env, stdout=f, stderr=subprocess.STDOUT
            ).returncode
        (logs / f"held4_{st['name']}.rc").write_text(f"{rc}\n")
        state = read_state()
        state[st["name"]] = rc
        write_state(state)
        if rc != 0:
            print(f"stage {st['name']} failed with exit {rc}")
            return rc
    if not a.dry_run:
        info = json.loads((H4 / "sample" / "sample.json").read_text())
        stk = json.loads((H4 / "held4_stack.json").read_text())
        print(
            "held4 built:",
            dict(
                seed=info["seed"],
                n_voxels=info["n_voxels"],
                n_grains=info["n_grains"],
                stack_matches_render=stk["stack_matches_render"],
            ),
        )
        if info["seed"] != 4 or not stk["stack_matches_render"]:
            return 3
    return 0


if __name__ == "__main__":
    sys.exit(main())
