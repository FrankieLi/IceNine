"""AR5 step 18: support for the new clean test sample ``held4`` in the script layer (no older script
is edited, no hash of any older freeze moves). Step 15's ``band_samples.py`` does this for ``held3``
only; this file does it for ``held3`` and ``held4`` with the same text replacements.

``held4`` is ``held`` (AR3 step 11: 250 Voronoi grains on the 24,570-voxel gradient grid, gradients
0 to 4 degrees, clean images, Q-max 8) built by the same script with the generator seed changed to
4: a new sample for a carry-over test, not an extra random seed for a comparison.

``boot()``: when ``RESEED_SAMPLE`` or ``AR3_SAMPLE`` is ``held3`` or ``held4``: import the AR3
modules as for ``held`` (their asserts accept it), then point them at the sample (``ar3_common.OUT``
so ``gradient_common`` derives its paths from the sample's folder, ``gradient_common.STACK``,
``reseed_common.SAMPLE``); both environment variables then say the sample. ``learnskip_common``
installs it in place of ``band_samples.boot`` (``band_samples.wanted`` too) before step 16's set-up
is imported.

CLI (the build stages of ``learnskip_build_held4.py`` use it)::

    AR3_SAMPLE=held4 uv run python scripts/autoresearch/ar5/learnskip_samples.py make|stack [args]
"""

import os
import sys
from pathlib import Path
from typing import Dict, List, Optional

HERE = Path(__file__).resolve().parent
AR3 = HERE.parent / "ar3"
ICE = HERE.parents[2]

SEEDS: Dict[str, int] = {"held3": 3, "held4": 4}
SEED_LINE = 'SEED = {"held": 1, "held2": 2}[C.SAMPLE_KIND]'
NAME_LINE = "K = C.SAMPLE_KIND  # AR5 step 12: the same note for held2"
SCRIPTS = {"make": "held_make_sample.py", "stack": "held_stack_note.py"}


def which_sample() -> Optional[str]:
    """``held3`` or ``held4`` when ``RESEED_SAMPLE`` or ``AR3_SAMPLE`` says so, else None."""
    for v in (os.environ.get("RESEED_SAMPLE"), os.environ.get("AR3_SAMPLE")):
        if v in SEEDS:
            return v
    return None


def wanted() -> bool:
    return which_sample() is not None


def boot() -> None:
    """Make the AR3 modules see ``held3`` or ``held4``; a no-op for other samples."""
    name = which_sample()
    if name is None:
        return
    sys.path.insert(0, str(AR3))
    os.environ["RESEED_SAMPLE"] = "held"
    os.environ["AR3_SAMPLE"] = "held"
    import ar3_common as C

    C.OUT = C.OUT.parent / name
    import gradient_common as G

    G.STACK = f"{name}_clean.npy"
    import reseed_common as RC

    RC.SAMPLE = name
    os.environ["RESEED_SAMPLE"] = name
    os.environ["AR3_SAMPLE"] = name


def patched_source(which: str, name: str) -> str:
    """The AR3 builder's source with the exact replacements for ``make`` or ``stack``."""
    patches: Dict[str, List[List[str]]] = {
        "make": [
            [SEED_LINE, f"SEED = {SEEDS[name]}  # AR5 step 18: {name} is held, seed {SEEDS[name]}"]
        ],
        "stack": [[NAME_LINE, f'K = "{name}"  # AR5 step 18: the same note for {name}']],
    }
    src = (AR3 / SCRIPTS[which]).read_text()
    for old, new in patches[which]:
        if src.count(old) != 1:
            raise SystemExit(f"learnskip_samples: {old!r} is not exactly once in {SCRIPTS[which]}")
        src = src.replace(old, new)
    return src


def run_patched(which: str, argv: List[str]) -> None:
    name = which_sample()
    if name is None:
        raise SystemExit("learnskip_samples: set AR3_SAMPLE to held3 or held4")
    boot()
    path = AR3 / SCRIPTS[which]
    code = compile(patched_source(which, name), str(path), "exec")
    sys.argv = [str(path)] + argv
    sys.path.insert(0, str(AR3))
    exec(code, {"__name__": "__main__", "__file__": str(path)})  # noqa: S102


if __name__ == "__main__":
    if len(sys.argv) < 2 or sys.argv[1] not in SCRIPTS:
        raise SystemExit("usage: learnskip_samples.py make|stack [args]")
    run_patched(sys.argv[1], sys.argv[2:])
