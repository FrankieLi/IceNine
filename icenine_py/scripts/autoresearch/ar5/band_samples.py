"""AR5 step 15: support for the new test sample ``held3`` in the script layer (no older script is
edited, no hash of any older freeze moves).

Step 12 added ``held2`` by editing four AR3 scripts (``ar3_common``, ``gradient_common``,
``reseed_common``, ``held_make_sample``, ``held_stack_note``). Those scripts are in older steps'
freeze lists, and this step may not edit older steps' scripts, so ``held3`` is made here at run
time instead:

``boot()``
    When ``RESEED_SAMPLE`` (a run) or ``AR3_SAMPLE`` (a build stage) is ``held3``: import the AR3
    modules as for ``held`` (their asserts accept it), then point them at ``held3``: ``ar3_common.
    OUT`` (so ``gradient_common`` derives ``GRAD``, ``SAMPLE``, ``TRUTH_MIC``, ``GRID_MIC`` and
    ``GRAINMAP_BIG`` from the ``held3`` folder), ``gradient_common.STACK`` (``held3_clean.npy``)
    and ``reseed_common.SAMPLE`` (``held3``, which names the run folders). Then both environment
    variables say ``held3`` for every module imported later. Every other sample: nothing happens.

``run_patched(script, replacements)``
    Runs an AR3 builder script as ``__main__`` after ``boot()``, with exact text replacements
    (each must occur exactly once or this raises): ``held_make_sample.py``'s generator seed line
    becomes ``SEED = 3`` (the plan: ``held`` built exactly as before, generator seed 3) and
    ``held_stack_note.py``'s sample name line becomes ``K = "held3"``. Nothing else changes.

CLI (the build stages of ``band_make_held3.py`` use it)::

    AR3_SAMPLE=held3 uv run python scripts/autoresearch/ar5/band_samples.py make|stack [args]
"""

import os
import sys
from pathlib import Path
from typing import Dict, List

HERE = Path(__file__).resolve().parent
AR3 = HERE.parent / "ar3"
ICE = HERE.parents[2]

SEED_LINE = 'SEED = {"held": 1, "held2": 2}[C.SAMPLE_KIND]'
NAME_LINE = "K = C.SAMPLE_KIND  # AR5 step 12: the same note for held2"
PATCHES: Dict[str, List[List[str]]] = {
    "make": [[SEED_LINE, "SEED = 3  # AR5 step 15: held3 is held with generator seed 3"]],
    "stack": [[NAME_LINE, 'K = "held3"  # AR5 step 15: the same note for held3']],
}
SCRIPTS = {"make": "held_make_sample.py", "stack": "held_stack_note.py"}


def wanted() -> bool:
    return "held3" in (os.environ.get("RESEED_SAMPLE"), os.environ.get("AR3_SAMPLE"))


def boot() -> None:
    """Make the AR3 modules see ``held3`` (see the module text); a no-op for other samples."""
    if not wanted():
        return
    sys.path.insert(0, str(AR3))
    os.environ["RESEED_SAMPLE"] = "held"
    os.environ["AR3_SAMPLE"] = "held"
    import ar3_common as C

    C.OUT = C.OUT.parent / "held3"
    import gradient_common as G

    G.STACK = "held3_clean.npy"
    import reseed_common as RC

    RC.SAMPLE = "held3"
    os.environ["RESEED_SAMPLE"] = "held3"
    os.environ["AR3_SAMPLE"] = "held3"


def patched_source(which: str) -> str:
    """The AR3 builder's source with the exact replacements for ``make`` or ``stack``."""
    src = (AR3 / SCRIPTS[which]).read_text()
    for old, new in PATCHES[which]:
        if src.count(old) != 1:
            raise SystemExit(f"band_samples: {old!r} is not exactly once in {SCRIPTS[which]}")
        src = src.replace(old, new)
    return src


def run_patched(which: str, argv: List[str]) -> None:
    boot()
    path = AR3 / SCRIPTS[which]
    code = compile(patched_source(which), str(path), "exec")
    sys.argv = [str(path)] + argv
    sys.path.insert(0, str(AR3))
    exec(code, {"__name__": "__main__", "__file__": str(path)})  # noqa: S102


if __name__ == "__main__":
    if len(sys.argv) < 2 or sys.argv[1] not in SCRIPTS:
        raise SystemExit("usage: band_samples.py make|stack [args]")
    run_patched(sys.argv[1], sys.argv[2:])
