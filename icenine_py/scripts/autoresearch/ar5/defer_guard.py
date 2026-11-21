"""AR5 step 22: the read guards on ``held2`` and ``held6`` (pure; no study imports, nothing is
opened here). Step 21's ``screen_guard`` with ``held6`` in the place of ``held5``.

Step 22 does not use ``held2`` at all. ``held5`` is a disclosed sample now. Guards:

1. ``refuse_sample(name)``: ``RESEED_SAMPLE=held2`` (or any name that is not a sample of this step)
   raises ``SystemExit(2)``. ``held6`` (the clean test sample) is refused until the file
   ``<cache>/logs/gate_passed`` exists (written by the implementer after the gate is judged).
2. ``guarded(path)``: a path with a part that contains ``held2`` (as given or resolved with links
   followed, missing parts kept) raises ``SystemExit(2)``. A path with a part that contains
   ``held6`` raises ``SystemExit(2)`` until the marker ``<cache>/logs/defer_runs.done`` exists (the
   run-all script writes it, holding ``exit=<code>``, after the last ``held6`` job, the repairs
   included). Nobody reads ``held6``'s logs or maps before then.
3. The one exception is a producer job (``SCREEN_PRODUCER=1``, set by the run-all script for the
   pass and screen jobs): it may open exactly the inputs it needs, by file name (``PRODUCER_OK``:
   the map before and after the pass, the pass's tries, the sample's grid and the image stack).
   Never the truth, a log or a result.

``install()`` puts these functions in the place of step 21's in ``screen_guard`` (same names), in
this process only; step 21's file is not edited. The main cache is ``DEFER_CACHE`` (the J5SD job's
own sub-cache sits inside it and has its own ``SCREEN_CACHE``).
"""

import os
import sys
from pathlib import Path
from typing import Any, Optional, Union

PathLike = Union[str, Path]

DISCLOSED = ("big", "grad", "mid", "held", "held3", "held4", "held5")
SAMPLES = DISCLOSED + ("held6",)
HELD_SAMPLES = ("held", "held3", "held4", "held5", "held6")
ICE = Path(__file__).resolve().parents[3]
RD = ICE / "benchmarks" / "autoresearch" / "ar5" / "defer"
MARKER = "logs/defer_runs.done"
GATE_PASSED = "logs/gate_passed"
PRODUCER_OK = (
    "before_J5.npy",
    "after_J5.npy",
    "tries.npz",
    "gradient_grid.mic",
    "held6_clean.npy",
)


def base_cache() -> Path:
    """The main cache of this step (the J5SD job's own sub-cache sits inside it)."""
    return Path(os.environ.get("DEFER_CACHE", os.environ.get("SCREEN_CACHE", str(RD / "cache"))))


def refuse(msg: str) -> None:
    """Say why on stderr and exit with code 2 (``SystemExit(2)``)."""
    print(msg, file=sys.stderr)
    raise SystemExit(2)


def refuse_sample(name: str, cache: Optional[Path] = None) -> None:
    """``SystemExit(2)`` unless ``name`` is a sample of step 22 (``held2`` in particular);
    ``held6`` only once ``<cache>/logs/gate_passed`` exists."""
    if "held2" in name:
        refuse("defer rule: held2 is not used by step 22 (RESEED_SAMPLE=held2 refused)")
    if name and name not in SAMPLES:
        refuse(f"defer rule: sample {name!r} is not one of {', '.join(SAMPLES)}")
    if name == "held6" and not ((cache or base_cache()) / GATE_PASSED).exists():
        refuse("defer rule: held6 is not built or run before the gate is judged (gate_passed)")


def _parts(p: PathLike) -> list:
    q = Path(p)
    out = list(q.parts)
    try:
        out += list(q.resolve().parts)
    except OSError:
        pass
    return out


def mentions_held2(p: PathLike) -> bool:
    return any("held2" in x for x in _parts(p))


def mentions_held6(p: PathLike, cache: Optional[Path] = None) -> bool:
    """A part of the path (below the cache folder, when it is under it) mentions ``held6``."""
    root = cache or base_cache()
    q = Path(p)
    for base in (root, root.resolve()):
        try:
            q = q.relative_to(base)
            break
        except ValueError:
            continue
    return any("held6" in x for x in _parts(q))


def held6_open(cache: Optional[Path] = None) -> bool:
    """True when the marker that says every ``held6`` run, pass and repair has ended exists."""
    return ((cache or base_cache()) / MARKER).exists()


def producer() -> bool:
    return os.environ.get("SCREEN_PRODUCER") == "1"


def guarded(p: PathLike, cache: Optional[Path] = None) -> Path:
    """The path as a ``Path``, or ``SystemExit(2)`` when any part of it mentions ``held2``, or
    ``held6`` before the marker of ended runs exists (the producer jobs' own inputs excepted)."""
    if mentions_held2(p):
        refuse(f"defer guard: nothing of held2 is opened by step 22 ({p})")
    if mentions_held6(p, cache) and not held6_open(cache):
        if not (producer() and Path(p).name in PRODUCER_OK):
            refuse(f"defer guard: held6's logs and maps are not read before its runs end ({p})")
    return Path(p)


def install(sg: Any) -> None:
    """Replace step 21's guard functions in the module ``sg`` (``screen_guard``) by these."""
    sg.SAMPLES = SAMPLES
    sg.HELD_SAMPLES = HELD_SAMPLES
    sg.DISCLOSED = DISCLOSED
    sg.refuse_sample = refuse_sample
    sg.guarded = guarded
    sg.mentions_held5 = lambda *a, **k: False  # held5 is disclosed in this step
    sg.held5_open = lambda *a, **k: True
