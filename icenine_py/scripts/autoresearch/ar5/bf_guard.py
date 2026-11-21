"""AR5 step 24: the read guards on ``held2``, ``fine7`` and ``held7`` (pure; no study imports,
nothing is opened here). Step 23's ``ball_guard`` with this step's names, markers and cache.

Step 24 does not use ``held2`` or ``fine7`` (``benchmarks/autoresearch/fine1``) at all. ``held7`` is
the clean test sample (built by this step); ``held6`` is disclosed. Guards:

1. ``refuse_sample(name)``: ``RESEED_SAMPLE`` with ``held2`` or ``fine7`` in it, or a name that is
   not a sample of this step, raises ``SystemExit(2)``. ``held7`` (the clean test sample) is
   refused until the file ``<cache>/logs/gate_passed`` exists (written by the implementer after
   g0, g1 and g2 are judged).
2. ``guarded(path)``: a path with a part that contains ``held2`` or ``fine7``, or a part equal to
   ``fine1`` (as given or resolved with links followed, missing parts kept) raises
   ``SystemExit(2)``. A path with a part that contains ``held7`` raises ``SystemExit(2)`` until
   the marker ``<cache>/logs/bf_runs.done`` exists (the run-all script writes it, holding
   ``exit=<code>``, after the last ``held7`` job, the passes and the repairs included). Nobody
   reads ``held7``'s logs or maps before then.
3. The one exception is a producer job (``SCREEN_PRODUCER=1``, set by the run-all script for the
   pass and screen jobs): it may open exactly the inputs it needs, by file name (``PRODUCER_OK``).
   Never the truth, a log or a result.

``defer_common`` (step 22's set-up) and ``ball_common`` (step 23's) call the guard functions of
their own steps; ``bf_common`` puts these functions in the place of theirs in those modules, in
this process only; no older file is edited.
"""

import os
import sys
from pathlib import Path
from typing import Any, Optional, Union

PathLike = Union[str, Path]

DISCLOSED = ("big", "grad", "mid", "held", "held3", "held4", "held5", "held6")
SAMPLES = DISCLOSED + ("held7",)
HELD_SAMPLES = ("held", "held3", "held4", "held5", "held6", "held7")
ICE = Path(__file__).resolve().parents[3]
RD = ICE / "benchmarks" / "autoresearch" / "ar5" / "bf"
MARKER = "logs/bf_runs.done"
GATE_PASSED = "logs/gate_passed"
PRODUCER_OK = (
    "before_J5DB.npy",
    "after_J5DB.npy",
    "before_J5.npy",
    "after_J5.npy",
    "tries.npz",
    "gradient_grid.mic",
    "held7_clean.npy",
)


def base_cache() -> Path:
    """The main cache of this step (the J5SDC job's own sub-cache sits inside it)."""
    return Path(
        os.environ.get(
            "BF_CACHE",
            os.environ.get("DEFER_CACHE", os.environ.get("SCREEN_CACHE", str(RD / "cache"))),
        )
    )


def refuse(msg: str) -> None:
    """Say why on stderr and exit with code 2 (``SystemExit(2)``)."""
    print(msg, file=sys.stderr)
    raise SystemExit(2)


def refuse_sample(name: str, cache: Optional[Path] = None) -> None:
    """``SystemExit(2)`` unless ``name`` is a sample of step 24 (``held2`` and ``fine7`` in
    particular); ``held7`` only once ``<cache>/logs/gate_passed`` exists."""
    if "held2" in name:
        refuse("bf rule: held2 is not used by step 24 (RESEED_SAMPLE=held2 refused)")
    if "fine7" in name:
        refuse("bf rule: fine7 is not read by step 24 (refused)")
    if name and name not in SAMPLES:
        refuse(f"bf rule: sample {name!r} is not one of {', '.join(SAMPLES)}")
    if name == "held7" and not ((cache or base_cache()) / GATE_PASSED).exists():
        refuse("bf rule: held7 is not built or run before the gate is judged (gate_passed)")


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


def mentions_fine7(p: PathLike) -> bool:
    return any(("fine7" in x or x == "fine1") for x in _parts(p))


def mentions_held7(p: PathLike, cache: Optional[Path] = None) -> bool:
    """A part of the path (below the cache folder, when it is under it) mentions ``held7``."""
    root = cache or base_cache()
    q = Path(p)
    for base in (root, root.resolve()):
        try:
            q = q.relative_to(base)
            break
        except ValueError:
            continue
    return any("held7" in x for x in _parts(q))


def held7_open(cache: Optional[Path] = None) -> bool:
    """True when the marker that says every ``held7`` run, pass and repair has ended exists."""
    return ((cache or base_cache()) / MARKER).exists()


def producer() -> bool:
    return os.environ.get("SCREEN_PRODUCER") == "1"


def guarded(p: PathLike, cache: Optional[Path] = None) -> Path:
    """The path as a ``Path``, or ``SystemExit(2)`` when any part of it mentions ``held2`` or
    ``fine7``, or ``held7`` before the marker of ended runs exists (the producer jobs' own inputs
    excepted)."""
    if mentions_held2(p):
        refuse(f"bf guard: nothing of held2 is opened by step 24 ({p})")
    if mentions_fine7(p):
        refuse(f"bf guard: nothing of fine7 is read by step 24 ({p})")
    if mentions_held7(p, cache) and not held7_open(cache):
        if not (producer() and Path(p).name in PRODUCER_OK):
            refuse(f"bf guard: held7's logs and maps are not read before its runs end ({p})")
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
