"""AR5 step 15: the read guards on ``held2`` (pure; no study imports, nothing is opened here).

Step 15 does not use ``held2`` at all (step 12's review opened its seed logs, and step 13's pass
result is not on record). Two guards:

1. ``refuse_sample(name)``: ``RESEED_SAMPLE=held2`` (or any name that is not one of the samples
   of this step: ``big``, ``grad``, ``mid``, ``held``, ``held3``) raises ``SystemExit(2)``.
2. ``guarded(path)``: a path with a part that contains ``held2`` (as given or resolved with links
   followed, missing parts kept: ``J_held2``, ``passes_held2``, ``pass_J_held2.log``,
   ``stats_held2.json``, ``held2_clean.npy``, ...) raises ``SystemExit(2)``. Every loader of this
   step goes through it.
"""

import sys
from pathlib import Path
from typing import Union

PathLike = Union[str, Path]

SAMPLES = ("big", "grad", "mid", "held", "held3")
HELD_SAMPLES = ("held", "held3")  # the ones behind this step's own freeze check


def refuse(msg: str) -> None:
    """Say why on stderr and exit with code 2 (``SystemExit(2)``)."""
    print(msg, file=sys.stderr)
    raise SystemExit(2)


def refuse_sample(name: str) -> None:
    """``SystemExit(2)`` unless ``name`` is a sample of step 15 (``held2`` in particular)."""
    if "held2" in name:
        refuse("band rule: held2 is not used by step 15 (RESEED_SAMPLE=held2 refused)")
    if name and name not in SAMPLES:
        refuse(f"band rule: sample {name!r} is not one of {', '.join(SAMPLES)}")


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


def guarded(p: PathLike) -> Path:
    """The path as a ``Path``, or ``SystemExit(2)`` when any part of it mentions ``held2``."""
    if mentions_held2(p):
        refuse(f"band guard: nothing of held2 is opened by step 15 ({p})")
    return Path(p)
