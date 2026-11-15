"""AR5 step 14: the two read guards (pure; no study imports except ``polish_freeze``, lazily).

1. Nothing under a ``passes_held2`` folder is opened, listed or statted until step 13's ``held2``
   result is on record (``step13_recorded=True``; phase A never sets it).
2. Nothing of ``held2``'s step-12 run folders (``H_held2``, ``G_held2``, ``A_held2``, ``S8_held2``,
   any ``*_held2`` folder of a cache) is opened before ``ar5/polish/frozen.json`` is committed,
   clean and agrees (``polish_freeze.problems() == ''``); this holds for the rescue logs in
   particular (``rescue_log.json.gz``, ``rescue.json``).

``guarded(path)`` returns the path as a ``Path`` or raises ``SystemExit`` (code 2); every loader of
step 14 goes through it. It looks at the path as given and as resolved (links followed, missing
parts kept), so a link into a ``passes_held2`` folder is refused too. No file is opened here.
"""

from pathlib import Path
from typing import Callable, Optional, Union

PathLike = Union[str, Path]

PASSES = "passes_held2"


def _parts(p: PathLike) -> list:
    q = Path(p)
    out = list(q.parts)
    try:
        out += list(q.resolve().parts)
    except OSError:
        pass
    return out


def in_passes_held2(p: PathLike) -> bool:
    return PASSES in _parts(p)


def in_held2_run(p: PathLike) -> bool:
    """A step-12 run folder of ``held2`` (a path part ``<arm>_held2``, not ``passes_held2``)."""
    return any(x.endswith("_held2") and x != PASSES for x in _parts(p))


def freeze_problems() -> str:
    import polish_freeze as PF

    return PF.problems()


def guarded(
    p: PathLike,
    step13_recorded: bool = False,
    problems: Optional[Callable[[], str]] = None,
) -> Path:
    """The path, or ``SystemExit(2)`` when it is under ``passes_held2`` (step 13's held2 result
    not recorded) or in a ``*_held2`` run folder before the freeze commit."""
    q = Path(p)
    if in_passes_held2(q) and not step13_recorded:
        raise SystemExit(f"polish guard: {PASSES} is not opened before step 13's held2 result")
    if in_held2_run(q):
        why = (problems or freeze_problems)()
        if why:
            raise SystemExit(f"polish guard: held2 is not read before the freeze commit ({why})")
    return q
