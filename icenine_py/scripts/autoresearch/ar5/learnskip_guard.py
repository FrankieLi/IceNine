"""AR5 step 18: the read guards on ``held2`` and ``held4`` (pure; no study imports, nothing is
opened here).

Step 18 does not use ``held2`` at all. Guards:

1. ``refuse_sample(name)``: ``RESEED_SAMPLE=held2`` (or any name that is not a sample of this step)
   raises ``SystemExit(2)``. ``held4`` (the clean test sample) is refused until the gate file
   ``ar5/learnskip/gate.json`` exists (the gate is evaluated on the disclosed samples first).
2. ``guarded(path)``: a path with a part that contains ``held2`` (as given or resolved with links
   followed, missing parts kept: ``J_held2``, ``passes_held2``, ``pass_J_held2.log``, ``stats_held2
   .json``, ``held2_clean.npy``, ...) raises ``SystemExit(2)``. A path with a part that contains
   ``held4`` (below the cache when it is under it) raises ``SystemExit(2)`` until the marker file
   ``<cache>/held4_runs_done.json`` exists (the plan: nobody reads ``held4``'s logs or maps before
   its runs end). The run-all script writes the marker after the last ``held4`` run and pass. The
   run scripts themselves write ``held4`` files through older code; this guard is for the readers
   of this step's layer (summary, stats, gate, final).
"""

import os
import sys
from pathlib import Path
from typing import Optional, Union

PathLike = Union[str, Path]

DISCLOSED = ("big", "grad", "mid", "held", "held3")
SAMPLES = DISCLOSED + ("held4",)
HELD_SAMPLES = ("held", "held3", "held4")  # the ones behind this step's own freeze check
ICE = Path(__file__).resolve().parents[3]
FS_DIR = ICE / "benchmarks" / "autoresearch" / "ar5" / "learnskip"
H4_MARKER = "held4_runs_done.json"


def cache_dir() -> Path:
    return Path(os.environ.get("LEARNSKIP_CACHE", str(FS_DIR / "cache")))


def refuse(msg: str) -> None:
    """Say why on stderr and exit with code 2 (``SystemExit(2)``)."""
    print(msg, file=sys.stderr)
    raise SystemExit(2)


def gate_file() -> Path:
    return FS_DIR / "gate.json"


def refuse_sample(name: str, gate: Optional[Path] = None) -> None:
    """``SystemExit(2)`` unless ``name`` is a sample of step 18 (``held2`` in particular);
    ``held4`` only once the gate file exists."""
    if "held2" in name:
        refuse("learnskip rule: held2 is not used by step 18 (RESEED_SAMPLE=held2 refused)")
    if name and name not in SAMPLES:
        refuse(f"learnskip rule: sample {name!r} is not one of {', '.join(SAMPLES)}")
    if name == "held4" and not (gate or gate_file()).exists():
        refuse(
            "learnskip rule: held4 is not built or run before the gate (ar5/learnskip/gate.json)"
        )


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


def mentions_held4(p: PathLike, cache: Optional[Path] = None) -> bool:
    """A part of the path (below the cache folder, when it is under it) mentions ``held4``."""
    root = cache or cache_dir()
    q = Path(p)
    for base in (root, root.resolve()):
        try:
            q = q.relative_to(base)
            break
        except ValueError:
            continue
    return any("held4" in x for x in _parts(q))


def held4_open(cache: Optional[Path] = None) -> bool:
    """True when the marker that says every ``held4`` run and pass has ended exists."""
    return ((cache or cache_dir()) / H4_MARKER).exists()


def guarded(p: PathLike, cache: Optional[Path] = None) -> Path:
    """The path as a ``Path``, or ``SystemExit(2)`` when any part of it mentions ``held2``, or
    ``held4`` before the marker of ended runs exists."""
    if mentions_held2(p):
        refuse(f"learnskip guard: nothing of held2 is opened by step 18 ({p})")
    if mentions_held4(p, cache) and not held4_open(cache):
        refuse(f"learnskip guard: held4's logs and maps are not read before its runs end ({p})")
    return Path(p)
