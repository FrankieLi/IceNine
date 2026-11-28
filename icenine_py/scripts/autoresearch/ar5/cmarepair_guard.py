"""AR5 step 20: the read guards on ``held2`` and ``held4`` (pure; no study imports, nothing is
opened here).

Step 20 does not use ``held2`` at all. Guards:

1. ``refuse_sample(name)``: ``RESEED_SAMPLE=held2`` (or any name that is not a sample of this step)
   raises ``SystemExit(2)``. ``held4`` (the clean test sample) is refused until the gate file
   ``ar5/cmarepair/gate.json`` exists (the gate is evaluated on the disclosed samples first).
2. ``guarded(path)``: a path with a part that contains ``held2`` (as given or resolved with links
   followed, missing parts kept) raises ``SystemExit(2)``. A path with a part that contains
   ``held4`` raises ``SystemExit(2)`` until the marker ``<cache>/logs/cmarepair_runs.done`` exists
   (the run-all script writes it, holding ``exit=<code>``, after the last ``held4`` job, repair
   included). Nobody reads ``held4``'s logs or maps before then.
3. The one exception is the repair job itself (``CMAREPAIR_PRODUCER=1``, set by the run-all script
   for that job only): it may open exactly the inputs it needs, by file name (``PRODUCER_OK``: the
   map before and after the pass, the pass's tries, the sample's grid and the image stack). Never
   the truth, a log or a result.
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
RD = ICE / "benchmarks" / "autoresearch" / "ar5" / "cmarepair"
MARKER = "logs/cmarepair_runs.done"
PRODUCER_OK = (
    "before_J5.npy",
    "after_J5.npy",
    "tries.npz",
    "gradient_grid.mic",
    "held4_clean.npy",
)


def cache_dir() -> Path:
    return Path(os.environ.get("CMAREPAIR_CACHE", str(RD / "cache")))


def refuse(msg: str) -> None:
    """Say why on stderr and exit with code 2 (``SystemExit(2)``)."""
    print(msg, file=sys.stderr)
    raise SystemExit(2)


def gate_file() -> Path:
    return RD / "gate.json"


def refuse_sample(name: str, gate: Optional[Path] = None) -> None:
    """``SystemExit(2)`` unless ``name`` is a sample of step 20 (``held2`` in particular);
    ``held4`` only once the gate file exists."""
    if "held2" in name:
        refuse("cmarepair rule: held2 is not used by step 20 (RESEED_SAMPLE=held2 refused)")
    if name and name not in SAMPLES:
        refuse(f"cmarepair rule: sample {name!r} is not one of {', '.join(SAMPLES)}")
    if name == "held4" and not (gate or gate_file()).exists():
        refuse("cmarepair rule: held4 is not built or run before the gate (cmarepair/gate.json)")


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
    """True when the marker that says every ``held4`` run, pass and repair has ended exists."""
    return ((cache or cache_dir()) / MARKER).exists()


def producer() -> bool:
    return os.environ.get("CMAREPAIR_PRODUCER") == "1"


def guarded(p: PathLike, cache: Optional[Path] = None) -> Path:
    """The path as a ``Path``, or ``SystemExit(2)`` when any part of it mentions ``held2``, or
    ``held4`` before the marker of ended runs exists (the repair job's own inputs excepted)."""
    if mentions_held2(p):
        refuse(f"cmarepair guard: nothing of held2 is opened by step 20 ({p})")
    if mentions_held4(p, cache) and not held4_open(cache):
        if not (producer() and Path(p).name in PRODUCER_OK):
            refuse(f"cmarepair guard: held4's logs and maps are not read before its runs end ({p})")
    return Path(p)
