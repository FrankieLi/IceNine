"""AR5 step 21: the read guards on ``held2`` and ``held5`` (pure; no study imports, nothing is
opened here). Step 20's ``cmarepair_guard`` with ``held5`` in the place of ``held4``.

Step 21 does not use ``held2`` at all. ``held4`` is a disclosed sample now (step 20 ended). Guards:

1. ``refuse_sample(name)``: ``RESEED_SAMPLE=held2`` (or any name that is not a sample of this step)
   raises ``SystemExit(2)``. ``held5`` (the clean test sample) is refused until the gate file
   ``ar5/screen/gate.json`` exists (the gate is evaluated on the disclosed samples first).
2. ``guarded(path)``: a path with a part that contains ``held2`` (as given or resolved with links
   followed, missing parts kept) raises ``SystemExit(2)``. A path with a part that contains
   ``held5`` raises ``SystemExit(2)`` until the marker ``<cache>/logs/screen_runs.done`` exists
   (the run-all script writes it, holding ``exit=<code>``, after the last ``held5`` job, the
   repairs included). Nobody reads ``held5``'s logs or maps before then.
3. The one exception is a producer job (``SCREEN_PRODUCER=1``, set by the run-all script for the
   ``J5C`` and ``J5S`` jobs only): it may open exactly the inputs it needs, by file name
   (``PRODUCER_OK``: the map before and after the pass, the pass's tries, the sample's grid and
   the image stack). Never the truth, a log or a result.
"""

import os
import sys
from pathlib import Path
from typing import Optional, Union

PathLike = Union[str, Path]

DISCLOSED = ("big", "grad", "mid", "held", "held3", "held4")
SAMPLES = DISCLOSED + ("held5",)
HELD_SAMPLES = ("held", "held3", "held4", "held5")
ICE = Path(__file__).resolve().parents[3]
RD = ICE / "benchmarks" / "autoresearch" / "ar5" / "screen"
MARKER = "logs/screen_runs.done"
PRODUCER_OK = (
    "before_J5.npy",
    "after_J5.npy",
    "tries.npz",
    "gradient_grid.mic",
    "held5_clean.npy",
)


def cache_dir() -> Path:
    return Path(os.environ.get("SCREEN_CACHE", str(RD / "cache")))


def refuse(msg: str) -> None:
    """Say why on stderr and exit with code 2 (``SystemExit(2)``)."""
    print(msg, file=sys.stderr)
    raise SystemExit(2)


def gate_file() -> Path:
    return RD / "gate.json"


def refuse_sample(name: str, gate: Optional[Path] = None) -> None:
    """``SystemExit(2)`` unless ``name`` is a sample of step 21 (``held2`` in particular);
    ``held5`` only once the gate file exists."""
    if "held2" in name:
        refuse("screen rule: held2 is not used by step 21 (RESEED_SAMPLE=held2 refused)")
    if name and name not in SAMPLES:
        refuse(f"screen rule: sample {name!r} is not one of {', '.join(SAMPLES)}")
    if name == "held5" and not (gate or gate_file()).exists():
        refuse("screen rule: held5 is not built or run before the gate (screen/gate.json)")


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


def mentions_held5(p: PathLike, cache: Optional[Path] = None) -> bool:
    """A part of the path (below the cache folder, when it is under it) mentions ``held5``."""
    root = cache or cache_dir()
    q = Path(p)
    for base in (root, root.resolve()):
        try:
            q = q.relative_to(base)
            break
        except ValueError:
            continue
    return any("held5" in x for x in _parts(q))


def held5_open(cache: Optional[Path] = None) -> bool:
    """True when the marker that says every ``held5`` run, pass and repair has ended exists."""
    return ((cache or cache_dir()) / MARKER).exists()


def producer() -> bool:
    return os.environ.get("SCREEN_PRODUCER") == "1"


def guarded(p: PathLike, cache: Optional[Path] = None) -> Path:
    """The path as a ``Path``, or ``SystemExit(2)`` when any part of it mentions ``held2``, or
    ``held5`` before the marker of ended runs exists (the producer jobs' own inputs excepted)."""
    if mentions_held2(p):
        refuse(f"screen guard: nothing of held2 is opened by step 21 ({p})")
    if mentions_held5(p, cache) and not held5_open(cache):
        if not (producer() and Path(p).name in PRODUCER_OK):
            refuse(f"screen guard: held5's logs and maps are not read before its runs end ({p})")
    return Path(p)
