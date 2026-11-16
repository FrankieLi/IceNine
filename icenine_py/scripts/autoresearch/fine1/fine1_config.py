"""fine1: the run configs of the four arms (pure text; the launchers and the run scripts use it).

Every arm's config is its base config with only these changed: ``SampleFilename`` (the fine grid
or the slice), ``InfileBasename`` (the fine images), ``MinSideLength`` and ``MaxInitSideLength``
0.0046875 (C++ splits any voxel above ``MinSideLength``; Python does not read them) and
``OutStructureBasename``. The C++ config is the Python arm ``PS`` config through
``cpptime/make_config.py`` (``convert``, imported unchanged: it drops the three Python-only
lines and makes the paths absolute). ``J5S`` adds the 12 library keys of AR5 step 21's
``lib_check.json``.

  PS   ``ar2/key_off_fused.config`` (BFSRevisitRefit 1, BFSReseed 1, FusedCost 1), grid search
       for every seed, stock MC local optimization
  CS   the C++ config of ``PS`` (mode p, mpirun -n 2: one server, one fitting client)
  CF   as CS, run with a fixed C++ build (an extra labelled run)
  J5S  ``ar2/key_on_fused.config`` (SeedByIndexing 1) plus the 12 keys below
"""

import json
import re
import sys
from pathlib import Path
from typing import Dict, List, Optional

HERE = Path(__file__).resolve().parent
ICE = HERE.parents[2]  # icenine_py
sys.path.insert(0, str(ICE / "scripts" / "autoresearch" / "cpptime"))
import make_config as MC  # noqa: E402  (cpptime, imported unchanged)

AR2 = ICE / "scripts" / "autoresearch" / "ar2"
BASE_PS = AR2 / "key_off_fused.config"
BASE_J5S = AR2 / "key_on_fused.config"
LIB_CHECK = ICE / "benchmarks" / "autoresearch" / "ar5" / "screen" / "lib_check.json"
SIDE = "0.0046875"
ARMS = ("CS", "J5S", "PS", "CF")  # the order the launchers run them
PYTHON_ARMS = ("J5S", "PS")
CPP_ARMS = ("CS", "CF")

J5S_KEYS: List[str] = [
    "BFSNeighborCMAFirst 1",
    "CMASigma0 0.05",
    "CMANeighborMaxEvals 60",
    "BoundaryPass 1",
    "SeedIndexCandidates 20 5",
    "SeedIndexStopHitRatio 0.95",
    "BFSNeighborCMANoFallback 1",
    "SeedIndexCMAFirst 200",
    "SeedRescueCMA 0.25 1000",
    "SeedIndexMCStallHit 500",
    "BoundaryPassRepair 1",
    "BoundaryPassScreen 1",
]
CHANGED_KEYS = (
    "SampleFilename",
    "InfileBasename",
    "MinSideLength",
    "MaxInitSideLength",
    "OutStructureBasename",
)


def _set(text: str, key: str, value: str) -> str:
    """Replace the value of the one line starting with ``key`` (keeps nothing else of the line)."""
    new, n = re.subn(rf"(?m)^{key}\b.*$", f"{key} {value}", text)
    assert n == 1, (key, n)
    return new


def python_config(arm: str, sample: Path, frames_basename: str, out_basename: str = "fine1") -> str:
    """The Python config text of ``PS`` or ``J5S``."""
    assert arm in PYTHON_ARMS, arm
    text = (BASE_J5S if arm == "J5S" else BASE_PS).read_text()
    text = _set(text, "SampleFilename", str(sample))
    text = _set(text, "InfileBasename", frames_basename)
    text = _set(text, "MinSideLength", SIDE)
    text = _set(text, "MaxInitSideLength", SIDE)
    text = _set(text, "OutStructureBasename", out_basename)
    if arm == "J5S":
        text = text.replace("# END OF FILE", "")
        text = text.rstrip("\n") + "\n\n# fine1 arm J5S: AR5 step 21's library keys\n"
        text += "\n".join(J5S_KEYS) + "\n\n# END OF FILE\n"
    return text


def cpp_config(sample: Path, frames_basename: str, out_basename: str = "fine1_cpp") -> str:
    """The C++ config text (arms CS and CF): the ``PS`` config through ``make_config.convert``."""
    text = MC.convert(BASE_PS.read_text(), sample, out_basename)
    text = _set(text, "InfileBasename", frames_basename)
    text = _set(text, "MinSideLength", SIDE)
    text = _set(text, "MaxInitSideLength", SIDE)
    return text


def keys_of(text: str) -> Dict[str, str]:
    """{key: rest of line} of the non-comment lines of a config text (a repeated key keeps its
    last line)."""
    out: Dict[str, str] = {}
    for ln in text.splitlines():
        s = ln.strip()
        if not s or s.startswith("#"):
            continue
        tok = s.split(None, 1)
        out[tok[0]] = tok[1].split("#")[0].strip() if len(tok) > 1 else ""
    return out


def changed_keys(base: str, new: str) -> List[str]:
    """Sorted keys whose value differs between two config texts, or that only one has."""
    a, b = keys_of(base), keys_of(new)
    return sorted(k for k in set(a) | set(b) if a.get(k) != b.get(k))


def lib_check_keys(path: Optional[Path] = None) -> List[str]:
    """The key list of step 21's ``lib_check.json`` (the library run of arm J5S)."""
    return list(json.loads((path or LIB_CHECK).read_text())["keys"])


def main() -> None:
    """``fine1_config.py --arm CS --grid slice --out FILE``: write one arm's config text."""
    import argparse

    import fine1_paths as P

    ap = argparse.ArgumentParser(description=main.__doc__)
    ap.add_argument("--arm", required=True, choices=ARMS)
    ap.add_argument("--grid", required=True, choices=("slice", "full"))
    ap.add_argument("--out", type=Path, required=True)
    a = ap.parse_args()
    grid = P.SLICE if a.grid == "slice" else P.GRID
    frames = str(P.FRAMES / "500Grains.sim")
    if a.arm in CPP_ARMS:
        text = cpp_config(grid, frames, "fine1_cpp")
    else:
        text = python_config(a.arm, grid, frames, f"fine1_{a.arm}")
    a.out.parent.mkdir(parents=True, exist_ok=True)
    a.out.write_text(text)
    print(f"wrote {a.out}")


if __name__ == "__main__":
    main()
