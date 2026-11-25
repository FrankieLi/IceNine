"""tex1: the run configs of the two arms (pure text; the launcher and the run script use it).

Every arm's config is its base config with only these changed: ``SampleFilename`` (the fine grid
or the slice), ``InfileBasename`` (the tex24 images), ``MinSideLength`` and ``MaxInitSideLength``
0.0046875 (C++ splits any voxel above ``MinSideLength``; Python does not read them),
``OutStructureBasename`` and ``SampleRadius 2.0`` (was 0.50: with 0.50 the C++ queue skipped the
voxels farther than 0.50 mm from its moving centre, about 10,000 of the 98,280 on ``fine7``; 2.0
skips none, see PLAN "Owner task of 24 Nov"; the Python run ignores the key).

  CSP  stock C++ in parallel: the PS config (``ar2/key_off_fused.config``) through
       ``cpptime/make_config.py`` (``convert``, imported unchanged: it drops the three Python-only
       lines and makes the paths absolute), mpirun -n 13: one server and 12 fitting clients
  RCP  the recipe of 24 Nov in library form, one worker: ``ar2/key_on_fused.config``
       (SeedByIndexing 1) plus the 12 keys of AR5 step 21's ``lib_check.json`` plus the four keys of
       steps 22, 28 and 29 (``SeedIndexMCDefer 1``, ``RayIndexMaxPairs 50000``,
       ``RayIndexFallbackFull 1``, ``SeedReuseKnown 1``)
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
SAMPLE_RADIUS = "2.0"
ARMS = ("CSP", "RCP")  # the order the launcher runs them
PYTHON_ARMS = ("RCP",)
CPP_ARMS = ("CSP",)
MPI_RANKS = {"CSP": 13}  # server + 12 fitting clients

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
STEP_KEYS: List[str] = [  # AR5 steps 22, 28, 29
    "SeedIndexMCDefer 1",
    "RayIndexMaxPairs 50000",
    "RayIndexFallbackFull 1",
    "SeedReuseKnown 1",
]
RCP_KEYS: List[str] = J5S_KEYS + STEP_KEYS
CHANGED_KEYS = (
    "SampleFilename",
    "InfileBasename",
    "MinSideLength",
    "MaxInitSideLength",
    "OutStructureBasename",
    "SampleRadius",
)


def _set(text: str, key: str, value: str) -> str:
    """Replace the value of the one line starting with ``key`` (keeps nothing else of the line)."""
    new, n = re.subn(rf"(?m)^{key}\b.*$", f"{key} {value}", text)
    assert n == 1, (key, n)
    return new


def python_config(arm: str, sample: Path, frames_basename: str, out_basename: str = "tex1") -> str:
    """The Python config text of ``RCP``."""
    assert arm in PYTHON_ARMS, arm
    text = BASE_J5S.read_text()
    text = _set(text, "SampleFilename", str(sample))
    text = _set(text, "InfileBasename", frames_basename)
    text = _set(text, "MinSideLength", SIDE)
    text = _set(text, "MaxInitSideLength", SIDE)
    text = _set(text, "OutStructureBasename", out_basename)
    text = _set(text, "SampleRadius", SAMPLE_RADIUS)
    text = text.replace("# END OF FILE", "")
    text = (
        text.rstrip("\n") + "\n\n# tex1 arm RCP: AR5 step 21's library keys and steps 22, 28, 29\n"
    )
    text += "\n".join(RCP_KEYS) + "\n\n# END OF FILE\n"
    return text


def cpp_config(sample: Path, frames_basename: str, out_basename: str = "tex1_cpp") -> str:
    """The C++ config text (arm CSP): the ``PS`` config through ``make_config.convert``."""
    text = MC.convert(BASE_PS.read_text(), sample, out_basename)
    text = _set(text, "InfileBasename", frames_basename)
    text = _set(text, "MinSideLength", SIDE)
    text = _set(text, "MaxInitSideLength", SIDE)
    text = _set(text, "SampleRadius", SAMPLE_RADIUS)
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
    """``tex1_config.py --arm CSP --grid slice --out FILE``: write one arm's config text."""
    import argparse

    import tex1_paths as P

    ap = argparse.ArgumentParser(description=main.__doc__)
    ap.add_argument("--arm", required=True, choices=ARMS)
    ap.add_argument("--grid", required=True, choices=("slice", "full"))
    ap.add_argument("--out", type=Path, required=True)
    a = ap.parse_args()
    grid = P.SLICE if a.grid == "slice" else P.GRID
    frames = str(P.FRAMES / "500Grains.sim")
    if a.arm in CPP_ARMS:
        text = cpp_config(grid, frames, "tex1_cpp")
    else:
        text = python_config(a.arm, grid, frames, f"tex1_{a.arm}")
    a.out.parent.mkdir(parents=True, exist_ok=True)
    a.out.write_text(text)
    print(f"wrote {a.out}")


if __name__ == "__main__":
    main()
