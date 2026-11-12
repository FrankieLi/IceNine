"""AR5 step 5 (task ar5e): shared set-up of the parity runs of the two library keys.

Library runs (``keys_run.py``) on the samples ``big`` (the 32-grain sample, ``AR3_SAMPLE=big``) and
``grad`` (the first gradient sample, ``AR3_SAMPLE`` unset), clean images, Q-max 8,
``key_on_fused.config`` plus the keys of the table below, ``default_rng(0)``, one worker each.

* ``S8L``  is step 10's ``S8L`` (parity check 1: the keys of S8, the two new keys absent);
* ``S10L`` adds ``BFSNeighborCMANoFallback 1`` (check 2, against the stored script run ``S10``);
* ``S11L`` adds ``SeedIndexCMAFirst 200`` (checks 3 and 4, against the stored script run ``S11``).

The stored runs and passes are untracked files of other worktrees, found from this worktree's
parent folder; they are read, never written. The library never sees the truth; only the check
script does. Results go to ``benchmarks/autoresearch/ar5/keys/`` (``cache/`` holds the large
files and is not committed). Import this module first: it reads ``AR3_SAMPLE`` before the AR3
common modules do.
"""

import gzip
import json
import os
import sys
from pathlib import Path
from typing import Any, Dict, List, Tuple

assert os.environ.get("AR3_SAMPLE", "gradient") in ("big", "gradient"), os.environ["AR3_SAMPLE"]
SAMPLE = "big" if os.environ.get("AR3_SAMPLE") == "big" else "grad"

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "ar3"))

import ar3_common as C  # noqa: E402
import numpy as np  # noqa: E402

ICE = C.ICE
WORKTREES = ICE.parent.parent
OUT = ICE / "benchmarks" / "autoresearch" / "ar5" / "keys"
CACHE = Path(os.environ.get("KEYS5_CACHE", str(OUT / "cache")))

_AR3 = "icenine_py/benchmarks/autoresearch/ar3"
_AR5 = "icenine_py/benchmarks/autoresearch/ar5"
_W8 = WORKTREES / f"agent-a635b2f904279bdc9/{_AR3}/reseed/cache"  # stored S8 (step 8)
_FIT = WORKTREES / f"agent-acc78ac2a483ed86f/{_AR5}/fit/cache"  # stored S10 (AR5 step 2)
_MIX = WORKTREES / f"agent-a5a3601b746680640/{_AR5}/mix/cache"  # stored S11 (AR5 step 3)

# the keys added to key_on_fused.config
S8_KEYS = [
    "BFSNeighborCMAFirst 1",
    "CMASigma0 0.05",
    "CMANeighborMaxEvals 60",
    "BoundaryPass 1",
    "SeedIndexCandidates 20 5",
    "SeedIndexStopHitRatio 0.95",
]
KEYS: Dict[str, List[str]] = {
    "S8L": S8_KEYS,
    "S10L": S8_KEYS + ["BFSNeighborCMANoFallback 1"],
    "S11L": S8_KEYS + ["BFSNeighborCMANoFallback 1", "SeedIndexCMAFirst 200"],
}
# stored run of the same kind and sample (folder with run.json, recon.mic.gz, records.json.gz)
STORED: Dict[Tuple[str, str], Path] = {
    ("S8L", "big"): _W8 / "S8_big",
    ("S10L", "big"): _FIT / "S10_big",
    ("S11L", "big"): _MIX / "S11_big",
    ("S11L", "grad"): _MIX / "S11_grad",
}
# stored pass of the run: (folder with tries.npz, map before as .npy, map after as .npy)
STORED_PASS: Dict[Tuple[str, str], Tuple[Path, Path, Path]] = {
    ("S8L", "big"): (_W8 / "passes_big/D4_S8", _W8 / "passes_big/before_S8.npy",
                     _W8 / "passes_big/after_S8.npy"),
    ("S10L", "big"): (_FIT / "passes_big/D4_S10", _FIT / "passes_big/before_S10.npy",
                      _FIT / "passes_big/after_S10.npy"),
    ("S11L", "big"): (_MIX / "passes_big/D4_S11", _MIX / "passes_big/before_S11.npy",
                      _MIX / "passes_big/after_S11.npy"),
    ("S11L", "grad"): (_MIX / "passes_grad/D4_S11", _MIX / "passes_grad/before_S11.npy",
                       _MIX / "passes_grad/after_S11.npy"),
}  # fmt: skip
# written down in the plan (big): evaluations before the pass, seeds, reseeds
EXPECT: Dict[Tuple[str, str], Dict[str, int]] = {
    ("S8L", "big"): dict(evals=7_201_863, seeds=37, reseeds=292, wrong_before=67),
    ("S10L", "big"): dict(evals=4_767_135, seeds=36, reseeds=297),
    ("S11L", "big"): dict(evals=2_879_705, seeds=35, reseeds=295, with_pass=4_015_342),
    ("S11L", "grad"): dict(with_pass=12_307_994),
}


def read_gz(p: Path) -> Any:
    with gzip.open(p, "rt") as f:
        return json.load(f)


def gz_bytes(path: Path) -> bytes:
    with gzip.open(path, "rb") as f:
        return f.read()


def is_time_key(k: str) -> bool:
    """Wall-clock counters (left out of every equality check)."""
    return "wall" in k or k.endswith("_s") or k in ("index_pairing_s",)


def counter_stats(stats: Dict[str, Any]) -> Dict[str, Any]:
    """The numeric counters of a run's stats without the time counters and the reseed list."""
    return {
        k: v
        for k, v in stats.items()
        if isinstance(v, (int, float)) and not isinstance(v, bool) and not is_time_key(k)
    }


def evals_total(st: Dict[str, Any]) -> int:
    """The run scripts' own total of a BFS (the pass is not in it)."""
    return int(
        sum(st[k] for k in ("n_evals_seed", "n_evals_neighbor", "n_evals_revisit"))
        + st.get("n_evals_reseed", 0)
        + st.get("n_evals_restart", 0)
    )


def orientations(mic: Any) -> np.ndarray:
    return np.stack([np.asarray(v.orientation, np.float64) for v in mic.voxels])
