"""AR3 step 10 (task ar3j): shared set-up of the validation runs of the library keys.

Runs ``D4L``, ``S5L`` and ``S8L`` use the sample ``big`` with clean images, Q-max 8,
``key_on_fused.config`` plus the keys of the table below, ``default_rng(0)``, one worker each.
Stored runs (untracked; in the worktrees of steps 7 and 8, found from this worktree's parent
folder) are read, never written. Results go to ``benchmarks/autoresearch/ar3/keys/``; ``cache/``
holds the large files (not committed).

Nothing here reads the truth except ``libkeys_score.py``; the library keys do not.
"""

import os

os.environ["AR3_SAMPLE"] = "big"  # before ar3_common is imported

import gzip  # noqa: E402
import json  # noqa: E402
from pathlib import Path  # noqa: E402
from typing import Any, Dict, List  # noqa: E402

import ar3_common as C  # noqa: E402
import numpy as np  # noqa: E402

WORKTREES = C.ICE.parent.parent  # sibling worktrees hold the untracked stored runs
OUT = C.ICE / "benchmarks" / "autoresearch" / "ar3" / "keys"
CACHE = Path(os.environ.get("KEYS_CACHE", str(OUT / "cache")))
_W7 = WORKTREES / "agent-a38ba6f5a25aa033e/icenine_py/benchmarks/autoresearch/ar3/big"
_W8 = WORKTREES / "agent-a635b2f904279bdc9/icenine_py/benchmarks/autoresearch/ar3/reseed/cache"
N = 24570

# keys added to key_on_fused.config (the plan's table)
D4L_KEYS = ["BFSNeighborCMAFirst 1", "CMASigma0 0.05", "CMANeighborMaxEvals 60", "BoundaryPass 1"]
KEYS: Dict[str, List[str]] = {
    "D4L": D4L_KEYS,
    "S5L": D4L_KEYS + ["SeedIndexCandidates 20 5"],
    "S8L": D4L_KEYS + ["SeedIndexCandidates 20 5", "SeedIndexStopHitRatio 0.95"],
}
# stored run of the same name: folder with run.json and recon.mic.gz (the map before the pass)
STORED = {
    "D4L": _W7 / "fine/runs/D4",
    "S5L": _W7 / "seed/runs/S5",
    "S8L": _W8 / "S8_big",
}
# stored pass of each map: folder with tries.npz, and the map before / after as .npy
STORED_PASS = {
    "D4L": (_W7 / "boundary/cache/passes/D4_D4", _W7 / "boundary/cache/passes/before_D4.npy",
            _W7 / "boundary/cache/passes/after_D4.npy"),
    "S5L": (_W7 / "boundary/cache/passes/D4_S5", _W7 / "boundary/cache/passes/before_S5.npy",
            _W7 / "boundary/cache/passes/after_S5.npy"),
    "S8L": (_W8 / "passes_big/D4_S8", _W8 / "passes_big/before_S8.npy",
            _W8 / "passes_big/after_S8.npy"),
}  # fmt: skip
# written down before the runs (plan, step 10): exact for the BFS before the pass
EXPECT = {
    "D4L": dict(seeds=38, reseeds=291, grid=8, evals=10_811_715, wrong_before=68,
                pass_tries=18_557, pass_evals=1_131_977, wrong_after=21),
    "S5L": dict(seeds=35, reseeds=314, grid=0, evals=11_867_509, wrong_before=71,
                pass_tries=18_563, pass_evals=1_132_343, wrong_after=30),
    "S8L": dict(seeds=37, reseeds=292, grid=0, evals=7_201_863, wrong_before=67,
                pass_tries=18_665, pass_evals=1_138_565, wrong_after=24),
}  # fmt: skip


def read_gz(p: Path) -> Any:
    with gzip.open(p, "rt") as f:
        return json.load(f)


def write_gz(obj: Any, p: Path) -> None:
    with gzip.open(p, "wt", compresslevel=6) as f:
        f.write(json.dumps(obj, separators=(",", ":"), default=float))


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
