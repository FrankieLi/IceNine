"""AR5 step 7 (task ar5g): where the stored passes are, and how to score a map.

The stored passes, run folders and maps are untracked files of other (kept) worktrees, found from
this worktree's parent folder; they are read, never written. Nothing here imports from the
unmerged step 4 and step 6 task branches; it only reads the files they left.

* ``mix`` (step 3, worktree ``agent-a5a3601b746680640``): run ``S11`` on ``big``, ``grad``, ``mid``.
* ``ver`` (step 6, worktree ``agent-a89458addd73332d9``): run ``S11n`` on ``mid``, ``grad``.
* ``held`` (step 4, worktrees ``agent-a814d0121e137b01c`` and ``agent-af8961874e131bda2``): run
  ``S11v2`` and the baseline ``S8_held``. Its entries are read only after ``frozen.json`` is
  committed (``gate_held.py``).

A stored pass folder has ``D4_<run>/tries.npz`` (voxel, cand, end, cost, hit, conf, evals per try;
``old_score`` (n, 3) = cost, hit ratio, confidence of every voxel's end before the pass; ``old_R``)
and ``pass.json``, and ``before_<run>.npy`` / ``after_<run>.npy`` (the maps, (n, 3, 3)).
Truth is used for scoring only (error = symmetry-reduced angle, cubic, to the truth ``.mic``).
"""

import hashlib
import json
import sys
from pathlib import Path
from typing import Any, Dict, List, Tuple

import numpy as np

HERE = Path(__file__).resolve().parent
ICE = HERE.parents[2]  # icenine_py
sys.path.insert(0, str(ICE / "scripts" / "common"))  # stats.py
WORKTREES = ICE.parent.parent
AR5 = "icenine_py/benchmarks/autoresearch/ar5"
AR3 = ICE / "benchmarks" / "autoresearch" / "ar3"
GATE = ICE / "benchmarks" / "autoresearch" / "ar5" / "gate"
MIX_TRACKED = ICE / "benchmarks" / "autoresearch" / "ar5" / "mix"
MIX_CACHE = WORKTREES / f"agent-a5a3601b746680640/{AR5}/mix/cache"
VER_CACHE = WORKTREES / f"agent-a89458addd73332d9/{AR5}/ver/cache"
N = 24570
WRONG = 1.0  # degrees

# (sample, run) -> (run folder with run.json, pass cache folder with passes_<sample>/)
DISCLOSED: List[Tuple[str, str]] = [
    ("mid", "S11n"),
    ("grad", "S11n"),
    ("big", "S11"),
    ("grad", "S11"),
    ("mid", "S11"),
]
# the plan's totals of the stored passes (tries only), written with the plan draft
PLAN_PASS_EVALS: Dict[Tuple[str, str], int] = {
    ("mid", "S11n"): 3_059_211,
    ("grad", "S11n"): 3_684_888,
    ("big", "S11"): 1_135_637,
    ("grad", "S11"): 3_685_803,
    ("mid", "S11"): 3_048_170,
}
PLAN_CUT_PCT: Dict[Tuple[str, str], float] = {
    ("mid", "S11n"): 68.9,
    ("grad", "S11n"): 83.3,
    ("big", "S11"): 75.9,
    ("grad", "S11"): 83.5,
    ("mid", "S11"): 69.2,
}


def cache_of(run: str) -> Path:
    return VER_CACHE if run == "S11n" else MIX_CACHE


def run_dir(sample: str, run: str) -> Path:
    return cache_of(run) / f"{run}_{sample}"


def pass_dir(sample: str, run: str) -> Path:
    return cache_of(run) / f"passes_{sample}"


def pass_files(sample: str, run: str) -> Dict[str, Path]:
    pd = pass_dir(sample, run)
    return dict(
        tries=pd / f"D4_{run}" / "tries.npz",
        pass_json=pd / f"D4_{run}" / "pass.json",
        before=pd / f"before_{run}.npy",
        after=pd / f"after_{run}.npy",
        run_json=run_dir(sample, run) / "run.json",
    )


def sha256(p: Path) -> str:
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for blk in iter(lambda: f.read(1 << 20), b""):
            h.update(blk)
    return h.hexdigest()


def truth_mic(sample: str) -> Path:
    if sample == "grad":
        return AR3 / "gradient" / "sample" / "gradient_truth.mic"
    return AR3 / sample / "gradient" / "sample" / "gradient_truth.mic"


class Truth:
    """Truth orientations of a sample and the symmetry-reduced error of a map."""

    def __init__(self, sample: str) -> None:
        from icenine.mic_file import MicFile
        from icenine.ray_indexing import proper_rotations
        from icenine.symmetry import create_cubic_symmetry

        self.sym = proper_rotations(create_cubic_symmetry(3.6))
        mic = MicFile.read(str(truth_mic(sample)))
        self.R = np.array([v.orientation for v in mic.voxels], dtype=np.float64)

    def err(self, R: np.ndarray) -> np.ndarray:
        from icenine.ray_indexing import misorientation_deg

        return np.asarray(misorientation_deg(R, self.R, self.sym), dtype=float)


def load_pass(sample: str, run: str) -> Dict[str, Any]:
    f = pass_files(sample, run)
    s = np.load(f["tries"])
    run_j = json.loads(f["run_json"].read_text())
    pj = json.loads(f["pass_json"].read_text())
    return dict(
        voxel=s["voxel"].astype(np.int64),
        hit=s["hit"],
        cost=s["cost"],
        evals=s["evals"].astype(np.int64),
        end=s["end"],
        old_score=s["old_score"],
        before=np.load(f["before"]).astype(np.float64),
        after=np.load(f["after"]).astype(np.float64),
        bfs_evals=int(run_j["evals_total"]),
        pass_json=pj,
    )


def s8_baseline(sample: str) -> Dict[str, Any]:
    """S8's errors (before and after its pass) and evaluation totals on a disclosed sample, from
    step 3's tracked ``part2_S11_<sample>*`` files (they hold the stored S8 run's numbers)."""
    e = np.load(MIX_TRACKED / f"part2_S11_{sample}_err.npz")
    j = json.loads((MIX_TRACKED / f"part2_S11_{sample}.json").read_text())
    return dict(
        err_before=np.asarray(e["s8_before"], dtype=float),
        err_after=np.asarray(e["s8_after"], dtype=float),
        stages=j["stages_s8"],
    )
