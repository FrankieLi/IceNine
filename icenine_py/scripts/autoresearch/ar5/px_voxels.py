"""AR5 step 27: the voxel lists (written once, hashed in the freeze).

  RESEED_SAMPLE=mid PX_PRE_FREEZE=1 uv run python scripts/autoresearch/ar5/px_voxels.py --write

``mid`` and ``held`` (t1): every stored seed step's voxel of step 22's J5D run (run order,
repeats kept in place). q1: 300 stored seed steps of each of the eight disclosed samples, chosen
with ``numpy.random.default_rng(27).choice(n_steps, 300, replace=False)`` and sorted, their
voxels. Only ``e1_log`` rows (voxel, reseed) of the disclosed samples are read.
"""

import json
import sys
from typing import Any, Dict, List

import px_common as PC
import px_rules as PR

import numpy as np  # noqa: E402

OUT = PC.RD / "voxels.json"


def build() -> Dict[str, Any]:
    out: Dict[str, Any] = dict(
        generator=PR.SEED_GENERATOR, n_per_sample=PR.N_PER_SAMPLE, q1={}, t1={}
    )
    for s in PR.Q1_SAMPLES:
        steps = PC.stored_steps(s)
        rng = np.random.default_rng(PR.SEED_GENERATOR)
        sel = sorted(int(i) for i in rng.choice(len(steps), PR.N_PER_SAMPLE, replace=False))
        out["q1"][s] = dict(n_steps=len(steps), steps=sel, voxels=[steps[i][0] for i in sel])
        if s in PR.FULL_SAMPLES:
            allsel: List[int] = list(range(len(steps)))
            out["t1"][s] = dict(n_steps=len(steps), voxels=[steps[i][0] for i in allsel])
    return out


if __name__ == "__main__":
    if "--write" not in sys.argv:
        raise SystemExit("use --write")
    d = build()
    OUT.write_text(json.dumps(d) + "\n")
    print({s: (v["n_steps"], len(v["voxels"])) for s, v in d["q1"].items()})
    print({s: (v["n_steps"], len(v["voxels"])) for s, v in d["t1"].items()})
    print("q1 voxels:", sum(len(v["voxels"]) for v in d["q1"].values()))
    print("t1 voxels:", sum(len(v["voxels"]) for v in d["t1"].values()))
