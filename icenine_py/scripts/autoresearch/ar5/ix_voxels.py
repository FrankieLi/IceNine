"""AR5 step 26: the exactness voxel lists (written once, hashed in the freeze).

  RESEED_SAMPLE=mid IX_PRE_FREEZE=1 uv run python scripts/autoresearch/ar5/ix_voxels.py --write

``mid`` and ``held``: every stored seed step's voxel of step 22's J5D run (run order, repeats
kept in place). ``held3`` to ``held6``, ``grad``, ``big``: 300 stored seed steps each, chosen with
``numpy.random.default_rng(26).choice(n_steps, 300, replace=False)`` and sorted, their voxels.
Only ``e1_log`` rows (voxel, reseed) of the disclosed samples are read; nothing else.
"""

import json
import sys
from typing import Any, Dict, List

import ix_common as IC
import ix_rules as IR

import numpy as np  # noqa: E402

OUT = IC.RD / "voxels.json"


def build() -> Dict[str, Any]:
    out: Dict[str, Any] = dict(
        generator=IR.SEED_GENERATOR, n_per_sample=IR.N_PER_SAMPLE, samples={}
    )
    for s in IR.FULL_SAMPLES + IR.SAMPLED:
        steps = IC.stored_steps(s)
        if s in IR.FULL_SAMPLES:
            sel: List[int] = list(range(len(steps)))
        else:
            rng = np.random.default_rng(IR.SEED_GENERATOR)
            sel = sorted(int(i) for i in rng.choice(len(steps), IR.N_PER_SAMPLE, replace=False))
        vox = [steps[i][0] for i in sel]
        out["samples"][s] = dict(n_steps=len(steps), steps=sel, voxels=vox)
    return out


if __name__ == "__main__":
    if "--write" not in sys.argv:
        raise SystemExit("use --write")
    d = build()
    OUT.write_text(json.dumps(d) + "\n")
    print(
        {
            s: (v["n_steps"], len(v["voxels"]), len(set(v["voxels"])))
            for s, v in d["samples"].items()
        }
    )
    print("total voxels:", sum(len(v["voxels"]) for v in d["samples"].values()))
