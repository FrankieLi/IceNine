"""AR3 step 7: per-voxel error of each run's map, before and after the pass, against the truth.

  AR3_SAMPLE=big uv run python scripts/autoresearch/ar3/big_errors.py

Error: the symmetry-reduced misorientation (degrees) to the voxel's own truth
(``gradient_common.Truth.err``); "wrong" is an error above ``gradient_common.WRONG`` (1 degree), as
in steps 5 and 6. Reads the maps in ``big/`` and the pass maps in ``big/boundary/cache/passes/``;
writes ``big/errors.npz`` (float32 errors of ``<kind>_before`` and ``<kind>_after``, the level and
grain of each voxel), which is all ``big_summary.py`` reads of the maps.
"""

import os

os.environ["AR3_SAMPLE"] = "big"

from typing import Dict  # noqa: E402

import ar3_common as C  # noqa: E402
import boundary_common as B  # noqa: E402
import gradient_common as G  # noqa: E402
import numpy as np  # noqa: E402
from big_pass import PASSES, load  # noqa: E402


def main() -> None:
    T = G.Truth()
    out: Dict[str, np.ndarray] = dict(level=T.level.astype(np.float32), grain=T.gmap)
    for kind in ("D", "D4", "S5", "P4R4"):
        R = load(kind)
        before = np.load(PASSES / f"before_{kind}.npy")
        assert np.array_equal(R, before)  # the pass started from this map
        out[f"{kind}_before"] = T.err(R).astype(np.float32)
        out[f"{kind}_after"] = T.err(np.load(PASSES / f"after_{kind}.npy")).astype(np.float32)
        print(
            kind,
            int((out[f"{kind}_before"] > G.WRONG).sum()),
            "->",
            int((out[f"{kind}_after"] > G.WRONG).sum()),
            flush=True,
        )
    np.savez_compressed(C.OUT / "errors.npz", **out)
    assert B.G is G


if __name__ == "__main__":
    main()
