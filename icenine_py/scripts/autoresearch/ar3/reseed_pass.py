"""AR3 step 8, part 3: the boundary pass (step 5, rule ``hit``) on the map of run S8.

  RESEED_SAMPLE=big  uv run python scripts/autoresearch/ar3/reseed_pass.py
  RESEED_SAMPLE=grad uv run python scripts/autoresearch/ar3/reseed_pass.py

The procedure is ``boundary_pass.py`` unchanged (called with the run's map as its input; one
generator ``default_rng(0)``, one worker), as ``seed_pass.py`` and ``big_pass.py`` call it; then
rule ``hit`` of ``boundary_lib.apply_rule`` gives the map after the pass. Reads
``reseed/cache/S8_<sample>/recon.mic.gz``; writes ``reseed/cache/passes_<sample>/``:
``before_S8.npy``, ``after_S8.npy`` and the pass folder ``D4_S8/`` (``tries.npz``, ``pass.json``).
"""

import argparse
import sys

import reseed_common as RC

import boundary_common as B  # noqa: E402
import boundary_lib as BL  # noqa: E402
import boundary_pass as BP  # noqa: E402
import numpy as np  # noqa: E402

PASSES = RC.CACHE / f"passes_{RC.SAMPLE}"
N = 24570


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--side-by-side", type=int, default=1)
    a = ap.parse_args()
    PASSES.mkdir(parents=True, exist_ok=True)
    B.RUNS = PASSES
    R, _ = B.G.read_map(RC.CACHE / f"S8_{RC.SAMPLE}" / "recon.mic.gz", N)
    npy = PASSES / "before_S8.npy"
    np.save(npy, R)
    sys.argv = [
        "boundary_pass",
        "--map",
        "D4",
        "--input",
        str(npy),
        "--tag",
        "S8",
        "--side-by-side",
        str(a.side_by_side),
    ]
    BP.main()
    t = BL.Tries(PASSES / "D4_S8")
    after, took = BL.apply_rule(t, R, "hit")
    np.save(PASSES / "after_S8.npy", after)
    print("after map saved; voxels that took a try:", int(took.sum()), flush=True)


if __name__ == "__main__":
    main()
