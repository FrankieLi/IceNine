"""AR3 step 6, part 2: the boundary pass (step 5, rule ``hit``) on the map of run S or K2.

  seed_pass.py --kind S

The procedure is ``boundary_pass.py`` unchanged (the same code, called with the run's map as its
input; its output folder is moved to ``seed/cache/passes/``), then rule ``hit`` of
``boundary_lib.apply_rule`` gives the map after the pass, written to
``seed/cache/passes/after_<kind>.npy``. Step 5's pass on D4 and on step 3's D were made the same
way (one generator ``default_rng(0)``, one worker).
"""

import argparse
import sys

import ar3_common as C  # noqa: F401
import boundary_common as B
import boundary_lib as BL
import boundary_pass as BP
import numpy as np
import seed_common as S

PASSES = S.CACHE / "passes"


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--kind", required=True, choices=["S", "K2"])
    ap.add_argument("--side-by-side", type=int, default=1)
    a = ap.parse_args()
    PASSES.mkdir(parents=True, exist_ok=True)
    B.RUNS = PASSES
    n = 24570
    R, _ = B.G.read_map(S.SEED / "runs" / a.kind / "recon.mic.gz", n)
    npy = PASSES / f"before_{a.kind}.npy"
    np.save(npy, R)
    sys.argv = [
        "boundary_pass",
        "--map",
        "D4",
        "--input",
        str(npy),
        "--tag",
        a.kind,
        "--side-by-side",
        str(a.side_by_side),
    ]
    BP.main()
    t = BL.Tries(PASSES / f"D4_{a.kind}")
    after, took = BL.apply_rule(t, R, "hit")
    np.save(PASSES / f"after_{a.kind}.npy", after)
    print("after map saved; voxels that took a try:", int(took.sum()), flush=True)


if __name__ == "__main__":
    main()
