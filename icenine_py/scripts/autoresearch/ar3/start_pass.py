"""AR3 step 9, part 3: the boundary pass (step 5, rule ``hit``) on the map of one whole-layer run.

  scripts/autoresearch/ar3/start_bg.sh <sample> pass_<kind> start_pass.py --kind D4|S5|S8|S9

The procedure is ``reseed_pass.py`` unchanged (``boundary_pass.py`` called with the run's map as
its input; one generator ``default_rng(0)``, one worker), with the run chosen by ``--kind``; then
rule ``hit`` of ``boundary_lib.apply_rule`` gives the map after the pass. Reads
``<RESEED_CACHE>/<kind>_<sample>/recon.mic.gz``; writes ``<RESEED_CACHE>/passes_<sample>/``:
``before_<kind>.npy``, ``after_<kind>.npy`` and the pass folder ``D4_<kind>/`` (``tries.npz``,
``pass.json``).
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
    ap.add_argument("--kind", required=True, choices=["D4", "S5", "S8", "S9"])
    ap.add_argument("--side-by-side", type=int, default=1)
    a = ap.parse_args()
    PASSES.mkdir(parents=True, exist_ok=True)
    B.RUNS = PASSES
    R, _ = B.G.read_map(RC.CACHE / f"{a.kind}_{RC.SAMPLE}" / "recon.mic.gz", N)
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
