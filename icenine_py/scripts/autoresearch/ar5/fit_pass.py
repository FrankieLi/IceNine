"""AR5 step 2, Part 2: the boundary pass (step 5, rule ``hit``) on the map of one run of this step.

  scripts/autoresearch/ar5/fit_bg.sh <sample> pass_<name>_<sample> fit_pass.py --name S10

The procedure is step 9's ``start_pass.py`` (itself step 8's ``reseed_pass.py``) unchanged:
``boundary_pass.py`` called with the run's map as its input, one generator ``default_rng(0)``, one
worker, every try a CMA fit (so the pass has no MC fallback to remove); then rule ``hit`` of
``boundary_lib.apply_rule`` gives the map after the pass. Only the run's name differs. Reads
``<cache>/<name>_<sample>/recon.mic.gz``; writes ``<cache>/passes_<sample>/``:
``before_<name>.npy``, ``after_<name>.npy`` and the folder ``D4_<name>/`` (``tries.npz``,
``pass.json``).
"""

import argparse
import sys

import fit_common as FC
import numpy as np

import boundary_common as B  # noqa: E402
import boundary_lib as BL  # noqa: E402
import boundary_pass as BP  # noqa: E402

PASSES = FC.CACHE / f"passes_{FC.SAMPLE}"
N = 24570


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--name", required=True)
    ap.add_argument("--side-by-side", type=int, default=1)
    a = ap.parse_args()
    PASSES.mkdir(parents=True, exist_ok=True)
    B.RUNS = PASSES
    R, _ = B.G.read_map(FC.CACHE / f"{a.name}_{FC.SAMPLE}" / "recon.mic.gz", N)
    npy = PASSES / f"before_{a.name}.npy"
    np.save(npy, R)
    sys.argv = [
        "boundary_pass",
        "--map",
        "D4",
        "--input",
        str(npy),
        "--tag",
        a.name,
        "--side-by-side",
        str(a.side_by_side),
    ]
    BP.main()
    t = BL.Tries(PASSES / f"D4_{a.name}")
    after, took = BL.apply_rule(t, R, "hit")
    np.save(PASSES / f"after_{a.name}.npy", after)
    print("after map saved; voxels that took a try:", int(took.sum()), flush=True)


if __name__ == "__main__":
    main()
