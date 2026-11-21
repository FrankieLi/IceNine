"""AR3 step 7: the boundary pass (step 5, rule ``hit``) on the map of one run on ``big``.

  uv run python scripts/autoresearch/ar3/big_pass.py --kind {D,D4,S5,P4R4} [--limit N]

The procedure is ``boundary_pass.py`` unchanged (called with the run's map as its input; one
generator ``default_rng(0)``, one worker), as ``seed_pass.py`` calls it for step 6; then rule
``hit`` of ``boundary_lib.apply_rule`` gives the map after the pass. The map of ``P4R4`` is the
painted map with the R4 fits put in (``boundary_common.load_map``). Writes
``big/boundary/cache/passes/``: ``before_<kind>.npy``, ``after_<kind>.npy`` and the pass folder
``D4_<kind>/`` (``tries.npz``, ``pass.json``). With ``--limit`` (a pilot), the folder is
``pilot_...`` and the after map is not written.
"""

import os

os.environ["AR3_SAMPLE"] = "big"

import argparse  # noqa: E402
import sys  # noqa: E402

import ar3_common as C  # noqa: E402,F401
import boundary_common as B  # noqa: E402
import boundary_lib as BL  # noqa: E402
import boundary_pass as BP  # noqa: E402
import numpy as np  # noqa: E402

PASSES = B.CACHE / "passes"
N = 24570


def load(kind: str) -> np.ndarray:
    if kind == "S5":
        return B.G.read_map(C.OUT / "seed" / "runs" / "S5" / "recon.mic.gz", N)[0]
    return B.load_map(kind, N)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--kind", required=True, choices=["D", "D4", "S5", "P4R4"])
    ap.add_argument("--side-by-side", type=int, default=1)
    ap.add_argument("--limit", type=int, default=0)
    a = ap.parse_args()
    PASSES.mkdir(parents=True, exist_ok=True)
    B.RUNS = PASSES
    R = load(a.kind)
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
    if a.limit:
        sys.argv += ["--limit", str(a.limit)]
    BP.main()
    if a.limit:
        return
    t = BL.Tries(PASSES / f"D4_{a.kind}")
    after, took = BL.apply_rule(t, R, "hit")
    np.save(PASSES / f"after_{a.kind}.npy", after)
    print("after map saved; voxels that took a try:", int(took.sum()), flush=True)


if __name__ == "__main__":
    main()
