"""Build the image stores a grain fit reads, from a tag's ASCII clean frames.

Writes (under ``Examples/Example2.ManyGrains/ScatteringData_PhaseD``, git-ignored):

* ``stacks/<tag>_clean.npy``: the uint8 binary stack (frame ``det * n_omega + omega``), the same
  layout ``write_binary_stack`` writes and ``ExperimentalData.from_binary_memmap`` reads.
* ``<tag>/intensity.npz``: the clean intensities as a sparse store: ``frame`` offsets (CSR over
  the flat frame index), ``row``, ``col`` (uint16) and ``val`` (float32).

Usage (from icenine_py/):
  uv run python scripts/grain_opt/build_store.py --tag grad [--variant clean]
"""

import argparse
import json
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
PHASE_D = ROOT / "Examples" / "Example2.ManyGrains" / "ScatteringData_PhaseD"
BASENAME = "500Grains.sim"
N_OMEGA, N_DET, NROW, NCOL = 180, 2, 2048, 2048


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--tag", required=True)
    ap.add_argument("--variant", default="clean")
    ap.add_argument("--out_tag", default=None, help="write stacks/<out_tag>_clean.npy and "
                    "<out_tag>/intensity_clean.npz instead (GO5m: a noisy variant read as a clean tag)")
    a = ap.parse_args()
    src = PHASE_D / a.tag / a.variant
    stack_path = PHASE_D / "stacks" / f"{a.tag}_{a.variant}.npy"
    store_path = PHASE_D / a.tag / f"intensity_{a.variant}.npz"
    if a.out_tag:
        stack_path = PHASE_D / "stacks" / f"{a.out_tag}_clean.npy"
        store_path = PHASE_D / a.out_tag / "intensity_clean.npz"
        store_path.parent.mkdir(parents=True, exist_ok=True)
    stack_path.parent.mkdir(parents=True, exist_ok=True)
    stack = np.lib.format.open_memmap(
        str(stack_path), mode="w+", dtype=np.uint8, shape=(N_OMEGA * N_DET, NROW, NCOL)
    )
    offs = [0]
    rows, cols, vals = [], [], []
    for d in range(N_DET):
        for o in range(N_OMEGA):
            fn = src / f"{BASENAME}{str(o).zfill(5)}.d{d}"
            arr = np.loadtxt(fn, delimiter=",", ndmin=2) if fn.stat().st_size else np.zeros((0, 3))
            if len(arr):
                arr = arr[arr[:, 2] > 0]
            j = arr[:, 0].astype(np.int64)
            k = arr[:, 1].astype(np.int64)
            frame = np.zeros((NROW, NCOL), dtype=np.uint8)
            frame[k, j] = 1
            stack[d * N_OMEGA + o] = frame
            order = np.lexsort((j, k))
            rows.append(k[order].astype(np.uint16))
            cols.append(j[order].astype(np.uint16))
            vals.append(arr[order, 2].astype(np.float32))
            offs.append(offs[-1] + len(order))
    stack.flush()
    del stack
    np.savez(
        store_path,
        frame=np.asarray(offs, dtype=np.int64),
        row=np.concatenate(rows),
        col=np.concatenate(cols),
        val=np.concatenate(vals),
    )
    info = dict(tag=a.tag, variant=a.variant, stack=str(stack_path), store=str(store_path),
                nonzero=int(offs[-1]), val_max=float(np.concatenate(vals).max()))
    print(json.dumps(info, indent=1))


if __name__ == "__main__":
    main()
