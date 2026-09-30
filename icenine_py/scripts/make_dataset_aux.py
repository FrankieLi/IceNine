#!/usr/bin/env python3
"""Build the aux table that goes with a toy orientation dataset (no windows are re-rendered).

Contents (padded to the dataset's peak count; multi-voxel tables have a leading voxel axis):
  nom_off      (.., M, 3) float32  exact nominal offsets: fractional centroid col/row (px) and the
                                   nominal omega relative to its frame centre (frames)
  pair_index   (.., M)    int64    partner entry of the same ray on the other detector, -1 if none
  det_idx      (.., M)    int64    detector of each entry (-1 for padding)
  frame_width_rad         float    |frame width| of the rotation scan (rad)
  voxel_indices / n_peaks_per_voxel copied from the dataset, used to check the aux matches it

Usage:
  cd icenine_py
  uv run python scripts/make_dataset_aux.py --data scripts/toy_orientation_stage3_multi_test.pt \
      --out scripts/toy_orientation_stage3_multi_aux.pt
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import torch

sys.path.insert(0, str(Path(__file__).parent))
from dataset_problems import iter_voxels, load_meta  # noqa: E402


def main():
    from icenine.orientation_eval import nominal_offsets, pair_index

    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--data", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    out = Path(args.out).resolve()
    data = load_meta(args.data)
    multi = bool(data.get("multi_voxel", False))
    M = int(data["n_peaks"])
    V = len(data["voxel_indices"]) if multi else 1
    nom = np.zeros((V, M, 3), dtype=np.float32)
    pidx = np.full((V, M), -1, dtype=np.int64)
    didx = np.full((V, M), -1, dtype=np.int64)
    fw = None
    for vox in iter_voxels(data):
        v, n = vox["v"], vox["n_pk"]
        nom[v, :n] = nominal_offsets(vox["obs"])
        pidx[v, :n] = pair_index(vox["problem"]["roi_list"])
        didx[v, :n] = vox["obs"].det_idx.numpy()
        fw = vox["obs"].frame_width_rad
        print(f"voxel {v}: {n} entries, {(pidx[v] >= 0).sum()} paired")
    aux = dict(
        nom_off=torch.from_numpy(nom if multi else nom[0]),
        pair_index=torch.from_numpy(pidx if multi else pidx[0]),
        det_idx=torch.from_numpy(didx if multi else didx[0]),
        frame_width_rad=float(fw),
        multi_voxel=multi,
        n_peaks=M,
        voxel_indices=data["voxel_indices"] if multi else torch.tensor([data["voxel_index"]]),
    )
    out.parent.mkdir(parents=True, exist_ok=True)
    torch.save(aux, out)
    print(f"saved {out}")


if __name__ == "__main__":
    main()
