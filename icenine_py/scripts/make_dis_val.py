"""Build the Huber-c validation set from the Step 4 distractor training set.

Samples 50 training samples for each of voxels 6, 11, 19, 24 (numpy seed 0) from
``toy_orientation_arch_dis_train.pt`` and writes them as a single-group dataset
(``magnitudes_deg`` = 0).  Used to select the Huber c of the robust GN baseline.
"""

import argparse

import numpy as np
import torch


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--train", default="scripts/toy_orientation_arch_dis_train.pt")
    ap.add_argument("--out", default="scripts/toy_orientation_arch_dis_val.pt")
    ap.add_argument("--voxels", type=int, nargs="+", default=[6, 11, 19, 24])
    ap.add_argument("--per-voxel", type=int, default=50)
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()

    tr = torch.load(args.train)
    vid = tr["voxel_id"].numpy()
    rng = np.random.default_rng(args.seed)
    idx = np.concatenate(
        [rng.choice(np.nonzero(vid == v)[0], args.per_voxel, replace=False) for v in args.voxels]
    )
    idx.sort()
    keep = ("windows", "dis_windows", "offsets_deg", "voxel_id")
    out = {k: v for k, v in tr.items() if k not in keep}
    for k in keep:
        out[k] = tr[k][idx].clone()
    nrm = tr["offsets_deg"][idx].double().norm(dim=1).numpy()
    print("offset norm percentiles (deg):", np.percentile(nrm, [0, 25, 50, 75, 100]))
    out["magnitudes_deg"] = torch.zeros(len(idx))  # single group
    torch.save(out, args.out)
    print(len(idx), "samples ->", args.out)


if __name__ == "__main__":
    main()
