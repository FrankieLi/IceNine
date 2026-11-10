"""AR2 step 2: do the same indexing and strict cost reproduce AR1's stored candidates?

Some voxels of the BFS pilot gave different tail candidates than AR1's stored ones. This script
indexes eight pilot voxels again in AR1's setting (the full sample, ``ar1_common.env``) and compares
the candidates and their strict costs with AR1's stored ones (``cand_pilot_base.npz`` and
``e2e_pilot_base_part0of1.json``). The comparison with the BFS setting (the region grid as the
sample, so a smaller hull for the ray pairing) is in ``key_summary.py`` (``ar1_comparison``).

  uv run python scripts/autoresearch/ar2/key_thread_check.py
Writes benchmarks/autoresearch/ar2/key/thread_check_default.json (the tag "one" is unused: a
single-thread variant was not run).
"""

import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[0] / "ar1"))
import ar1_common as AC  # noqa: E402

from icenine import ray_indexing as RI  # noqa: E402

OUT = AC.ICE / "benchmarks" / "autoresearch" / "ar2" / "key"
VOXELS = (24281, 9061, 5319, 19480, 2032, 983, 4206, 11294)


def main() -> None:
    tag = "one" if "--one" in sys.argv else "default"
    import torch

    if tag == "one":
        torch.set_num_threads(1)
    z = np.load(AC.OUT / "cand_pilot_base.npz")
    rows = json.loads((AC.OUT / "e2e_pilot_base_part0of1.json").read_text())["rows"]
    ar1 = {int(r["voxel"]): r for r in rows if not r.get("no_candidates")}
    p = RI.IndexParams()
    AC.init("full")
    out = []
    for v in VOXELS:
        e = AC.env("full", v)
        sc = e["sc"]
        crystal = RI.crystal_from_smooth_cost(sc, AC.symmetry_of(e), e["phase"], p.q_max)
        res, att = RI.index_voxel(sc, e["verts"], crystal, p, 20)
        R0 = z[f"v{v}"]
        same_n = res.R.shape == R0.shape
        dmax = float(np.abs(res.R - R0).max()) if same_n else None
        lf = e["lf"]
        costs = [float(lf.evaluate(res.R[j], e["verts"], e["phase"]).cost) for j in range(5)]
        out.append(
            dict(
                voxel=v,
                n_att=int(len(att)),
                n_cand=int(len(res.R)),
                n_cand_ar1=int(len(R0)),
                max_abs_diff_R=dmax,
                same_strict_first5=bool(np.array_equal(costs, ar1[v]["strict_costs"][:5])),
                n_first5_same=int(sum(a == b for a, b in zip(costs, ar1[v]["strict_costs"][:5]))),
                torch_threads=torch.get_num_threads(),
            )
        )
        print(out[-1], flush=True)
    (OUT / f"thread_check_{tag}.json").write_text(json.dumps(out, indent=1) + "\n")


if __name__ == "__main__":
    main()
