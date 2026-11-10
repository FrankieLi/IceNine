"""AR2 step 2: the library path equals the script path, on real voxels, in the BFS pilot's setting.

Builds the pilot's reconstruction setup (the region grid as the sample, as ``pilot_bfs.py``), and
for a few voxels of the stored key-on run's ``index_log.json`` runs (a) ``index_seed`` (the library
path) and (b) the path of ``holdout_run.py`` / ``ar1/e2e.py`` written out here: index, strict cost
on the first 5 candidates, ``refine_from_candidates`` (rng seed 0) on the 2 lowest, the lowest final
cost. They must agree exactly, and with the entry stored in the run's log.

  uv run python scripts/autoresearch/ar2/key_path_check.py [N]
Writes benchmarks/autoresearch/ar2/key/path_check.json.
"""

import json
import os
import re
import sys
from pathlib import Path
from typing import Any, Dict, List

import numpy as np

HERE = Path(__file__).resolve().parent
ICE = HERE.parents[2]
EX = ICE.parent / "Examples" / "Example2.ManyGrains"
OUT = ICE / "benchmarks" / "autoresearch" / "ar2" / "key"
REGION = ICE / "benchmarks" / "phase_d_pilot" / "region_grid.mic"


def main() -> None:
    n_vox = int(sys.argv[1]) if len(sys.argv) > 1 else 6
    cfg_text = (HERE / "key_on.config").read_text()
    cfg_text, k = re.subn(r"(?m)^SampleFilename\s+\S+", f"SampleFilename {REGION}", cfg_text)
    assert k == 1
    tmp = OUT / "cache" / "path_check.config"
    tmp.parent.mkdir(parents=True, exist_ok=True)
    tmp.write_text(cfg_text)
    os.chdir(EX)
    import torch

    torch.set_num_threads(1)
    from icenine import ray_indexing as RI
    from icenine.config_file import ConfigFile
    from icenine.experimental_data import ExperimentalData
    from icenine.orientation_search import SearchCandidate
    from icenine.reconstructor import (
        AdaptiveVoxelReconstructor,
        _get_voxel_vertices,
        setup_reconstruction,
    )

    cfg = ConfigFile.from_file(str(tmp))
    stack = EX / "ScatteringData_PhaseD" / "stacks" / "clean.npy"
    setup = setup_reconstruction(cfg, exp_data=ExperimentalData.from_binary_memmap(stack, 180, 2))
    rec = AdaptiveVoxelReconstructor(setup)
    mic = setup.sample.get_mic()
    log = json.loads((OUT / "runs" / "key_on" / "index_log.json").read_text())
    pick = log[:: max(1, len(log) // n_vox)][:n_vox]
    rows: List[Dict[str, Any]] = []
    for e in pick:
        r = int(e["voxel"])
        voxel = mic.voxels[r]
        verts = _get_voxel_vertices(voxel)
        phase = voxel.phase
        lib = rec.index_seed(verts, phase)
        # the script path, written out
        sc, crystal, p, _ = rec._index_context(phase)
        res, _att = RI.index_voxel(sc, verts, crystal, RI.IndexParams(), 20)
        lf = rec._make_local_cost_fn()
        costs = np.array(
            [float(lf.evaluate(res.R[j], verts, phase).cost) for j in range(5)], dtype=np.float64
        )
        order = [int(j) for j in np.argsort(costs, kind="stable")[:2]]
        fins = []
        for j in order:
            c = rec.refine_from_candidates(
                [SearchCandidate(orientation=res.R[j].copy(), cost=float(costs[j]))],
                verts,
                phase,
                rng=np.random.default_rng(0),
            )
            fins.append(
                dict(
                    cand=j,
                    cost=float(c.cost),
                    R=np.asarray(c.orientation, np.float64),
                    evals=int(rec.last_eval_counts[1]),
                )
            )
        best = min(fins, key=lambda f: f["cost"])
        same_script = (
            np.array_equal(lib.candidates, res.R)
            and lib.strict_costs == [float(c) for c in costs]
            and lib.refined == order
            and lib.finisher_costs == [f["cost"] for f in fins]
            and lib.finisher_evals == [f["evals"] for f in fins]
            and lib.taken == best["cand"]
            and np.array_equal(lib.orientation, best["R"])
            and lib.evals == 5 + sum(f["evals"] for f in fins)
        )
        same_log = (
            lib.strict_costs == e["strict_costs"]
            and lib.refined == e["refined"]
            and lib.finisher_costs == e["finisher_costs"]
            and lib.finisher_evals == e["finisher_evals"]
            and lib.evals == e["evals"]
            and lib.taken == e["taken"]
        )
        rows.append(
            dict(voxel=r, n_att=lib.n_att, evals=lib.evals, same_as_script=bool(same_script),
                 same_as_stored_log=bool(same_log))  # fmt: skip
        )
        print(rows[-1], flush=True)
    summ = dict(
        n=len(rows),
        same_as_script=int(sum(x["same_as_script"] for x in rows)),
        same_as_stored_log=int(sum(x["same_as_stored_log"] for x in rows)),
        rows=rows,
    )
    (OUT / "path_check.json").write_text(json.dumps(summ, indent=1) + "\n")
    print(summ["same_as_script"], "of", summ["n"], "identical to the script path;")
    print(summ["same_as_stored_log"], "of", summ["n"], "identical to the stored log")


if __name__ == "__main__":
    main()
