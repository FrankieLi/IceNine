"""AR5 step 9, phase B: the library pass on the stored map of ``Gl`` gives the script pass's tries?

  RESEED_SAMPLE=mid uv run python scripts/autoresearch/ar5/rescue_lib_pass.py

As ``keys_pass.py`` of AR5 step 5, for ``mid``: reads the map the script pass of ``Gl`` read
(``passes_mid/before_Gl.npy``: the stored ``.mic`` with its rounded angles; the library run
``LibG`` stored the same bytes before its pass) into the BFS's map in memory, every voxel FITTED,
and applies ``icenine.boundary_pass.run_boundary_pass`` (the function the library run calls).
Compares its tries (voxel, candidate, cost, hit ratio, evaluations, end) and evaluation total with
``passes_mid/D4_Gl/tries.npz`` and the changed voxels with ``after_Gl.npy``. The pass inside the
library run started from the unrounded map in memory; this check says whether the library pass
and the script pass agree when they start from the same map. Writes ``rescue/lib_pass.json``.
"""

import json
import os
import sys
import time
from pathlib import Path
from typing import Any, Dict

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "ar3"))

import reseed_common as RC  # noqa: E402  (first: RESEED_SAMPLE -> AR3_SAMPLE)

import ar3_common as C  # noqa: E402
import gradient_common as G  # noqa: E402
import numpy as np  # noqa: E402

RES = C.ICE / "benchmarks" / "autoresearch" / "ar5" / "rescue"
CACHE = Path(os.environ.get("RESCUE_CACHE", str(RES / "cache")))


def main() -> None:
    from icenine import boundary_pass as LBP
    from icenine.mic_file import ReconstructionState

    assert RC.SAMPLE == "mid"
    pdir = CACHE / "passes_mid"
    out = CACHE / "LibPass_mid"
    if out.exists():
        raise SystemExit(f"refusing to start: {out} exists")
    out.mkdir(parents=True)
    label = C.preflight_label(out)
    setup, bfs, _t = C.build(C.CONFIG, out, sample=G.GRID_MIC, stack_name=G.STACK)
    mic = setup.sample.get_mic()
    R = np.load(pdir / "before_Gl.npy").astype(np.float64)
    assert len(R) == len(mic.voxels)
    for v, r in zip(mic.voxels, R):
        v.orientation = r.copy()
        v.reconstruction_id = ReconstructionState.FITTED
    bfs.stats = bfs._new_stats()
    t0 = time.time()
    res = LBP.run_boundary_pass(bfs, mic)
    wall = time.time() - t0
    s = np.load(pdir / "D4_Gl" / "tries.npz")
    n = len(res.tries)
    lib = dict(
        voxel=np.array([t["voxel"] for t in res.tries]),
        cand=np.array([t["cand"] for t in res.tries]),
        end=np.array([t["end"] for t in res.tries], np.float32).reshape(-1, 3, 3),
        cost=np.array([t["cost"] for t in res.tries]),
        hit=np.array([t["hit"] for t in res.tries]),
        evals=np.array([t["evals"] for t in res.tries]),
    )
    info: Dict[str, Any] = dict(
        sample="mid",
        n_tries=dict(library=n, script=int(len(s["voxel"]))),
        evals=dict(library=int(res.n_evals_tries), script=int(s["evals"].sum())),
        wall_s=wall,
        timing_label=label,
        wall_label="contended; not a timing result",
    )
    if n == len(s["voxel"]):
        for f in ("voxel", "cand", "evals", "cost", "hit"):
            info[f"differ_{f}"] = int((lib[f] != s[f]).sum())
        info["differ_end"] = int((~(lib["end"] == s["end"]).reshape(n, -1).all(1)).sum())
    after = np.load(pdir / "after_Gl.npy")
    got = np.stack([np.asarray(v.orientation, np.float64) for v in mic.voxels])
    changed = set(np.where((got != R).reshape(len(R), -1).any(1))[0].tolist())
    stored = set(np.where((after != R).reshape(len(R), -1).any(1))[0].tolist())
    info["voxels_changed"] = dict(library=len(changed), script=len(stored))
    info["same_changed_voxels"] = changed == stored
    info["map_after_equal"] = bool(np.array_equal(got, after))
    info["exact"] = bool(
        info["n_tries"]["library"] == info["n_tries"]["script"]
        and info["evals"]["library"] == info["evals"]["script"]
        and all(info.get(f"differ_{f}", 1) == 0 for f in ("voxel", "cand", "evals", "cost", "hit"))
        and info["same_changed_voxels"]
    )
    (RES / "lib_pass.json").write_text(json.dumps(info, indent=1) + "\n")
    print(json.dumps(info, indent=1), flush=True)
    os._exit(0)


if __name__ == "__main__":
    main()
