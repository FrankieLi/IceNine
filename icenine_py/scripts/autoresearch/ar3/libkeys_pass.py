"""AR3 step 10, pass check (1): the library pass on a stored map gives the stored tries.

  libkeys_pass.py --kind S8L

Reads the map the script read (``before_S8.npy``: the run's ``.mic`` with its rounded angles) and
puts it into the BFS's map in memory, with every voxel FITTED. ``icenine.boundary_pass.
run_boundary_pass`` is applied to it. Exact: its tries (voxel, candidate, end as float32, cost, hit
ratio, confidence, evaluations), its evaluation total and the map after rule ``hit`` equal the
stored ``tries.npz`` and ``after_<map>.npy`` (the script stored the ends as float32, so the
library's float64 ends are compared after the same cast; ``after`` only at the voxels that took a
try, and with the stored float32 values). Writes ``keys/pass_check_<kind>.json``.
"""

import argparse
import json
import time
from pathlib import Path
from typing import Any, Dict

import libkeys_common as K
import numpy as np

import ar3_common as C  # noqa: E402
import gradient_common as G  # noqa: E402


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--kind", required=True, choices=sorted(K.KEYS))
    a = ap.parse_args()
    from icenine import boundary_pass as LBP
    from icenine.mic_file import ReconstructionState

    pdir, before_npy, after_npy = K.STORED_PASS[a.kind]
    out = K.CACHE / f"pass_check_{a.kind}"
    out.mkdir(parents=True, exist_ok=True)
    label = C.preflight_label(out)
    setup, bfs, _t = C.build(C.CONFIG, out, sample=G.GRID_MIC, stack_name=G.STACK)
    mic = setup.sample.get_mic()
    R = np.load(before_npy).astype(np.float64)
    assert len(R) == len(mic.voxels) == K.N
    for v, r in zip(mic.voxels, R):
        v.orientation = r.copy()
        v.reconstruction_id = ReconstructionState.FITTED
    bfs.stats = bfs._new_stats()
    t0 = time.time()
    res = LBP.run_boundary_pass(bfs, mic)
    wall = time.time() - t0
    s = np.load(pdir / "tries.npz")
    n = len(res.tries)
    lib = dict(
        voxel=np.array([t["voxel"] for t in res.tries]),
        cand=np.array([t["cand"] for t in res.tries]),
        end=np.array([t["end"] for t in res.tries], np.float32).reshape(-1, 3, 3),
        cost=np.array([t["cost"] for t in res.tries]),
        hit=np.array([t["hit"] for t in res.tries]),
        conf=np.array([t["conf"] for t in res.tries]),
        evals=np.array([t["evals"] for t in res.tries]),
    )
    info: Dict[str, Any] = dict(
        kind=a.kind,
        n_tries=(n, int(len(s["voxel"]))),
        evals=(int(res.n_evals_tries), int(s["evals"].sum())),
        score_evals_library=int(res.n_evals_score),  # boundary voxels only; the script scored all
        boundary_voxels=(len(res.boundary_voxels), int(len(np.unique(s["voxel"])))),
        wall_s=wall,
        timing_label=label,
    )
    if n == len(s["voxel"]):
        for f in ("voxel", "cand", "evals", "cost", "hit", "conf"):
            bad = np.where(lib[f] != s[f])[0]
            info[f"differ_{f}"] = int(len(bad))
            if len(bad):
                info[f"first_{f}"] = dict(
                    try_index=int(bad[0]), voxel=int(lib["voxel"][bad[0]]),
                    new=float(lib[f][bad[0]]), stored=float(s[f][bad[0]]),
                )  # fmt: skip
        bad_end = np.where(~(lib["end"] == s["end"]).reshape(n, -1).all(1))[0]
        info["differ_end"] = int(len(bad_end))
        # old-end scores of the boundary voxels
        bv = np.array(res.boundary_voxels)
        ours = np.array([res.old_score[v] for v in res.boundary_voxels])
        info["differ_old_score"] = int((ours != s["old_score"][bv]).any(1).sum())
    after = np.load(after_npy)
    got = K.orientations(mic)
    changed = np.where((got != R).reshape(len(R), -1).any(1))[0]
    stored_changed = np.where((after != R).reshape(len(R), -1).any(1))[0]
    info["voxels_changed"] = (int(len(changed)), int(len(stored_changed)))
    same_set = set(changed.tolist()) == set(stored_changed.tolist())
    info["same_changed_voxels"] = bool(same_set)
    if same_set:
        # the script's map holds the float32-rounded end made orthonormal (boundary_lib.Tries:
        # SVD of the float32 end); the library keeps its float64 end. Compare after the same steps.
        def orth(m: np.ndarray) -> np.ndarray:
            u, _s, vt = np.linalg.svd(m)
            return u @ vt

        stored_like = orth(got[changed].astype(np.float32).astype(np.float64))
        diff = [int(v) for v, m in zip(changed, stored_like) if not (m == after[v]).all()]
        info["changed_voxels_with_a_different_end_after_the_script_rounding"] = len(diff)
        info["first_differing_voxel"] = diff[:5]
        info["max_abs_diff_library_vs_stored_map_at_changed_voxels"] = float(
            np.abs(got[changed] - after[changed]).max() if len(changed) else 0.0
        )
        unchanged_equal = (
            np.delete(got, changed, axis=0) == np.delete(after, changed, axis=0)
        ).all()
        info["unchanged_voxels_equal"] = bool(unchanged_equal)
    (K.OUT / f"pass_check_{a.kind}.json").write_text(json.dumps(info, indent=1) + "\n")
    print(json.dumps(info, indent=1))


if __name__ == "__main__":
    main()
