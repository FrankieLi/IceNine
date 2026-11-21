"""AR5 step 1, part 3c: the vote count of the unexplained rays (H3) at flagged voxels.

  RESEED_SAMPLE=big|grad|mid uv run python scripts/autoresearch/ar5/stop_part3_votes.py [--limit N]

For every voxel flagged at the picked threshold (hit ratio of the map after the pass under it) and
for every 20th unflagged boundary voxel by index: take the rays attached to the voxel
(``SmoothCost.attached_rays``), drop those the voxel's orientation after the pass explains
(``explained_mask``), index the rest (``index_rays``) and keep the highest vote count of the
candidates (0 when there are none). No library function does this; the steps are the library's.
Writes ``cache/votes_<sample>.json``; ``stop_part3_summary.py`` writes the slim result. Truth is
used only to label a voxel wrong or right (error above 1 degree).
"""

import argparse
import json
import time
from typing import Any, Dict, List

import numpy as np
import stop_common as SC
import stop_part3_lib as L
import stop_rules as SR

import reseed_common as RC  # noqa: E402


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--limit", type=int, default=0)
    a = ap.parse_args()
    thr = json.loads((SC.AR5 / "part3_pick.json").read_text())["threshold"]
    M = L.load(SC.SAMPLE)
    flagged = M["hit"] < thr - SR.TOL
    sel_un = SR.sample_unflagged(flagged, M["boundary"], every=20)
    vox = np.concatenate([np.where(flagged)[0], sel_un])
    if a.limit:
        vox = np.concatenate([np.where(flagged)[0][: a.limit], sel_un[: a.limit]])
    print(SC.SAMPLE, "flagged", int(flagged.sum()), "unflagged sampled", len(sel_un), flush=True)
    from icenine import ray_indexing as RI

    env = RC.make_env()
    rec = env.bfs.reconstructor
    after = np.load(SC.PASS_DIR / "after_S8.npy")
    rows: List[Dict[str, Any]] = []
    t0 = time.time()
    for i, v in enumerate(vox):
        v = int(v)
        vertices = env.hooks.vertices(v)
        phase = env.mic.voxels[v].phase
        sc, crystal, params, _ = rec._index_context(phase)
        att = sc.attached_rays(vertices)
        g = RI.scattering_vectors(sc, att)
        ex = RI.explained_mask(after[v], g, crystal, params.ang_tol_deg, params.len_tol)
        gu = g[~ex]
        res = RI.index_rays(gu, crystal, params, 20)
        count = int(np.max(res.votes)) if len(res.votes) else 0
        rows.append(
            dict(
                voxel=v,
                flagged=bool(flagged[v]),
                wrong=bool(M["wrong"][v]),
                n_attached=int(len(att)),
                n_unexplained=int(len(gu)),
                max_votes=count,
                n_cand=int(len(res.votes)),
            )
        )
        if i % 200 == 0:
            print(i, len(vox), f"{time.time() - t0:.0f} s", flush=True)
    name = f"votes_{SC.SAMPLE}" + (f"_limit{a.limit}" if a.limit else "")
    (SC.AR5 / "cache" / f"{name}.json").write_text(
        json.dumps(
            dict(
                sample=SC.SAMPLE,
                threshold=thr,
                rows=rows,
                wall_s=time.time() - t0,
                timing_label="contended (another agent on the machine); not a timing result",
            )
        )
        + "\n"
    )
    print("done", len(rows), f"{time.time() - t0:.0f} s", flush=True)


if __name__ == "__main__":
    main()
