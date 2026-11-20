"""AR5 step 26, phase 0: where the old indexer's seconds go on 200 seed voxels of ``mid``.

  RESEED_SAMPLE=mid IX_PRE_FREEZE=1 uv run python scripts/autoresearch/ar5/ix_profile.py [--label L]

200 voxels: a fixed stride over the stored seed steps of ``mid`` (J5D's ``e1_log``, run order;
first 200 distinct voxels of that stride). Per voxel the old ``index_voxel`` stays untouched; the
seconds are split by ``attach`` (``sc.attached_rays``), ``g`` (scattering vectors), ``vote``,
``cluster``, ``refine`` and ``merge`` (``Candidates.seconds`` of ``index_rays``). One thread.
"""

import argparse
import json
import time
from typing import Any, Dict, List

import ix_common as IC

import numpy as np  # noqa: E402


def pick_voxels(steps: Any, n: int = 200) -> List[int]:
    stride = max(1, len(steps) // n)
    out: List[int] = []
    for v, _ in steps[::stride]:
        if v not in out:
            out.append(v)
        if len(out) == n:
            break
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--label", default="contended")
    ap.add_argument("--n", type=int, default=200)
    a = ap.parse_args()
    s = IC.SAMPLE
    steps = IC.stored_steps(s)
    vox = pick_voxels(steps, a.n)
    S = IC.setup()
    env, rec = S["env"], S["rec"]
    from icenine import ray_indexing as RI

    tot: Dict[str, float] = dict(attach=0.0, g=0.0, vote=0.0, cluster=0.0, refine=0.0, merge=0.0)
    n_votes = n_rays = n_pairs = 0
    t_all = 0.0
    rows = []
    for v in vox:
        verts = env.hooks.vertices(v)
        ph = int(env.mic.voxels[v].phase)
        sc, crystal, params, _ = rec._index_context(ph)
        t0 = time.perf_counter()
        att = sc.attached_rays(verts)
        t1 = time.perf_counter()
        g = RI.scattering_vectors(sc, att)
        t2 = time.perf_counter()
        res = RI.index_rays(g, crystal, params, RI.SEED_TOP_K)
        t3 = time.perf_counter()
        sec = res.seconds
        tot["attach"] += t1 - t0
        tot["g"] += t2 - t1
        for k in ("vote", "cluster", "refine", "merge"):
            tot[k] += sec.get(k, 0.0)
        t_all += t3 - t0
        n_votes += res.n_votes
        n_rays += res.n_rays
        n_pairs += res.n_pairs
        rows.append([v, len(att), res.n_votes, len(res.R), t3 - t0])
    out: Dict[str, Any] = dict(
        sample=s,
        n_voxels=len(vox),
        stride=len(steps) // a.n,
        label=a.label,
        seconds=tot,
        total_s=t_all,
        shares={k: x / t_all for k, x in tot.items()},
        per_voxel_median_s=float(np.median([r[4] for r in rows])),
        mean_rays=n_rays / len(vox),
        mean_votes=n_votes / len(vox),
        mean_pairs=n_pairs / len(vox),
        rows=rows,
    )
    IC.RD.mkdir(parents=True, exist_ok=True)
    (IC.RD / f"profile_{s}_{a.label}.json").write_text(json.dumps(out, indent=1) + "\n")
    print(json.dumps({k: v for k, v in out.items() if k != "rows"}, indent=1))


if __name__ == "__main__":
    main()
