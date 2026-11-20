"""AR5 step 27, gate t1 (and the old side of c1): paired indexing seconds on the stored seed voxels.

  RESEED_SAMPLE=<mid|held> OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 uv run python \
      scripts/autoresearch/ar5/px_t1.py [--limit N] [--tag T]

The voxels are those of ``ar5/px/voxels.json`` key ``t1`` (frozen): every stored seed step's voxel
of step 22's J5D run. Blocks of 100 voxels: in every block the old setting (``max_pairs``
600,000) and the new (50,000) each index the same voxels once, with ``index_voxel`` (the old
indexer in both), the order alternating from block to block (old first in even blocks). Per block
the seconds, the summed ray pairs and votes of each setting are kept. The preflight of the start
and of the end is KEPT (never required: t1 is judged with any load), with the load average at the
start of every block. No evaluation is made; no truth is read. Writes ``<px>/t1_<sample>.json``.
"""

import argparse
import json
import os
import sys
import time
from pathlib import Path
from typing import Any, Dict, List

import px_common as PC
import px_rules as PR

sys.path.insert(0, str(PC.ICE / "scripts" / "common"))
import preflight as PF  # noqa: E402


def sanitize(x: Any) -> Any:
    """Home paths in the preflight's process commands are replaced (no absolute home path)."""
    s = json.dumps(x)
    s = s.replace(str(Path.home()), "~")
    return json.loads(s)


def verdict(info: Dict[str, Any]) -> str:
    try:
        PF.require_quiet(info=info)
        return "quiet"
    except PF.MachineBusyError as e:
        return f"refused: {e}"


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--tag", default="")
    a = ap.parse_args()
    s = PC.SAMPLE
    start = PF.preflight()
    vox_all = json.loads((PC.RD / "voxels.json").read_text())["t1"][s]["voxels"]
    vox = vox_all[: a.limit] if a.limit else vox_all
    S = PC.setup()
    env, rec = S["env"], S["rec"]
    from icenine import ray_indexing as RI

    import dataclasses

    tot = {k: dict(s=0.0, pairs=0, votes=0) for k in ("old", "new")}
    blocks: List[Dict[str, Any]] = []
    t_start = time.time()
    for b0 in range(0, len(vox), PR.BLOCK):
        blk = vox[b0 : b0 + PR.BLOCK]
        ctx = []
        for v in blk:
            verts = env.hooks.vertices(v)
            ph = int(env.mic.voxels[v].phase)
            sc, crystal, params, _ = rec._index_context(ph)
            ctx.append((sc, verts, crystal, params))
        load0 = os.getloadavg()[0]
        row: Dict[str, Any] = dict(b0=b0, load1_at_start=load0)
        order = ("old", "new") if (b0 // PR.BLOCK) % 2 == 0 else ("new", "old")
        for arm in order:
            mp = PR.OLD_MAX_PAIRS if arm == "old" else PR.NEW_MAX_PAIRS
            pairs = votes = 0
            t0 = time.perf_counter()
            for sc, verts, crystal, params in ctx:
                p = dataclasses.replace(params, max_pairs=mp)
                res, _att = RI.index_voxel(sc, verts, crystal, p, RI.SEED_TOP_K)
                pairs += int(res.n_pairs)
                votes += int(res.n_votes)
            dt = time.perf_counter() - t0
            row[arm] = dict(s=dt, pairs=pairs, votes=votes)
            tot[arm]["s"] += dt
            tot[arm]["pairs"] += pairs
            tot[arm]["votes"] += votes
        blocks.append(row)
        done = b0 + len(blk)
        print(
            f"{s}: {done}/{len(vox)} voxels, old {tot['old']['s']:.0f}s new {tot['new']['s']:.0f}s "
            f"({100 * (1 - tot['new']['s'] / tot['old']['s']):.1f}% less), "
            f"{time.time() - t_start:.0f}s, load {load0:.1f}",
            flush=True,
        )
    end = PF.preflight()
    out: Dict[str, Any] = dict(
        sample=s,
        tag=a.tag,
        n_voxels=len(vox),
        old=tot["old"],
        new=tot["new"],
        less_pct_seconds=100.0 * (1.0 - tot["new"]["s"] / tot["old"]["s"]),
        less_pct_pairs=100.0 * (1.0 - tot["new"]["pairs"] / tot["old"]["pairs"]),
        less_pct_votes=100.0 * (1.0 - tot["new"]["votes"] / tot["old"]["votes"]),
        blocks=blocks,
        verdict_start=verdict(start),
        verdict_end=verdict(end),
        timing_label="contended" if "refused" in verdict(start) + verdict(end) else "single-worker",
        preflight_start=sanitize(start),
        preflight_end=sanitize(end),
        wall_s=time.time() - t_start,
    )
    name = f"t1_{s}{'_' + a.tag if a.tag else ''}{'_limit' if a.limit else ''}.json"
    (PC.RD / name).write_text(json.dumps(out, indent=1) + "\n")
    print(
        json.dumps(
            {
                k: v
                for k, v in out.items()
                if k not in ("blocks", "preflight_start", "preflight_end")
            },
            indent=1,
        )
    )


if __name__ == "__main__":
    main()
