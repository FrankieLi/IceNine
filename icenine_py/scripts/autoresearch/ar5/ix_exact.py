"""AR5 step 26, gates e1 and w1: the old and the new indexer on the exactness voxels of a sample.

  RESEED_SAMPLE=<sample> OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 uv run python \
      scripts/autoresearch/ar5/ix_exact.py [--limit N] [--resume] [--require-quiet]

The voxels are those of ``ar5/ix/voxels.json`` (frozen). Blocks of 100 voxels: in every block the
old ``index_voxel`` and the new ``index_voxel_fast`` each run once on the same voxels, the order
alternating from block to block (old first in even blocks). Per voxel the outputs are compared
bit for bit (tier 1) and, when they differ, by the tier 2 rule (same count, order, attached rays,
rotation entries within 1e-12). Per block the seconds of each arm (attaching the rays, voting,
clustering, refining and merging; one thread) are summed. Writes ``<ix>/exact_<sample>.json``
(summary, per-block seconds, the voxels that differ) and the preflight of the start and the end
(``--require-quiet``: refuse to start unless the preflight passes; a run without it is labelled
"contended" and is outside gate w1). No evaluation is made; no truth is read.
"""

import argparse
import json
import sys
import time
from pathlib import Path
from typing import Any, Dict, List

import ix_common as IC
import ix_rules as IR

import numpy as np  # noqa: E402

sys.path.insert(0, str(IC.ICE / "scripts" / "common"))
import preflight as PF  # noqa: E402


def sanitize(x: Any) -> Any:
    """Home paths in the preflight's process commands are replaced (no absolute home path)."""
    s = json.dumps(x)
    s = s.replace(str(Path.home()), "~")
    return json.loads(s)


def tier1(a: Any, b: Any) -> bool:
    return bool(
        a.R.shape == b.R.shape
        and np.array_equal(a.R, b.R)
        and np.array_equal(a.explained, b.explained)
        and np.array_equal(a.votes, b.votes)
        and (a.n_rays, a.n_entries, a.n_pairs, a.n_votes)
        == (b.n_rays, b.n_entries, b.n_pairs, b.n_votes)
    )


def tier2(a: Any, b: Any, att_a: Any, att_b: Any) -> bool:
    return bool(
        a.R.shape == b.R.shape
        and np.array_equal(a.explained, b.explained)
        and np.array_equal(a.votes, b.votes)
        and np.array_equal(att_a, att_b)
        and (a.n_rays, a.n_entries, a.n_pairs, a.n_votes)
        == (b.n_rays, b.n_entries, b.n_pairs, b.n_votes)
        and (a.R.size == 0 or float(np.abs(a.R - b.R).max()) <= 1e-12)
    )


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--require-quiet", action="store_true")
    ap.add_argument("--tag", default="")
    a = ap.parse_args()
    s = IC.SAMPLE
    start: Dict[str, Any] = PF.preflight()
    quiet_start = True
    try:
        PF.require_quiet(info=start)
    except PF.MachineBusyError as e:
        quiet_start = False
        if a.require_quiet:
            print("refusing to start:", e, file=sys.stderr)
            raise SystemExit(2)
    vox_all = json.loads((IC.RD / "voxels.json").read_text())["samples"][s]["voxels"]
    vox = vox_all[: a.limit] if a.limit else vox_all
    S = IC.setup()
    env, rec = S["env"], S["rec"]
    from icenine import ray_indexing as RI

    t_old = t_new = 0.0
    sec_old: Dict[str, float] = {}
    sec_new: Dict[str, float] = {}
    blocks: List[List[float]] = []
    differ: List[Dict[str, Any]] = []
    n_t1 = n_t2 = 0
    t_start = time.time()
    for b0 in range(0, len(vox), IR.BLOCK):
        blk = vox[b0 : b0 + IR.BLOCK]
        ctx = []
        for v in blk:
            verts = env.hooks.vertices(v)
            ph = int(env.mic.voxels[v].phase)
            sc, crystal, params, _ = rec._index_context(ph)
            ctx.append((v, sc, verts, crystal, params))
        res: Dict[str, List[Any]] = {"old": [], "new": []}
        secs = {"old": 0.0, "new": 0.0}
        order = ("old", "new") if (b0 // IR.BLOCK) % 2 == 0 else ("new", "old")
        for arm in order:
            fn = RI.index_voxel if arm == "old" else RI.index_voxel_fast
            t0 = time.perf_counter()
            for v, sc, verts, crystal, params in ctx:
                res[arm].append(fn(sc, verts, crystal, params, RI.SEED_TOP_K))
            secs[arm] = time.perf_counter() - t0
        t_old += secs["old"]
        t_new += secs["new"]
        blocks.append([b0, secs["old"], secs["new"]])
        for i, (v, *_rest) in enumerate(ctx):
            (ra, att_a), (rb, att_b) = res["old"][i], res["new"][i]
            for k, x in ra.seconds.items():
                sec_old[k] = sec_old.get(k, 0.0) + x
            for k, x in rb.seconds.items():
                sec_new[k] = sec_new.get(k, 0.0) + x
            ok1 = tier1(ra, rb) and bool(np.array_equal(att_a, att_b))
            ok2 = ok1 or tier2(ra, rb, att_a, att_b)
            n_t1 += int(ok1)
            n_t2 += int(ok2)
            if not ok1:
                differ.append(
                    dict(
                        voxel=int(v),
                        tier2_ok=bool(ok2),
                        n_old=int(len(ra.R)),
                        n_new=int(len(rb.R)),
                        max_dR=(
                            float(np.abs(ra.R - rb.R).max())
                            if ra.R.shape == rb.R.shape and ra.R.size
                            else None
                        ),
                    )
                )
        done = b0 + len(blk)
        print(
            f"{s}: {done}/{len(vox)} voxels, tier1 {n_t1}, differ {len(differ)}, old {t_old:.0f}s "
            f"new {t_new:.0f}s ({100 * (1 - t_new / t_old):.1f}% less), "
            f"{time.time() - t_start:.0f}s",
            flush=True,
        )
    end = PF.preflight()
    quiet_end = True
    try:
        PF.require_quiet(info=end)
    except PF.MachineBusyError:
        quiet_end = False
    label = "single-worker" if (quiet_start and quiet_end) else "contended"
    out: Dict[str, Any] = dict(
        sample=s,
        tag=a.tag,
        n_voxels=len(vox),
        tier1_identical=n_t1,
        tier2_ok=n_t2,
        n_differ=len(differ),
        differ=differ,
        old_s=t_old,
        new_s=t_new,
        less_pct=100.0 * (1.0 - t_new / t_old) if t_old else None,
        old_parts_s=sec_old,
        new_parts_s=sec_new,
        blocks=blocks,
        timing_label=label,
        quiet_start=quiet_start,
        quiet_end=quiet_end,
        preflight_start=sanitize(start),
        preflight_end=sanitize(end),
        wall_s=time.time() - t_start,
    )
    name = f"exact_{s}{'_' + a.tag if a.tag else ''}{'_limit' if a.limit else ''}.json"
    (IC.RD / name).write_text(json.dumps(out, indent=1) + "\n")
    print(
        json.dumps(
            {
                k: v
                for k, v in out.items()
                if k not in ("differ", "blocks", "preflight_start", "preflight_end")
            },
            indent=1,
        )
    )


if __name__ == "__main__":
    main()
