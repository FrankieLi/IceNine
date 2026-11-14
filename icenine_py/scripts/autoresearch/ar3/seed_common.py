"""AR3 step 6 (task ar3f): fallback voxels of the stored gradient-sample runs, replay helpers.

A fallback voxel is a voxel where a stored run made a grid fallback (the seed step returned no
candidate or its end failed the seed test, and the caller ran the full grid search). The runs keep
no per-seed log, so the voxel is found from what they do keep:

  D4, D (BFS):   a reseed in ``run.json`` ``stats.reseeds`` whose ``full_search_evaluations`` is
                 above 100,000 (a grid search costs about 340,000 evaluations, an index seed about
                 2,600; nothing lies between them except the 20 or 63 reseeds that cost 5,000 to
                 200,000, which are index seeds with long local optimizations); a first seed is a
                 record with ``source == "seed"``, ``reseeded`` false and ``n_evals`` above 100,000.
  P, P4 (paint): a row of the paint log (``paint.json.gz`` ``log``) with ``seed == "grid"``.

The checks against the run's own counters (``index_fallback_seed_test``) are in ``counts()``.
"""

import gzip
import json
from pathlib import Path
from typing import Any, Dict, List

import ar3_common as C
import boundary_common as B
import gradient_common as G
import numpy as np

SEED = C.OUT / "seed"
CACHE = SEED / "cache"
WORKTREES = C.ICE.parent.parent  # sibling worktrees hold untracked run files of steps 3 and 4
# Transient input (not in git): step 3's D ``records.json`` in an agent worktree. It is read only
# by the part 0 listing (``seed_part01.py``, and ``seed_extra.py`` through ``counts()``) and is
# gone when that worktree is cleaned up. The summary scripts without ``--slim`` / ``--pack`` read
# committed files only.
D_RECORDS = (
    WORKTREES / "agent-ae68d32fb972a8ad4/icenine_py/benchmarks/autoresearch/ar3/gradient/runs/D"
)
BIG = 100_000  # evaluations: a full search is above this, an index seed is far below
RUNS = ("D4", "D", "P4", "P")
MAP_OF = {"D4": "D4", "D": "D", "P4": "P4R4", "P": "PR"}  # run -> stored map of boundary_common
TEST = 0.80  # the seed test (``MinAccelerationThreshold`` of the run's config)


def _rj(p: Path) -> Dict[str, Any]:
    return json.loads(p.read_text())


def bfs_run(name: str) -> Dict[str, Any]:
    """Records, reseed list and counters of a BFS run (D4 or D)."""
    base = C.OUT / "fine" / "runs" / "D4" if name == "D4" else G.GRAD / "runs" / "D"
    st = _rj(base / "run.json")["stats"]
    if name == "D4":
        rec = G.load_gz(base / "records.json.gz")
    else:
        rec = json.loads((D_RECORDS / "records.json").read_text())
    return dict(stats=st, records=rec)


def paint_run(name: str) -> Dict[str, Any]:
    base = G.GRAD / "runs" / "P" if name == "P" else C.OUT / "fine" / "runs" / "P4"
    return G.load_gz(base / "paint.json.gz")


def fallbacks_bfs(name: str) -> List[Dict[str, Any]]:
    r = bfs_run(name)
    st, rec = r["stats"], r["records"]
    rows: List[Dict[str, Any]] = []
    for x in st["reseeds"]:
        if x["full_search_evaluations"] > BIG:
            rows.append(
                dict(
                    run=name,
                    voxel=int(x["voxel"]),
                    phase="reseed",
                    cluster=int(x["cluster_size"]),
                    full_evals=int(x["full_search_evaluations"]),
                    accepted=bool(x["accepted"]),
                )
            )
    for k, x in rec.items():
        if x["source"] == "seed" and not x["reseeded"] and x["n_evals"] > BIG:
            rows.append(
                dict(
                    run=name,
                    voxel=int(k),
                    phase="seed",
                    cluster=0,
                    full_evals=int(x["n_evals"]),
                    accepted=True,
                )
            )
    return rows


def fallbacks_paint(name: str) -> List[Dict[str, Any]]:
    pj = paint_run(name)
    rows: List[Dict[str, Any]] = []
    for x in pj["log"]:
        if x.get("seed") == "grid":
            rows.append(
                dict(
                    run=name,
                    voxel=int(x["voxel"]),
                    phase="paint",
                    cluster=0,
                    full_evals=0,
                    accepted=True,
                    index_evals=int(x.get("index_evals", 0)),
                    index_hit=x.get("index_hit"),
                )
            )
    return rows


def all_fallbacks() -> List[Dict[str, Any]]:
    rows = fallbacks_bfs("D4") + fallbacks_bfs("D") + fallbacks_paint("P4") + fallbacks_paint("P")
    return rows


def controls_pool(name: str = "D4") -> List[Dict[str, Any]]:
    """The index seeds of a BFS run that were used (the seed step passed the test): every reseed
    below the big-evaluation line and every first seed that is not reseeded and below it."""
    r = bfs_run(name)
    st, rec = r["stats"], r["records"]
    pool = [
        dict(
            voxel=int(x["voxel"]),
            phase="reseed",
            cluster=int(x["cluster_size"]),
            stored_evals=int(x["full_search_evaluations"]),
        )
        for x in st["reseeds"]
        if x["full_search_evaluations"] <= BIG
    ]
    pool += [
        dict(voxel=int(k), phase="seed", cluster=0, stored_evals=int(x["n_evals"]))
        for k, x in rec.items()
        if x["source"] == "seed" and not x["reseeded"] and x["n_evals"] <= BIG
    ]
    return pool


def draw_controls(n: int = 200) -> List[Dict[str, Any]]:
    pool = controls_pool("D4")
    rng = np.random.default_rng(0)
    idx = rng.choice(len(pool), size=n, replace=False)
    return [pool[int(i)] for i in sorted(idx)]


def counts() -> Dict[str, Any]:
    out: Dict[str, Any] = {}
    for name in ("D4", "D"):
        st = bfs_run(name)["stats"]
        rows = fallbacks_bfs(name)
        out[name] = dict(
            listed=len(rows),
            listed_reseed=sum(r["phase"] == "reseed" for r in rows),
            listed_seed=sum(r["phase"] == "seed" for r in rows),
            counter=int(st["index_fallback_seed_test"]),
            counter_no_candidate=int(st.get("index_fallback_no_candidate", 0)),
            grid_evals_counter=int(st["index_fallback_grid_evals"]),
            n_reseeds=int(st["n_reseeds"]),
            n_seeds=int(st["n_seeds"]),
        )
    for name in ("P4", "P"):
        pj = paint_run(name)
        rows = fallbacks_paint(name)
        out[name] = dict(
            listed=len(rows),
            seeds=len(pj["entries"]),
            index_rows=sum(1 for x in pj["log"] if x.get("seed") == "index"),
            seed_grid_evals=int(pj["evals"].get("seed_grid", 0)),
        )
    return out


def err1(env: B.Env, R: np.ndarray, v: int) -> float:
    """Error of one orientation against the truth of voxel ``v`` (analysis only)."""
    return float(env.mis(np.asarray(R, np.float64)[None], env.Rtrue[v][None])[0])


def read_gz(path: Path) -> Any:
    with gzip.open(path, "rt") as f:
        return json.load(f)
