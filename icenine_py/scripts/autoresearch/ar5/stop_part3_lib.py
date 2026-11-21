"""AR5 step 1, part 3: the map of S8 after the boundary pass, with hit ratio, error and source.

The hit ratio of a voxel after the pass comes from the stored scores: its old end's, or the end of
the try the pass took (``stop_rules.choose_try`` with every try kept; checked against the stored
map after the pass by ``stop_part1_starts.py``). The source of a voxel is that of its record in S8
(``records.json.gz``): seed, neighbor, revisit, or reseed (a seed that was reseeded). A boundary
voxel is a voxel with at least one try. The error is for scoring only.
"""

import json
from pathlib import Path
from typing import Any, Dict, List

import numpy as np
import stop_part1_lib as P
import stop_rules as SR

SOURCES = ("seed", "neighbor", "revisit", "reseed")
STORED = P.AR5 / "cache" / "stored"


def load(sample: str) -> Dict[str, Any]:
    import gzip

    S = P.load(sample)
    voxel = S["voxel"]
    old_hit = S["old_hit"].astype(np.float64)
    old_cost = S["old_cost"].astype(np.float64)
    ch = SR.choose_try(
        voxel,
        old_hit,
        old_cost,
        S["end_hit"].astype(np.float64),
        S["end_cost"].astype(np.float64),
        np.ones(len(voxel), bool),
    )
    hit = SR.hit_after(old_hit, voxel, S["end_hit"].astype(np.float64), ch)
    n = len(old_hit)
    boundary = np.zeros(n, bool)
    boundary[voxel] = True
    with gzip.open(STORED / sample / f"S8_{sample}" / "records.json.gz", "rt") as f:
        rec = json.load(f)
    src = np.array(["none"] * n, dtype=object)
    for k, r in rec.items():
        s = r["source"]
        src[int(k)] = "reseed" if (s == "seed" and r["reseeded"]) else s
    err_after = S["err_after"].astype(np.float64)
    return dict(
        hit=hit,
        wrong=err_after > P.WRONG,
        err=err_after,
        boundary=boundary,
        source=src,
        n=n,
        old_hit=old_hit,
        choice=ch,
    )


def table(M: Dict[str, Any], thresholds: List[float] = list(SR.THRESHOLDS)) -> Dict[str, Any]:
    """Per threshold: all voxels, by source, boundary voxels alone."""
    out: Dict[str, Any] = {}
    for t in thresholds:
        row = dict(all=SR.flag_counts(M["hit"], M["wrong"], t))
        for s in SOURCES:
            m = M["source"] == s
            row[s] = SR.flag_counts(M["hit"][m], M["wrong"][m], t)
        b = M["boundary"]
        row["boundary"] = SR.flag_counts(M["hit"][b], M["wrong"][b], t)
        out[f"{t:.2f}"] = row
    return out
