"""AR5 step 16: readers of the stored records of the four disclosed samples (read only).

Arm ``J``'s stored run (step 14: ``e1_log``, ``mc_traces``, ``run.json``, ``polish_log``) and
step 15's recording rows (the stall hit ratios of every 500- and 1000-stall of every seed-step MC
refinement of ``J``) are linked into step 16's cache by ``high500_link_cache.sh`` (per file) and
read here through ``high500_guard.guarded``. Nothing here opens ``held2`` or ``held3``.
"""

import json
import re
from typing import Any, Dict, List

import high500_common as BC  # noqa: E402  (first)

import high500_rules as HR  # noqa: E402
import numpy as np  # noqa: E402
import stallhit_rules as SR  # noqa: E402

TOL = 1e-9


def stored_dir(s: str) -> Any:
    return BC.HCACHE / f"J_{s}"


def load(s: str) -> Dict[str, Any]:
    d = stored_dir(s)
    return dict(
        e1=BC.load_gz(d / "e1_log.json.gz"),
        traces=BC.load_gz(d / "mc_traces.json.gz"),
        run=BC.load_json(d / "run.json"),
        stats=BC.load_json(BC.HD.parent / "polish" / f"stats_{s}.json"),
    )


def counts(s: str, d: Dict[str, Any]) -> Dict[str, Any]:
    """The stored counts against the plan's: ``J`` total with the pass, the MC evaluations of the
    seed steps, and the number of steps and traces."""
    mc = int(sum(int(f[5]) for e in d["e1"] for f in e["fits"] if f[2] == "M"))
    got = dict(
        total=int(d["stats"]["arms"]["J"]["stages"]["total"]),
        mc=mc,
        steps=len(d["e1"]),
        traces=len(d["traces"]),
        mc_fits=int(sum(1 for e in d["e1"] for f in e["fits"] if f[2] == "M")),
    )
    exp = HR.EXPECT[s]
    ok = got["total"] == exp["total"] and got["mc"] == exp["mc"] and got["traces"] == got["mc_fits"]
    return dict(sample=s, expected=exp, got=got, ok=bool(ok))


def same_R(a: Any, b: Any) -> bool:
    x = np.asarray(a, np.float64).ravel()
    y = np.asarray(b, np.float64).ravel()
    return bool(np.max(np.abs(np.round(x, 12) - np.round(y, 12))) <= TOL)


def mc_steps(d: Dict[str, Any]) -> List[Dict[str, Any]]:
    """The stored seed steps with an MC fit: step number (1-based over all steps), the entry and
    the stored traces by candidate ``j``."""
    by: Dict[int, Dict[int, Dict[str, Any]]] = {}
    for t in d["traces"]:
        by.setdefault(int(t["step"]), {})[int(t["j"])] = t
    out = []
    for si, e in enumerate(d["e1"], start=1):
        if any(f[2] == "M" for f in e["fits"]):
            out.append(dict(step=si, voxel=int(e["voxel"]), entry=e, traces=by.get(si, {})))
    return out


def read_records(s: str) -> List[Dict[str, Any]]:
    """The recording rows of sample ``s`` (all shards), sorted by step."""
    rx = re.compile(rf"^record_{s}_s\d+of\d+\.jsonl$")
    rows: List[Dict[str, Any]] = []
    for p in sorted((BC.HCACHE / "record").glob("record_*.jsonl")):
        if rx.match(p.name):
            for line in BC.HG.guarded(p).read_text().splitlines():
                if line.strip():
                    rows.append(json.loads(line))
    rows.sort(key=lambda r: int(r["step"]))
    return rows


def j_taken(e: Dict[str, Any]) -> Dict[str, Any]:
    """Arm ``J``'s taken end of one stored seed step (``stallhit_replay.stored_taken``)."""
    fits = [dict(pos=f[0], j=f[1], kind=f[2], hit=f[3], cost=f[4], evals=f[5]) for f in e["fits"]]
    w = SR.walk_step(fits, {})
    p = w["taken"]
    which = w["kept"][p]["which"]
    row = [f for f in e["fits"] if f[0] == p and f[2] == which][0]
    return dict(
        pos=int(p),
        which=which,
        hit=float(w["kept"][p]["hit"]),
        R=np.asarray(row[6], np.float64),
        stopped=bool(w["stopped"]),
    )


def j_fallbacks(e1: List[Dict[str, Any]], plog: List[Dict[str, Any]]) -> Dict[int, Dict[str, Any]]:
    """Step number -> the polish-log row of its fallback call (a failed seed test), matched in
    order (voxels must agree; a mismatch raises)."""
    att = [r for r in plog if r["attempted"]]
    failed = [
        si
        for si, e in enumerate(e1, start=1)
        if e["fits"] and j_taken(e)["hit"] < SR.SEED_TEST and e["fin_hit"]
    ]
    if len(att) != len(failed):
        raise SystemExit(f"fallback rows {len(att)} != failed seed tests {len(failed)}")
    out = {}
    for si, r in zip(failed, att):
        if int(r["voxel"]) != int(e1[si - 1]["voxel"]):
            raise SystemExit(f"fallback row voxel differs at step {si}")
        out[si] = r
    return out
