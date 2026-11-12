"""AR5 step 5 (fixes round): the contents of the boundary pass inside the library runs.

  KEYS5_CACHE=<folder of the library runs> uv run python \\
      scripts/autoresearch/ar5/keys_pass_contents.py

Reads stored files only; runs nothing. For each library run (``<KEYS5_CACHE>/<kind>_<sample>/
pass_tries.npz``, written by ``keys_run.py``) it compares the tries of the pass inside the run
with the stored script pass (``tries.npz`` of ``keys_common.STORED_PASS``), try by try in order:
the voxel and candidate, the cost, the hit ratio, the evaluations and the end. The angle between
two ends is the rotation angle of ``A^T B``, with no symmetry reduction, computed as
``atan2(|skew part|, (trace - 1) / 2)`` so that small angles are accurate. It also compares the
sets of voxels the two passes changed: the library set is the run's ``changed`` list; the stored set
is the voxels whose orientation differs between ``before_<name>.npy`` and ``after_<name>.npy``.
Both sets are checked against a second count (the run's ``boundary_changed`` counter and the
stored pass rebuilt with ``icenine.boundary_pass.keep_hit`` from the stored scores and tries).

Writes ``benchmarks/autoresearch/ar5/keys/pass_contents.json`` (paths in it are relative to the
worktrees folder) and ``tables.md`` (block ``pass_contents``, copied into the findings page by
``scripts/dev/sync_doc_tables.py``). ``--tables-only`` rewrites ``tables.md`` from the JSON.
"""

import argparse
import json
import sys
from pathlib import Path
from typing import Any, Dict, List, Set

import keys_common as K
import numpy as np

from icenine.boundary_pass import keep_hit

RUNS = [("S8L", "big"), ("S10L", "big"), ("S11L", "big"), ("S11L", "grad")]


def rel(p: Path) -> str:
    try:
        return str(p.resolve().relative_to(K.WORKTREES.resolve()))
    except ValueError:
        return p.name


def angle_deg(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    """Rotation angle in degrees of ``a[i]^T b[i]`` (no symmetry reduction)."""
    m = np.einsum("nji,njk->nik", a, b)
    tr = np.einsum("nii->n", m)
    s = np.stack([m[:, 2, 1] - m[:, 1, 2], m[:, 0, 2] - m[:, 2, 0], m[:, 1, 0] - m[:, 0, 1]], 1)
    return np.degrees(np.arctan2(np.linalg.norm(s, axis=1) / 2.0, (tr - 1.0) / 2.0))


def stored_changed_by_rule(b: Any) -> Set[int]:
    """The stored pass rebuilt from its scores and tries with rule ``hit``."""
    out: Set[int] = set()
    voxel = b["voxel"]
    old = b["old_score"]
    starts = np.flatnonzero(np.r_[True, voxel[1:] != voxel[:-1]])
    ends = np.r_[starts[1:], len(voxel)]
    for s, e in zip(starts, ends):
        v = int(voxel[s])
        tries = [(float(b["cost"][i]), float(b["hit"][i])) for i in range(s, e)]
        if keep_hit((float(old[v, 0]), float(old[v, 1])), tries) > 0:
            out.add(v)
    return out


def compare(kind: str, sample: str) -> Dict[str, Any]:
    lib_dir = K.CACHE / f"{kind}_{sample}"
    pdir, before_npy, after_npy = K.STORED_PASS[(kind, sample)]
    a, b = np.load(lib_dir / "pass_tries.npz"), np.load(pdir / "tries.npz")
    run = json.loads((lib_dir / "run.json").read_text())
    res: Dict[str, Any] = dict(
        library=rel(lib_dir / "pass_tries.npz"),
        stored=rel(pdir / "tries.npz"),
        tries=[int(len(a["voxel"])), int(len(b["voxel"]))],
        evals=[int(a["evals"].sum()), int(b["evals"].sum())],
    )
    lib_changed = {int(v) for v in a["changed"]}
    before = np.load(before_npy).astype(np.float64)
    after = np.load(after_npy).astype(np.float64)
    st_changed = {int(v) for v in np.flatnonzero((before != after).reshape(len(before), -1).any(1))}
    by_rule = stored_changed_by_rule(b)
    res["changed"] = dict(
        library=len(lib_changed),
        library_counter=int(run["stats"]["boundary_changed"]),
        stored=len(st_changed),
        stored_rebuilt_by_rule=len(by_rule),
        stored_rule_equals_map=by_rule == st_changed,
        both=len(lib_changed & st_changed),
        only_library=len(lib_changed - st_changed),
        only_stored=len(st_changed - lib_changed),
    )
    if len(a["voxel"]) != len(b["voxel"]):
        res["same_length"] = False
        return res
    n = len(a["voxel"])
    other_cand = np.flatnonzero((a["voxel"] != b["voxel"]) | (a["cand"] != b["cand"]))
    ang = angle_deg(a["end"].astype(np.float64), b["end"].astype(np.float64))
    entry = np.abs(a["end"].astype(np.float64) - b["end"].astype(np.float64)).reshape(n, -1)
    worst = int(np.argmax(ang))
    rows: List[Dict[str, int]] = [
        dict(
            try_index=int(i),
            voxel=int(a["voxel"][i]),
            cand_library=int(a["cand"][i]),
            cand_stored=int(b["cand"][i]),
        )
        for i in other_cand
    ]
    res.update(
        same_length=True,
        tries_with_another_candidate=len(other_cand),
        another_candidate=rows,
        cost_differs=int((a["cost"] != b["cost"]).sum()),
        hit_ratio_differs=int((a["hit"] != b["hit"]).sum()),
        evals_differ=int((a["evals"] != b["evals"]).sum()),
        end_exactly_equal=int((a["end"] == b["end"]).reshape(n, -1).all(1).sum()),
        end_angle_max_deg=float(ang.max()),
        end_angle_max_at=dict(try_index=worst, voxel=int(a["voxel"][worst])),
        end_angle_over_0_01_deg=int((ang > 0.01).sum()),
        end_angle_over_0_1_deg=int((ang > 0.1).sum()),
        end_angle_median_deg=float(np.median(ang)),
        end_entry_max_abs_difference=float(entry.max()),
    )
    return res


def pair(a: int, b: int) -> str:
    return f"{a:,} / {b:,}"


def write_doc_tables(data: Dict[str, Any]) -> None:
    """Write ``tables.md`` (block ``pass_contents``) from the JSON, for the findings page."""
    sys.path.insert(0, str(K.ICE / "scripts" / "common"))
    from doc_tables import markdown_table, write_tables

    sample = {"big": "big", "grad": "first gradient"}
    rows = []
    for name, r in data["runs"].items():
        kind, s = name.split("_")
        c = r["changed"]
        row: Dict[str, Any] = {
            "Run": kind,
            "Sample": sample[s],
            "Tries (library / stored)": pair(*r["tries"]),
            "Changed voxels (library / stored)": pair(c["library"], c["stored"]),
            "Changed in one run only (library / stored)": pair(c["only_library"], c["only_stored"]),
        }
        if r["same_length"]:
            row.update(
                {
                    "Tries with another candidate": r["tries_with_another_candidate"],
                    "Costs differ": r["cost_differs"],
                    "Largest angle between ends (deg)": r["end_angle_max_deg"],
                    "Ends over 0.01 deg apart": r["end_angle_over_0_01_deg"],
                }
            )
        rows.append(row)
    cols = [
        "Run",
        "Sample",
        "Tries (library / stored)",
        "Tries with another candidate",
        "Costs differ",
        "Largest angle between ends (deg)",
        "Ends over 0.01 deg apart",
        "Changed voxels (library / stored)",
        "Changed in one run only (library / stored)",
    ]
    fmt = {
        "Tries with another candidate": ",d",
        "Costs differ": ",d",
        "Ends over 0.01 deg apart": ",d",
        "Largest angle between ends (deg)": ".2f",
    }
    write_tables(K.OUT / "tables.md", {"pass_contents": markdown_table(rows, cols, fmt)})


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--tables-only", action="store_true", help="rewrite tables.md from pass_contents.json"
    )
    path = K.OUT / "pass_contents.json"
    if ap.parse_args().tables_only:
        write_doc_tables(json.loads(path.read_text()))
        return
    out: Dict[str, Any] = dict(
        note=(
            "Pass inside the library run against the stored script pass, try by try. Read from "
            "stored files only. Angles: rotation angle of A^T B, no symmetry reduction."
        ),
        runs={f"{k}_{s}": compare(k, s) for k, s in RUNS},
    )
    K.OUT.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(out, indent=1) + "\n")
    write_doc_tables(out)
    for name, r in out["runs"].items():
        short = {k: v for k, v in r.items() if k not in ("another_candidate", "library", "stored")}
        print(name, json.dumps(short))


if __name__ == "__main__":
    main()
