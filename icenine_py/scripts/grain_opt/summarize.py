"""Tables of grain-fit runs: per-voxel error against truth for the start, stage A, stage B and the
AR3 step 4 CMA per-voxel fit (``err_C``, only voxels that run fitted: painted, non-seed voxels).

Rows: all voxels of the fitted groups, the voxels with a C record (same-voxel comparison), interior
and boundary voxels, and per gradient level. Columns: n, median, 95th percentile, above 0.1 deg,
wrong (above 1 deg), with Wilson 95% intervals; mean strict cost.

Usage (from icenine_py/): uv run python scripts/grain_opt/summarize.py A_sink AB_sink ...
"""

import argparse
import json
from typing import Any, Dict, List

import common as CM  # noqa: I001
import numpy as np


def table(rows: List[Dict[str, Any]], keys: List[str]) -> List[str]:
    out = ["| group | stage | n | median | p95 | above 0.1 | wrong | strict |", "|---|---|---|---|---|---|---|---|"]

    def add(name: str, sel: List[Dict[str, Any]]) -> None:
        for k in keys:
            vals = [r for r in sel if f"err_{k}" in r]
            if not vals:
                continue
            e = np.array([r[f"err_{k}"] for r in vals])
            st = CM.err_stats(e)
            sc = np.array([r.get(f"strict_{k}", np.nan) for r in vals])
            out.append(
                f"| {name} | {k} | {st['n']} | {st['median']:.4f} | {st['p95']:.4f} | {st['gt01']} "
                f"| {st['wrong']} | {np.nanmean(sc):.4f} |"
            )

    add("all", rows)
    withC = [r for r in rows if "err_C" in r]
    add("voxels with a C fit", withC)
    add("interior (C voxels)", [r for r in withC if not r["boundary"]])
    add("boundary (C voxels)", [r for r in withC if r["boundary"]])
    for lv in sorted({r["level"] for r in rows}):
        add(f"level {lv} (C voxels)", [r for r in withC if r["level"] == lv])
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("runs", nargs="+")
    a = ap.parse_args()
    md: List[str] = []
    summ: Dict[str, Any] = {}
    for name in a.runs:
        d = CM.OUT / "runs" / name
        rows = CM.G.load_gz(d / "voxels.json.gz")
        run = json.loads((d / "run.json").read_text())
        keys = ["start"] + [k for k in ("A_raw", "A", "B_raw", "B") if f"err_{k}" in rows[0]] + ["C"]
        md += [f"## {name}", "", f"args: `{json.dumps(run['args'])}`", ""]
        md += table(rows, keys) + [""]
        acc = {s: sum(1 for g in run["grains"] if g.get(s, {}).get("accepted")) for s in ("A", "B")}
        ts = {s: sum(g.get(s, {}).get("t_s", 0.0) for g in run["grains"]) for s in ("A", "B")}
        ev = {s: sum(g.get(s, {}).get("evals", 0) for g in run["grains"]) for s in ("A", "B")}
        md += [
            f"grains: {len(run['grains'])}; accepted A {acc['A']}, B {acc['B']}; "
            f"time A {ts['A']:.0f} s, B {ts['B']:.0f} s; objective evaluations A {ev['A']}, B {ev['B']}; "
            f"timing: {run.get('timing', '')}",
            "",
        ]
        summ[name] = dict(accepted=acc, t_s=ts, evals=ev)
    (CM.OUT / "tables.md").write_text("\n".join(md) + "\n")
    print("\n".join(md))


if __name__ == "__main__":
    main()
