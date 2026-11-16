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
    ap.add_argument("--out", default="")
    a = ap.parse_args()
    md: List[str] = []
    summ: Dict[str, Any] = {}
    for name in a.runs:
        d = CM.OUT / "runs" / name
        rows = CM.G.load_gz(d / "voxels.json.gz")
        run = json.loads((d / "run.json").read_text())
        keys = [k[4:] for k in rows[0] if k.startswith("err_") and k != "err_C"] + ["C"]
        md += [f"## {name}", "", f"args: `{json.dumps(run['args'])}`", ""]
        md += table(rows, keys) + [""]
        G = run["grains"]
        if "stages" in G[0]:  # run_pipeline.py
            nst = len(G[0]["stages"])
            acc = [sum(1 for g in G if g["stages"][i]["accepted"]) for i in range(nst)]
            nv = np.array([g["n_vox"] for g in G], float)
            rpv = float((np.array([g["renders_per_voxel"] for g in G]) * nv).sum() / nv.sum())
            tsec = float(sum(s["t_s"] for g in G for s in g["stages"]))
            md += [f"grains: {len(G)}; stages {run['stages']}; accepted per stage {acc}; "
                   f"renders per voxel {rpv:.1f}; stage wall time {tsec:.1f} s "
                   f"({tsec / nv.sum() * 1000:.1f} ms per voxel); timing: {run.get('timing', '')}", ""]
            last = f"s{nst - 1}"
            withC = [r for r in rows if "err_C" in r]
            e = np.array([r[f"err_{last}"] for r in withC])
            summ[name] = dict(stages=run["stages"], renders_per_voxel=rpv, wall_s=tsec,
                              ms_per_voxel=tsec / nv.sum() * 1000, accepted=acc, n=len(e),
                              k_gt01=int((e > 0.1).sum()), median=float(np.median(e)),
                              p95=float(np.percentile(e, 95)), k_wrong=int((e > 1).sum()),
                              n_all=len(rows),
                              k_gt01_all=int((np.array([r[f"err_{last}"] for r in rows]) > 0.1).sum()))
        else:
            acc = {s: sum(1 for g in G if g.get(s, {}).get("accepted")) for s in ("A", "B")}
            md += [f"grains: {len(G)}; accepted A {acc['A']}, B {acc['B']}", ""]
            summ[name] = dict(accepted=acc)
    ap_out = CM.OUT / (a.out or "tables")
    ap_out.parent.mkdir(parents=True, exist_ok=True)
    (ap_out.parent / f"{ap_out.name}.md").write_text("\n".join(md) + "\n")
    (ap_out.parent / f"{ap_out.name}.json").write_text(json.dumps(summ, indent=1) + "\n")
    for k, v in summ.items():
        if "renders_per_voxel" in v:
            print(f"{k:16s} r/vox {v['renders_per_voxel']:6.1f} ms/vox {v['ms_per_voxel']:6.1f} "
                  f"C-voxels n {v['n']} gt0.1 {v['k_gt01']} med {v['median']:.4f} p95 {v['p95']:.4f} "
                  f"wrong {v['k_wrong']} | all n {v['n_all']} gt0.1 {v['k_gt01_all']} acc {v['accepted']}")


if __name__ == "__main__":
    main()
