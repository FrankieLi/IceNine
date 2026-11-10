"""AR1 stress test: recall when the voxel has fewer rays (grains with few spots).

Keeps every m-th attached ray (a fixed stride, no random draw) and indexes the rest, on the
pilot's full-search voxels. The true grain then explains about 1/m of its usual rays. Also a
parameter check at m = 1 (tolerances half and double of the defaults) so the tolerance choice is
visible.

  uv run python scripts/autoresearch/ar1/sparse.py
Writes benchmarks/autoresearch/ar1/sparse_pilot.json.
"""

import json
import sys
from dataclasses import replace
from pathlib import Path
from typing import Any, Dict, List

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import ar1_common as AC  # noqa: E402

from icenine import ray_indexing as RI  # noqa: E402

STRIDES = [1, 2, 3, 4, 6, 8]
VARIANTS = {
    "default": {},
    "ang 0.35": dict(ang_tol_deg=0.35),
    "ang 1.4": dict(ang_tol_deg=1.4),
    "pair 0.5": dict(pair_tol_deg=0.5),
    "pair 2": dict(pair_tol_deg=2.0),
    "cluster 0.5": dict(cluster_deg=0.5),
    "cluster 2": dict(cluster_deg=2.0),
}


def recall_row(res: RI.Candidates, Rt: np.ndarray, sym: np.ndarray) -> Dict[str, Any]:
    d = RI.misorientation_deg(res.R, Rt[None], sym) if len(res.R) else np.zeros(0)
    return {f"k{K}": bool(len(d[:K]) and d[:K].min() < 2.0) for K in (1, 5, 20)}


def main() -> None:
    stored = AC.pilot_stored()
    AC.init("full")
    base = RI.IndexParams()
    data: Dict[str, Dict[str, List[Any]]] = {}
    cases = [(f"stride {m}", {}, m) for m in STRIDES] + [
        (name, kw, 1) for name, kw in VARIANTS.items() if kw
    ]
    for v in sorted(stored):
        e = AC.env("full", v)
        sc = e["sc"]
        crystal = RI.crystal_from_smooth_cost(sc, AC.symmetry_of(e), e["phase"], base.q_max)
        att = sc.attached_rays(e["verts"])
        g_all = RI.scattering_vectors(sc, att)
        for name, kw, m in cases:
            p = replace(base, **kw)
            g = g_all[::m]
            res = RI.index_rays(g, crystal, p, 20)
            row = recall_row(res, e["R_true"], crystal.sym)
            row["n_true"] = int(
                RI.explained_mask(e["R_true"], g, crystal, p.ang_tol_deg, p.len_tol).sum()
            )
            row["seconds"] = res.seconds.get("total", 0.0)
            data.setdefault(name, {}).setdefault("rows", []).append(row)
        print(v, flush=True)
    out: Dict[str, Any] = {}
    for name, d in data.items():
        rows = d["rows"]
        out[name] = dict(
            n=len(rows),
            **{k: int(sum(r[k] for r in rows)) for k in ("k1", "k5", "k20")},
            median_true_rays=float(np.median([r["n_true"] for r in rows])),
            median_seconds=float(np.median([r["seconds"] for r in rows])),
        )
        print(name, out[name])
    AC.OUT.mkdir(parents=True, exist_ok=True)
    (AC.OUT / "sparse_pilot.json").write_text(json.dumps(out, indent=1) + "\n")


if __name__ == "__main__":
    main()
