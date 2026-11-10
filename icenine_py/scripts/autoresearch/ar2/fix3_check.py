"""Review fix 3 of AR2 step 1: does the robust ``misorientation_deg`` change any stored result?

The library function is used inside the indexer (merging near-duplicate candidates at 1 degree)
and by the AR1 and AR2 scripts for distances to the truth. Per set, this script re-indexes every
voxel with the current code and compares the candidates with the stored ones, then recomputes the
candidate distances to the truth with the current function and compares them with the stored
distances (computed with the old ``arccos`` form).

  uv run python scripts/autoresearch/ar2/fix3_check.py holdout   # AR2's 300 holdout voxels
  uv run python scripts/autoresearch/ar2/fix3_check.py pilot     # AR1's 69 pilot voxels
  uv run python scripts/autoresearch/ar2/fix3_check.py e0        # AR1's 200 E0 voxels
Writes benchmarks/autoresearch/ar2/holdout/fix3_check_<set>.json.
"""

import argparse
import json
import sys
from pathlib import Path
from typing import Any, Dict, List, Tuple

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "ar1"))
import ar1_common as AC  # noqa: E402

from icenine import ray_indexing as RI  # noqa: E402

OUT = AC.ICE / "benchmarks" / "autoresearch" / "ar2" / "holdout"
TOP_K = 20
TOLS = (0.5, 1.0, 2.0, 5.0)


def stored(what: str) -> Tuple[str, Dict[int, np.ndarray], Dict[int, List[float]]]:
    """Kind for ``AC.init``, stored candidates and stored distances per voxel."""
    cand: Dict[int, np.ndarray] = {}
    dist: Dict[int, List[float]] = {}
    if what == "holdout":
        for f in sorted(OUT.glob("cand_part*.npz")):
            z = np.load(f)
            cand.update({int(k[1:]): z[k] for k in z.files})
        for f in sorted(OUT.glob("run_part*.json")):
            for r in json.loads(f.read_text())["rows"]:
                dist[r["voxel"]] = r.get("dist", [])
        return "full", cand, dist
    rows, c = AC.load_recall(what)
    cand = {int(k[1:]): v for k, v in c.items()}
    dist = {r["voxel"]: r["dist"] for r in rows}
    return ("e0" if what == "e0" else "full"), cand, dist


def ar1_medians(
    what: str, old: Dict[int, List[float]], new: Dict[int, np.ndarray]
) -> Dict[str, Any]:
    """AR1's candidate-distance medians (table ``ar1_dist``, as AR1 summary.py computes them) and
    the no-finisher ``med_err`` of its end-to-end table, with the old and the new distances. The
    finisher arms' errors cannot be recomputed: e2e.py did not store the final orientations."""

    def dist_stats(D: List[np.ndarray]) -> Dict[str, float]:
        best20 = np.array([d[:20].min() if len(d) else np.inf for d in D])
        top1 = np.array([d[0] if len(d) else np.inf for d in D])
        ok = best20 < 2.0
        return dict(
            best_med=float(np.median(best20[ok])),
            best_p90=float(np.percentile(best20[ok], 90)),
            best_max=float(best20[ok].max()),
            top1_med=float(np.median(top1[top1 < 2.0])),
        )

    vs = sorted(old)
    res: Dict[str, Any] = dict(
        old=dist_stats([np.asarray(old[v]) for v in vs]), new=dist_stats([new[v] for v in vs])
    )
    tag = "e0_base_part" if what == "e0" else "pilot_base_part"
    files = sorted(AC.OUT.glob(f"e2e_{tag}*.json"))
    rows = [r for f in files for r in json.loads(f.read_text())["rows"]]
    rows = [r for r in rows if not r.get("no_candidates")]
    for K in ("1", "5", "20"):
        o = np.array([r["err_pick"][K] for r in rows])
        n = np.array([float(new[r["voxel"]][r["pick"][K]]) for r in rows])
        res[f"no_finisher_first{K}_med_err"] = dict(
            old=float(np.median(o)), new=float(np.median(n))
        )
    return res


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("what", choices=["holdout", "pilot", "e0"])
    a = ap.parse_args()
    kind, cand, dist_old = stored(a.what)
    p = RI.IndexParams()
    AC.init(kind)
    per: List[Dict[str, Any]] = []
    d_old_all: List[float] = []
    d_new_all: List[float] = []
    new_by_voxel: Dict[int, np.ndarray] = {}
    for v in sorted(cand):
        e = AC.env(kind, v)
        sc = e["sc"]
        crystal = RI.crystal_from_smooth_cost(sc, AC.symmetry_of(e), e["phase"], p.q_max)
        res, _ = RI.index_voxel(sc, e["verts"], crystal, p, TOP_K)
        same = res.R.shape == cand[v].shape
        diff = float(np.abs(res.R - cand[v]).max()) if same and len(res.R) else 0.0
        d_new = (
            RI.misorientation_deg(cand[v], e["R_true"][None], crystal.sym)
            if len(cand[v])
            else np.zeros(0)
        )
        d_old = np.asarray(dist_old[v], dtype=np.float64)
        per.append(dict(voxel=v, same_shape=same, max_abs_diff=diff, n=int(len(res.R))))
        d_old_all += d_old.tolist()
        d_new_all += d_new.tolist()
        new_by_voxel[v] = d_new
        print(v, same, f"{diff:.2e}", flush=True)
    do, dn = np.array(d_old_all), np.array(d_new_all)
    small = do < 1.0
    out = dict(
        set=a.what,
        voxels=len(per),
        candidates_identical=int(sum(r["same_shape"] and r["max_abs_diff"] == 0.0 for r in per)),
        candidates_max_abs_diff=max(r["max_abs_diff"] for r in per),
        candidates_same_shape=int(sum(r["same_shape"] for r in per)),
        distances=int(len(do)),
        dist_max_abs_change=float(np.abs(dn - do).max()),
        dist_max_abs_change_below_1deg=float(np.abs(dn - do)[small].max()) if small.any() else 0.0,
        dist_max_abs_change_above_1deg=(
            float(np.abs(dn - do)[~small].max()) if (~small).any() else 0.0
        ),
        dist_exact_zero_old=int((do == 0.0).sum()),
        dist_exact_zero_new=int((dn == 0.0).sum()),
        threshold_flips={f"{t:g}": int(((do < t) != (dn < t)).sum()) for t in TOLS},
        per_voxel=per,
    )
    if a.what != "holdout":
        out["ar1_medians"] = ar1_medians(a.what, dist_old, new_by_voxel)
    (OUT / f"fix3_check_{a.what}.json").write_text(json.dumps(out, indent=1) + "\n")
    print({k: v for k, v in out.items() if k != "per_voxel"})


if __name__ == "__main__":
    main()
