"""AR2 step 1: the finisher arms (best 1, best 2, best 3) on the holdout voxels.

Per voxel: attach rays, index (first 20 candidates), recall of the truth, strict cost on the first 5
candidates, the finisher (``refine_from_candidates``, rng seed 0, as AR1's e2e.py) on the three
lowest-strict-cost of those 5, one run each. Arm "best N" keeps the lowest final cost of the first N
runs; its evaluations are 5 strict + the N finisher runs. The voxel sets are read from voxels.json
(fixed and committed before this ran).

  uv run python scripts/autoresearch/ar2/holdout_run.py run --part 0/3
  uv run python scripts/autoresearch/ar2/holdout_run.py check      # AR1 candidates reproduce?
Writes benchmarks/autoresearch/ar2/holdout/run_part<i>of<n>.json and cand_part<i>of<n>.npz.
"""

import argparse
import json
import sys
import time
from pathlib import Path
from typing import Any, Dict, List

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "ar1"))
import ar1_common as AC  # noqa: E402

from icenine import ray_indexing as RI  # noqa: E402

OUT = AC.ICE / "benchmarks" / "autoresearch" / "ar2" / "holdout"
TOP_K = 20
K_STRICT = 5
N_ARMS = (1, 2, 3)


def load_sets() -> Dict[int, str]:
    d = json.loads((OUT / "voxels.json").read_text())
    return {v: name for name, s in d["sets"].items() for v in s["voxels"]}


def check() -> None:
    """Re-index 3 pilot voxels and compare with AR1's stored candidates."""
    z = np.load(AC.OUT / "cand_pilot_base.npz")
    keys = sorted(z.files)
    pick = [keys[0], keys[len(keys) // 2], keys[-1]]
    p = RI.IndexParams()
    AC.init("full")
    res_rows = []
    for k in pick:
        v = int(k[1:])
        e = AC.env("full", v)
        sc = e["sc"]
        crystal = RI.crystal_from_smooth_cost(sc, AC.symmetry_of(e), e["phase"], p.q_max)
        res, _ = RI.index_voxel(sc, e["verts"], crystal, p, TOP_K)
        same_shape = res.R.shape == z[k].shape
        dmax = float(np.abs(res.R - z[k]).max()) if same_shape else float("nan")
        res_rows.append(dict(voxel=v, n=len(res.R), n_stored=len(z[k]), max_abs_diff=dmax))
        print(res_rows[-1], flush=True)
    (OUT / "ar1_reproduction_check.json").write_text(json.dumps(res_rows, indent=1) + "\n")


def run(part: str) -> None:
    sets = load_sets()
    voxels = sorted(sets)
    pi, pn = (int(x) for x in part.split("/"))
    voxels = voxels[pi::pn]
    p = RI.IndexParams()
    AC.init("full")
    from icenine.reconstructor import SearchCandidate

    rows: List[Dict[str, Any]] = []
    cands: Dict[str, np.ndarray] = {}
    stem = f"part{pi}of{pn}"
    for v in voxels:
        e = AC.env("full", v)
        sc, rec, lf = e["sc"], e["rec"], e["lf"]
        verts, phase, Rt = e["verts"], e["phase"], e["R_true"]
        crystal = RI.crystal_from_smooth_cost(sc, AC.symmetry_of(e), e["phase"], p.q_max)
        sym = RI.proper_rotations(AC.symmetry_of(e))
        t0 = time.perf_counter()
        res, att = RI.index_voxel(sc, verts, crystal, p, TOP_K)
        index_s = time.perf_counter() - t0
        g = RI.scattering_vectors(sc, att)
        n_true = int(RI.explained_mask(Rt, g, crystal, p.ang_tol_deg, p.len_tol).sum())
        R = res.R
        row: Dict[str, Any] = dict(
            voxel=v,
            set=sets[v],
            n_att=int(len(att)),
            n_true_explained=n_true,
            n_cand=int(len(R)),
            index_s=index_s,
        )
        cands[f"v{v}"] = R
        if len(R) == 0:
            row["no_candidates"] = True
            rows.append(row)
            continue
        dist = RI.misorientation_deg(R, Rt[None], sym)
        k5 = min(K_STRICT, len(R))
        costs = np.array(
            [float(lf.evaluate(R[j], verts, phase).cost) for j in range(k5)], dtype=np.float64
        )
        order = [int(j) for j in np.argsort(costs, kind="stable")[:3]]
        fins: List[Dict[str, Any]] = []
        for j in order:
            t1 = time.perf_counter()
            r = rec.refine_from_candidates(
                [SearchCandidate(orientation=R[j].copy(), cost=float(costs[j]))],
                verts,
                phase,
                rng=np.random.default_rng(0),
            )
            Rf = np.asarray(r.orientation, dtype=np.float64)
            fins.append(
                dict(
                    cand=j,
                    cost_before=float(costs[j]),
                    cost_after=float(r.cost),
                    evals=int(rec.last_eval_counts[1]),
                    seconds=time.perf_counter() - t1,
                    err=float(RI.misorientation_deg(Rf[None], Rt[None], sym)[0]),
                    R_final=Rf.reshape(-1).tolist(),
                )
            )
        arms = {}
        for n in N_ARMS:
            use = fins[: min(n, len(fins))]
            best = min(use, key=lambda f: f["cost_after"])
            arms[str(n)] = dict(
                err=best["err"],
                wrong=bool(best["err"] > 1.0),
                taken_cand=best["cand"],
                evals=int(k5 + sum(f["evals"] for f in use)),
                n_refined=len(use),
            )
        row.update(
            dist=[float(d) for d in dist],
            strict_costs=[float(c) for c in costs],
            explained=[int(x) for x in res.explained],
            finishers=fins,
            arms=arms,
        )
        rows.append(row)
        print(
            v,
            sets[v],
            f"att={row['n_att']} true_expl={n_true} best20={dist.min():.2f}",
            " ".join(f"b{n}={arms[str(n)]['err']:.3f}/{arms[str(n)]['evals']}" for n in N_ARMS),
            f"idx={index_s:.1f}s",
            flush=True,
        )
        (OUT / f"run_{stem}.json").write_text(json.dumps(dict(rows=rows, label="contended")) + "\n")
    (OUT / f"run_{stem}.json").write_text(json.dumps(dict(rows=rows, label="contended")) + "\n")
    np.savez_compressed(OUT / f"cand_{stem}.npz", **cands)
    (OUT / f"run_{stem}.done").write_text("done\n")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("what", choices=["run", "check"])
    ap.add_argument("--part", default="0/1")
    a = ap.parse_args()
    check() if a.what == "check" else run(a.part)
