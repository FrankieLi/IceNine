"""AR1 steps 2-5: recall of the ray-indexing proposer on E0 and on the pilot's full-search voxels.

Per voxel: the attached rays, how many the true orientation explains, the candidates (top K),
their distance to the truth (symmetry-reduced), the time, and for a miss the evidence for its
reason (rays explained by the truth, votes near the truth, the best cluster's votes).

  uv run python scripts/autoresearch/ar1/recall.py e0 [--n 200] [--tag base]
  uv run python scripts/autoresearch/ar1/recall.py pilot [--tag base]
Writes benchmarks/autoresearch/ar1/recall_<kind>_<tag>.json and cand_<kind>_<tag>.npz.
"""

import argparse
import json
import sys
import time
from dataclasses import asdict
from pathlib import Path
from typing import Any, Dict, List

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import ar1_common as AC  # noqa: E402

from icenine import ray_indexing as RI  # noqa: E402

TOP_K = 20


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("kind", choices=["e0", "pilot"])
    ap.add_argument("--n", type=int, default=0)
    ap.add_argument("--tag", default="base")
    ap.add_argument("--set", nargs="*", default=[], help="IndexParams overrides key=value")
    a = ap.parse_args()
    p = RI.IndexParams()
    for kv in a.set:
        k, v = kv.split("=")
        setattr(p, k, type(getattr(p, k))(v))
    kind = "e0" if a.kind == "e0" else "full"
    stored = AC.e0_stored() if a.kind == "e0" else AC.pilot_stored()
    voxels = sorted(stored)
    if a.n:
        voxels = voxels[: a.n]
    AC.init(kind)
    rows: List[Dict[str, Any]] = []
    cand_R: Dict[str, np.ndarray] = {}
    for v in voxels:
        e = AC.env(kind, v)
        sc = e["sc"]
        crystal = RI.crystal_from_smooth_cost(sc, AC.symmetry_of(e), e["phase"], p.q_max)
        sym = crystal.sym
        t0 = time.perf_counter()
        res, att = RI.index_voxel(sc, e["verts"], crystal, p, TOP_K, keep_votes=True)
        wall = time.perf_counter() - t0
        g = RI.scattering_vectors(sc, att)
        Rt = e["R_true"]
        n_true = int(RI.explained_mask(Rt, g, crystal, p.ang_tol_deg, p.len_tol).sum())
        if len(res.R):
            dist = RI.misorientation_deg(res.R, Rt[None], sym)
        else:
            dist = np.zeros(0)
        if res.vote_R is not None and len(res.vote_R):
            vd = RI.misorientation_deg(res.vote_R, Rt[None], sym)
            votes_truth = int((vd < p.cluster_deg).sum())
        else:
            votes_truth = 0
        # strict refinement of the truth for the "explained if refined" comparison
        _, n_true_ref = RI.refine_rotation(Rt, g, crystal, p)
        rows.append(
            dict(
                voxel=v,
                stored_err=stored[v]["err"],
                stored_evals=stored[v]["evals"],
                n_att=int(len(att)),
                n_true_explained=n_true,
                n_true_explained_refined=n_true_ref,
                n_votes=int(res.n_votes),
                n_pairs=int(res.n_pairs),
                votes_at_truth=votes_truth,
                dist=[float(x) for x in dist],
                explained=[int(x) for x in res.explained],
                cluster_votes=[int(x) for x in res.votes],
                seconds={k: float(x) for k, x in res.seconds.items()},
                wall_s=wall,
            )
        )
        cand_R[f"v{v}"] = res.R
        r = rows[-1]
        print(
            v,
            f"att={r['n_att']} true_expl={n_true} votes={res.n_votes} at_truth={votes_truth}",
            f"best={min(dist) if len(dist) else float('nan'):.2f}",
            f"rank1={dist[0] if len(dist) else float('nan'):.2f}",
            f"{wall:.2f}s",
            flush=True,
        )
    AC.OUT.mkdir(parents=True, exist_ok=True)
    stem = f"{a.kind}_{a.tag}"
    (AC.OUT / f"recall_{stem}.json").write_text(
        json.dumps(dict(params=asdict(p), rows=rows)) + "\n"
    )
    np.savez_compressed(AC.OUT / f"cand_{stem}.npz", **cand_R)


if __name__ == "__main__":
    main()
