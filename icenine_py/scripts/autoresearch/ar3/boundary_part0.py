"""AR3 step 5, part 0: can the pass work at all? Records and single strict scores, no search.

  uv run python scripts/autoresearch/ar3/boundary_part0.py

For each of the six stored maps: the number of boundary voxels and of candidates (the cost of the
pass before it runs), and for every wrong voxel (error above 1 degree) whether it is a boundary
voxel, whether a candidate lies within 1 degree of its true orientation, and one strict score of
the true orientation and of the voxel's end (cost, hit ratio, confidence). The truth is used here
for analysis and for the gate only. Writes ``boundary/part0.json``.

Gate (plan): over the newly wrong voxels of the three step 4 maps (wrong in the step 4 map, right in
the map it was compared with), if fewer than half have both a candidate within 1 degree of the
truth and a truth that scores better than the end by at least one of the three measures, stop.
"""

import json
import time
from typing import Any, Dict, List

import boundary_common as B
import numpy as np


def main() -> None:
    B.BND.mkdir(parents=True, exist_ok=True)
    out: Dict[str, Any] = dict(maps={}, rows={})
    envs: Dict[str, B.Env] = {}
    errs: Dict[str, np.ndarray] = {}
    maps: Dict[str, np.ndarray] = {}
    for name, (ek, _desc, _cmp) in B.MAPS.items():
        if ek not in envs:
            t0 = time.time()
            envs[ek] = B.Env(ek)
            print("env", ek, f"{time.time() - t0:.0f}s", flush=True)
        env = envs[ek]
        maps[name] = B.load_map(name, env.n)
        errs[name] = env.err(maps[name])
    for name, (ek, _desc, cmp_) in B.MAPS.items():
        env = envs[ek]
        R, e = maps[name], errs[name]
        isb, cands = env.boundary(R)
        wrong = e > B.WRONG
        newly = wrong & ~(errs[cmp_] > B.WRONG) if cmp_ else np.zeros(env.n, bool)
        n_c = np.array([len(cands.get(v, [])) for v in range(env.n)])
        st = dict(
            n_voxels=env.n,
            n_wrong=int(wrong.sum()),
            n_newly_wrong=int(newly.sum()),
            n_boundary=int(isb.sum()),
            n_candidates=int(n_c.sum()),
            candidates_by_count={
                str(k): int((n_c == k).sum()) for k in range(1, int(n_c.max()) + 1)
            },
            max_evals_of_pass=int(n_c.sum()) * (B.RF.CMA_MAX_EVALS + 1),
            wrong_that_are_boundary=int((wrong & isb).sum()),
        )
        rows: List[Dict[str, Any]] = []
        t0 = time.time()
        for v in np.where(wrong)[0]:
            v = int(v)
            cs = cands.get(v, [])
            c_err = (
                env.mis(R[cs], np.repeat(env.Rtrue[v][None], len(cs), 0)) if cs else np.array([])
            )
            tr = env.score3(env.Rtrue[v], v)
            en = env.score3(R[v], v)
            rows.append(
                dict(
                    voxel=v,
                    err=float(e[v]),
                    level=None if env.level is None else float(env.level[v]),
                    grain=int(env.gmap[v]),
                    isA=bool(env.isA[v]),
                    newly=bool(newly[v]),
                    boundary=bool(isb[v]),
                    n_cand=len(cs),
                    cand_err_min=float(c_err.min()) if len(cs) else None,
                    cand_within_1=bool(len(cs) and c_err.min() <= 1.0),
                    truth=dict(cost=tr[0], hit=tr[1], conf=tr[2]),
                    end=dict(cost=en[0], hit=en[1], conf=en[2]),
                    better_cost=bool(tr[0] < en[0]),
                    better_hit=bool(tr[1] > en[1]),
                    better_conf=bool(tr[2] > en[2]),
                )
            )
        print(
            name,
            st["n_boundary"],
            "boundary",
            st["n_candidates"],
            "candidates",
            len(rows),
            "wrong scored",
            f"{time.time() - t0:.0f}s",
            flush=True,
        )
        out["maps"][name] = st
        out["rows"][name] = rows
    # ---- gate
    tot = ok = 0
    gate: Dict[str, Any] = {}
    for name in B.STEP4:
        rs = [r for r in out["rows"][name] if r["newly"]]
        good = [
            r
            for r in rs
            if r["cand_within_1"] and (r["better_cost"] or r["better_hit"] or r["better_conf"])
        ]
        gate[name] = dict(newly_wrong=len(rs), both=len(good))
        tot += len(rs)
        ok += len(good)
    gate["all"] = dict(newly_wrong=tot, both=ok, share=ok / tot if tot else None)
    gate["passes"] = bool(tot and ok >= tot / 2)
    out["gate"] = gate
    (B.BND / "part0.json").write_text(json.dumps(out, default=float) + "\n")
    print("gate", json.dumps(gate), flush=True)


if __name__ == "__main__":
    main()
