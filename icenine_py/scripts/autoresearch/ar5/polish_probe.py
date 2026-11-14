"""AR5 step 14, phase A: the probe replay of rule ``G2`` from the stored rescue ends.

  RESEED_SAMPLE=<mid|grad|held|held2|big> uv run python scripts/autoresearch/ar5/polish_probe.py \
      [--forms G2a,G2b] [--limit N]

1. Stored counts (``polish_cases.counts``): failed rescues, their grid evaluations and the run total
   with the pass must equal the plan's; else exit 1 (stop rule).
2. Prediction 1: every stored rescue end of every failed rescue is re-scored (strict cost and hit
   ratio of one overlap evaluation); the stored cost and hit ratio must agree to 1e-9; else exit 1.
3. ``G2`` of each form from the rescue ends in ``G1``'s order, first end at 0.95 taken
   (``polish_rules.polish_choose``). Per try: start, end hit ratio, cost, evaluations, error to the
   truth (symmetry-reduced; scoring only), angle to the grid end. Per case: accepted or not,
   evaluations of the probe, the grid evaluations it would save, the grid end's hit ratio and error.

Writes ``polish/probe_<form>_<sample>.json`` and ``polish/pred1_<sample>.json``. ``held2`` runs only
after the freeze commit (``polish_common`` refuses it before). It replays the whole changed path at
the voxel (``G2`` acts only after a failed rescue); growth after the voxel is not replayed: a cost
estimate and a voxel check, not the result.
"""

import argparse
import json
import os
import sys
import time
from typing import Any, Dict, List

import polish_common as PC  # noqa: E402  (first)

import numpy as np  # noqa: E402
import polish_cases as CS  # noqa: E402
import polish_rules as PR  # noqa: E402
import reseed_common as RC  # noqa: E402

TOL = 1e-9


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--forms", default="G2a,G2b")
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--tag", default="")
    a = ap.parse_args()
    s = RC.SAMPLE
    forms = a.forms.split(",") if a.forms else []
    PC.PO.mkdir(parents=True, exist_ok=True)
    c = CS.counts(s)
    print(json.dumps(c), flush=True)
    if not c["ok"]:
        (PC.PO / f"counts_{s}.json").write_text(json.dumps(c, indent=1) + "\n")
        print("STOP: a stored count does not reproduce", file=sys.stderr)
        sys.exit(1)
    cases = CS.failed_cases(s)
    if a.limit:
        cases = cases[: a.limit]
    ge = CS.grid_ends(s, [int(x["voxel"]) for x in cases])
    env = RC.make_env()
    rec = env.hooks.rec
    T = env.T
    orig_lo = rec.local_optimization
    ratios = env.bfs._ratios
    t0 = time.time()

    # ---- prediction 1: re-score every stored rescue end
    rows1: List[Dict[str, Any]] = []
    worst = 0.0
    for cs in cases:
        v = int(cs["voxel"])
        verts = env.hooks.vertices(v)
        ph = int(env.mic.voxels[v].phase)
        for t in cs["tries"]:
            info = rec.evaluate_overlap(np.asarray(t["R"], np.float64).reshape(3, 3), verts, ph)
            dc = abs(float(info.cost) - float(t["cost"]))
            dh = abs(float(ratios(info)[0]) - float(t["hit"]))
            worst = max(worst, dc, dh)
            rows1.append(dict(voxel=v, k=t["k"], d_cost=dc, d_hit=dh))
    p1 = dict(
        sample=s,
        ends=len(rows1),
        worst_abs_diff=worst,
        tol=TOL,
        ok=bool(worst <= TOL),
        rows=rows1,
    )
    (PC.PO / f"pred1_{s}{a.tag}.json").write_text(json.dumps(p1, indent=1) + "\n")
    print(f"prediction 1: {len(rows1)} ends, worst diff {worst:.3g}, ok={p1['ok']}", flush=True)
    if not p1["ok"]:
        print("STOP: prediction 1 fails", file=sys.stderr)
        sys.exit(1)

    # ---- G2 probe
    for form in forms:
        out_rows: List[Dict[str, Any]] = []
        for cs in cases:
            v = int(cs["voxel"])
            verts = env.hooks.vertices(v)
            ph = int(env.mic.voxels[v].phase)
            g = ge[v]
            g_err = float(T.err(g["R"][None], [v])[0])
            ends = [np.asarray(t["R"], np.float64).reshape(3, 3) for t in cs["tries"]]

            def two(j: int) -> Dict[str, Any]:
                r = PC.g2_try(form, rec, orig_lo, ratios, ends[j], verts, ph)
                R = np.asarray(r["R"], np.float64)
                err = float(T.err(R[None], [v])[0])
                ang = float(env.mis(R[None], g["R"][None])[0])
                ang_start = float(env.mis(R[None], ends[j][None])[0])
                return dict(
                    hit=r["hit"],
                    cost=r["cost"],
                    evals=r["evals"],
                    grid_evals=r["grid_evals"],
                    err=err,
                    angle_to_grid_end=ang,
                    angle_to_rescue_end=ang_start,
                    R=R,
                )

            ch = PR.polish_choose(len(ends), two)
            tk = ch["tries"][ch["taken"]] if ch["accepted"] else None
            rescue_err = [float(T.err(e[None], [v])[0]) for e in ends]
            out_rows.append(
                dict(
                    sample=s,
                    voxel=v,
                    form=form,
                    n_ends=len(ends),
                    accepted=ch["accepted"],
                    n_tries=len(ch["tries"]),
                    probe_evals=ch["evals"],
                    grid_evals=int(cs["grid_evals"]),
                    saved=int(cs["grid_evals"]) if ch["accepted"] else 0,
                    grid_end_hit=g["hit"],
                    grid_end_err=g_err,
                    grid_end_reached_095=bool(g["hit"] >= PR.TAU),
                    rescue_end_errs=rescue_err,
                    taken_err=None if tk is None else tk["err"],
                    taken_hit=None if tk is None else tk["hit"],
                    taken_angle_to_grid_end=None if tk is None else tk["angle_to_grid_end"],
                    wrong_accepted=bool(tk is not None and tk["err"] > 1.0),
                    tries=[
                        {
                            k: (x.ravel().tolist() if hasattr(x, "ravel") else x)
                            for k, x in t.items()
                        }
                        for t in ch["tries"]
                    ],
                )
            )
            print(
                f"{s} {form} voxel {v}: accepted={ch['accepted']} tries={len(ch['tries'])} "
                f"evals={ch['evals']} {time.time() - t0:.0f}s",
                flush=True,
            )
        (PC.PO / f"probe_{form}_{s}{a.tag}.json").write_text(
            json.dumps(
                dict(
                    sample=s,
                    form=form,
                    rule=PR.RULE,
                    wall_s=time.time() - t0,
                    wall_label="contended; not a timing result",
                    rows=out_rows,
                ),
                indent=1,
            )
            + "\n"
        )
    os._exit(0)


if __name__ == "__main__":
    main()
