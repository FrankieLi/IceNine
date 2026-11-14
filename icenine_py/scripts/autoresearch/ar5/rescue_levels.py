"""AR5 step 9, phase A step 5: the per-level record of the 23 disclosed grid fallbacks.

  RESEED_SAMPLE=<grad|mid> uv run python scripts/autoresearch/ar5/rescue_levels.py [--limit N]

Runs the grid fallback again on each failed seed step of arm ``A`` (``rescue/cases_<sample>.json``)
with the draws arm ``A`` gave it (aligned generator of (voxel, "grid", 0)) and records, through the
reconstructor's own ``recorder`` hook (it draws nothing; results are bit-identical with it off),
what the stored runs do not hold, per level: the discrete stage (FZ orientations x local grid,
evaluations, candidates, their scores), the quick MC (ranked candidates, costs, how many are
kept), the level best, and after the last level FindOptimal's candidates, the variance step and
the final end. Truth errors of the level bests and of the final end are added for scoring only.
The grid fallback's evaluations (global + local) are summed per case; the totals must equal arm
``A``'s stored 3,726,723 (``grad``) / 4,007,241 (``mid``) when all cases of a sample are run.
Per-level evaluation counts of the quick MC are not separable (the library keeps one local
counter); the discrete stage's are exact (printed). Writes the small summary
``cache/levels_<sample>.json`` (untracked) and the arrays ``cache/levels/<sample>_<voxel>.npz``
(untracked). Data for a later learned step; nothing here feeds a decision of this step.
"""

import argparse
import contextlib
import io
import json
import os
import re
import time
from typing import Any, Dict, List

import rescue_common as RCM  # noqa: E402  (first)

import fit_common as FC  # noqa: E402
import numpy as np  # noqa: E402
import order_draws as OD  # noqa: E402
import reseed_common as RC  # noqa: E402

DISC = re.compile(r"Level (\d+): discrete search \((\d+) FZ × (\d+) local = (\d+) evals")
EXPECT = {"grad": 3_726_723, "mid": 4_007_241}


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--limit", type=int, default=0)
    a = ap.parse_args()
    s = FC.SAMPLE
    assert s in ("grad", "mid")
    cases = json.loads((RCM.RESCUE / f"cases_{s}.json").read_text())
    if a.limit:
        cases = cases[: a.limit]
    env = RC.make_env()
    rec = env.hooks.rec
    T = env.T
    out_dir = RCM.RCACHE / "levels"
    out_dir.mkdir(parents=True, exist_ok=True)
    events: List[Any] = []
    rec.recorder = lambda ev, data: events.append((ev, data))
    rows: List[Dict[str, Any]] = []
    t0 = time.time()
    for c in cases:
        v = int(c["voxel"])
        verts = env.hooks.vertices(v)
        ph = int(env.mic.voxels[v].phase)
        events.clear()
        buf = io.StringIO()
        with contextlib.redirect_stdout(buf):
            res = rec.reconstruct_voxel(
                voxel_vertices=verts,
                phase_index=ph,
                rng=OD.make_rng(OD.voxel_key(verts), "grid", 0),
            )
        g, loc, tot = rec.last_eval_counts
        disc = {
            int(m.group(1)): dict(fz=int(m.group(2)), local=int(m.group(3)), evals=int(m.group(4)))
            for m in DISC.finditer(buf.getvalue())
        }
        arrays: Dict[str, Any] = {}
        levels: List[Dict[str, Any]] = []
        for ev, d in events:
            if ev == "discrete":
                L = int(d["level"])
                arrays[f"disc_R_{L}"] = np.asarray(d["R"], np.float32)
                arrays[f"disc_score_{L}"] = np.asarray(d["score"], np.float32)
                levels.append(
                    dict(
                        level=L,
                        n_q_max=float(d["n_q_max"]),
                        diameter_deg=float(np.degrees(d["diameter"])),
                        n_candidates=int(len(d["score"])),
                        best_discrete_score=float(np.min(d["score"])),
                        **disc.get(L, {}),
                    )
                )
            elif ev == "quick_mc":
                L = levels[-1]["level"]
                arrays[f"qmc_R_{L}"] = np.asarray(d["R"], np.float32)
                arrays[f"qmc_cost_{L}"] = np.asarray(d["cost"], np.float32)
                levels[-1].update(
                    n_keep=int(d["n_keep"]),
                    best_qmc_cost=float(np.min(d["cost"])),
                    best_qmc_err=float(
                        T.err(np.asarray(d["R"])[[int(np.argmin(d["cost"]))]], [v])[0]
                    ),
                )
        fo = [d for ev, d in events if ev == "find_candidate"]
        final = [d for ev, d in events if ev == "final"]
        row = dict(
            sample=s,
            voxel=v,
            step=c["step"],
            global_evals=int(g),
            local_evals=int(loc),
            grid_evals=int(tot),
            levels=levels,
            find_candidates=len(fo),
            find_best_cost=float(min(d["cost"] for d in fo)) if fo else None,
            final_cost=float(res.cost),
            final_err=float(T.err(np.asarray(res.orientation)[None], [v])[0]),
            grid_end_err_in_run_a=float(c["grid_end_err"]),
        )
        if final:
            arrays["final_R"] = np.asarray(final[0]["R"], np.float64)
        np.savez_compressed(out_dir / f"{s}_{v}.npz", **arrays)
        rows.append(row)
        print(
            f"{s} voxel {v}: evals={tot} final_err={row['final_err']:.4f} {time.time()-t0:.0f}s",
            flush=True,
        )
    total = int(sum(r["grid_evals"] for r in rows))
    summ = dict(
        sample=s,
        cases=len(rows),
        grid_evals_total=total,
        stored_total_in_arm_A=EXPECT[s],
        reproduces_stored_total=bool(total == EXPECT[s]) if not a.limit else None,
        wall_s=time.time() - t0,
        wall_label="contended; not a timing result",
        rows=rows,
    )
    name = f"levels_{s}" + (f"_limit{a.limit}" if a.limit else "")
    (RCM.RCACHE / f"{name}.json").write_text(json.dumps(summ, indent=1) + "\n")
    os._exit(0)


if __name__ == "__main__":
    main()
