"""AR5 step 12, phase A step 4: the gate replay of rule PH(3000) along the whole changed path.

  RESEED_SAMPLE=<big|grad|mid> uv run python scripts/autoresearch/ar5/stallhit_replay.py \
      [--k 3000] [--shard I/N] [--limit M] [--voxels V1,V2] [--tag NAME]

For every seed step of arm ``G`` (step 9) that ran an MC refinement it runs the step live, with the
running code of the whole-layer runs (``order_run.OrderPick.e1_step`` through ``stallhit_run.
StallHitPick``: interleaved order, 20 costed, 5 refined) and rule ``PH(k)`` on the MC refinements,
and compares every fit with the stored one (``e1_log.json.gz``). The seed step is deterministic and
independent of the rest of the run (every fit has a fresh generator), so a replay in isolation is
the same step. Per fit: ``same`` (hit ratio, cost, end equal; evaluations equal to the stored ones
plus one for each stall check that went on) or ``after_last`` (the same end, cost and hit ratio,
fewer evaluations: stopped after its last gain) or ``moved`` (any other difference, expected only
where the stored trace has a stall before a later gain). The class the stored trace allows
(``stallhit_rules.stall_events``; a stop needs a stall) is checked per fit (``pred_mismatch``). The
hit ratio of every stall check is stored. Along the path of a step whose kept ends changed: the seed
test (0.80) on the new taken end; if it fails, rule ``G1``'s rescue (``rescue_rules.rescue_choose``,
``rescue_common.cma_rescue``) and, if nothing is rescued, the grid search with the draws arm ``G``
gave it (generator of (voxel, "grid", 0)). A step whose fits are the same or ``after_last`` has the
path of ``G`` (its stored rescue and grid evaluations are used). Truth (symmetry-reduced angle) is
for scoring only. Writes ``stallhit/replay_<sample>[_sI]<tag>.json`` (one row per step). ``held``
and ``held2`` are not accepted.
"""

import argparse
import contextlib
import io
import json
import os
import time
from typing import Any, Dict, List, Optional

import stallhit_common as SC  # noqa: E402  (first)

import numpy as np  # noqa: E402
import order_draws as OD  # noqa: E402
import order_run as ORUN  # noqa: E402
import stallhit_rules as SR  # noqa: E402
import stallhit_run as SRUN  # noqa: E402
import stallhit_whatif as SW  # noqa: E402
import reseed_common as RC  # noqa: E402
import rescue_common as RCM  # noqa: E402
import rescue_rules as RR  # noqa: E402

TOL_R = 1e-9
TOL_X = 1e-12


def same_R(a: Any, b: Any) -> bool:
    x = np.asarray(a, np.float64).ravel()
    y = np.asarray(b, np.float64).ravel()
    return bool(np.max(np.abs(np.round(x, 12) - np.round(y, 12))) <= TOL_R)


def classify(live: List[Any], stored: List[Any], went_on: int = 0) -> str:
    """``live`` and ``stored`` are fit rows ``[pos, j, kind, hit, cost, evals, R]``; ``went_on`` is
    the number of stall checks of the live fit that did not stop it (one counted evaluation each).
    """
    if list(live[:3]) != list(stored[:3]):
        return "bad_structure"
    eq_h = abs(live[3] - stored[3]) <= TOL_X
    eq_c = abs(live[4] - stored[4]) <= TOL_X
    eq_r = same_R(live[6], stored[6])
    if eq_h and eq_c and eq_r:
        return "same" if int(live[5]) == int(stored[5]) + went_on else "after_last"
    return "moved"


def allowed_classes(trace: List[List[float]], stored_fit: List[Any], k: int) -> List[str]:
    """The classes the stored trace allows for a live MC fit under ``PH(k)``. No stall: ``same``.
    One stall, at the last gain: its best is G's own end, so ``after_last`` if G's end is at 0.95 or
    more (``same`` when the stop saves nothing), else ``same``. A stall before a later gain: any
    class (the hit ratio of that best is not on record)."""
    ev = SR.stall_events(trace, int(stored_fit[5]), k)
    if not ev:
        return ["same"]
    if len(ev) == 1 and ev[0]["last"]:
        return ["same", "after_last"] if float(stored_fit[3]) >= SR.TAU else ["same"]
    return ["same", "after_last", "moved"]


def stored_taken(fits: List[Dict[str, Any]], e: Dict[str, Any]) -> Dict[str, Any]:
    w = SR.walk_step(fits, {})
    p, which = w["taken"], w["kept"][w["taken"]]["which"]
    row = [f for f in e["fits"] if f[0] == p and f[2] == which][0]
    return dict(pos=p, which=which, hit=float(w["kept"][p]["hit"]), R=row[6], stopped=w["stopped"])


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--k", type=int, default=SR.K)
    ap.add_argument("--shard", default="0/1")
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--voxels", default="")
    ap.add_argument("--tag", default="")
    a = ap.parse_args()
    s = RC.SAMPLE
    assert s in ("big", "grad", "mid"), "the replay is for the disclosed samples"
    si, sn = (int(x) for x in a.shard.split("/"))
    d = SW.load(s)
    sts = SW.steps_with_mc(d)
    rlog: Dict[int, Dict[str, Any]] = {}
    for r in RC.read_gz(SC.SCACHE / f"G_{s}" / "rescue_log.json.gz"):
        if r["attempted"]:
            assert r["voxel"] not in rlog
            rlog[int(r["voxel"])] = r
    mine = [st for i, st in enumerate(sts) if i % sn == si]
    if a.voxels:
        want = [int(x) for x in a.voxels.split(",")]
        mine = [st for st in sts if st["voxel"] in want]
    if a.limit:
        mine = mine[: a.limit]

    env = RC.make_env()
    rec = env.hooks.rec
    T = env.T
    from icenine import ray_indexing as RI

    sym = RI.proper_rotations(env.setup.exp_setup.get_sample_symmetry())
    SRUN.STALLHIT["k"] = a.k
    pick = SRUN.StallHitPick(env.bfs, env.mic, rec.local_optimization, sym, use_tries=False)
    ORUN.ORD["order"] = "A"
    ratios = env.bfs._ratios
    rows: List[Dict[str, Any]] = []
    t0 = time.time()
    for st in mine:
        v = int(st["voxel"])
        e = st["entry"]
        verts = env.hooks.vertices(v)
        ph = int(env.mic.voxels[v].phase)
        pick.cur, pick.reseed = v, bool(e["reseed"])
        n_log = len(pick.e1_log)
        with contextlib.redirect_stdout(io.StringIO()):
            out = pick.e1_step(verts, ph, n_cost=20, n_refine=5)
        ne = pick.e1_log[n_log]
        live_tr = {int(x["j"]): x for x in ORUN.TRACES}
        ORUN.TRACES.clear()
        checks_by_j = {j: x["checks"] for j, x in live_tr.items()}
        classes: List[str] = []
        pred_mismatch = 0
        first_moved: Optional[int] = None
        bad = False
        for i, (lf, sf) in enumerate(zip(ne["fits"], e["fits"])):
            if first_moved is not None:
                break
            went_on = int(live_tr[int(lf[1])]["continued"]) if lf[2] == "M" else 0
            c = classify(lf, sf, went_on)
            classes.append(c)
            if c == "bad_structure":
                bad = True
            if c == "moved":
                first_moved = i
            if sf[2] == "M":
                if c not in allowed_classes(st["traces"][sf[0]]["trace"], sf, a.k):
                    pred_mismatch += 1
            elif c != "same":
                bad = True  # a CMA fit never changes
        same_len = len(ne["fits"]) == len(e["fits"])
        st_taken = stored_taken(st["fits"], e)
        n_taken_pos = out.refined.index(out.taken)
        new_taken_R = np.asarray(out.orientation, np.float64)
        new_hit = float(ne["fin_hit"][n_taken_pos])
        path_same = first_moved is None and same_len and not bad
        g_row = rlog.get(v)
        g_failed = st_taken["hit"] < SR.SEED_TEST
        g_fb_evals = 0 if g_row is None else int(g_row["rescue_evals"] + g_row["grid_evals"])
        g_total = int(e["evals"]) + g_fb_evals
        err_g: Optional[float]
        if not g_failed:
            err_g = float(T.err(np.asarray(st_taken["R"], np.float64).reshape(3, 3)[None], [v])[0])
        elif g_row is not None and g_row["rescued"]:
            Rg = np.asarray(g_row["tries"][-1]["R"], np.float64).reshape(3, 3)
            err_g = float(T.err(Rg[None], [v])[0])
        else:
            npz = SC.SCACHE / "levels" / f"{s}_{v}.npz"
            err_g = None
            if npz.exists():
                Rg = np.asarray(np.load(npz)["final_R"], np.float64)
                err_g = float(T.err(Rg[None], [v])[0])
        row: Dict[str, Any] = dict(
            sample=s,
            step=int(st["step"]),
            voxel=v,
            classes=classes,
            checks=[[int(j), c] for j, c in sorted(checks_by_j.items()) if c],
            pred_mismatch=pred_mismatch,
            bad=bool(bad),
            path_same=bool(path_same),
            g_evals_e1=int(e["evals"]),
            new_evals_e1=int(ne["evals"]),
            g_fallback_evals=g_fb_evals,
            g_total=g_total,
            g_stopped=bool(st_taken["stopped"]),
            new_stopped=bool(ne["stop"] >= 0),
            g_hit=st_taken["hit"],
            new_hit=new_hit,
            g_failed=bool(g_failed),
            err_g=err_g,
        )
        taken_same = same_R(new_taken_R, st_taken["R"])
        row["taken_same"] = bool(taken_same)
        row["mis_to_g_deg"] = (
            0.0
            if taken_same
            else float(
                env.mis(
                    new_taken_R[None], np.asarray(st_taken["R"], np.float64).reshape(3, 3)[None]
                )[0]
            )
        )
        if path_same:
            row.update(
                path="same",
                new_total=int(ne["evals"]) + g_fb_evals,
                new_failed=bool(g_failed),
                err_new=err_g,
                new_fallback_evals=g_fb_evals,
            )
        else:
            fb = 0
            new_failed = new_hit < SR.SEED_TEST
            path = "seed"
            R_final = new_taken_R
            if new_failed:
                kept = [
                    dict(R=R, cost=float(c)) for R, c in zip(out.finisher_R, out.finisher_costs)
                ]

                def one(k: int) -> Dict[str, Any]:
                    r = RCM.cma_rescue(rec, pick.orig_lo, ratios, kept[k]["R"], verts, ph)
                    return dict(hit=r["hit"], cost=r["cost"], evals=r["evals"], R=r["R"])

                res = RR.rescue_choose([k["cost"] for k in kept], one)
                fb = int(res["evals"])
                if res["rescued"]:
                    path = "rescued"
                    R_final = np.asarray(res["tries"][res["taken"]]["R"], np.float64)
                else:
                    path = "grid"
                    with contextlib.redirect_stdout(io.StringIO()):
                        gres = rec.reconstruct_voxel(
                            voxel_vertices=verts,
                            phase_index=ph,
                            rng=OD.make_rng(OD.voxel_key(verts), "grid", 0),
                        )
                    fb += int(rec.last_eval_counts[2])
                    R_final = np.asarray(gres.orientation, np.float64)
            row.update(
                path=path,
                new_total=int(ne["evals"]) + fb,
                new_failed=bool(new_failed),
                err_new=float(T.err(R_final[None], [v])[0]),
                new_fallback_evals=fb,
                final_hit=float(ratios(rec.evaluate_overlap(R_final, verts, ph))[0]),
            )
        rows.append(row)
        print(
            f"{s} step {st['step']} voxel {v}: {classes} path={row['path']} "
            f"evals {g_total}->{row['new_total']} {time.time() - t0:.0f}s",
            flush=True,
        )
    name = f"replay_{s}" + (f"_s{si}of{sn}" if sn > 1 else "") + a.tag
    (SC.SH / f"{name}.json").write_text(
        json.dumps(
            dict(
                sample=s,
                k=a.k,
                rule=SR.RULE,
                shard=a.shard,
                n_steps_replayed=len(rows),
                wall_s=time.time() - t0,
                wall_label="contended; not a timing result",
                rows=rows,
            ),
            separators=(",", ":"),
        )
        + "\n"
    )
    os._exit(0)


if __name__ == "__main__":
    main()
