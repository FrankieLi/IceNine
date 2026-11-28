"""AR5 step 24, phase A step 4 (gate g0): each stored seed step of step 22's live ``J5D`` run, run
again alone in the script layer from the same voxel and image, with rule ``BF``.

  RESEED_SAMPLE=mid SCREEN_CACHE=<cache> uv run python \
      scripts/autoresearch/ar5/bf_replay.py [--limit N] [--resume] [--judge]

A seed step is deterministic and independent of the rest of a run (every refinement has its own
generator), so a replay in isolation is the same step. For every stored step of
``<cache>/J5D_<s>/e1_log.json.gz`` (run order) it runs ``bf_pick.BFPick`` (rule ``BF``) and
compares (prediction 1):

* a step that ends in the ball pass equals step 23's stored replay row (``replay_on_<s>``): the
  stop, the taken end, the evaluations and the ball stage of every refinement;
* a step that falls back equals ``J5D``'s stored step: every fit (hit ratio, cost, end orientation,
  evaluations) and the stop, the taken end; its evaluations equal the ball pass's plus the
  refinements' evaluations beyond their first stages (step 23's row and the stored traces);
* the first stage of every fallback refinement equals the kept one (checked inside the pick; a
  mismatch stops the run), no check and no pause fires in the ball pass, no lost end.

One row a step goes to ``<cache>/replay_bf_<s>/rows.jsonl`` (``--resume`` skips the steps already
written); with ``--judge`` the verdict is written to ``ar5/bf/g0.json`` (and the summary to
``<cache>/replay_bf_<s>/summary.json``). The sample is disclosed; no truth is read. Wall seconds
are labelled contended.
"""

import argparse
import contextlib
import io
import json
import sys
import time
from pathlib import Path
from typing import Any, Dict, List

import bf_common as BC  # noqa: E402  (first)

import ball_pick as BP  # noqa: E402
import ball_replay as BRP  # noqa: E402
import bf_pick as BFP  # noqa: E402
import bf_rules as BRU  # noqa: E402
import defer_pick as DP  # noqa: E402
import high500_pick as HP  # noqa: E402
import high500_rules as HR  # noqa: E402
import numpy as np  # noqa: E402
import order_run as ORUN  # noqa: E402
import reseed_common as RC  # noqa: E402

RD = BC.RD
KEYS_BALL = (
    "stop",
    "evals",
    "costed",
    "tried",
    "n_cma",
    "taken",
    "refinements",
    "first_stage_all_equal",
    "lost_end",
    "gained_end",
    "new_end_pos",
    "j5d_stop",
    "j5d_evals",
    "reseed",
    "voxel",
)


def run_live(pick: Any, env: Any, e: Dict[str, Any]) -> Dict[str, Any]:
    v = int(e["voxel"])
    verts = env.hooks.vertices(v)
    ph = int(env.mic.voxels[v].phase)
    pick.cur, pick.reseed = v, bool(e["reseed"])
    n_log, n_tr, n_ball = len(pick.e1_log), len(ORUN.TRACES), len(BP.LOG)
    n_feat = len(BFP.FEATURES)
    with contextlib.redirect_stdout(io.StringIO()):
        pick.e1_step(verts, ph, n_cost=20, n_refine=5)
    ne = pick.e1_log[n_log]
    out = dict(
        entry=ne, traces=ORUN.TRACES[n_tr:], ball=BP.LOG[n_ball:], n_feat=len(BFP.FEATURES) - n_feat
    )
    del ORUN.TRACES[n_tr:]
    del BP.LOG[n_ball:]
    del pick.e1_log[n_log:]
    return out


def compare_fallback(e: Dict[str, Any], ne: Dict[str, Any]) -> Dict[str, Any]:
    """The live fallback step against J5D's stored step."""
    fits_ok = bool(
        len(ne["fits"]) == len(e["fits"])
        and all(BRP.same_fit(x, y) for x, y in zip(ne["fits"], e["fits"]))
    )
    tk_new, tk_old = BRP.taken_row(ne), BRP.taken_row(e)
    return dict(
        fits_equal=fits_ok,
        stop_equal=bool(int(ne["stop"]) == int(e["stop"])),
        taken_equal=bool(
            list(tk_new[:3]) == list(tk_old[:3])
            and abs(tk_new[3] - tk_old[3]) <= BRP.TOL
            and BRP.same_R(tk_new[6], tk_old[6])
        ),
    )


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--resume", action="store_true")
    ap.add_argument("--judge", action="store_true")
    ap.add_argument("--fallback-only", action="store_true")
    a = ap.parse_args()
    s = BC.SAMPLE
    cache = Path(BC.RCACHE)
    out = cache / f"replay_bf_{s}{'_pilot' if a.limit else ''}"
    out.mkdir(parents=True, exist_ok=True)
    rows_p = out / "rows.jsonl"
    done_idx = set()
    if rows_p.exists():
        if not a.resume:
            raise SystemExit(f"{rows_p} exists (use --resume)")
        done_idx = {int(json.loads(x)["index"]) for x in rows_p.read_text().splitlines() if x}
    entries = BC.load_gz(cache / f"J5D_{s}" / "e1_log.json.gz")
    traces = BC.load_gz(cache / f"J5D_{s}" / "mc_traces.json.gz")
    by = BRP.group_traces(entries, traces)
    cb_rows = {
        int(r["index"]): r
        for r in (
            json.loads(x)
            for x in BC.DG.guarded(cache / f"replay_on_{s}" / "rows.jsonl").read_text().splitlines()
            if x
        )
    }
    todo = [i for i in range(len(entries)) if i not in done_idx and entries[i]["fits"]]
    if a.fallback_only:  # the pilot: the steps step 23's ball pass did not end
        todo = [i for i in todo if int(cb_rows[i]["stop"]) < 0]
    if a.limit:
        todo = todo[: a.limit]
    env = RC.make_env()
    rec = env.hooks.rec
    from icenine import ray_indexing as RI

    sym = RI.proper_rotations(env.setup.exp_setup.get_sample_symmetry())
    HP.H5.update(k=HR.K, aux=())
    ORUN.ORD["order"] = "A"
    BFP.BF["on"] = True
    BP.BALL["on"] = True
    DP.DEFER["on"] = True
    pick = BFP.BFPick(env.bfs, env.mic, rec.local_optimization, sym, use_tries=False)
    t0 = time.time()
    n_done = 0
    ev_new = ev_old = 0
    with open(rows_p, "a") as f:
        for i in todo:
            e = entries[i]
            live = run_live(pick, env, e)
            ne = live["entry"]
            ball_step = int(ne["stop"]) >= 0 and not ne["bf"]["fallback"]
            fallback = bool(ne["bf"]["fallback"])
            base = BRP.step_row(i, e, by.get(i, {}), dict(live, traces=[]), True, sym)
            cb = cb_rows[i]
            row: Dict[str, Any] = dict(
                index=i,
                voxel=int(e["voxel"]),
                reseed=bool(e["reseed"]),
                kind="fallback" if fallback else ("ball" if ball_step else "none"),
                stop=int(ne["stop"]),
                evals=int(ne["evals"]),
                j5d_stop=int(e["stop"]),
                j5d_evals=int(e["evals"]),
                cb_stop=int(cb["stop"]),
                cb_evals=int(cb["evals"]),
                lost_end=BRU.lost_end(int(e["stop"]), int(ne["stop"])),
            )
            if fallback:
                stage1 = {int(b["j"]): int(b["stage1"]) for b in live["ball"]}
                expect = int(cb["evals"]) + sum(
                    int(rf["j5d_evals"]) - int(rf["stage1"])
                    for rf in cb.get("refinements", [])
                    if rf["j5d_evals"] is not None
                )
                row.update(compare_fallback(e, ne))
                row.update(expected_evals=expect, evals_equal=bool(expect == int(ne["evals"])))
                row["cb_ball_equal"] = bool(
                    int(cb["stop"]) < 0
                    and [
                        (r["j"], r["stage1"], r["proposals"], r["gains"], r["halvings"], r["end"])
                        for r in base["refinements"]
                    ][: len(cb.get("refinements", []))]
                    == [
                        (r["j"], r["stage1"], r["proposals"], r["gains"], r["halvings"], r["end"])
                        for r in cb.get("refinements", [])
                    ]
                )
                del stage1
            else:
                diff = [
                    k
                    for k in KEYS_BALL
                    if k in cb and (base.get(k) != cb.get(k) if k != "evals" else False)
                ]
                if int(ne["evals"]) != int(cb["evals"]):
                    diff.append("evals")
                if abs(base["taken_moved_deg"] - cb["taken_moved_deg"]) > 1e-9:
                    diff.append("taken_moved_deg")
                row.update(equals_step23=bool(not diff), diff_keys=diff)
            f.write(json.dumps(row) + "\n")
            f.flush()
            n_done += 1
            ev_new += row["evals"]
            ev_old += row["j5d_evals"]
            if n_done % 25 == 0 or n_done == len(todo):
                print(
                    f"{s} bf: {n_done}/{len(todo)} steps, evals {ev_new:,} vs J5D {ev_old:,}, "
                    f"{time.time() - t0:.0f}s (contended)",
                    flush=True,
                )
    summary = dict(
        sample=s,
        steps_run=n_done,
        evals=ev_new,
        j5d_evals=ev_old,
        ball_stats=BP.STATS,
        bf_stats=BFP.STATS,
        wall_s=time.time() - t0,
        wall_label="contended; not a timing result",
    )
    (out / "summary.json").write_text(json.dumps(summary, indent=1) + "\n")
    (out / "replay.done").write_text("0\n")
    if a.judge and not a.limit:
        return judge(s, rows_p, summary)
    return 0


def judge(s: str, rows_p: Path, summary: Dict[str, Any]) -> int:
    """Prediction 1 on the written rows; the verdict goes to ``ar5/bf/g0.json``."""
    rows: List[Dict[str, Any]] = [json.loads(x) for x in rows_p.read_text().splitlines() if x]
    ball = [r for r in rows if r["kind"] == "ball"]
    fb = [r for r in rows if r["kind"] == "fallback"]
    plan_index = 2_099_764  # the plan's mid row: BF index (prediction 1)
    v = dict(
        sample=s,
        steps=len(rows),
        ball_steps=len(ball),
        fallback_steps=len(fb),
        ball_equal_step23=sum(bool(r["equals_step23"]) for r in ball),
        fallback_fits_equal=sum(bool(r["fits_equal"]) for r in fb),
        fallback_stop_equal=sum(bool(r["stop_equal"]) for r in fb),
        fallback_taken_equal=sum(bool(r["taken_equal"]) for r in fb),
        fallback_evals_equal=sum(bool(r["evals_equal"]) for r in fb),
        fallback_ball_equal_step23=sum(bool(r["cb_ball_equal"]) for r in fb),
        index_stage=sum(r["evals"] for r in rows),
        index_stage_plan=plan_index,
        lost_ends=sum(bool(r["lost_end"]) for r in rows),
        checks_in_ball_pass=int(BP.STATS["checks"]),
        pauses_in_ball_pass=int(BP.STATS["pauses"]),
        first_stage_mismatch=int(BFP.STATS["stage1_mismatch"]),
        first_stage_checked=int(BFP.STATS["stage1_checked"]),
    )
    v["g0_passes"] = bool(
        v["ball_steps"] == 1453
        and v["fallback_steps"] == 173
        and v["ball_equal_step23"] == v["ball_steps"]
        and v["fallback_fits_equal"] == v["fallback_steps"]
        and v["fallback_stop_equal"] == v["fallback_steps"]
        and v["fallback_taken_equal"] == v["fallback_steps"]
        and v["fallback_evals_equal"] == v["fallback_steps"]
        and v["index_stage"] == v["index_stage_plan"]
        and v["lost_ends"] == 0
        and v["checks_in_ball_pass"] == 0
        and v["pauses_in_ball_pass"] == 0
        and v["first_stage_mismatch"] == 0
    )
    v["stats_note"] = "steps with no candidates are not replayed; wall: " + summary["wall_label"]
    (RD / "g0.json").write_text(json.dumps(v, indent=1) + "\n")
    print(json.dumps(v, indent=1))
    return 0 if v["g0_passes"] else 1


if __name__ == "__main__":
    sys.exit(main())
