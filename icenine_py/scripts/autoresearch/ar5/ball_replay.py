"""AR5 step 23, phase A steps 4 and 5 (gates g0 and g1): each stored seed step of step 22's live
``J5D`` run, run again alone in the script layer from the same voxel and image.

  RESEED_SAMPLE=<sample> SCREEN_CACHE=<cache> uv run python \
      scripts/autoresearch/ar5/ball_replay.py --mode on|off [--limit N] [--resume]

A seed step is deterministic and independent of the rest of a run (every refinement has its own
generator; ``defer_g0.py`` relies on the same fact), so a replay in isolation is the same step. For
every stored step of ``<cache>/J5D_<s>/e1_log.json.gz`` (run order) it runs the live step with the
running code of the whole-layer run (``ball_pick.BallPick``):

* ``--mode off``: rule ``CB`` off. The step is ``J5D``'s: every fit equals the stored one (hit
  ratio, cost, end orientation, evaluations), the stop position and the step's evaluations equal
  (the first half of g0, on ``mid``).
* ``--mode on``: rule ``CB`` on. The first stage of every refinement equals ``J5D``'s stored one,
  evaluation by evaluation (the gains of the trace up to the end of the first stage; the second
  half of g0), and the step's evaluations, stop, candidates tried, taken end, and each
  refinement's second stage (proposals, gains, halvings, end) are written.

One row a step is appended to ``<cache>/replay_<mode>_<s>/rows.jsonl`` (``--resume`` skips the steps
already written); at the end ``summary.json`` and ``replay.done`` (exit code). A stall check or
pause under ``CB`` raises ``BallCheckFired`` (stop rule). No truth is read; the samples are
disclosed (``held6`` included). Wall seconds are labelled contended.
"""

import argparse
import contextlib
import io
import json
import sys
import time
from pathlib import Path
from typing import Any, Dict, List

import ball_common as BC  # noqa: E402  (first)

import ball_pick as BP  # noqa: E402
import ball_rules as BR  # noqa: E402
import defer_pick as DP  # noqa: E402
import high500_pick as HP  # noqa: E402
import high500_rules as HR  # noqa: E402
import numpy as np  # noqa: E402
import order_run as ORUN  # noqa: E402
import reseed_common as RC  # noqa: E402

TOL = 1e-12
TOL_R = 1e-9
MOVED_DEG = 0.05


def same_R(a: Any, b: Any) -> bool:
    x = np.asarray(a, np.float64).ravel()
    y = np.asarray(b, np.float64).ravel()
    return bool(np.max(np.abs(np.round(x, 12) - np.round(y, 12))) <= TOL_R)


def same_fit(live: List[Any], stored: List[Any]) -> bool:
    """``[pos, j, kind, hit, cost, evals, R]``: equal hit ratio and cost to 1e-12, the same end to
    1e-9 and the same evaluations."""
    return bool(
        list(live[:3]) == list(stored[:3])
        and abs(live[3] - stored[3]) <= TOL
        and abs(live[4] - stored[4]) <= TOL
        and same_R(live[6], stored[6])
        and int(live[5]) == int(stored[5])
    )


def taken_row(entry: Dict[str, Any]) -> Any:
    """The taken end of a logged step: the lowest-cost fit row of the stopped position, or the
    lowest-cost row when no end reaches 0.95."""
    fits = entry["fits"]
    stop = entry["stop"]
    if stop >= 0:
        return min([f for f in fits if f[0] == stop], key=lambda f: f[4])
    return min(fits, key=lambda f: f[4])


def misorientation_deg(R1: Any, R2: Any, sym: Any) -> float:
    """The smallest angle (degrees) between ``R1`` and ``R2``, over the crystal's rotations."""
    a = np.asarray(R1, np.float64).reshape(3, 3)
    b = np.asarray(R2, np.float64).reshape(3, 3)
    best = 180.0
    for o in sym:
        m = a @ (np.asarray(o, np.float64) @ b).T
        c = np.clip((np.trace(m) - 1.0) / 2.0, -1.0, 1.0)
        best = min(best, float(np.degrees(np.arccos(c))))
    return best


def group_traces(
    entries: List[Dict[str, Any]], traces: List[Dict[str, Any]]
) -> Dict[int, Dict[int, Dict[str, Any]]]:
    """The stored MC traces by step index (0-based) and candidate ``j``: the row's ``step`` is the
    run's 1-based seed-step counter. A step and a candidate has one row."""
    by: Dict[int, Dict[int, Dict[str, Any]]] = {}
    for t in traces:
        i = int(t["step"]) - 1
        assert 0 <= i < len(entries) and int(entries[i]["voxel"]) == int(t["voxel"]), (
            i,
            t["voxel"],
        )
        assert int(t["j"]) not in by.setdefault(i, {}), "two stored rows for one refinement"
        by[i][int(t["j"])] = t
    return by


def run_live(pick: Any, env: Any, e: Dict[str, Any]) -> Dict[str, Any]:
    v = int(e["voxel"])
    verts = env.hooks.vertices(v)
    ph = int(env.mic.voxels[v].phase)
    pick.cur, pick.reseed = v, bool(e["reseed"])
    n_log, n_tr, n_ball = len(pick.e1_log), len(ORUN.TRACES), len(BP.LOG)
    with contextlib.redirect_stdout(io.StringIO()):
        pick.e1_step(verts, ph, n_cost=20, n_refine=5)
    ne = pick.e1_log[n_log]
    tr = ORUN.TRACES[n_tr:]
    ball = BP.LOG[n_ball:]
    del ORUN.TRACES[n_tr:]
    del BP.LOG[n_ball:]
    del pick.e1_log[n_log:]
    return dict(entry=ne, traces=tr, ball=ball)


def step_row(
    i: int,
    e: Dict[str, Any],
    stored_tr: Dict[int, Dict[str, Any]],
    live: Dict[str, Any],
    on: bool,
    sym: Any,
) -> Dict[str, Any]:
    ne = live["entry"]
    tk_new, tk_old = taken_row(ne), taken_row(e)
    row: Dict[str, Any] = dict(
        index=i,
        voxel=int(e["voxel"]),
        reseed=bool(e["reseed"]),
        j5d_stop=int(e["stop"]),
        stop=int(ne["stop"]),
        j5d_evals=int(e["evals"]),
        evals=int(ne["evals"]),
        costed=int(ne["costed"]),
        j5d_tried=len({f[0] for f in e["fits"]}),
        tried=len({f[0] for f in ne["fits"]}),
        n_cma=sum(1 for f in ne["fits"] if f[2] == "C"),
        j5d_paused=len(e["defer"]["paused"]),
        taken_moved_deg=misorientation_deg(tk_new[6], tk_old[6], sym),
        lost_end=BR.lost_end(int(e["stop"]), int(ne["stop"])),
        gained_end=bool(int(e["stop"]) < 0 and int(ne["stop"]) >= 0),
        new_end_pos=int(ne["stop"]),
    )
    if not on:
        row["fits_equal"] = bool(
            len(ne["fits"]) == len(e["fits"])
            and all(same_fit(x, y) for x, y in zip(ne["fits"], e["fits"]))
        )
        row["pauses"] = len(ne["defer"]["paused"])
        return row
    refs: List[Dict[str, Any]] = []
    all_equal = True
    n_cmp = 0
    for b in live["ball"]:
        t = stored_tr.get(int(b["j"]))
        eq = None
        if t is not None:
            eq = BR.first_stage_equal(int(b["stage1"]), b["trace"], t["trace"])
            all_equal &= eq
            n_cmp += 1
        sp = BR.split_refinement(int(b["evals"]), int(b["stage1"]), int(b["proposals"]))
        refs.append(
            dict(
                j=int(b["j"]),
                evals=int(b["evals"]),
                stage1=sp["stage1"],
                stage2=sp["stage2"],
                final=sp["final"],
                proposals=int(b["proposals"]),
                gains=int(b["gains"]),
                halvings=int(b["halvings"]),
                end=b["end"],
                r_final_deg=float(b["r_final_deg"]),
                first_stage_equal=eq,
                j5d_evals=None if t is None else int(t["evals"]),
            )
        )
    row.update(
        refinements=refs,
        first_stage_all_equal=bool(all_equal),
        n_first_stage_compared=n_cmp,
        checks=int(sum(len(t["checks"]) for t in live["traces"])),
        j5d_taken=[int(tk_old[0]), tk_old[2], float(tk_old[3])],
        taken=[int(tk_new[0]), tk_new[2], float(tk_new[3])],
    )
    return row


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--mode", required=True, choices=("on", "off"))
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--resume", action="store_true")
    a = ap.parse_args()
    s = BC.SAMPLE
    on = a.mode == "on"
    cache = BC.RCACHE
    out = Path(cache) / f"replay_{a.mode}_{s}{'_pilot' if a.limit else ''}"
    out.mkdir(parents=True, exist_ok=True)
    rows_p = out / "rows.jsonl"
    done_idx = set()
    if rows_p.exists():
        if not a.resume:
            raise SystemExit(f"{rows_p} exists (use --resume)")
        done_idx = {int(json.loads(x)["index"]) for x in rows_p.read_text().splitlines() if x}
    entries = BC.load_gz(cache / f"J5D_{s}" / "e1_log.json.gz")
    traces = BC.load_gz(cache / f"J5D_{s}" / "mc_traces.json.gz")
    by = group_traces(entries, traces)
    todo = [i for i in range(len(entries)) if i not in done_idx and entries[i]["fits"]]
    if a.limit:
        todo = todo[: a.limit]
    env = RC.make_env()
    rec = env.hooks.rec
    from icenine import ray_indexing as RI

    sym = RI.proper_rotations(env.setup.exp_setup.get_sample_symmetry())
    HP.H5.update(k=HR.K, aux=())
    ORUN.ORD["order"] = "A"
    BP.BALL["on"] = on
    DP.DEFER["on"] = True
    pick = BP.BallPick(env.bfs, env.mic, rec.local_optimization, sym, use_tries=False)
    t0 = time.time()
    n_done = 0
    ev_new = ev_old = 0
    with open(rows_p, "a") as f:
        for i in todo:
            e = entries[i]
            live = run_live(pick, env, e)
            row = step_row(i, e, by.get(i, {}), live, on, sym)
            f.write(json.dumps(row) + "\n")
            f.flush()
            n_done += 1
            ev_new += row["evals"]
            ev_old += row["j5d_evals"]
            if n_done % 25 == 0 or n_done == len(todo):
                print(
                    f"{s} {a.mode}: {n_done}/{len(todo)} steps, evals {ev_new:,} vs J5D "
                    f"{ev_old:,}, {time.time() - t0:.0f}s (contended)",
                    flush=True,
                )
    (out / "summary.json").write_text(
        json.dumps(
            dict(
                sample=s,
                mode=a.mode,
                steps_run=n_done,
                evals=ev_new,
                j5d_evals=ev_old,
                stats=BP.STATS,
                wall_s=time.time() - t0,
                wall_label="contended; not a timing result",
            ),
            indent=1,
        )
        + "\n"
    )
    (out / "replay.done").write_text("0\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
