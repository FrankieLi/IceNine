"""AR5 step 22, phase A step 4 (gate g0): the live seed step of arm ``J5D`` against arm ``J5``'s
stored seed steps on ``mid``, in isolation.

  RESEED_SAMPLE=mid uv run python scripts/autoresearch/ar5/defer_g0.py [--per-type 8]

A seed step is deterministic and independent of the rest of a run (every refinement has its own
generator; ``stallhit_replay.py`` relies on the same fact), so a replay in isolation is the same
step. For a chosen set of ``mid``'s stored seed steps (the first ``--per-type`` of each type A, C, N
and none, in run order; ``defer_rules.walk_step`` gives the type) it runs the live step twice with
the running code of the whole-layer run (``defer_pick.DeferPick``):

* **DEF off** (nothing paused): every fit equals the stored one (hit ratio, cost, end orientation,
  evaluations): ends and evaluations are ``J5``'s (the first half of g0).
* **DEF on**: type N and none steps equal ``J5``'s stored fits (ends, evaluations: a resumed
  refinement equals ``J5``'s unpaused one, evaluation by evaluation: its trace of gains, its
  checks and its evaluations equal the stored ``mc_traces`` row's); type A steps have ``J5``'s
  taken end and ``J5``'s evaluations minus the paused tails (exactly the replay's saving); type C
  steps are listed with the old and the new taken ends. The re-run's evaluations (not counted) are
  listed apart.

Writes ``ar5/defer/g0.json`` (one row per step and mode and the verdict). Exit 1 when g0 fails.
Neither the truth nor held6 is read; the sample is the disclosed ``mid``.
"""

import argparse
import contextlib
import io
import json
import sys
import time
from typing import Any, Dict, List

import defer_common as DC  # noqa: E402  (first)

import defer_pick as DP  # noqa: E402
import defer_rules as DR  # noqa: E402
import high500_pick as HP  # noqa: E402
import high500_rules as HR  # noqa: E402
import numpy as np  # noqa: E402
import order_run as ORUN  # noqa: E402
import reseed_common as RC  # noqa: E402

TOL = 1e-12
TOL_R = 1e-9


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


def load_steps(cache: Any, s: str) -> List[Dict[str, Any]]:
    """The stored seed steps with their MC traces by candidate position and their walk."""
    entries = DC.load_gz(cache / f"J5_{s}" / "e1_log.json.gz")
    traces = iter(DC.load_gz(cache / f"J5_{s}" / "mc_traces.json.gz"))
    out: List[Dict[str, Any]] = []
    for i, x in enumerate(entries):
        by: Dict[int, Dict[str, Any]] = {}
        for f in x["fits"]:
            if f[2] == "M":
                t = next(traces)
                assert t["voxel"] == x["voxel"] and t["evals"] == f[5]
                by[int(f[0])] = t
        w = DR.walk_step(x, by)
        out.append(dict(index=i, entry=x, traces=by, walk=w))
    return out


def choose(steps: List[Dict[str, Any]], per_type: int) -> List[Dict[str, Any]]:
    got: Dict[str, List[Dict[str, Any]]] = {"A": [], "C": [], "N": [], "-": []}
    for st in steps:
        if not st["traces"]:
            continue  # no MC refinement in the step: nothing to pause
        lst = got[st["walk"]["type"]]
        if len(lst) < per_type:
            lst.append(st)
    return sorted([x for v in got.values() for x in v], key=lambda x: x["index"])


def run_live(pick: Any, env: Any, st: Dict[str, Any], on: bool) -> Dict[str, Any]:
    e = st["entry"]
    v = int(e["voxel"])
    verts = env.hooks.vertices(v)
    ph = int(env.mic.voxels[v].phase)
    pick.cur, pick.reseed = v, bool(e["reseed"])
    DP.DEFER["on"] = on
    n_log, n_tr = len(pick.e1_log), len(ORUN.TRACES)
    rerun0 = DP.STATS["rerun_evals"]
    with contextlib.redirect_stdout(io.StringIO()):
        pick.e1_step(verts, ph, n_cost=20, n_refine=5)
    ne = pick.e1_log[n_log]
    live_tr = ORUN.TRACES[n_tr:]
    del ORUN.TRACES[n_tr:]
    return dict(entry=ne, traces=live_tr, rerun_evals=DP.STATS["rerun_evals"] - rerun0)


def taken_row(entry: Dict[str, Any]) -> Any:
    """The taken end of a logged step: the first fit row of the stopped position's kept end (the
    lowest cost of its CMA and MC fits), or the lowest-cost end when no end reaches 0.95."""
    fits = entry["fits"]
    stop = entry["stop"]
    if stop >= 0:
        rows = [f for f in fits if f[0] == stop]
        return min(rows, key=lambda f: f[4])
    return min(fits, key=lambda f: f[4])


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--per-type", type=int, default=8)
    a = ap.parse_args()
    s = DC.SAMPLE
    if s != "mid":
        raise SystemExit("defer_g0.py: g0 is judged on mid (set RESEED_SAMPLE=mid)")
    steps = load_steps(DC.RCACHE, s)
    mine = choose(steps, a.per_type)
    env = RC.make_env()
    rec = env.hooks.rec
    from icenine import ray_indexing as RI

    sym = RI.proper_rotations(env.setup.exp_setup.get_sample_symmetry())
    HP.H5.update(k=HR.K, aux=())
    ORUN.ORD["order"] = "A"
    pick = DP.DeferPick(env.bfs, env.mic, rec.local_optimization, sym, use_tries=False)
    rows: List[Dict[str, Any]] = []
    t0 = time.time()
    for st in mine:
        e, w = st["entry"], st["walk"]
        row: Dict[str, Any] = dict(
            index=st["index"], voxel=int(e["voxel"]), type=w["type"], stop_j5=int(e["stop"])
        )
        # ---- DEF off: J5's step, fit by fit
        off = run_live(pick, env, st, False)
        ne = off["entry"]
        fits_equal = len(ne["fits"]) == len(e["fits"]) and all(
            same_fit(x, y) for x, y in zip(ne["fits"], e["fits"])
        )
        row["off_fits_equal"] = bool(fits_equal)
        row["off_evals_equal"] = int(ne["evals"]) == int(e["evals"])
        row["off_stop_equal"] = int(ne["stop"]) == int(e["stop"])
        row["off_pauses"] = len(ne["defer"]["paused"])
        # ---- DEF on
        on = run_live(pick, env, st, True)
        ne = on["entry"]
        d = ne["defer"]
        row.update(
            on_evals=int(ne["evals"]),
            j5_evals=int(e["evals"]),
            on_stop=int(ne["stop"]),
            paused=d["paused"],
            resumed=d["resumed"],
            found=d["found"],
            rerun_evals=int(on["rerun_evals"]),
            expected_evals=w["expected_evals"],
        )
        tk_new, tk_old = taken_row(ne), taken_row(e)
        row["taken_same"] = bool(same_R(tk_new[6], tk_old[6]))
        if w["type"] in ("A", "N", "-"):
            row["on_evals_ok"] = int(ne["evals"]) == int(w["expected_evals"])
            row["on_end_ok"] = row["taken_same"] and int(ne["stop"]) == int(e["stop"])
        else:  # type C: the end may move; listed with the old and the new taken end
            row["on_evals_ok"] = True
            row["on_end_ok"] = True
            row["c_moved"] = not (row["taken_same"] and int(ne["stop"]) == int(e["stop"]))
            row["c_old_taken"] = [int(tk_old[0]), tk_old[2], float(tk_old[3])]
            row["c_new_taken"] = [int(tk_new[0]), tk_new[2], float(tk_new[3])]
        # a resumed refinement equals J5's unpaused one, evaluation by evaluation
        res_ok = True
        n_res = 0
        by_j = {int(t["j"]): t for t in on["traces"] if t.get("resumed")}
        for p in d["resumed"]:
            st_tr = st["traces"][p]
            lv = by_j.get(int(e["fits"][[f[0] for f in e["fits"]].index(p)][1]))
            n_res += 1
            if lv is None:
                res_ok = False
                continue
            res_ok &= bool(
                lv["evals"] == st_tr["evals"]
                and lv["trace"] == st_tr["trace"]
                and [list(c) for c in lv["checks"]] == [list(c) for c in st_tr["checks"]]
            )
        row["resumed_equal_j5"] = bool(res_ok)
        row["n_resumed_checked"] = n_res
        rows.append(row)
        print(
            f"step {st['index']} voxel {row['voxel']} {w['type']}: off equal "
            f"{row['off_fits_equal']} | on evals {row['on_evals']} vs J5 {row['j5_evals']} "
            f"(expected {row['expected_evals']}) paused {d['paused']} resumed {d['resumed']} "
            f"found {d['found']} {time.time() - t0:.0f}s",
            flush=True,
        )
    by_type = {t: sum(1 for r in rows if r["type"] == t) for t in ("A", "C", "N", "-")}
    ok_off = all(r["off_fits_equal"] and r["off_evals_equal"] and r["off_stop_equal"] for r in rows)
    ok_off &= all(r["off_pauses"] == 0 for r in rows)
    ok_on = all(r["on_evals_ok"] and r["on_end_ok"] and r["resumed_equal_j5"] for r in rows)
    ok_rerun = DP.STATS["rerun_checked"] == DP.STATS["resumes"]
    verdict = bool(ok_off and ok_on and ok_rerun and rows)
    out = dict(
        sample=s,
        steps_by_type=by_type,
        n_steps=len(rows),
        def_off_equals_j5=bool(ok_off),
        def_on_a_n_none_as_expected=bool(ok_on),
        every_resume_rerun_checked=bool(ok_rerun),
        stats=DP.STATS,
        g0_passes=verdict,
        wall_s=time.time() - t0,
        wall_label="contended; not a timing result",
        rows=rows,
    )
    (DC.RD / "g0.json").write_text(json.dumps(out, indent=1) + "\n")
    print(json.dumps({k: v for k, v in out.items() if k != "rows"}, indent=1))
    return 0 if verdict else 1


if __name__ == "__main__":
    sys.exit(main())
