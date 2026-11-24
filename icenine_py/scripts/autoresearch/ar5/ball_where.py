"""AR5 step 23, information only (the step failed its gate g1): where were the 0.95 ends that rule
``CB`` loses?

  RESEED_SAMPLE=<sample> SCREEN_CACHE=<cache> uv run python scripts/autoresearch/ar5/ball_where.py

For every replayed step that loses its 0.95 end (``J5D`` had one, ``CB`` has none; from
``replay_on_<s>/rows.jsonl``) the step is run again under ``CB`` (as ``ball_replay.py``), keeping
the first stage's best orientation of every refinement (the start of the second stage; the call of
``ball_rules.ball_stage`` is wrapped here, no pinned script is touched). For the refinement that
gave ``J5D`` its end (the MC fit of the stopped position, from the stored ``e1_log`` row) the stored
trace gives the evaluation of its LAST GAIN, and with the first stage's length ``s1`` of the same
candidate in the ``CB`` run:

* ``stage1``: last gain at evaluation ``s1`` or earlier (the end was already in the first stage);
* ``stage2_first_block``: after ``s1`` and up to ``s1 + 11`` (the start evaluation and the first
  block of 10 proposals: before the first restart; a gain here continues the block's best);
* ``stage2_later``: after that (a restart, or a continuation after a gain, of the second stage).

The distance is the angle (degrees, smallest over the crystal's rotations) from that end to the
first stage's best, in bins 0-0.1, 0.1-0.33, 0.33-0.64, over 0.64. Wide restarts reach out to
about 0.32, 0.64 and 1.11 degrees; the ball starts at 0.329. Disclosed samples, for information.
Writes ``<cache>/where_<s>.json``.
"""

import json
import sys
from pathlib import Path
from typing import Any, Dict, List

import ball_common as BC  # noqa: E402  (first)

import ball_pick as BP  # noqa: E402
import ball_replay as BRP  # noqa: E402
import ball_rules as BR  # noqa: E402
import defer_pick as DP  # noqa: E402
import high500_pick as HP  # noqa: E402
import high500_rules as HR  # noqa: E402
import numpy as np  # noqa: E402
import order_run as ORUN  # noqa: E402
import reseed_common as RC  # noqa: E402

BINS = (0.1, 0.33, 0.64)
CAPTURED: List[np.ndarray] = []
_real = BR.ball_stage


def _wrapped(evaluate: Any, R0: Any, *a: Any, **k: Any) -> Any:
    CAPTURED.append(np.array(R0, dtype=np.float64, copy=True))
    return _real(evaluate, R0, *a, **k)


def bin_of(d: float) -> str:
    lo = 0.0
    for hi in BINS:
        if d < hi:
            return f"{lo:g}-{hi:g}"
        lo = hi
    return f"over {BINS[-1]:g}"


def main() -> int:
    s = BC.SAMPLE
    cache = Path(BC.RCACHE)
    rows = [
        json.loads(x)
        for x in (cache / f"replay_on_{s}" / "rows.jsonl").read_text().splitlines()
        if x
    ]
    lost = [r for r in rows if r["lost_end"]]
    entries = BC.load_gz(cache / f"J5D_{s}" / "e1_log.json.gz")
    by = BRP.group_traces(entries, BC.load_gz(cache / f"J5D_{s}" / "mc_traces.json.gz"))
    BR.ball_stage = _wrapped  # type: ignore[assignment]
    env = RC.make_env()
    rec = env.hooks.rec
    from icenine import ray_indexing as RI

    sym = RI.proper_rotations(env.setup.exp_setup.get_sample_symmetry())
    HP.H5.update(k=HR.K, aux=())
    ORUN.ORD["order"] = "A"
    BP.BALL["on"] = True
    DP.DEFER["on"] = True
    pick = BP.BallPick(env.bfs, env.mic, rec.local_optimization, sym, use_tries=False)
    out: List[Dict[str, Any]] = []
    for r in lost:
        i = int(r["index"])
        e = entries[i]
        del CAPTURED[:]
        n0 = len(BP.LOG)
        live = BRP.run_live(pick, env, e)
        assert len(CAPTURED) == len(live["ball"])
        p = int(e["stop"])
        f = next(x for x in e["fits"] if x[0] == p and x[2] == "M" and x[3] >= 0.95)
        j = int(f[1])
        t = by[i][j]
        last = int(t["trace"][-1][0])
        k = next((q for q, b in enumerate(live["ball"]) if int(b["j"]) == j), None)
        row: Dict[str, Any] = dict(
            index=i,
            voxel=int(e["voxel"]),
            reseed=bool(e["reseed"]),
            j5d_stop=p,
            j=j,
            last_gain=last,
        )
        if k is None:
            row.update(zone="candidate_not_refined_by_CB")
        else:
            s1 = int(live["ball"][k]["stage1"])
            dist = BRP.misorientation_deg(f[6], CAPTURED[k], sym)
            zone = (
                "stage1"
                if last <= s1
                else ("stage2_first_block" if last <= s1 + 11 else "stage2_later")
            )
            row.update(stage1=s1, zone=zone, dist_deg=dist, dist_bin=bin_of(dist))
        out.append(row)
        del BP.LOG[n0:]
    (cache / f"where_{s}.json").write_text(json.dumps(dict(sample=s, rows=out), indent=1) + "\n")
    print(s, "lost steps rerun:", len(out))
    return 0


if __name__ == "__main__":
    sys.exit(main())
