"""AR5 step 15, phase A step 3: the recording replay of the seed steps' MC refinements of arm ``J``.

  RESEED_SAMPLE=<big|grad|mid|held> uv run python scripts/autoresearch/ar5/band_record.py \
      [--shard I/N] [--limit M] [--voxels V1,V2] [--tag NAME]

For every seed step of the stored arm ``J`` (step 14, linked read only by ``band_link_cache.sh``)
that ran an MC refinement it runs the step live, with the running code of the whole-layer runs
(``band_pick.BandPick``: interleaved order, 20 costed, 5 refined) in the recording mode: the acting
rule is arm ``J``'s own ``PH(3000)``, so the step must equal the stored one, and, besides, the hit
ratio of the running best is recorded at every 500 and every 1000 stall (never acting). The seed
step is deterministic and independent of the rest of the run (every fit has a fresh generator), so
a replay in isolation is the same step (the same argument as ``stallhit_replay.py``).

Prediction 1 (checked here, per fit): hit ratio, cost, orientation and evaluations equal the stored
ones (to 1e-9; evaluations exactly), the stored gain trace equal, and the stored ``PH(3000)`` stall
checks equal (evaluation, hit ratio to 1e-9, stopped). The stored counts of the plan's stop rule
(``J`` total with the pass, MC evaluations, steps, traces) are checked first; a difference exits 1.

Writes ``<cache>/record/record_<s>_s<I>of<N><tag>.jsonl`` (one row per step, appended as it goes: a
killed run is resumed by running the same command again) and, at the end of the shard,
``record_<s>_s<I>of<N><tag>.done``. Nothing printed or written to a tracked file here contains a
stall hit ratio.
"""

import argparse
import contextlib
import io
import json
import os
import sys
import time
from typing import Any, Dict, List, Optional

import band_common as BC  # noqa: E402  (first)

import band_pick as BP  # noqa: E402
import band_rules as BR  # noqa: E402
import numpy as np  # noqa: E402
import order_run as ORUN  # noqa: E402
import reseed_common as RC  # noqa: E402

TOL = 1e-9
EXPECT: Dict[str, Dict[str, int]] = {
    # J total with the pass, MC evaluations (plan, "What the stored record says" and stop rule)
    "held": dict(total=24_457_773, mc=11_381_609),
    "mid": dict(total=11_295_158, mc=3_795_130),
    "grad": dict(total=9_175_867, mc=1_078_188),
    "big": dict(total=3_635_256, mc=424_450),
}


def stored_dir(s: str) -> Any:
    return BC.BCACHE / f"J_{s}"


def load(s: str) -> Dict[str, Any]:
    d = stored_dir(s)
    return dict(
        e1=BC.load_gz(d / "e1_log.json.gz"),
        traces=BC.load_gz(d / "mc_traces.json.gz"),
        run=BC.load_json(d / "run.json"),
        stats=BC.load_json(BC.BD.parent / "polish" / f"stats_{s}.json"),
    )


def counts(s: str, d: Dict[str, Any]) -> Dict[str, Any]:
    """The stored counts against the plan's: ``J`` total with the pass, the MC evaluations of the
    seed steps, and the number of steps and traces."""
    mc = int(sum(int(f[5]) for e in d["e1"] for f in e["fits"] if f[2] == "M"))
    got = dict(
        total=int(d["stats"]["arms"]["J"]["stages"]["total"]),
        mc=mc,
        steps=len(d["e1"]),
        traces=len(d["traces"]),
        mc_fits=int(sum(1 for e in d["e1"] for f in e["fits"] if f[2] == "M")),
    )
    exp = EXPECT[s]
    ok = got["total"] == exp["total"] and got["mc"] == exp["mc"] and got["traces"] == got["mc_fits"]
    return dict(sample=s, expected=exp, got=got, ok=bool(ok))


def same_R(a: Any, b: Any) -> bool:
    x = np.asarray(a, np.float64).ravel()
    y = np.asarray(b, np.float64).ravel()
    return bool(np.max(np.abs(np.round(x, 12) - np.round(y, 12))) <= TOL)


def fit_diff(live: List[Any], stored: List[Any]) -> Dict[str, Any]:
    """Compare one fit row ``[pos, j, kind, hit, cost, evals, R]`` with its stored one."""
    same_struct = list(live[:3]) == list(stored[:3])
    return dict(
        struct=bool(same_struct),
        d_hit=abs(float(live[3]) - float(stored[3])),
        d_cost=abs(float(live[4]) - float(stored[4])),
        evals_equal=int(live[5]) == int(stored[5]),
        R_equal=same_R(live[6], stored[6]),
    )


def checks_equal(live: List[List[Any]], stored: List[List[Any]]) -> bool:
    if len(live) != len(stored):
        return False
    for a, b in zip(live, stored):
        if (
            int(a[0]) != int(b[0])
            or abs(float(a[1]) - float(b[1])) > TOL
            or bool(a[2]) != bool(b[2])
        ):
            return False
    return True


def trace_equal(live: List[List[float]], stored: List[List[float]]) -> bool:
    if len(live) != len(stored):
        return False
    return all(
        int(a[0]) == int(b[0]) and abs(float(a[1]) - float(b[1])) <= TOL
        for a, b in zip(live, stored)
    )


def mc_steps(d: Dict[str, Any]) -> List[Dict[str, Any]]:
    """The stored seed steps with an MC fit: step number (1-based over all steps), the entry and
    the stored traces by candidate ``j``."""
    by: Dict[int, Dict[int, Dict[str, Any]]] = {}
    for t in d["traces"]:
        by.setdefault(int(t["step"]), {})[int(t["j"])] = t
    out = []
    for si, e in enumerate(d["e1"], start=1):
        if any(f[2] == "M" for f in e["fits"]):
            out.append(dict(step=si, voxel=int(e["voxel"]), entry=e, traces=by.get(si, {})))
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--shard", default="0/1")
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--voxels", default="")
    ap.add_argument("--tag", default="")
    a = ap.parse_args()
    s = RC.SAMPLE
    assert s in EXPECT, "the recording replay is for the four disclosed samples"
    si_, sn = (int(x) for x in a.shard.split("/"))
    d = load(s)
    c = counts(s, d)
    print(json.dumps(c), flush=True)
    outdir = BC.BCACHE / "record"
    outdir.mkdir(parents=True, exist_ok=True)
    if not c["ok"]:
        (outdir / f"counts_{s}.json").write_text(json.dumps(c, indent=1) + "\n")
        print("STOP: a stored count does not reproduce", file=sys.stderr)
        sys.exit(1)
    steps = mc_steps(d)
    mine = [st for i, st in enumerate(steps) if i % sn == si_]
    if a.voxels:
        want = [int(x) for x in a.voxels.split(",")]
        mine = [st for st in steps if st["voxel"] in want]
    if a.limit:
        mine = mine[: a.limit]
    name = f"record_{s}_s{si_}of{sn}{a.tag}"
    part = outdir / f"{name}.jsonl"
    done_steps = set()
    if part.exists():
        for line in part.read_text().splitlines():
            if line.strip():
                done_steps.add(int(json.loads(line)["step"]))
    mine = [st for st in mine if st["step"] not in done_steps]
    print(f"{s} shard {a.shard}: {len(mine)} steps to replay ({len(done_steps)} done)", flush=True)

    env = RC.make_env()
    rec = env.hooks.rec
    from icenine import ray_indexing as RI

    sym = RI.proper_rotations(env.setup.exp_setup.get_sample_symmetry())
    BP.BAND.update(k=BR.PH_K, low=None, aux=BR.RECORD_KS)
    pick = BP.BandPick(env.bfs, env.mic, rec.local_optimization, sym, use_tries=False)
    ORUN.ORD["order"] = "A"
    t0 = time.time()
    bad_steps = 0
    with open(part, "a") as out:
        for st in mine:
            v = int(st["voxel"])
            e = st["entry"]
            verts = env.hooks.vertices(v)
            ph = int(env.mic.voxels[v].phase)
            pick.cur, pick.reseed = v, bool(e["reseed"])
            n_log = len(pick.e1_log)
            with contextlib.redirect_stdout(io.StringIO()):
                pick.e1_step(verts, ph, n_cost=20, n_refine=5)
            ne = pick.e1_log[n_log]
            live_tr = {int(x["j"]): x for x in ORUN.TRACES}
            ORUN.TRACES.clear()
            struct_ok = len(ne["fits"]) == len(e["fits"])
            diffs: List[Dict[str, Any]] = []
            refs: List[Dict[str, Any]] = []
            for lf, sf in zip(ne["fits"], e["fits"]):
                df = fit_diff(lf, sf)
                diffs.append(df)
                if sf[2] != "M":
                    continue
                j = int(sf[1])
                lt = live_tr[j]
                stt = st["traces"][j]
                tr_ok = trace_equal(lt["trace"], stt["trace"])
                ck_ok = checks_equal(lt["checks"], stt.get("checks", []))
                final_cost = float(lt["trace"][-1][1]) if lt["trace"] else float("inf")
                refs.append(
                    dict(
                        pos=int(sf[0]),
                        j=j,
                        evals=int(sf[5]),
                        evals_live=int(lt["evals"]),
                        cont_j=int(lt["continued"]),
                        cut_j=bool(lt["cut"]),
                        final_cost=final_cost,
                        end_hit=float(lf[3]),
                        end_cost=float(lf[4]),
                        trace_ok=bool(tr_ok),
                        checks_ok=bool(ck_ok),
                        events=lt["events"],
                    )
                )
            ok = (
                struct_ok
                and all(
                    x["struct"]
                    and x["d_hit"] <= TOL
                    and x["d_cost"] <= TOL
                    and x["evals_equal"]
                    and x["R_equal"]
                    for x in diffs
                )
                and all(r["trace_ok"] and r["checks_ok"] for r in refs)
                and int(ne["evals"]) == int(e["evals"])
            )
            if not ok:
                bad_steps += 1
            row = dict(
                sample=s,
                step=int(st["step"]),
                voxel=v,
                reseed=bool(e["reseed"]),
                evals_stored=int(e["evals"]),
                evals_live=int(ne["evals"]),
                p1_ok=bool(ok),
                worst_d=float(max([x["d_hit"] for x in diffs] + [x["d_cost"] for x in diffs])),
                refinements=refs,
            )
            out.write(json.dumps(row, separators=(",", ":")) + "\n")
            out.flush()
            print(
                f"{s} step {st['step']} voxel {v}: p1_ok={ok} evals {e['evals']} "
                f"{time.time() - t0:.0f}s",
                flush=True,
            )
    (outdir / f"{name}.done").write_text(
        json.dumps(dict(steps=len(mine), bad_steps=bad_steps, wall_s=time.time() - t0)) + "\n"
    )
    print(f"done: {len(mine)} steps, {bad_steps} not equal to the stored step", flush=True)
    os._exit(1 if bad_steps else 0)


if __name__ == "__main__":
    main()
