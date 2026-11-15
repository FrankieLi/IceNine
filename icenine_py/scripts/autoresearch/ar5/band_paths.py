"""AR5 step 15, phase A step 5 (gate, live part): the changed paths of rule ``PB(k)``.

  RESEED_SAMPLE=<big|grad|mid|held> uv run python scripts/autoresearch/ar5/band_paths.py \
      --form PB1000|PB500 [--shard I/N] [--limit M] [--tag NAME]

Needs the recording replay (``band_record.py``) of the sample and ``ar5/band/frozen.json``
committed. From the recorded stall hit ratios, ``band_rules.refine_cuts`` walks every seed step of
arm ``J`` under ``PB(k)`` by arithmetic. A step with no refinement that stops before a later gain
has the path of ``J`` (the same ends, costs and hit ratios; fewer or more evaluations): it is
priced by arithmetic in ``band_gate.py``. A step with a ``moved`` refinement (its kept end differs
from ``J``'s) is run live here, with the running code (``band_pick.BandPick``, mode ``PB(k)``: 20
costed, 5 refined), and followed along the whole changed path: the seed test (0.80) on the new
taken end; if it fails, rule ``G1``'s rescue, then ``G2b`` (step 14's frozen form), then the grid
search with the draws arm ``J`` gave it (generator of (voxel, "grid", 0)). The truth (symmetry-
reduced angle) is for scoring only. Growth after the voxel is not replayed: a cost estimate and a
voxel check, not the result.

Writes ``<cache>/paths/paths_<form>_<s>_s<I>of<N><tag>.jsonl`` (one row per replayed step, appended
as it goes; the same command resumes) and ``.done``.
"""

import argparse
import contextlib
import io
import json
import os
import re
import sys
import time
from typing import Any, Dict, List, Optional

import band_common as BC  # noqa: E402  (first)

import band_pick as BP  # noqa: E402
import band_record as BREC  # noqa: E402
import band_rules as BR  # noqa: E402
import numpy as np  # noqa: E402
import order_draws as OD  # noqa: E402
import order_run as ORUN  # noqa: E402
import polish_rules as PR  # noqa: E402
import reseed_common as RC  # noqa: E402
import rescue_common as RCM  # noqa: E402
import rescue_rules as RR  # noqa: E402
import stallhit_rules as SR  # noqa: E402

N_VOXELS = 24570


def read_records(s: str) -> List[Dict[str, Any]]:
    """The recording rows of sample ``s`` (all shards, no tagged pilot files), sorted by step."""
    rx = re.compile(rf"^record_{s}_s\d+of\d+\.jsonl$")
    rows: List[Dict[str, Any]] = []
    for p in sorted((BC.BCACHE / "record").glob("record_*.jsonl")):
        if rx.match(p.name):
            for line in BC.BG.guarded(p).read_text().splitlines():
                if line.strip():
                    rows.append(json.loads(line))
    rows.sort(key=lambda r: int(r["step"]))
    return rows


def j_taken(e: Dict[str, Any]) -> Dict[str, Any]:
    """Arm ``J``'s taken end of one stored seed step (``stallhit_replay.stored_taken``)."""
    fits = [dict(pos=f[0], j=f[1], kind=f[2], hit=f[3], cost=f[4], evals=f[5]) for f in e["fits"]]
    w = SR.walk_step(fits, {})
    p = w["taken"]
    which = w["kept"][p]["which"]
    row = [f for f in e["fits"] if f[0] == p and f[2] == which][0]
    return dict(
        pos=int(p),
        which=which,
        hit=float(w["kept"][p]["hit"]),
        R=np.asarray(row[6], np.float64),
        stopped=bool(w["stopped"]),
    )


def j_fallbacks(e1: List[Dict[str, Any]], plog: List[Dict[str, Any]]) -> Dict[int, Dict[str, Any]]:
    """Step number -> the polish-log row of its fallback call (a failed seed test), matched in
    order (voxels must agree; a mismatch raises)."""
    att = [r for r in plog if r["attempted"]]
    failed = [
        si
        for si, e in enumerate(e1, start=1)
        if e["fits"] and j_taken(e)["hit"] < SR.SEED_TEST and e["fin_hit"]
    ]
    if len(att) != len(failed):
        raise SystemExit(f"fallback rows {len(att)} != failed seed tests {len(failed)}")
    out = {}
    for si, r in zip(failed, att):
        if int(r["voxel"]) != int(e1[si - 1]["voxel"]):
            raise SystemExit(f"fallback row voxel differs at step {si}")
        out[si] = r
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--form", required=True, choices=list(BR.FORMS))
    ap.add_argument("--shard", default="0/1")
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--tag", default="")
    ap.add_argument("--force-steps", default="", help="code pilot: step numbers, no recording read")
    a = ap.parse_args()
    if a.force_steps and not a.tag:
        raise SystemExit("--force-steps needs a --tag (its rows are never read by the gate)")
    s = RC.SAMPLE
    assert s in BREC.EXPECT
    k = BR.FORMS[a.form]
    si_, sn = (int(x) for x in a.shard.split("/"))
    d = BREC.load(s)
    plog = BC.load_gz(BREC.stored_dir(s) / "polish_log.json.gz")
    fb = j_fallbacks(d["e1"], plog)
    moved: List[Dict[str, Any]] = []
    if a.force_steps:
        # a code pilot: replay these steps live whatever the recording says (no recording read)
        for x in a.force_steps.split(","):
            moved.append(dict(step=int(x), voxel=int(d["e1"][int(x) - 1]["voxel"])))
    else:
        recs = read_records(s)
        if not recs:
            raise SystemExit("no recording rows: run band_record.py first")
        bad = [r["step"] for r in recs if not r["p1_ok"]]
        if bad:
            raise SystemExit(f"prediction 1 failed at steps {bad[:5]}: stop rule")
        for r in recs:
            cuts = BR.refine_cuts(r["refinements"], k)
            if any(c["kind"] == "moved" for c in cuts.values()):
                moved.append(r)
    mine = [r for i, r in enumerate(moved) if i % sn == si_]
    if a.limit:
        mine = mine[: a.limit]
    outdir = BC.BCACHE / "paths"
    outdir.mkdir(parents=True, exist_ok=True)
    name = f"paths_{a.form}_{s}_s{si_}of{sn}{a.tag}"
    part = outdir / f"{name}.jsonl"
    done = set()
    if part.exists():
        for line in part.read_text().splitlines():
            if line.strip():
                done.add(int(json.loads(line)["step"]))
    mine = [r for r in mine if int(r["step"]) not in done]
    print(f"{s} {a.form}: {len(moved)} moved steps, {len(mine)} to replay here", flush=True)

    env = RC.make_env()
    rec = env.hooks.rec
    T = env.T
    from icenine import ray_indexing as RI
    import gradient_common as G

    sym = RI.proper_rotations(env.setup.exp_setup.get_sample_symmetry())
    BP.BAND.update(k=k, low=BR.LOW, aux=())
    pick = BP.BandPick(env.bfs, env.mic, rec.local_optimization, sym, use_tries=False)
    ORUN.ORD["order"] = "A"
    ratios = env.bfs._ratios
    R_map, _ = G.read_map(BC.BG.guarded(BREC.stored_dir(s) / "recon.mic.gz"), N_VOXELS)
    g2 = BC.g2_form()
    t0 = time.time()
    with open(part, "a") as out:
        for r in mine:
            si = int(r["step"])
            v = int(r["voxel"])
            e = d["e1"][si - 1]
            assert int(e["voxel"]) == v
            verts = env.hooks.vertices(v)
            ph = int(env.mic.voxels[v].phase)
            pick.cur, pick.reseed = v, bool(e["reseed"])
            n_log = len(pick.e1_log)
            with contextlib.redirect_stdout(io.StringIO()):
                o = pick.e1_step(verts, ph, n_cost=20, n_refine=5)
            ne = pick.e1_log[n_log]
            live_tr = {int(x["j"]): x for x in ORUN.TRACES}
            ORUN.TRACES.clear()
            jt = j_taken(e)
            jrow = fb.get(si)
            j_failed = jt["hit"] < SR.SEED_TEST
            j_fb = (
                0
                if jrow is None
                else int(jrow["rescue_evals"]) + int(jrow["polish_evals"]) + int(jrow["grid_evals"])
            )
            if not j_failed:
                err_j = float(T.err(jt["R"].reshape(3, 3)[None], [v])[0])
                j_path = "seed"
            elif jrow["rescued"]:
                err_j = float(
                    T.err(np.asarray(jrow["tries"][-1]["R"], np.float64).reshape(3, 3)[None], [v])[
                        0
                    ]
                )
                j_path = "rescued"
            elif jrow["polished"]:
                Rg = np.asarray(jrow["polish_tries"][-1]["R"], np.float64).reshape(3, 3)
                err_j = float(T.err(Rg[None], [v])[0])
                j_path = "polished"
            else:
                err_j = float(T.err(np.asarray(R_map[v], np.float64)[None], [v])[0])
                j_path = "grid"
            n_taken_pos = o.refined.index(o.taken)
            new_R = np.asarray(o.orientation, np.float64)
            new_hit = float(ne["fin_hit"][n_taken_pos])
            taken_same = bool(BREC.same_R(new_R, jt["R"]))
            mis = 0.0 if taken_same else float(env.mis(new_R[None], jt["R"].reshape(3, 3)[None])[0])
            new_failed = new_hit < SR.SEED_TEST
            fb_n = 0
            path = "seed"
            R_final = new_R
            g2_acc = False
            if new_failed:
                kept = [dict(R=R, cost=float(c)) for R, c in zip(o.finisher_R, o.finisher_costs)]

                def one(j: int) -> Dict[str, Any]:
                    x = RCM.cma_rescue(rec, pick.orig_lo, ratios, kept[j]["R"], verts, ph)
                    return dict(hit=x["hit"], cost=x["cost"], evals=x["evals"], R=x["R"])

                res = RR.rescue_choose([x["cost"] for x in kept], one)
                fb_n = int(res["evals"])
                if res["rescued"]:
                    path = "rescued"
                    R_final = np.asarray(res["tries"][res["taken"]]["R"], np.float64)
                else:
                    ends = [t["R"] for t in res["tries"]]

                    def two(j: int) -> Dict[str, Any]:
                        return PR.g2_try(g2, rec, pick.orig_lo, ratios, ends[j], verts, ph)

                    pc = PR.polish_choose(len(ends), two)
                    fb_n += int(pc["evals"])
                    if pc["accepted"]:
                        path = "polished"
                        g2_acc = True
                        R_final = np.asarray(pc["tries"][pc["taken"]]["R"], np.float64)
                    else:
                        path = "grid"
                        with contextlib.redirect_stdout(io.StringIO()):
                            gres = rec.reconstruct_voxel(
                                voxel_vertices=verts,
                                phase_index=ph,
                                rng=OD.make_rng(OD.voxel_key(verts), "grid", 0),
                            )
                        fb_n += int(rec.last_eval_counts[2])
                        R_final = np.asarray(gres.orientation, np.float64)
            err_new = float(T.err(R_final.reshape(3, 3)[None], [v])[0])
            lows = [j for j, x in live_tr.items() if x["cut"] and x["stop_kind"] == "low"]
            row = dict(
                sample=s,
                form=a.form,
                step=si,
                voxel=v,
                j_path=j_path,
                path=path,
                j_stop=int(e["stop"]),
                new_stop=int(ne["stop"]),
                stopped_low=len(lows),
                j_hit=jt["hit"],
                new_hit=new_hit,
                j_failed=bool(j_failed),
                new_failed=bool(new_failed),
                taken_same=taken_same,
                mis_to_j_deg=mis,
                j_e1_evals=int(e["evals"]),
                new_e1_evals=int(ne["evals"]),
                j_fallback_evals=j_fb,
                new_fallback_evals=fb_n,
                j_total=int(e["evals"]) + j_fb,
                new_total=int(ne["evals"]) + fb_n,
                g2_accepted=bool(g2_acc),
                err_j=err_j,
                err_new=err_new,
            )
            out.write(json.dumps(row, separators=(",", ":")) + "\n")
            out.flush()
            print(
                f"{s} {a.form} step {si} voxel {v}: path={path} taken_same={taken_same} "
                f"evals {row['j_total']}->{row['new_total']} {time.time() - t0:.0f}s",
                flush=True,
            )
    (outdir / f"{name}.done").write_text(
        json.dumps(dict(steps=len(mine), wall_s=time.time() - t0)) + "\n"
    )
    os._exit(0)


if __name__ == "__main__":
    main()
