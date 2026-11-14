"""AR3 step 6, parts 0 and 1: replay the seed step and the remedies on the fallback voxels, on 200
control index seeds and (remedy 5) on the small-cluster reseeds of D4.

  uv run python scripts/autoresearch/ar3/seed_part01.py --workers 5 [--limit N]

Writes ``seed/cache/part01.jsonl`` (one line per task; a rerun skips tasks already there) and the
run notes ``seed/cache/part01_info.json``. ``seed_summary.py`` reads the lines and writes the slim
committed files. One thread per worker (set in ``ar3_common``).
"""

import argparse
import json
import multiprocessing as mp
import time
from typing import Any, Dict, List

import ar3_common as C  # noqa: F401  (sets the thread counts before numpy / torch load)
import boundary_common as B
import numpy as np
import seed_common as S

ENV: Dict[str, Any] = {}


def init() -> None:
    import seed_replay as SR

    t0 = time.time()
    env = B.Env("grad")
    n = env.n
    maps = {m: B.load_map(m, n) for m in ("D4", "D", "P4R4", "PR")}
    ENV.update(env=env, rp=SR.Replayer(env), maps=maps, n=n)
    ENV["build_s"] = time.time() - t0


def run_map(run: str) -> str:
    return S.MAP_OF[run]


def work(task: Dict[str, Any]) -> Dict[str, Any]:
    env, rp, maps = ENV["env"], ENV["rp"], ENV["maps"]
    v = int(task["voxel"])
    t0 = time.time()
    if task["kind"] == "small":
        Rm = maps["D4"]
        t = rp.tries(v, Rm)
        end = Rm[v]
        d = (
            float(env.mis(np.asarray(t["R"], np.float64)[None], end[None])[0])
            if t["R"] is not None
            else None
        )
        return dict(
            kind="small",
            voxel=v,
            cluster=task["cluster"],
            stored_hit=task["stored_hit"],
            stored_evals=task["stored_evals"],
            n_starts=t["n_starts"],
            evals=t["evals"],
            hit=t["hit"],
            err=t["err"],
            d_to_stored=d,
            stored_err=S.err1(env, end, v),
            wall_s=time.time() - t0,
        )
    stored: Dict[str, Dict[str, Any]] = {}
    use: Dict[str, np.ndarray] = {}
    for run, info in task["runs"].items():
        if run in ("D4", "D"):
            R = maps[run][v]
        else:
            R = np.asarray(info["R"], np.float64)
        stored[run] = dict(R=R, hit=float(info["hit"]))
        use[run] = maps[run_map(run)] if run not in ("D4", "D") else maps[run]
    out = rp.voxel(v, use, stored)
    out["kind"] = task["kind"]
    out["runs"] = {r: {k: i[k] for k in i if k != "R"} for r, i in task["runs"].items()}
    return out


def build_tasks() -> List[Dict[str, Any]]:
    d4 = S.bfs_run("D4")
    dd = S.bfs_run("D")
    recs = {"D4": d4["records"], "D": dd["records"]}
    paint = {"P4": S.paint_run("P4"), "P": S.paint_run("P")}
    grid_entry = {
        r: {int(e["voxel"]): e for e in paint[r]["entries"] if e.get("kind") == "grid"}
        for r in paint
    }
    by: Dict[int, Dict[str, Any]] = {}
    for f in S.all_fallbacks():
        v, run = f["voxel"], f["run"]
        if run in recs:
            info = dict(
                hit=recs[run][str(v)]["hit_ratio"],
                phase=f["phase"],
                cluster=f["cluster"],
                full_evals=f["full_evals"],
            )
        else:
            e = grid_entry[run][v]
            info = dict(
                hit=e["hit"],
                R=e["R"],
                phase="paint",
                cluster=0,
                full_evals=0,
                index_evals=f.get("index_evals", 0),
            )
        by.setdefault(v, dict(kind="fallback", voxel=v, runs={}))["runs"][run] = info
    tasks = list(by.values())
    fb_voxels = set(by)
    for c in S.draw_controls(200):
        v = c["voxel"]
        t = dict(
            kind="control",
            voxel=v,
            runs={
                "D4": dict(
                    hit=recs["D4"][str(v)]["hit_ratio"],
                    phase=c["phase"],
                    cluster=c["cluster"],
                    stored_evals=c["stored_evals"],
                )
            },
        )
        t["also_fallback"] = v in fb_voxels
        tasks.append(t)
    for x in d4["stats"]["reseeds"]:
        if x["cluster_size"] <= 2:
            v = int(x["voxel"])
            tasks.append(
                dict(
                    kind="small",
                    voxel=v,
                    cluster=int(x["cluster_size"]),
                    stored_hit=float(recs["D4"][str(v)]["hit_ratio"]),
                    stored_evals=int(x["full_search_evaluations"]),
                    fallback=bool(x["full_search_evaluations"] > S.BIG),
                )
            )
    return tasks


def key(t: Dict[str, Any]) -> str:
    return f"{t['kind']}:{t['voxel']}"


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=5)
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--kinds", default="fallback,control,small")
    a = ap.parse_args()
    S.CACHE.mkdir(parents=True, exist_ok=True)
    path = S.CACHE / ("part01_pilot.jsonl" if a.limit else "part01.jsonl")
    done = set()
    if path.exists():
        for line in path.read_text().splitlines():
            r = json.loads(line)
            done.add(f"{r['kind']}:{r['voxel']}")
    tasks = [t for t in build_tasks() if t["kind"] in a.kinds.split(",") and key(t) not in done]
    if a.limit:
        tasks = tasks[: a.limit]
    print(f"{len(tasks)} tasks, {len(done)} done", flush=True)
    t0 = time.time()
    with mp.get_context("spawn").Pool(a.workers, initializer=init) as pool:
        with path.open("a") as f:
            for i, r in enumerate(pool.imap_unordered(work, tasks), 1):
                f.write(json.dumps(r, default=float) + "\n")
                f.flush()
                if i % 10 == 0 or i == len(tasks):
                    print(f"{i}/{len(tasks)} {time.time() - t0:.0f}s", flush=True)
    (S.CACHE / "part01_info.json").write_text(
        json.dumps(dict(wall_s=time.time() - t0, workers=a.workers, n_tasks=len(tasks))) + "\n"
    )
    print("done", time.time() - t0, flush=True)


if __name__ == "__main__":
    main()
