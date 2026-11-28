"""AR3 step 8, part 3 pilot: the E1 loop against the library's seed step and against part 2.

  RESEED_SAMPLE=big uv run python scripts/autoresearch/ar3/reseed_s8_check.py --n 12

For the first ``n`` stored S5 seed steps (spread over the voxels): (a) with a stop test no
refinement can pass (tau 2.0) the E1 loop must give the library's ``index_seed(20, 5)`` result:
same evaluations, same orientation; (b) with the chosen tau its evaluations must equal the ones
part 2 computed from the replay rows (``reseed_e1.e1_step``). Writes
``reseed/s8_check_<sample>.json``.
"""

import argparse
import json

import numpy as np
import reseed_common as RC
import reseed_e1 as E1
import reseed_s8 as S8


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=12)
    a = ap.parse_args()
    run = "S5L" if RC.SAMPLE == "big" else "S5"
    rows = E1.load(RC.SAMPLE, run)
    step = max(1, len(rows) // a.n)
    rows = rows[::step][: a.n]
    env = RC.make_env()
    rec = env.hooks.rec
    pick = S8.E1Pick(env.bfs, env.mic, rec.local_optimization, env.hooks.sym, use_tries=False)
    res = []
    for r in rows:
        v = r["voxel"]
        verts = env.hooks.vertices(v)
        ph = int(env.mic.voxels[v].phase)
        pick.cur, pick.reseed = v, bool(r["reseed"])
        lib = rec.index_seed(verts, ph, n_cost=20, n_refine=5)
        mine_all = pick.e1_step(verts, ph, n_cost=20, n_refine=5, tau=2.0)
        mine = pick.e1_step(verts, ph, n_cost=20, n_refine=5)
        want = E1.e1_step(r, S8.TAU)
        res.append(
            dict(
                voxel=v,
                lib_evals=int(lib.evals),
                all_evals=int(mine_all.evals),
                same_orientation=bool(np.array_equal(lib.orientation, mine_all.orientation)),
                same_cost=bool(lib.cost == mine_all.cost),
                e1_evals=int(mine.evals),
                part2_evals=int(want["evals_e1"]),
            )
        )
        print(res[-1], flush=True)
    ok = all(
        x["lib_evals"] == x["all_evals"]
        and x["same_orientation"]
        and x["same_cost"]
        and x["e1_evals"] == x["part2_evals"]
        for x in res
    )
    out = dict(sample=RC.SAMPLE, tau=S8.TAU, n=len(res), all_equal=ok, rows=res)
    (RC.OUT / f"s8_check_{RC.SAMPLE}.json").write_text(json.dumps(out, indent=1) + "\n")
    print("all equal:", ok)


if __name__ == "__main__":
    main()
