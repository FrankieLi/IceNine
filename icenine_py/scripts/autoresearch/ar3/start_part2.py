"""AR3 step 9, part 2, remedy F2 on stored cases: one sample, one set (pick or test).

  scripts/autoresearch/ar3/start_bg.sh big f2_pick_big start_part2.py --set pick
  scripts/autoresearch/ar3/start_bg.sh grad f2_test_grad start_part2.py --set test

``--set test`` refuses to run unless ``start/f_pick.json`` (the pick, written by
``start_part2_pick.py`` from the pick sets only) is committed in git: the order pick, then test, is
shown by the commit history.

F2 (the plan's table): before a rejected fit is accepted as rejected, one more local optimization
from its end (CMA, first step 0.05, 0.1 or 0.2 degrees, at most 60 evaluations plus the final
one), judged by the same relative test (hit ratio / max(flood best before, hit ratio) > 0.9). It
runs on every rejected fit of the run (the BFS does not know the class), so its evaluations on
every rejected fit count. If the new end has the higher hit ratio it replaces the rejected fit's
end (this rule was written here, before any F2 result was read).

Stored cases (the logged run of ``S5``): for every rejected fit (relative test) of a voxel in a true
grain of the set: the three settings' hit ratio at the end, whether the relative test passes, the
error at the end, the evaluations. A reseed cluster of the set (set by the true grain of the
reseed voxel) *vanishes* under a setting when the last rejected fit of every one of its voxels
(before the reseed) passes. The target of F2 is the seed-step evaluations of the reseed clusters
(``full_search_evaluations`` of the reseed list of ``run.json``); the saving is the seed steps of
the clusters that vanish, less F2's own evaluations on all rejected fits of the set. Cases: the
voxels of a vanishing cluster at their first cluster; *fixed* is a voxel whose final map entry is
wrong (error above 1 degree) and F2's accepted end is right; *broken* is a voxel whose final map
entry is right and F2's end is wrong. (The stored final entry includes the reseed and later fits;
what F2 would have changed downstream is not measured.)

Writes ``start/part2_<set>_<sample>.json``.
"""

import argparse
import json
import subprocess
from typing import Any, Dict, List

import numpy as np
import reseed_common as RC
import start_part1 as P1

import boundary_common as B  # noqa: E402
import ar3_common as C  # noqa: E402

SIGMAS = (0.05, 0.1, 0.2)
START = P1.START


def require_pick_committed() -> None:
    rel = "icenine_py/benchmarks/autoresearch/ar3/start/f_pick.json"
    root = C.ICE.parent
    tracked = subprocess.run(
        ["git", "-C", str(root), "ls-files", "--error-unmatch", rel], capture_output=True
    )
    dirty = subprocess.run(
        ["git", "-C", str(root), "status", "--porcelain", "--", rel],
        capture_output=True,
        text=True,
    ).stdout.strip()
    assert tracked.returncode == 0 and not dirty, "commit start/f_pick.json before the test set"


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--set", required=True, choices=["pick", "test"])
    a = ap.parse_args()
    if a.set == "test":
        require_pick_committed()
    want_even = a.set == "pick"
    d = P1.LOGGED[RC.SAMPLE]
    f = dict(np.load(START / f"feat_{RC.SAMPLE}.npz"))
    z = np.load(d / "fits.npz")
    Re_all = B.orthonormal(z["R"][:, 9:].astype(np.float64).reshape(-1, 3, 3))
    env = RC.make_env()
    n = len(f["voxel"])
    in_set = (f["grain"] % 2 == 0) == want_even
    best_a = np.maximum(f["hit_end"], f["best_before"])
    rejected = (best_a > 0) & (f["hit_end"] / np.where(best_a > 0, best_a, 1) <= 0.9)
    todo = np.nonzero(in_set & rejected)[0]
    print(f"{a.set}: {int(in_set.sum())} fits in the set, {len(todo)} rejected", flush=True)
    res: Dict[float, Dict[int, Dict[str, Any]]] = {s: {} for s in SIGMAS}

    def f2_row(i: int) -> None:
        """F2 on rejected fit ``i`` for the three settings (memoized)."""
        v = int(f["voxel"][i])
        for s in SIGMAS:
            if i in res[s]:
                continue
            ef = P1.extra_fit(env, v, Re_all[i], s)
            h2 = ef["hit"]
            use_new = h2 > float(f["hit_end"][i])
            hit = h2 if use_new else float(f["hit_end"][i])
            best = float(f["best_before"][i])
            ba = max(hit, best)
            res[s][i] = dict(
                evals=ef["evals"],
                hit=hit,
                passes=bool(ba > 0 and hit / ba > 0.9),
                err=(
                    float(env.T.err(ef["end"][None], np.array([v]))[0])
                    if use_new
                    else float(f["err_end"][i])
                ),
            )

    for j, i in enumerate(todo):
        f2_row(int(i))
        if j % 1000 == 0:
            print(f"  F2 fits {j}/{len(todo)}", flush=True)
    in_set_rows = set(int(i) for i in todo)
    # clusters
    clusters = RC.read_gz(d / "clusters.json.gz")
    floods = RC.read_gz(d / "floods.json.gz")
    reseed_floods = [x for x in floods if x["reseed"]]
    run = json.loads((d / "run.json").read_text())
    resl = run["stats"]["reseeds"]
    assert len(resl) == len(clusters) == len(reseed_floods)
    first_row: Dict[int, int] = {}
    for i, fl in enumerate(f["flood"]):
        first_row.setdefault(int(fl), i)
    R_final, _ = B.G.read_map(d / "recon.mic.gz", env.n)
    err_final = env.T.err(R_final, np.arange(env.n))
    gm = np.asarray(env.gmap)
    seen: set = set()
    cl_out: List[Dict[str, Any]] = []
    for k, (c, fl) in enumerate(zip(clusters, reseed_floods)):
        assert resl[k]["voxel"] == fl["voxel"]
        if (gm[fl["voxel"]] % 2 == 0) != want_even:
            for m in c["members"]:
                seen.add(int(m))
            continue
        lim = first_row.get(fl["flood"], n)
        last: Dict[int, int] = {}
        for i in range(lim):
            last[int(f["voxel"][i])] = i
        mem = []
        for m in c["members"]:
            m = int(m)
            first = m not in seen
            seen.add(m)
            i = last.get(m)
            mem.append((m, i, first))
        row: Dict[str, Any] = dict(
            k=k, seed_evals=int(resl[k]["full_search_evaluations"]), size=len(mem)
        )
        for _m, i, _fi in mem:
            if i is not None:
                f2_row(int(i))  # a member in a grain of the other set: its fit is needed too
        for s in SIGMAS:
            ok = all(i is not None and res[s][int(i)]["passes"] for _, i, _ in mem)
            row[f"vanishes_{s}"] = bool(ok)
            fixed = broken = both_wrong = 0
            if ok:
                for m, i, first in mem:
                    if not first:
                        continue
                    right_f = bool(err_final[m] <= RC.RIGHT)
                    right_2 = bool(res[s][int(i)]["err"] <= RC.RIGHT)
                    fixed += int((not right_f) and right_2)
                    broken += int(right_f and (not right_2))
                    both_wrong += int((not right_f) and (not right_2))
            row[f"fixed_{s}"], row[f"broken_{s}"], row[f"both_wrong_{s}"] = (
                fixed,
                broken,
                both_wrong,
            )
        cl_out.append(row)
    out: Dict[str, Any] = dict(
        sample=RC.SAMPLE,
        set=a.set,
        fits_in_set=int(in_set.sum()),
        rejected_fits=len(todo),
        clusters=len(cl_out),
        seed_evals_clusters=int(sum(x["seed_evals"] for x in cl_out)),
        settings={},
    )
    for s in SIGMAS:
        van = [x for x in cl_out if x[f"vanishes_{s}"]]
        own = [r for i, r in res[s].items() if i in in_set_rows]  # the set's rejected fits
        ev = int(sum(r["evals"] for r in own))
        out["settings"][str(s)] = dict(
            f2_evals=ev,
            rejected_pass=int(sum(1 for r in own if r["passes"])),
            clusters_vanish=len(van),
            seed_evals_vanish=int(sum(x["seed_evals"] for x in van)),
            fixed=int(sum(x[f"fixed_{s}"] for x in cl_out)),
            broken=int(sum(x[f"broken_{s}"] for x in cl_out)),
            both_wrong=int(sum(x[f"both_wrong_{s}"] for x in cl_out)),
            err_median_new=float(np.median([r["err"] for r in own])) if own else None,
        )
    out["err_median_stored_end"] = float(np.median(f["err_end"][todo])) if len(todo) else None
    (START / f"part2_{a.set}_{RC.SAMPLE}_f2.json").write_text(json.dumps(out, indent=1) + "\n")
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
