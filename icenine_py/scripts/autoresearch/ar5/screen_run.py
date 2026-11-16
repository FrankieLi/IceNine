"""AR5 step 21: arms ``J5C`` and ``J5S`` on one stored map before the pass.

  RESEED_SAMPLE=<sample> SCREEN_CACHE=<dir> uv run python scripts/autoresearch/ar5/screen_run.py \
      --arm J5C [--source stored|full] | --arm J5S  [--limit N]

``J5C`` (step 20's arm): J5's stored pass tries (``--source stored``: ``D4_J5/tries.npz``, ends
orthonormalized, as step 20 read them; ``--source full``: the gate's full-budget re-run,
``full_<sample>/tries.npz``), the state after the pass, step 19's contested flag and step 20's
repair. ``J5S``: J5's map before the pass, the pass again with screen tries (``screen_pass``), the
state after it, the wide flag and the same repair. Neither reads the truth: the neighbours come
from the sample's grid. Writes ``<cache>/<arm>_<sample>/``: ``after.npy`` (the map after the
repair), ``fits.npz`` (one row per fit), ``repair.json`` (counts: the pass's evaluations, the
repair's, flagged voxels split into contested and extra, the extra voxels listed, fits, changed
voxels, wall time with the label "contended", not a timing result); J5S also ``tries.npz`` and
``after_pass.npy``, J5C ``after_pass.npy``. ``--limit N`` takes the first N flagged voxels (J5C) or
boundary voxels (J5S), a smoke test (folder suffix ``_smoke``). The freeze check runs first for
every sample (the stop rule: nothing runs before ``frozen.json`` is committed).
"""

import argparse
import json
import os
import signal
import time
from typing import Any, Dict

import screen_common as SC  # noqa: E402  (first: guards, freeze for the held samples, cache)

SC.SF.check_or_exit()  # nothing runs before the freeze is committed, for every sample
import numpy as np  # noqa: E402
import screen_data as SD  # noqa: E402
import screen_pass as SP  # noqa: E402
import screen_rules as SR  # noqa: E402


def j5c_state(d: Dict[str, Any], source: str) -> Dict[str, Any]:
    """J5C: the state after the pass and the contested flag from J5's tries (stored or re-run)."""
    y = d["sym"]
    if source == "stored":
        st = d["stored"]
        t_voxel, t_end = st["voxel"], SD.orthonormal(st["end"])
        t_cost, t_hit, evals = st["cost"], st["hit"], st["evals"]
    else:
        z = np.load(SC.SG.guarded(SC.RCACHE / f"full_{SC.SAMPLE}" / "tries.npz"))
        t_voxel, t_end = z["voxel"].astype(np.int64), z["end"].astype(np.float64)
        t_cost, t_hit, evals = z["cost"], z["hit"], z["evals"].astype(np.int64)
    R, cost, hit = SR.after_pass_state(d["R_before"], d["old_score"], t_voxel, t_end, t_cost, t_hit)
    other = SR.other_end_hit(R, d["R_before"], d["old_score"][:, 1], t_voxel, t_end, t_hit, y)
    flag, ratio, margin = SR.contested(hit, other, d["nbrs"])
    flagged = [int(v) for v in np.where(flag)[0]]
    cands = SR.repair_candidates(R, d["nbrs"], y, flagged)
    return dict(
        R=R,
        cost=cost,
        hit=hit,
        flagged=flagged,
        cands=cands,
        extra=[],
        n_contested=len(flagged),
        pass_evals=int(np.sum(evals)),
        n_tries=int(len(t_voxel)),
    )


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--arm", required=True, choices=("J5C", "J5S"))
    ap.add_argument("--source", default="stored", choices=("stored", "full"))
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--cap-min", type=float, default=0.0, help="stop rule: exit 3 at the cap")
    a = ap.parse_args()
    s = SC.SAMPLE
    if not s:
        raise SystemExit("screen_run.py: set RESEED_SAMPLE")
    if a.arm == "J5S" and a.source != "stored":
        raise SystemExit("--source is for J5C only")
    out = SC.RCACHE / (f"{a.arm}_{s}" + ("_smoke" if a.limit else ""))
    if (out / "repair.json").exists():
        raise SystemExit(f"refusing: {out / 'repair.json'} exists")
    out.mkdir(parents=True, exist_ok=True)
    if a.cap_min > 0:  # a run at its cap writes capped.json and exits 3 (the stop rule)

        def _cap(_sig: int, _frm: Any) -> None:
            (out / "capped.json").write_text(json.dumps(dict(cap_min=a.cap_min)) + "\n")
            os._exit(3)

        signal.signal(signal.SIGALRM, _cap)
        signal.alarm(int(a.cap_min * 60))
    d = SD.load_inputs(s, truth=False)
    y = d["sym"]
    import boundary_common as B  # noqa: E402  (after screen_common: the sample is set)

    env = B.Env("grad")
    assert env.n == len(d["R_before"]), (env.n, len(d["R_before"]))
    for v in range(env.n):  # the reconstructor's neighbour relation is the one the flag uses
        assert sorted(int(j) for j in env.nbr[v]) == d["nbrs"][v], v
    t_pass0 = time.time()
    extra: Dict[str, Any] = {}
    if a.arm == "J5C":
        stt = j5c_state(d, a.source)
        flagged = stt["flagged"][: a.limit] if a.limit else stt["flagged"]
        R, cost, hit, cands = stt["R"], stt["cost"], stt["hit"], stt["cands"]
        pass_evals, n_tries, n_contested, extra_vox = (
            stt["pass_evals"],
            stt["n_tries"],
            stt["n_contested"],
            [],
        )
    else:
        res = SP.run_pass(env, d["R_before"], SR.SCREEN_MAX_EVALS, limit=a.limit, log=print)
        R, cost, hit = SR.after_pass_state(
            d["R_before"], d["old_score"], res["voxel"], res["end"], res["cost"], res["hit"]
        )
        fl = SR.screen_flags(
            R,
            cost,
            hit,
            d["old_score"][:, 1],
            d["R_before"],
            res["voxel"],
            res["end"],
            res["hit"],
            d["nbrs"],
            y,
        )
        flagged, cands = fl["flagged"], fl["cands"]
        pass_evals = SR.pass_evals_of(res["evals"])
        n_tries = int(len(res["voxel"]))
        n_contested = int(fl["contested"].sum())
        extra_vox = [int(v) for v in np.where(fl["extra"])[0]]
        extra = dict(contested_mask=fl["contested"], wide_mask=fl["wide"])
        SP.save_pass(
            out,
            res,
            dict(
                sample=s,
                mode="screen",
                budget=SR.SCREEN_MAX_EVALS,
                limit=a.limit,
                n_tries=n_tries,
                evals=pass_evals,
                wall_s=res["wall_s"],
                timing_label="contended (not a timing result)",
            ),
        )
    wall_pass = time.time() - t_pass0
    np.save(out / "after_pass.npy", R)
    rec = env.hooks.rec
    kw = SR.fit_kwargs()

    def fit(v: int, R0: np.ndarray, rng: np.random.Generator) -> Any:
        r = rec.local_optimization(
            voxel_vertices=env.hooks.vertices(v),
            phase_index=env.mic.voxels[v].phase,
            initial_orientation=np.asarray(R0, dtype=np.float64),
            rng=rng,
            **kw,
        )
        h, _conf = env.bfs._ratios(r.overlap_info)
        return (
            np.asarray(r.orientation, dtype=np.float64),
            float(r.cost),
            float(h),
            int(rec.last_local_optimization_evals),
        )

    def score(v: int, Rv: np.ndarray) -> Any:
        c, h, _f = env.score3(Rv, v)
        return c, h

    t0 = time.time()
    rep = SR.repair(R, cost, hit, flagged, cands, fit, score=score)
    wall = time.time() - t0
    np.save(out / "after.npy", rep.R)
    if extra:
        np.savez_compressed(out / "flags.npz", **extra)
    np.savez_compressed(
        out / "fits.npz",
        voxel=np.array([r["voxel"] for r in rep.fits], np.int32),
        start=np.array([-1 if r["start"] is None else r["start"] for r in rep.fits], np.int32),
        end=np.array([r["end"] for r in rep.fits], np.float64).reshape(-1, 3, 3),
        cost=np.array([r["cost"] for r in rep.fits]),
        hit=np.array([r["hit"] for r in rep.fits]),
        evals=np.array([r["evals"] for r in rep.fits], np.int32),
    )
    ev = np.array([r["evals"] for r in rep.fits])
    stored_pass_evals = int(d["stored"]["evals"].sum())
    info = dict(
        sample=s,
        arm=a.arm,
        source=a.source if a.arm == "J5C" else "screen",
        limit=a.limit,
        n_tries=n_tries,
        pass_evals=pass_evals,
        stored_pass_evals=stored_pass_evals,
        n_flagged=len(rep.flagged),
        n_contested=n_contested,
        n_extra=len(extra_vox),
        extra_voxels=extra_vox,
        n_fits=len(rep.fits),
        n_changed=len(rep.changed),
        changed=rep.changed,
        n_evals=int(rep.n_evals),
        n_evals_fits=int(rep.n_evals_fits),
        n_evals_score=int(rep.n_evals_score),
        median_evals_per_fit=float(np.median(ev)) if len(ev) else 0.0,
        wall_pass_s=wall_pass,
        wall_repair_s=wall,
        timing_label="contended (not a timing result)",
        optimizer=SR.CR.REPAIR_OPTIMIZER,
        cma_sigma0_deg=SR.CR.CMA_SIGMA0_DEG,
        cma_max_evals=SR.REPAIR_MAX_EVALS,
    )
    (out / "repair.json").write_text(json.dumps(info, indent=1) + "\n")
    print(
        "done", json.dumps({k: v for k, v in info.items() if k not in ("changed", "extra_voxels")})
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
