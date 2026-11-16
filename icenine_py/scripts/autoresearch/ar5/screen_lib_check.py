"""AR5 step 21, phase B: one library run on ``mid`` with the key ``BoundaryPassRepair 1`` and
``BoundaryPassScreen 1`` on, and a check that the library's screen pass and wide repair equal the
script rule of arm ``J5S`` (``screen_rules``) on the same inputs. As in step 20 the library run
cannot equal the script arm (its in-run pass reads the unrounded map), so the check is on this
run's own inputs, bit for bit.

  RESEED_SAMPLE=mid SCREEN_CACHE=<cache> uv run python \
      scripts/autoresearch/ar5/screen_lib_check.py [--cap-min 240]

The run is ``LibJ5`` of step 16 (``rescue_lib_run.S11L`` keys, ``SeedRescueCMA 0.25 1000``,
``SeedIndexMCStallHit 500``: the library form of the script arm ``J5l``) plus ``BoundaryPassRepair
1``. It is NOT the script arm ``J5`` (which has rule ``G2b`` and the pass with aligned draws, both
script-layer), so its map and counts cannot equal ``J5C``'s. What is compared, on this run's own
inputs, is the repair itself:

1. state: every voxel's strict score is recomputed on the map after the pass and compared with the
   library's state (voxel ``cost`` and ``overlap_ratio``, boundary voxels from the pass's re-score);
2. the script rule (``screen_rules``: the screen flag and step 20 repair, the CMA local
   optimization with the
   per-fit generator) is applied to the state rebuilt the script way (re-scored old ends, the
   pass's tries) and compared with the library's repair: flagged voxels, the order of the local
   optimizations (voxel, start), the evaluations of each, the total evaluations, the changed
   voxels and the final map (exact, byte for byte, as float64).

Writes ``<cache>/LibJ5S_mid/`` and ``benchmarks/autoresearch/ar5/screen/lib_check.json``.
Truth of ``mid`` (a disclosed sample) is used for the accuracy line only. Wall times are contended.
"""

import argparse
import json
import os
import sys
import threading
import time
from pathlib import Path
from typing import Any, Dict

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import screen_common as RCM  # noqa: E402  (first: guards, cache, stale freeze checks off)

import ar3_common as C  # noqa: E402
import screen_rules as SR  # noqa: E402
import gradient_common as G  # noqa: E402
import numpy as np  # noqa: E402
import rescue_lib_run as RLR  # noqa: E402
import reseed_common as RC  # noqa: E402
import stallhit_cache_guard as SCG  # noqa: E402

KEYS = RLR.S11L + [
    "SeedRescueCMA 0.25 1000",
    "SeedIndexMCStallHit 500",
    "BoundaryPassRepair 1",
    "BoundaryPassScreen 1",
]
NAME = "LibJ5S"
CACHE = RCM.RCACHE
OUTJ = RCM.RD / "lib_check.json"


def compare(bfs: Any, mic: Any, res: Any) -> Dict[str, Any]:
    """The checks 1 and 2 of the module text; every number is a count or a maximum."""
    from icenine.ray_indexing import misorientation_deg, proper_rotations
    from icenine.reconstructor import _get_voxel_vertices

    rep = res.repair
    n = len(mic.voxels)
    y = proper_rotations(bfs.setup.exp_setup.get_sample_symmetry())
    rec = bfs.reconstructor
    cost_fn = rec._make_local_cost_fn()
    nbrs = [
        sorted(
            int(j)
            for j in mic.get_neighbors(v, radius=2.0 * mic.voxels[v].side_length)
            if int(j) != v
        )
        for v in range(n)
    ]
    # the map before the pass: the state before the repair, with the pass's changes undone
    Rb = rep.R_pre.copy()
    for e in bfs.boundary_log:
        Rb[e["voxel"]] = e["old"]
    # 1. strict re-score of every voxel on the map after the pass
    old_score = np.zeros((n, 3))
    after_score = np.zeros((n, 3))
    for v in range(n):
        vox = mic.voxels[v]
        verts = _get_voxel_vertices(vox)
        for arr, R in ((old_score, Rb[v]), (after_score, rep.R_pre[v])):
            info = cost_fn.evaluate(np.asarray(R, np.float64), verts, vox.phase)
            h, c = bfs._ratios(info)
            arr[v] = (float(info.cost), float(h), float(c))
    d_hit = np.abs(after_score[:, 1] - rep.hit_pre)
    d_cost = np.abs(after_score[:, 0] - rep.cost_pre)
    state = dict(
        max_abs_hit_diff=float(d_hit.max()),
        max_abs_cost_diff=float(d_cost.max()),
        n_hit_differs_over_1e9=int((d_hit > 1e-9).sum()),
        n_cost_differs_over_1e9=int((d_cost > 1e-9).sum()),
    )
    # 2. the script rule on the state rebuilt the script way
    tv = np.array([t["voxel"] for t in res.tries], np.int64)
    te = np.array([t["end"] for t in res.tries], np.float64).reshape(-1, 3, 3)
    tc = np.array([t["cost"] for t in res.tries])
    th = np.array([t["hit"] for t in res.tries])
    R, cost, hit = SR.after_pass_state(Rb, old_score, tv, te, tc, th)
    fl = SR.screen_flags(R, cost, hit, old_score[:, 1], Rb, tv, te, th, nbrs, y)
    flagged, cands = fl["flagged"], fl["cands"]
    verts_of = {}

    def fit(v: int, R0: np.ndarray, rng: np.random.Generator) -> Any:
        vox = mic.voxels[v]
        if v not in verts_of:
            verts_of[v] = _get_voxel_vertices(vox)
        r = rec.local_optimization(
            voxel_vertices=verts_of[v],
            phase_index=vox.phase,
            initial_orientation=np.asarray(R0, np.float64),
            rng=rng,
            **SR.fit_kwargs(),
        )
        h, _c = bfs._ratios(r.overlap_info)
        return (
            np.asarray(r.orientation, np.float64),
            float(r.cost),
            float(h),
            int(rec.last_local_optimization_evals),
        )

    script = SR.repair(R, cost, hit, flagged, cands, fit)
    lib_final = np.stack([np.asarray(v.orientation, np.float64) for v in mic.voxels])
    same = dict(
        extra=bool([int(v) for v in np.where(fl["extra"])[0]] == rep.extra),
        n_extra=[int(fl["extra"].sum()), len(rep.extra)],
        pass_tries_all_2_evals=bool(all(t["evals"] == 2 for t in res.tries)),
        pass_ends_are_starts=bool(all(np.array_equal(t["end"], Rb[t["cand"]]) for t in res.tries)),
        flagged=bool(flagged == rep.flagged),
        n_flagged=[len(flagged), len(rep.flagged)],
        fit_order=bool(
            [(f["voxel"], f["start"]) for f in script.fits]
            == [(f["voxel"], f["start"]) for f in rep.fits]
        ),
        n_fits=[len(script.fits), len(rep.fits)],
        fit_evals=bool([f["evals"] for f in script.fits] == [f["evals"] for f in rep.fits]),
        n_evals=[int(script.n_evals), int(rep.n_evals)],
        changed=bool(script.changed == rep.changed),
        n_changed=[len(script.changed), len(rep.changed)],
        final_map_bytes_equal=bool(script.R.tobytes() == lib_final.tobytes()),
        final_map_max_abs_diff=float(np.abs(script.R - lib_final).max()),
        fit_ends_bytes_equal=bool(
            np.array([f["end"] for f in script.fits]).tobytes()
            == np.array([f["end"] for f in rep.fits]).tobytes()
        ),
    )
    same["all_equal"] = bool(
        same["extra"]
        and same["pass_tries_all_2_evals"]
        and same["pass_ends_are_starts"]
        and same["flagged"]
        and same["fit_order"]
        and same["fit_evals"]
        and same["n_evals"][0] == same["n_evals"][1]
        and same["changed"]
        and same["final_map_bytes_equal"]
        and same["fit_ends_bytes_equal"]
    )
    # accuracy line (truth of mid, scoring only): wrong before and after the repair
    T = G.Truth()
    e_pre = np.asarray(misorientation_deg(rep.R_pre, T.R, y), float)
    e_post = np.asarray(misorientation_deg(lib_final, T.R, y), float)
    wp, wq = e_pre > G.WRONG, e_post > G.WRONG
    acc = dict(
        wrong_after_pass=int(wp.sum()),
        wrong_after_repair=int(wq.sum()),
        fixed=int((wp & ~wq).sum()),
        broken=int((~wp & wq).sum()),
    )
    return dict(state=state, script_vs_library=same, accuracy_library_run=acc)


def run(out: Path, cap_min: float) -> None:
    from icenine import boundary_pass as LBP

    label = C.preflight_label(out)
    setup, bfs, t_setup = C.build(
        C.CONFIG, out, sample=G.GRID_MIC, stack_name=G.STACK, extra="\n".join(KEYS)
    )
    mic = setup.sample.get_mic()
    p = setup.search_params
    assert p.seed_by_indexing and p.fused_cost and p.local_optimizer == "mc"
    assert (
        p.boundary_pass
        and p.boundary_pass_repair
        and p.boundary_pass_screen
        and p.boundary_pass_skip_hit == 0.0
    )
    assert p.seed_index_mc_stall_hit == 500
    captured: Dict[str, Any] = {}
    orig_pass = LBP.run_boundary_pass

    def capture(bfs_: Any, mic_: Any, **kw: Any) -> Any:
        captured["stats_before"] = {k: v for k, v in bfs_.stats.items() if k != "reseeds"}
        res = orig_pass(bfs_, mic_, **kw)
        captured["result"] = res
        captured["kw"] = kw
        return res

    LBP.run_boundary_pass = capture
    t0 = time.time()

    def watchdog() -> None:
        time.sleep(cap_min * 60.0)
        (out / "capped.json").write_text(json.dumps(dict(cap_min=cap_min)) + "\n")
        print(f"CAP of {cap_min} min reached; stopped", flush=True)
        os._exit(3)

    threading.Thread(target=watchdog, daemon=True).start()
    order = bfs.reconstruct_sample(output_mic=str(out / "recon.mic"), rng=np.random.default_rng(0))
    wall = time.time() - t0
    res = captured["result"]
    assert captured["kw"].get("repair") is True and captured["kw"].get("screen") is True, captured[
        "kw"
    ]
    chk = compare(bfs, mic, res)
    st = bfs.stats
    info = dict(
        name=NAME,
        sample=RC.SAMPLE,
        keys=KEYS,
        n_processed=len(order),
        wall_total_s=wall,
        setup_s=t_setup,
        timing_label=label,
        wall_label="contended (other jobs on the machine); not a timing result",
        library_stats={k: v for k, v in st.items() if k.startswith(("n_evals_boundary", "repair"))},
        n_boundary_voxels=len(res.boundary_voxels),
        n_pass_tries=len(res.tries),
        n_repair_flagged=len(res.repair.flagged),
        n_repair_fits=len(res.repair.fits),
        n_repair_evals=int(res.repair.n_evals),
        context_script_J5S_mid=dict(flagged=1187, evals=759_275, pass_evals=100_058),
        check=chk,
    )
    OUTJ.write_text(json.dumps(info, indent=1, default=float) + "\n")
    (out / "run.done").write_text("done\n")
    print(NAME, RC.SAMPLE, "done", wall, "s; all_equal:", chk["script_vs_library"]["all_equal"])


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--cap-min", type=float, default=240.0)
    a = ap.parse_args()
    assert RC.SAMPLE == "mid", "the library check runs on mid only"
    out = CACHE / f"{NAME}_{RC.SAMPLE}"
    SCG.require_inside(CACHE, [out])
    if out.exists():
        raise SystemExit(f"refusing to start: {out} exists")
    out.mkdir(parents=True)
    run(out, a.cap_min)
    os._exit(0)


if __name__ == "__main__":
    main()
