"""AR5 step 7, Phase B: the library boundary pass on a stored map, with the key on or off.

  uv run python scripts/autoresearch/ar5/gate_live.py --sample mid --run S11n --h 0      # parity
  uv run python scripts/autoresearch/ar5/gate_live.py --sample mid --run S11n --h 0.95   # live

Reads the map before the pass of the stored run (``before_<run>.npy``, the run's ``.mic`` with its
rounded angles) into the BFS's map in memory, every voxel FITTED, and applies
``icenine.boundary_pass.run_boundary_pass`` with ``skip_hit=h`` (as step 5's ``keys_pass.py`` does
for the key-off pass). With ``--h 0`` the check is exact: every try (voxel, candidate, end, cost,
hit ratio, confidence, evaluations) and the map after the pass equal the stored ones (stop rule
otherwise). With ``--h`` > 0 it compares the live pass with the replay of ``gate_replay.py``: the
skipped voxels, the evaluations (stop rule if they differ), the wrong voxels, fixed / broken
against the full pass and against S8, the median error, the run total. Truth is used for scoring
only. The pass has one generator, so the tries after a skipped voxel need not equal the stored
ones; the evaluation count (61 per try) does not depend on it. Results go to
``benchmarks/autoresearch/ar5/gate/live_<sample>_<run>_h<h>.json`` (tracked, slim) and the maps to
``gate/cache/`` (not committed). Wall times are contended and are not timing results.
"""

import argparse
import json
import os
import sys
import time
from typing import Any, Dict

import gate_common as GC

_ap = argparse.ArgumentParser()
_ap.add_argument("--sample", required=True, choices=["big", "grad", "mid", "held"])
_ap.add_argument("--run", required=True)
_ap.add_argument("--h", type=float, required=True)
ARGS = _ap.parse_args()
if ARGS.sample != "grad":  # ar3_common reads the sample kind when it is imported
    os.environ["AR3_SAMPLE"] = ARGS.sample
sys.path.insert(0, str(GC.ICE / "scripts" / "autoresearch" / "ar3"))

import ar3_common as C  # noqa: E402
import gate_replay as RP  # noqa: E402
import gate_rules as GR  # noqa: E402
import gradient_common as G  # noqa: E402
import numpy as np  # noqa: E402


def main() -> None:
    from icenine import boundary_pass as LBP
    from icenine.mic_file import ReconstructionState

    sample, run, h = ARGS.sample, ARGS.run, ARGS.h
    tag = f"live_{sample}_{run}_h{h:g}"
    out = GC.GATE / "cache" / tag
    if out.exists():
        raise SystemExit(f"refusing to start: {out} exists")
    out.mkdir(parents=True)
    label = C.preflight_label(out)
    P = GC.load_pass(sample, run)
    stored = np.load(GC.pass_files(sample, run)["tries"])
    setup, bfs, _t = C.build(C.CONFIG, out, sample=G.GRID_MIC, stack_name=G.STACK)
    mic = setup.sample.get_mic()
    R = P["before"]
    assert len(R) == len(mic.voxels) == GC.N
    for v, r in zip(mic.voxels, R):
        v.orientation = r.copy()
        v.reconstruction_id = ReconstructionState.FITTED
    bfs.stats = bfs._new_stats()
    t0 = time.time()
    res = LBP.run_boundary_pass(bfs, mic, skip_hit=h) if h > 0 else LBP.run_boundary_pass(bfs, mic)
    wall = time.time() - t0
    live = np.stack([np.asarray(v.orientation, np.float64) for v in mic.voxels])
    np.save(out / "after_live.npy", live)

    old_hit, old_cost = P["old_score"][:, 1], P["old_score"][:, 0]
    rep = GR.replay(
        old_hit, old_cost, P["voxel"], P["hit"], P["cost"], P["evals"], h if h > 0 else 2.0
    )
    n = len(res.tries)
    lv = np.array([t["voxel"] for t in res.tries])
    info: Dict[str, Any] = dict(
        sample=sample,
        run=run,
        h=h,
        wall_s=wall,
        timing_label=label,
        live_tries=n,
        replay_tries=rep["tries_kept"],
        live_pass_evals=int(res.n_evals_tries),
        replay_pass_evals=int(rep["evals"]),
        stored_full_pass_evals=int(rep["evals_full"]),
        live_score_evals=int(res.n_evals_score),
        stats={
            k: v for k, v in bfs.stats.items() if k.startswith(("boundary", "n_evals_boundary"))
        },
        boundary_voxels=len(res.boundary_voxels),
        live_skipped=len(res.skipped),
        replay_skipped=int(rep["voxels_skipped"]),
        same_skipped_voxels=bool(set(res.skipped) == set(np.flatnonzero(rep["skipped"]).tolist())),
        evals_equal_replay=bool(res.n_evals_tries == rep["evals"] and n == rep["tries_kept"]),
        own_hit_max_abs_diff=float(
            max(abs(res.old_score[v][1] - old_hit[v]) for v in res.old_score)
        ),
    )
    # tries against the stored tries of the same voxels (exact when nothing is skipped)
    kept = rep["kept"]
    if n == int(kept.sum()):
        for f in ("voxel", "cand", "evals", "cost", "hit", "conf"):
            lib = np.array([t[f] for t in res.tries])
            info[f"differ_{f}"] = int((lib != stored[f][kept]).sum())
        end = np.array([t["end"] for t in res.tries], np.float32).reshape(-1, 3, 3)
        info["differ_end"] = int((~(end == stored["end"][kept]).reshape(n, -1).all(1)).sum())
    info["voxels_changed"] = int(len(res.changed))
    changed_stored = np.flatnonzero((P["after"] != P["before"]).reshape(GC.N, -1).any(1))
    info["stored_voxels_changed"] = int(len(changed_stored))
    if h == 0:
        info["map_equals_stored_after"] = bool(np.array_equal(live, P["after"]))
        info["same_changed_voxels"] = bool(set(res.changed) == set(changed_stored.tolist()))
        info["tries_equal_stored"] = bool(
            n == len(stored["voxel"])
            and all(info.get(f"differ_{f}") == 0 for f in ("voxel", "cand", "evals", "cost"))
            and all(info.get(f"differ_{f}") == 0 for f in ("hit", "conf", "end"))
        )
    # scores
    T = GC.Truth(sample)
    s8 = GC.s8_baseline(sample) if sample != "held" else None
    if sample == "held":
        s8d = GC.HELD_S8 / "passes_held"
        s8 = GC.s8_from_maps(
            T,
            GC.HELD_S8 / "S8_held" / "run.json",
            s8d / "D4_S8" / "pass.json",
            s8d / "after_S8.npy",
        )
    e_live = T.err(live)
    e_before, e_full = T.err(P["before"]), T.err(P["after"])
    w_before, w_full, w_live = (e > GC.WRONG for e in (e_before, e_full, e_live))
    w_s8 = s8["err_after"] > GC.WRONG
    fixed, broken = GR.fixed_broken(w_full, w_live)
    fx8, br8 = GR.fixed_broken(w_s8, w_live)
    # the replay map's errors at the same h, for the "within 2" comparison
    e_rep = T.err(RP.replay_map(P, rep))
    w_rep = e_rep > GC.WRONG
    fr, br = GR.fixed_broken(w_rep, w_live)
    side = RP.side(e_live)
    total = P["bfs_evals"] + int(res.n_evals_tries)
    info.update(
        wrong_before=int(w_before.sum()),
        wrong_full=int(w_full.sum()),
        wrong_replay=int(w_rep.sum()),
        live=side,
        median_full=float(np.median(e_full)),
        median_replay=float(np.median(e_rep)),
        fixed_vs_full=fixed,
        broken_vs_full=broken,
        fixed_vs_s8=fx8,
        broken_vs_s8=br8,
        fixed_vs_replay=fr,
        broken_vs_replay=br,
        total_with_pass=total,
        total_with_full_pass=P["bfs_evals"] + int(rep["evals_full"]),
        s8_total=int(s8["stages"]["total"]),
        s8_median=float(np.median(s8["err_after"])),
        bfs_evals=P["bfs_evals"],
        cut=1.0 - res.n_evals_tries / rep["evals_full"],
    )
    info["rule"] = GR.gate_rule(
        cut=info["cut"],
        wrong_before=info["wrong_before"],
        wrong_full=info["wrong_full"],
        fixed=fixed,
        broken=broken,
        med_full=info["median_full"],
        med_skip=side["median_err_deg"],
        total_skip=total,
        total_s8=info["s8_total"],
        fixed_s8=fx8,
        broken_s8=br8,
        med_s8=info["s8_median"],
    )
    (GC.GATE / f"{tag}.json").write_text(json.dumps(info, indent=1, default=float) + "\n")
    print(json.dumps({k: v for k, v in info.items() if k != "stats"}, indent=1, default=float))


if __name__ == "__main__":
    main()
