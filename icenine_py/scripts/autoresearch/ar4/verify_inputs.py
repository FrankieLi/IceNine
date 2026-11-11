"""AR4 step 1: the pairs of one sample and the numbers of one evaluation of each.

  RESEED_SAMPLE=big|grad uv run python scripts/autoresearch/ar4/verify_inputs.py [--limit N]

(``mid`` is run only by ``verify_score_mid.py`` after the thresholds are committed.)

Pairs (plan, "Pairs"), kinds in ``KINDS``:
  after / before   every voxel with its S8 orientation after the boundary pass; before it too where
                   it differs from the one after (more than 1e-6 in the matrix)
  try_start        at every boundary voxel, each distinct neighbour orientation (the stored tries of
                   the pass: voxel and neighbour of ``starts_<sample>.npz``, the neighbour's
                   orientation from the map before the pass)
  seed_end         the refined end of every stored seed step of S8 (a replay with the library's
                   ``index_seed``: 20 costed, 5 refined, stop at hit 0.95; checked against the
                   stored log's finisher evaluations and hit ratios)
  seed_cand        the costed candidates (first 20 of the indexer's list) of every 10th seed step
  turn_01/03/05    the truth turned by 0.1, 0.3, 0.5 degree about a random axis at every 20th voxel

Per pair, one evaluation of the library's strict cost (fused, as the stored runs) with the per-peak
pixel counts of the same overlap data (``verify_common.Evaluator.peak_rows``), the number of
attached rays the orientation explains, and the error against the truth (the label's source;
never an input). The vote count among the unexplained rays costs 50 to 200 ms a pair (index_rays),
over the budget for about 10^5 pairs; it is not computed (``best_votes`` = -1; deviation reported).

Writes ``cache/inputs_<sample>.npz`` and ``cache/seeds_<sample>.json`` (untracked) and the slim
``inputs_<sample>.json`` (tracked): counts, evaluations, the sum-back check.
"""

import argparse
import json
import time
from typing import Any, Dict, List

import numpy as np
import verify_common as VC
import verify_lib as VL

KINDS = (
    "after",
    "before",
    "try_start",
    "seed_end",
    "seed_cand",
    "turn_01",
    "turn_03",
    "turn_05",
)
TURNS = {"turn_01": 0.1, "turn_03": 0.3, "turn_05": 0.5}
SEED_N_COST, SEED_N_REFINE, SEED_STOP = 20, 5, 0.95
SEED_TEST = 0.80


def replay_seed_steps(E: VC.Evaluator, limit: int) -> List[Dict[str, Any]]:
    """Replay every seed step of the stored S8 log; keep the ends, costed candidates, checks."""
    log = VC.RC.read_gz(VC.S8_DIR / "e1_log.json.gz")
    if limit:
        log = log[:limit]
    out = []
    t0 = time.time()
    for i, e in enumerate(log):
        v = int(e["voxel"])
        phase = E.mic.voxels[v].phase
        r = E.rec.index_seed(E.env.hooks.vertices(v), phase, SEED_N_COST, SEED_N_REFINE, SEED_STOP)
        hits = []
        for R in r.finisher_R:
            info = E.rec.evaluate_overlap(R, E.env.hooks.vertices(v), phase)
            hits.append(float(E.env.bfs._ratios(info)[0]))
        stored_hit = e["fin_hit"]
        passed = bool(e["stop"] >= 0 or (stored_hit and max(stored_hit) >= SEED_TEST))
        out.append(
            dict(
                step=i,
                voxel=v,
                reseed=bool(e["reseed"]),
                stored_stop=int(e["stop"]),
                stored_fin_evals=[int(x) for x in e["fin_evals"]],
                stored_fin_hit=[float(x) for x in e["fin_hit"]],
                replay_fin_evals=[int(x) for x in r.finisher_evals],
                replay_fin_hit=hits,
                evals_match=[int(x) for x in r.finisher_evals] == [int(x) for x in e["fin_evals"]],
                hit_match=bool(
                    len(hits) == len(stored_hit) and np.allclose(hits, stored_hit, atol=1e-9)
                ),
                passed_seed_test=passed,
                taken_pos=int(r.refined.index(r.taken)) if r.taken in r.refined else -1,
                stopped_at=int(r.stopped_at),
                ends=[np.asarray(x).tolist() for x in r.finisher_R],
                cand=(
                    [np.asarray(x).tolist() for x in (r.candidates[:SEED_N_COST])]
                    if r.candidates is not None and i % 10 == 0
                    else []
                ),
                stored_evals=int(e["evals"]),
            )
        )
        if i % 100 == 0:
            print("seed step", i, len(log), f"{time.time() - t0:.0f} s", flush=True)
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--limit", type=int, default=0)
    a = ap.parse_args()
    s = VC.SAMPLE
    assert s in VC.TRAIN or s == VC.TEST
    t_start = time.time()
    E = VC.Evaluator()
    before = np.load(VC.PASS_DIR / "before_S8.npy")
    after = np.load(VC.PASS_DIR / "after_S8.npy")
    tries = np.load(VC.ICE / "benchmarks" / "autoresearch" / "ar5" / f"starts_{s}.npz")
    VC.CACHE.mkdir(parents=True, exist_ok=True)
    seeds_path = VC.CACHE / f"seeds_{s}.json"
    suffix = f"_limit{a.limit}" if a.limit else ""
    if seeds_path.exists() and not a.limit:
        seeds = json.loads(seeds_path.read_text())["steps"]
        print("seed steps read from", seeds_path.name, len(seeds), flush=True)
    else:
        seeds = replay_seed_steps(E, a.limit)
        (VC.CACHE / f"seeds_{s}{suffix}.json").write_text(json.dumps(dict(steps=seeds)) + "\n")
    t_seeds = time.time()

    # ---- the task list: (voxel, orientation, kind, aux, step)
    tasks: List[Any] = []
    n = VC.N
    vox = np.arange(n)
    if a.limit:
        vox = vox[:: n // 600]
    for v in vox:
        tasks.append((int(v), after[v], 0, 0, -1))
        if not np.allclose(before[v], after[v], atol=1e-6):
            tasks.append((int(v), before[v], 1, 0, -1))
    tv, tc = tries["voxel"], tries["cand"]
    if a.limit:
        keep = np.isin(tv, vox)
        tv, tc = tv[keep], tc[keep]
    for v, c in zip(tv, tc):
        tasks.append((int(v), before[int(c)], 2, int(c), -1))
    for st in seeds:
        for pos, R in enumerate(st["ends"]):
            tasks.append((st["voxel"], np.asarray(R), 3, pos, st["step"]))
        for j, R in enumerate(st["cand"]):
            tasks.append((st["voxel"], np.asarray(R), 4, j, st["step"]))
    truth_R = E.env.Rtrue
    axes = VL.random_axes(3 * len(range(0, n, 20)), seed=0)
    for q, v in enumerate(range(0, n, 20)):
        if a.limit and v not in set(vox.tolist()):
            continue
        for m, deg in enumerate((0.1, 0.3, 0.5)):
            R = VL.rotate_about_axis(truth_R[v], deg, axes[3 * q + m])
            tasks.append((int(v), R, 5 + m, 0, -1))
    tasks.sort(key=lambda t: (t[0], t[2]))
    print(SAMPLE_LINE(s, len(tasks)), flush=True)

    # ---- one evaluation per pair
    P = len(tasks)
    voxel = np.array([t[0] for t in tasks], np.int32)
    kind = np.array([t[2] for t in tasks], np.int8)
    aux = np.array([t[3] for t in tasks], np.int32)
    step = np.array([t[4] for t in tasks], np.int32)
    R_all = np.stack([np.asarray(t[1], np.float64).reshape(9) for t in tasks])
    hit = np.zeros(P)
    cost = np.zeros(P)
    conf = np.zeros(P)
    n_expl = np.zeros(P)
    d_po = np.zeros(P)
    d_pod = np.zeros(P)
    d_stage = np.zeros(P)
    d_hit = np.zeros(P)
    counts = np.zeros((P, 3), np.int64)
    n_peaks = np.zeros(P, np.int32)
    rows: List[np.ndarray] = []
    t0 = time.time()
    for i, t in enumerate(tasks):
        v, R = t[0], np.asarray(t[1], np.float64)
        r = E.peak_rows(R, v)
        hit[i], cost[i], conf[i] = r["hit"], r["cost"], r["conf"]
        counts[i] = r["counts"]
        rw = r["rows"]
        rows.append(rw)
        n_peaks[i] = len(rw)
        if "counts_rows" in r:
            po, pod, _ = r["counts_rows"]
            d_po[i] = po - r["counts"][0]
            d_pod[i] = pod - r["counts"][1]
            d_stage[i] = abs(po - r["counts_stage_d"][0]) + abs(pod - r["counts_stage_d"][1])
            hs = float((rw[:, 3] * rw[:, 4]).sum() / max(rw[:, 4].sum(), 1))
            d_hit[i] = hs - r["hit"]
        n_expl[i] = E.explained_and_votes_cheap(R, v)
        if i % 5000 == 0:
            print(i, P, f"{time.time() - t0:.0f} s", flush=True)
    t_end = time.time()
    flat = np.concatenate(rows) if rows else np.zeros((0, len(VC.PEAK_COLS)), np.float32)
    offsets = np.concatenate([[0], np.cumsum([len(x) for x in rows])]).astype(np.int64)
    # ---- the error against the truth (labels only), chunked
    err = np.zeros(P)
    for lo in range(0, P, 4000):
        sl = slice(lo, lo + 4000)
        err[sl] = np.asarray(
            E.env.T.err(R_all[sl].reshape(-1, 3, 3), voxel[sl].astype(np.int64)), dtype=float
        )
    shares = [x[:, 3] for x in rows]
    summ = np.stack(
        [VL.summary_vector(hit[i], cost[i], conf[i], shares[i], n_expl[i], -1.0) for i in range(P)]
    )
    np.savez_compressed(
        VC.CACHE / f"inputs_{s}{suffix}.npz",
        voxel=voxel,
        kind=kind,
        aux=aux,
        step=step,
        R=R_all,
        err=err,
        summary=summ,
        n_peaks=n_peaks,
        counts=counts,
        peaks=flat,
        offsets=offsets,
    )
    seed_ok = [x["evals_match"] and x["hit_match"] for x in seeds]
    info = dict(
        sample=s,
        n_pairs=int(P),
        n_evaluations=int(P),
        library_calls_per_pair=2,
        by_kind={k: int((kind == j).sum()) for j, k in enumerate(KINDS)},
        n_right=int(VL.right_label(err).sum()),
        n_seed_steps=len(seeds),
        seed_replay_evals_and_hits_match=int(sum(seed_ok)),
        seed_replay_mismatch=int(len(seed_ok) - sum(seed_ok)),
        sum_back=dict(
            n_pairs_checked=int((n_peaks > 0).sum()),
            max_abs_pixel_overlap_diff_vs_fused=float(np.abs(d_po).max()),
            max_abs_pixel_on_detector_diff_vs_fused=float(np.abs(d_pod).max()),
            max_abs_hit_ratio_diff_vs_fused=float(np.abs(d_hit).max()),
            n_pairs_pixel_counts_differ_from_fused=int(((d_po != 0) | (d_pod != 0)).sum()),
            max_abs_pixel_diff_vs_stage_d_aggregate=float(d_stage.max()),
        ),
        wall_seed_replay_s=t_seeds - t_start,
        wall_pairs_s=t_end - t0,
        timing_label="contended (another agent on the machine); not a timing result",
    )
    (VC.AR4 / f"inputs_{s}{suffix}.json").write_text(json.dumps(info, indent=1) + "\n")
    print(json.dumps(info), flush=True)


def SAMPLE_LINE(s: str, n: int) -> str:
    return f"{s}: {n} pairs"


if __name__ == "__main__":
    main()
