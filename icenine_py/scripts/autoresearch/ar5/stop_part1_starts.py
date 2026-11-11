"""AR5 step 1, part 1a: the start score of every stored try of S8's boundary pass.

  RESEED_SAMPLE=big|grad|mid uv run python scripts/autoresearch/ar5/stop_part1_starts.py

A try starts from a neighbour's orientation (the map before the pass). Its *start score* is the
strict cost and hit ratio of that orientation at the voxel: one evaluation per try
(``boundary_common.Env.score3``, the score the pass itself uses for the voxels' old ends). The
script also checks that the stored tries, applied with rule ``hit`` and no test, give the stored
map after the pass, and scores each try's end against the truth (scoring only: no rule reads it).

Writes ``benchmarks/autoresearch/ar5/starts_<sample>.npz`` (slim, tracked): per try voxel, cand,
start_hit, start_cost, end_hit, end_cost, evals, end_err; per voxel old_hit, old_cost, err_before,
err_after. Timings are contended (another agent runs on the machine).
"""

import json
import time

import numpy as np
import stop_common as SC
import stop_rules as SR

import boundary_common as B  # noqa: E402
import boundary_lib as BL  # noqa: E402
import reseed_common as RC  # noqa: E402


def main() -> None:
    s = SC.SAMPLE
    t0 = time.time()
    t = BL.Tries(SC.TRIES_DIR)
    before = np.load(SC.PASS_DIR / "before_S8.npy")
    after = np.load(SC.PASS_DIR / "after_S8.npy")
    env: B.Env = RC.make_env()
    assert env.n == SC.N == len(before)
    # the pass as stored: rule hit, no test, equals the stored map after the pass
    keep_all = np.ones(len(t.voxel), bool)
    ch = SR.choose_try(t.voxel, t.old_score[:, 1], t.old_score[:, 0], t.hit, t.cost, keep_all)
    again = before.copy()
    took = ch >= 0
    again[took] = t.end[ch[took]]
    same = bool(np.allclose(again, after, atol=1e-6))
    assert same, "stored tries with rule hit do not give the stored map after the pass"
    # start scores: one evaluation per try
    sh = np.zeros(len(t.voxel))
    sc = np.zeros(len(t.voxel))
    for i, (v, c) in enumerate(zip(t.voxel, t.cand)):
        cost, hit, _ = env.score3(before[c], int(v))
        sc[i], sh[i] = cost, hit
    t1 = time.time()
    # truth, for scoring only
    end_err = np.asarray(env.T.err(t.end, t.voxel), dtype=float)
    e = np.load(SC.errors_path(s))
    eb, ea = e["S8_before"], e["S8_after"]
    assert np.allclose(env.err(before), eb, atol=1e-4) and np.allclose(
        env.err(after), ea, atol=1e-4
    )
    np.savez_compressed(
        SC.AR5 / f"starts_{s}.npz",
        voxel=t.voxel.astype(np.int32),
        cand=t.cand.astype(np.int32),
        start_hit=sh.astype(np.float32),
        start_cost=sc.astype(np.float32),
        end_hit=t.hit.astype(np.float32),
        end_cost=t.cost.astype(np.float32),
        evals=t.evals.astype(np.int32),
        end_err=end_err.astype(np.float32),
        old_hit=t.old_score[:, 1].astype(np.float32),
        old_cost=t.old_score[:, 0].astype(np.float32),
        err_before=eb.astype(np.float32),
        err_after=ea.astype(np.float32),
    )
    info = dict(
        sample=s,
        n_tries=int(len(t.voxel)),
        pass_evals=int(t.evals.sum()),
        start_evals=int(len(t.voxel)),
        wrong_before=int((eb > SC.WRONG).sum()),
        wrong_after=int((ea > SC.WRONG).sum()),
        stored_map_reproduced=same,
        wall_total_s=time.time() - t0,
        wall_start_scores_s=t1 - t0,
        timing_label="contended (another agent on the machine); not a timing result",
    )
    (SC.AR5 / f"starts_{s}.json").write_text(json.dumps(info, indent=1) + "\n")
    print(json.dumps(info), flush=True)


if __name__ == "__main__":
    main()
