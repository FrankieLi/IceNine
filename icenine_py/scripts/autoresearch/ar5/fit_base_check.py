"""AR5 step 2, Part 2: the base check. Does the script layer with both changes off give stored S8?

  scripts/autoresearch/ar5/fit_bg.sh big base_check fit_base_check.py

Compares ``<cache>/S8x_big/`` (``fit_run.py --name S8x``, arm M, MC fallback on) with the stored S8
run of ``big``: evaluations by stage and in total (7,201,863 stored), the seed-step log (voxel,
reseed, candidates, costed, evaluations, refinement evaluations, hit ratios, stop), the wrapper
counts, the final map (orientations as written to the .mic file) and the wrong count (voxels more
than 1 degree from the truth; 67 stored). Writes
``benchmarks/autoresearch/ar5/fit/base_check.json``.
"""

import json

import fit_common as FC
import numpy as np
import reseed_common as RC

import gradient_common as G  # noqa: E402

N = 24570


def main() -> None:
    assert FC.SAMPLE == "big"
    new_d = FC.CACHE / "S8x_big"
    old_d = FC.STORED_S8["big"]
    new = json.loads((new_d / "run.json").read_text())
    old = json.loads((old_d / "run.json").read_text())
    keys = ["n_evals_seed", "n_evals_reseed", "n_evals_neighbor", "n_evals_revisit", "index_evals"]
    stage = {
        k: [new["stats"][k], old["stats"][k], new["stats"][k] == old["stats"][k]] for k in keys
    }
    wk = ["wrapper_calls", "picked_fell", "fallback_runs", "fallback_won", "evals_fallback"]
    wrap = {
        k: [new["wrapper"][k], old["wrapper"][k], new["wrapper"][k] == old["wrapper"][k]]
        for k in wk
    }
    ln, lo = RC.read_gz(new_d / "e1_log.json.gz"), RC.read_gz(old_d / "e1_log.json.gz")
    log_same = ln == lo
    Rn, _ = G.read_map(new_d / "recon.mic.gz", N)
    Ro, _ = G.read_map(old_d / "recon.mic.gz", N)
    T = G.Truth()
    en, eo = T.err(Rn), T.err(Ro)
    out = dict(
        evals_total=[
            new["evals_total"],
            old["evals_total"],
            new["evals_total"] == old["evals_total"],
        ],
        stages=stage,
        wrapper=wrap,
        seed_log_identical=bool(log_same),
        seed_steps=[len(ln), len(lo)],
        map_max_abs_diff=float(np.abs(Rn - Ro).max()),
        map_identical=bool(np.array_equal(Rn, Ro)),
        wrong_new=int((en > G.WRONG).sum()),
        wrong_stored=int((eo > G.WRONG).sum()),
        error_arrays_identical=bool(np.array_equal(en, eo)),
    )
    out["reproduces"] = bool(
        out["evals_total"][2]
        and all(v[2] for v in stage.values())
        and out["seed_log_identical"]
        and out["map_identical"]
        and out["wrong_new"] == out["wrong_stored"] == 67
    )
    (FC.FIT / "base_check.json").write_text(json.dumps(out, indent=1) + "\n")
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
