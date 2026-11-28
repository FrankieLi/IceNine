"""AR5 step 2, Part 0: the stored numbers of run S8 (seed refinements and MC fallbacks).

  scripts/autoresearch/ar5/fit_bg.sh <sample> part0_<sample> fit_part0.py

Reads, for one sample, the stored S8 run (``run.json``, ``e1_log.json.gz``) and its boundary pass
(``passes_<sample>/D4_S8/pass.json``) and writes
``benchmarks/autoresearch/ar5/fit/part0_<sample>.json``: the seed steps, the refinements and
their evaluations, the costing evaluations, the MC fallbacks of
the neighbour fits (runs, evaluations, runs that ended lower than the CMA end, ``fallback_won``),
the total with the pass, and the shares of the total. Nothing is run.
"""

import json

import fit_common as FC
import reseed_common as RC


def main() -> None:
    s = FC.SAMPLE
    d = FC.STORED_S8[s]
    run = json.loads((d / "run.json").read_text())
    log = RC.read_gz(d / "e1_log.json.gz")
    pj = json.loads((FC.STORED_PASS[s] / "D4_S8" / "pass.json").read_text())
    st, w = run["stats"], run["wrapper"]
    n_ref = sum(len(e["fin_evals"]) for e in log)
    ev_ref = sum(sum(e["fin_evals"]) for e in log)
    ev_cost = sum(e["costed"] for e in log)
    ev_seed = sum(e["evals"] for e in log if "evals" in e)
    grid = int(st.get("index_fallback_grid_evals", 0))
    pre = int(run["evals_total"])
    total = pre + int(pj["evals"])
    ev_fb = int(w["evals_fallback"])
    out = dict(
        sample=s,
        seed_steps=len(log),
        seed_steps_stopped_early=sum(1 for e in log if e["stop"] >= 0),
        seed_refinements=n_ref,
        seed_refinement_evals=ev_ref,
        seed_refinement_evals_mean=ev_ref / max(1, n_ref),
        seed_costing_evals=ev_cost,
        seed_step_evals=ev_seed,
        seed_step_evals_check=int(run["pick"]["evals"]),
        grid_fallback_evals=grid,
        grid_fallbacks=int(st.get("index_fallback_seed_test", 0)),
        neighbour_fit_calls=int(w["wrapper_calls"]),
        mc_fallback_runs=int(w["fallback_runs"]),
        mc_fallback_share=w["fallback_runs"] / w["wrapper_calls"],
        mc_fallback_evals=ev_fb,
        mc_fallback_evals_mean=ev_fb / max(1, w["fallback_runs"]),
        mc_fallback_won=int(w["fallback_won"]),
        evals_before_pass=pre,
        pass_evals=int(pj["evals"]),
        evals_with_pass=total,
        share_seed_refinements=ev_ref / total,
        share_mc_fallback=ev_fb / total,
        share_both=(ev_ref + ev_fb) / total,
        stages=dict(
            seed=int(st["n_evals_seed"]),
            reseed=int(st.get("n_evals_reseed", 0)),
            neighbour=int(st["n_evals_neighbor"]),
            revisit=int(st["n_evals_revisit"]),
        ),
    )
    (FC.FIT / f"part0_{s}.json").write_text(json.dumps(out, indent=1) + "\n")
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
