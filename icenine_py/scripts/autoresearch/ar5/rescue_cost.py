"""AR5 step 9, phase A step 1: the stored arm ``A`` runs reproduce step 8's totals; cost check.

  RESEED_SAMPLE=<big|grad|mid> uv run python scripts/autoresearch/ar5/rescue_cost.py

Disclosed samples only (``held`` is refused here). Reads the linked arm ``A`` run and pass of the
sample (``rescue/cache/A_<sample>``, ``passes_<sample>/D4_A``; read only) and checks the totals the
plan's stop rule names: the grid-fallback evaluations (3,726,723 ``grad``, 4,007,241 ``mid``) and
the totals with the pass (12,189,777 ``grad``, 15,645,472 ``mid``, 3,807,747 ``big``).
Finds the failed seed steps of arm ``A`` (the taken end's hit ratio under 0.80; stop rule of the
seed step: none of the kept ends reached 0.95) from ``e1_log.json.gz`` and, for ``grad`` and
``mid``, runs the cost check (prediction 1): the strict cost of the grid fallback's end (the
voxel's orientation in ``recon.mic.gz``, which a later neighbour fit may have moved by a hair; the
fallback's own end is not stored) against the cost of each kept end of the step (stored in the
log). One strict evaluation per case (23 in all). Writes ``rescue/stored_<sample>.json``;
the cases (voxel, kept ends, truth errors) are written for the replay to
``rescue/cases_<sample>.json``. Truth is used for scoring only.
"""

import json
import os
import sys
from typing import Any, Dict, List

import rescue_common as RCM  # noqa: E402  (first)

import fit_common as FC  # noqa: E402
import numpy as np  # noqa: E402
import order_logs as OL  # noqa: E402
import reseed_common as RC  # noqa: E402

import gradient_common as G  # noqa: E402

EXPECT = {
    "grad": dict(fallback_evals=3_726_723, total=12_189_777, fallbacks=11),
    "mid": dict(fallback_evals=4_007_241, total=15_645_472, fallbacks=12),
    "big": dict(fallback_evals=0, total=3_807_747, fallbacks=0),
}
N = 24570


def failed_steps(log: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
    """Steps of arm A that failed the seed test (no kept end at 0.95, taken hit under 0.80)."""
    out = []
    for i, e in enumerate(log):
        if not e["fits"]:
            continue
        t = OL.taken_end(e)
        if e["stop"] < 0 and t["hit"] < 0.80:
            out.append(dict(step=i, entry=e, kept=OL.kept_ends(e)))
    return out


def main() -> None:
    s = FC.SAMPLE
    if s == "held":
        raise SystemExit("rescue_cost.py is for the disclosed samples; held is not read here")
    rd = RCM.RCACHE / f"A_{s}"
    run = json.loads((rd / "run.json").read_text())
    pj = json.loads((RCM.RCACHE / f"passes_{s}" / "D4_A" / "pass.json").read_text())
    st = run["stats"]
    log = RC.read_gz(rd / "e1_log.json.gz")
    fs = failed_steps(log)
    exp = EXPECT[s]
    got = dict(
        fallback_evals=int(st.get("index_fallback_grid_evals", 0)),
        total=int(run["evals_total"]) + int(pj["evals"]),
        fallbacks=int(st.get("index_fallback_seed_test", 0)),
        failed_steps_in_log=len(fs),
        no_candidate_steps=sum(1 for e in log if not e["fits"]),
    )
    ok = all(got[k] == exp[k] for k in exp) and got["failed_steps_in_log"] == exp["fallbacks"]
    out: Dict[str, Any] = dict(sample=s, expected=exp, got=got, reproduces=bool(ok))
    print(json.dumps(out, indent=1))
    if not ok:
        (RCM.RESCUE / f"stored_{s}.json").write_text(json.dumps(out, indent=1) + "\n")
        print("STOP: a stored total does not reproduce", file=sys.stderr)
        sys.exit(1)
    cases: List[Dict[str, Any]] = []
    if fs:
        env = RC.make_env()
        rec = env.hooks.rec
        T = env.T
        R_map, _ = G.read_map(rd / "recon.mic.gz", N)
        records = RC.read_gz(rd / "records.json.gz")
        for f in fs:
            e = f["entry"]
            v = int(e["voxel"])
            verts = env.hooks.vertices(v)
            ph = int(env.mic.voxels[v].phase)
            strict_cf = rec._index_context(ph)[3]
            c_grid = float(strict_cf.evaluate(np.asarray(R_map[v], np.float64), verts, ph).cost)
            kept = f["kept"]
            for k in kept:
                k["err"] = float(T.err(np.asarray(k["R"], np.float64).reshape(1, 3, 3), [v])[0])
            cases.append(
                dict(
                    sample=s,
                    step=f["step"],
                    voxel=v,
                    reseed=bool(e["reseed"]),
                    kept=kept,
                    grid_end_cost=c_grid,
                    grid_end_hit=float(records[str(v)]["hit_ratio"]),
                    grid_end_err=float(T.err(R_map[v][None], [v])[0]),
                    grid_cheaper_than_all_kept=bool(c_grid < min(k["cost"] for k in kept)),
                    min_kept_cost=float(min(k["cost"] for k in kept)),
                    n_kept=len(kept),
                )
            )
    out["cases"] = len(cases)
    out["grid_end_cheaper_than_every_kept_end"] = sum(
        1 for c in cases if c["grid_cheaper_than_all_kept"]
    )
    out["strict_evaluations"] = len(cases)
    (RCM.RESCUE / f"stored_{s}.json").write_text(json.dumps(out, indent=1) + "\n")
    (RCM.RESCUE / f"cases_{s}.json").write_text(json.dumps(cases, indent=1) + "\n")
    print(json.dumps({k: v for k, v in out.items() if k != "got"}, indent=1))
    os._exit(0)


if __name__ == "__main__":
    main()
