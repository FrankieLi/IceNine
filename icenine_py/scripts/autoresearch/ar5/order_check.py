"""AR5 step 8, Phase A step 4: do the running orders give what the stored replay gives?

  RESEED_SAMPLE=mid uv run python scripts/autoresearch/ar5/order_check.py [--n 20]

Takes 20 stored seed steps of ``mid`` that have candidates (10 where the two orders spend a
different number of evaluations, 10 where they spend the same; evenly spaced in step order), runs
each step with the running code of the whole-layer runs (``order_run.OrderPick.e1_step``) in both
orders, and compares with the replay: every fit that ran (CMA end or MC end of a candidate: hit
ratio, cost, evaluations) against the stored per-candidate ends, and the step (evaluations, taken
candidate, hit ratio of the taken end, MC refinements, stopped early) against
``mix_rules.mixed_step``
(interleaved) and ``order_rules.sweep_step`` (sweep). An MC end the replay did not compute (its
arm had stopped earlier) is counted as ``not_stored``, not compared. Truth is used for the error
of the taken end only. Writes ``benchmarks/autoresearch/ar5/order/check.json``.
"""

import argparse
import json
from typing import Any, Dict, List

import order_common as OC  # noqa: E402  (first)

import fit_part1_lib as L  # noqa: E402
import mix_rules as MR  # noqa: E402
import numpy as np  # noqa: E402
import order_rules as OR  # noqa: E402
import order_run as ORUN  # noqa: E402
import reseed_common as RC  # noqa: E402


def spread(xs: List[Any], k: int) -> List[Any]:
    idx = np.linspace(0, len(xs) - 1, k).round().astype(int)
    return [xs[i] for i in idx]


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=20)
    a = ap.parse_args()
    assert OC.SAMPLE == "mid"
    from icenine import ray_indexing as RI

    rows = [r for r in L.load("mid") if r["arms"]]
    diff = [
        r
        for r in rows
        if MR.mixed_step(r["costed"], r["arms"]["C200"], r["arms"]["M"])["evals"]
        != OR.sweep_step(r["costed"], r["arms"]["C200"], r["arms"]["M"])["evals"]
    ]
    same = [r for r in rows if r not in diff]
    chosen = spread(diff, a.n // 2) + spread(same, a.n - a.n // 2)
    env = RC.make_env()
    rec = env.hooks.rec
    sym = RI.proper_rotations(env.setup.exp_setup.get_sample_symmetry())
    pick = ORUN.OrderPick(env.bfs, env.mic, rec.local_optimization, sym, use_tries=False)
    out_rows: List[Dict[str, Any]] = []
    for row in chosen:
        v = int(row["voxel"])
        verts = env.hooks.vertices(v)
        ph = int(env.mic.voxels[v].phase)
        res_row: Dict[str, Any] = dict(step=row["step"], voxel=v)
        for order, step_fn in (("A", MR.mixed_step), ("Y", OR.sweep_step)):
            ORUN.ORD["order"] = order
            pick.cur, pick.reseed = v, False
            n_log = len(pick.e1_log)
            out = pick.e1_step(verts, ph, n_cost=20, n_refine=5)
            e = pick.e1_log[n_log]
            exp = step_fn(row["costed"], row["arms"]["C200"], row["arms"]["M"])
            bad: List[str] = []
            not_stored = 0
            for pos, j, kind, hit, cost, ev, _R in e["fits"]:
                arm = row["arms"]["C200" if kind == "C" else "M"]
                if pos >= len(arm["fin_hit"]):
                    not_stored += 1
                    continue
                if (
                    abs(hit - arm["fin_hit"][pos]) > 1e-12
                    or abs(cost - arm["fin_cost"][pos]) > 1e-12
                    or ev != arm["fin_evals"][pos]
                ):
                    bad.append(f"fit {kind}{pos}")
            live = dict(
                evals=int(e["evals"]),
                taken=int(out.refined.index(out.taken)) if out.taken is not None else None,
                stopped=bool(e["stop"] >= 0),
                mc_ran=int(sum(1 for f in e["fits"] if f[2] == "M")),
            )
            live["hit"] = float(e["fin_hit"][live["taken"]]) if live["taken"] is not None else 0.0
            if not exp["missing"]:
                for k in ("evals", "taken", "stopped", "mc_ran"):
                    if live[k] != exp[k]:
                        bad.append(f"step {k}: live {live[k]} stored {exp[k]}")
                if abs(live["hit"] - exp["hit"]) > 1e-12:
                    bad.append("step hit")
            err = float(env.T.err(np.asarray(out.orientation)[None], [v])[0])
            if not exp["missing"] and abs(err - exp["err"]) > 1e-9:
                bad.append(f"error {err} vs {exp['err']}")
            res_row[order] = dict(
                live=live,
                stored={k: exp[k] for k in ("evals", "taken", "stopped", "mc_ran", "hit")},
                n_fits=len(e["fits"]),
                not_stored=not_stored,
                missing_in_replay=bool(exp["missing"]),
                bad=bad,
                same=not bad,
            )
        out_rows.append(res_row)
        print(row["step"], v, {o: res_row[o]["same"] for o in ("A", "Y")}, flush=True)
    summ = dict(
        sample="mid",
        n_steps=len(out_rows),
        n_orders_same={o: sum(r[o]["same"] for r in out_rows) for o in ("A", "Y")},
        n_fits_compared=sum(r[o]["n_fits"] - r[o]["not_stored"] for r in out_rows for o in "AY"),
        n_fits_not_stored=sum(r[o]["not_stored"] for r in out_rows for o in "AY"),
        n_steps_missing=sum(r[o]["missing_in_replay"] for r in out_rows for o in "AY"),
        all_same=all(r[o]["same"] for r in out_rows for o in "AY"),
        evals_live={o: sum(r[o]["live"]["evals"] for r in out_rows) for o in ("A", "Y")},
        evals_stored={o: sum(r[o]["stored"]["evals"] for r in out_rows) for o in ("A", "Y")},
    )
    (OC.ORDER / "check.json").write_text(
        json.dumps(dict(summary=summ, rows=out_rows), indent=1) + "\n"
    )
    print(json.dumps(summ, indent=1))


if __name__ == "__main__":
    main()
