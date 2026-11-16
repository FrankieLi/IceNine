"""AR5 step 22, phase A step 5 (gates g2 and g3): the live ``J5D`` run on ``mid`` against the
replay, and ``J5SD`` against ``J5S`` on ``mid`` and ``held``.

  RESEED_SAMPLE=mid uv run python scripts/autoresearch/ar5/defer_gate.py

Reads the runs of ``defer_run_all.sh --stage gate`` (``<cache>/J5D_<s>`` and the screen arm's
``<cache>/j5sd/J5S_<s>``) and step 21's stored ``J5S`` maps (read only). The truth of the
disclosed samples is used for scoring only. Nothing of ``held6`` is read.

* **g2** (``mid``): the seed steps of the live ``J5D`` run equal the replay in every seed step up to
  the first type C step in run order: the same voxel, the same stop position, the taken end equal
  to ``J5``'s and the evaluations equal to ``J5``'s minus the paused tails (``walk_step``).
* **g3** (``mid`` and ``held``): ``J5SD`` against ``J5S``: at least 7% fewer evaluations
  (``J5D`` search + the screen pass + the wide repair, against step 21's ``J5S`` total), broken
  minus fixed at most 3 (wrong: over 1 degree from the truth, symmetry-reduced; counted by voxel),
  median error over all voxels at most 0.001 degree higher, grid searches at most 2 more than
  ``J5``'s.

Writes ``ar5/defer/gate.json`` and ``gate_tables.md`` (blocks for ``sync_doc_tables.py``). Exit 1
when a gate fails (the plan's stop rule: no tuning, no threshold, no ``held6``). The file
``gate_passed`` that opens ``held6`` is NOT written here.
"""

import json
import sys
from pathlib import Path
from typing import Any, Dict, List

import defer_common as DC  # noqa: E402  (first)

import defer_g0 as G0  # noqa: E402
import defer_rules as DR  # noqa: E402
import numpy as np  # noqa: E402
import screen_data as SD  # noqa: E402
import screen_rules as SR  # noqa: E402

sys.path.insert(0, str(DC.ICE / "scripts" / "common"))
import doc_tables as DT  # noqa: E402
import stats as ST  # noqa: E402


def g2(cache: Path) -> Dict[str, Any]:
    steps = G0.load_steps(cache, "mid")
    live = DC.load_gz(cache / "J5D_mid" / "e1_log.json.gz")
    first_c = next((st["index"] for st in steps if st["walk"]["type"] == "C"), len(steps))
    bad: List[str] = []
    n = 0
    for st in steps[:first_c]:
        i = st["index"]
        if i >= len(live):
            bad.append(f"step {i}: the live run has fewer steps")
            break
        e, ne, w = st["entry"], live[i], st["walk"]
        n += 1
        if int(ne["voxel"]) != int(e["voxel"]) or bool(ne["reseed"]) != bool(e["reseed"]):
            bad.append(f"step {i}: voxel {ne['voxel']} vs {e['voxel']}")
        elif int(ne["stop"]) != int(e["stop"]):
            bad.append(f"step {i}: stop {ne['stop']} vs {e['stop']}")
        elif int(ne["evals"]) != int(w["expected_evals"]):
            bad.append(f"step {i}: evals {ne['evals']} vs expected {w['expected_evals']}")
        elif not G0.same_R(G0.taken_row(ne)[6], G0.taken_row(e)[6]):
            bad.append(f"step {i}: taken end differs")
    return dict(
        ok=bool(not bad and n == first_c and first_c > 0),
        steps_compared=n,
        first_c_step=first_c,
        live_steps=len(live),
        j5_steps=len(steps),
        differences=bad[:20],
        n_differences=len(bad),
    )


def g3(cache: Path, s: str, total_j5s: int) -> Dict[str, Any]:
    run = cache / f"J5D_{s}"
    arm = cache / "j5sd" / f"J5S_{s}"
    st = json.loads(DC.DG.guarded(run / "run.json").read_text())["stats"]
    search = int(
        st["n_evals_seed"] + st["n_evals_reseed"] + st["n_evals_neighbor"] + st["n_evals_revisit"]
    )
    rep = json.loads(DC.DG.guarded(arm / "repair.json").read_text())
    total = search + int(rep["pass_evals"]) + int(rep["n_evals"])
    y = SD.sym()
    truth = SD.read_truth(s)
    R_t = truth["R_t"]
    new = np.load(DC.DG.guarded(arm / "after.npy"))
    old = np.load(DC.DG.guarded(cache / f"stored_J5S_{s}" / "after.npy"))
    w_old, w_new = SR.wrong_mask(old, R_t, y), SR.wrong_mask(new, R_t, y)
    fixed, broken = SR.fixed_broken(w_old, w_new)
    e_old = np.asarray(SR.misorientation_deg(old, R_t, y), float)
    e_new = np.asarray(SR.misorientation_deg(new, R_t, y), float)
    pol_d = json.loads(DC.DG.guarded(run / "polish.json").read_text())
    pol_j = json.loads(DC.DG.guarded(cache / f"J5_{s}" / "polish.json").read_text())
    grid_d, grid_j = int(pol_d["grid_fallbacks_left"]), int(pol_j["grid_fallbacks_left"])
    g = DR.GATE
    fewer = 1.0 - total / total_j5s
    n = len(w_new)
    lo, hi = ST.wilson(int(w_new.sum()), n)
    row = dict(
        sample=s,
        j5s_total=int(total_j5s),
        j5sd_total=int(total),
        j5d_search=search,
        screen_pass_evals=int(rep["pass_evals"]),
        repair_evals=int(rep["n_evals"]),
        fewer=float(fewer),
        wrong_j5s=int(w_old.sum()),
        wrong_j5sd=int(w_new.sum()),
        wrong_j5sd_wilson=[float(lo), float(hi)],
        fixed=int(fixed),
        broken=int(broken),
        broken_minus_fixed=int(broken - fixed),
        median_err_j5s=float(np.median(e_old)),
        median_err_j5sd=float(np.median(e_new)),
        median_err_higher=float(np.median(e_new) - np.median(e_old)),
        grid_j5=grid_j,
        grid_j5d=grid_d,
        grid_more=grid_d - grid_j,
    )
    row["ok_fewer"] = bool(fewer >= g["g3_fewer_min"])
    row["ok_broken"] = bool(row["broken_minus_fixed"] <= g["g3_broken_minus_fixed_max"])
    row["ok_median"] = bool(row["median_err_higher"] <= g["g3_median_error_higher_max_deg"])
    row["ok_grid"] = bool(row["grid_more"] <= g["g3_grid_more_max"])
    row["ok"] = bool(row["ok_fewer"] and row["ok_broken"] and row["ok_median"] and row["ok_grid"])
    return row


def main() -> int:
    cache = DC.RCACHE
    tab = json.loads((DC.RD / "replay.json").read_text())
    totals = {r["sample"]: int(r["total"]) for r in tab["rows"]}
    out: Dict[str, Any] = dict(
        rule=DR.GATE, g0=json.loads((DC.RD / "g0.json").read_text())["g0_passes"]
    )
    out["g1"] = bool(tab["reproduces_plan"])
    out["g2"] = g2(cache)
    out["g3"] = {s: g3(cache, s, totals[s]) for s in DR.GATE["g3_samples"]}
    out["g3_ok"] = all(r["ok"] for r in out["g3"].values())
    out["passed"] = bool(out["g0"] and out["g1"] and out["g2"]["ok"] and out["g3_ok"])
    (DC.RD / "gate.json").write_text(json.dumps(out, indent=1) + "\n")
    cols = [
        "sample",
        "j5s_total",
        "j5sd_total",
        "fewer",
        "wrong_j5s",
        "wrong_j5sd",
        "fixed",
        "broken",
        "median_err_higher",
        "grid_more",
        "ok",
    ]
    fmt = dict(j5s_total=",d", j5sd_total=",d", fewer=".4f", median_err_higher=".5f")
    DT.write_tables(
        DC.RD / "gate_tables.md",
        dict(ar5_defer_gate_g3=DT.markdown_table(list(out["g3"].values()), cols, fmt)),
    )
    print(json.dumps(out, indent=1))
    return 0 if out["passed"] else 1


if __name__ == "__main__":
    sys.exit(main())
