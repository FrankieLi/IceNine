"""AR3 step 11: add the two rows of run S8 on ``held`` (before and after the pass) to the pareto
rows file, in the format of the S8 rows, and the sample entry ``held_24570``.

  uv run python scripts/autoresearch/ar3/held_pareto_rows.py

Reads ``benchmarks/autoresearch/ar3/held/summary.json`` only. Rows with ``step`` "AR3 step 11" are
replaced, so the script can be run again. A run at its cap gets no row (it is listed under
``excluded``).
"""

import json
from pathlib import Path
from typing import Any, Dict, List

HERE = Path(__file__).resolve().parent
ICE = HERE.parents[2]
SUMMARY = ICE / "benchmarks" / "autoresearch" / "ar3" / "held" / "summary.json"
ROWS = ICE / "docs" / "autoresearch" / "pareto" / "pareto_rows.json"
STEP = "AR3 step 11"
SAMPLE_ID = "held_24570"
N = 24570


def commas(n: int) -> str:
    return f"{n:,}"


def sample_entry(d: Dict[str, Any]) -> Dict[str, Any]:
    nl = d["sample"]["n_by_level"]
    levels = ", ".join(f"{nl[k]}" for k in ("0.0", "0.2", "0.5", "1.0", "2.0", "4.0"))
    return {
        "id": SAMPLE_ID,
        "label": "Whole layer, 250 grains with orientation gradients (sample held): 24,570 voxels",
        "n_voxels": N,
        "wrong_definition": (
            "more than 1 degree from the true orientation, after crystal symmetry "
            "(icenine_py/docs/autoresearch/findings_ar3_gradient.md:14)"
        ),
        "note": (
            "AR3 step 11. 250 grains from a Voronoi grain map on the same 24,570-voxel grid, built "
            "like mid with generator seed 1; orientations drawn so that adjacent grains differ by "
            "at least 5 degrees; orientation gradients at six levels (0, 0.2, 0.5, 1, 2 and 4 "
            f"degrees on {levels} grains); clean images. Held out: no setting of any method was "
            "chosen on this sample. In AR3 step 11 only run S8 was run on it; runs of later steps "
            "are in their own rows (AR5 step 4, D-20261111-N18)."
        ),
    }


def rows_for(d: Dict[str, Any]) -> List[Dict[str, Any]]:
    st = d["held"]["stages"]
    nfb = d["held"]["grid_fallbacks"]
    flag = (
        "run S8 on the held-out sample held: run S5 whose seed step stops at the first refinement "
        f"with a hit ratio of at least 0.95; {nfb} grid-search fallbacks; no setting was chosen "
        "on this sample; this is the baseline for later whole-layer tests on held; second run, "
        "cap 480 minutes (D-20261111-N17); the first run reached its 120 minute cap"
    )
    inc_b = f"seed step {commas(st['seed_step'])}, "
    if st["grid_fallback"]:
        inc_b += f"grid-search fallbacks {commas(st['grid_fallback'])}, "
    inc_b += (
        f"neighbour local optimizations {commas(st['neighbour'])}, "
        f"revisit local optimizations {commas(st['revisit'])}"
    )
    out: List[Dict[str, Any]] = []
    for after in (False, True):
        side = d["held"]["after" if after else "before"]
        total = st["total"] if after else st["before_pass"]
        out.append(
            {
                "sample": SAMPLE_ID,
                "step": STEP,
                "method": "BFS indexed, best 5 of 20 seeds, stop at hit ratio 0.95, CMA"
                + (" + boundary pass" if after else "")
                + " (run S8)",
                "family": "BFS with indexed seeds",
                "boundary_pass": after,
                "evals_total": total,
                "evals_per_voxel": round(total / N, 2),
                "per_voxel_derived": True,
                "wrong": side["wrong"],
                "n": N,
                "interval_pct": [
                    round(100 * side["wilson_lo"], 2),
                    round(100 * side["wilson_hi"], 2),
                ],
                "estimated": False,
                "what_evals_include": inc_b
                + (f", boundary pass {commas(st['pass_'])}" if after else ""),
                "source_json": (
                    "icenine_py/benchmarks/autoresearch/ar3/held/summary.json :: "
                    f"held.{{stages, {'after' if after else 'before'}}}"
                ),
                "source_page": "icenine_py/docs/autoresearch/findings_ar3_held.md",
                "flag": flag + ("" if after else "; derived: evals_total = total minus the pass"),
            }
        )
    return out


def main() -> None:
    d = json.loads(SUMMARY.read_text())
    doc = json.loads(ROWS.read_text())
    doc["rows"] = [r for r in doc["rows"] if r.get("step") != STEP]
    doc["excluded"] = [e for e in doc["excluded"] if not e["what"].startswith(STEP)]
    doc["samples"] = [s for s in doc["samples"] if s["id"] != SAMPLE_ID]
    if d.get("unfinished"):
        doc["excluded"].append(
            dict(
                what=f"{STEP} run S8 on held",
                why="reached the 120 minute cap; no finished map and no pass",
            )
        )
        print("unfinished: no rows")
    else:
        doc["samples"].append(sample_entry(d))
        doc["rows"] += rows_for(d)
        print("rows added: S8 held before and after the pass")
    ROWS.write_text(json.dumps(doc, indent=1, ensure_ascii=False) + "\n")


if __name__ == "__main__":
    main()
