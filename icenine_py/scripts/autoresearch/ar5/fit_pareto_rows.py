"""AR5 step 2: add the rows of the finished S10 runs to the pareto rows file.

  uv run python scripts/autoresearch/ar5/fit_pareto_rows.py

The file is ``docs/autoresearch/pareto/pareto_rows.json``. Two rows per finished run and sample
(before the pass: ``boundary_pass`` false, the evaluations without the pass; after: with the
pass), in the format of the S8 rows. Rows with ``step`` "AR5 step 2" are replaced, so the script
can be run again. Reads ``part2_<run>_<sample>.json`` only; an unfinished run gets no row (it is
listed under ``excluded`` instead).
"""

import json
from pathlib import Path
from typing import Any, Dict, List

HERE = Path(__file__).resolve().parent
ICE = HERE.parents[2]
FIT = ICE / "benchmarks" / "autoresearch" / "ar5" / "fit"
ROWS = ICE / "docs" / "autoresearch" / "pareto" / "pareto_rows.json"
STEP = "AR5 step 2"
SAMPLE_ID = {"big": "big_24570", "grad": "gradient_24570", "mid": "mid_24570"}
SAMPLE_TEXT = {
    "big": "the sample with larger grains",
    "grad": "the first gradient sample",
    "mid": "mid",
}
METHOD = {
    "S10": "BFS indexed, best 5 of 20 seeds, stop at hit ratio 0.95, CMA fits without the MC "
    "fallback",
    "S10a": "BFS indexed, best 5 of 20 seeds, stop at hit ratio 0.95, CMA seed refinement",
    "S10b": "BFS indexed, best 5 of 20 seeds, stop at hit ratio 0.95, CMA fits without the MC "
    "fallback",
}


def commas(n: int) -> str:
    return f"{n:,}"


def rows_for(run: str, sample: str, d: Dict[str, Any]) -> List[Dict[str, Any]]:
    st = d["stages"]
    n = 24570
    out: List[Dict[str, Any]] = []
    inc_b = f"seed step {commas(st['seed_step'])}, "
    if st["grid_fallback"]:
        inc_b += f"grid-search fallbacks {commas(st['grid_fallback'])}, "
    inc_b += (
        f"neighbour local optimizations {commas(st['neighbour'])}, "
        f"revisit local optimizations {commas(st['revisit'])}"
    )
    flag = (
        f"run {run}: run S8 with no MC fallback after a CMA neighbour or revisit fit (a fit "
        f"returns what its CMA run returned) (on {SAMPLE_TEXT[sample]}); "
        f"{d['grid_fallbacks']} grid-search fallbacks; no seed-refinement arm met the plan's "
        "1% rule, so the seed step is S8's; script layer, no library key; "
        "the change was chosen from stored counts before the run"
    )
    for after in (False, True):
        side = d["after" if after else "before"]
        total = st["total"] if after else st["before_pass"]
        inc = inc_b + (f", boundary pass {commas(st['pass_'])}" if after else "")
        out.append(
            {
                "sample": SAMPLE_ID[sample],
                "step": STEP,
                "method": METHOD[run] + (" + boundary pass" if after else "") + f" (run {run})",
                "family": "BFS with indexed seeds",
                "boundary_pass": after,
                "evals_total": total,
                "evals_per_voxel": round(total / n, 2),
                "per_voxel_derived": True,
                "wrong": side["wrong"],
                "n": n,
                "interval_pct": [
                    round(100 * side["wilson_lo"], 2),
                    round(100 * side["wilson_hi"], 2),
                ],
                "estimated": False,
                "what_evals_include": inc,
                "source_json": (
                    f"icenine_py/benchmarks/autoresearch/ar5/fit/part2_{run}_{sample}.json :: "
                    f"stages, {'after' if after else 'before'}"
                ),
                "source_page": "icenine_py/docs/autoresearch/findings_ar5_fits.md",
                "flag": flag + ("" if after else "; derived: evals_total = total minus the pass"),
            }
        )
    return out


def main() -> None:
    doc = json.loads(ROWS.read_text())
    doc["rows"] = [r for r in doc["rows"] if r.get("step") != STEP]
    doc["excluded"] = [e for e in doc["excluded"] if not e["what"].startswith(STEP)]
    added: List[str] = []
    for run in ("S10", "S10a", "S10b"):
        for sample in ("big", "grad", "mid"):
            p = FIT / f"part2_{run}_{sample}.json"
            if not p.exists():
                continue
            d = json.loads(p.read_text())
            if d.get("unfinished"):
                doc["excluded"].append(
                    dict(
                        what=f"{STEP} run {run} on {SAMPLE_TEXT[sample]}",
                        why="reached the 90 minute cap; no finished map and no pass",
                    )
                )
                continue
            doc["rows"] += rows_for(run, sample, d)
            added.append(f"{run} {sample}")
    ROWS.write_text(json.dumps(doc, indent=1, ensure_ascii=False) + "\n")
    print("rows added for:", ", ".join(added) or "none")


if __name__ == "__main__":
    main()
