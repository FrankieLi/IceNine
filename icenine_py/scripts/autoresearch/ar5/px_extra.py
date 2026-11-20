"""AR5 step 27 review fixes: counts that ``px_summary`` adds to ``summary.json`` (no run, no truth).

  uv run python scripts/autoresearch/ar5/px_extra.py

Reads the stored and new J5D logs (``polish_log``, ``polish.json``, ``run.json``) and
``q1_held.json``
and writes ``ar5/px/extra.json``: the grid-search voxels by arm, the seed-test failures and rescues,
the evaluation accounting by stage, the ms per evaluation by stage, and the q1 detail of two held
voxels. Descriptive only: no cause is claimed.
"""

import json
import sys
from pathlib import Path
from typing import Any, Dict, List

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import px_score as PS  # noqa: E402

C, RD = PS.CACHE, PS.RD


def runs(s: str) -> Dict[str, Path]:
    return dict(
        stored=C / f"J5D_{s}", old=C / "live_old" / f"J5D_{s}", new=C / "live_new" / f"J5D_{s}"
    )


def main() -> None:
    out: Dict[str, Any] = {}
    for s in ("mid", "held"):
        d: Dict[str, Any] = {}
        sets: Dict[str, List[int]] = {}
        for k, p in runs(s).items():
            lp = p / "polish_log.json.gz"
            if not lp.exists():  # the stored run's log is not linked; the old rerun equals it
                lp = runs(s)["old"] / "polish_log.json.gz"
            log = PS.read_gz(lp)
            pol = PS.BS.read_json(p / "polish.json")
            st = PS.BS.read_json(p / "run.json")["stats"]
            grid = sorted(int(r["voxel"]) for r in log if r["grid_evals"] > 0)
            fail = sorted(int(r["voxel"]) for r in log)
            sets[k] = grid
            d[k] = dict(
                grid_voxels=grid,
                grid_evals=int(pol["grid_evals"]),
                seed_test_failures=int(st["index_fallback_seed_test"]),
                rescues_attempted=int(pol["rescues_attempted"]),
                rescued=int(pol["rescued"]),
                polished=int(pol["polished"]),
                failing_voxels=fail,
                ms_per_eval=dict(
                    neighbour=1000 * st["wall_neighbor_s"] / st["n_evals_neighbor"],
                    revisit=1000 * st["wall_revisit_s"] / st["n_evals_revisit"],
                ),
            )
        a, b = set(sets["old"]), set(sets["new"])
        d["grid_overlap_old_new"] = dict(
            shared=sorted(a & b), old_only=sorted(a - b), new_only=sorted(b - a)
        )
        fo, fn = set(d["old"]["failing_voxels"]), set(d["new"]["failing_voxels"])
        d["failing_overlap_old_new"] = dict(shared=len(fo & fn), old=len(fo), new=len(fn))
        out[s] = d
    sp = json.loads((RD / "stage_split.json").read_text())
    extra = {s: sp[f"{s}_new"]["total"] - sp[f"{s}_stored"]["total"] for s in ("mid", "held")}
    gridx = {s: sp[f"{s}_new"]["grid"] - sp[f"{s}_stored"]["grid"] for s in ("mid", "held")}
    out["accounting"] = dict(
        extra_evals=extra,
        extra_grid_evals=gridx,
        pooled_extra=sum(extra.values()),
        pooled_extra_grid=sum(gridx.values()),
        held_grid_share=gridx["held"] / extra["held"],
        pooled_grid_share=sum(gridx.values()) / sum(extra.values()),
    )
    q = json.loads((RD / "q1_held.json").read_text())
    rows = {r["voxel"]: r for r in q["rows"]}
    out["q1_held_voxels"] = {
        str(v): dict(base=rows[v]["base"], new=rows[v]["new"]) for v in (19825, 12171) if v in rows
    }
    (RD / "extra.json").write_text(json.dumps(out, indent=1) + "\n")
    print(
        json.dumps(
            {k: v for k, v in out.items() if k in ("accounting", "q1_held_voxels")}, indent=1
        )[:1800]
    )
    for s in ("mid", "held"):
        d = out[s]
        print(
            s,
            d["grid_overlap_old_new"],
            d["failing_overlap_old_new"],
            {
                k: (
                    d[k]["seed_test_failures"],
                    d[k]["rescued"],
                    d[k]["polished"],
                    d[k]["ms_per_eval"],
                )
                for k in ("stored", "old", "new")
            },
        )


if __name__ == "__main__":
    main()
