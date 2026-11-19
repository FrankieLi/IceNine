"""AR5 step 24, information only, after the gate failed (decision D-20261119-N4): the per-stage
evaluations of ``J5SD`` (step 22's stored arm) and ``J5SDB`` (live) on ``mid`` and ``held``, with
the grid searches and seed-test failures, so the reader can see where the shortfall on ``mid`` sits.

  uv run python scripts/autoresearch/ar5/bf_info.py

Reads the runs of ``bf_run_all.sh --stage gate`` and step 22's stored arms (links); arithmetic only.
Writes ``ar5/bf/info.json`` and ``info_tables.md``. It does not change any verdict.
"""

import json
import sys
from pathlib import Path
from typing import Any, Dict

sys.path.insert(0, str(Path(__file__).resolve().parent))
import bf_score as BS  # noqa: E402
import doc_tables as DT  # noqa: E402

RD = BS.ICE / "benchmarks" / "autoresearch" / "ar5" / "bf"
KEYS = [
    "index",
    "rescue",
    "polish",
    "grid",
    "bfs_neighbour",
    "bfs_revisit",
    "screen_pass",
    "repair",
]


def main() -> int:
    rows: Dict[str, Any] = {}
    tab = []
    tot = []
    for s in ("mid", "held"):
        r = BS.one(s)
        rows[s] = r
        for arm, key, st in (
            ("J5SD", "j5sd", r["stages_j5sd"]),
            ("J5SDB", "j5sdb", r["stages"]),
        ):
            d = dict(sample=s, arm=arm, total=r["totals"][key], **{k: st[k] for k in KEYS})
            tab.append(d)
        diff = dict(sample=s, arm="J5SDB - J5SD", total=r["totals"]["j5sdb"] - r["totals"]["j5sd"])
        for k in KEYS:
            diff[k] = r["stages"][k] - r["stages_j5sd"][k]
        tab.append(diff)
        tot.append(
            dict(
                sample=s,
                fewer=100 * r["fewer"]["j5sd"],
                wrong=f"{r['wrong']['j5sd']['n']} / {r['wrong']['j5sdb']['n']}",
                fixed_broken=f"{r['pair']['j5sd']['fixed']} / {r['pair']['j5sd']['broken']}",
                median_higher=r["pair"]["j5sd"]["median_higher"],
                grid=f"{r['grid']['j5d']} / {r['grid']['j5db']}",
                seed_test=f"{r['seed_test']['j5d']} / {r['seed_test']['j5db']}",
                seeds_steps=f"{r['bf']['steps']}",
            )
        )
    fmt = {k: ",d" for k in ["total"] + KEYS}
    DT.write_tables(
        RD / "info_tables.md",
        dict(
            ar5_bf_info_stages=DT.markdown_table(tab, ["sample", "arm", "total"] + KEYS, fmt),
            ar5_bf_info_counts=DT.markdown_table(
                tot,
                [
                    "sample",
                    "fewer",
                    "wrong",
                    "fixed_broken",
                    "median_higher",
                    "grid",
                    "seed_test",
                    "seeds_steps",
                ],
                dict(fewer=".2f", median_higher=".5f"),
            ),
        ),
    )
    (RD / "info.json").write_text(json.dumps(rows, indent=1, default=float) + "\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
