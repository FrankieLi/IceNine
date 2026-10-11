"""Merge the level-0 tables (trace, six voxels, E0 on / rule, cost offset) into one marker file for
scripts/dev/sync_doc_tables.py: benchmarks/phase_d_seed_diag/level0/doc_tables_all.md."""

import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ICE = HERE.parents[1]
sys.path.insert(0, str(ICE / "scripts" / "common"))
import doc_tables as DT  # noqa: E402

OUT = ICE / "benchmarks" / "phase_d_seed_diag" / "level0"


def main() -> None:
    parts = []
    for f in [
        "tables.md",
        "six_voxels_tables.md",
        "e0/tables_on.md",
        "e0/tables_rule.md",
        "e0/tables_cpp.md",
        "e0/tables_debug.md",
        "e0/tables_debug_levels.md",
    ]:
        parts.append((OUT / f).read_text().rstrip("\n") + "\n")
    c = json.loads((OUT / "cost_offset_check.json").read_text())
    rows = [
        dict(rule=k, **{x: v[x] for x in ("n_within_1e6", "n_within_1e3", "mean_signed",
                                          "median_signed", "n_python_lower", "n_python_higher",
                                          "corr")})
        for k, v in (("off (library default)", c["rule_off"]), ("on", c["rule_on"]))
    ]  # fmt: skip
    t = DT.markdown_table(
        rows,
        ["rule", "n_within_1e6", "n_within_1e3", "mean_signed", "median_signed",
         "n_python_lower", "n_python_higher", "corr"],
        formats={"mean_signed": ".5f", "median_signed": ".2e", "corr": ".5f"},
    )  # fmt: skip
    parts.append(f"<!-- table:level0_cost_offset -->\n{t}\n<!-- /table:level0_cost_offset -->\n")
    (OUT / "doc_tables_all.md").write_text("\n".join(parts))
    print("ok")


if __name__ == "__main__":
    main()
