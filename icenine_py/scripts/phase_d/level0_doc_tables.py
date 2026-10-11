"""Merge the level-0 tables (trace, six voxels, E0 on / rule, cost offset) into one marker file for
scripts/dev/sync_doc_tables.py: benchmarks/phase_d_seed_diag/level0/doc_tables_all.md."""

import json

import numpy as np
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
        "e0/tables_fairness.md",
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
    cc = json.loads((OUT / "e0" / "cost_check.json").read_text())
    crow = []
    for msin in (0.3, 0.0):
        for rule in (False, True):
            sel = [x for x in cc if x["py_min_sin_eta"] == msin and x["py_rule"] == rule]
            dq = np.array([abs(x["py_quality"] - x["cpp_quality"]) for x in sel])
            dc = np.array([abs(x["py_confidence"] - x["cpp_confidence"]) for x in sel])
            same = sum(
                x["py_peak_overlap"] == x["cpp_peak_overlap"]
                and x["py_peak_on_det"] == x["cpp_peak_on_det"]
                for x in sel
            )
            crow.append(dict(py_min_sin_eta=msin, py_rule="on" if rule else "off", n=len(sel),
                             median_abs_dq=float(np.median(dq)), max_abs_dq=float(dq.max()),
                             median_abs_dconf=float(np.median(dc)), max_abs_dconf=float(dc.max()),
                             same_peak_counts=same))  # fmt: skip
    t = DT.markdown_table(
        crow, list(crow[0].keys()),
        formats={"median_abs_dq": ".4f", "max_abs_dq": ".4f", "median_abs_dconf": ".4f",
                 "max_abs_dconf": ".4f", "py_min_sin_eta": ".1f"},
    )  # fmt: skip
    parts.append(f"<!-- table:e0_cost_check -->\n{t}\n<!-- /table:e0_cost_check -->\n")
    (OUT / "doc_tables_all.md").write_text("\n".join(parts))
    print("ok")


if __name__ == "__main__":
    main()
