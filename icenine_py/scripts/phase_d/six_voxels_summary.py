"""Six interior voxels, seed 1, single worker: opt-in C++ overlap rule on (tag cpprule) next to the
stored Python-before (tag mcport) and the C++ columns of mc_parity (tables.md). Writes
benchmarks/phase_d_seed_diag/level0/six_voxels.json and six_voxels_tables.md."""

import json
import re
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ICE = HERE.parents[1]
sys.path.insert(0, str(ICE / "scripts" / "common"))
import doc_tables as DT  # noqa: E402
import gzio as GZ  # noqa: E402

SD = ICE / "benchmarks" / "phase_d_seed_diag"
OUT = SD / "level0"
VOX = [15901, 17034, 18558, 19889, 22867, 5242]


def main() -> None:
    blk = (SD / "mc_parity" / "tables.md").read_text()
    blk = blk.split("<!-- table:mc_parity_stages -->")[1].split("<!-- /table:mc_parity_stages -->")[
        0
    ]
    cpp = {}
    for line in blk.splitlines():
        c = [x.strip().replace(",", "") for x in line.strip("|").split("|")]
        if c and c[0].isdigit():
            cpp[int(c[0])] = dict(total=int(c[8]), l0=int(c[10]), allc=int(c[12]))
    par = json.loads((SD / "mc_parity" / "summary.json").read_text())["voxels"]
    rows, summ = [], {}
    for v in VOX:
        a = GZ.load_json(SD / "runs" / f"cpprule_full_v{v}.json")
        b = GZ.load_json(SD / "runs" / f"mcport_full_v{v}.json")
        n_all = sum(x["n_returned"] for x in a["levels"])
        r = dict(
            voxel=v,
            l0_cpp=cpp[v]["l0"], l0_before=b["levels"][0]["n_returned"],
            l0_rule=a["levels"][0]["n_returned"],
            all_cpp=cpp[v]["allc"], all_before=sum(x["n_returned"] for x in b["levels"]),
            all_rule=n_all,
            evals_cpp=cpp[v]["total"], evals_before=b["evals_total"], evals_rule=a["evals_total"],
            wall_cpp_s=sum(par[str(v)]["cpp_adap_sec"]) / len(par[str(v)]["cpp_adap_sec"]),
            wall_before_s=b["wall_s"], wall_rule_s=a["wall_s"],
            label_rule=a["label"], err_rule=a["err_deg"], err_before=b["err_deg"],
        )  # fmt: skip
        r["evals_ratio_before"] = r["evals_before"] / r["evals_cpp"]
        r["evals_ratio_rule"] = r["evals_rule"] / r["evals_cpp"]
        summ[v] = r
        rows.append(r)
    (OUT / "six_voxels.json").write_text(json.dumps(summ, indent=1, default=float) + "\n")
    t1 = DT.markdown_table(
        rows,
        ["voxel", "l0_cpp", "l0_before", "l0_rule", "all_cpp", "all_before", "all_rule"],
        formats={
            k: ",d" for k in ["l0_cpp", "l0_before", "l0_rule", "all_cpp", "all_before", "all_rule"]
        },
    )
    t2 = DT.markdown_table(
        rows,
        ["voxel", "evals_cpp", "evals_before", "evals_rule", "evals_ratio_before",
         "evals_ratio_rule", "wall_cpp_s", "wall_before_s", "wall_rule_s", "label_rule",
         "err_before", "err_rule"],
        formats={"evals_cpp": ",d", "evals_before": ",d", "evals_rule": ",d",
                 "evals_ratio_before": ".2f", "evals_ratio_rule": ".2f", "wall_cpp_s": ".1f",
                 "wall_before_s": ".0f", "wall_rule_s": ".0f", "err_before": ".4f",
                 "err_rule": ".4f"},
    )  # fmt: skip
    DT.write_tables(
        OUT / "six_voxels_tables.md", {"level0_six_candidates": t1, "level0_six_evals": t2}
    )
    print((OUT / "six_voxels_tables.md").read_text())


if __name__ == "__main__":
    main()
