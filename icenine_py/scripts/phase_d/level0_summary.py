"""Summary tables for the level-0 candidate trace (C++ vs Python, voxel 15901 plus level-0 counts of
the other five interior voxels). Reads benchmarks/phase_d_seed_diag/level0/*.json and the stored
mc_parity table; writes summary.json and tables.md (doc_tables blocks).

  uv run python scripts/phase_d/level0_summary.py
"""

import json
import re
from pathlib import Path
from typing import Any, Dict, List

HERE = Path(__file__).resolve().parent
ICE = HERE.parents[1]
import sys  # noqa: E402

sys.path.insert(0, str(ICE / "scripts" / "common"))
import doc_tables as DT  # noqa: E402
import gzio as GZ  # noqa: E402

OUT = ICE / "benchmarks" / "phase_d_seed_diag" / "level0"
PAR = ICE / "benchmarks" / "phase_d_seed_diag" / "mc_parity"
VOXELS = [15901, 17034, 18558, 19889, 22867, 5242]


def load(name: str) -> Any:
    return GZ.load_json(OUT / name)


def main() -> None:
    summ: Dict[str, Any] = {}
    tables: Dict[str, str] = {}
    cpp = load("cpp_v15901.json")
    variants = [
        ("library default", "", "py_v15901_L0.json", "compare_v15901.json"),
        (
            "C++ overlap rule",
            "_cpprule",
            "py_cpprule_v15901_L0.json",
            "compare_cpprule_v15901.json",
        ),
        (
            "C++ rule + centroid",
            "_cpprule_centroid",
            "py_cpprule_centroid_v15901_L0.json",
            "compare_cpprule_centroid_v15901.json",
        ),
    ]
    # --- S4 settings (library default vs C++, all levels run in the full default run) ---
    full = load("py_v15901.json")
    rows = []
    for lv, (c, p) in enumerate(zip(cpp["levels"], full["levels"])):
        rows.append(
            dict(
                level=lv,
                n_fz_cpp=c["n_fz"],
                n_fz_py=p["n_fz"],
                n_local_cpp=c["n_local"],
                n_local_py=p["n_local"],
                nqmax_cpp=c["max_q"],
                n_recip_cpp=c["n_recip"],
                n_recip_py=p["n_recip"],
                diam_deg_cpp=c["diameter_rad"] * 57.29577951308232,
                diam_deg_py=p["diameter_rad"] * 57.29577951308232,
                n_sym_cpp=c["n_sym"],
                n_sym_py=p["n_sym"],
            )
        )
    summ["S4"] = rows
    tables["level0_settings"] = DT.markdown_table(
        rows,
        list(rows[0].keys()),
        formats=dict(nqmax_cpp=".0f", diam_deg_cpp=".3f", diam_deg_py=".3f"),
    )
    # --- stages S1, S2, S3 at level 0, voxel 15901, per Python variant ---
    srows = []
    for label, _tag, _pyf, cmpf in variants:
        c = load(cmpf)
        s1, s2, s3 = c["S1_screen"], c["S2_recost"], c["S3_filter"]
        s4 = c["S4_levels"][0]
        srows.append(
            dict(
                python=label,
                screened_pairs_cpp=s1["total_pass_cpp"],
                screened_pairs_py=s1["total_pass_py"],
                fz_identical_screen=s1["fz_identical_pass_set"],
                peak_overlap_equal=s1["peak_overlap_equal_frac"],
                recost_median_diff=s2["median"],
                recost_within_1e6=s2["frac_abs_lt_1e_6"],
                cliques_cpp=s4["n_cliques"][0],
                cliques_py=s4["n_cliques"][1],
                kept_cpp=s3["total_kept_cpp"],
                kept_py=s3["total_kept_py"],
            )
        )
    summ["stages_15901"] = srows
    tables["level0_stages"] = DT.markdown_table(
        srows,
        list(srows[0].keys()),
        formats=dict(peak_overlap_equal=".3f", recost_median_diff=".2g", recost_within_1e6=".3f"),
    )
    # --- S3 on C++'s own inputs ---
    s3 = load("s3_on_cpp_inputs_v15901.json")
    summ["S3_on_cpp_inputs"] = {k: v for k, v in s3.items() if k != "examples"}
    tables["level0_s3_cpp_inputs"] = DT.markdown_table(
        [summ["S3_on_cpp_inputs"]],
        ["voxel", "n_cliques_with_ge2", "same_kept_set", "same_kept_order"],
    )
    # --- per-level candidate counts, voxel 15901 ---
    lrows = []
    cfull = load("py_cpprule_centroid_v15901.json")
    for lv, c in enumerate(cpp["levels"]):
        d = full["levels"][lv] if lv < len(full["levels"]) else {}
        e = cfull["levels"][lv] if lv < len(cfull["levels"]) else {}
        lrows.append(
            dict(
                level=lv,
                n_fz_cpp=c["n_fz"],
                cand_cpp=c["n_candidates"],
                n_fz_py=d.get("n_fz"),
                cand_py=d.get("n_candidates"),
                n_fz_probe=e.get("n_fz"),
                cand_probe=e.get("n_candidates"),
            )
        )
    summ["levels_15901"] = lrows
    tables["level0_levels"] = DT.markdown_table(lrows, list(lrows[0].keys()))
    # --- level-0 count over the six interior voxels ---
    stage = (PAR / "tables.md").read_text()
    blk = stage.split("<!-- table:mc_parity_stages -->")[1].split(
        "<!-- /table:mc_parity_stages -->"
    )[0]
    stored: Dict[int, Dict[str, int]] = {}
    for line in blk.splitlines():
        cells = [x.strip() for x in line.strip("|").split("|")]
        if cells and cells[0].isdigit():
            stored[int(cells[0])] = dict(
                py=int(cells[9].replace(",", "")), cpp=int(cells[10].replace(",", ""))
            )
    vrows = []
    for v in VOXELS:
        f = OUT / f"py_cpprule_centroid_v{v}_L0.json"
        if v == 15901:
            f = OUT / "py_cpprule_centroid_v15901.json"
        probe = GZ.load_json(f)["levels"][0]["n_candidates"] if GZ.exists(f) else None
        vrows.append(
            dict(
                voxel=v,
                cpp=stored[v]["cpp"],
                py_default=stored[v]["py"],
                py_probe=probe,
                ratio_default=stored[v]["py"] / stored[v]["cpp"],
                ratio_probe=None if probe is None else probe / stored[v]["cpp"],
            )
        )
    summ["voxels_level0"] = vrows
    tables["level0_voxels"] = DT.markdown_table(
        vrows,
        list(vrows[0].keys()),
        formats=dict(
            cpp=",d", py_default=",d", py_probe=",d", ratio_default=".3f", ratio_probe=".3f"
        ),
    )
    # --- per-peak evidence from the C++ dump (3 FZ points, j = 0, adaptive level 0) ---
    L = (OUT / "cpp_peak_dump_fz0-2.txt").read_text().splitlines()
    cls: Dict[str, List[int]] = {}
    for i, line in enumerate(L):
        if not line.startswith("L0Y omega"):
            continue
        w = [x for x in L[max(0, i - 3) : i] if x.startswith("L0W")]
        if len(w) < 2:
            key = "not on a detector (no projection hit)"
        else:
            br = [int(re.search(r"bright (\d+)", x).group(1)) > 0 for x in w]  # type: ignore
            lit1 = int(re.search(r" in (\d+)", w[1]).group(1)) > 0  # type: ignore
            key = (
                f"det0 bright={br[0]}, det1 lit={lit1}, det1 bright={br[1]}"
                if br[0] or lit1
                else "neither detector lit"
            )
        ov = int(re.search(r"peakOv (\d+)", line).group(1))  # type: ignore
        cls.setdefault(key, []).append(ov)
    prow = [
        dict(pattern=k, n_peaks=len(v), cpp_peak_overlap_sum=sum(v)) for k, v in sorted(cls.items())
    ]
    summ["peak_patterns"] = prow
    tables["level0_peak_patterns"] = DT.markdown_table(
        prow, ["pattern", "n_peaks", "cpp_peak_overlap_sum"]
    )
    (OUT / "summary.json").write_text(json.dumps(summ, indent=1, default=float) + "\n")
    DT.write_tables(OUT / "tables.md", tables)
    print((OUT / "tables.md").read_text())


if __name__ == "__main__":
    main()
