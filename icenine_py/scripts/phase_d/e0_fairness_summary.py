"""Step 4 summary: Python (min_sin_eta 0.3 as E0, and 0), C++ rule off / on, and C++ with three
random seeds, on the 400 E0 cases; Python-on against C++ on the voxels the switch broke.
No significance tests: counts and paired splits (fixed / broken / both wrong) by voxel."""

import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ICE = HERE.parents[1]
sys.path.insert(0, str(ICE / "scripts" / "findoptimal_robustness"))
sys.path.insert(0, str(ICE / "scripts" / "common"))
import common as C  # noqa: E402
import doc_tables as DT  # noqa: E402

CACHE = HERE / "cache"
BASE = ICE / "scripts" / "cma_finisher" / "cache" / "e0_run_mcfaithful"
OUT = ICE / "benchmarks" / "phase_d_seed_diag" / "level0" / "e0"
SEEDS = ["", "_s1", "_s2"]  # "" = the default mt19937 seed 5489


def main() -> None:
    names = sorted(p.name for p in (CACHE / "on").glob("v*.npz"))
    Rt = {n: np.load(CACHE / "on" / n)["R_true"] for n in names}

    def err_py(path: Path, key: str, n: str) -> float:
        return float(C.err_deg(np.load(path)[key][None], Rt[n][None])[0])

    wrong = {}  # arm -> {name: bool}
    for n in names:
        for arm, path, key in (
            ("py_off", BASE / n, "mc_R"),
            ("py_on", CACHE / "on" / n, "R"),
            ("py_off_ms0", CACHE / "off_ms0" / n, "R"),
            ("py_on_ms0", CACHE / "on_ms0" / n, "R"),
        ):
            wrong.setdefault(arm, {})[n] = err_py(path, key, n) > C.WRONG_DEG
        for sd in SEEDS:
            cf = CACHE / f"e0_cpp{sd}" / n.replace(".npz", ".json")
            d = json.loads(cf.read_text())
            R = np.array(d["CV"][1:]).reshape(1, 3, 3)
            wrong.setdefault(f"cpp{sd or '_s0'}", {})[n] = (
                float(C.err_deg(R, Rt[n][None])[0]) > C.WRONG_DEG
            )
    arms = list(wrong)
    rows, summ = [], {}
    for var in C.VARIANTS:
        sel = [n for n in names if n.endswith(f"_{var}.npz")]
        label = "realistic" if var == "all" else var
        row = dict(variant=label, n=len(sel))
        for a in arms:
            row[a] = int(sum(wrong[a][n] for n in sel))
        cpp_arms = [f"cpp{sd or '_s0'}" for sd in SEEDS]
        row["cpp_mean"] = float(np.mean([row[a] for a in cpp_arms]))
        row["cpp_wrong_all_seeds"] = int(sum(all(wrong[a][n] for a in cpp_arms) for n in sel))
        row["cpp_wrong_any_seed"] = int(sum(any(wrong[a][n] for a in cpp_arms) for n in sel))
        rows.append(row)
        summ[label] = row
    # voxels the switch broke, for both min_sin_eta settings, and C++ on them
    broke_rows = []
    for tag, off, on in (
        ("min_sin_eta 0.3", "py_off", "py_on"),
        ("min_sin_eta 0", "py_off_ms0", "py_on_ms0"),
    ):
        for var in C.VARIANTS:
            sel = [n for n in names if n.endswith(f"_{var}.npz")]
            newly = [n for n in sel if wrong[on][n] and not wrong[off][n]]
            r = dict(setting=tag, variant="realistic" if var == "all" else var, n=len(sel),
                     off_wrong=sum(wrong[off][n] for n in sel),
                     on_wrong=sum(wrong[on][n] for n in sel),
                     broken=len(newly), fixed=sum(wrong[off][n] and not wrong[on][n] for n in sel),
                     both_wrong=sum(wrong[off][n] and wrong[on][n] for n in sel))  # fmt: skip
            for sd in SEEDS:
                a = f"cpp{sd or '_s0'}"
                r[f"{a}_wrong_on_broken"] = sum(wrong[a][n] for n in newly)
            broke_rows.append(r)
    summ["broken"] = broke_rows
    (OUT / "fairness_summary.json").write_text(json.dumps(summ, indent=1, default=float))
    t1 = DT.markdown_table(rows, list(rows[0].keys()), formats={"cpp_mean": ".1f"})
    t2 = DT.markdown_table(broke_rows, list(broke_rows[0].keys()))
    DT.write_tables(
        OUT / "tables_fairness.md", {"e0_fairness_counts": t1, "e0_fairness_broken": t2}
    )
    print(t1)
    print(t2)


if __name__ == "__main__":
    main()
