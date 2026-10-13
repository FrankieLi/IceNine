"""C++ right / wrong on the E0 cases (e0_cpp.py output) next to the Python switch off (stored
faithful-MC run) and on (e0_cpprule.py `on`), paired by voxel (fixed / broken / both wrong)."""

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
import stats as S  # noqa: E402

OUT = ICE / "benchmarks" / "phase_d_seed_diag" / "level0" / "e0"
BASE = ICE / "scripts" / "cma_finisher" / "cache" / "e0_run_mcfaithful"


def main() -> None:
    rows, summ, per_case = [], {}, {}
    for var in C.VARIANTS:
        wc, wo, wn, ids = [], [], [], []
        for f in sorted((HERE / "cache" / "on").glob(f"v*_{var}.npz")):
            cf = HERE / "cache" / "e0_cpp" / f.name.replace(".npz", ".json")
            d = json.loads(cf.read_text())
            z, b = np.load(f), np.load(BASE / f.name)
            Rt = z["R_true"][None]
            Rc = np.array(d["CV"][1:]).reshape(1, 3, 3)
            ec = C.err_deg(Rc, Rt)[0]
            eo = C.err_deg(b["mc_R"][None], Rt)[0]
            en = C.err_deg(z["R"][None], Rt)[0]
            wc.append(ec > C.WRONG_DEG)
            wo.append(eo > C.WRONG_DEG)
            wn.append(en > C.WRONG_DEG)
            ids.append(f.stem)
            per_case[f.stem] = dict(cpp_err=float(ec), py_off_err=float(eo), py_on_err=float(en),
                                    cpp_cost=float(d["CV"][0]))  # fmt: skip
        wc, wo, wn = map(np.array, (wc, wo, wn))
        n = len(wc)
        newly = wn & ~wo
        label = "realistic" if var == "all" else var
        s = dict(
            n=n, cpp_wrong=int(wc.sum()), py_off_wrong=int(wo.sum()), py_on_wrong=int(wn.sum()),
            cpp_vs_off=dict(zip(("cpp_only_wrong", "off_only_wrong", "both_wrong"),
                                (int((wc & ~wo).sum()), int((~wc & wo).sum()),
                                 int((wc & wo).sum())))),
            cpp_vs_on=dict(zip(("cpp_only_wrong", "on_only_wrong", "both_wrong"),
                               (int((wc & ~wn).sum()), int((~wc & wn).sum()),
                                int((wc & wn).sum())))),
            newly_wrong_with_switch=int(newly.sum()),
            cpp_wrong_among_newly=int((wc & newly).sum()),
        )  # fmt: skip
        summ[var] = s
        lo_c, hi_c = S.wilson(s["cpp_wrong"], n)
        rows.append(dict(variant=label, n=n,
                         py_off=f"{s['py_off_wrong']}", py_on=f"{s['py_on_wrong']}",
                         cpp=f"{s['cpp_wrong']} [{lo_c:.3f}, {hi_c:.3f}]",
                         newly_wrong=s["newly_wrong_with_switch"],
                         cpp_wrong_among_newly=s["cpp_wrong_among_newly"],
                         cpp_only_vs_off=s["cpp_vs_off"]["cpp_only_wrong"],
                         off_only_vs_cpp=s["cpp_vs_off"]["off_only_wrong"],
                         cpp_only_vs_on=s["cpp_vs_on"]["cpp_only_wrong"],
                         on_only_vs_cpp=s["cpp_vs_on"]["on_only_wrong"]))  # fmt: skip
    (OUT / "cpp_e0_summary.json").write_text(
        json.dumps(dict(summary=summ, cases=per_case), indent=1)
    )
    t = DT.markdown_table(rows, list(rows[0].keys()))
    DT.write_tables(OUT / "tables_cpp.md", {"e0_cpp_vs_python": t})
    print(t)


if __name__ == "__main__":
    main()
