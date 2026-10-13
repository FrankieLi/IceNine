#!/usr/bin/env python3
"""B3 rerun with the C++-faithful MC and VarianceMinimizing: before / after for the rows that use
them (mc_deployed, mc_april, vm_small, default finisher), with the unchanged non-MC rows of the
original B3 cache as the reference.

  bench.py run --cache run_mcfaithful --methods mc_deployed,mc_april,vm_small --workers 10
  uv run python scripts/finisher_bench/mc_faithful_summary.py

Reads cache/run (original B3) and cache/run_mcfaithful (rerun; same cases, same seeds, same
images), writes benchmarks/finisher_bench/mc_faithful/{summary.json,tables.md}. Error, wrong and
the pairing rules are those of finisher_summary.py.
"""

import json
import sys
from pathlib import Path
from typing import Any, Dict, List, Tuple

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import finisher_summary as FS  # noqa: E402
from doc_tables import markdown_table, write_tables  # noqa: E402
from stats import mcnemar_exact, paired_discordant, win_rate  # noqa: E402

OUT = HERE.parents[1] / "benchmarks" / "finisher_bench" / "mc_faithful"
SETS = {"T5-H3": ("H3",), "T5-H0": ("H0",), "SW": ("SW",)}
ROWS: List[Tuple[str, str]] = [
    ("mc_deployed", "natural"),
    ("mc_april", "250"),
    ("mc_april", "1000"),
    ("mc_april", "2600"),
    ("vm_small", "250"),
    ("vm_small", "1000"),
    ("vm_small", "2600"),
    ("vm_small", "10000"),
    ("finisher", "natural"),
]
REF_ROWS: List[Tuple[str, str]] = [("nm", "250"), ("cma_02", "250"), ("cma_02", "1000")]


def row(d: Dict[str, Any], key: str, budget: str) -> Dict[str, np.ndarray]:
    Rt, ct = d["R_true"], d["cost_true"]
    if key == "finisher":
        return dict(
            err=FS.angles(d["fin_R"], Rt), gap=d["fin_cost"] - ct, used=d["fin_evals"].astype(float)
        )
    mi = d["methods"].index(key)
    k = len(d["ckpts"]) - 1 if budget == "natural" else d["ckpts"].index(int(budget))
    return dict(
        err=FS.angles(d["res_R"][:, :, mi, k], Rt),
        gap=d["res_cost"][:, :, mi, k] - ct,
        used=d["res_used"][:, :, mi, k].astype(float),
    )


def cell(e: np.ndarray) -> str:
    n = len(e)
    lt, wr = 100 * (e < 0.02).sum() / n, 100 * (e > FS.WRONG_CUT).sum() / n
    return f"{np.median(e):.4f} / {lt:.0f}% / {wr:.1f}%"


def stats(e: np.ndarray, u: np.ndarray) -> Dict[str, Any]:
    n = len(e)
    return dict(
        n=n, err_median=float(np.median(e)), evals_median=float(np.median(u)),
        lt002=FS.frac(int((e < 0.02).sum()), n), wrong=FS.frac(int((e > FS.WRONG_CUT).sum()), n),
    )  # fmt: skip


def main() -> None:
    old = FS.load(FS.CACHE_DIR / "run")
    new = FS.load(FS.CACHE_DIR / "run_mcfaithful")
    for k in ("kind", "vox", "ri", "j"):
        assert np.array_equal(old[k], new[k]), f"case order differs: {k}"
    assert np.array_equal(old["R_start"], new["R_start"])
    S: Dict[str, Any] = dict(n_cases=int(len(old["kind"])), sets={})
    tables: Dict[str, str] = {}
    for vi, var in enumerate(FS.VARIANTS):
        for sname, kinds in SETS.items():
            mask = np.isin(old["kind"], kinds)
            vox = old["vox"][mask]
            lines: List[Dict[str, Any]] = []
            rec: Dict[str, Any] = dict(n=int(mask.sum()), n_voxels=int(len(np.unique(vox))))
            for key, b in ROWS + REF_ROWS:
                ro = row(old, key, b)
                eo, uo = ro["err"][vi][mask], ro["used"][vi][mask]
                is_ref = (key, b) in REF_ROWS
                if is_ref:
                    lines.append(dict(method=f"{FS.rowname(key, b)} (unchanged)",
                                      before=cell(eo), after=cell(eo), ev_before=np.median(uo),
                                      ev_after=np.median(uo)))  # fmt: skip
                    rec[f"{key}@{b}"] = dict(before=stats(eo, uo))
                    continue
                rn = row(new, key, b)
                en, un = rn["err"][vi][mask], rn["used"][vi][mask]
                wr, tie, _ = win_rate(en, eo, FS.TIE_DEG)
                bw, cw = paired_discordant(en > FS.WRONG_CUT, eo > FS.WRONG_CUT)
                b2, c2 = paired_discordant(en >= 0.02, eo >= 0.02)
                rec[f"{key}@{b}"] = dict(
                    before=stats(eo, uo), after=stats(en, un),
                    after_vs_before=dict(
                        win_rate_after=wr, tie_frac=tie, sign=FS.sign_test(en - eo),
                        voxel_sign=FS.voxel_sign_test(en - eo, vox),
                        mcnemar_wrong=dict(only_after_wrong=bw, only_before_wrong=cw,
                                           p=mcnemar_exact(bw, cw)),
                        mcnemar_lt002=dict(only_after_fails=b2, only_before_fails=c2,
                                           p=mcnemar_exact(b2, c2)),
                    ),
                )  # fmt: skip
                lines.append(dict(method=FS.rowname(key, b), before=cell(eo), after=cell(en),
                                  ev_before=np.median(uo), ev_after=np.median(un)))  # fmt: skip
            S["sets"][f"{sname}|{var}"] = rec
            tables[f"mcf_{sname}_{var}"] = markdown_table(
                lines, ["method", "before", "after", "ev_before", "ev_after"],
                formats={"ev_before": ".0f", "ev_after": ".0f"},
            )  # fmt: skip
    # how the deployed MC ended (before / after): no improvement over the start, ended before its
    # 200-step budget (restarts exhausted), converged
    dep: List[Dict[str, Any]] = []
    for vi, var in enumerate(FS.VARIANTS):
        for sname, kinds in SETS.items():
            mask = np.isin(old["kind"], kinds)
            n = int(mask.sum())
            r: Dict[str, Any] = dict(set=f"{sname} {var}", n=n)
            for tag, d in (("before", old), ("after", new)):
                mi = d["methods"].index("mc_deployed")
                k = len(d["ckpts"]) - 1
                imp = d["res_cost"][vi, :, mi, k] < d["cost_start"][vi] - 1e-12
                early = d["total_evals"][vi, :, mi] < 200
                ni, ee = (~imp) & mask, early & mask
                fn, fe = FS.frac(int(ni.sum()), n), FS.frac(int(ee.sum()), n)
                r[f"noimp_{tag}"] = FS.fmt_frac(fn)
                r[f"early_{tag}"] = FS.fmt_frac(fe)
                r[f"both_{tag}"] = int((ni & ee).sum())
                r[f"evals_med_{tag}"] = float(np.median(d["total_evals"][vi, mask, mi]))
            dep.append(r)
    S["deployed_mc"] = dep
    tables["mcf_deployed"] = markdown_table(
        dep, ["set", "n", "noimp_before", "noimp_after", "early_before", "early_after",
              "both_before", "both_after", "evals_med_before", "evals_med_after"],
        formats={"evals_med_before": ".0f", "evals_med_after": ".0f"},
    )  # fmt: skip
    # bit-identity of the default finisher to the stored T5 result is not expected any more
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "summary.json").write_text(json.dumps(S, indent=1, default=float))
    write_tables(OUT / "tables.md", tables)
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
