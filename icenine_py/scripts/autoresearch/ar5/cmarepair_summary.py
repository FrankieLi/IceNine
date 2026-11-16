"""AR5 step 20, phase A step 5 (after the held4 marker): validity, statistics, decision rules,
predictions, tables, plot and pareto rows of arm J5C (J5 plus the CMA repair).

  uv run python scripts/autoresearch/ar5/cmarepair_summary.py

Reads the five disclosed samples' stored numbers (``stats_<s>.json`` of step 16, ``stats_mid.json``
of step 18 for J5L, this step's ``gate.json`` rows, step 19's ``summary.json`` for J5R) and, for
``held4``, its runs (``S8_held4``, ``J5_held4``, ``passes_held4``, ``J5C_held4``), readable only
once ``cache/logs/cmarepair_runs.done`` exists (``cmarepair_guard``). Truth is used for scoring
only. Writes ``stats_held4.json``, ``summary.json``, ``tables.md``, ``pareto_rows.json`` and
``ar5_cmarepair.png`` in ``benchmarks/autoresearch/ar5/cmarepair``. No p-values.

Rules are ``cmarepair_rules.JUDGE`` as frozen. Evaluations: search = the run before the pass, pass
= the boundary pass, repair = the CMA fits; total = their sum. "Right" = at most 1 degree from the
truth (the house threshold).
"""

import hashlib
import json
import subprocess
import sys
from pathlib import Path
from typing import Any, Dict, List

import cmarepair_data as RD
import cmarepair_freeze as CF
import cmarepair_guard as CG
import cmarepair_rules as CR
import numpy as np

HERE = Path(__file__).resolve().parent
ICE = HERE.parents[2]
sys.path.insert(0, str(ICE / "scripts" / "common"))
import doc_tables as DT  # noqa: E402
import stats as ST  # noqa: E402

OUT = CG.RD
C = CG.cache_dir()
N = 24570
J = CR.JUDGE
ORDER = ("held4",) + CR.SAMPLES


def stage_evals(run: Dict[str, Any], pass_evals: int) -> Dict[str, int]:
    st = run["stats"]
    d = dict(
        seed_step=int(run["pick"]["evals"]),
        grid_fallback=int(st.get("index_fallback_grid_evals", 0)),
        neighbour=int(st["n_evals_neighbor"]),
        revisit=int(st["n_evals_revisit"]),
        pass_=int(pass_evals),
    )
    assert d["seed_step"] + d["grid_fallback"] + d["neighbour"] + d["revisit"] == int(
        run["evals_total"]
    )
    d["before_pass"] = int(run["evals_total"])
    d["total"] = d["before_pass"] + d["pass_"]
    return d


def side(err: np.ndarray) -> Dict[str, Any]:
    k = int((err > CR.WRONG_DEG).sum())
    lo, hi = ST.wilson(k, len(err))
    return dict(wrong=k, wilson_lo=lo, wilson_hi=hi, median_err_deg=float(np.median(err)))


# ---- 1. validity -------------------------------------------------------------------------------


def validity() -> Dict[str, Any]:
    marker = C / CG.MARKER
    v: Dict[str, Any] = dict(marker=marker.read_text().splitlines()[0])
    for arm in ("S8", "J5"):
        rd = C / f"{arm}_held4"
        run = json.loads(CG.guarded(rd / "run.json").read_text())
        v[f"{arm}_n_processed"] = int(run["n_processed"])
        v[f"{arm}_capped"] = (rd / "capped.json").exists()
        v[f"{arm}_pass_json"] = (C / "passes_held4" / f"D4_{arm}" / "pass.json").exists()
        v[f"{arm}_done"] = (C / "logs" / f"{arm}_held4.done").read_text().strip()
    rep = json.loads(CG.guarded(C / "J5C_held4" / "repair.json").read_text())
    v["repair_done"] = (C / "logs" / "repair_held4.done").read_text().strip()
    v["repair_fits_all_121"] = bool(rep["median_evals_per_fit"] == CR.EVALS_PER_FIT)
    v["freeze_check"] = CF.problems() or "ok"
    t = subprocess.run(
        ["git", "log", "-1", "--format=%ct", "--", str(CF.FROZEN)],
        cwd=str(ICE),
        capture_output=True,
        text=True,
    ).stdout.strip()
    froz = float(t)
    files: List[Path] = []
    for p in (
        C / "S8_held4",
        C / "J5_held4",
        C / "passes_held4",
        C / "J5C_held4",
        ICE / "benchmarks" / "autoresearch" / "ar3" / "held4",
    ):
        files += [f for f in p.rglob("*") if f.is_file()]
    mt = [f.stat().st_mtime for f in files]
    v["freeze_commit_epoch"] = froz
    v["n_held4_output_files"] = len(files)
    v["oldest_held4_output_minus_freeze_s"] = float(min(mt) - froz)
    v["freeze_older_than_every_held4_output"] = bool(min(mt) > froz)
    fz = json.loads(CF.FROZEN.read_text())
    v["frozen_hash_count"] = len(fz["files"])
    v["frozen_sha_of_frozen_json"] = hashlib.sha256(CF.FROZEN.read_bytes()).hexdigest()[:16]
    v["ok"] = bool(
        v["marker"] == "exit=0"
        and v["S8_n_processed"] == N
        and v["J5_n_processed"] == N
        and not v["S8_capped"]
        and not v["J5_capped"]
        and v["S8_pass_json"]
        and v["J5_pass_json"]
        and v["S8_done"] == "run=0 pass=0"
        and v["J5_done"] == "run=0 pass=0"
        and v["repair_done"] == "exit=0"
        and v["freeze_check"] == "ok"
        and v["freeze_older_than_every_held4_output"]
    )
    return v


# ---- 2. statistics -----------------------------------------------------------------------------


def held4_rows() -> Dict[str, Any]:
    y = RD.sym()
    tr = RD.read_truth("held4")
    Rt = tr["R_t"]
    arms: Dict[str, Any] = {}
    for a in ("S8", "J5"):
        run = json.loads(CG.guarded(C / f"{a}_held4" / "run.json").read_text())
        pj = json.loads(CG.guarded(C / "passes_held4" / f"D4_{a}" / "pass.json").read_text())
        R = np.load(CG.guarded(C / "passes_held4" / f"after_{a}.npy")).astype(np.float64)
        err = np.asarray(CR.misorientation_deg(R, Rt, y), float)
        arms[a] = dict(stages=stage_evals(run, int(pj["evals"])), after=side(err), err=err, R=R)
    rep = json.loads(CG.guarded(C / "J5C_held4" / "repair.json").read_text())
    Rc = np.load(CG.guarded(C / "J5C_held4" / "after.npy"))
    errc = np.asarray(CR.misorientation_deg(Rc, Rt, y), float)
    w = {k: arms[k]["err"] > CR.WRONG_DEG for k in arms}
    wc = errc > CR.WRONG_DEG
    f_c_j5, b_c_j5 = CR.fixed_broken(w["J5"], wc)
    f_c_s8, b_c_s8 = CR.fixed_broken(w["S8"], wc)
    f_j5_s8, b_j5_s8 = CR.fixed_broken(w["S8"], w["J5"])
    j5t = arms["J5"]["stages"]["total"]
    row = dict(
        sample="held4",
        flagged=rep["n_flagged"],
        fits=rep["n_fits"],
        changed=rep["n_changed"],
        repair_evals=rep["n_evals"],
        median_evals_per_fit=rep["median_evals_per_fit"],
        s8_search=arms["S8"]["stages"]["before_pass"],
        s8_pass=arms["S8"]["stages"]["pass_"],
        s8_total=arms["S8"]["stages"]["total"],
        j5_search=arms["J5"]["stages"]["before_pass"],
        j5_pass=arms["J5"]["stages"]["pass_"],
        j5_total=j5t,
        wrong_S8=arms["S8"]["after"]["wrong"],
        wrong_J5=arms["J5"]["after"]["wrong"],
        wrong_J5C=int(wc.sum()),
        wilson_S8=[arms["S8"]["after"]["wilson_lo"], arms["S8"]["after"]["wilson_hi"]],
        wilson_J5=[arms["J5"]["after"]["wilson_lo"], arms["J5"]["after"]["wilson_hi"]],
        wilson_J5C=list(ST.wilson(int(wc.sum()), N)),
        median_S8=arms["S8"]["after"]["median_err_deg"],
        median_J5=arms["J5"]["after"]["median_err_deg"],
        median_J5C=float(np.median(errc)),
        J5C_vs_J5=dict(fixed=f_c_j5, broken=b_c_j5),
        J5C_vs_S8=dict(fixed=f_c_s8, broken=b_c_s8, net=b_c_s8 - f_c_s8),
        J5_vs_S8=dict(fixed=f_j5_s8, broken=b_j5_s8, net=b_j5_s8 - f_j5_s8),
        both_wrong_J5C_vs_J5=int((w["J5"] & wc).sum()),
        share_pct=100.0 * rep["n_evals"] / j5t,
    )
    return row


def disclosed_rows() -> Dict[str, Any]:
    gate = json.loads((CG.RD / "gate.json").read_text())["rows"]
    rows: Dict[str, Any] = {}
    for s in CR.SAMPLES:
        g = gate[s]
        st = json.loads((RD.STATS16 / f"stats_{s}.json").read_text())["arms"]
        rows[s] = dict(
            sample=s,
            flagged=g["flagged"],
            fits=g["fits"],
            changed=g["changed"],
            repair_evals=g["repair_evals"],
            median_evals_per_fit=g["median_evals_per_fit"],
            s8_search=st["S8"]["stages"]["before_pass"],
            s8_pass=st["S8"]["stages"]["pass_"],
            s8_total=st["S8"]["stages"]["total"],
            j5_search=st["J5"]["stages"]["before_pass"],
            j5_pass=st["J5"]["stages"]["pass_"],
            j5_total=g["j5_total"],
            wrong_S8=g["wrong_S8"],
            wrong_J5=g["wrong_J5"],
            wrong_J5C=g["wrong_J5C"],
            wilson_S8=[st["S8"]["after"]["wilson_lo"], st["S8"]["after"]["wilson_hi"]],
            wilson_J5=[st["J5"]["after"]["wilson_lo"], st["J5"]["after"]["wilson_hi"]],
            wilson_J5C=list(ST.wilson(g["wrong_J5C"], N)),
            median_S8=st["S8"]["after"]["median_err_deg"],
            median_J5=g["median_err_J5"],
            median_J5C=g["median_err_J5C"],
            J5C_vs_J5=dict(fixed=g["fixed_vs_J5"], broken=g["broken_vs_J5"]),
            J5C_vs_S8=g["J5C_vs_S8"],
            J5_vs_S8=g["J5_vs_S8"],
            share_pct=g["share_pct"],
            in_sample=True,
        )
    return rows


def judge(r: Dict[str, Any]) -> Dict[str, Any]:
    s = r["sample"]
    cap = J["max_share_pct_grad"] if s == "grad" else J["max_share_pct"]
    tot = r["j5_total"] + r["repair_evals"]
    cut = 100.0 * (1 - tot / r["s8_total"])
    out = dict(
        j5c_total=tot,
        j5c_search=r["j5_search"],
        j5c_pass=r["j5_pass"],
        r1=bool(r["share_pct"] <= cap),
        r1_value=r["share_pct"],
        r1_cap=cap,
        r2=bool(r["J5C_vs_J5"]["broken"] - r["J5C_vs_J5"]["fixed"] <= J["max_net_broken_vs_j5"]),
        r2_value=r["J5C_vs_J5"]["broken"] - r["J5C_vs_J5"]["fixed"],
        r3=bool(r["median_J5C"] - r["median_J5"] <= J["max_median_rise_vs_j5_deg"]),
        r3_value=r["median_J5C"] - r["median_J5"],
        cut_vs_s8_pct=cut,
        r4_cut=bool(cut >= J["min_cut_vs_s8_pct"]),
        r4_net=bool(r["J5C_vs_S8"]["net"] <= J["max_net_broken_vs_s8"]),
        r4_median=bool(r["median_J5C"] - r["median_S8"] <= J["max_median_rise_vs_s8_deg"]),
        r4_median_value=r["median_J5C"] - r["median_S8"],
    )
    out["r4"] = bool(out["r4_cut"] and out["r4_net"] and out["r4_median"])
    out["pass"] = bool(out["r1"] and out["r2"] and out["r3"] and out["r4"])
    return out


def predictions(h: Dict[str, Any]) -> Dict[str, Any]:
    p2 = dict(
        flagged=(150 <= h["flagged"] <= 250),
        cost=(0.4 <= h["share_pct"] <= 1.0),
        fixed=(8 <= h["J5C_vs_J5"]["fixed"] <= 25),
        broken=(0 <= h["J5C_vs_J5"]["broken"] <= 2),
        median=(h["median_J5C"] - h["median_J5"] <= 0.0001),
    )
    return dict(
        p2=p2,
        p2_all=all(p2.values()),
        p3_j5=(-10 <= h["J5_vs_S8"]["net"] <= 5),
        p3_j5c=(-30 <= h["J5C_vs_S8"]["net"] <= -5),
    )


# ---- 3. plot, tables, pareto -------------------------------------------------------------------


def plot(rows: Dict[str, Any], extra: Dict[str, Any]) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axs = plt.subplots(2, 3, figsize=(13, 7))
    for ax, s in zip(axs.ravel(), ORDER):
        r = rows[s]
        pts = [
            ("S8", r["s8_total"], r["wrong_S8"], "tab:gray", "o"),
            ("J5", r["j5_total"], r["wrong_J5"], "tab:blue", "s"),
            ("J5C", r["j5_total"] + r["repair_evals"], r["wrong_J5C"], "tab:red", "D"),
        ]
        if s in extra.get("J5R", {}):
            e = extra["J5R"][s]
            pts.append(("J5R (step 19)", e[0], e[1], "tab:green", "^"))
        if s == "mid":
            for lab, (x, yv), col in (
                ("J5L", extra["J5L"], "tab:cyan"),
                ("J5L+C", extra["J5LC"], "tab:orange"),
            ):
                pts.append((lab, x, yv, col, "x"))
        for lab, x, yv, col, m in pts:
            ax.scatter([x / 1e6], [yv], c=col, marker=m, s=60, label=lab)
        ax.set_title(s + ("  (clean test)" if s == "held4" else "  (in-sample)"))
        ax.set_xlabel("evaluations with pass, millions")
        ax.set_ylabel("wrong voxels after pass")
        ax.grid(alpha=0.3)
    axs[0, 0].legend(fontsize=7)
    fig.suptitle(
        "AR5 step 20: evaluations against wrong voxels (more than 1 degree), 24,570 voxels"
    )
    fig.tight_layout()
    fig.savefig(OUT / "ar5_cmarepair.png", dpi=110)


def tables(rows: Dict[str, Any], jd: Dict[str, Any]) -> Dict[str, str]:
    L = [
        "| sample | S8 total | J5 search / pass / total | repair | J5C total | share of J5 "
        "| wrong S8 / J5 / J5C | median error over all 24,570 voxels, S8 / J5 / J5C (deg) |",
        "|---|---|---|---|---|---|---|---|",
    ]
    for s in ORDER:
        r = rows[s]
        L.append(
            f"| {s}{'' if s == 'held4' else ' (in-sample)'} | {r['s8_total']:,} | "
            f"{r['j5_search']:,} / {r['j5_pass']:,} / {r['j5_total']:,} | {r['repair_evals']:,} | "
            f"{jd[s]['j5c_total']:,} | {r['share_pct']:.2f}% | {r['wrong_S8']} / {r['wrong_J5']} / "
            f"{r['wrong_J5C']} | {r['median_S8']:.5f} / {r['median_J5']:.5f} / "
            f"{r['median_J5C']:.5f} |"
        )
    out = {"evals": "\n".join(L) + "\n"}
    L = [
        "| sample | flagged / fits | J5C vs J5 fixed / broken | J5C vs S8 fixed / broken (net) "
        "| J5 vs S8 fixed / broken (net) | r1 | r2 | r3 | r4 | pass |",
        "|---|---|---|---|---|---|---|---|---|---|",
    ]
    for s in ORDER:
        r, j = rows[s], jd[s]
        L.append(
            f"| {s} | {r['flagged']} / {r['fits']} | {r['J5C_vs_J5']['fixed']} / "
            f"{r['J5C_vs_J5']['broken']} | {r['J5C_vs_S8']['fixed']} / {r['J5C_vs_S8']['broken']} "
            f"({r['J5C_vs_S8']['net']:+d}) | {r['J5_vs_S8']['fixed']} / {r['J5_vs_S8']['broken']} "
            f"({r['J5_vs_S8']['net']:+d}) | {'yes' if j['r1'] else 'NO'} | "
            f"{'yes' if j['r2'] else 'NO'} | {'yes' if j['r3'] else 'NO'} | "
            f"{'yes' if j['r4'] else 'NO'} | {'yes' if j['pass'] else 'NO'} |"
        )
    out["pairs"] = "\n".join(L) + "\n"
    return out


def main() -> int:
    if not CG.held4_open():
        raise SystemExit("held4 marker missing: not reading held4")
    if "--tables-only" in sys.argv:  # phase B: rewrite tables.md (marked blocks) from summary.json
        d = json.loads((OUT / "summary.json").read_text())
        DT.write_tables(OUT / "tables.md", tables(d["rows"], d["judge"]))
        print((OUT / "tables.md").read_text())
        return 0
    val = validity()
    rows = {"held4": held4_rows(), **disclosed_rows()}
    (OUT / "stats_held4.json").write_text(json.dumps(rows["held4"], indent=1) + "\n")
    jd = {s: judge(rows[s]) for s in ORDER}
    pred = predictions(rows["held4"])
    mid_ok = rows["mid"]["J5C_vs_S8"]["net"] <= J["mid_net_vs_s8_max"]
    step_pass = bool(all(jd[s]["pass"] for s in ORDER) and mid_ok and val["ok"])
    # extra points for the plot: step 19's J5R, step 18's J5L and its repair
    s19 = json.loads((ICE / "benchmarks/autoresearch/ar5/repair/summary.json").read_text())
    extra: Dict[str, Any] = dict(
        J5R={
            s: (s19[s]["j5_total"] + s19[s]["repair_evals"], s19[s]["wrong_J5R"])
            for s in CR.SAMPLES
        }
    )
    gate = json.loads((CG.RD / "gate.json").read_text())["rows"]["mid_J5L"]
    extra["J5L"] = (CR.J5L["total"], CR.J5L["wrong"])
    extra["J5LC"] = (CR.J5L["total"] + gate["repair_evals"], gate["wrong_J5C"])
    plot(rows, extra)
    DT.write_tables(OUT / "tables.md", tables(rows, jd))
    pareto = []
    for s in ORDER:
        r = rows[s]
        for arm, tot, wr in (
            ("S8", r["s8_total"], r["wrong_S8"]),
            ("J5", r["j5_total"], r["wrong_J5"]),
            ("J5C", jd[s]["j5c_total"], r["wrong_J5C"]),
        ):
            lo, hi = ST.wilson(wr, N)
            pareto.append(
                dict(
                    sample=s,
                    in_sample=(s != "held4"),
                    arm=arm,
                    evals_total=int(tot),
                    evals_per_voxel=round(tot / N, 2),
                    wrong=int(wr),
                    fraction_right=round(1 - wr / N, 6),
                    wrong_wilson95=[round(lo, 6), round(hi, 6)],
                )
            )
    (OUT / "pareto_rows.json").write_text(json.dumps(pareto, indent=1) + "\n")
    out = dict(
        validity=val,
        rows=rows,
        judge=jd,
        mid_net_vs_s8_ok=bool(mid_ok),
        predictions=pred,
        step_pass=step_pass,
        j5l_mid=dict(
            repair_evals=gate["repair_evals"],
            wrong_J5C=gate["wrong_J5C"],
            fixed_vs_J5L=gate["fixed_vs_J5"],
            broken_vs_J5L=gate["broken_vs_J5"],
        ),
        note="in-sample rows: the flag and this form were chosen on them; held4 is the clean test",
    )
    (OUT / "summary.json").write_text(json.dumps(out, indent=1, default=float) + "\n")
    print((OUT / "tables.md").read_text())
    print("validity ok:", val["ok"], "| step pass:", step_pass, "| predictions:", json.dumps(pred))
    return 0


if __name__ == "__main__":
    sys.exit(main())
