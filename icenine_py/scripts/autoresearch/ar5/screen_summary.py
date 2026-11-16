"""AR5 step 21, phase A step 5: validity, statistics, predictions, plot, tables and pareto rows for
arm J5S on held5 (the clean test) and the six disclosed samples (in-sample).

  uv run python scripts/autoresearch/ar5/screen_summary.py

Run only after ``logs/screen_runs.done`` exists (``screen_guard`` refuses every held5 path before).
Writes ``benchmarks/autoresearch/ar5/screen/``: ``summary.json``, ``tables.md`` (doc_tables blocks),
``pareto_rows.json`` and ``ar5_screen.png``. No p-values: counts, paired fixed / broken by voxel,
Wilson 95% intervals. Medians are over all voxels of the sample (24,570 on held5).
"""

import hashlib
import json
import subprocess
import sys
from pathlib import Path
from typing import Any, Dict, List

import screen_guard as SG
import screen_rules as SR
from screen_freeze import FROZEN, ICE, problems

import numpy as np  # noqa: E402

sys.path.insert(0, str(ICE / "scripts" / "common"))
import doc_tables as DT  # noqa: E402
import stats as ST  # noqa: E402
import screen_data as SD  # noqa: E402

OUT = SG.RD
C = SG.cache_dir()
ORDER = ("held5",) + SR.SAMPLES
WRONG = SR.WRONG_DEG
J = SR.JUDGE


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


def jl(p: Path) -> Any:
    return json.loads(SG.guarded(p).read_text())


def ld(p: Path) -> np.ndarray:
    return np.load(SG.guarded(p)).astype(np.float64)


# ---- 1. validity -------------------------------------------------------------------------------


def validity() -> Dict[str, Any]:
    marker = C / SG.MARKER
    v: Dict[str, Any] = dict(marker=marker.read_text().splitlines()[0])
    marker_mtime = marker.stat().st_mtime
    n_vox = 24570
    for arm in ("S8", "J5"):
        rd = C / f"{arm}_held5"
        run = jl(rd / "run.json")
        v[f"{arm}_n_processed"] = int(run["n_processed"])
        v[f"{arm}_capped"] = (rd / "capped.json").exists()
        v[f"{arm}_pass_json"] = (C / "passes_held5" / f"D4_{arm}" / "pass.json").exists()
        v[f"{arm}_done"] = (C / "logs" / f"{arm}_held5.done").read_text().strip()
    for arm in ("J5C", "J5S"):
        rd = C / f"{arm}_held5"
        rep = jl(rd / "repair.json")
        v[f"{arm}_done"] = (C / "logs" / f"{arm}_held5.done").read_text().strip()
        v[f"{arm}_capped"] = (rd / "capped.json").exists()
        ev = np.load(SG.guarded(rd / "fits.npz"))["evals"]
        v[f"{arm}_all_fits_121"] = bool((ev == SR.EVALS_PER_FIT).all())
        v[f"{arm}_evals_equal_fit_sum"] = bool(
            int(ev.sum()) + rep["n_evals_score"] == rep["n_evals"]
            and rep["n_evals_fits"] == int(ev.sum())
        )
    sc = jl(C / "J5S_held5" / "repair.json")
    tr = np.load(SG.guarded(C / "J5S_held5" / "tries.npz"))
    v["J5S_tries_all_2_evals"] = bool((tr["evals"] == SR.EVALS_PER_SCREEN).all())
    v["J5S_pass_evals_equal_2_per_try"] = bool(
        sc["pass_evals"] == 2 * sc["n_tries"] == int(tr["evals"].sum())
    )
    v["J5S_extra_voxels_listed"] = bool(len(sc["extra_voxels"]) == sc["n_extra"])
    v["J5S_flagged_is_contested_plus_extra"] = bool(
        sc["n_flagged"] == sc["n_contested"] + sc["n_extra"]
    )
    v["J5C_J5S_same_J5_run"] = bool(
        jl(C / "J5C_held5" / "repair.json")["stored_pass_evals"] == sc["stored_pass_evals"]
    )
    # the arm reads no truth: the loader of the two arms is called with truth=False
    src = (Path(__file__).resolve().parent / "screen_run.py").read_text()
    v["arms_load_inputs_without_truth"] = bool(
        "load_inputs(s, truth=False)" in src and "read_truth" not in src
    )
    v["freeze_check"] = problems() or "ok"
    t = subprocess.run(
        ["git", "log", "-1", "--format=%ct", "--", str(FROZEN)],
        cwd=str(ICE),
        capture_output=True,
        text=True,
    ).stdout.strip()
    froz = float(t)
    files: List[Path] = []
    for p in (
        C / "S8_held5",
        C / "J5_held5",
        C / "passes_held5",
        C / "J5C_held5",
        C / "J5S_held5",
        ICE / "benchmarks" / "autoresearch" / "ar3" / "held5",
    ):
        files += [f for f in p.rglob("*") if f.is_file()]
    mt = [f.stat().st_mtime for f in files]
    v["freeze_commit_epoch"] = froz
    v["n_held5_output_files"] = len(files)
    v["oldest_held5_output_minus_freeze_s"] = float(min(mt) - froz)
    v["freeze_older_than_every_held5_output"] = bool(min(mt) > froz)
    # informational (this summary and any reader after the marker raise these): the access times
    # of the held5 maps were also taken by hand BEFORE any held5 read; see the report
    # access times of the held5 maps: none after the marker (macOS may not update atime on every
    # read: this is consistent with, not proof of, no read before the marker)
    maps = [f for f in files if f.suffix in (".npy", ".npz", ".mic", ".gz")]
    late = [f.name for f in maps if f.stat().st_atime > marker_mtime + 1]
    v["n_map_files_checked"] = len(maps)
    v["n_map_files_accessed_after_marker"] = len(late)
    v["after_S8_atime_minus_marker_s"] = float(
        (C / "passes_held5" / "after_S8.npy").stat().st_atime - marker_mtime
    )
    fz = json.loads(FROZEN.read_text())
    v["frozen_hash_count"] = len(fz["files"]) + len(fz["library"])
    v["frozen_sha_of_frozen_json"] = hashlib.sha256(FROZEN.read_bytes()).hexdigest()[:16]
    v["ok"] = bool(
        v["marker"] == "exit=0"
        and v["S8_n_processed"] == n_vox
        and v["J5_n_processed"] == n_vox
        and not any(v[f"{a}_capped"] for a in ("S8", "J5", "J5C", "J5S"))
        and v["S8_pass_json"]
        and v["J5_pass_json"]
        and v["S8_done"] == "run=0 pass=0"
        and v["J5_done"] == "run=0 pass=0"
        and v["J5C_done"] == "exit=0"
        and v["J5S_done"] == "exit=0"
        and v["J5C_all_fits_121"]
        and v["J5S_all_fits_121"]
        and v["J5C_evals_equal_fit_sum"]
        and v["J5S_evals_equal_fit_sum"]
        and v["J5S_tries_all_2_evals"]
        and v["J5S_pass_evals_equal_2_per_try"]
        and v["J5S_extra_voxels_listed"]
        and v["J5S_flagged_is_contested_plus_extra"]
        and v["J5C_J5S_same_J5_run"]
        and v["arms_load_inputs_without_truth"]
        and v["freeze_check"] == "ok"
        and v["freeze_older_than_every_held5_output"]
    )
    return v


# ---- 2. statistics -----------------------------------------------------------------------------


def side(err: np.ndarray) -> Dict[str, Any]:
    k = int((err > WRONG).sum())
    lo, hi = ST.wilson(k, len(err))
    return dict(wrong=k, wilson=[float(lo), float(hi)], median=float(np.median(err)))


def pair(w_old: np.ndarray, w_new: np.ndarray) -> Dict[str, int]:
    f, b = SR.fixed_broken(w_old, w_new)
    return dict(fixed=int(f), broken=int(b), net=int(b - f), both_wrong=int((w_old & w_new).sum()))


def row(s: str) -> Dict[str, Any]:
    y = SD.sym()
    tr = SD.read_truth(s)
    Rt = tr["R_t"]
    n = len(Rt)
    if s == "held5":
        p = C / "passes_held5"
        st: Dict[str, Any] = {}
        for a in ("S8", "J5"):
            st[a] = stage_evals(
                jl(C / f"{a}_held5" / "run.json"), int(jl(p / f"D4_{a}" / "pass.json")["evals"])
            )
        c_path = C / "J5C_held5"
    else:
        p = C / f"passes_{s}"
        if s == "held4":
            f = jl(SD.STATS["held4"])
            st = dict(
                S8=dict(before_pass=f["s8_search"], pass_=f["s8_pass"], total=f["s8_total"]),
                J5=dict(before_pass=f["j5_search"], pass_=f["j5_pass"], total=f["j5_total"]),
            )
        else:
            a = jl(SD.STATS16 / f"stats_{s}.json")["arms"]
            st = {k: a[k]["stages"] for k in ("S8", "J5")}
        c_path = C / f"stored_J5C_{s}"
    rc = jl(c_path / "repair.json")
    rs = jl(C / f"J5S_{s}" / "repair.json")
    maps = dict(
        S8=ld(p / "after_S8.npy"),
        J5=ld(p / "after_J5.npy"),
        J5C=ld(c_path / "after.npy"),
        J5S=ld(C / f"J5S_{s}" / "after.npy"),
    )
    err = {k: np.asarray(SR.misorientation_deg(R, Rt, y), float) for k, R in maps.items()}
    w = {k: e > WRONG for k, e in err.items()}
    j5t, s8t = int(st["J5"]["total"]), int(st["S8"]["total"])
    j5c_t = j5t + int(rc["n_evals"])
    j5s_t = SR.total_evals(j5t, rs["stored_pass_evals"], rs["pass_evals"], rs["n_evals"])
    sides = {k: side(e) for k, e in err.items()}
    pairs = {
        "J5S_vs_J5C": pair(w["J5C"], w["J5S"]),
        "J5S_vs_J5": pair(w["J5"], w["J5S"]),
        "J5S_vs_S8": pair(w["S8"], w["J5S"]),
        "J5C_vs_J5": pair(w["J5"], w["J5C"]),
        "J5C_vs_S8": pair(w["S8"], w["J5C"]),
        "J5_vs_S8": pair(w["S8"], w["J5"]),
    }
    judge = SR.judge(
        j5s_t,
        j5c_t,
        s8t,
        pairs["J5S_vs_J5C"]["net"],
        pairs["J5S_vs_S8"]["net"],
        sides["J5S"]["median"],
        sides["J5C"]["median"],
        sides["S8"]["median"],
    )
    return dict(
        sample=s,
        in_sample=s != "held5",
        n_voxels=n,
        evals=dict(
            S8=dict(search=st["S8"]["before_pass"], pass_=st["S8"]["pass_"], total=s8t),
            J5=dict(search=st["J5"]["before_pass"], pass_=st["J5"]["pass_"], total=j5t),
            J5C=dict(
                search=st["J5"]["before_pass"],
                pass_=st["J5"]["pass_"],
                repair=int(rc["n_evals"]),
                total=j5c_t,
            ),
            J5S=dict(
                search=st["J5"]["before_pass"],
                pass_=int(rs["pass_evals"]),
                repair=int(rs["n_evals"]),
                total=j5s_t,
            ),
        ),
        flagged=dict(
            J5C=int(rc["n_flagged"]),
            J5S=int(rs["n_flagged"]),
            J5S_contested=int(rs["n_contested"]),
            J5S_extra=int(rs["n_extra"]),
        ),
        fits=dict(J5C=int(rc["n_fits"]), J5S=int(rs["n_fits"])),
        changed=dict(J5C=int(rc["n_changed"]), J5S=int(rs["n_changed"])),
        saving_vs_J5C_pct=SR.saving_pct(j5c_t, j5s_t),
        saving_vs_J5_pct=SR.saving_pct(j5t, j5s_t),
        cut_vs_S8_pct=SR.saving_pct(s8t, j5s_t),
        repair_share_of_J5_pct=100.0 * rs["n_evals"] / j5t,
        pass_share_of_J5_pct=100.0 * rs["pass_evals"] / j5t,
        J5C_repair_share_of_J5_pct=100.0 * rc["n_evals"] / j5t,
        arms=sides,
        pairs=pairs,
        median_drop_vs_J5C_deg=sides["J5C"]["median"] - sides["J5S"]["median"],
        judge=judge,
        stop_line=dict(
            saving_under_8=bool(judge["saving_vs_j5c_pct"] < 8.0),
            net_over_2=bool(pairs["J5S_vs_J5C"]["net"] > 2),
            repair_over_10_pct_of_J5=bool(100.0 * rs["n_evals"] / j5t > 10.0),
        ),
    )


def predictions(r: Dict[str, Any]) -> Dict[str, Any]:
    p2 = dict(
        saving_11_18=bool(11.0 <= r["saving_vs_J5C_pct"] <= 18.0),
        flagged_1500_2200=bool(1500 <= r["flagged"]["J5S"] <= 2200),
        repair_5_8_pct=bool(5.0 <= r["repair_share_of_J5_pct"] <= 8.0),
        fixed_2_12=bool(2 <= r["pairs"]["J5S_vs_J5C"]["fixed"] <= 12),
        broken_0_4=bool(0 <= r["pairs"]["J5S_vs_J5C"]["broken"] <= 4),
        median_lower_0p001_0p004=bool(0.001 <= r["median_drop_vs_J5C_deg"] <= 0.004),
    )
    net_c = r["pairs"]["J5C_vs_S8"]["net"]
    net_s = r["pairs"]["J5S_vs_S8"]["net"]
    p3 = dict(
        j5c_net_in_m35_m5=bool(-35 <= net_c <= -5),
        j5s_within_3_of_j5c=bool(abs(net_s - net_c) <= 3),
        j5c_net=net_c,
        j5s_net=net_s,
    )
    return dict(
        p2=p2,
        p2_all=all(p2.values()),
        p3=p3,
        p3_all=bool(p3["j5c_net_in_m35_m5"] and p3["j5s_within_3_of_j5c"]),
    )


# ---- 3. the g0 note: why the replayed J5C map differs from the stored one ----------------------


def g0_note(s: str) -> Dict[str, Any]:
    """Step 20 read J5's tries as float32 and orthonormalized them; the gate's replay used the
    re-run's float64 ends. Compare the two states after the pass and the two final maps."""
    d = SD.load_inputs(s, truth=True)
    y = d["sym"]
    Rt = d["R_t"]
    st = d["stored"]
    full = np.load(SG.guarded(C / f"full_{s}" / "tries.npz"))
    args = (d["R_before"], d["old_score"])
    Ra, _c, _h = SR.after_pass_state(
        *args, st["voxel"], SD.orthonormal(st["end"]), st["cost"], st["hit"]
    )
    Rb_, _c2, _h2 = SR.after_pass_state(
        *args, full["voxel"], full["end"].astype(np.float64), full["cost"], full["hit"]
    )
    start_diff = np.abs(Ra - Rb_).reshape(len(Ra), -1).max(1) > 0
    a = ld(C / f"stored_J5C_{s}" / "after.npy")
    b = ld(C / f"J5C_{s}" / "after.npy")
    fin = np.abs(a - b).reshape(len(a), -1).max(1)
    diff = np.where(fin > 1e-6)[0]
    nb = d["nbrs"]
    # a fit at voxel v starts from the snapshot at v and at its candidates (neighbours): a start
    # that differs between the two sources can send the hard-cost CMA down another path
    touched = [bool(start_diff[v] or any(start_diff[j] for j in nb[v])) for v in diff]
    mis = np.asarray(SR.misorientation_deg(a[diff], b[diff], y), float) if len(diff) else []
    wa = np.asarray(SR.misorientation_deg(a, Rt, y), float) > WRONG
    wb = np.asarray(SR.misorientation_deg(b, Rt, y), float) > WRONG
    return dict(
        sample=s,
        n_start_differs_after_pass=int(start_diff.sum()),
        max_start_diff=float(np.abs(Ra - Rb_).max()),
        n_final_differ_over_1e6=int(len(diff)),
        n_of_those_with_start_diff_at_voxel_or_neighbour=int(sum(touched)),
        max_final_misorientation_deg=float(np.max(mis)) if len(diff) else 0.0,
        n_wrong_label_changes=int((wa != wb).sum()),
    )


# ---- 4. tables, plot, pareto ------------------------------------------------------------------


def tables(rows: Dict[str, Any]) -> Dict[str, str]:
    t: Dict[str, str] = {}
    ev = []
    for s in ORDER:
        r = rows[s]
        e = r["evals"]
        ev.append(
            dict(
                sample=s + ("" if s == "held5" else " (in-sample)"),
                S8=e["S8"]["total"],
                J5_search=e["J5"]["search"],
                J5_pass=e["J5"]["pass_"],
                J5=e["J5"]["total"],
                J5C_repair=e["J5C"]["repair"],
                J5C=e["J5C"]["total"],
                J5S_pass=e["J5S"]["pass_"],
                J5S_repair=e["J5S"]["repair"],
                J5S=e["J5S"]["total"],
                vs_J5C=r["saving_vs_J5C_pct"],
                vs_S8=r["cut_vs_S8_pct"],
            )
        )
    cols = list(ev[0])
    fm = {c: ",d" for c in cols if c not in ("sample", "vs_J5C", "vs_S8")}
    fm.update(vs_J5C=".1f", vs_S8=".1f")
    t["evals"] = DT.markdown_table(ev, cols, fm)
    fl = []
    for s in ORDER:
        r = rows[s]
        fl.append(
            dict(
                sample=s,
                J5C_flagged=r["flagged"]["J5C"],
                J5S_flagged=r["flagged"]["J5S"],
                contested=r["flagged"]["J5S_contested"],
                wide_extra=r["flagged"]["J5S_extra"],
                J5C_fits=r["fits"]["J5C"],
                J5S_fits=r["fits"]["J5S"],
                J5C_changed=r["changed"]["J5C"],
                J5S_changed=r["changed"]["J5S"],
                J5S_repair_pct_of_J5=r["repair_share_of_J5_pct"],
                J5S_pass_pct_of_J5=r["pass_share_of_J5_pct"],
            )
        )
    cols = list(fl[0])
    fm = {c: ",d" for c in cols if c != "sample" and "pct" not in c}
    fm.update(J5S_repair_pct_of_J5=".2f", J5S_pass_pct_of_J5=".2f")
    t["flags"] = DT.markdown_table(fl, cols, fm)
    wr = []
    for s in ORDER:
        r = rows[s]
        a = r["arms"]
        p = r["pairs"]
        wr.append(
            dict(
                sample=s,
                S8=a["S8"]["wrong"],
                J5=a["J5"]["wrong"],
                J5C=a["J5C"]["wrong"],
                J5S=a["J5S"]["wrong"],
                J5S_wilson="[%.3f%%, %.3f%%]"
                % (100 * a["J5S"]["wilson"][0], 100 * a["J5S"]["wilson"][1]),
                vs_J5C="%d / %d" % (p["J5S_vs_J5C"]["fixed"], p["J5S_vs_J5C"]["broken"]),
                vs_J5="%d / %d" % (p["J5S_vs_J5"]["fixed"], p["J5S_vs_J5"]["broken"]),
                vs_S8="%d / %d" % (p["J5S_vs_S8"]["fixed"], p["J5S_vs_S8"]["broken"]),
                J5C_vs_S8="%d / %d" % (p["J5C_vs_S8"]["fixed"], p["J5C_vs_S8"]["broken"]),
                net_J5C_J5S_vs_S8="%+d, %+d" % (p["J5C_vs_S8"]["net"], p["J5S_vs_S8"]["net"]),
            )
        )
    t["wrong"] = DT.markdown_table(wr, list(wr[0]))
    md = []
    for s in ORDER:
        r = rows[s]
        a = r["arms"]
        md.append(
            dict(
                sample=s,
                n=r["n_voxels"],
                S8=a["S8"]["median"],
                J5=a["J5"]["median"],
                J5C=a["J5C"]["median"],
                J5S=a["J5S"]["median"],
                drop_vs_J5C=r["median_drop_vs_J5C_deg"],
            )
        )
    fm = {c: ".5f" for c in ("S8", "J5", "J5C", "J5S", "drop_vs_J5C")}
    fm["n"] = ",d"
    t["median"] = DT.markdown_table(md, list(md[0]), fm)
    jr = []
    for s in ORDER:
        r = rows[s]
        j = r["judge"]
        jr.append(
            dict(
                sample=s,
                rule1_saving_vs_J5C=("pass " if j["ok1"] else "FAIL ")
                + "%.1f%%" % j["saving_vs_j5c_pct"],
                rule2_net_vs_J5C=("pass " if j["ok2"] else "FAIL ")
                + "%+d" % r["pairs"]["J5S_vs_J5C"]["net"],
                rule3_median=("pass " if j["ok3"] else "FAIL ")
                + "%+.5f" % -r["median_drop_vs_J5C_deg"],
                rule4_house=("pass " if j["ok4"] else "FAIL ")
                + "%.1f%%, %+d" % (j["cut_vs_s8_pct"], r["pairs"]["J5S_vs_S8"]["net"]),
                J5S_passes="yes" if j["ok"] else "no",
            )
        )
    t["rules"] = DT.markdown_table(jr, list(jr[0]))
    return t


def plot(rows: Dict[str, Any]) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axs = plt.subplots(2, 4, figsize=(15, 6.5))
    col = dict(S8="tab:gray", J5="tab:blue", J5C="tab:orange", J5S="tab:red")
    for ax, s in zip(axs.ravel(), ORDER):
        r = rows[s]
        for a in ("S8", "J5", "J5C", "J5S"):
            x = r["evals"][a]["total"] / 1e6
            yv = r["arms"][a]["wrong"]
            lo, hi = r["arms"][a]["wilson"]
            n = r["n_voxels"]
            ax.errorbar(
                [x],
                [yv],
                yerr=[[yv - lo * n], [hi * n - yv]],
                fmt="o",
                color=col[a],
                capsize=3,
                label=a,
            )
            ax.annotate(a, (x, yv), textcoords="offset points", xytext=(4, 4), fontsize=8)
        ax.set_title(s + (" (clean test)" if s == "held5" else " (in-sample)"), fontsize=10)
        ax.set_xlabel("evaluations in the run (millions)")
        ax.set_ylabel("wrong voxels (> 1 degree)")
        ax.grid(alpha=0.3)
    axs.ravel()[-1].axis("off")
    axs.ravel()[-1].text(
        0.0,
        0.5,
        "Left is cheaper, lower is better.\nBars: Wilson 95%.\nJ5S: screen tries + wide repair.",
        fontsize=9,
    )
    fig.suptitle("AR5 step 21: run evaluations against wrong voxels, after the repair")
    fig.tight_layout()
    fig.savefig(OUT / "ar5_screen.png", dpi=110)


def pareto(rows: Dict[str, Any]) -> List[Dict[str, Any]]:
    out = []
    for s in ORDER:
        r = rows[s]
        arms = ("S8", "J5", "J5C", "J5S") if s == "held5" else ("J5S",)
        for a in arms:
            e = r["evals"][a]["total"]
            k = r["arms"][a]["wrong"]
            lo, hi = r["arms"][a]["wilson"]
            n = r["n_voxels"]
            out.append(
                dict(
                    sample=s,
                    in_sample=s != "held5",
                    arm=a,
                    evals_total=int(e),
                    evals_per_voxel=round(e / n, 2),
                    wrong=int(k),
                    n=n,
                    fraction_right=round(1 - k / n, 6),
                    wrong_wilson95=[round(lo, 6), round(hi, 6)],
                )
            )
    return out


def main() -> int:
    if not (C / SG.MARKER).exists():
        raise SystemExit("the marker of ended held5 runs does not exist: nothing is read")
    val = validity()
    rows = {s: row(s) for s in ORDER}
    pred = predictions(rows["held5"])
    note = {s: g0_note(s) for s in SR.SAMPLES}
    verdict = dict(
        held5_pass=rows["held5"]["judge"]["ok"],
        disclosed_pass={s: rows[s]["judge"]["ok"] for s in SR.SAMPLES},
        step_passes=bool(all(rows[s]["judge"]["ok"] for s in ORDER)),
        stop_the_line=rows["held5"]["stop_line"],
    )
    summ = dict(validity=val, rows=rows, predictions=pred, g0_note=note, verdict=verdict)
    (OUT / "summary.json").write_text(json.dumps(summ, indent=1) + "\n")
    plot(rows)
    DT.write_tables(OUT / "tables.md", tables(rows))
    (OUT / "pareto_rows.json").write_text(json.dumps(pareto(rows), indent=1) + "\n")
    print(json.dumps(verdict, indent=1))
    print("validity ok:", val["ok"])
    return 0


if __name__ == "__main__":
    sys.exit(main())
