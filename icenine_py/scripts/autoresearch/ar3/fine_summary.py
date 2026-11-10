"""AR3 step 4: summary of part 1 (fits from the stored painted starts with a finer box).

  uv run python scripts/autoresearch/ar3/fine_summary.py

Reads ``fine/runs/<arm>/refine.json.gz`` (arms B06, B03, B015, L, C), ``fine/runs/check/``,
step 3's ``gradient/runs/P`` (the stored paint, ``refine_R.json.gz`` = arm B5, ``run.json``) and
``gradient/reading0.json.gz`` (strict cost of each voxel's true orientation; analysis only).
Writes ``fine/summary.json``, ``fine/tables.md`` and ``fine/ar3_fine.png``. No p-values. Counts by
voxel with Wilson 95% intervals. A pick is made on the pick set A (grains with an even id) and
reported on the test set B (odd id). A difference is a gain or a loss only when intervals do not
overlap.
"""

import json
from pathlib import Path
from typing import Any, Dict, List, Tuple

import ar3_common as C
import gradient_common as G
import numpy as np
import run_gradient as RG

import doc_tables as DT  # noqa: E402
import stats as ST  # noqa: E402

FINE = C.OUT / "fine"
FRUNS = FINE / "runs"
LEVELS = G.LEVELS
ARMS_FULL = ["B5", "B06", "B03", "C"]  # arms run on all painted voxels
ARMS_SUB = ["B015", "L"]  # arms run on the first N painted voxels of the order
ARM_LABEL = {
    "start": "start (no fit)",
    "B5": "B5 stock box 5 deg",
    "B06": "B06 box 5/8 deg",
    "B03": "B03 box 5/16 deg",
    "B015": "B015 box 5/32 deg",
    "L": "L ladder 5/32, 5/8, 5",
    "C": "C CMA, step 0.05 deg",
}
ERR_EDGES = [0.0, 0.05, 0.1, 0.2, 0.5, 1e9]


def wil(k: int, n: int) -> Dict[str, Any]:
    lo, hi = ST.wilson(int(k), int(n)) if n else (float("nan"), float("nan"))
    return dict(k=int(k), n=int(n), lo=float(lo), hi=float(hi))


def wil_str(d: Dict[str, Any]) -> str:
    if not d["n"]:
        return "0 of 0"
    nd = 2 if d["hi"] < 0.01 else 1
    return f"{d['k']} of {d['n']} ({100 * d['lo']:.{nd}f} to {100 * d['hi']:.{nd}f}%)"


def overlap(a: Dict[str, Any], b: Dict[str, Any]) -> bool:
    return not (a["hi"] < b["lo"] or b["hi"] < a["lo"])


def lvl_name(L: float) -> str:
    return f"level {L:g}"


def bin_name(lo: float, hi: float) -> str:
    return f"{lo:g} to {hi:g}" if hi < 1e8 else f"{lo:g} and above"


class Arm:
    """One arm's fits: voxels, end orientations, costs, evaluations, start error and end error."""

    def __init__(self, name: str, rows: List[Dict[str, Any]], D: "Data") -> None:
        self.name = name
        self.rows = rows
        self.v = np.array([r["voxel"] for r in rows], dtype=int)
        self.end = np.array([r["end"] for r in rows], dtype=np.float64)
        self.end_cost = np.array([r["end_cost"] for r in rows])
        self.start_cost = np.array([r["start_cost"] for r in rows])
        self.evals = np.array([r["evals"] for r in rows], dtype=np.int64)
        self.e_start = D.eP[self.v]
        self.e_end = D.T.err(self.end, self.v)
        self.moved = D.T.mis(D.starts[self.v], self.end)
        self.lev = D.lev[self.v]
        self.isA = D.isA[self.v]
        self.cost_true = D.cost_true[self.v]

    def mask(self, which: str) -> np.ndarray:
        if which == "A":
            return self.isA
        if which == "B":
            return ~self.isA
        return np.ones(len(self.v), bool)


class Data:
    def __init__(self) -> None:
        self.T = G.Truth()
        T = self.T
        self.n = T.n
        self.lev = T.level
        self.isA = (T.gmap % 2) == 0
        P = RG.RUNS / "P"
        res = RG.load_paint(P)
        self.starts = np.zeros((self.n, 3, 3))
        for v in range(self.n):
            _R, src = res.final_orientation(v)
            if src == "grain" and res.stage.get(v) != "seed":
                g = res.painted_grain(v)
                assert g is not None
                self.starts[v] = np.asarray(res.entries[g].R, np.float64)
        self.eP = np.zeros(self.n)
        ok = np.abs(self.starts).sum(axis=(1, 2)) > 0
        self.eP[ok] = T.err(self.starts[ok], np.where(ok)[0])
        self.runP = json.loads((P / "run.json").read_text())
        r0 = np.array(G.load_gz(G.GRAD / "reading0.json.gz")["rows"])
        self.cost_true = np.full(self.n, np.nan)
        self.cost_true[r0[:, 0].astype(int)] = r0[:, 3]
        assert not np.isnan(self.cost_true).any()
        self.arms: Dict[str, Arm] = {}
        self.arms["B5"] = Arm("B5", G.load_gz(P / "refine_R.json.gz"), self)
        for a in ARMS_FULL[1:] + ARMS_SUB:
            f = FRUNS / a / "refine.json.gz"
            if f.exists():
                self.arms[a] = Arm(a, G.load_gz(f), self)
        self.order = self.arms["B5"].v  # the voxel order of the fits
        self.check_b5 = json.loads((FRUNS / "check" / "check_b5.json").read_text())
        self.runs = {
            a: json.loads((FRUNS / a / "run.json").read_text())
            for a in self.arms
            if a != "B5" and (FRUNS / a / "run.json").exists()
        }
        self.trace_check = {
            a: json.loads((FRUNS / a / "trace_check.json").read_text())
            for a in self.arms
            if (FRUNS / a / "trace_check.json").exists()
        }
        for a, arm in self.arms.items():
            if a in ARMS_FULL:
                assert len(arm.v) == len(self.order) and (arm.v == self.order).all()
            else:
                assert (arm.v == self.order[: len(arm.v)]).all()

    def sub(self, arm: Arm, other: Arm) -> np.ndarray:
        """Indices of ``other`` rows whose voxels are in ``arm`` (same order prefix)."""
        return np.arange(len(arm.v))


def stats_of(e_end: np.ndarray, e_start: np.ndarray, a: Arm, m: np.ndarray) -> Dict[str, Any]:
    n = int(m.sum())
    ee, es = e_end[m], e_start[m]
    ev = a.evals[m]
    lowered = int((a.end_cost[m] < a.start_cost[m]).sum())
    le_true = int((a.end_cost[m] <= a.cost_true[m]).sum())
    return dict(
        n=n,
        closer=wil(int((ee < es - G.HURT).sum()), n),
        hurt=wil(int((ee > es + G.HURT).sum()), n),
        median=float(np.median(ee)) if n else float("nan"),
        p95=float(np.percentile(ee, 95)) if n else float("nan"),
        gt01=wil(int((ee > 0.1).sum()), n),
        wrong=wil(int((ee > G.WRONG).sum()), n),
        lowered=wil(lowered, n),
        le_true=wil(le_true, n),
        evals_median=float(np.median(ev)) if n else float("nan"),
        evals_mean=float(ev.mean()) if n else float("nan"),
        evals_p95=float(np.percentile(ev, 95)) if n else float("nan"),
        evals_max=int(ev.max()) if n else 0,
        evals_total=int(ev.sum()),
    )


def start_stats(a: Arm, m: np.ndarray) -> Dict[str, Any]:
    es = a.e_start[m]
    n = int(m.sum())
    return dict(
        n=n,
        median=float(np.median(es)),
        p95=float(np.percentile(es, 95)),
        gt01=wil(int((es > 0.1).sum()), n),
        wrong=wil(int((es > G.WRONG).sum()), n),
        le_true=wil(int((a.start_cost[m] <= a.cost_true[m]).sum()), n),
    )


def restrict(D: Data, arm: Arm, n_first: int) -> Arm:
    """``arm`` cut to its first ``n_first`` fits (used to compare on a common voxel set)."""
    return Arm(arm.name, arm.rows[:n_first], D)


def main_table(D: Data, S: Dict[str, Any], T: Dict[str, str], which: str) -> None:
    """All arms on one voxel set. ``which``: 'full' (arms run on every painted voxel) or 'sub'
    (the first N voxels, every arm cut to them)."""
    arms = {a: D.arms[a] for a in D.arms}
    if which == "full":
        names = ARMS_FULL
        n_first = len(D.order)
    else:
        names = ["B5", "B06", "B03", "C", "B015", "L"]
        names = [a for a in names if a in arms]
        n_first = min(len(arms[a].v) for a in ARMS_SUB if a in arms)
    out: Dict[str, Any] = {}
    rows = []
    cut = {a: restrict(D, arms[a], n_first) for a in names}
    s0 = start_stats(cut[names[0]], np.ones(n_first, bool))
    out["start"] = s0
    rows.append(
        dict(
            arm=ARM_LABEL["start"],
            n=n_first,
            closer="",
            hurt="",
            med=s0["median"],
            p95=s0["p95"],
            gt01=wil_str(s0["gt01"]),
            wrong=wil_str(s0["wrong"]),
            lowered="",
            le_true=wil_str(s0["le_true"]),
            ev_med="",
            ev_tot="",
        )
    )
    for a in names:
        st = stats_of(cut[a].e_end, cut[a].e_start, cut[a], np.ones(n_first, bool))
        out[a] = st
        rows.append(
            dict(
                arm=ARM_LABEL[a],
                n=n_first,
                closer=wil_str(st["closer"]),
                hurt=wil_str(st["hurt"]),
                med=st["median"],
                p95=st["p95"],
                gt01=wil_str(st["gt01"]),
                wrong=wil_str(st["wrong"]),
                lowered=wil_str(st["lowered"]),
                le_true=wil_str(st["le_true"]),
                ev_med=st["evals_median"],
                ev_tot=st["evals_total"],
            )
        )
    S["main_" + which] = out
    T["ar3d_main_" + which] = DT.markdown_table(
        rows,
        [
            "arm",
            "n",
            "closer",
            "hurt",
            "med",
            "p95",
            "gt01",
            "wrong",
            "lowered",
            "le_true",
            "ev_med",
            "ev_tot",
        ],
        {"n": ",", "med": ".4f", "p95": ".4f", "ev_med": ",.0f", "ev_tot": ","},
    )


def evals_table(D: Data, S: Dict[str, Any], T: Dict[str, str]) -> None:
    """Evaluations per fit: how long the tail is."""
    rows = []
    S["evals_dist"] = {}
    for a, arm in D.arms.items():
        ev = arm.evals
        d = dict(
            n=len(ev),
            median=float(np.median(ev)),
            mean=float(ev.mean()),
            p95=float(np.percentile(ev, 95)),
            max=int(ev.max()),
            total=int(ev.sum()),
        )
        S["evals_dist"][a] = d
        rows.append(dict(arm=ARM_LABEL[a], n=d["n"], **{k: d[k] for k in d if k != "n"}))
    T["ar3d_evals"] = DT.markdown_table(
        rows,
        ["arm", "n", "median", "mean", "p95", "max", "total"],
        {"n": ",", "median": ",.0f", "mean": ",.0f", "p95": ",.0f", "max": ",", "total": ","},
    )


def level_tables(D: Data, S: Dict[str, Any], T: Dict[str, str]) -> None:
    """Share above 0.1 deg after the fit, by level, for the arms run on every voxel."""
    rows = []
    S["by_level"] = {}
    groups: List[Tuple[str, Any]] = [("all", None)] + [(lvl_name(L), L) for L in LEVELS]
    for name, L in groups:
        row: Dict[str, Any] = dict(level=name)
        d: Dict[str, Any] = {}
        m0 = np.ones(len(D.order), bool) if L is None else (D.arms["B5"].lev == L)
        row["n"] = int(m0.sum())
        s = wil(int((D.arms["B5"].e_start[m0] > 0.1).sum()), int(m0.sum()))
        d["start"] = s
        row["start"] = wil_str(s)
        for a in ARMS_FULL:
            arm = D.arms[a]
            w = wil(int((arm.e_end[m0] > 0.1).sum()), int(m0.sum()))
            d[a] = w
            row[a] = wil_str(w)
        S["by_level"][name] = d
        rows.append(row)
    T["ar3d_level"] = DT.markdown_table(rows, ["level", "n", "start"] + ARMS_FULL, {"n": ","})


def bin_tables(D: Data, S: Dict[str, Any], T: Dict[str, str]) -> None:
    """Closer / hurt by the size of the start's error, arms run on every voxel."""
    rows = []
    S["by_start_error"] = {}
    b5 = D.arms["B5"]
    for lo, hi in zip(ERR_EDGES[:-1], ERR_EDGES[1:]):
        m = (b5.e_start >= lo) & (b5.e_start < hi)
        n = int(m.sum())
        if not n:
            continue
        name = bin_name(lo, hi)
        row: Dict[str, Any] = dict(err=name, n=n)
        d: Dict[str, Any] = {}
        for a in ARMS_FULL:
            arm = D.arms[a]
            c = wil(int((arm.e_end[m] < arm.e_start[m] - G.HURT).sum()), n)
            h = int((arm.e_end[m] > arm.e_start[m] + G.HURT).sum())
            d[a] = dict(closer=c, hurt=h)
            row[a] = f"{c['k']} (hurt {h})"
        S["by_start_error"][name] = d
        rows.append(row)
    T["ar3d_bins"] = DT.markdown_table(rows, ["err", "n"] + ARMS_FULL, {"n": ","})


def set_tables(D: Data, S: Dict[str, Any], T: Dict[str, str]) -> Dict[str, Any]:
    """The pick on A (arms run on every voxel) and its report on B; the gate."""
    pick: Dict[str, Any] = {}
    A = {a: D.arms[a] for a in ARMS_FULL}
    stA = {a: stats_of(A[a].e_end, A[a].e_start, A[a], A[a].mask("A")) for a in A}
    best = min(ARMS_FULL, key=lambda a: stA[a]["gt01"]["k"] / stA[a]["gt01"]["n"])
    ties = [a for a in ARMS_FULL if overlap(stA[a]["gt01"], stA[best]["gt01"])]
    chosen = min(ties, key=lambda a: stA[a]["evals_total"])
    pick.update(
        lowest_share_arm=best,
        overlapping_arms=ties,
        chosen=chosen,
        rule=(
            "on A: the arm with the lowest share above 0.1 deg; among arms whose Wilson intervals "
            "overlap with it, the arm with the fewest total evaluations"
        ),
    )
    rows = []
    stB: Dict[str, Any] = {}
    sB = start_stats(A["B5"], A["B5"].mask("B"))
    for a in ARMS_FULL:
        stB[a] = stats_of(A[a].e_end, A[a].e_start, A[a], A[a].mask("B"))
    S["set_A"] = stA
    S["set_B"] = dict(start=sB, **stB)
    sA0 = start_stats(A["B5"], A["B5"].mask("A"))
    S["set_A"]["start"] = sA0
    for name, st0, st in (("A", sA0, stA), ("B", sB, stB)):
        rows.append(
            dict(
                set=name,
                arm="start (no fit)",
                n=st0["n"],
                med=st0["median"],
                gt01=wil_str(st0["gt01"]),
                wrong=wil_str(st0["wrong"]),
                ev_tot="",
                mark="",
            )
        )
        for a in ARMS_FULL:
            s = st[a]
            mark = ""
            if name == "A":
                mark = "picked" if a == chosen else ("in overlap" if a in ties else "")
            else:
                gain_s = s["gt01"]["hi"] < st0["gt01"]["lo"]
                gain_b5 = s["gt01"]["hi"] < stB["B5"]["gt01"]["lo"]
                mark = ("gain over start" if gain_s else "no gain over start") + (
                    "; gain over B5" if (gain_b5 and a != "B5") else ""
                )
            rows.append(
                dict(
                    set=name,
                    arm=ARM_LABEL[a],
                    n=s["n"],
                    med=s["median"],
                    gt01=wil_str(s["gt01"]),
                    wrong=wil_str(s["wrong"]),
                    ev_tot=s["evals_total"],
                    mark=mark,
                )
            )
    T["ar3d_sets"] = DT.markdown_table(
        rows,
        ["set", "arm", "n", "med", "gt01", "wrong", "ev_tot", "mark"],
        {"n": ",", "med": ".4f", "ev_tot": ","},
    )
    gains: Dict[str, Any] = {}
    for a in ARMS_FULL:
        gains[a] = dict(
            over_start=bool(stB[a]["gt01"]["hi"] < sB["gt01"]["lo"]),
            over_B5=bool(a != "B5" and stB[a]["gt01"]["hi"] < stB["B5"]["gt01"]["lo"]),
        )
    pick["gains_on_B"] = gains
    # the subset arms (B015, L): on their own first-N voxels, B half against the start there
    sub: Dict[str, Any] = {}
    for a in ARMS_SUB:
        if a not in D.arms:
            continue
        arm = D.arms[a]
        mB = arm.mask("B")
        st = stats_of(arm.e_end, arm.e_start, arm, mB)
        s0 = start_stats(arm, mB)
        sub[a] = dict(
            n=int(mB.sum()),
            end_gt01=st["gt01"],
            start_gt01=s0["gt01"],
            gain_over_start=bool(st["gt01"]["hi"] < s0["gt01"]["lo"]),
        )
    pick["subset_arms_on_B"] = sub
    pick["gate_open"] = bool(
        any(g["over_start"] for g in gains.values())
        or any(s["gain_over_start"] for s in sub.values())
    )
    S["pick"] = pick
    return pick


def ladder_table(D: Data, S: Dict[str, Any], T: Dict[str, str]) -> None:
    if "L" not in D.arms:
        return
    arm = D.arms["L"]
    n = len(arm.rows)
    from collections import Counter

    last = Counter(r["last_rung"] for r in arm.rows)
    fell = [r["end_cost"] < r["start_cost"] for r in arm.rows]
    rows = []
    S["ladder"] = {}
    for box in (5.0 / 32, 5.0 / 8, 5.0):
        k = int(last.get(box, 0))
        ev = [r["evals"] for r in arm.rows if r["last_rung"] == box]
        S["ladder"][f"{box:.5f}"] = dict(
            ended_here=k,
            evals_median=float(np.median(ev)) if ev else float("nan"),
        )
        rows.append(
            dict(
                rung=f"box {box:g} deg",
                ended=wil_str(wil(k, n)),
                evals_med=float(np.median(ev)) if ev else float("nan"),
            )
        )
    S["ladder"]["cost_fell"] = wil(int(sum(fell)), n)
    S["ladder"]["n"] = n
    T["ar3d_ladder"] = DT.markdown_table(
        rows, ["rung", "ended", "evals_med"], {"evals_med": ",.0f"}
    )


def figure(D: Data, S: Dict[str, Any]) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(1, 3, figsize=(16, 4.6))
    cols = {"B5": "k", "B06": "tab:blue", "B03": "tab:cyan", "C": "tab:red"}
    for a in ARMS_FULL:
        arm = D.arms[a]
        e = np.sort(arm.e_end)
        ax[0].plot(e, np.arange(1, len(e) + 1) / len(e), label=ARM_LABEL[a], color=cols[a])
    e = np.sort(D.arms["B5"].e_start)
    ax[0].plot(e, np.arange(1, len(e) + 1) / len(e), "--", color="gray", label=ARM_LABEL["start"])
    ax[0].set_xscale("log")
    ax[0].set_xlim(1e-3, 3)
    ax[0].axvline(0.1, color="gray", lw=0.5, ls=":")
    ax[0].set_xlabel("error to the truth (deg)")
    ax[0].set_ylabel("share of painted voxels at or below")
    ax[0].set_title("Error after the fit, all painted voxels")
    ax[0].legend(fontsize=7, loc="lower right")
    # share above 0.1 against evaluations per fit
    mf = S["main_full"]
    xs, ys, lo, hi, nm = [], [], [], [], []
    for a in ARMS_FULL:
        g = mf[a]["gt01"]
        xs.append(mf[a]["evals_mean"])
        ys.append(g["k"] / g["n"])
        lo.append(g["k"] / g["n"] - g["lo"])
        hi.append(g["hi"] - g["k"] / g["n"])
        nm.append(a)
    ax[1].errorbar(xs, ys, yerr=[lo, hi], fmt="o", color="k", capsize=3)
    for x, y, n_ in zip(xs, ys, nm):
        ax[1].annotate(n_, (x, y), textcoords="offset points", xytext=(5, 5), fontsize=8)
    s0 = mf["start"]["gt01"]
    ax[1].axhline(s0["k"] / s0["n"], color="gray", ls="--", lw=0.8)
    ax[1].set_xscale("log")
    ax[1].set_xlabel("mean evaluations per fit")
    ax[1].set_ylabel("share above 0.1 deg after the fit")
    ax[1].set_title("Cost against error (dashed: the start)")
    # error by level
    w = 0.18
    for i, a in enumerate(["start"] + ARMS_FULL):
        med = []
        for L in LEVELS:
            m = D.arms["B5"].lev == L
            arm = D.arms["B5" if a == "start" else a]
            e = arm.e_start[m] if a == "start" else arm.e_end[m]
            med.append(float(np.percentile(e, 95)))
        ax[2].bar(
            np.arange(len(LEVELS)) + (i - 2) * w,
            med,
            w,
            label=ARM_LABEL[a],
            color=("gray" if a == "start" else cols[a]),
        )
    ax[2].set_xticks(range(len(LEVELS)))
    ax[2].set_xticklabels([f"{L:g}" for L in LEVELS])
    ax[2].set_xlabel("grain level (deg)")
    ax[2].set_ylabel("95th percentile of error (deg)")
    ax[2].set_title("95th percentile by level")
    ax[2].legend(fontsize=6)
    fig.tight_layout()
    fig.savefig(FINE / "ar3_fine.png", dpi=110)


def main() -> None:
    D = Data()
    S: Dict[str, Any] = dict(
        n_voxels=D.n,
        n_painted=len(D.order),
        truth_sha256=D.T.sha,
        check_b5=D.check_b5,
        trace_check=D.trace_check,
        runs=D.runs,
        evals_R_stored=D.runP["evals_R"],
        evals_B5_summed=int(D.arms["B5"].evals.sum()),
    )
    assert S["evals_R_stored"] == S["evals_B5_summed"]
    T: Dict[str, str] = {}
    main_table(D, S, T, "full")
    if all(a in D.arms for a in ARMS_SUB):
        main_table(D, S, T, "sub")
    evals_table(D, S, T)
    level_tables(D, S, T)
    bin_tables(D, S, T)
    set_tables(D, S, T)
    ladder_table(D, S, T)
    if (FRUNS / "C_trace" / "refine2.json.gz").exists():
        import fine_part2 as P2

        assert S["pick"]["chosen"] == "C", "part 2 is written for the picked arm C"
        S["trace_check_C_trace"] = json.loads((FRUNS / "C_trace" / "trace_check.json").read_text())
        P2.part2(D, S, T, "C_trace")
    figure(D, S)
    (FINE / "summary.json").write_text(json.dumps(S, indent=1, default=float) + "\n")
    DT.write_tables(FINE / "tables.md", T)
    print("wrote summary.json, tables.md, ar3_fine.png")


if __name__ == "__main__":
    main()
