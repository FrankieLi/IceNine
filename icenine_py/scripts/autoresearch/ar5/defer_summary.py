"""AR5 step 22, phase A step 7: the numbers, the plot and the verdict of arm ``J5SD`` on the seven
disclosed samples and the clean test sample ``held6``, by the frozen decision rules.

  uv run python scripts/autoresearch/ar5/defer_summary.py

Reads (read only) the runs of ``defer_run_all.sh`` in ``benchmarks/autoresearch/ar5/defer/cache``
(``J5D_<s>``, ``j5sd/J5S_<s>``; on ``held6`` also ``S8_held6``, ``J5_held6``, ``J5S_held6`` and
their passes), the stored arms of the disclosed samples (links: ``J5_<s>``, ``stored_passes_<s>``,
``stored_J5S_<s>``), step 21's ``tables.md`` (the S8 and J5 totals), and the truth maps. ``held6``
is opened only after ``logs/defer_runs.done`` exists (``defer_guard``). Writes ``summary.json``,
``tables.md`` (blocks for ``sync_doc_tables.py``) and the figure ``ar5_defer_evals_wrong.png``.
No p-values: counts, Wilson intervals (``scripts/common/stats.py``) and paired fixed / broken by
voxel (wrong: over 1 degree from the truth, symmetry-reduced). The decision rules are
``defer_rules.JUDGE`` (frozen); nothing is tuned here.
"""

import gzip
import json
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import numpy as np

HERE = Path(__file__).resolve().parent
ICE = HERE.parents[2]
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(ICE / "scripts" / "common"))
import defer_guard as DG  # noqa: E402
import defer_rules as DR  # noqa: E402
import doc_tables as DT  # noqa: E402
import screen_rules as SR  # noqa: E402
import stats as ST  # noqa: E402

RD = ICE / "benchmarks" / "autoresearch" / "ar5" / "defer"
CACHE = RD / "cache"
S21_TABLES = ICE / "benchmarks" / "autoresearch" / "ar5" / "screen" / "tables.md"
AR3 = ICE / "benchmarks" / "autoresearch" / "ar3"
SAMPLES = DR.DISCLOSED + (DR.NEW_SAMPLE,)
ORDER = ("held5", "held4", "held3", "held", "mid", "grad", "big", "held6")


# ---- pure helpers (tested) -----------------------------------------------------------------------
def judge_vs_j5s(
    sample: str, fewer: float, broken_minus_fixed: int, median_higher: float
) -> Dict[str, Any]:
    """Rules (1) to (3) of the plan against ``J5S``."""
    j = DR.JUDGE
    # held6 is a held type: the plan's 7% for "the held types and mid"
    r1 = bool(fewer >= j["r1_fewer_min"].get(sample, j["r1_fewer_min"]["held"]))
    r2 = bool(broken_minus_fixed <= j["r2_broken_minus_fixed_max"])
    r3 = bool(median_higher <= j["r3_median_error_higher_max_deg"])
    return dict(r1=r1, r2=r2, r3=r3, ok=bool(r1 and r2 and r3))


def judge_vs_s8(fewer: float, broken_minus_fixed: int, median_higher: float) -> Dict[str, Any]:
    """Rule (4), the house rule against ``S8``."""
    h = DR.JUDGE["r4_vs_s8"]
    a = bool(fewer >= h["fewer"])
    b = bool(broken_minus_fixed <= h["broken_minus_fixed_max"])
    c = bool(median_higher <= h["median_error_higher_max_deg"])
    return dict(fewer=a, broken=b, median=c, ok=bool(a and b and c))


def step_verdict(rows: Dict[str, Dict[str, Any]]) -> Dict[str, Any]:
    """The step passes only if ``J5SD`` passes on ``held6`` and on all seven disclosed samples."""
    per = {s: bool(rows[s]["vs_j5s"]["ok"] and rows[s]["vs_s8"]["ok"]) for s in rows}
    need = set(DR.DISCLOSED) | {DR.NEW_SAMPLE}
    return dict(per_sample=per, step_passes=bool(need <= set(per) and all(per[s] for s in need)))


def in_range(x: float, lo: float, hi: float) -> bool:
    return bool(lo <= x <= hi)


# ---- readers -------------------------------------------------------------------------------------
def read_json(p: Path) -> Any:
    return json.loads(DG.guarded(p).read_text())


def read_gz(p: Path) -> Any:
    with gzip.open(DG.guarded(p), "rt") as f:
        return json.load(f)


def table_rows() -> Dict[str, Dict[str, int]]:
    """Step 21's table (first row of a sample): S8, J5_search, J5_pass, J5, J5S totals."""
    rows: Dict[str, Dict[str, int]] = {}
    for line in DG.guarded(S21_TABLES).read_text().splitlines():
        c = [x.strip() for x in line.strip("|").split("|")]
        if (
            len(c) == 12
            and c[0]
            and c[0].split()[0] in DR.DISCLOSED
            and c[0].split()[0] not in rows
        ):

            def n(s: str) -> int:
                return int(s.replace(",", ""))

            rows[c[0].split()[0]] = dict(
                s8=n(c[1]), j5_search=n(c[2]), j5_pass=n(c[3]), j5=n(c[4]), j5s=n(c[9])
            )
    return rows


def truth(s: str) -> np.ndarray:
    from icenine.mic_file import MicFile

    d = AR3 / "gradient" / "sample" if s == "grad" else AR3 / s / "gradient" / "sample"
    mic = MicFile.read(str(DG.guarded(d / "gradient_truth.mic")))
    return np.array([v.orientation for v in mic.voxels], dtype=np.float64)


def search_evals(st: Dict[str, Any]) -> int:
    return int(
        st["n_evals_seed"] + st["n_evals_reseed"] + st["n_evals_neighbor"] + st["n_evals_revisit"]
    )


def stages(run: Path, pol: Dict[str, Any], pass_ev: int, rep_ev: int) -> Dict[str, int]:
    st = read_json(run / "run.json")["stats"]
    return dict(
        index=int(st["index_evals"]),
        rescue=int(pol["rescue_evals"]),
        polish=int(pol["polish_evals"]),
        grid=int(pol["grid_evals"]),
        bfs_neighbour=int(st["n_evals_neighbor"]),
        bfs_revisit=int(st["n_evals_revisit"]),
        screen_pass=int(pass_ev),
        repair=int(rep_ev),
    )


def c_steps(s: str, run_j5: Path, run_d: Path) -> Dict[str, Any]:
    """Type counts of the replay of ``J5``'s stored seed steps, and the live ``J5D`` steps that
    match
    a type C step of the replay by voxel and reseed flag (only a step with a unique match): how many
    of them ended at another candidate. After the first C step the run order differs, so this is
    a partial matching, said so."""
    ent = read_gz(run_j5 / "e1_log.json.gz")
    trs = iter(read_gz(run_j5 / "mc_traces.json.gz"))
    walks: List[Dict[str, Any]] = []
    for x in ent:
        by: Dict[int, Dict[str, Any]] = {}
        for f in x["fits"]:
            if f[2] == "M":
                by[int(f[0])] = next(trs)
        walks.append(DR.walk_step(x, by))
    live = read_gz(run_d / "e1_log.json.gz")
    key: Dict[Tuple[int, bool], List[int]] = {}
    for i, x in enumerate(ent):
        key.setdefault((int(x["voxel"]), bool(x["reseed"])), []).append(i)
    n_c = moved = 0
    for e in live:
        ii = key.get((int(e["voxel"]), bool(e["reseed"])), [])
        if len(ii) == 1 and walks[ii[0]]["type"] == "C":
            n_c += 1
            moved += int(int(e["stop"]) != int(ent[ii[0]]["stop"]))
    types = {t: sum(1 for w in walks if w["type"] == t) for t in ("A", "C", "N", "-")}
    return dict(
        replay_types=types,
        j5_steps=len(ent),
        live_steps=len(live),
        c_matched=n_c,
        c_other_candidate=moved,
        c_other_share=(moved / n_c if n_c else None),
    )


def one(s: str, tab: Dict[str, Dict[str, int]], y: np.ndarray) -> Dict[str, Any]:
    new = s == DR.NEW_SAMPLE
    run_d = CACHE / f"J5D_{s}"
    arm = CACHE / "j5sd" / f"J5S_{s}"
    run_j5 = CACHE / f"J5_{s}"
    ps = CACHE / (f"passes_{s}" if new else f"stored_passes_{s}")
    R_t = truth(s)
    st_d = read_json(run_d / "run.json")["stats"]
    pol_d = read_json(run_d / "polish.json")
    pol_j = read_json(run_j5 / "polish.json")
    rep = read_json(arm / "repair.json")
    j5sd = search_evals(st_d) + int(rep["pass_evals"]) + int(rep["n_evals"])
    if new:
        st5 = read_json(run_j5 / "run.json")["stats"]
        j5 = search_evals(st5) + int(read_json(ps / "D4_J5" / "pass.json")["evals"])
        st8 = read_json(CACHE / "S8_held6" / "run.json")["stats"]
        s8 = search_evals(st8) + int(read_json(ps / "D4_S8" / "pass.json")["evals"])
        rep_s = read_json(CACHE / "J5S_held6" / "repair.json")
        j5s = search_evals(st5) + int(rep_s["pass_evals"]) + int(rep_s["n_evals"])
        m_j5s = np.load(DG.guarded(CACHE / "J5S_held6" / "after.npy"))
    else:
        t = tab[s]
        j5, s8, j5s = t["j5"], t["s8"], t["j5s"]
        m_j5s = np.load(DG.guarded(CACHE / f"stored_J5S_{s}" / "after.npy"))
    m_j5 = np.load(DG.guarded(ps / "after_J5.npy"))
    m_s8 = np.load(DG.guarded(ps / "after_S8.npy"))
    m_d = np.load(DG.guarded(arm / "after.npy"))
    w = {
        k: SR.wrong_mask(m, R_t, y) for k, m in dict(j5=m_j5, s8=m_s8, j5s=m_j5s, j5sd=m_d).items()
    }
    err = {
        k: np.asarray(SR.misorientation_deg(m, R_t, y), float)
        for k, m in dict(j5=m_j5, s8=m_s8, j5s=m_j5s, j5sd=m_d).items()
    }
    n = len(w["j5sd"])
    wrong = {}
    for k in w:
        lo, hi = ST.wilson(int(w[k].sum()), n)
        wrong[k] = dict(n=int(w[k].sum()), wilson=[float(lo), float(hi)])
    pair = {}
    for k in ("j5s", "j5", "s8"):
        f, b = SR.fixed_broken(w[k], w["j5sd"])
        pair[k] = dict(
            fixed=int(f),
            broken=int(b),
            net=int(b - f),
            median_higher=float(np.median(err["j5sd"]) - np.median(err[k])),
        )
    med = {k: float(np.median(v)) for k, v in err.items()}
    tot = dict(s8=int(s8), j5=int(j5), j5s=int(j5s), j5sd=int(j5sd))
    fewer = {k: 1.0 - j5sd / tot[k] for k in ("j5s", "j5", "s8")}
    df = read_json(run_d / "defer.json")
    pj = read_json(run_j5 / "polish.json")
    row: Dict[str, Any] = dict(
        sample=s,
        totals=tot,
        fewer=fewer,
        wrong=wrong,
        median_err_deg=med,
        pair=pair,
        stages=stages(run_d, pol_d, rep["pass_evals"], rep["n_evals"]),
        j5d_search=search_evals(st_d),
        seed_step_evals=int(st_d["n_evals_seed"] + st_d["n_evals_reseed"]),
        j5_seed_step_evals=None,
        defer={k: df[k] for k in DEFER_KEYS},
        rerun_evals_listed_apart=int(df["rerun_evals"]),
        grid=dict(j5=int(pj["grid_fallbacks_left"]), j5d=int(pol_d["grid_fallbacks_left"])),
        rescues=dict(
            j5_attempted=int(pj["rescues_attempted"]),
            j5_rescued=int(pj["rescued"]),
            j5d_attempted=int(pol_d["rescues_attempted"]),
            j5d_rescued=int(pol_d["rescued"]),
        ),
        c=c_steps(s, run_j5, run_d),
    )
    row["j5_seed_step_evals"] = int(
        read_json(run_j5 / "run.json")["stats"]["n_evals_seed"]
        + read_json(run_j5 / "run.json")["stats"]["n_evals_reseed"]
    )
    row["vs_j5s"] = judge_vs_j5s(s, fewer["j5s"], pair["j5s"]["net"], pair["j5s"]["median_higher"])
    row["vs_s8"] = judge_vs_s8(fewer["s8"], pair["s8"]["net"], pair["s8"]["median_higher"])
    return row


DEFER_KEYS = (
    "steps",
    "steps_with_pause",
    "pauses",
    "resumes",
    "found_tried_after_pause",
    "found_in_resume",
    "none_all_resumed",
    "rerun_checked",
)


def predictions(rows: Dict[str, Dict[str, Any]]) -> Dict[str, Any]:
    """Prediction 2 on the disclosed samples and prediction 3 on held6, each part judged."""
    p2: Dict[str, Any] = {}
    for s in DR.DISCLOSED:
        r, t = rows[s], DR.PLAN_TABLE[s]
        sav = 100 * r["fewer"]["j5s"]
        sh = r["c"]["c_other_share"]
        p2[s] = dict(
            saving_pct=sav,
            bounds_widened=[t["lo"] - 1.5, t["hi"] + 1.5],
            saving_inside=in_range(sav, t["lo"] - 1.5, t["hi"] + 1.5),
            c_other_share=sh,
            c_other_inside=None if sh is None else in_range(sh, 0.30, 0.70),
            net_inside=in_range(r["pair"]["j5s"]["net"], -3, 3),
            grid_within_2=abs(r["grid"]["j5d"] - r["grid"]["j5"]) <= 2,
        )
    h = rows[DR.NEW_SAMPLE]
    p3 = dict(
        saving_pct=100 * h["fewer"]["j5s"],
        saving_inside=in_range(100 * h["fewer"]["j5s"], 8.0, 14.0),
        net=h["pair"]["j5s"]["net"],
        net_inside=in_range(h["pair"]["j5s"]["net"], -3, 3),
        median_higher_deg=h["pair"]["j5s"]["median_higher"],
        median_inside=abs(h["pair"]["j5s"]["median_higher"]) <= 0.0005,
        house_rule=h["vs_s8"]["ok"],
    )
    return dict(p2=p2, p3=p3)


def figure(rows: Dict[str, Dict[str, Any]], path: Path) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axs = plt.subplots(2, 4, figsize=(15, 6.5))
    for ax, s in zip(axs.ravel(), ORDER):
        r = rows[s]
        for lab, k, ek, col, mk in (
            ("S8", "s8", "s8", "tab:gray", "s"),
            ("J5S", "j5s", "j5s", "tab:blue", "o"),
            ("J5SD", "j5sd", "j5sd", "tab:red", "D"),
        ):
            ax.scatter(r["totals"][ek] / 1e6, r["wrong"][k]["n"], c=col, marker=mk, s=60, label=lab)
            lo, hi = r["wrong"][k]["wilson"]
            n = 24570
            ax.vlines(r["totals"][ek] / 1e6, lo * n, hi * n, color=col, lw=1)
        ax.set_title(f"{s}{' (clean test)' if s == DR.NEW_SAMPLE else ''}")
        ax.set_xlabel("evaluations (millions)")
        ax.set_ylabel("wrong voxels")
        ax.set_ylim(bottom=0)
    axs[0, 0].legend()
    fig.suptitle("Run evaluations against wrong voxels (over 1 degree), Wilson 95% bars")
    fig.tight_layout()
    fig.savefig(path, dpi=110)
    plt.close(fig)


def tables(rows: Dict[str, Dict[str, Any]], pred: Dict[str, Any]) -> Dict[str, str]:
    def row_ev(s: str) -> Dict[str, Any]:
        r = rows[s]
        t = r["totals"]
        return dict(
            sample=s,
            s8=t["s8"],
            j5=t["j5"],
            j5s=t["j5s"],
            j5sd=t["j5sd"],
            vs_j5s=100 * r["fewer"]["j5s"],
            vs_j5=100 * r["fewer"]["j5"],
            vs_s8=100 * r["fewer"]["s8"],
        )

    def row_acc(s: str) -> Dict[str, Any]:
        r = rows[s]
        w = r["wrong"]
        lo, hi = w["j5sd"]["wilson"]
        return dict(
            sample=s,
            w_s8=w["s8"]["n"],
            w_j5=w["j5"]["n"],
            w_j5s=w["j5s"]["n"],
            w_j5sd=w["j5sd"]["n"],
            wilson=f"{100 * lo:.2f}-{100 * hi:.2f}%",
            fb_j5s=f"{r['pair']['j5s']['fixed']}/{r['pair']['j5s']['broken']}",
            fb_j5=f"{r['pair']['j5']['fixed']}/{r['pair']['j5']['broken']}",
            fb_s8=f"{r['pair']['s8']['fixed']}/{r['pair']['s8']['broken']}",
            med_j5s=r["median_err_deg"]["j5s"],
            med_j5sd=r["median_err_deg"]["j5sd"],
        )

    def row_dec(s: str) -> Dict[str, Any]:
        r = rows[s]
        a, b = r["vs_j5s"], r["vs_s8"]
        return dict(
            sample=s,
            fewer=100 * r["fewer"]["j5s"],
            r1=a["r1"],
            net=r["pair"]["j5s"]["net"],
            r2=a["r2"],
            med=r["pair"]["j5s"]["median_higher"],
            r3=a["r3"],
            house=b["ok"],
            passes=bool(a["ok"] and b["ok"]),
        )

    def row_def(s: str) -> Dict[str, Any]:
        r = rows[s]
        d, c = r["defer"], r["c"]
        return dict(
            sample=s,
            seed=r["seed_step_evals"],
            j5_seed=r["j5_seed_step_evals"],
            steps=d["steps"],
            with_pause=d["steps_with_pause"],
            pauses=d["pauses"],
            resumes=d["resumes"],
            tried=d["found_tried_after_pause"],
            resume=d["found_in_resume"],
            none=d["none_all_resumed"],
            rerun=r["rerun_evals_listed_apart"],
            c_matched=c["c_matched"],
            c_moved=c["c_other_candidate"],
            grid=f"{r['grid']['j5']}/{r['grid']['j5d']}",
            rescued=f"{r['rescues']['j5_rescued']}/{r['rescues']['j5d_rescued']}",
        )

    def row_stage(s: str) -> Dict[str, Any]:
        d = dict(rows[s]["stages"])
        d["sample"] = s
        return d

    o = ORDER
    return dict(
        ar5_defer_evals=DT.markdown_table(
            [row_ev(s) for s in o],
            ["sample", "s8", "j5", "j5s", "j5sd", "vs_j5s", "vs_j5", "vs_s8"],
            dict(s8=",d", j5=",d", j5s=",d", j5sd=",d", vs_j5s=".2f", vs_j5=".2f", vs_s8=".2f"),
        ),
        ar5_defer_accuracy=DT.markdown_table(
            [row_acc(s) for s in o],
            [
                "sample",
                "w_s8",
                "w_j5",
                "w_j5s",
                "w_j5sd",
                "wilson",
                "fb_j5s",
                "fb_j5",
                "fb_s8",
                "med_j5s",
                "med_j5sd",
            ],
            dict(med_j5s=".5f", med_j5sd=".5f"),
        ),
        ar5_defer_decision=DT.markdown_table(
            [row_dec(s) for s in o],
            ["sample", "fewer", "r1", "net", "r2", "med", "r3", "house", "passes"],
            dict(fewer=".2f", med=".5f"),
        ),
        ar5_defer_steps=DT.markdown_table(
            [row_def(s) for s in o],
            [
                "sample",
                "seed",
                "j5_seed",
                "steps",
                "with_pause",
                "pauses",
                "resumes",
                "tried",
                "resume",
                "none",
                "rerun",
                "c_matched",
                "c_moved",
                "grid",
                "rescued",
            ],
            dict(seed=",d", j5_seed=",d", rerun=",d"),
        ),
        ar5_defer_stages=DT.markdown_table(
            [row_stage(s) for s in o],
            [
                "sample",
                "index",
                "rescue",
                "polish",
                "grid",
                "bfs_neighbour",
                "bfs_revisit",
                "screen_pass",
                "repair",
            ],
            {
                k: ",d"
                for k in (
                    "index",
                    "rescue",
                    "polish",
                    "grid",
                    "bfs_neighbour",
                    "bfs_revisit",
                    "screen_pass",
                    "repair",
                )
            },
        ),
    )


def main() -> int:
    if not (CACHE / DG.MARKER).exists():
        print("held6's runs have not ended (no defer_runs.done): nothing is read")
        return 2
    marker = (CACHE / DG.MARKER).read_text().splitlines()[0]
    if marker != "exit=0":
        print("the runs did not end cleanly:", marker)
        return 2
    from icenine.ray_indexing import proper_rotations
    from icenine.symmetry import create_cubic_symmetry

    y = proper_rotations(create_cubic_symmetry(3.61))
    tab = table_rows()
    rows = {s: one(s, tab, y) for s in ORDER}
    verdict = step_verdict(rows)
    pred = predictions(rows)
    out = dict(rows=rows, verdict=verdict, predictions=pred, rule=DR.JUDGE)
    (RD / "summary.json").write_text(json.dumps(out, indent=1, default=float) + "\n")
    DT.write_tables(RD / "tables.md", tables(rows, pred))
    figure(rows, RD / "ar5_defer_evals_wrong.png")
    for s in ORDER:
        r = rows[s]
        print(
            f"{s:6s} J5SD {r['totals']['j5sd']:>11,} vs J5S {r['totals']['j5s']:>11,} "
            f"({100 * r['fewer']['j5s']:.2f}% fewer) net {r['pair']['j5s']['net']:+d} "
            f"med {r['pair']['j5s']['median_higher']:+.5f} | vs S8 {100 * r['fewer']['s8']:.1f}% "
            f"net {r['pair']['s8']['net']:+d} med {r['pair']['s8']['median_higher']:+.5f} | "
            f"vs_j5s {r['vs_j5s']['ok']} house {r['vs_s8']['ok']}"
        )
    print("step passes:", verdict["step_passes"])
    return 0


if __name__ == "__main__":
    sys.exit(main())
