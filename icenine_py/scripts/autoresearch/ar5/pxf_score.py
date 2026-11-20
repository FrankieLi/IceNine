"""AR5 step 28: scoring of finished runs against the frozen gates (reads maps and the truth for
scoring only; nothing is tuned here). Not frozen: it applies the frozen rules to finished runs.

  uv run python scripts/autoresearch/ar5/pxf_score.py live|h1

``live``: for ``mid`` and ``held``, the new J5SD (``live_pxf``) against stored J5SD (step 22;
the old arm is step 27's old rerun, equal to it): evaluations (search + screen pass + wide repair),
wrong voxels (over 1 degree from the truth, symmetry reduced; Wilson 95%), fixed / broken by voxel,
median error, grid searches, wall seconds with the preflight label and the load, indexer work from
the wrappers' logs (fallback indexings included), the fallback counts, and the replay check (a
fallback step at a voxel seeded once in the old run equals the stored old step: evaluations, fits,
hit ratio). ``h1``: the new J5SD on ``held7`` against ``S8``, plus the plain J5SD (reported, not a
gate); its logs and maps open only after ``logs/bf_runs.done`` (reads go through
``bf_guard.guarded``).
"""

import gzip
import json
import os
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional

import numpy as np

HERE = Path(__file__).resolve().parent
ICE = HERE.parents[2]
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(ICE / "scripts" / "common"))


def _default_cache() -> Path:
    parts = ICE.parts
    if ".claude" in parts and "worktrees" in parts:
        main = Path(*parts[: parts.index(".claude")])
        return main / "icenine_py" / "benchmarks" / "autoresearch" / "ar5" / "pxf" / "cache"
    return ICE / "benchmarks" / "autoresearch" / "ar5" / "pxf" / "cache"


os.environ.setdefault("PXF_CACHE", str(_default_cache()))
for _k in ("SCREEN_CACHE", "DEFER_CACHE", "BALL_CACHE", "BF_CACHE"):
    os.environ.setdefault(_k, os.environ["PXF_CACHE"])

import bf_guard as BG  # noqa: E402
import bf_score as BS  # noqa: E402
import pxf_rules as PR  # noqa: E402
import screen_rules as SR  # noqa: E402
import stats as ST  # noqa: E402

CACHE = Path(os.environ["PXF_CACHE"])
RD = ICE / "benchmarks" / "autoresearch" / "ar5" / "pxf"
PXD = ICE / "benchmarks" / "autoresearch" / "ar5" / "px"  # step 27's tracked results
SCORE27 = PXD / "score_live.json"
T1_27 = PXD / "score_t1.json"


def read_gz(p: Path) -> Any:
    with gzip.open(BG.guarded(p), "rt") as f:
        return json.load(f)


def wilson(k: int, n: int) -> List[float]:
    lo, hi = ST.wilson(int(k), int(n))
    return [float(lo), float(hi)]


def index_work(c: Path, s: str) -> Dict[str, Any]:
    """Summed wrapper logs of the run, the boundary pass and the screen/repair arm."""
    out: Dict[str, Any] = {}
    tot = dict(calls=0, pairs=0, votes=0, seconds=0.0)
    for part in ("run", "pass", "arm"):
        p = c / f"pxf_index_{part}_{s}.json.gz"
        if not BG.guarded(p).exists():
            out[part] = None
            continue
        d = read_gz(p)
        out[part] = d["summed"]
        for k in tot:
            tot[k] += d["summed"][k]
    out["total"] = tot
    return out


def job_times(job: str) -> Dict[str, Any]:
    lg = CACHE / "logs"
    rd = lambda n: (lg / n).read_text().strip() if (lg / n).exists() else None  # noqa: E731
    t0, t1 = rd(f"{job}.t0"), rd(f"{job}.t1")
    load: List[float] = []
    lf = lg / f"{job}.load"
    if lf.exists():
        for ln in lf.read_text().splitlines():
            parts = ln.replace("{", " ").replace("}", " ").split()
            if len(parts) >= 3:
                load.append(float(parts[1]))
    return dict(
        label=rd(f"{job}.label"),
        wall_s=(int(t1) - int(t0)) if t0 and t1 else None,
        load_samples=len(load),
        load_mean=float(np.mean(load)) if load else None,
        load_max=float(np.max(load)) if load else None,
    )


def maps_scores(maps: Dict[str, np.ndarray], R_t: np.ndarray) -> Dict[str, Any]:
    y = BS.sym()
    w = {k: SR.wrong_mask(m, R_t, y) for k, m in maps.items()}
    err = {k: np.asarray(SR.misorientation_deg(m, R_t, y), float) for k, m in maps.items()}
    n = len(next(iter(w.values())))
    return dict(
        n=n,
        w=w,
        err=err,
        wrong={k: dict(n=int(w[k].sum()), wilson=wilson(int(w[k].sum()), n)) for k in w},
        median_err_deg={k: float(np.median(v)) for k, v in err.items()},
    )


def pair(sc: Dict[str, Any], base: str, new: str) -> Dict[str, Any]:
    f, b = SR.fixed_broken(sc["w"][base], sc["w"][new])
    return dict(
        fixed=int(f),
        broken=int(b),
        net=int(b - f),
        median_higher=float(np.median(sc["err"][new]) - np.median(sc["err"][base])),
    )


def step_counts(c: Path, s: str) -> Dict[str, Any]:
    """The seed steps of the run from ``pxf_steps_<s>``: steps, failures, fallbacks, outcomes."""
    rows = read_gz(c / f"pxf_steps_{s}.json.gz")["rows"]
    fb = [r for r in rows if r["fallback"]]
    fail_after = [r for r in rows if (r["fallback"] and not r["passed_second"])]
    fail_nofb = [r for r in rows if (not r["fallback"] and r["failed_first"])]
    return dict(
        seed_steps=len(rows),
        failed_first=sum(1 for r in rows if r["failed_first"]),
        failed_first_subsampled=sum(1 for r in rows if r["failed_first"] and r["subsampled"]),
        fallbacks=len(fb),
        fallbacks_passed=sum(1 for r in fb if r["passed_second"]),
        fallbacks_failed=len(fail_after),
        failed_not_subsampled=len(fail_nofb),
        to_rescue=len(fail_after) + len(fail_nofb),
        first_step_evals=int(sum(r["first"]["evals"] for r in rows)),
        fallback_step_evals=int(sum(r["second"]["evals"] for r in fb)),
        fallback_pairs=int(sum(r["second"]["n_pairs"] for r in fb)),
        fallback_votes=int(sum(r["second"]["n_votes"] for r in fb)),
        fallback_index_s=float(sum(r["second"]["index_s"] for r in fb)),
        first_index_s=float(sum(r["first"]["index_s"] for r in rows)),
    )


def replay_check(c: Path, s: str, old_run: Optional[Path] = None) -> Dict[str, Any]:
    """A fallback step at a voxel seeded once in the old run (stored J5D) equals the stored old
    step: evaluations, candidate count, fits (orientation, hit, cost, evaluations), hit ratio."""
    old = read_gz((old_run or CACHE / f"J5D_{s}") / "e1_log.json.gz")
    by_old: Dict[int, List[Dict[str, Any]]] = {}
    for x in old:
        by_old.setdefault(int(x["voxel"]), []).append(x)
    new = read_gz(c / f"J5D_{s}" / "e1_log.json.gz")
    fb_rows: Dict[int, List[Dict[str, Any]]] = {}
    for x in new:
        if x.get("pxf_stage") == "fallback":
            fb_rows.setdefault(int(x["voxel"]), []).append(x)
    checked = same = 0
    diffs: List[Dict[str, Any]] = []
    for v, xs in sorted(fb_rows.items()):
        if len(by_old.get(v, [])) != 1 or len(xs) != 1:
            continue
        o, n = by_old[v][0], xs[0]
        checked += 1
        keys = ("evals", "n_cand", "fits", "fin_hit", "fin_evals")
        bad = [k for k in keys if o.get(k) != n.get(k)]
        if bad:
            diffs.append(
                dict(voxel=v, differs=bad, old_evals=o.get("evals"), new_evals=n.get("evals"))
            )
        else:
            same += 1
    return dict(
        fallback_voxels=len(fb_rows),
        seeded_once_in_old=checked,
        equal=same,
        differ=len(diffs),
        diffs=diffs[:20],
        ok=bool(not diffs),
    )


def stored_old(s: str) -> Optional[Dict[str, Any]]:
    """Step 27's old run of ``s`` (tracked ``score_live.json``): totals, wall, indexer work."""
    if not SCORE27.exists():
        return None
    r = json.loads(SCORE27.read_text())[s]
    return dict(
        total=r["totals"]["old"],
        job=r["job"]["old"],
        index_work=r["index_work"]["old"],
        grid_fallbacks_left=r["grid_fallbacks_left"]["old"],
    )


def live(s: str) -> Dict[str, Any]:
    R_t = BS.truth(s)
    stored_run, stored_arm = CACHE / f"J5D_{s}", CACHE / f"stored_J5SD_{s}"
    c = CACHE / "live_pxf"
    new_run, new_arm = c / f"J5D_{s}", c / "j5sd" / f"J5S_{s}"
    tot = dict(stored=BS.run_total(stored_run, stored_arm), new=BS.run_total(new_run, new_arm))
    maps = dict(
        stored=np.load(BG.guarded(stored_arm / "after.npy")),
        new=np.load(BG.guarded(new_arm / "after.npy")),
    )
    pol = dict(
        stored=BS.read_json(stored_run / "polish.json"), new=BS.read_json(new_run / "polish.json")
    )
    st = dict(
        stored=BS.read_json(stored_run / "run.json")["stats"],
        new=BS.read_json(new_run / "run.json")["stats"],
    )
    sc = maps_scores(maps, R_t)
    old = stored_old(s)
    row: Dict[str, Any] = dict(
        sample=s,
        totals=tot,
        stored_total_matches_step22=bool(tot["stored"]["total"] == PR.STORED_TOTALS[s]),
        wrong=sc["wrong"],
        median_err_deg=sc["median_err_deg"],
        new_vs_stored=pair(sc, "stored", "new"),
        grid_fallbacks_left={k: int(v["grid_fallbacks_left"]) for k, v in pol.items()},
        index_fallback_seed_test={k: int(v["index_fallback_seed_test"]) for k, v in st.items()},
        n_seeds={k: int(v["index_runs"]) for k, v in st.items()},
        run_index_s={k: float(v["index_s"]) for k, v in st.items()},
        run_wall_s=float(BS.read_json(new_run / "run.json")["wall_total_s"]),
        job=dict(new=job_times(f"live_pxf_{s}"), old=old["job"] if old else None),
        index_work=dict(new=index_work(c, s), old=old["index_work"] if old else None),
        old_run=old,
        steps=step_counts(c, s),
        replay=replay_check(c, s),
    )
    pa = row["index_work"]["new"]
    row["pass_and_arm_index_calls"] = int(
        sum((pa[k] or dict(calls=0))["calls"] for k in ("pass", "arm"))
    )
    row["c3"] = PR.judge_c3(
        row["new_vs_stored"]["net"],
        row["new_vs_stored"]["median_higher"],
        row["grid_fallbacks_left"]["new"] - row["grid_fallbacks_left"]["stored"],
    )
    ow, nw = (old["job"]["wall_s"] if old else None), row["job"]["new"]["wall_s"]
    row["wall_fewer"] = PR.fewer(nw, ow) if ow and nw else None
    return row


def live_all() -> Dict[str, Any]:
    rows = {s: live(s) for s in PR.FULL_SAMPLES}
    t1 = json.loads(T1_27.read_text())["pooled"]
    npairs = sum(r["index_work"]["new"]["total"]["pairs"] for r in rows.values())
    nvotes = sum(r["index_work"]["new"]["total"]["votes"] for r in rows.values())
    tot = sum(r["totals"]["new"]["total"] for r in rows.values())
    nw = sum(r["job"]["new"]["wall_s"] for r in rows.values())
    judge = dict(
        c1=PR.judge_c1(npairs, nvotes, int(t1["old_pairs"]), int(t1["old_votes"])),
        c2=PR.judge_c2(tot),
        c3={s: r["c3"] for s, r in rows.items()},
        t2=dict(PR.judge_t2(nw), new_wall_s=int(nw), old_wall_s=PR.OLD_WALL_S),
        replay_ok=all(r["replay"]["ok"] for r in rows.values()),
        pass_arm_calls_zero=all(r["pass_and_arm_index_calls"] == 0 for r in rows.values()),
    )
    judge["c1"]["new_pairs"], judge["c1"]["new_votes"] = int(npairs), int(nvotes)
    judge["c1"]["old_pairs"], judge["c1"]["old_votes"] = int(t1["old_pairs"]), int(t1["old_votes"])
    judge["all_ok"] = bool(
        judge["c1"]["ok"]
        and judge["c2"]["ok"]
        and all(v["ok"] for v in judge["c3"].values())
        and judge["replay_ok"]
    )
    return dict(rows=rows, judge=judge)


def h1() -> Dict[str, Any]:
    s = "held7"
    c, co = CACHE / "live_pxf", CACHE / "live_old"
    R_t = BS.truth(s)
    arms = dict(
        pxf=(c / f"J5D_{s}", c / "j5sd" / f"J5S_{s}"),
        old=(co / f"J5D_{s}", co / "j5sd" / f"J5S_{s}"),
    )
    tot = {k: BS.run_total(r, m) for k, (r, m) in arms.items()}
    s8 = BS.search_evals(BS.read_json(c / f"S8_{s}" / "run.json")["stats"]) + int(
        BS.read_json(c / "passes_held7" / "D4_S8" / "pass.json")["evals"]
    )
    maps = {k: np.load(BG.guarded(m / "after.npy")) for k, (r, m) in arms.items()}
    maps["s8"] = np.load(BG.guarded(c / "passes_held7" / "after_S8.npy"))
    sc = maps_scores(maps, R_t)
    pr = pair(sc, "s8", "pxf")
    fe = PR.fewer(tot["pxf"]["total"], s8)
    pol = {k: BS.read_json(r / "polish.json") for k, (r, m) in arms.items()}
    row = dict(
        sample=s,
        totals=dict(s8=int(s8), pxf=int(tot["pxf"]["total"]), old=int(tot["old"]["total"])),
        stages={k: BS.stage_split(r, tot[k]) for k, (r, m) in arms.items()},
        wrong=sc["wrong"],
        median_err_deg=sc["median_err_deg"],
        pxf_vs_s8=pr,
        old_vs_s8=pair(sc, "s8", "old"),
        pxf_vs_old=pair(sc, "old", "pxf"),
        fewer=fe,
        fewer_old=PR.fewer(tot["old"]["total"], s8),
        grid_fallbacks_left={k: int(v["grid_fallbacks_left"]) for k, v in pol.items()},
        job=dict(
            pxf=job_times("live_pxf_held7"),
            old=job_times("live_old_held7"),
            s8=job_times("S8_held7"),
        ),
        index_work=dict(pxf=index_work(c, s), old=index_work(co, s)),
        steps=step_counts(c, s),
        replay=replay_check(c, s, co / f"J5D_{s}"),
        h1=PR.judge_h1(fe, pr["net"], pr["median_higher"]),
    )
    return row


if __name__ == "__main__":
    what = sys.argv[1] if len(sys.argv) > 1 else ""
    if what == "live":
        out: Any = live_all()
    elif what == "h1":
        out = h1()
    else:
        raise SystemExit("usage: pxf_score.py live|h1")
    (RD / f"score_{what}.json").write_text(json.dumps(out, indent=1) + "\n")
    print(json.dumps(out, indent=1)[:6000])
