"""AR5 step 29: scoring of finished runs against the frozen gates (reads maps and the truth for
scoring only; nothing is tuned here). Not frozen: it applies the frozen rules to finished runs.

  uv run python scripts/autoresearch/ar5/ru_score.py live [--write-gate]
  uv run python scripts/autoresearch/ar5/ru_score.py held8

``live``: for ``mid`` and ``held``, the new J5SD (rules PXF and RU; ``<cache>/live_ru``) against
step 28's stored PXF J5SD runs (the main checkout's ``ar5/pxf/cache/live_pxf``; read only, not
rerun): evaluations (search + screen pass + wide repair), wrong voxels (over 1 degree from the
truth, symmetry reduced; Wilson 95%), fixed / broken by voxel, median error, grid searches, wall
seconds with the preflight label and the load, the try counts and the validity checks. r1 to r3
are judged; ``--write-gate`` writes ``<cache>/logs/gate_passed`` only when all three pass.
``held8``: S8, the PXF J5SD and the new J5SD (h1, h2); its logs and maps open only after
``logs/bf_runs.done`` (reads go through ``bf_guard.guarded``).
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

import ru_common as RC  # noqa: E402  (first: freeze, guards, cache)

import bf_guard as BG  # noqa: E402
import bf_score as BS  # noqa: E402
import ru_rules as RR  # noqa: E402
import screen_rules as SR  # noqa: E402
import stats as ST  # noqa: E402

CACHE = Path(os.environ["RU_CACHE"])
RD = ICE / "benchmarks" / "autoresearch" / "ar5" / "ru"


def _old_cache() -> Path:
    parts = ICE.parts
    if ".claude" in parts and "worktrees" in parts:
        main = Path(*parts[: parts.index(".claude")])
    else:
        main = ICE.parent
    return main / "icenine_py" / "benchmarks" / "autoresearch" / "ar5" / "pxf" / "cache"


OLD = Path(os.environ.get("RU_OLD_CACHE", str(_old_cache())))


def read_gz(p: Path) -> Any:
    with gzip.open(BG.guarded(p), "rt") as f:
        return json.load(f)


def wilson(k: int, n: int) -> List[float]:
    lo, hi = ST.wilson(int(k), int(n))
    return [float(lo), float(hi)]


def job_times(job: str, logs: Optional[Path] = None) -> Dict[str, Any]:
    lg = logs or (CACHE / "logs")
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
    both = int((sc["w"][base] & sc["w"][new]).sum())
    return dict(
        fixed=int(f),
        broken=int(b),
        both_wrong=both,
        net=int(b - f),
        median_higher=float(np.median(sc["err"][new]) - np.median(sc["err"][base])),
    )


def index_s(c: Path, s: str, part: str = "run") -> Dict[str, Any]:
    p = c / f"pxf_index_{part}_{s}.json.gz"
    if not BG.guarded(p).exists():
        return dict(calls=0, pairs=0, votes=0, seconds=0.0)
    d = read_gz(p)["summed"]
    return {k: d[k] for k in ("calls", "pairs", "votes", "seconds")}


def try_counts(c: Path, s: str) -> Dict[str, Any]:
    """Rule RU's rows of the run (``ru_steps_<s>``) and the validity checks against the PXF
    wrapper's own rows and the indexer's calls."""
    rows = read_gz(c / f"ru_steps_{s}.json.gz")["rows"]
    pxf = read_gz(c / f"pxf_steps_{s}.json.gz")["rows"]
    idx = index_s(c, s)
    tried = [r for r in rows if r["tried"]]
    ok = [r for r in rows if r["success"]]
    sec = [r for r in ok if r["second"]]
    fails = [r for r in tried if not r["success"]]
    fb = sum(1 for r in pxf if r["fallback"])
    n = len(rows)
    low = [r for r in ok if r["fits"][-1]["hit"] < RR.ACCEPT]
    return dict(
        seed_steps=n,
        reseeds=sum(1 for r in rows if r["reseed"]),
        no_start=sum(1 for r in rows if not r["tried"]),
        tried=len(tried),
        success=len(ok),
        success_first_fit=len(ok) - len(sec),
        success_second_fit=len(sec),
        success_share=float(len(ok) / n) if n else None,
        success_wilson=wilson(len(ok), n) if n else None,
        failed_tries=len(fails),
        try_evals=int(sum(r["evals"] for r in tried)),
        try_evals_in_successes=int(sum(r["evals"] for r in ok)),
        try_evals_in_failures=int(sum(r["evals"] for r in fails)),
        try_seconds=float(sum(r["try_s"] for r in tried)),
        pxf_steps_after_try=len(pxf),
        pxf_fallbacks=fb,
        index_calls=int(idx["calls"]),
        index_seconds=float(idx["seconds"]),
        # validity: calls = steps - successes + second tries (the PXF fallbacks)
        check_pxf_steps_equal_steps_minus_successes=bool(len(pxf) == n - len(ok)),
        check_calls_equal=bool(idx["calls"] == n - len(ok) + fb),
        check_no_success_under_095=bool(not low),
    )


def stored_wall(job: str) -> Dict[str, Any]:
    return job_times(job, OLD / "logs")


def old_dirs(s: str) -> Dict[str, Path]:
    return dict(run=OLD / "live_pxf" / f"J5D_{s}", arm=OLD / "live_pxf" / "j5sd" / f"J5S_{s}")


def live(s: str) -> Dict[str, Any]:
    R_t = BS.truth(s)
    od = old_dirs(s)
    c = CACHE / "live_ru"
    nrun, narm = c / f"J5D_{s}", c / "j5sd" / f"J5S_{s}"
    tot = dict(old=BS.run_total(od["run"], od["arm"]), new=BS.run_total(nrun, narm))
    maps = dict(old=np.load(od["arm"] / "after.npy"), new=np.load(narm / "after.npy"))
    pol = dict(old=BS.read_json(od["run"] / "polish.json"), new=BS.read_json(nrun / "polish.json"))
    sc = maps_scores(maps, R_t)
    pr = pair(sc, "old", "new")
    stage = dict(old=BS.stage_split(od["run"], tot["old"]), new=BS.stage_split(nrun, tot["new"]))
    ow, nw = stored_wall(f"live_pxf_{s}"), job_times(f"live_ru_{s}")
    row: Dict[str, Any] = dict(
        sample=s,
        totals=tot,
        old_total_reproduces=bool(tot["old"]["total"] == RR.OLD_TOTALS[s]),
        wrong=sc["wrong"],
        median_err_deg=sc["median_err_deg"],
        new_vs_old=pr,
        grid_fallbacks_left={k: int(v["grid_fallbacks_left"]) for k, v in pol.items()},
        stages=stage,
        job=dict(new=nw, old=ow),
        run_wall_s=float(BS.read_json(nrun / "run.json")["wall_total_s"]),
        old_run_wall_s=float(BS.read_json(od["run"] / "run.json")["wall_total_s"]),
        tries=try_counts(c, s),
        index_old=index_s(OLD / "live_pxf", s),
        fewer_evals=RR.fewer(tot["new"]["total"], tot["old"]["total"]),
    )
    row["r2"] = RR.judge_r2(
        pr["net"],
        pr["median_higher"],
        row["grid_fallbacks_left"]["new"] - row["grid_fallbacks_left"]["old"],
    )
    return row


def live_all() -> Dict[str, Any]:
    rows = {s: live(s) for s in RR.FULL_SAMPLES}
    tot = sum(r["totals"]["new"]["total"] for r in rows.values())
    nw = sum(r["job"]["new"]["wall_s"] for r in rows.values())
    ow = sum(r["job"]["old"]["wall_s"] for r in rows.values())
    judge = dict(
        r1=RR.judge_r1(tot),
        r2={s: r["r2"] for s, r in rows.items()},
        r3=dict(
            RR.judge_r3(nw),
            new_wall_s=int(nw),
            old_wall_s_frozen=RR.OLD_WALL_S,
            old_wall_s_from_logs=int(ow),
            labels={s: r["job"]["new"]["label"] for s, r in rows.items()},
        ),
        old_totals_reproduce=all(r["old_total_reproduces"] for r in rows.values()),
        checks_ok=all(
            r["tries"]["check_pxf_steps_equal_steps_minus_successes"]
            and r["tries"]["check_calls_equal"]
            and r["tries"]["check_no_success_under_095"]
            for r in rows.values()
        ),
    )
    judge["all_ok"] = bool(
        judge["r1"]["ok"]
        and judge["r3"]["ok"]
        and all(v["ok"] for v in judge["r2"].values())
        and judge["old_totals_reproduce"]
        and judge["checks_ok"]
    )
    return dict(rows=rows, judge=judge)


def held8() -> Dict[str, Any]:
    s = "held8"
    c, cr = CACHE / "live_pxf", CACHE / "live_ru"
    R_t = BS.truth(s)
    arms = dict(
        pxf=(c / f"J5D_{s}", c / "j5sd" / f"J5S_{s}"),
        ru=(cr / f"J5D_{s}", cr / "j5sd" / f"J5S_{s}"),
    )
    tot = {k: BS.run_total(r, m) for k, (r, m) in arms.items()}
    s8 = BS.search_evals(BS.read_json(c / f"S8_{s}" / "run.json")["stats"]) + int(
        BS.read_json(c / f"passes_{s}" / "D4_S8" / "pass.json")["evals"]
    )
    maps = {k: np.load(BG.guarded(m / "after.npy")) for k, (r, m) in arms.items()}
    maps["s8"] = np.load(BG.guarded(c / f"passes_{s}" / "after_S8.npy"))
    sc = maps_scores(maps, R_t)
    ru_s8, ru_pxf, pxf_s8 = pair(sc, "s8", "ru"), pair(sc, "pxf", "ru"), pair(sc, "s8", "pxf")
    pol = {k: BS.read_json(r / "polish.json") for k, (r, m) in arms.items()}
    f1, f2 = RR.fewer(tot["ru"]["total"], s8), RR.fewer(tot["ru"]["total"], tot["pxf"]["total"])
    return dict(
        sample=s,
        totals=dict(s8=int(s8), pxf=int(tot["pxf"]["total"]), ru=int(tot["ru"]["total"])),
        stages={k: BS.stage_split(r, tot[k]) for k, (r, m) in arms.items()},
        wrong=sc["wrong"],
        median_err_deg=sc["median_err_deg"],
        ru_vs_s8=ru_s8,
        ru_vs_pxf=ru_pxf,
        pxf_vs_s8=pxf_s8,
        fewer_vs_s8=f1,
        fewer_vs_pxf=f2,
        grid_fallbacks_left={k: int(v["grid_fallbacks_left"]) for k, v in pol.items()},
        job=dict(
            s8=job_times("S8_held8"),
            pxf=job_times("live_pxf_held8"),
            ru=job_times("live_ru_held8"),
        ),
        tries=try_counts(cr, s),
        h1=RR.judge_h1(f1, ru_s8["net"], ru_s8["median_higher"]),
        h2=RR.judge_h2(f2, ru_pxf["net"], ru_pxf["median_higher"]),
    )


if __name__ == "__main__":
    what = sys.argv[1] if len(sys.argv) > 1 else ""
    if what == "live":
        out: Any = live_all()
        if "--write-gate" in sys.argv and out["judge"]["all_ok"]:
            (CACHE / "logs").mkdir(parents=True, exist_ok=True)
            (CACHE / "logs" / "gate_passed").write_text("r1, r2, r3 passed (ru_score.py live)\n")
    elif what == "held8":
        out = held8()
    else:
        raise SystemExit("usage: ru_score.py live [--write-gate] | held8")
    RD.mkdir(parents=True, exist_ok=True)
    (RD / f"score_{what}.json").write_text(json.dumps(out, indent=1) + "\n")
    print(json.dumps(out, indent=1)[:7000])
