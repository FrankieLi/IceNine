"""AR5 step 27: scoring of finished runs against the frozen gates (reads maps and the truth for
scoring only; nothing is tuned here). Not frozen: it applies the frozen rules to finished runs.

  uv run python scripts/autoresearch/ar5/px_score.py q1|t1|live|h1

``q1``: reads ``q1_<sample>.json`` (8 samples). ``t1``: reads ``t1_<sample>.json``. ``live``:
for ``mid`` and ``held``, the new J5SD (``live_new``) and the old J5SD run again (``live_old``)
against stored J5SD (step 22) : evaluations (search + screen pass + wide repair), wrong voxels
(over 1 degree from the truth, symmetry reduced; Wilson 95%), fixed / broken by voxel, median
error, grid searches, wall seconds with the preflight label and the load, indexer work from the
wrappers' logs. ``h1``: the new J5SD on ``held7`` against ``S8`` (its logs and maps open only
after ``logs/bf_runs.done``: reads go through ``bf_guard.guarded``).
"""

import gzip
import json
import os
import sys
from pathlib import Path
from typing import Any, Dict, List

import numpy as np

HERE = Path(__file__).resolve().parent
ICE = HERE.parents[2]
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(ICE / "scripts" / "common"))


def _default_cache() -> Path:
    parts = ICE.parts
    if ".claude" in parts and "worktrees" in parts:
        main = Path(*parts[: parts.index(".claude")])
        return main / "icenine_py" / "benchmarks" / "autoresearch" / "ar5" / "px" / "cache"
    return ICE / "benchmarks" / "autoresearch" / "ar5" / "px" / "cache"


os.environ.setdefault("PX_CACHE", str(_default_cache()))
for _k in ("SCREEN_CACHE", "DEFER_CACHE", "BALL_CACHE", "BF_CACHE"):
    os.environ.setdefault(_k, os.environ["PX_CACHE"])

import bf_guard as BG  # noqa: E402
import bf_score as BS  # noqa: E402
import px_rules as PR  # noqa: E402
import screen_rules as SR  # noqa: E402
import stats as ST  # noqa: E402

CACHE = Path(os.environ["PX_CACHE"])
RD = ICE / "benchmarks" / "autoresearch" / "ar5" / "px"
DEFER_SUMMARY = BS.DEFER_SUMMARY


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
        p = c / f"px_index_{part}_{s}.json.gz"
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


def live(s: str) -> Dict[str, Any]:
    R_t = BS.truth(s)
    stored_run, stored_arm = CACHE / f"J5D_{s}", CACHE / f"stored_J5SD_{s}"
    arms = {
        "old": (CACHE / "live_old" / f"J5D_{s}", CACHE / "live_old" / "j5sd" / f"J5S_{s}"),
        "new": (CACHE / "live_new" / f"J5D_{s}", CACHE / "live_new" / "j5sd" / f"J5S_{s}"),
    }
    tot = {"stored": BS.run_total(stored_run, stored_arm)}
    maps = {"stored": np.load(BG.guarded(stored_arm / "after.npy"))}
    pol = {"stored": BS.read_json(stored_run / "polish.json")}
    st = {"stored": BS.read_json(stored_run / "run.json")["stats"]}
    for a, (r, m) in arms.items():
        tot[a] = BS.run_total(r, m)
        maps[a] = np.load(BG.guarded(m / "after.npy"))
        pol[a] = BS.read_json(r / "polish.json")
        st[a] = BS.read_json(r / "run.json")["stats"]
    sc = maps_scores(maps, R_t)
    row: Dict[str, Any] = dict(
        sample=s,
        totals={k: v for k, v in tot.items()},
        stored_total_matches_step22=bool(tot["stored"]["total"] == PR.STORED_TOTALS[s]),
        old_total_equals_stored=bool(tot["old"]["total"] == tot["stored"]["total"]),
        old_map_equals_stored=bool(np.array_equal(maps["old"], maps["stored"])),
        wrong=sc["wrong"],
        median_err_deg=sc["median_err_deg"],
        new_vs_stored=pair(sc, "stored", "new"),
        new_vs_old=pair(sc, "old", "new"),
        old_vs_stored=pair(sc, "stored", "old"),
        grid_fallbacks_left={k: int(v["grid_fallbacks_left"]) for k, v in pol.items()},
        index_fallback_seed_test={k: int(v["index_fallback_seed_test"]) for k, v in st.items()},
        n_seeds={k: int(v["index_runs"]) for k, v in st.items()},
        run_index_s={k: float(v["index_s"]) for k, v in st.items()},
        run_wall_s={
            a: float(BS.read_json(r / "run.json")["wall_total_s"]) for a, (r, m) in arms.items()
        },
        job=dict(old=job_times(f"live_old_{s}"), new=job_times(f"live_new_{s}")),
        index_work=dict(
            old=index_work(CACHE / "live_old", s), new=index_work(CACHE / "live_new", s)
        ),
    )
    row["c3"] = PR.judge_c3(
        row["new_vs_stored"]["net"],
        row["new_vs_stored"]["median_higher"],
        row["grid_fallbacks_left"]["new"] - row["grid_fallbacks_left"]["stored"],
    )
    ow, nw = row["job"]["old"]["wall_s"], row["job"]["new"]["wall_s"]
    row["wall_fewer"] = PR.fewer(nw, ow) if ow and nw else None
    return row


def h1() -> Dict[str, Any]:
    s = "held7"
    c = CACHE / "live_new"
    R_t = BS.truth(s)
    run, arm = c / f"J5D_{s}", c / "j5sd" / f"J5S_{s}"
    tot = BS.run_total(run, arm)
    s8 = BS.search_evals(BS.read_json(c / f"S8_{s}" / "run.json")["stats"]) + int(
        BS.read_json(c / "passes_held7" / "D4_S8" / "pass.json")["evals"]
    )
    maps = dict(
        j5sd_new=np.load(BG.guarded(arm / "after.npy")),
        s8=np.load(BG.guarded(c / "passes_held7" / "after_S8.npy")),
    )
    sc = maps_scores(maps, R_t)
    pr = pair(sc, "s8", "j5sd_new")
    fe = PR.fewer(tot["total"], s8)
    pol = BS.read_json(run / "polish.json")
    row = dict(
        sample=s,
        totals=dict(s8=int(s8), j5sd_new=int(tot["total"])),
        stages=BS.stage_split(run, tot),
        wrong=sc["wrong"],
        median_err_deg=sc["median_err_deg"],
        new_vs_s8=pr,
        fewer=fe,
        grid_fallbacks_left=int(pol["grid_fallbacks_left"]),
        job=dict(new=job_times("live_new_held7"), s8=job_times("S8_held7")),
        index_work=index_work(c, s),
        h1=PR.judge_h1(fe, pr["net"], pr["median_higher"]),
    )
    return row


def q1_all() -> Dict[str, Any]:
    rows = {}
    tl = tw = tg = 0
    for s in PR.Q1_SAMPLES:
        d = json.loads((RD / f"q1_{s}.json").read_text())
        rows[s] = {k: d[k] for k in ("n", "was_in_top20", "lost", "gained", "wall_s")}
        rk = {"base": [], "new": []}
        for r in d["rows"]:
            for k in rk:
                rk[k].append(r[k]["rank"])
        for k in rk:
            a = np.array(rk[k])
            rows[s][f"top1_{k}"] = int((a == 0).sum())
            rows[s][f"top5_{k}"] = int(((a >= 0) & (a < 5)).sum())
            rows[s][f"top20_{k}"] = int((a >= 0).sum())
        rows[s]["pairs_base"] = int(sum(r["base"]["n_pairs"] for r in d["rows"]))
        rows[s]["pairs_new"] = int(sum(r["new"]["n_pairs"] for r in d["rows"]))
        rows[s]["votes_base"] = int(sum(r["base"]["n_votes"] for r in d["rows"]))
        rows[s]["votes_new"] = int(sum(r["new"]["n_votes"] for r in d["rows"]))
        rows[s]["s_base"] = float(sum(r["base"]["s"] for r in d["rows"]))
        rows[s]["s_new"] = float(sum(r["new"]["s"] for r in d["rows"]))
        tl += d["lost"]
        tw += d["was_in_top20"]
        tg += d["gained"]
    lo, hi = ST.wilson(tl, tw)
    return dict(rows=rows, judge=PR.judge_q1(tl, tw, tg), lost_wilson=[float(lo), float(hi)])


def t1_all() -> Dict[str, Any]:
    rows = {}
    so = sn = po = pn = vo = vn = 0.0
    for s in PR.FULL_SAMPLES:
        d = json.loads((RD / f"t1_{s}.json").read_text())
        rows[s] = {
            k: d[k]
            for k in (
                "n_voxels",
                "old",
                "new",
                "less_pct_seconds",
                "less_pct_pairs",
                "less_pct_votes",
                "verdict_start",
                "verdict_end",
                "timing_label",
                "wall_s",
            )
        }
        so += d["old"]["s"]
        sn += d["new"]["s"]
        po += d["old"]["pairs"]
        pn += d["new"]["pairs"]
        vo += d["old"]["votes"]
        vn += d["new"]["votes"]
    return dict(
        rows=rows,
        pooled=dict(
            old_s=so,
            new_s=sn,
            old_pairs=int(po),
            new_pairs=int(pn),
            old_votes=int(vo),
            new_votes=int(vn),
        ),
        judge=PR.judge_t(PR.fewer(sn, so), 0.50),
    )


if __name__ == "__main__":
    what = sys.argv[1] if len(sys.argv) > 1 else ""
    if what == "q1":
        out = q1_all()
    elif what == "t1":
        out = t1_all()
    elif what == "live":
        out = {s: live(s) for s in PR.FULL_SAMPLES}
    elif what == "h1":
        out = h1()
    else:
        raise SystemExit("usage: px_score.py q1|t1|live|h1")
    (RD / f"score_{what}.json").write_text(json.dumps(out, indent=1) + "\n")
    print(json.dumps(out, indent=1)[:6000])
