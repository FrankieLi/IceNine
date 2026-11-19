"""AR5 step 23: scoring helpers shared by ``ball_gate.py`` and ``ball_summary.py`` (reads maps and
the truth for scoring only; nothing is tuned here). Not frozen: it applies the frozen rules to
finished runs.

``one(s)``: the numbers of arm ``J5SDC`` on sample ``s`` against ``J5SD`` (step 22's stored arm) and
``S8``: evaluations (search + screen pass + wide repair), wrong voxels (over 1 degree from the
truth, symmetry-reduced; Wilson 95%), fixed / broken by voxel, median error over all voxels, grid
searches, stage split, second-stage counts. On ``held7`` the baselines are this step's own runs
(``S8_held7``, ``J5D_held7`` and ``j5sd/J5S_held7``); on the disclosed samples they are step 22's
stored arms (links ``J5D_<s>``, ``stored_J5SD_<s>``, ``stored_S8_<s>``) and step 22's tracked
``summary.json`` (the S8 totals). Reads go through ``ball_guard.guarded`` (``held7`` stays closed
until ``logs/ball_runs.done``).
"""

import gzip
import json
import sys
from pathlib import Path
from typing import Any, Dict

import numpy as np

HERE = Path(__file__).resolve().parent
ICE = HERE.parents[2]
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(ICE / "scripts" / "common"))
import ball_guard as BG  # noqa: E402
import ball_rules as BR  # noqa: E402
import screen_rules as SR  # noqa: E402
import stats as ST  # noqa: E402

CACHE = BG.base_cache()
AR3 = ICE / "benchmarks" / "autoresearch" / "ar3"
DEFER_SUMMARY = ICE / "benchmarks" / "autoresearch" / "ar5" / "defer" / "summary.json"


def read_json(p: Path) -> Any:
    return json.loads(BG.guarded(p).read_text())


def read_gz(p: Path) -> Any:
    with gzip.open(BG.guarded(p), "rt") as f:
        return json.load(f)


def sym() -> Any:
    from icenine.ray_indexing import proper_rotations
    from icenine.symmetry import create_cubic_symmetry

    return proper_rotations(create_cubic_symmetry(3.61))


def truth(s: str) -> np.ndarray:
    from icenine.mic_file import MicFile

    d = AR3 / "gradient" / "sample" if s == "grad" else AR3 / s / "gradient" / "sample"
    mic = MicFile.read(str(BG.guarded(d / "gradient_truth.mic")))
    return np.array([v.orientation for v in mic.voxels], dtype=np.float64)


def search_evals(st: Dict[str, Any]) -> int:
    return int(
        st["n_evals_seed"] + st["n_evals_reseed"] + st["n_evals_neighbor"] + st["n_evals_revisit"]
    )


def run_total(run: Path, arm: Path) -> Dict[str, int]:
    """The search of a run and the screen pass and wide repair of its arm."""
    st = read_json(run / "run.json")["stats"]
    rep = read_json(arm / "repair.json")
    s = search_evals(st)
    return dict(
        search=s,
        screen_pass=int(rep["pass_evals"]),
        repair=int(rep["n_evals"]),
        total=int(s + rep["pass_evals"] + rep["n_evals"]),
    )


def stage_split(run: Path, tot: Dict[str, int]) -> Dict[str, int]:
    st = read_json(run / "run.json")["stats"]
    pol = read_json(run / "polish.json")
    return dict(
        index=int(st["index_evals"]),
        rescue=int(pol["rescue_evals"]),
        polish=int(pol["polish_evals"]),
        grid=int(pol["grid_evals"]),
        bfs_neighbour=int(st["n_evals_neighbor"]),
        bfs_revisit=int(st["n_evals_revisit"]),
        screen_pass=tot["screen_pass"],
        repair=tot["repair"],
    )


def sample_paths(s: str) -> Dict[str, Path]:
    """Where the baselines of ``s`` live (this step's own on ``held7``)."""
    if s == "held7":
        return dict(
            base_run=CACHE / "J5D_held7",
            base_arm=CACHE / "j5sd" / "J5S_held7",
            base_map=CACHE / "j5sd" / "J5S_held7" / "after.npy",
            s8_map=CACHE / "passes_held7" / "after_S8.npy",
        )
    return dict(
        base_run=CACHE / f"J5D_{s}",
        base_arm=CACHE / f"stored_J5SD_{s}",
        base_map=CACHE / f"stored_J5SD_{s}" / "after.npy",
        s8_map=CACHE / f"stored_S8_{s}" / "after.npy",
    )


def s8_total(s: str) -> int:
    if s == "held7":
        st = read_json(CACHE / "S8_held7" / "run.json")["stats"]
        return search_evals(st) + int(
            read_json(CACHE / "passes_held7" / "D4_S8" / "pass.json")["evals"]
        )
    return int(read_json(DEFER_SUMMARY)["rows"][s]["totals"]["s8"])


def one(s: str) -> Dict[str, Any]:
    p = sample_paths(s)
    run_c = CACHE / f"J5DC_{s}"
    arm_c = CACHE / "j5sdc" / f"J5S_{s}"
    y = sym()
    R_t = truth(s)
    tot_c = run_total(run_c, arm_c)
    tot_b = run_total(p["base_run"], p["base_arm"])
    s8 = s8_total(s)
    m_c = np.load(BG.guarded(arm_c / "after.npy"))
    m_b = np.load(BG.guarded(p["base_map"]))
    m_8 = np.load(BG.guarded(p["s8_map"]))
    maps = dict(j5sdc=m_c, j5sd=m_b, s8=m_8)
    w = {k: SR.wrong_mask(m, R_t, y) for k, m in maps.items()}
    err = {k: np.asarray(SR.misorientation_deg(m, R_t, y), float) for k, m in maps.items()}
    n = len(w["j5sdc"])
    wrong = {}
    for k in w:
        lo, hi = ST.wilson(int(w[k].sum()), n)
        wrong[k] = dict(n=int(w[k].sum()), wilson=[float(lo), float(hi)])
    pair = {}
    for k in ("j5sd", "s8"):
        f, b = SR.fixed_broken(w[k], w["j5sdc"])
        pair[k] = dict(
            fixed=int(f),
            broken=int(b),
            net=int(b - f),
            median_higher=float(np.median(err["j5sdc"]) - np.median(err[k])),
        )
    tot = dict(s8=int(s8), j5sd=tot_b["total"], j5sdc=tot_c["total"])
    fewer = {k: BR.fewer(tot["j5sdc"], tot[k]) for k in ("j5sd", "s8")}
    pol_c = read_json(run_c / "polish.json")
    pol_b = read_json(p["base_run"] / "polish.json")
    ball = read_json(run_c / "ball.json")
    hi5 = read_json(run_c / "high500.json")
    dfr = read_json(run_c / "defer.json")
    s22 = None
    if s != "held7":
        s22 = int(read_json(DEFER_SUMMARY)["rows"][s]["totals"]["j5sd"])
    row: Dict[str, Any] = dict(
        sample=s,
        j5sd_total_matches_step22=(None if s22 is None else bool(s22 == tot["j5sd"])),
        totals=tot,
        fewer=fewer,
        wrong=wrong,
        median_err_deg={k: float(np.median(v)) for k, v in err.items()},
        pair=pair,
        stages=stage_split(run_c, tot_c),
        stages_j5sd=stage_split(p["base_run"], tot_b),
        grid=dict(j5d=int(pol_b["grid_fallbacks_left"]), j5dc=int(pol_c["grid_fallbacks_left"])),
        rescues=dict(
            j5d_attempted=int(pol_b["rescues_attempted"]),
            j5d_rescued=int(pol_b["rescued"]),
            j5dc_attempted=int(pol_c["rescues_attempted"]),
            j5dc_rescued=int(pol_c["rescued"]),
        ),
        ball={k: ball[k] for k in ball if k not in ("rule", "name")},
        checks=int(hi5.get("checks_total", 0)),
        pauses=int(dfr.get("pauses", 0)),
        vs_j5sd=BR.judge_vs_j5sd(
            s, fewer["j5sd"], pair["j5sd"]["net"], pair["j5sd"]["median_higher"]
        ),
        vs_s8=BR.judge_vs_s8(fewer["s8"], pair["s8"]["net"], pair["s8"]["median_higher"]),
    )
    return row
