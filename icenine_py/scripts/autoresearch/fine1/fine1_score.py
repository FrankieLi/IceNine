"""fine1 phase D: score the six whole-layer arms against the fine7 truth (light, one thread).

  uv run python scripts/autoresearch/fine1/fine1_score.py [--arms CS,CF,...]

Per arm, read the finished map (``fine1_cpp.mic.LBFS`` of the C++ arms, ``recon.mic.gz`` of the
Python arms), match every output voxel to the fine grid row of the same direction by left vertex
(tolerance 1e-3 mm, exactly one row, as ``score_cpp.py``; the fine grid's nearest same-direction
rows are 0.0047 mm apart), and take the cubic-symmetry misorientation to the truth ``.mic``.
Wrong: over 1 degree; a grid row with no output voxel, or with a zero orientation, counts as wrong
and is listed apart ("missing"). Wilson 95% interval (``common/stats.py``); median, 90th, 99th
percentile of the error over voxels that have an output. Evaluations: Python arms from
``run.json`` (all, and global searches apart); C++ arms the sum of the ``[cost evaluations]
client N: X`` lines of the build of branch ``feature/eval-reduction-cpptrace`` (client side cost
evaluations; the server evaluates nothing) in the launcher's log. Global searches of the C++ arms:
``First Voxel Required`` lines of the client logs that exist (parse_logs.py, unchanged). Pairs:
CSP against CS and CFP against CF (map identical or how many voxels differ, and the paired
fixed / broken wrong counts), plus the other pairs of interest. No p-values.

Writes ``benchmarks/autoresearch/fine1/score.json`` and ``score_tables.md``; per-voxel errors and
orientations go to ``<cache>/score/`` (not tracked).
"""

import argparse
import gzip
import json
import shutil
import sys
import tempfile
from pathlib import Path
from typing import Any, Dict, Optional

import numpy as np
from scipy.spatial import cKDTree

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1] / "common"))
sys.path.insert(0, str(HERE.parents[0] / "cpptime"))
import doc_tables as DT  # noqa: E402
import fine1_paths as P  # noqa: E402
import parse_logs as PL  # noqa: E402
import stats as ST  # noqa: E402

ARMS = ("CS", "CF", "CSP", "CFP", "J5S", "PS")
WORKERS = {"CS": 1, "CF": 1, "CSP": 12, "CFP": 12, "J5S": 1, "PS": 1}
WRONG_DEG = 1.0
TOL_MM = 1e-3
RUNS = P.RUNS / "whole"
OUT = P.CACHE / "score"


def read_map(arm: str) -> Any:
    from icenine.mic_file import MicFile

    run = RUNS / arm
    if arm in ("CS", "CF", "CSP", "CFP"):
        return MicFile.read(str(run / "fine1_cpp.mic.LBFS"))
    with tempfile.TemporaryDirectory() as td:
        tmp = Path(td) / "recon.mic"
        with gzip.open(run / "recon.mic.gz", "rb") as fi, open(tmp, "wb") as fo:
            shutil.copyfileobj(fi, fo)
        return MicFile.read(str(tmp))


def grid_rows() -> Any:
    a = np.loadtxt(str(P.GRID), skiprows=1, ndmin=2)
    return a[:, :2], a[:, 3].astype(int)


def match(
    out_xy: np.ndarray, out_dir: np.ndarray, g_xy: np.ndarray, g_dir: np.ndarray
) -> np.ndarray:
    """Grid row of every output voxel (exactly one candidate of the same direction)."""
    idx = np.full(len(out_xy), -1, dtype=np.int64)
    for code in (1, 2):
        rows = np.nonzero(g_dir == code)[0]
        sel = np.nonzero(out_dir == code)[0]
        if len(sel) == 0:
            continue
        d, j = cKDTree(g_xy[rows]).query(out_xy[sel], k=1)
        assert float(d.max()) < TOL_MM, f"an output voxel is {d.max()} mm from every grid row"
        idx[sel] = rows[j]
    assert (idx >= 0).all() and len(set(idx.tolist())) == len(idx), "match is not one to one"
    return idx


def score_arm(arm: str, truth_R: np.ndarray, g_xy: np.ndarray, g_dir: np.ndarray, sym: Any) -> Any:
    from icenine.ray_indexing import misorientation_deg

    mic = read_map(arm)
    n = len(g_xy)
    xy = np.array([v.position[:2] for v in mic.voxels], dtype=np.float64)
    up = np.array([1 if v.points_up else 2 for v in mic.voxels])
    R = np.array([v.orientation for v in mic.voxels], dtype=np.float64)
    idx = match(xy, up, g_xy, g_dir)
    R_all = np.zeros((n, 3, 3))
    has = np.zeros(n, dtype=bool)
    R_all[idx] = R
    has[idx] = np.abs(np.linalg.det(R)) > 1e-6
    err = np.full(n, np.inf)
    err[has] = np.asarray(misorientation_deg(R_all[has], truth_R[has], sym), dtype=float)
    wrong = ~has | (err > WRONG_DEG)
    lo, hi = ST.wilson(int(wrong.sum()), n)
    e = err[has]
    return (
        dict(
            n_voxels=n,
            n_output_voxels=int(len(idx)),
            n_missing=int((~has).sum()),
            n_wrong=int(wrong.sum()),
            n_wrong_with_output=int((has & wrong).sum()),
            wrong_fraction=float(wrong.mean()),
            wilson95=[lo, hi],
            median_err_deg=float(np.median(e)),
            p90_err_deg=float(np.percentile(e, 90)),
            p99_err_deg=float(np.percentile(e, 99)),
            n_over_0p05=int((e > 0.05).sum()),
        ),
        R_all,
        has,
        wrong,
        err,
    )


def eval_lines(arm: str) -> Any:
    """(sum, number of lines) of the client cost-evaluation lines of a C++ arm, from the tracked
    record ``cost_evaluation_lines.txt`` (the launcher logs of CSP and CFP were overwritten by a
    test; see the header of that file)."""
    ev = n = 0
    sec = ""
    for ln in (P.BENCH / "cost_evaluation_lines.txt").read_text().splitlines():
        if ln.startswith("# ") and ln.split()[1] in ("CS", "CF", "CSP", "CFP"):
            sec = ln.split()[1]
        m = PL.re.search(r"\[cost evaluations\] client (\d+):\s*(\d+)", ln)
        if m and sec == arm:
            ev += int(m.group(2))
            n += 1
    return ev, n


def cpp_counts(arm: str) -> Dict[str, Any]:
    ev, n_ev = eval_lines(arm)
    run = RUNS / arm
    clients = sorted(run.glob("ParallelReconLog_*.log"))
    gs = 0
    local = 0
    gsec = 0.0
    logs = []
    for f in clients:
        if f.name.endswith("_0.log"):
            continue
        r = PL.parse_client_log(f.read_text())
        gs += r["global_searches"]
        local += r["local_optimizations"]
        gsec += r["global_s_sum"]
        logs.append(f.name)
    t = PL.parse_time_file((run / "time.txt").read_text())
    return dict(
        evals_all=ev,
        evals_lines=n_ev,
        evals_global=None,
        global_searches_in_logs=gs,
        client_logs_found=logs,
        local_optimizations_in_logs=local,
        global_search_seconds_in_logs=gsec,
        wall_s=t["real"],
        user_s=t["user"],
        sys_s=t["sys"],
    )


def py_counts(arm: str) -> Dict[str, Any]:
    d = json.loads((RUNS / arm / "run.json").read_text())
    return dict(
        evals_all=d["evals_all"],
        evals_global=d["evals_global"],
        evals_by_key=d["evals_by_key"],
        global_searches_in_logs=d["global_searches"],
        n_seeds=d["n_seeds"],
        n_reseeds=d["n_reseeds"],
        wall_s=d["wall_total_s"],
        wall_global_s=d["wall_global_s"],
        wall_other_s=d["wall_other_s"],
        n_without_orientation=d["n_without_orientation"],
    )


def pair_stats(d: np.ndarray, wa: np.ndarray, wb: np.ndarray) -> Dict[str, Any]:
    """Counts of a map pair over the voxels BOTH maps have: ``d`` (m,) misorientation between the
    two maps, ``wa`` / ``wb`` (m,) wrong flags. A voxel one map lacks is not in any count."""
    return dict(
        n_shared=int(len(d)),
        n_same_1e3deg=int((d < 1e-3).sum()),
        n_differ_over_0p05=int((d > 0.05).sum()),
        n_differ_over_1=int((d > 1.0).sum()),
        median_diff_deg=float(np.median(d)) if len(d) else float("nan"),
        wrong_both=int((wa & wb).sum()),
        wrong_only_a=int((wa & ~wb).sum()),
        wrong_only_b=int((~wa & wb).sum()),
        wrong_neither=int((~wa & ~wb).sum()),
    )


def pair(a: str, b: str, maps: Dict[str, Any]) -> Dict[str, Any]:
    from icenine.ray_indexing import misorientation_deg, proper_rotations
    from icenine.symmetry import create_cubic_symmetry

    sym = proper_rotations(create_cubic_symmetry(3.6))
    both = maps[a]["has"] & maps[b]["has"]
    d = np.asarray(misorientation_deg(maps[a]["R"][both], maps[b]["R"][both], sym), dtype=float)
    out = pair_stats(d, maps[a]["wrong"][both], maps[b]["wrong"][both])
    out.update(
        a=a,
        b=b,
        n_lacking_in_a=int((~maps[a]["has"]).sum()),
        n_lacking_in_b=int((~maps[b]["has"]).sum()),
        note="counts over the voxels both maps have; wrong_only_a = fixed by b, only_b = broken",
    )
    return out


def time_seconds(arm: str) -> Dict[str, float]:
    """real, user, sys of ``/usr/bin/time -p`` for the launcher command of the arm."""
    t = PL.parse_time_file((P.LOGS / f"whole_{arm}.time").read_text())
    return dict(real=t["real"], user=t["user"], sys=t["sys"])


def derived(arms: Dict[str, Any]) -> Dict[str, Any]:
    """Hours, CPU hours and the ratios the docs quote (one source: ``/usr/bin/time -p`` real)."""
    out: Dict[str, Any] = {}
    for k, d in arms.items():
        t = time_seconds(k)
        out[k] = dict(
            wall_h=t["real"] / 3600.0,
            cpu_h=(t["user"] + t["sys"]) / 3600.0,
            cores_used=(t["user"] + t["sys"]) / t["real"],
            evals_millions=d["evals_all"] / 1e6,
        )
    r = lambda x, y: float(arms[x]["wall_s"] / arms[y]["wall_s"])  # noqa: E731
    out["ratios"] = dict(
        wall_CS_over_CSP=r("CS", "CSP"),
        wall_CF_over_CFP=r("CF", "CFP"),
        wall_PS_over_J5S=r("PS", "J5S"),
        evals_PS_over_J5S=float(arms["PS"]["evals_all"] / arms["J5S"]["evals_all"]),
        evals_CS_over_CSP=float(arms["CS"]["evals_all"] / arms["CSP"]["evals_all"]),
    )
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--arms", default=",".join(ARMS))
    a = ap.parse_args()
    from icenine.mic_file import MicFile
    from icenine.ray_indexing import proper_rotations
    from icenine.symmetry import create_cubic_symmetry

    sym = proper_rotations(create_cubic_symmetry(3.6))
    truth = MicFile.read(str(P.TRUTH))
    truth_R = np.array([v.orientation for v in truth.voxels], dtype=np.float64)
    g_xy, g_dir = grid_rows()
    assert len(truth_R) == len(g_xy) == P.N_FINE
    OUT.mkdir(parents=True, exist_ok=True)
    arms: Dict[str, Any] = {}
    maps: Dict[str, Any] = {}
    for arm in a.arms.split(","):
        s, R, has, wrong, err = score_arm(arm, truth_R, g_xy, g_dir, sym)
        c = cpp_counts(arm) if arm in ("CS", "CF", "CSP", "CFP") else py_counts(arm)
        c["wall_inprocess_s"] = c.get("wall_s") if arm in ("J5S", "PS") else None
        c["wall_s"] = time_seconds(arm)["real"]  # one source for every arm
        arms[arm] = dict(workers=WORKERS[arm], **s, **c)
        maps[arm] = dict(R=R, has=has, wrong=wrong)
        np.save(OUT / f"{arm}_err.npy", err)
        print(arm, s["n_wrong"], s["wilson95"], s["median_err_deg"], c["evals_all"], flush=True)
    common = np.ones(P.N_FINE, dtype=bool)
    for k in maps:
        common &= maps[k]["has"]
    for k, d in arms.items():
        w = maps[k]["wrong"][common]
        d["common_n"] = int(common.sum())
        d["common_wrong"] = int(w.sum())
        d["common_wilson95"] = list(ST.wilson(int(w.sum()), int(common.sum())))
        d["common_share"] = d["common_wrong"] / d["common_n"]
        e = np.load(OUT / f"{k}_err.npy")[common]
        d["common_median_err_deg"] = float(np.median(e))
    pairs = [
        pair(x, y, maps)
        for x, y in (
            ("CS", "CSP"),
            ("CF", "CFP"),
            ("CS", "CF"),
            ("CS", "PS"),
            ("CS", "J5S"),
            ("PS", "J5S"),
            ("CSP", "CFP"),
        )
        if x in maps and y in maps
    ]
    der = derived(arms) if len(arms) == len(ARMS) else {}
    out = dict(
        label="contended: every arm's quiet check failed (a cpptools process on one core). "
        "Wall seconds are `/usr/bin/time -p` real of the launcher command for every arm "
        "(wall_inprocess_s is the Python run's own clock). CSP and CFP ran 13 MPI ranks on a "
        "12-core machine (8 performance, 4 efficiency) with cpptools running: oversubscribed.",
        arms=arms,
        pairs=pairs,
        derived=der,
    )
    (P.BENCH / "score.json").write_text(json.dumps(out, indent=1, default=float) + "\n")
    write_tables(arms, pairs, der)
    return 0


def write_tables(arms: Dict[str, Any], pairs: Any, der: Dict[str, Any]) -> None:
    """``score_tables.md``: marker-delimited tables (``doc_tables``) the docs copy from."""
    n_c = next(iter(arms.values()))["common_n"]
    rows = []
    for k, d in arms.items():
        lo, hi = d["common_wilson95"]
        rows.append(
            dict(
                arm=k,
                workers="1 worker" if d["workers"] == 1 else "12 workers, 13 ranks on 12 cores",
                common_wrong=d["common_wrong"],
                share=f"{100 * d['common_share']:.3f}% ({100 * lo:.3f}, {100 * hi:.3f})",
                evals=d["evals_all"],
                wall_s=d["wall_s"],
                wall_h=f"{d['wall_s'] / 3600:.1f}",
            )
        )
    t: Dict[str, str] = {}
    t["fine1_whole_common"] = DT.markdown_table(
        rows,
        ["arm", "workers", "common_wrong", "share", "evals", "wall_s", "wall_h"],
        formats={"evals": ",d", "wall_s": ",.0f", "common_wrong": ",d"},
    )
    rows2 = []
    for k, d in arms.items():
        lo, hi = d["wilson95"]
        rows2.append(
            dict(
                arm=k,
                wrong=d["n_wrong"],
                no_output=d["n_missing"],
                wrong_with_output=d["n_wrong_with_output"],
                share=f"{100 * d['wrong_fraction']:.3f}% ({100 * lo:.3f}, {100 * hi:.3f})",
                median=d["median_err_deg"],
            )
        )
    t["fine1_whole_all"] = DT.markdown_table(
        rows2,
        ["arm", "wrong", "no_output", "wrong_with_output", "share", "median"],
        formats={"wrong": ",d", "no_output": ",d", "wrong_with_output": ",d", "median": ".4f"},
    )
    rows3 = [
        dict(
            pair=f"{p['a']} vs {p['b']}",
            shared=p["n_shared"],
            same=p["n_same_1e3deg"],
            over_005=p["n_differ_over_0p05"],
            over_1=p["n_differ_over_1"],
            median=p["median_diff_deg"],
            both=p["wrong_both"],
            only_a=p["wrong_only_a"],
            only_b=p["wrong_only_b"],
        )
        for p in pairs
    ]
    t["fine1_whole_pairs"] = DT.markdown_table(
        rows3,
        ["pair", "shared", "same", "over_005", "over_1", "median", "both", "only_a", "only_b"],
        formats={
            k: ",d" for k in ("shared", "same", "over_005", "over_1", "both", "only_a", "only_b")
        }
        | {"median": ".3f"},
    )
    if der:
        rows4 = [
            dict(
                arm=k,
                wall_h=der[k]["wall_h"],
                cpu_h=der[k]["cpu_h"],
                cores=der[k]["cores_used"],
                evals_m=der[k]["evals_millions"],
            )
            for k in arms
        ]
        t["fine1_whole_cpu"] = DT.markdown_table(
            rows4,
            ["arm", "wall_h", "cpu_h", "cores", "evals_m"],
            formats={"wall_h": ".1f", "cpu_h": ".1f", "cores": ".1f", "evals_m": ",.1f"},
        )
        pj = json.loads((P.BENCH / "projection.json").read_text())["arms"]
        rows5 = [
            dict(
                arm=k,
                low=pj[k]["projection"]["T_low_h"],
                high=pj[k]["projection"]["T_high_h"],
                measured=der[k]["wall_h"],
            )
            for k in ("CS", "CF", "PS", "J5S")
        ]
        t["fine1_projection_vs_measured"] = DT.markdown_table(
            rows5,
            ["arm", "low", "high", "measured"],
            formats={k: ".1f" for k in ("low", "high", "measured")},
        )
        rr = der["ratios"]
        t["fine1_ratios"] = DT.markdown_table(
            [
                dict(quantity="CS wall / CSP wall (oversubscribed)", value=rr["wall_CS_over_CSP"]),
                dict(quantity="CF wall / CFP wall (oversubscribed)", value=rr["wall_CF_over_CFP"]),
                dict(quantity="PS wall / J5S wall", value=rr["wall_PS_over_J5S"]),
                dict(quantity="PS evaluations / J5S evaluations", value=rr["evals_PS_over_J5S"]),
                dict(quantity="CS evaluations / CSP evaluations", value=rr["evals_CS_over_CSP"]),
            ],
            ["quantity", "value"],
            formats={"value": ".2f"},
        )
    DT.write_tables(P.BENCH / "score_tables.md", t)
    print(f"{n_c} common voxels")
    print((P.BENCH / "score_tables.md").read_text())


if __name__ == "__main__":
    sys.exit(main())
