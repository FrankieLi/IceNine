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
sys.path.insert(0, str(P_ := HERE.parents[1] / "common"))
sys.path.insert(0, str(HERE.parents[0] / "cpptime"))
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


def cpp_counts(arm: str) -> Dict[str, Any]:
    ev = 0
    n_ev = 0
    for ln in (P.LOGS / f"whole_{arm}.log").read_text(errors="replace").splitlines():
        m = PL.re.search(r"\[cost evaluations\] client (\d+):\s*(\d+)", ln)
        if m:
            ev += int(m.group(2))
            n_ev += 1
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


def pair(a: str, b: str, maps: Dict[str, Any]) -> Dict[str, Any]:
    from icenine.ray_indexing import misorientation_deg, proper_rotations
    from icenine.symmetry import create_cubic_symmetry

    sym = proper_rotations(create_cubic_symmetry(3.6))
    Ra, ha, wa = maps[a]["R"], maps[a]["has"], maps[a]["wrong"]
    Rb, hb, wb = maps[b]["R"], maps[b]["has"], maps[b]["wrong"]
    both = ha & hb
    d = np.full(len(Ra), np.inf)
    d[both] = np.asarray(misorientation_deg(Ra[both], Rb[both], sym), dtype=float)
    return dict(
        a=a,
        b=b,
        n_identical_1e6=int((d < 1e-6).sum()),
        n_differ_1e6=int((d >= 1e-6).sum()),
        n_differ_over_0p05=int((d > 0.05).sum()),
        n_differ_over_1=int((d > 1.0).sum()),
        maps_identical=bool((d < 1e-6).all()),
        wrong_both=int((wa & wb).sum()),
        wrong_only_a=int((wa & ~wb).sum()),
        wrong_only_b=int((~wa & wb).sum()),
        wrong_neither=int((~wa & ~wb).sum()),
        note=f"fixed by {b} against {a}: wrong_only_a; broken by {b}: wrong_only_b",
    )


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
        arms[arm] = dict(workers=WORKERS[arm], **s, **c)
        maps[arm] = dict(R=R, has=has, wrong=wrong)
        np.save(OUT / f"{arm}_err.npy", err)
        print(arm, s["n_wrong"], s["wilson95"], s["median_err_deg"], c["evals_all"], flush=True)
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
    # the voxels every arm has an output for (the C++ arms leave an edge band out)
    common = np.ones(P.N_FINE, dtype=bool)
    for k in maps:
        common &= maps[k]["has"]
    for k, d in arms.items():
        w = maps[k]["wrong"][common]
        d["common_n"] = int(common.sum())
        d["common_wrong"] = int(w.sum())
        d["common_wilson95"] = list(ST.wilson(int(w.sum()), int(common.sum())))
    out = dict(
        label="contended: every arm's quiet check failed (a cpptools process on one core); "
        "wall seconds are `/usr/bin/time -p` real of the whole launcher command",
        arms=arms,
        pairs=pairs,
    )
    (P.BENCH / "score.json").write_text(json.dumps(out, indent=1, default=float) + "\n")
    L = [
        "Reading line: one row per arm on the whole fine layer (98,280 voxels, clean images). "
        "Wrong is over 1 degree from the truth, cubic symmetry; a voxel with no output counts as "
        "wrong. All times contended.",
        "",
        "| arm | workers | evaluations | wrong | share wrong (95% Wilson) | missing "
        "| median error (deg) | wall (s) |",
        "|---|---|---|---|---|---|---|---|",
    ]
    for k, d in arms.items():
        lo, hi = d["wilson95"]
        L.append(
            f"| {k} | {d['workers']} | {d['evals_all']:,} | {d['n_wrong']} | "
            f"{100 * d['wrong_fraction']:.3f}% ({100 * lo:.3f}, {100 * hi:.3f}) | "
            f"{d['n_missing']} | {d['median_err_deg']:.4f} | {d['wall_s']:.0f} |"
        )
    n_c = next(iter(arms.values()))["common_n"]
    L += [
        "",
        f"Reading line: the same count on the {n_c} voxels that every arm has an output for.",
        "",
        "| arm | wrong on the common voxels | share (95% Wilson) |",
        "|---|---|---|",
    ]
    for k, d in arms.items():
        lo, hi = d["common_wilson95"]
        L.append(
            f"| {k} | {d['common_wrong']} | "
            f"{100 * d['common_wrong'] / n_c:.3f}% ({100 * lo:.3f}, {100 * hi:.3f}) |"
        )
    L += [
        "",
        "Reading line: voxel-by-voxel comparison of two maps; `differ` counts voxels whose two "
        "orientations are more than 1e-6 degree apart.",
        "",
        "| pair | identical | differ | differ over 0.05 deg | differ over 1 deg "
        "| wrong in both | only in first | only in second |",
        "|---|---|---|---|---|---|---|---|",
    ]
    for p in pairs:
        L.append(
            f"| {p['a']} vs {p['b']} | {p['n_identical_1e6']} | {p['n_differ_1e6']} | "
            f"{p['n_differ_over_0p05']} | {p['n_differ_over_1']} | {p['wrong_both']} | "
            f"{p['wrong_only_a']} | {p['wrong_only_b']} |"
        )
    (P.BENCH / "score_tables.md").write_text("\n".join(L) + "\n")
    print("\n".join(L))
    return 0


if __name__ == "__main__":
    sys.exit(main())
