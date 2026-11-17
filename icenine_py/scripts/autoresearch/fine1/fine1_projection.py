"""fine1 phase B: read the slice runs (counts and seconds only) and project the whole layer.

  uv run python scripts/autoresearch/fine1/fine1_projection.py

Reads, per arm, ``<cache>/runs/slice/<arm>/``: for the Python arms ``run.json`` (counts and seconds
of the BFS; no map, no record); for the C++ arms the client log's global searches, local
optimizations and timers (``cpptime/parse_logs.py``, imported unchanged) and the evaluation count
line ``[cost evaluations] client N: X`` that the build of branch ``feature/eval-reduction-cpptrace``
prints to the run's stdout (the launcher's ``<cache>/logs/slice_<arm>.log``). The grain map is read
for grain counts only (grains in the layer, grains touching the slice). No orientation, wrong count
or truth is read.

The frozen formula: T = t_g x G + t_r x (98,280 / 2,000); t_g the slice's mean seconds per global
search, G the projected global searches, t_r the slice's seconds outside global searches. Low G =
slice searches x (grains in the layer / grains touching the slice); high G = slice searches x
49.14. Cap = 1.5 x the high projection. Writes ``benchmarks/autoresearch/fine1/projection.json``
and ``projection_tables.md``.
"""

import json
import re
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[0] / "cpptime"))
import fine1_paths as P  # noqa: E402
import parse_logs as PL  # noqa: E402  (cpptime, unchanged)

SCALE = P.N_FINE / P.SLICE_N  # 49.14
EVAL_LINE = re.compile(r"\[cost evaluations\] client \d+:\s*(\d+)")


def parse_eval_count(text: str) -> Optional[int]:
    """Sum of the ``[cost evaluations] client N: X`` lines of a C++ run's stdout, or None."""
    vals = [int(m.group(1)) for m in EVAL_LINE.finditer(text)]
    return int(sum(vals)) if vals else None


def project(
    t_g: float, n_slice_global: int, t_r: float, n_grains: int, n_touch: int
) -> Dict[str, Any]:
    """The frozen projection: hours at the low and the high G, and the cap (1.5 x high)."""
    g_low = n_slice_global * n_grains / n_touch
    g_high = n_slice_global * SCALE
    t_low = (t_g * g_low + t_r * SCALE) / 3600.0
    t_high = (t_g * g_high + t_r * SCALE) / 3600.0
    return dict(
        t_g_s=t_g,
        slice_global_searches=n_slice_global,
        t_r_s=t_r,
        G_low=g_low,
        G_high=g_high,
        T_low_h=t_low,
        T_high_h=t_high,
        cap_h=1.5 * t_high,
        passes_48h_low=bool(t_low > 48.0),
        passes_48h_high=bool(t_high > 48.0),
    )


def grain_counts() -> Dict[str, int]:
    gmap = np.load(P.GRAINMAP)["grain"]
    rows = np.array(json.loads(P.SLICE_JSON.read_text())["rows"])
    return dict(n_grains=int(len(np.unique(gmap))), n_touch=int(len(np.unique(gmap[rows]))))


def read_python(run: Path) -> Dict[str, Any]:
    d = json.loads((run / "run.json").read_text())
    return dict(
        kind="python",
        wall_total_s=d["wall_total_s"],
        global_searches=d["global_searches"],
        n_seeds=d["n_seeds"],
        n_reseeds=d["n_reseeds"],
        wall_global_s=d["wall_global_s"],
        wall_other_s=d["wall_other_s"],
        evals_all=d["evals_all"],
        evals_global=d["evals_global"],
        evals_by_key=d["evals_by_key"],
        n_without_orientation=d["n_without_orientation"],
        timing_label=d["timing_label"],
    )


def read_cpp(run: Path, arm: str, log_dir: Path) -> Dict[str, Any]:
    r = PL.parse_run(run, run / "time.txt")
    ev_log = log_dir / f"slice_{arm}.log"
    ev = parse_eval_count(ev_log.read_text(errors="replace")) if ev_log.exists() else None
    return dict(
        kind="cpp",
        wall_total_s=float(r["time_real"]),
        client_total_s=r["total_processing_s"],
        global_searches=r["global_searches"],
        rejected_seeds=r["rejected_seeds"],
        local_optimizations=r["local_optimizations"],
        wall_global_s=r["global_s_sum"],
        wall_other_s=float(r["total_processing_s"]) - r["global_s_sum"],
        read_s=r.get("read_s"),
        evals_all=ev,
        eval_count_line="present" if ev is not None else "missing",
    )


def main() -> int:
    runs = P.RUNS / "slice"
    gc = grain_counts()
    arms: Dict[str, Any] = {}
    for arm in ("CS", "J5S", "PS", "CF"):
        d = read_cpp(runs / arm, arm, P.LOGS) if arm in ("CS", "CF") else read_python(runs / arm)
        d["marker"] = (runs / f"{arm}.done").read_text().strip()
        n = d["global_searches"]
        d["projection"] = project(d["wall_global_s"] / n, n, d["wall_other_s"], **gc)
        arms[arm] = d
    # the first PS try (stopped at the launcher's default 1 h cap): counts and seconds only
    first: Dict[str, Any] = {}
    prog = sorted(runs.glob("interrupted_*/PS/progress.jsonl"))
    if prog:
        lines = [json.loads(x) for x in prog[-1].read_text().splitlines() if x.strip()]
        last = lines[-1]
        st = last["stats"]
        first = dict(
            progress_lines=len(lines),
            elapsed_s_at_last_line=last["elapsed_s"],
            n_records_at_last_line=last["n_records"],
            n_seeds=st.get("n_seeds"),
            n_reseeds=st.get("n_reseeds"),
            marker=(runs / "PS.cap1h.marker").read_text().strip(),
        )
    out = dict(
        formula="T = t_g x G + t_r x 98,280 / 2,000",
        scale=SCALE,
        grains=gc,
        label="contended: the timing preflight refused (a cpptools process at 100% CPU); CS also "
        "shared the machine with a C++ rerun, the other arms only with cpptools",
        arms=arms,
        ps_first_try_capped_at_1h=first,
    )
    (P.BENCH / "projection.json").write_text(json.dumps(out, indent=1, default=float) + "\n")
    L: List[str] = [
        "Reading line: per arm, from its slice run (2,000 voxels, contended): global searches, "
        "seconds per global search t_g, seconds outside them t_r, the projected searches G (low, "
        "high), hours T for the whole layer (low, high) and the cap (1.5 x high).",
        "",
        "| arm | slice s | evaluations | global searches | t_g (s) | t_r (s) | G low | G high "
        "| T low (h) | T high (h) | cap (h) |",
        "|---|---|---|---|---|---|---|---|---|---|---|",
    ]
    for arm, d in arms.items():
        p = d["projection"]
        ev = "missing" if d["evals_all"] is None else f"{d['evals_all']:,}"
        L.append(
            f"| {arm} | {d['wall_total_s']:.0f} | {ev} | {d['global_searches']} | "
            f"{p['t_g_s']:.2f} | {p['t_r_s']:.1f} | {p['G_low']:.0f} | {p['G_high']:.0f} | "
            f"{p['T_low_h']:.1f} | {p['T_high_h']:.1f} | {p['cap_h']:.1f} |"
        )
    (P.BENCH / "projection_tables.md").write_text("\n".join(L) + "\n")
    print("\n".join(L))
    print("grains", gc, "first PS try", first)
    return 0


if __name__ == "__main__":
    sys.exit(main())
