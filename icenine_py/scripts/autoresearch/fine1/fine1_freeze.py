"""fine1: the freeze. ``frozen.json`` is committed before any method output of ``fine7`` exists.

  uv run python scripts/autoresearch/fine1/fine1_freeze.py --write [--cs-binary PATH]
  uv run python scripts/autoresearch/fine1/fine1_freeze.py --check [--full] [--outputs]

``--write`` records the seed, the grid rule, the cap rule, the arms with their exact settings, the
slice, the projection formula, the decision rules and the sha256 of the grid, the truth, the slice
grid, the stack, every image frame, every fine1 script, the scripts of earlier steps the runs
import and every library file (``icenine/*.py``). ``--check`` exits 2 unless ``frozen.json`` is
committed, unchanged since, and every recorded hash agrees (``--full`` hashes all 360 frames too;
the stack is always hashed). ``--outputs`` also exits 2 if any file under the runs folder is older
than the freeze commit (a method output that existed before the freeze). No earlier step's frozen
file is read or edited.
"""

import argparse
import json
import subprocess
import sys
import time
from pathlib import Path
from typing import Any, Dict, List

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import fine1_config as CF  # noqa: E402
import fine1_lib as F  # noqa: E402
import fine1_paths as P  # noqa: E402

ICE = P.ICE
OWN_SUFFIXES = (".py", ".sh", ".cpp")
EARLIER = [  # scripts of earlier steps that the fine1 runs and launchers import or call, unchanged
    "scripts/autoresearch/cpptime/make_config.py",
    "scripts/autoresearch/cpptime/guard.py",
    "scripts/autoresearch/cpptime/parse_logs.py",
    "scripts/autoresearch/cpptime/score_cpp.py",
    "scripts/autoresearch/item14/quiet.py",
    "scripts/autoresearch/ar2/key_off_fused.config",
    "scripts/autoresearch/ar2/key_on_fused.config",
    "scripts/autoresearch/ar3/ar3_common.py",
    "scripts/phase_d/render_full.py",
    "../scripts/dev/timing_preflight.py",
]
ARM_TEXT = {
    "PS": "stock Python: grid search for every seed, stock MC local optimization "
    "(ar2/key_off_fused.config: BFSRevisitRefit 1, BFSReseed 1, FusedCost 1)",
    "CS": "stock C++: mode p, mpirun -n 2 (one server, one fitting client), the PS config through "
    "cpptime/make_config.py",
    "J5S": "latest recipe, from library keys: ar2/key_on_fused.config plus "
    + "; ".join(CF.J5S_KEYS),
    "CF": "CS with part 1's fixed C++ build, an extra labelled run (only if a fixed binary is "
    "given)",
    "CSP": "CS with all cores as fitting workers (owner decision D-20261120-N6): the same binary "
    "and the same config text as CS, mpirun -n 13 (one server, 12 clients); run alone, later",
    "CFP": "CF with all cores as fitting workers (D-20261120-N6): as CSP with the corrected binary",
}
DECISION_RULES = [
    "If any arm's high projection passes 48 hours, stop and ask the owner before starting it.",
    "Caps: 1.5 times each arm's high projection.",
    "No arm is changed after the slice.",
    "Reported: the three (or four) points as measured; no winner rule, this is a measurement.",
    "Wrong: over 1 degree from the truth, cubic symmetry; a voxel with no output orientation "
    "counts as wrong and is listed apart.",
]
PROJECTION = (
    "T = t_g x G + t_r x (98,280 / 2,000); t_g the slice's mean seconds per global search, G the "
    "projected global searches, t_r the slice's seconds outside global searches. Low G = slice "
    "searches x (grains in the layer / grains touching the slice), from grain-map counts only; "
    "high G = slice searches x 49.14."
)
STOP_RULE = (
    "stop and report if the gradient gate fails, if the grid check fails, if a run starts before "
    "fine1/frozen.json is committed, if any log or map of fine7 is read before all end markers "
    "exist, if a projection passes 48 hours (ask the owner), if a run hits its cap or its log "
    "stalls for 30 minutes, or if a timed run's quiet check fails (relabel as contended and "
    "report; no rerun)."
)


def own_scripts() -> List[Path]:
    return sorted(p for p in HERE.iterdir() if p.suffix in OWN_SUFFIXES)


def script_hashes() -> Dict[str, str]:
    out = {f"scripts/autoresearch/fine1/{p.name}": P.sha256(p) for p in own_scripts()}
    for rel in EARLIER:
        out[rel] = P.sha256(ICE / rel)
    return out


def library_hashes() -> Dict[str, str]:
    return {f"icenine/{p.name}": P.sha256(p) for p in sorted((ICE / "icenine").glob("*.py"))}


def data_hashes(full: bool) -> Dict[str, Any]:
    out: Dict[str, Any] = dict(
        grid=P.sha256(P.GRID),
        truth=P.sha256(P.TRUTH),
        slice_grid=P.sha256(P.SLICE),
        stack=P.sha256(P.STACK),
    )
    if full:
        out["frames"] = {p.name: P.sha256(p) for p in sorted(P.FRAMES.iterdir())}
    return out


def write(cs_binary: str, refreeze_of: str = "", refreeze_note: str = "") -> None:
    gc = json.loads(P.GRID_CHECK.read_text())
    gr = json.loads(P.GRADIENT_CHECK.read_text())
    info = json.loads(P.SAMPLE_JSON.read_text())
    assert gc["ok"] and gr["gate"]["passed"], "grid check or gradient gate not passed"
    d: Dict[str, Any] = dict(
        task="fine1",
        sample="fine7",
        seed=P.SEED,
        generator_note="a new sample (as held3 to held6), not an extra random seed for a "
        "comparison",
        grid_rule="every triangle of SimInput/rand_500grains_1mm_neworient_s0_grid.mic split into "
        "four as C++ SubdivideTriangles does (fine1_lib.split_once): generation 7, side "
        "0.0046875 mm, 98,280 voxels, same area; C++ MinSideLength and MaxInitSideLength 0.0046875 "
        "find generation 7 so nothing is split again",
        n_voxels=P.N_FINE,
        n_grains=info["n_grains"],
        cap_rule=f"after a grain's gradient is drawn, if its largest neighbour-pair "
        f"misorientation (edge-sharing voxels of one grain) is {F.GRADIENT_CAP_DEG} degree or "
        "more, "
        "its rotation vectors are scaled down so that it is that value (fine1_lib."
        "capped_gradient_orientations); capped grains listed in sample/caps.json",
        n_capped=info["n_capped"],
        gate=dict(
            pair_max_below_deg=F.GATE_PAIR_MAX_DEG, spread_at_least_deg=F.GATE_SPREAD_MIN_DEG
        ),
        gate_result=gr["gate"],
        arms=ARM_TEXT,
        arm_order=["CS", "J5S", "PS", "CF"],
        config_changes=list(CF.CHANGED_KEYS),
        images="render_full.py --tag fine7 --variants clean --mic <truth>, Q-max 8, clean",
        slice=dict(
            n=P.SLICE_N,
            rule=json.loads(P.SLICE_JSON.read_text())["rule"],
            grid_sha256=P.sha256(P.SLICE),
        ),
        projection=PROJECTION,
        decision_rules=DECISION_RULES,
        stop_rule=STOP_RULE,
        proposed_cs_binary=dict(
            name=Path(cs_binary).name if cs_binary else None,
            sha256=P.sha256(Path(cs_binary)) if cs_binary else None,
            note="the build of branch feature/eval-reduction-cpptrace at commit f47bf032 (it "
            "prints the client cost-evaluation count; its map is byte-equal to the stored "
            "binary's on the pilot), kept in the main checkout under icenine_py/benchmarks/"
            "autoresearch/cppchk/cache/; the path is an argument of the launcher",
        ),
        sha256=data_hashes(full=True),
        scripts=script_hashes(),
        library=library_hashes(),
        frozen_at=time.strftime("%Y-%m-%dT%H:%M:%S%z"),
        note="no method (PS, CS, J5S, CF) was run on fine7 or its slice, and no log or map of "
        "fine7 exists, before this file was committed",
    )
    if refreeze_of:
        old = json.loads(P.FROZEN.read_text())
        d["proposed_cs_binary"] = old["proposed_cs_binary"]
        hist = list(old.get("refreeze_history", []))
        if "refreeze" in old:
            hist.append(old["refreeze"])
        hist.append(dict(previous_freeze_commit=refreeze_of, note=refreeze_note))
        d["refreeze_history"] = hist
        d["refreeze"] = hist[-1]
    P.FROZEN.write_text(json.dumps(d, indent=1) + "\n")
    print("wrote", P.FROZEN)


def commit_time(first: bool = False) -> int:
    """Commit time of frozen.json: the latest commit, or with ``first`` the very first one (an
    amendment re-commits the file; outputs of arms run under an earlier freeze stay valid)."""
    out = subprocess.run(
        ["git", "log", "--format=%ct", "--", str(P.FROZEN)],
        cwd=str(ICE),
        capture_output=True,
        text=True,
    ).stdout.split()
    if not out:
        return 0
    return int(out[-1] if first else out[0])


def problems(full: bool, outputs: bool) -> str:
    """'' when frozen.json is committed, clean and agrees with the files; else the reason."""
    if not P.FROZEN.exists():
        return "fine1 frozen.json does not exist"
    if not commit_time():
        return "fine1 frozen.json is not committed"
    st = subprocess.run(
        ["git", "status", "--porcelain", "--", str(P.FROZEN)],
        cwd=str(ICE),
        capture_output=True,
        text=True,
    ).stdout.strip()
    if st:
        return "fine1 frozen.json differs from its commit"
    saved = json.loads(P.FROZEN.read_text())
    for name, cur in (("scripts", script_hashes()), ("library", library_hashes())):
        for k, v in cur.items():
            if saved[name].get(k) != v:
                return f"{k} differs from the frozen hash"
        for k in saved[name]:
            if k not in cur:
                return f"{k} is frozen but missing"
    cur_d = data_hashes(full)
    for k, v in cur_d.items():
        if k == "frames":
            if v != saved["sha256"]["frames"]:
                return "an image frame differs from the frozen hash"
        elif saved["sha256"][k] != v:
            return f"{k} differs from the frozen hash"
    if outputs and (P.RUNS / "whole").exists():  # the slice runs predate a re-freeze, by design
        t = commit_time(first=True)
        for p in (P.RUNS / "whole").rglob("*"):
            if p.is_file() and int(p.stat().st_mtime) < t:
                return f"{p} is older than the freeze commit"
    return ""


def check_or_exit(full: bool = False, outputs: bool = False) -> None:
    p = problems(full, outputs)
    if p:
        print("freeze check failed:", p, file=sys.stderr)
        raise SystemExit(2)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--write", action="store_true")
    ap.add_argument("--check", action="store_true")
    ap.add_argument("--full", action="store_true", help="--check: hash every image frame")
    ap.add_argument("--outputs", action="store_true", help="--check: no output older than freeze")
    ap.add_argument("--cs-binary", default="", help="--write: the C++ binary proposed for CS")
    ap.add_argument("--refreeze-of", default="", help="--write: commit of the freeze replaced")
    ap.add_argument("--refreeze-note", default="", help="--write: why and what changed")
    a = ap.parse_args()
    if a.write:
        write(a.cs_binary, a.refreeze_of, a.refreeze_note)
    elif a.check:
        check_or_exit(a.full, a.outputs)
        print("freeze check passed")
    else:
        raise SystemExit("use --write or --check")


if __name__ == "__main__":
    main()
