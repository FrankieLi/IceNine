"""tex1: the freeze. ``frozen.json`` is committed before any method output of ``tex24`` exists.

  uv run python scripts/autoresearch/tex1/tex1_freeze.py --write [--cs-binary PATH]
  uv run python scripts/autoresearch/tex1/tex1_freeze.py --check [--full] [--outputs]

``--write`` records the seed, the texture rules, the arms with their exact settings, the slice,
the decision rules and the sha256 of the grid, the truth, the slice grid, the stack, every image
frame, every tex1 script, the scripts of earlier steps the runs import and every library file
(``icenine/*.py``). ``--check`` exits 2 unless ``frozen.json`` is committed, unchanged since, and
every recorded hash agrees (``--full`` hashes all 360 frames too; the stack is always hashed).
``--outputs`` also exits 2 if any file under the runs folder is older than the freeze commit (a
method output that existed before the freeze). No earlier step's frozen file is edited; fine1's
is read for the grid hash only.
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
import tex1_config as CF  # noqa: E402
import tex1_lib as F  # noqa: E402
import tex1_paths as P  # noqa: E402

ICE = P.ICE
OWN_SUFFIXES = (".py", ".sh", ".cpp")
EARLIER = [  # scripts of earlier steps that the tex1 runs and launchers import or call, unchanged
    "scripts/autoresearch/cpptime/make_config.py",
    "scripts/autoresearch/cpptime/guard.py",
    "scripts/autoresearch/cpptime/parse_logs.py",
    "scripts/autoresearch/cpptime/score_cpp.py",
    "scripts/autoresearch/cpptime/r4_limit.txt",
    "scripts/autoresearch/item14/quiet.py",
    "scripts/autoresearch/ar2/key_off_fused.config",
    "scripts/autoresearch/ar2/key_on_fused.config",
    "scripts/autoresearch/ar3/ar3_common.py",
    "scripts/autoresearch/ar3/held_sample_lib.py",
    "scripts/autoresearch/ar1/misses.py",
    "scripts/phase_d/grains.py",
    "scripts/phase_d/render_full.py",
    "scripts/plot_mic.py",
    "../scripts/dev/timing_preflight.py",
]
ARM_TEXT = {
    "CSP": "stock C++ in parallel (owner decision D-20261120-N6): the ar2/key_off_fused.config "
    "config through cpptime/make_config.py with SampleRadius 2.0, mpirun -n 13 (one server, 12 "
    "fitting clients), the build with the count line (cppchk IceNine_unfixed_f47bf032); run "
    "alone on all cores",
    "RCP": "the recipe of 24 Nov in library form, one worker: ar2/key_on_fused.config plus "
    + "; ".join(CF.RCP_KEYS),
}
DECISION_RULES = [
    "Two arms only (owner task of 2026-11-24): stock C++ in parallel (CSP) and the recipe (RCP).",
    "Caps: CSP 9 hours, RCP 16 hours.",
    "No arm is changed after the slice pilot; the pilot is a check of the 16 keys together.",
    "Reported: the two points as measured; no winner rule, this is a measurement.",
    "Wrong: over 1 degree from the voxel's own truth, cubic symmetry; a voxel with no output "
    "orientation counts as wrong and is listed apart.",
]
STOP_RULE = (
    "stop and report if the sample gate fails, if a run starts before tex1/frozen.json is "
    "committed, if any map or log of tex24 is read before both end markers exist (slice counts "
    "and seconds excepted), if a run hits its cap or its log stalls 30 minutes, if CSP leaves any "
    "voxel without output (report, no rerun), or if held2 is opened or listed or a "
    "pareto_rows.json read."
)


def own_scripts() -> List[Path]:
    return sorted(p for p in HERE.iterdir() if p.suffix in OWN_SUFFIXES)


def script_hashes() -> Dict[str, str]:
    out = {f"scripts/autoresearch/tex1/{p.name}": P.sha256(p) for p in own_scripts()}
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
    chk = json.loads((P.BENCH / "sample" / "sample_check.json").read_text())
    info = json.loads(P.SAMPLE_JSON.read_text())
    assert chk["gate"]["passed"], "the sample gate has not passed"
    assert P.sha256(P.TRUTH) == info["truth_sha256"], "truth differs from sample.json"
    d: Dict[str, Any] = dict(
        task="tex1",
        sample="tex24",
        seed=P.SEED,
        generator_note="a new sample (as held3 to held6), not an extra random seed for a "
        "comparison",
        grid="fine7_grid.mic of fine1 (98,280 voxels, side 0.0046875 mm; sha256 checked against "
        "fine1/frozen.json)",
        n_voxels=P.N_FINE,
        n_grains=info["n_grains"],
        texture_rule=f"{info['n_noise_grains']} noise grains (default_rng([24, 1])) and "
        f"{info['n_gradient_grains']} gradient grains (default_rng([24, 2])), drawn "
        f"independently; noise: per voxel a rotation vector uniform in the {F.NOISE_BALL_DEG} "
        f"degree ball; gradient: rotation vector linear plus {F.CURVATURE} curvature in the "
        f"position, largest edge-pair misorientation t_g uniform on [{F.GRADIENT_MIN_DEG}, "
        f"{F.GRADIENT_CAP_DEG}) degree (tex1_lib.textured_orientations)",
        n_by_class=info["n_by_class"],
        gate_result=chk["gate"],
        arms=ARM_TEXT,
        arm_order=["CSP", "RCP"],
        config_changes=list(CF.CHANGED_KEYS),
        images="render_full.py --tag tex24 --variants clean --mic <truth>, Q-max 8, clean",
        slice=dict(
            n=P.SLICE_N,
            rule=json.loads(P.SLICE_JSON.read_text())["rule"],
            grid_sha256=P.sha256(P.SLICE),
        ),
        decision_rules=DECISION_RULES,
        stop_rule=STOP_RULE,
        proposed_cs_binary=dict(
            name=Path(cs_binary).name if cs_binary else None,
            sha256=P.sha256(Path(cs_binary)) if cs_binary else None,
            note="the build of branch feature/eval-reduction-cpptrace at commit f47bf032 (develop "
            "plus prints only; it prints one cost-evaluation count line per client), kept in "
            "the main checkout under icenine_py/benchmarks/autoresearch/cppchk/cache/; the path "
            "is an argument of the launcher",
        ),
        sha256=data_hashes(full=True),
        sample_records=dict(
            sample_json=P.sha256(P.SAMPLE_JSON),
            grainmap=P.sha256(P.GRAINMAP),
            sample_check=P.sha256(P.BENCH / "sample" / "sample_check.json"),
        ),
        scripts=script_hashes(),
        library=library_hashes(),
        frozen_at=time.strftime("%Y-%m-%dT%H:%M:%S%z"),
        note="no method (CSP, RCP) was run on tex24, and no log or map of tex24 exists, before "
        "this file was committed; the RCP slice pilot runs after it",
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
        return "tex1 frozen.json does not exist"
    if not commit_time():
        return "tex1 frozen.json is not committed"
    st = subprocess.run(
        ["git", "status", "--porcelain", "--", str(P.FROZEN)],
        cwd=str(ICE),
        capture_output=True,
        text=True,
    ).stdout.strip()
    if st:
        return "tex1 frozen.json differs from its commit"
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
    if outputs and (P.RUNS / "whole").exists():
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
