#!/usr/bin/env python3
"""C++ on the E0 cases (per-voxel images of scripts/findoptimal_robustness, Q-max 8 config).

For each (voxel, variant) the stored pixel keys are written as C++ `.d` frames ("col, row,
intensity", one file per frame and detector), a one-voxel .mic is cut from the ManyGrains mic, and
the scratch debug build (patches 1-3, `CTRACE=1`) runs the adaptive reconstruction of that voxel.
The log gives the per-level candidate lists (CD, CQ, CK), the FindOptimal candidates (CF) and the
result (CR, CV). Output: scripts/phase_d/cache/e0_cpp/<voxel>_<variant>.json (gitignored).

  E0_CPP_BIN=<scratch>/IceNine E0_CPP_SCRATCH=<dir> \\
    uv run python scripts/phase_d/e0_cpp.py run [--workers 8] [--only-file LIST] [--limit N]
  uv run python scripts/phase_d/e0_cpp.py summary
"""

import argparse
import json
import multiprocessing as mp
import os
import re
import subprocess
import sys
import time
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import numpy as np

HERE = Path(__file__).resolve().parent
ICE = HERE.parents[1]
EX = ICE.parent / "Examples" / "Example2.ManyGrains"
sys.path.insert(0, str(ICE / "scripts" / "findoptimal_robustness"))
sys.path.insert(0, str(ICE / "scripts" / "common"))
SEED = os.environ.get("CSEED", "")  # patch 4b: mt19937 seed of the C++ MC (default 5489)
CACHE = HERE / "cache" / ("e0_cpp" + (f"_s{SEED}" if SEED else ""))
IMG = ICE / "scripts" / "findoptimal_robustness" / "cache" / "images"
MIC = EX / "SimInput" / "rand_500grains_1mm_neworient_s0.mic"
TEMPLATE = ICE / "benchmarks" / "phase_d_seed_diag" / "level0" / "e0" / "cpp_template.config"
H = W = 2048
N_DET = 2


def scratch() -> Path:
    return Path(os.environ["E0_CPP_SCRATCH"])


def write_case(vidx: int, variant: str, root: Path) -> Path:
    """Frames, one-voxel mic and config for a case; returns the config path."""
    keys = np.load(IMG / f"v{vidx}_{variant}.npz")["keys"].astype(np.int64)
    hw = H * W
    img, pos = keys // hw, keys % hw
    frame, det = img // N_DET, img % N_DET
    row, col = pos // W, pos % W
    d = root / "frames"
    d.mkdir(parents=True, exist_ok=True)
    for f in range(180):
        for k in range(N_DET):
            m = (frame == f) & (det == k)
            with open(d / f"500Grains.sim{f:05d}.d{k}", "w") as fh:
                for c, r in zip(col[m], row[m]):
                    fh.write(f"{c}, {r}, 1000.0\n")
    lines = MIC.read_text().splitlines()
    vox = lines[vidx + 1].split()
    vox[6:9] = ["0.000000"] * 3
    (root / "one.mic").write_text(lines[0] + "\n" + "\t".join(vox) + "\n")
    cfg = TEMPLATE.read_text()
    cfg = cfg.replace("@FRAMES@", f"{d}/500Grains.sim").replace("@MIC@", f"{root}/one.mic")
    cfg = cfg.replace("@OUT@", f"{root}/out")
    (root / "run.config").write_text(cfg)
    return root / "run.config"


def parse_log(text: str) -> Dict[str, Any]:
    """CD/CQ/CK/CF/CR/CV lines of the adaptive part of the log."""
    out: Dict[str, Any] = dict(CD={}, CQ={}, CK={}, CF=[], CR=None, CV=None)
    started = False
    for ln in text.splitlines():
        if ln.startswith("Normal Reconstruction"):
            started = True
        if not started:
            continue
        t = ln.split()
        if not t:
            continue
        if t[0] in ("CD", "CQ"):
            lv = int(t[1])
            out[t[0]].setdefault(lv, []).append([float(x) for x in t[3:]])  # cost + 9
        elif t[0] == "CK":
            out["CK"][int(t[1])] = int(t[2])
        elif t[0] == "CF":
            i = t.index("in")
            out["CF"].append([float(x) for x in t[2:i]] + [float(x) for x in t[i + 1 :]])
        elif t[0] == "CR":
            out["CR"] = [float(x) for x in t[1:2]] + [float(x) for x in t[4:13]]
            out["converged"] = t[3] == "1"
        elif t[0] == "CV":
            out["CV"] = [float(x) for x in t[1:11]]
    return out


def run_case(item: Tuple[int, str]) -> str:
    vidx, variant = item
    res = CACHE / f"v{vidx}_{variant}.json"
    if res.exists():
        return f"skip {vidx} {variant}"
    root = scratch() / f"case_{vidx}_{variant}"
    t0 = time.time()
    cfg = write_case(vidx, variant, root)
    env = dict(os.environ, CTRACE="1")
    p = subprocess.run(
        [os.environ["E0_CPP_BIN"], "r", str(cfg)],
        cwd=EX, env=env, capture_output=True, text=True, check=False,
    )  # fmt: skip
    d = parse_log(p.stdout)
    d.update(voxel=vidx, variant=variant, wall_s=time.time() - t0, returncode=p.returncode)
    m = re.search(r"adap_sec=([\d.]+) adap_evals=(\d+)", p.stdout)
    if m:
        d["adap_sec"], d["adap_evals"] = float(m.group(1)), int(m.group(2))
    CACHE.mkdir(parents=True, exist_ok=True)
    res.write_text(json.dumps(d))
    # frames are 360 small files per case: drop them after the run
    import shutil

    shutil.rmtree(root, ignore_errors=True)
    return f"{vidx} {variant} {d['wall_s']:.0f}s rc={p.returncode} CV={'yes' if d['CV'] else 'NO'}"


def cases(only: Optional[str], limit: int) -> List[Tuple[int, str]]:
    import common as C

    info = dict(np.load(C.OUT_DIR / "voxels.npz"))
    vox = [int(v) for v in info["voxel_indices"]][: C.N_VOXELS]
    out = [(v, var) for v in vox for var in C.VARIANTS
           if (C.CACHE_DIR / "e0" / f"v{v}_{var}.npz").exists()]  # fmt: skip
    if only:
        keep = {tuple(x.split(":")) for x in Path(only).read_text().split()}
        out = [(v, var) for v, var in out if (str(v), var) in keep]
    return out[:limit] if limit else out


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    ap.add_argument("cmd", choices=["run"])
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--only-file", default=None, help="file of voxel:variant entries")
    ap.add_argument("--limit", type=int, default=0)
    a = ap.parse_args()
    items = [
        c for c in cases(a.only_file, a.limit) if not (CACHE / f"v{c[0]}_{c[1]}.json").exists()
    ]
    print(f"{len(items)} cases, {a.workers} workers", flush=True)
    t0 = time.time()
    with mp.get_context("spawn").Pool(a.workers) as pool:
        for k, r in enumerate(pool.imap_unordered(run_case, items), 1):
            print(f"  [{k}/{len(items)}] {r} ({time.time() - t0:.0f}s)", flush=True)


if __name__ == "__main__":
    main()
