"""Task A step 0: write the fixed orientation lists for the C++ parity driver.

Each record is (kind, id, centre 3x3 float32, local 3x3 float32); the driver evaluates the product
local @ centre (float32, as C++ does) with three arms. Kind 0 = a single orientation (local =
identity); kind 1 = an (FZ centre, local-grid matrix) pair of the level-0 screen.

Per voxel of 15901, 22867, 5242:
  truth, and 0.1 / 1 / 3 deg off the truth (3 random axes each)     -> 10 records
  random orientations (uniform over SO(3); 200 over the three voxels) -> 67 / 67 / 66
  random level-0 pairs (FZ point, local-grid index j), 300           -> kind 1
  offset sweep (random centre, local rotation of 0.05 .. 5 deg), 8 x 25 -> kind 1
For 15901 only (the stored level-0 trace has per-FZ-point pass sets): 500 pairs where the stored
C++ and Python (C++ rule, centroid) trace disagree (pass status first, then screen overlap count).

  uv run python scripts/eval_reduction/parity_inputs.py --out DIR
Writes DIR/in_<voxel>.bin and DIR/in_<voxel>.json (labels).
"""

import argparse
import json
import struct
import sys
from pathlib import Path
from typing import Any, Dict, List

import numpy as np
from scipy.spatial.transform import Rotation

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import parity_common as PC  # noqa: E402

TRACE = PC.ICE / "benchmarks" / "phase_d_seed_diag" / "level0"
RNG_SEED = 12345  # fixes the input list only (not an experiment seed)
SWEEP_DEG = [0.05, 0.1, 0.2, 0.5, 1.0, 2.0, 3.0, 5.0]


def rot_about(axis: np.ndarray, deg: float) -> np.ndarray:
    return Rotation.from_rotvec(axis / np.linalg.norm(axis) * np.radians(deg)).as_matrix()


def differing_pairs(n_want: int) -> List[Dict[str, Any]]:
    cpp = json.load(open(TRACE / "cpp_v15901.json"))["fz"]
    py = json.load(open(TRACE / "py_cpprule_centroid_v15901.json"))["fz"]
    pass_diff: List[Dict[str, Any]] = []
    cnt_diff: List[Dict[str, Any]] = []
    for c, p in zip(cpp, py):
        assert c["fz"] == p["fz"]
        for j, (a, b) in enumerate(zip(c["peak_overlap"], p["peak_overlap"])):
            if (a > 0) != (b > 0):
                pass_diff.append(dict(fz=c["fz"], j=j, cpp=a, py=b, why="pass"))
            elif a != b:
                cnt_diff.append(dict(fz=c["fz"], j=j, cpp=a, py=b, why="count"))
    rng = np.random.default_rng(RNG_SEED)
    out = pass_diff[:n_want]
    if len(out) < n_want:
        k = n_want - len(out)
        pick = rng.permutation(len(cnt_diff))[:k]
        out += [cnt_diff[i] for i in sorted(pick)]
    print(
        f"differing pairs available: pass {len(pass_diff)}, count-only {len(cnt_diff)};"
        f" using {len(out)}"
    )
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("--n-random", type=int, default=200)
    ap.add_argument("--n-pairs-random", type=int, default=300)
    ap.add_argument("--n-pairs-differ", type=int, default=500)
    ap.add_argument("--n-sweep", type=int, default=25)
    a = ap.parse_args()
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    ctx = PC.PyContext()
    fz = np.asarray(ctx.setup.fz_orientations, dtype=np.float64)
    lg = ctx.local_grid()
    eye = np.eye(3)
    rng = np.random.default_rng(RNG_SEED)
    rng_sweep = np.random.default_rng(RNG_SEED + 1)
    n_vox = len(PC.PARITY_VOXELS)
    for vi, v in enumerate(PC.PARITY_VOXELS):
        recs: List[Dict[str, Any]] = []

        def add(kind: int, label: str, centre: np.ndarray, local: np.ndarray, **meta: Any) -> None:
            recs.append(
                dict(
                    id=len(recs),
                    kind=kind,
                    label=label,
                    centre=f32(centre),
                    local=f32(local),
                    **meta,
                )
            )

        def f32(x: np.ndarray) -> np.ndarray:
            return np.asarray(x, dtype=np.float32)

        R = ctx.truth(v)
        add(0, "truth", R, eye)
        for deg in (0.1, 1.0, 3.0):
            for k in range(3):
                ax = rng.normal(size=3)
                add(0, f"off_{deg}", rot_about(ax, deg) @ R, eye, deg=deg)
        n = a.n_random // n_vox + (1 if vi < a.n_random % n_vox else 0)
        for k in range(n):
            add(
                0,
                "random",
                Rotation.random(random_state=int(rng.integers(1 << 30))).as_matrix(),
                eye,
            )
        for k in range(a.n_pairs_random):
            i, j = int(rng.integers(len(fz))), int(rng.integers(len(lg)))
            add(1, "pair_random", fz[i], lg[j], fz=i, j=j)
        # offset sweep for the interpolation: a random orientation and a local rotation of theta
        # degrees about a random axis (separate generator, so the lists above do not change)
        for theta in SWEEP_DEG:
            for k in range(a.n_sweep):
                c = Rotation.random(random_state=int(rng_sweep.integers(1 << 30))).as_matrix()
                add(1, "sweep", c, rot_about(rng_sweep.normal(size=3), theta), theta=theta)
        if v == 15901:
            for d in differing_pairs(a.n_pairs_differ):
                add(1, "pair_differ_" + d["why"], fz[d["fz"]], lg[d["j"]], **d)
        with open(out / f"in_{v}.bin", "wb") as fh:
            fh.write(struct.pack("<i", len(recs)))
            for r in recs:
                fh.write(struct.pack("<ii", r["kind"], r["id"]))
                fh.write(r["centre"].astype("<f4").tobytes())
                fh.write(r["local"].astype("<f4").tobytes())
        meta = [{k: x for k, x in r.items() if k not in ("centre", "local")} for r in recs]
        for r, m in zip(recs, meta):
            m["centre"] = r["centre"].tolist()
            m["local"] = r["local"].tolist()
        (out / f"in_{v}.json").write_text(json.dumps(meta))
        print(f"voxel {v}: {len(recs)} records")


if __name__ == "__main__":
    main()
