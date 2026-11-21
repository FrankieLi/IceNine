"""A few frames of ``stacks/clean.npy`` against the ASCII files C++ reads (cpptime).

  uv run python scripts/autoresearch/cpptime/stacks_check.py [--frames 0,17,90,179] [--out F]

Python's pilot runs read the memory-mapped stack; C++ reads the 360 ASCII files. For each
(omega, detector) frame the set of lit pixels must be the same. The ASCII line is
``<x>, <y>, <intensity>`` and is compared with ``stack[det * 180 + omega][y, x] > 0`` (and with
``[x, y]``: the one that matches is reported). Exit 0 if every frame matches, else 1.
"""

import argparse
import json
import sys
from pathlib import Path
from typing import Dict, Set, Tuple

import numpy as np

ICE = Path(__file__).resolve().parents[3]
DATA = ICE.parent / "Examples" / "Example2.ManyGrains" / "ScatteringData_PhaseD"
N_OMEGA = 180


def ascii_pixels(path: Path) -> Set[Tuple[int, int]]:
    """Set of (first, second) coordinates of the lit pixels of an ASCII frame file."""
    out: Set[Tuple[int, int]] = set()
    for ln in path.read_text().splitlines():
        t = ln.replace(",", " ").split()
        if len(t) >= 3 and float(t[2]) > 0:
            out.add((int(t[0]), int(t[1])))
    return out


def stack_pixels(frame: np.ndarray, swap: bool) -> Set[Tuple[int, int]]:
    """Lit pixels of a frame as (x, y) pairs; ``swap`` reads the frame as [x, y]."""
    a, b = np.nonzero(frame)
    pairs = zip(a.tolist(), b.tolist()) if swap else zip(b.tolist(), a.tolist())
    return set(pairs)


def check(frames: list, stack_path: Path, data: Path = DATA) -> Dict[str, object]:
    stack = np.load(str(stack_path), mmap_mode="r")
    res: Dict[str, object] = {}
    for f in frames:
        omega = int(f)
        for det in (0, 1):
            txt = data / "full" / "clean" / f"500Grains.sim{omega:05d}.d{det}"
            a = ascii_pixels(txt)
            fr = np.asarray(stack[det * N_OMEGA + omega])
            res[f"omega{omega}_det{det}"] = dict(
                n_ascii=len(a),
                n_stack=int((fr > 0).sum()),
                equal_xy=a == stack_pixels(fr, swap=False),
                equal_yx=a == stack_pixels(fr, swap=True),
            )
    return res


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--frames", default="0,17,90,179", help="omega indices, both detectors each")
    ap.add_argument("--stack", type=Path, default=DATA / "stacks" / "clean.npy")
    ap.add_argument("--out", type=Path, default=None)
    a = ap.parse_args()
    res = check([int(x) for x in a.frames.split(",")], a.stack)
    text = json.dumps(res, indent=1)
    if a.out:
        a.out.write_text(text + "\n")
    print(text)
    ok = all(bool(v["equal_xy"]) or bool(v["equal_yx"]) for v in res.values())  # type: ignore
    print("OK" if ok else "MISMATCH")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
