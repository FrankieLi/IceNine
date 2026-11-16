"""fine1: pack the clean frames of ``fine7`` into the uint8 stack and check it equals the frames.

  uv run python scripts/autoresearch/fine1/fine1_stack.py [--build]

The stack is built as for ``held`` (``held_stack_note.py``): ``write_binary_stack(<frames>,
'500Grains.sim', 'd', 5, 180, 2, 2048, 2048, <stack>)``, a 0/1 mask, frame ``det * 180 + omega``,
pixel ``[k, j]`` for the ASCII line ``j, k, intensity``. The check, for all 360 frames: the set of
lit pixels of the ASCII file (what C++ reads) equals the set of non-zero pixels of the stack (what
Python reads), and the stack's total equals the render record's ``lit_clean_total``. Writes
``benchmarks/autoresearch/fine1/stack_check.json``; exits 1 on any difference.
"""

import argparse
import json
import sys
from pathlib import Path
from typing import Any, Dict

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import fine1_paths as P  # noqa: E402

RECORD = P.ICE / "scripts" / "phase_d" / "results" / "render_fine7_nodrop.json"


def ascii_lit(path: Path) -> np.ndarray:
    """Lit pixels (n, 2) of an ASCII frame as (k, j) = (row, column), intensity above zero."""
    a = np.loadtxt(str(path), delimiter=",", usecols=(0, 1, 2), ndmin=2)
    if a.size == 0:
        return np.zeros((0, 2), dtype=np.int64)
    a = a[a[:, 2] > 0]
    return np.stack([a[:, 1], a[:, 0]], axis=1).astype(np.int64)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--build", action="store_true", help="build the stack if it is missing")
    a = ap.parse_args()
    built = False
    if not P.STACK.exists():
        if not a.build:
            raise SystemExit(f"{P.STACK} is missing; run with --build")
        from icenine.experimental_data import write_binary_stack

        write_binary_stack(P.FRAMES, "500Grains.sim", "d", 5, 180, 2, 2048, 2048, P.STACK)
        built = True
    st = np.load(str(P.STACK), mmap_mode="r")
    assert st.shape == (360, 2048, 2048), st.shape
    bad = []
    total_ascii = 0
    total_stack = 0
    for det in (0, 1):
        for om in range(180):
            fr = np.asarray(st[det * 180 + om])
            lit = ascii_lit(P.FRAMES / f"500Grains.sim{om:05d}.d{det}")
            nz = np.argwhere(fr > 0)
            total_ascii += len(lit)
            total_stack += len(nz)
            same = len(lit) == len(nz) and np.array_equal(
                np.unique(lit[:, 0] * 2048 + lit[:, 1]), nz[:, 0] * 2048 + nz[:, 1]
            )
            if not same:
                bad.append(dict(det=det, omega=om, n_ascii=len(lit), n_stack=len(nz)))
    rec = json.loads(RECORD.read_text())
    out: Dict[str, Any] = dict(
        stack=str(P.STACK.relative_to(P.BENCH)),
        built_now=built,
        stack_shape=list(st.shape),
        n_frames=360,
        n_frames_different=len(bad),
        different=bad[:10],
        lit_pixels_ascii=int(total_ascii),
        lit_pixels_stack=int(total_stack),
        render_lit_clean_total=int(rec["lit_clean_total"]),
        render_spots_total=int(rec["spots_total"]),
        render_wall_s=float(rec["wall_s"]),
        render_workers=int(rec["workers"]),
        render_timing_label=rec["timing_label"],
        stack_sha256=P.sha256(P.STACK),
    )
    out["ok"] = bool(
        not bad
        and out["lit_pixels_ascii"] == out["lit_pixels_stack"] == out["render_lit_clean_total"]
    )
    P.STACK_JSON.write_text(json.dumps(out, indent=1) + "\n")
    print(json.dumps(out, indent=1))
    return 0 if out["ok"] else 1


if __name__ == "__main__":
    sys.exit(main())
