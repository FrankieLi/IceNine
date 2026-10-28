"""S2 Steps 1-2 review fix 4: what ``-ffp-contract=off`` changes in the key-off C code.

On clang arm64 the flag leaves ``_fused_eval.c``'s object code unchanged (its multiply-adds are
explicit ``fmaf`` calls under a no-contraction pragma) but removes clang's default in-expression
contraction from ``_rasterize.c``: the multiply-adds of the Sutherland-Hodgman intersections
(``p0.y + slope * (b - p0.x)``). This script compares ``triangle_overlap`` of two builds of the
``_rasterize`` module on random triangles, most of which cross the image border (where the clip
computes intersections), and classifies each differing triangle against the Python reference
(``ImageData._sutherland_hodgman_clip`` + ``_scanline_fill``, the hard-mode path of
``get_triangle_overlap_property`` without C).

Usage (from ``icenine_py/``; the two ``.so`` files in separate directories)::

    uv run python scripts/speed_parity/s2/contract_fuzz.py OLD.so NEW.so --out OUT.json

OLD: built without the flag (for example at commit 37c1e47, or with ``extra_compile_args``
removed from ``setup.py``); NEW: the current build.
"""

import argparse
import importlib.machinery
import importlib.util
import json
from pathlib import Path
from types import ModuleType
from typing import Any, Dict, List, Tuple

import numpy as np

from icenine.image_data import ImageData

SEED = 20261027
ROWS, COLS = 64, 80
N = 200_000
# (label, low, high, half-size of the vertex offsets in px); "integer" floors to the same vertices
ARMS = [
    ("small_border", -6.0, 86.0, 6.0),
    ("large_border", -40.0, 120.0, 60.0),
    ("integer_border", -6.0, 86.0, 6.0),
]


def load(path: str, name: str) -> ModuleType:
    loader = importlib.machinery.ExtensionFileLoader(name, path)
    spec = importlib.util.spec_from_file_location(name, path, loader=loader)
    assert spec is not None
    mod = importlib.util.module_from_spec(spec)
    loader.exec_module(mod)
    return mod


def python_reference(img: np.ndarray, tri: np.ndarray) -> Tuple[int, int]:
    """(overlap, lit) of the Python reference, from the double coordinates."""
    rows, cols = img.shape

    def trunc(v: float) -> float:
        return -1.0 if v < 0 else float(int(v))

    poly = [(trunc(x), trunc(y)) for x, y in tri]
    clipped = ImageData._sutherland_hodgman_clip(poly, 0.0, float(cols - 1), 0.0, float(rows - 1))
    if len(clipped) < 3:
        return (0, 0)
    pixels = ImageData._scanline_fill([(round(x), round(y)) for x, y in clipped])
    lit = overlap = 0
    for c, r in pixels or []:
        if 0 <= c < cols and 0 <= r < rows:
            lit += 1
            overlap += int(img[r, c] > 0)
    return (overlap, lit)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("old_so")
    ap.add_argument("new_so")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    old = load(args.old_so, "old._rasterize")
    new = load(args.new_so, "new._rasterize")
    rng = np.random.default_rng(SEED)
    results: Dict[str, Any] = {
        "meta": dict(seed=SEED, rows=ROWS, cols=COLS, n=N, image="uniform, 30% bright")
    }
    for label, lo, hi, scale in ARMS:
        img = (rng.random((ROWS, COLS)) < 0.3).astype(np.float32)
        base = rng.uniform(lo, hi, size=(N, 2))
        tris = base[:, None, :] + rng.uniform(-scale, scale, size=(N, 3, 2))
        if label == "integer_border":
            tris = np.round(tris) + 0.5
        crossing = differ = 0
        agrees = {"old": 0, "new": 0, "neither": 0}
        examples: List[Dict[str, Any]] = []
        for t in tris:
            flat = [float(v) for v in t.reshape(-1)]
            a = tuple(old.triangle_overlap(img, *flat))
            b = tuple(new.triangle_overlap(img, *flat))
            x, y = t[:, 0], t[:, 1]
            if x.min() < 0 or y.min() < 0 or x.max() >= COLS - 1 or y.max() >= ROWS - 1:
                crossing += 1
            if a != b:
                differ += 1
                ref = python_reference(img, t)
                agrees["old" if a == ref else "new" if b == ref else "neither"] += 1
                if len(examples) < 5:
                    examples.append(dict(tri=t.tolist(), old=a, new=b, python_reference=ref))
        results[label] = dict(
            n=N,
            crossing_border=crossing,
            differ=differ,
            python_reference_agrees_with=agrees,
            examples=examples,
        )
        print(f"{label:15s} n={N} crossing {crossing} differ {differ} reference agrees {agrees}")
    Path(args.out).write_text(json.dumps(results, indent=1) + "\n")


if __name__ == "__main__":
    main()
