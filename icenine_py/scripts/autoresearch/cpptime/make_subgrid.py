"""Write a small connected sub-grid of the pilot region (cpptime short checks).

  uv run python scripts/autoresearch/cpptime/make_subgrid.py --n 30 --centre-row 0 \
      --out DIR/sub.mic --rows-out DIR/sub_rows.npy

Takes the ``n`` rows of ``region_grid.mic`` nearest to the centre row (stable sort by distance,
ties by row), keeps the header line and the row order of the pilot grid, and saves the row
indices (for ``score_cpp.py --rows``).
"""

import argparse
from pathlib import Path
from typing import Tuple

import numpy as np

ICE = Path(__file__).resolve().parents[3]
REGION_GRID = ICE / "benchmarks" / "phase_d_pilot" / "region_grid.mic"


def nearest_rows(xy: np.ndarray, centre: int, n: int) -> np.ndarray:
    """Sorted row indices of the ``n`` points nearest to point ``centre`` (itself included)."""
    d = np.hypot(xy[:, 0] - xy[centre, 0], xy[:, 1] - xy[centre, 1])
    return np.sort(np.argsort(d, kind="stable")[:n])


def make(grid: Path, centre: int, n: int) -> Tuple[str, np.ndarray]:
    lines = grid.read_text().splitlines()
    body = [ln for ln in lines[1:] if ln.strip()]
    xy = np.array([[float(t) for t in ln.split()[:2]] for ln in body])
    rows = nearest_rows(xy, centre, n)
    return "\n".join([lines[0]] + [body[i] for i in rows]) + "\n", rows


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--n", type=int, default=30)
    ap.add_argument("--centre-row", type=int, default=0)
    ap.add_argument("--grid", type=Path, default=REGION_GRID)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--rows-out", type=Path, required=True)
    a = ap.parse_args()
    text, rows = make(a.grid, a.centre_row, a.n)
    a.out.parent.mkdir(parents=True, exist_ok=True)
    a.out.write_text(text)
    np.save(a.rows_out, rows)
    print(f"wrote {a.out} ({len(rows)} voxels)")


if __name__ == "__main__":
    main()
