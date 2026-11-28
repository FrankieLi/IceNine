"""fine1: the generation-7 grid of the fine layer, and its slice.

  uv run python scripts/autoresearch/fine1/fine1_make_grid.py

Splits every triangle of ``SimInput/rand_500grains_1mm_neworient_s0_grid.mic`` (generation 6,
24,570 voxels) into four, the way C++ does (``fine1_lib.split_once``), and writes
``<cache>/grid/fine7_grid.mic`` in the text format of the source grid (tab separated, ten
columns, orientation zero): 98,280 voxels of generation 7, side 0.0046875 mm. Then the checks:
the count, the area, every child inside its parent, every child row equal (within 1e-6 mm) to the
rows of a literal port of C++ ``SubdivideTriangles`` for all 24,570 parents, the file read back by
the Python ``MicFile``, and C++'s generation rule (``SetMinimumResolution``, side stored as a
32-bit float) finding generation 7 so that nothing is split again for MinSideLength 0.0046875.
Writes the slice grid (the 2,000 voxels nearest the layer centre by the left vertex, as
``render_full.py --n-voxels``) and the tracked records ``grid_check.json`` and ``slice.json``.
Exits 1 if a check fails.
"""

import json
import sys
from pathlib import Path
from typing import Any, Dict

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import fine1_lib as F  # noqa: E402
import fine1_paths as P  # noqa: E402


def read_rows(path: Path) -> Dict[str, Any]:
    lines = path.read_text().splitlines()
    body = [ln for ln in lines[1:] if ln.strip()]
    tok = [ln.split() for ln in body]
    return dict(
        header=lines[0].strip(),
        side0=float(lines[0]),
        x=np.array([float(t[0]) for t in tok]),
        y=np.array([float(t[1]) for t in tok]),
        z=np.array([float(t[2]) for t in tok]),
        direction=np.array([int(t[3]) for t in tok]),
        gen=np.array([int(t[4]) for t in tok]),
        phase=np.array([int(t[5]) for t in tok]),
        conf=np.array([float(t[9]) if len(t) > 9 else 1.0 for t in tok]),
        n_tokens=sorted({len(t) for t in tok}),
    )


def write_rows(
    path: Path,
    header: str,
    xy: np.ndarray,
    z: np.ndarray,
    up: np.ndarray,
    gen: np.ndarray,
    phase: np.ndarray,
    conf: np.ndarray,
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    out = [header]
    for i in range(len(xy)):
        out.append(
            f"{xy[i, 0]:.9f}\t{xy[i, 1]:.9f}\t{z[i]:.6f}\t{F.UP if up[i] else F.DOWN}\t"
            f"{int(gen[i])}\t{int(phase[i])}\t0.000000\t0.000000\t0.000000\t{conf[i]:.6f}"
        )
    path.write_text("\n".join(out) + "\n")


def cpp_gen(init_side: float, min_side: float) -> int:
    """C++ ``SetMinimumResolution``: ceil(log(init / min) / log 2) with 32-bit floats."""
    r = np.float32(np.float32(init_side) / np.float32(min_side))
    return int(np.ceil(np.log(np.float64(r)) / np.log(2.0)))


def main() -> int:
    src = read_rows(P.SRC_GRID)
    n = len(src["x"])
    side0 = src["side0"]
    gens = set(src["gen"].tolist())
    assert gens == {6} and n == 24570, (gens, n)
    pside = side0 / 2.0**6
    pxy = np.stack([src["x"], src["y"]], axis=1)
    pup = src["direction"] == F.UP
    cxy, cup, cside = F.split_once(pxy, pup, pside)
    rep = np.repeat(np.arange(n), 4)
    write_rows(
        P.GRID,
        src["header"],
        cxy,
        src["z"][rep],
        cup,
        np.full(4 * n, 7),
        src["phase"][rep],
        src["conf"][rep],
    )
    chk: Dict[str, Any] = dict(F.grid_checks(pxy, pup, pside, cxy, cup, cside))
    # the file as written, read back by the Python reader
    from icenine.mic_file import MicFile

    back = MicFile.read(str(P.GRID))
    bxy = np.array([v.position[:2] for v in back.voxels])
    bup = np.array([v.points_up for v in back.voxels])
    chk["python_reader_n"] = len(back.voxels)
    chk["python_reader_sides_mm"] = sorted({round(float(v.side_length), 12) for v in back.voxels})
    chk["python_reader_all_generation_7"] = bool(all(v.generation == 7 for v in back.voxels))
    chk["file_vs_arrays_max_abs_mm"] = float(np.abs(bxy - cxy).max())
    chk["file_direction_codes_equal"] = bool(np.array_equal(bup, cup))
    # a literal port of C++ SubdivideTriangles for every parent
    worst = 0.0
    bad_dir = 0
    for i in range(n):
        kids = F.cpp_subdivide(pxy[i, 0], pxy[i, 1], bool(pup[i]), pside, 6, 7)
        assert len(kids) == 4
        for k, (kx, ky, ku) in enumerate(kids):
            worst = max(worst, abs(kx - bxy[4 * i + k, 0]), abs(ky - bxy[4 * i + k, 1]))
            bad_dir += int(ku != bool(bup[4 * i + k]))
    chk["cpp_port_all_parents"] = dict(
        n_parents=n,
        max_abs_diff_mm=float(worst),
        n_direction_mismatch=bad_dir,
        tolerance_mm=1e-6,
    )
    chk["cpp_generation_rule"] = dict(
        min_side_length=P.SIDE_FINE_MM,
        generation_found_for_header_side=cpp_gen(side0, P.SIDE_FINE_MM),
        generation_of_rows=7,
        nothing_split_again=bool(cpp_gen(side0, P.SIDE_FINE_MM) == 7),
        generation_found_for_source_grid=cpp_gen(side0, 0.009375),
        float_type="32-bit (XDM++/libXDM/Types.h)",
    )
    chk["n_voxels_expected"] = P.N_FINE
    chk["sha256_source_grid"] = P.sha256(P.SRC_GRID)
    chk["sha256_grid"] = P.sha256(P.GRID)
    chk["ok"] = bool(
        chk["ok"]
        and chk["n_child"] == P.N_FINE
        and chk["python_reader_n"] == P.N_FINE
        and chk["python_reader_all_generation_7"]
        and chk["file_vs_arrays_max_abs_mm"] < 1e-6
        and chk["file_direction_codes_equal"]
        and worst < 1e-6
        and bad_dir == 0
        and chk["cpp_generation_rule"]["nothing_split_again"]
    )
    # the slice: the 2,000 fine voxels nearest the centre (left vertex), as render_full --n-voxels
    raw = np.loadtxt(P.GRID, skiprows=1)
    rows = F.nearest_centre_rows(raw[:, :2], P.SLICE_N)
    write_rows(
        P.SLICE,
        src["header"],
        cxy[rows],
        src["z"][rep][rows],
        cup[rows],
        np.full(len(rows), 7),
        src["phase"][rep][rows],
        src["conf"][rep][rows],
    )
    sl: Dict[str, Any] = dict(
        n=int(len(rows)),
        rule="np.sort(np.argsort(|left vertex (x, y)|)[:2000]) over the rows of fine7_grid.mic",
        rows=[int(r) for r in rows],
        sha256_slice_grid=P.sha256(P.SLICE),
    )
    P.BENCH.mkdir(parents=True, exist_ok=True)
    P.GRID_CHECK.write_text(json.dumps(chk, indent=1) + "\n")
    P.SLICE_JSON.write_text(json.dumps(sl) + "\n")
    print(json.dumps({k: v for k, v in chk.items() if k != "python_reader_sides_mm"}, indent=1))
    print("slice voxels", len(rows), "grid ok:", chk["ok"])
    return 0 if chk["ok"] else 1


if __name__ == "__main__":
    sys.exit(main())
