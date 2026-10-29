"""Task G2: helpers of the BFS pilot runs (scripts/phase_d/pilot_bfs.py does the runs).

subregion  write a sub-region of the pilot region grid (the N voxels nearest the region's
           centre, stable sort, ties by file index) for short pilots and the key-off check
compare    key-off check: two run folders of pilot_bfs.py are bit-identical (recon.mic,
           records.json without the wall times, run.json stats without the wall times)

uv run python scripts/eval_reduction/g2_bfs.py subregion --n 150 --out <grid.mic>
uv run python scripts/eval_reduction/g2_bfs.py compare <run dir A> <run dir B> [--json OUT]
"""

import argparse
import json
from pathlib import Path
from typing import Any, Dict, List

import numpy as np

ICE = Path(__file__).resolve().parents[2]
REGION = ICE / "benchmarks" / "phase_d_pilot"
CENTRE_XY = (0.17, 0.0)  # region.json centre_xy_mm


def subregion(n: int, out: Path) -> None:
    lines = (REGION / "region_grid.mic").read_text().splitlines()
    head, rows = lines[0], [r for r in lines[1:] if r.strip()]
    xy = np.array([[float(t) for t in r.split()[:2]] for r in rows])
    # triangle centroid: the first two columns are a vertex; the region was selected by centroid,
    # but for a compact patch the vertex distance is enough (the patch is not compared to it)
    d = np.hypot(xy[:, 0] - CENTRE_XY[0], xy[:, 1] - CENTRE_XY[1])
    keep = np.sort(np.argsort(d, kind="stable")[:n])
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text("\n".join([head] + [rows[i] for i in keep]) + "\n")
    print(f"wrote {out}: {len(keep)} of {len(rows)} voxels")


def _records(p: Path) -> Dict[str, Any]:
    rec = json.loads((p / "records.json").read_text())
    for r in rec.values():
        r.pop("wall_s", None)
    return rec


def _stats(p: Path) -> Dict[str, Any]:
    st = json.loads((p / "run.json").read_text())["stats"]
    out = {k: v for k, v in st.items() if "wall" not in k and k != "reseeds"}
    out["reseeds"] = [
        {k: v for k, v in e.items() if "seconds" not in k} for e in st.get("reseeds", [])
    ]
    return out


def compare(a: Path, b: Path) -> Dict[str, Any]:
    same_mic = (a / "recon.mic").read_bytes() == (b / "recon.mic").read_bytes()
    ra, rb = _records(a), _records(b)
    sa, sb = _stats(a), _stats(b)
    diff_rec: List[str] = [k for k in sorted(set(ra) | set(rb)) if ra.get(k) != rb.get(k)]
    diff_st = [k for k in sorted(set(sa) | set(sb)) if sa.get(k) != sb.get(k)]
    res = dict(
        a=str(a.name),
        b=str(b.name),
        recon_mic_identical=bool(same_mic),
        n_records=len(ra),
        records_differing=diff_rec,
        stats_differing=diff_st,
        n_voxels=sa.get("n_voxels"),
        n_seeds=sa.get("n_seeds"),
        n_reseeds=sa.get("n_reseeds"),
        n_neighbor_fits=sa.get("n_neighbor_fits"),
        n_revisit_attempted=sa.get("n_revisit_attempted"),
        identical=bool(same_mic and not diff_rec and not diff_st),
    )
    return res


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd", choices=["subregion", "compare"])
    ap.add_argument("dirs", nargs="*")
    ap.add_argument("--n", type=int, default=150)
    ap.add_argument("--out", default=None)
    ap.add_argument("--json", default=None)
    a = ap.parse_args()
    if a.cmd == "subregion":
        subregion(a.n, Path(a.out))
    else:
        res = compare(Path(a.dirs[0]), Path(a.dirs[1]))
        print(json.dumps(res, indent=1))
        if a.json:
            Path(a.json).write_text(json.dumps(res, indent=1) + "\n")


if __name__ == "__main__":
    main()
