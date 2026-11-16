"""fine1: the gradient check, a hard gate before the freeze.

  uv run python scripts/autoresearch/fine1/fine1_gradient_check.py

Reads the written truth back from the ``.mic`` (not from the builder's arrays) and, over all
neighbour pairs of one grain (two voxels sharing an edge), measures the misorientation with the
cubic symmetry operators: the count, the maximum, the median, the 90th, 99th and 99.9th
percentiles and a histogram; and per grain: the spread (largest misorientation between a voxel
and the grain's centre orientation), its drawn level, its voxels and its largest neighbour step.
Gate: the maximum is under 0.5 degree and the largest grain spread is at least 3 degrees.
Writes ``benchmarks/autoresearch/fine1/gradient_check.json`` and ``gradient_tables.md``; exits 1
if the gate fails.
"""

import json
import sys
from pathlib import Path
from typing import Any, Dict, List

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import fine1_lib as F  # noqa: E402
import fine1_paths as P  # noqa: E402


def main() -> int:
    from icenine.mic_file import MicFile
    from icenine.ray_indexing import misorientation_deg, proper_rotations
    from icenine.symmetry import create_cubic_symmetry

    mic = MicFile.read(str(P.TRUTH))
    R = np.array([v.orientation for v in mic.voxels], dtype=np.float64)
    xy = np.array([v.position[:2] for v in mic.voxels], dtype=np.float64)
    up = np.array([v.points_up for v in mic.voxels], dtype=bool)
    side = float(mic.voxels[0].side_length)
    assert len(R) == P.N_FINE and all(abs(v.side_length - side) < 1e-12 for v in mic.voxels)
    info = json.loads(P.SAMPLE_JSON.read_text())
    gmap = np.load(P.GRAINMAP)["grain"]
    assert info["truth_sha256"] == P.sha256(P.TRUTH), "truth file differs from sample.json"
    sym = proper_rotations(create_cubic_symmetry(3.6))
    centres = F.triangle_vertices(xy, up, side).mean(axis=1)
    pairs = F.edge_pairs(centres, side)
    same = gmap[pairs[:, 0]] == gmap[pairs[:, 1]]
    inner = pairs[same]
    ang = np.asarray(misorientation_deg(R[inner[:, 0]], R[inner[:, 1]], sym), dtype=float)
    # per grain
    grains: List[Dict[str, Any]] = []
    gid_of_pair = gmap[inner[:, 0]]
    capped = {c["id"] for c in json.loads(P.CAPS_JSON.read_text())["capped"]}
    for g in info["grains"]:
        gi = int(g["id"])
        sel = np.nonzero(gmap == gi)[0]
        Rc = np.broadcast_to(np.array(g["R_centre"]), (len(sel), 3, 3))
        spread = float(np.max(misorientation_deg(R[sel], Rc, sym)))
        a = ang[gid_of_pair == gi]
        grains.append(
            dict(
                id=gi,
                level_deg=float(g["level_deg"]),
                n_voxels=int(len(sel)),
                spread_deg=spread,
                n_pairs=int(len(a)),
                max_step_deg=float(a.max()) if len(a) else 0.0,
                capped=gi in capped,
            )
        )
    spreads = np.array([g["spread_deg"] for g in grains])
    gt = F.gate(len(ang), float(ang.max()), float(spreads.max()))
    big = [g for g in grains if g["n_voxels"] > 1]
    out: Dict[str, Any] = dict(
        gate=gt,
        n_edge_pairs_all=int(len(pairs)),
        n_pairs_between_grains=int((~same).sum()),
        n_pairs_inside_grains=int(len(ang)),
        pair_max_deg=float(ang.max()),
        pair_median_deg=float(np.median(ang)),
        pair_percentiles_deg=F.percentiles(ang),
        pair_histogram=F.histogram(ang),
        n_pairs_at_or_over_0p45=int((ang >= 0.45).sum()),
        spread_max_deg=float(spreads.max()),
        spread_median_deg=float(np.median(spreads)),
        spread_median_of_grains_over_level0_deg=float(
            np.median([g["spread_deg"] for g in grains if g["level_deg"] > 0 and g["n_voxels"] > 1])
        ),
        spread_by_level={
            str(L): [
                float(np.min([g["spread_deg"] for g in big if g["level_deg"] == L])),
                float(np.median([g["spread_deg"] for g in big if g["level_deg"] == L])),
                float(np.max([g["spread_deg"] for g in big if g["level_deg"] == L])),
            ]
            for L in P.LEVELS
            if any(g["level_deg"] == L for g in big)
        },
        n_grains=len(grains),
        n_grains_capped=int(sum(g["capped"] for g in grains)),
        truth_sha256=info["truth_sha256"],
        grains=grains,
    )
    P.GRADIENT_CHECK.write_text(json.dumps(out, indent=1) + "\n")
    # a plain table file: summary, histogram, grains
    pc = out["pair_percentiles_deg"]
    lines = [
        "Reading line: neighbour pairs are two voxels of one grain that share an edge; the "
        "misorientation is in degrees, with cubic symmetry, read back from the written truth.",
        "",
        "| quantity | value |",
        "|---|---|",
        f"| neighbour pairs inside grains | {out['n_pairs_inside_grains']} |",
        f"| maximum | {out['pair_max_deg']:.4f} |",
        f"| median | {out['pair_median_deg']:.4f} |",
        f"| 90th percentile | {pc['p90']:.4f} |",
        f"| 99th percentile | {pc['p99']:.4f} |",
        f"| 99.9th percentile | {pc['p99.9']:.4f} |",
        f"| largest grain spread | {out['spread_max_deg']:.3f} |",
        f"| median grain spread (all {len(grains)} grains) | {out['spread_median_deg']:.3f} |",
        f"| grains capped | {out['n_grains_capped']} |",
        f"| gate passed | {gt['passed']} |",
        "",
        "| bin (degree) | pairs |",
        "|---|---|",
    ]
    for b in out["pair_histogram"]:
        hi = "and over" if b["hi"] > 1e6 else f"to {b['hi']:g}"
        lines.append(f"| {b['lo']:g} {hi} | {b['count']} |")
    lines += [
        "",
        "| grain | level | voxels | spread | largest step | capped |",
        "|---|---|---|---|---|---|",
    ]
    for g in grains:
        lines.append(
            f"| {g['id']} | {g['level_deg']:g} | {g['n_voxels']} | {g['spread_deg']:.3f} | "
            f"{g['max_step_deg']:.4f} | {'yes' if g['capped'] else ''} |"
        )
    (P.BENCH / "gradient_tables.md").write_text("\n".join(lines) + "\n")
    print(json.dumps({k: v for k, v in out.items() if k not in ("grains", "pair_histogram")}))
    print("GATE", "PASSED" if gt["passed"] else "FAILED")
    return 0 if gt["passed"] else 1


if __name__ == "__main__":
    sys.exit(main())
