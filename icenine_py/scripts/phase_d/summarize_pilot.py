"""Collect the scored Phase D pilot runs into benchmarks/phase_d_pilot/{summary.json,tables.md}.

Usage (from icenine_py/): uv run python scripts/phase_d/summarize_pilot.py --runs DIR
(DIR holds <arm>_<variant>/score.json and voxels.npz written by score_bfs.py).
Timings are contended (several concurrent single-process runs). P-values are grain-clustered.
"""

import argparse
import json
import sys
from pathlib import Path
from typing import Any, Dict, List

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(ROOT / "scripts" / "common"))

import score_bfs as S  # noqa: E402
from doc_tables import markdown_table, write_tables  # noqa: E402
from stats import wilson  # noqa: E402

RUNS = ["mc_clean", "cma_clean", "cma_noretry_clean", "mc_realistic_q16", "cma_realistic_q16"]
PAIRS = [
    ("mc_clean", "cma_clean"),
    ("cma_clean", "cma_noretry_clean"),
    ("mc_realistic_q16", "cma_realistic_q16"),
]


def sample_counts() -> Dict[str, int]:
    """Voxels, grains and connected grain pieces of the full sample, read from the sample files.

    Pieces are counted under the BFS neighbour radius (2 sides between left-vertex positions, as
    `MicFile.get_neighbors` in the BFS), the count the seed-count oracle used (497, MH "Phase D
    seed-cost diagnosis")."""
    from scipy.sparse import coo_matrix
    from scipy.sparse.csgraph import connected_components
    from scipy.spatial import cKDTree

    ex = ROOT.parent / "Examples" / "Example2.ManyGrains" / "SimInput"
    grain = np.load(ex / "rand_500grains_1mm_neworient_s0_grainmap.npy")
    lines = (ex / "rand_500grains_1mm_neworient_s0_grid.mic").read_text().splitlines()
    pos = np.array([[float(t) for t in ln.split()[:2]] for ln in lines[1:] if ln.strip()])
    side = json.loads((ROOT / "benchmarks" / "phase_d_pilot" / "region.json").read_text())[
        "side_length_mm"
    ]
    pr = cKDTree(pos).query_pairs(2.0 * side, output_type="ndarray")
    pr = pr[grain[pr[:, 0]] == grain[pr[:, 1]]]
    n = len(grain)
    comp = coo_matrix((np.ones(len(pr)), (pr[:, 0], pr[:, 1])), shape=(n, n))
    return {
        "n_voxels": int(n),
        "n_grains": int(len(np.unique(grain))),
        "n_pieces_bfs_radius": int(connected_components(comp, directed=False)[0]),
    }


def pct(x: Any) -> Any:
    return None if x is None else 100.0 * x


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--runs", required=True)
    ap.add_argument("--out", default=str(ROOT / "benchmarks" / "phase_d_pilot"))
    a = ap.parse_args()
    runs = Path(a.runs)
    sc: Dict[str, Any] = {}
    vox: Dict[str, Any] = {}
    for r in RUNS:
        if (runs / r / "score.json").exists():
            sc[r] = json.loads((runs / r / "score.json").read_text())
            vox[r] = np.load(runs / r / "voxels.npz")
    pairs = {}
    for x, y in PAIRS:
        if x in vox and y in vox:
            pairs[f"{x} vs {y}"] = S.clustered_paired_p(
                vox[x]["wrong"], vox[y]["wrong"], vox[x]["grain"]
            )
    # accepted = not left unresolved; "accepted wrong" = the BFS accepted a wrong orientation
    acc_wrong = {r: vox[r]["wrong"] & (vox[r]["source"] != "unresolved") for r in vox}
    pairs_acc = {}
    for x, y in PAIRS:
        if x in vox and y in vox:
            pairs_acc[f"{x} vs {y}"] = S.clustered_paired_p(
                acc_wrong[x], acc_wrong[y], vox[x]["grain"]
            )
    diag: Dict[str, Any] = {}
    cnt = sample_counts()
    n_full, n_grain_full, n_pieces, n_reg = (
        cnt["n_voxels"],
        cnt["n_grains"],
        cnt["n_pieces_bfs_radius"],
        2000,
    )
    full_grain = np.load(
        ROOT.parent
        / "Examples"
        / "Example2.ManyGrains"
        / "SimInput"
        / "rand_500grains_1mm_neworient_s0_grainmap.npy"
    )
    full_size = np.bincount(full_grain)
    for r, v in vox.items():
        g, src, wrong = v["grain"], v["source"], v["wrong"]
        gs = np.unique(g)
        seeded = set(g[src == "seed"].tolist())
        rejected_seed = set(g[v["seed_rejected"]].tolist())
        lost = [x for x in gs if not (~wrong[g == x]).any()]
        found = [x for x in gs if x not in set(lost)]
        size = {x: int((g == x).sum()) for x in gs}
        whole = [x for x in gs if size[x] == full_size[x]]
        trunc = [x for x in gs if size[x] < full_size[x]]
        bnd = v["boundary"]
        lo_, hi_ = wilson(len(lost), len(gs))
        unres = src == "unresolved"
        in_lost = np.isin(g, lost)
        c = sc[r]["cost"]
        nb_rev_h = (c["wall_neighbor_s"] + c["wall_revisit_s"]) / 3600.0
        seed_h = c["wall_seed_s"] / 3600.0
        diag[r] = {
            "n_grains": int(len(gs)),
            "grains_with_an_accepted_seed": len(seeded),
            "lost_grains": len(lost),
            "lost_grains_with_an_accepted_seed": len(set(lost) & seeded),
            "lost_grains_with_a_rejected_seed": len(set(lost) & rejected_seed),
            "lost_grains_no_full_search": len(set(lost) - rejected_seed - seeded),
            "lost_grains_rate": len(lost) / len(gs),
            "lost_grains_wilson_lo": lo_,
            "lost_grains_wilson_hi": hi_,
            "lost_median_region_voxels": float(np.median([size[x] for x in lost])),
            "found_median_region_voxels": float(np.median([size[x] for x in found])),
            "lost_boundary_voxel_share": float(bnd[np.isin(g, lost)].mean()),
            "found_boundary_voxel_share": float(bnd[~np.isin(g, lost)].mean()),
            "grains_wholly_in_region": len(whole),
            "wholly_in_region_lost": len(set(whole) & set(lost)),
            "grains_truncated_by_region": len(trunc),
            "truncated_lost": len(set(trunc) & set(lost)),
            "accepted_wrong_rate": float(acc_wrong[r].sum() / max((~unres).sum(), 1)),
            "voxels_in_lost_grains": int(in_lost.sum()),
            "unresolved_in_lost_grains": int((unres & in_lost).sum()),
            "unresolved_total": int(unres.sum()),
            "accepted_wrong": int(acc_wrong[r].sum()),
            "accepted_wrong_cross_boundary_like": int((acc_wrong[r] & v["cross"]).sum()),
            "accepted_total": int((~unres).sum()),
            "proj_full_h_voxel_linear": sc[r]["cost"]["wall_total_s"] / 3600.0 * n_full / n_reg,
            "proj_full_h_every_piece_seeded": n_pieces * c["wall_per_seed_s"] / 3600.0
            + nb_rev_h * n_full / n_reg,
            "proj_full_h_seeds_per_grain": seed_h * n_grain_full / len(gs)
            + nb_rev_h * n_full / n_reg,
        }
    summary = {
        "region": json.loads((ROOT / "benchmarks" / "phase_d_pilot" / "region.json").read_text()),
        "runs": sc,
        "paired_grain_clustered": pairs,
        "paired_grain_clustered_accepted_wrong": pairs_acc,
        "diagnostics": diag,
        "projection_note": (
            "full-sample projection (contended, single run): voxel-linear scales the whole "
            f"region wall by {n_full}/{n_reg}; seeds-per-grain scales the seed time by "
            f"{n_grain_full}/60 grains and the neighbour plus revisit time by {n_full}/{n_reg}; "
            f"every-piece-seeded charges {n_pieces} grain pieces (connected under the BFS "
            "neighbour radius of 2 sides) one seed each at the run's mean seed time (the cost if "
            "every unreached piece got a full search; neighbour and revisit time unchanged). "
            "Assumes the region is typical."
        ),
        "sample_counts": cnt,
        "seed_count_oracle": {
            "predicted_swallowed_grains": 140,
            "of_grains": 497,
            "source": "MIGRATION_HISTORY, Phase D seed-cost diagnosis, Seed count bullet",
        },
        "all_arms_share_one_seed_order": "every run uses default_rng(0) and shuffles the voxels "
        "before any other draw, so the five runs start from the same seed order",
        "timing_label": sc[next(iter(sc))]["cost"]["timing_label"],
    }
    (Path(a.out) / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")

    def row(r: str) -> Dict[str, Any]:
        s = sc[r]
        o = s["overall"]
        return {
            "run": r,
            "n": o["n"],
            "wrong": o["wrong"],
            "rate": pct(o["wrong_rate"]),
            "wilson": f"{pct(o['wilson_lo']):.1f}-{pct(o['wilson_hi']):.1f}",
            "grain_boot": f"{pct(o['cluster_boot_lo']):.1f}-{pct(o['cluster_boot_hi']):.1f}",
            "median_right": o["median_right_err_deg"],
            "unres": s["provenance"]["unresolved_final"],
            "found": s["grains"]["found"],
            "partial": s["grains"]["partial"],
            "lost": s["grains"]["lost"],
            "frag": s["grains"]["fragmented_found"],
            "wall_h": s["cost"]["wall_total_s"] / 3600.0,
        }

    t1 = (
        markdown_table(
            [row(r) for r in sc],
            [
                "run",
                "n",
                "wrong",
                "rate",
                "wilson",
                "grain_boot",
                "median_right",
                "unres",
                "found",
                "partial",
                "lost",
                "frag",
                "wall_h",
            ],
            formats={"rate": ".2f", "median_right": ".4f", "wall_h": ".2f"},
        )
        .replace("| rate |", "| wrong % |")
        .replace("| wilson |", "| Wilson 95% (%) |")
        .replace("| grain_boot |", "| grain-bootstrap 95% (%) |")
        .replace("| median_right |", "| median right err (deg) |")
        .replace("| unres |", "| unresolved |")
        .replace("| wall_h |", "| wall (h, contended) |")
        .replace("| frag |", "| fragmented (found grains) |")
    )

    def split_rows(key_a: str, key_b: str) -> List[Dict[str, Any]]:
        out = []
        for r in sc:
            for key in (key_a, key_b):
                o = sc[r][key]
                out.append(
                    {
                        "run": r,
                        "subset": key,
                        "n": o["n"],
                        "wrong": o["wrong"],
                        "rate": pct(o["wrong_rate"]),
                        "median_right": o["median_right_err_deg"],
                    }
                )
        return out

    fmt = {"rate": ".2f", "median_right": ".4f"}
    cols = ["run", "subset", "n", "wrong", "rate", "median_right"]
    t2 = markdown_table(split_rows("boundary", "interior"), cols, fmt)
    t3 = markdown_table(split_rows("near_axis", "far_axis"), cols, fmt)

    src_rows = []
    for r in sc:
        for nm, o in sc[r]["by_source"].items():
            src_rows.append(
                {
                    "run": r,
                    "source": nm,
                    "n": o["n"],
                    "wrong": o["wrong"],
                    "rate": pct(o["wrong_rate"]),
                    "median_right": o["median_right_err_deg"],
                }
            )
    t4 = markdown_table(src_rows, ["run", "source", "n", "wrong", "rate", "median_right"], fmt)

    def cost_row(r: str) -> Dict[str, Any]:
        c, p = sc[r]["cost"], sc[r]["provenance"]
        return {
            "run": r,
            "seeds": p["seeds"],
            "seed_rej": p["seed_rejected"],
            "nb_fits": p["neighbor_fits"],
            "nb_first": p["neighbor_accepted_first"],
            "retry": f"{p['retry_accepted']}/{p['retry_attempted']}",
            "revisit": f"{p['revisit_accepted']}/{p['revisit_attempted']}",
            "ev_seed": c["evals_seed"],
            "ev_nb": c["evals_neighbor"],
            "ev_rev": c["evals_revisit"],
            "w_seed": c["wall_seed_s"] / 60,
            "w_nb": c["wall_neighbor_s"] / 60,
            "w_rev": c["wall_revisit_s"] / 60,
            "s_per_seed": c["wall_per_seed_s"],
            "s_per_nb": c["wall_per_neighbor_fit_s"],
        }

    t5 = (
        markdown_table(
            [cost_row(r) for r in sc],
            [
                "run",
                "seeds",
                "seed_rej",
                "nb_fits",
                "nb_first",
                "retry",
                "revisit",
                "ev_seed",
                "ev_nb",
                "ev_rev",
                "w_seed",
                "w_nb",
                "w_rev",
                "s_per_seed",
                "s_per_nb",
            ],
            formats={
                "w_seed": ".1f",
                "w_nb": ".1f",
                "w_rev": ".1f",
                "s_per_seed": ".1f",
                "s_per_nb": ".3f",
                "ev_seed": ",d",
                "ev_nb": ",d",
                "ev_rev": ",d",
            },
        )
        .replace("| retry |", "| retry acc/att |")
        .replace("| revisit |", "| revisit acc/att |")
        .replace("| w_seed |", "| seed (min) |")
        .replace("| w_nb |", "| neighbour (min) |")
        .replace("| w_rev |", "| revisit (min) |")
        .replace("| s_per_seed |", "| s / seed |")
        .replace("| s_per_nb |", "| s / neighbour fit |")
    )

    def unres_rows() -> List[Dict[str, Any]]:
        out = []
        for r in sc:
            u = sc[r]["unresolved"]
            e = u["errors_deg"]
            out.append(
                {
                    "run": r,
                    "unresolved": u["n"],
                    "wrong": u["wrong"],
                    "median_err": float(np.median(e)) if e else None,
                    "max_err": max(e) if e else None,
                    "cross": sc[r]["wrong_cross_boundary_like"],
                    "wrong_all": sc[r]["wrong_total"],
                    "wrong_boundary": sc[r]["wrong_boundary"],
                }
            )
        return out

    t6 = markdown_table(
        unres_rows(),
        [
            "run",
            "unresolved",
            "wrong",
            "median_err",
            "max_err",
            "wrong_all",
            "wrong_boundary",
            "cross",
        ],
        {"median_err": ".3f", "max_err": ".3f"},
    )
    t7 = markdown_table(
        [{"pair": k, **{kk: v for kk, v in d.items()}} for k, d in pairs.items()],
        ["pair", "wrong_a", "wrong_b", "n_clusters", "n_clusters_differing", "p"],
        {"p": ".3g"},
    )
    d_rows = [{"run": r, **d} for r, d in diag.items()]
    t8 = markdown_table(
        d_rows,
        [
            "run",
            "n_grains",
            "grains_with_an_accepted_seed",
            "lost_grains",
            "lost_grains_no_full_search",
            "lost_grains_with_a_rejected_seed",
            "lost_grains_with_an_accepted_seed",
            "voxels_in_lost_grains",
            "unresolved_in_lost_grains",
            "unresolved_total",
            "accepted_total",
            "accepted_wrong",
            "accepted_wrong_cross_boundary_like",
        ],
    )
    t11 = markdown_table(
        d_rows,
        [
            "run",
            "lost_grains_rate",
            "lost_grains_wilson_lo",
            "lost_grains_wilson_hi",
            "lost_median_region_voxels",
            "found_median_region_voxels",
            "lost_boundary_voxel_share",
            "found_boundary_voxel_share",
            "grains_wholly_in_region",
            "wholly_in_region_lost",
            "grains_truncated_by_region",
            "truncated_lost",
            "accepted_wrong_rate",
        ],
        {
            "lost_grains_rate": ".3f",
            "lost_grains_wilson_lo": ".3f",
            "lost_grains_wilson_hi": ".3f",
            "lost_median_region_voxels": ".1f",
            "found_median_region_voxels": ".1f",
            "lost_boundary_voxel_share": ".3f",
            "found_boundary_voxel_share": ".3f",
            "accepted_wrong_rate": ".4f",
        },
    )
    t9 = markdown_table(
        [{"run": r, **d} for r, d in diag.items()],
        [
            "run",
            "proj_full_h_voxel_linear",
            "proj_full_h_seeds_per_grain",
            "proj_full_h_every_piece_seeded",
        ],
        {
            "proj_full_h_voxel_linear": ".1f",
            "proj_full_h_seeds_per_grain": ".1f",
            "proj_full_h_every_piece_seeded": ".1f",
        },
    )
    t10 = markdown_table(
        [{"pair": k, **d} for k, d in pairs_acc.items()],
        ["pair", "wrong_a", "wrong_b", "n_clusters", "n_clusters_differing", "p"],
        {"p": ".3g"},
    )
    write_tables(
        Path(a.out) / "tables.md",
        {
            "phase_d_pilot_main": t1,
            "phase_d_pilot_boundary": t2,
            "phase_d_pilot_axis": t3,
            "phase_d_pilot_source": t4,
            "phase_d_pilot_cost": t5,
            "phase_d_pilot_unresolved": t6,
            "phase_d_pilot_pairs": t7,
            "phase_d_pilot_diag": t8,
            "phase_d_pilot_lost_grains": t11,
            "phase_d_pilot_projection": t9,
            "phase_d_pilot_pairs_accepted": t10,
        },
    )
    print(t1)


if __name__ == "__main__":
    main()
