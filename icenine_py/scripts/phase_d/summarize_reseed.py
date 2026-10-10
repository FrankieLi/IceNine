"""Compare the BFS-reseed re-pilot with the Phase D pilot (same region, images and seed order).

Usage (from icenine_py/):
  uv run python scripts/phase_d/summarize_reseed.py [--runs benchmarks/phase_d_pilot_reseed/runs]
Reads <runs>/<arm>_<variant>/{run.json,records.json,score.json,voxels.npz} (score_bfs.py) and the
pilot's benchmarks/phase_d_pilot/runs; writes benchmarks/phase_d_pilot_reseed/{summary.json,
tables.md}. Timings are contended (4 concurrent single-process runs). P-values are
grain-clustered sign-flip tests on per-grain differences of the wrong counts.

Pilot counterpart of each re-pilot run: mc_clean -> mc_clean, cma_noretry_clean ->
cma_noretry_clean, mc_realistic_q16 -> mc_realistic_q16, and cma_noretry_realistic_q16 ->
cma_realistic_q16 (the pilot had no cma_noretry run on realistic_q16; its cma arm had the retry).
"""

import argparse
import json
import sys
from pathlib import Path
from typing import Any, Dict, List

import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
from scipy.spatial import cKDTree

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(ROOT / "scripts" / "common"))

import score_bfs as S  # noqa: E402
from doc_tables import markdown_table, write_tables  # noqa: E402

PAIRS = {
    "mc_clean": "mc_clean",
    "cma_noretry_clean": "cma_noretry_clean",
    "mc_realistic_q16": "mc_realistic_q16",
    "cma_noretry_realistic_q16": "cma_realistic_q16",
}


def unresolved_clusters(
    run_dir: Path, recs: Dict[int, Any], grain: np.ndarray, lost: set
) -> Dict[str, Any]:
    """Connected clusters (radius 2 x side) of the voxels that end unresolved, and how many hold
    a rejected ordinary seed (those are final by the rule) or a rejected reseed voxel."""
    from icenine.mic_file import MicFile

    mic = MicFile.read(str(run_dir / "recon.mic"))
    pos = np.array([v.position[:2] for v in mic.voxels], dtype=float)
    side = mic.voxels[0].side_length
    unres = np.array([i for i, r in recs.items() if r["source"] == "unresolved"], dtype=int)
    n = len(pos)
    pr = cKDTree(pos).query_pairs(2.0 * side, output_type="ndarray")
    keep = np.isin(pr[:, 0], unres) & np.isin(pr[:, 1], unres)
    adj = coo_matrix((np.ones(int(keep.sum())), (pr[keep, 0], pr[keep, 1])), shape=(n, n))
    lab = connected_components(adj, directed=False)[1]
    seed_rej = {i for i, r in recs.items() if r["seed_rejected"] and not r.get("reseeded")}
    reseed_rej = {i for i, r in recs.items() if r["seed_rejected"] and r.get("reseeded")}
    sizes: Dict[int, int] = {}
    for i in unres:
        sizes[int(lab[i])] = sizes.get(int(lab[i]), 0) + 1
    rows = []
    for c, sz in sorted(sizes.items(), key=lambda t: -t[1]):
        members = set(unres[lab[unres] == c].tolist())
        rows.append(
            {
                "size": sz,
                "has_rejected_ordinary_seed": bool(members & seed_rej),
                "has_rejected_reseed": bool(members & reseed_rej),
                "grains": len({int(grain[i]) for i in members}),
                "lost_grains_inside": len({int(grain[i]) for i in members} & lost),
            }
        )
    return {"n_clusters": len(rows), "clusters": rows}


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--runs", default=str(ROOT / "benchmarks" / "phase_d_pilot_reseed" / "runs"))
    ap.add_argument("--pilot", default=str(ROOT / "benchmarks" / "phase_d_pilot" / "runs"))
    ap.add_argument(
        "--v1", default=None, help="summary.json of the first (per-cluster finality) reseed run"
    )
    ap.add_argument("--prefix", default="phase_d_reseed", help="table name prefix")
    ap.add_argument("--out", default=str(ROOT / "benchmarks" / "phase_d_pilot_reseed"))
    a = ap.parse_args()
    runs, pilot, out = Path(a.runs), Path(a.pilot), Path(a.out)

    v1 = json.loads(Path(a.v1).read_text())["runs"] if a.v1 else {}
    summary: Dict[str, Any] = {"runs": {}, "pairs": {}}
    cmp_rows: List[Dict[str, Any]] = []
    cnt_rows: List[Dict[str, Any]] = []
    log_rows: List[Dict[str, Any]] = []
    for new, old in PAIRS.items():
        if not (runs / new / "score.json").exists():
            continue
        sn = json.loads((runs / new / "score.json").read_text())
        so = json.loads((pilot / old / "score.json").read_text())
        vn, vo = np.load(runs / new / "voxels.npz"), np.load(pilot / old / "voxels.npz")
        run = json.loads((runs / new / "run.json").read_text())
        recs = {int(k): v for k, v in json.loads((runs / new / "records.json").read_text()).items()}
        st = run["stats"]
        assert run["bfs_reseed"] is True
        assert (vn["grain"] == vo["grain"]).all()
        acc_wrong = lambda v: v["wrong"] & (v["source"] != "unresolved")  # noqa: E731
        pair = S.clustered_paired_p(vo["wrong"], vn["wrong"], vo["grain"])
        pair_acc = S.clustered_paired_p(acc_wrong(vo), acc_wrong(vn), vo["grain"])
        # lost grains that the reseed recovered, and grains lost only in the re-pilot
        lost_o, lost_n = set(so["grains"]["lost_ids"]), set(sn["grains"]["lost_ids"])
        reseeded = [i for i, r in recs.items() if r.get("reseeded")]
        reseed_log = st["reseeds"]
        n_acc = sum(1 for e in reseed_log if e["accepted"])
        n_rej = len(reseed_log) - n_acc
        sizes = [e["cluster_size"] for e in reseed_log]
        wall = st["wall_total_s"]
        ucl = unresolved_clusters(runs / new, recs, vn["grain"], lost_n)
        ev_tot = (
            st["n_evals_seed"]
            + st["n_evals_reseed"]
            + st["n_evals_neighbor"]
            + st["n_evals_revisit"]
        )
        ev_reseed_all = sum(e["evaluations"] for e in reseed_log)
        unres_n = int((vn["source"] == "unresolved").sum())
        # why the remaining lost grains are lost: a failed full search in the grain, or none at all
        searched = {int(vn["grain"][i]) for i, r in recs.items() if r["seed_rejected"]}
        accepted_full = {
            int(vn["grain"][i])
            for i, r in recs.items()
            if r["source"] == "seed" and not r["seed_rejected"]
        }
        lost_failed = len(lost_n & searched)
        lost_nosearch = len(lost_n - searched - accepted_full)
        unres_o = int((vo["source"] == "unresolved").sum())
        summary["runs"][new] = {
            "pilot_counterpart": old,
            "wrong_pilot": int(vo["wrong"].sum()),
            "wrong_reseed": int(vn["wrong"].sum()),
            "unresolved_pilot": unres_o,
            "unresolved_reseed": unres_n,
            "lost_pilot": so["grains"]["lost"],
            "lost_reseed": sn["grains"]["lost"],
            "lost_pilot_ids": sorted(lost_o),
            "lost_reseed_ids": sorted(lost_n),
            "lost_in_pilot_not_lost_now": len(lost_o - lost_n),
            "lost_now_not_lost_in_pilot": len(lost_n - lost_o),
            "accepted_wrong_pilot": int(acc_wrong(vo).sum()),
            "accepted_wrong_reseed": int(acc_wrong(vn).sum()),
            "n_reseeded_voxels": len(reseeded),
            "lost_with_a_failed_full_search": lost_failed,
            "lost_without_any_full_search": lost_nosearch,
            "lost_with_an_accepted_full_search": len(lost_n & accepted_full),
            "failed_voxels_total": sum(1 for r in recs.values() if r["seed_rejected"]),
            "v1": v1.get(new),
            "unresolved_clusters": ucl,
            "evals_total": ev_tot,
            "evals_reseed_with_expansion": ev_reseed_all,
            "reseed_share_of_evals": ev_reseed_all / ev_tot,
            "reseeds": reseed_log,
            "wall_total_s": wall,
            "wall_reseed_total_s": st["wall_reseed_total_s"],
            "wall_reseed_full_search_s": st["wall_reseed_s"],
            "reseed_share_of_wall": st["wall_reseed_total_s"] / wall,
            "reseed_full_search_share_of_wall": st["wall_reseed_s"] / wall,
            "n_seeds": st["n_seeds"],
            "n_seed_rejected": st["n_seed_rejected"],
            "n_evals_seed": st["n_evals_seed"],
            "n_evals_reseed": st["n_evals_reseed"],
            "timing_label": run["timing_label"],
            "torch_threads": run["torch_threads"],
            "overall_reseed": sn["overall"],
            "overall_pilot": so["overall"],
            "grains_reseed": {k: v for k, v in sn["grains"].items() if not k.endswith("_ids")},
            "grains_pilot": {k: v for k, v in so["grains"].items() if not k.endswith("_ids")},
        }
        summary["pairs"][new] = {"vs": old, "all_voxels": pair, "accepted_wrong": pair_acc}
        cmp_rows.append(
            {
                "run": new,
                "pilot": old,
                "wrong_pilot": int(vo["wrong"].sum()),
                "wrong_reseed": int(vn["wrong"].sum()),
                "unres_pilot": unres_o,
                "unres_reseed": unres_n,
                "lost_pilot": so["grains"]["lost"],
                "lost_reseed": sn["grains"]["lost"],
                "wrong_v1": v1[new]["wrong_reseed"] if new in v1 else None,
                "unres_v1": v1[new]["unresolved_reseed"] if new in v1 else None,
                "lost_v1": v1[new]["lost_reseed"] if new in v1 else None,
                "acc_wrong_pilot": int(acc_wrong(vo).sum()),
                "acc_wrong_reseed": int(acc_wrong(vn).sum()),
                "p": pair["p"],
                "n_clusters_differing": pair["n_clusters_differing"],
            }
        )
        cnt_rows.append(
            {
                "run": new,
                "seeds": st["n_seeds"],
                "seed_rej": st["n_seed_rejected"],
                "reseeds": len(reseed_log),
                "accepted": n_acc,
                "rejected": n_rej,
                "median_cluster": float(np.median(sizes)) if sizes else None,
                "lost_failed": lost_failed,
                "lost_nosearch": lost_nosearch,
                "final_voxels": st["n_reseed_final_voxels"],
                "reseed_min": st["wall_reseed_total_s"] / 60.0,
                "full_search_min": st["wall_reseed_s"] / 60.0,
                "share_pct": 100.0 * st["wall_reseed_total_s"] / wall,
                "wall_h": wall / 3600.0,
                "wall_pilot_h": so["cost"]["wall_total_s"] / 3600.0,
                "ev_reseed": ev_reseed_all,
                "ev_share_pct": 100.0 * ev_reseed_all / ev_tot,
                "unres_clusters": ucl["n_clusters"],
                "largest_unres": ucl["clusters"][0]["size"] if ucl["clusters"] else 0,
                "largest_grains": ucl["clusters"][0]["grains"] if ucl["clusters"] else 0,
                "largest_lost_inside": (
                    ucl["clusters"][0]["lost_grains_inside"] if ucl["clusters"] else 0
                ),
                "largest_has_rej_seed": (
                    ucl["clusters"][0]["has_rejected_ordinary_seed"] if ucl["clusters"] else False
                ),
            }
        )
        for k, e in enumerate(reseed_log):
            log_rows.append(
                {
                    "run": new,
                    "k": k + 1,
                    "cluster": e["cluster_size"],
                    "voxel": e["voxel"],
                    "accepted": "yes" if e["accepted"] else "no",
                    "left_refit": e["cluster_voxels_still_refit"],
                    "seconds": e["seconds"],
                    "evals": e["evaluations"],
                }
            )
    (out / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    t1 = markdown_table(
        cmp_rows,
        [
            "run",
            "pilot",
            "wrong_pilot",
            "wrong_reseed",
            "unres_pilot",
            "unres_reseed",
            "lost_pilot",
            "lost_v1",
            "lost_reseed",
            "wrong_v1",
            "unres_v1",
            "acc_wrong_pilot",
            "acc_wrong_reseed",
            "n_clusters_differing",
            "p",
        ],
        {"p": ".3g"},
    )
    t2 = markdown_table(
        cnt_rows,
        [
            "run",
            "seeds",
            "seed_rej",
            "reseeds",
            "accepted",
            "rejected",
            "median_cluster",
            "final_voxels",
            "lost_failed",
            "lost_nosearch",
            "reseed_min",
            "full_search_min",
            "share_pct",
            "ev_share_pct",
            "unres_clusters",
            "largest_unres",
            "largest_grains",
            "largest_lost_inside",
            "largest_has_rej_seed",
            "wall_h",
            "wall_pilot_h",
            "ev_reseed",
        ],
        {
            "median_cluster": ".1f",
            "reseed_min": ".1f",
            "full_search_min": ".1f",
            "share_pct": ".1f",
            "ev_share_pct": ".1f",
            "wall_h": ".2f",
            "wall_pilot_h": ".2f",
            "ev_reseed": ",d",
        },
    )
    t3 = markdown_table(
        log_rows,
        ["run", "k", "cluster", "voxel", "accepted", "left_refit", "seconds", "evals"],
        {"seconds": ".0f", "evals": ",d"},
    )
    three = [
        {
            "run": r["run"],
            "wrong_pilot": r["wrong_pilot"],
            "wrong_v1": r["wrong_v1"],
            "wrong_v2": r["wrong_reseed"],
            "unres_pilot": r["unres_pilot"],
            "unres_v1": r["unres_v1"],
            "unres_v2": r["unres_reseed"],
            "lost_pilot": r["lost_pilot"],
            "lost_v1": r["lost_v1"],
            "lost_v2": r["lost_reseed"],
        }
        for r in cmp_rows
    ]
    t0 = markdown_table(
        three,
        [
            "run",
            "wrong_pilot",
            "wrong_v1",
            "wrong_v2",
            "unres_pilot",
            "unres_v1",
            "unres_v2",
            "lost_pilot",
            "lost_v1",
            "lost_v2",
        ],
    )
    write_tables(
        out / "tables.md",
        {
            f"{a.prefix}_three": t0,
            f"{a.prefix}_compare": t1,
            f"{a.prefix}_counts": t2,
            f"{a.prefix}_log": t3,
        },
    )
    print(t1)
    print(t2)


if __name__ == "__main__":
    main()
