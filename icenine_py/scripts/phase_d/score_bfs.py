"""Score one Phase D pilot BFS run against the truth (symmetry-reduced misorientation).

Pure helpers (tested in tests/test_phase_d_scoring.py) plus a CLI:
  uv run python scripts/phase_d/score_bfs.py --run DIR/<arm>_<variant> [--out score.json]

Per run it writes score.json (numbers) and voxels.npz (per-voxel error, source, flags) next to
the run. Wrong means a symmetry-reduced misorientation above 1 degree. Interval conventions:
`wilson` treats voxels as independent (optimistic, voxels of a grain are correlated);
`cluster_boot` is a percentile bootstrap over grains. P-values between arms are grain-clustered
(sign-flip test on per-grain differences of the wrong counts, `clustered_paired_p`).
"""

import argparse
import json
import sys
from pathlib import Path
from typing import Any, Dict, List, Sequence

import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
from scipy.spatial import cKDTree

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(ROOT / "scripts" / "common"))

import grains as G  # noqa: E402
from stats import wilson  # noqa: E402

WRONG_DEG = 1.0
SOURCES = ["seed", "neighbor", "neighbor_retry", "revisit", "restart", "unresolved"]


def pair_misorientation_deg(qa: np.ndarray, qb: np.ndarray, sym: np.ndarray) -> np.ndarray:
    """Symmetry-reduced misorientation (deg) of row i of qa with row i of qb; (N,4) [w,x,y,z]."""
    conj = qa * np.array([1.0, -1.0, -1.0, -1.0])
    rel = G._qmul(conj, qb)  # (N, 4)
    prod = G._qmul(rel[:, None, :], np.asarray(sym)[None, :, :])  # (N, S, 4)
    w = np.abs(prod[..., 0]).max(axis=1)
    return np.degrees(2.0 * np.arccos(np.clip(w, 0.0, 1.0)))


def mats_to_quats(mats: np.ndarray) -> np.ndarray:
    from icenine.sampling import matrix_to_quaternion

    return np.array([matrix_to_quaternion(np.asarray(m, dtype=np.float64)) for m in mats])


def wrong_summary(err_deg: np.ndarray, thr: float = WRONG_DEG) -> Dict[str, Any]:
    """Count, rate, Wilson interval and the median error of the right (<= thr) answers."""
    err = np.asarray(err_deg, dtype=float)
    n = int(len(err))
    k = int((err > thr).sum())
    right = err[err <= thr]
    out: Dict[str, Any] = {"n": n, "wrong": k}
    if n:
        lo, hi = wilson(k, n)
        out.update(wrong_rate=k / n, wilson_lo=lo, wilson_hi=hi)
    else:
        out.update(wrong_rate=None, wilson_lo=None, wilson_hi=None)
    out["n_right"] = int(len(right))
    out["median_right_err_deg"] = float(np.median(right)) if len(right) else None
    out["median_err_deg"] = float(np.median(err)) if n else None
    return out


def cluster_bootstrap_rate(
    wrong: np.ndarray, cluster: np.ndarray, n_boot: int = 4000, seed: int = 0
) -> Dict[str, float]:
    """Percentile 95% interval of the wrong rate, resampling clusters (grains) with replacement."""
    wrong = np.asarray(wrong, dtype=float)
    ids, inv = np.unique(cluster, return_inverse=True)
    k = np.bincount(inv, weights=wrong, minlength=len(ids))
    n = np.bincount(inv, minlength=len(ids)).astype(float)
    rng = np.random.default_rng(seed)
    pick = rng.integers(0, len(ids), size=(n_boot, len(ids)))
    rates = k[pick].sum(axis=1) / n[pick].sum(axis=1)
    return {
        "cluster_boot_lo": float(np.percentile(rates, 2.5)),
        "cluster_boot_hi": float(np.percentile(rates, 97.5)),
    }


def clustered_paired_p(
    wrong_a: np.ndarray,
    wrong_b: np.ndarray,
    cluster: np.ndarray,
    n_perm: int = 20000,
    seed: int = 0,
) -> Dict[str, Any]:
    """Two-sided sign-flip test on per-cluster differences (wrong_a - wrong_b summed per grain).

    Voxels of one grain are not independent, so the unit is the grain. Exact enumeration when
    there are at most 20 non-zero clusters (2^20 patterns, done in chunks), else n_perm random
    sign flips (p = (1+hits)/(1+n), so its floor is 1/(1+n_perm))."""
    wa, wb = np.asarray(wrong_a, dtype=float), np.asarray(wrong_b, dtype=float)
    ids, inv = np.unique(cluster, return_inverse=True)
    d = np.bincount(inv, weights=wa - wb, minlength=len(ids))
    nz = d[d != 0]
    obs = abs(float(nz.sum()))
    out: Dict[str, Any] = {
        "wrong_a": int(wa.sum()),
        "wrong_b": int(wb.sum()),
        "n_clusters": int(len(ids)),
        "n_clusters_differing": int(len(nz)),
    }
    if len(nz) == 0:
        out["p"] = 1.0
        return out
    if len(nz) <= 20:
        hits, total = 0, 2 ** len(nz)
        for lo in range(0, total, 1 << 16):
            ids_ = np.arange(lo, min(lo + (1 << 16), total))
            signs = ((ids_[:, None] >> np.arange(len(nz))[None, :]) & 1) * 2 - 1
            hits += int((np.abs(signs @ nz) >= obs - 1e-9).sum())
        out["p"] = hits / total
        out["exact"] = True
    else:
        rng = np.random.default_rng(seed)
        signs = rng.integers(0, 2, size=(n_perm, len(nz))) * 2 - 1
        stat = np.abs(signs @ nz)
        out["p"] = float((1 + (stat >= obs - 1e-9).sum()) / (1 + n_perm))
    return out


def grain_status(
    err: np.ndarray, grain: np.ndarray, quats: np.ndarray, sym: np.ndarray, thr: float = WRONG_DEG
) -> Dict[str, Any]:
    """Per-grain outcome in the scored voxel set.

    found: more than half of the grain's voxels are right; partial: some but not more than half;
    lost: none right. fragmented: the reconstructed orientations of the grain's voxels form two or
    more single-linkage clusters at `thr` degrees (so a grain can be both found and fragmented).
    """
    found = partial = lost = frag = 0
    lost_ids: List[int] = []
    frag_ids: List[int] = []
    found_ids: List[int] = []
    for g in np.unique(grain):
        m = grain == g
        n_right = int((err[m] <= thr).sum())
        if n_right == 0:
            lost += 1
            lost_ids.append(int(g))
        elif n_right * 2 > m.sum():
            found += 1
            found_ids.append(int(g))
        else:
            partial += 1
        q = quats[m]
        if len(q) > 1:
            mm = G.misorientation_matrix_deg(q, q, sym)
            n_comp, _ = connected_components(coo_matrix(mm <= thr), directed=False)
            if n_comp >= 2:
                frag += 1
                frag_ids.append(int(g))
    return {
        "n_grains": int(len(np.unique(grain))),
        "found": found,
        "partial": partial,
        "lost": lost,
        "fragmented": frag,
        "fragmented_found": len(set(frag_ids) & set(found_ids)),
        "lost_ids": lost_ids,
        "fragmented_ids": frag_ids,
    }


def cross_boundary_wrong(
    err: np.ndarray,
    grain: np.ndarray,
    pos: np.ndarray,
    quats: np.ndarray,
    truth_quats: np.ndarray,
    sym: np.ndarray,
    radius: float,
    thr: float = WRONG_DEG,
) -> np.ndarray:
    """Per voxel: wrong AND its orientation matches (<= thr) the truth of a neighbouring voxel
    (within `radius`) that belongs to a different grain, i.e. an error that looks like a
    neighbouring grain's orientation (a propagated cross-boundary error)."""
    flag = np.zeros(len(err), dtype=bool)
    tree = cKDTree(pos)
    for i in np.nonzero(err > thr)[0]:
        nb = [j for j in tree.query_ball_point(pos[i], radius) if grain[j] != grain[i]]
        if not nb:
            continue
        d = pair_misorientation_deg(
            np.repeat(quats[i : i + 1], len(nb), axis=0), truth_quats[nb], sym
        )
        flag[i] = bool((d <= thr).any())
    return flag


def by_group(err: np.ndarray, labels: Sequence[str], names: Sequence[str]) -> Dict[str, Any]:
    lab = np.asarray(labels)
    return {nm: wrong_summary(err[lab == nm]) for nm in names if (lab == nm).any()}


def score_run(run_dir: Path) -> Dict[str, Any]:
    from icenine.mic_file import MicFile
    from icenine.orientation_search import get_symmetry_quaternions
    from icenine.symmetry import create_cubic_symmetry

    sym = get_symmetry_quaternions(create_cubic_symmetry(1.0))
    bench = ROOT / "benchmarks" / "phase_d_pilot"
    ex = ROOT.parent / "Examples" / "Example2.ManyGrains" / "SimInput"
    idx = np.load(bench / "region_index.npy")
    region = json.loads((bench / "region.json").read_text())
    side = region["side_length_mm"]
    truth = MicFile.read(str(ex / "rand_500grains_1mm_neworient_s0.mic"))
    full_grain = np.load(ex / "rand_500grains_1mm_neworient_s0_grainmap.npy")
    recon = MicFile.read(str(run_dir / "recon.mic"))
    run = json.loads((run_dir / "run.json").read_text())
    records = {int(k): v for k, v in json.loads((run_dir / "records.json").read_text()).items()}
    assert len(recon.voxels) == len(idx)

    tv = [truth.voxels[i] for i in idx]
    pos = np.array([v.position[:2] for v in tv], dtype=float)
    rpos = np.array([v.position[:2] for v in recon.voxels], dtype=float)
    assert np.abs(pos - rpos).max() < 1e-4, "reconstructed voxel order differs from the region"
    q_true = mats_to_quats(np.array([v.orientation for v in tv]))
    q_rec = mats_to_quats(np.array([v.orientation for v in recon.voxels]))
    err = pair_misorientation_deg(q_rec, q_true, sym)
    grain = full_grain[idx]

    # boundary flag from the full sample (a different-grain voxel within 1.01 side of the centroid)
    from make_region import boundary_flags, centroids_from_lines  # noqa: E402

    lines = (ex / "rand_500grains_1mm_neworient_s0_grid.mic").read_text().splitlines()
    raw = np.array([[float(t) for t in ln.split()[:6]] for ln in lines[1:] if ln.strip()])
    cent, _ = centroids_from_lines(raw, float(lines[0].split()[0]))
    bnd = boundary_flags(cent, full_grain, side)[idx]
    r = np.hypot(cent[idx, 0], cent[idx, 1])
    near = r < region["near_axis_r_mm"]

    n = len(idx)
    visited = np.array([i in records for i in range(n)])
    src = np.array([records[i]["source"] if i in records else "not_visited" for i in range(n)])
    wrong = err > WRONG_DEG
    cross = cross_boundary_wrong(err, grain, cent[idx], q_rec, q_true, sym, 1.01 * side)

    scored = visited  # smoke runs score the voxels that BFS reached
    s = {
        "run": run_dir.name,
        "n_region": n,
        "n_scored": int(scored.sum()),
        "overall": wrong_summary(err[scored]),
        "boundary": wrong_summary(err[scored & bnd]),
        "interior": wrong_summary(err[scored & ~bnd]),
        "near_axis": wrong_summary(err[scored & near]),
        "far_axis": wrong_summary(err[scored & ~near]),
        "grains": grain_status(err[scored], grain[scored], q_rec[scored], sym),
        "by_source": by_group(err[scored], src[scored], SOURCES),
        "wrong_cross_boundary_like": int((cross & scored).sum()),
        "wrong_total": int((wrong & scored).sum()),
        "wrong_boundary": int((wrong & scored & bnd).sum()),
        "wrong_near_axis": int((wrong & scored & near).sum()),
        "unresolved": {
            "n": int((src == "unresolved").sum()),
            "errors_deg": [float(x) for x in np.sort(err[src == "unresolved"])],
            "wrong": int(((src == "unresolved") & wrong).sum()),
        },
    }
    s["overall"].update(cluster_bootstrap_rate(wrong[scored], grain[scored]))
    sts = run["stats"]
    s["provenance"] = {
        "seeds": sts["n_seeds"],
        "seed_rejected": sts["n_seed_rejected"],
        "neighbor_fits": sts["n_neighbor_fits"],
        "neighbor_accepted_first": sts["n_neighbor_accepted_first"],
        "retry_attempted": sts["n_retry_attempted_neighbor"] + sts["n_retry_attempted_revisit"],
        "retry_accepted": sts["n_retry_accepted_neighbor"] + sts["n_retry_accepted_revisit"],
        "revisit_attempted": sts["n_revisit_attempted"],
        "revisit_accepted": sts["n_revisit_accepted"],
        "revisit_rejected": sts["n_revisit_rejected"],
        "revisit_capped": sts["n_revisit_capped"],
        "unresolved_final": sts["n_unresolved"],
        "source_counts": {nm: int((src == nm).sum()) for nm in SOURCES + ["not_visited"]},
    }
    s["cost"] = {
        "evals_seed": sts["n_evals_seed"],
        "evals_neighbor": sts["n_evals_neighbor"],
        "evals_revisit": sts["n_evals_revisit"],
        "wall_seed_s": sts["wall_seed_s"],
        "wall_neighbor_s": sts["wall_neighbor_s"],
        "wall_revisit_s": sts["wall_revisit_s"],
        "wall_total_s": sts["wall_total_s"],
        "wall_per_seed_s": sts["wall_seed_s"] / max(sts["n_seeds"], 1),
        "wall_per_neighbor_fit_s": sts["wall_neighbor_s"] / max(sts["n_neighbor_fits"], 1),
        "timing_label": run.get("timing_label"),
    }
    s["options"] = {
        k: run[k]
        for k in (
            "local_optimizer",
            "cma_neighbor_max_evals",
            "cma_retry_sigma0_deg",
            "bfs_revisit_refit",
            "bfs_revisit_max",
            "bfs_rng_seed",
            "max_voxels",
        )
    }
    np.savez(
        run_dir / "voxels.npz",
        seed_rejected=np.array(
            [bool(records[i]["seed_rejected"]) if i in records else False for i in range(n)]
        ),
        region_pos=idx,
        err_deg=err,
        grain=grain,
        source=src,
        boundary=bnd,
        near_axis=near,
        wrong=wrong,
        scored=scored,
        cross=cross,
        n_evals=np.array([records[i]["n_evals"] if i in records else 0 for i in range(n)]),
        wall_s=np.array([records[i]["wall_s"] if i in records else 0.0 for i in range(n)]),
    )
    return s


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--run", required=True)
    ap.add_argument(
        "--timing-label",
        default="contended: 5 concurrent single-process runs on a 12-core machine, plus 2 A/B "
        "processes for about 3 min, load average 3.0 at launch (preflight.json)",
    )
    a = ap.parse_args()
    run_dir = Path(a.run)
    s = score_run(run_dir)
    s["cost"]["timing_label"] = a.timing_label
    (run_dir / "score.json").write_text(json.dumps(s, indent=2) + "\n")
    print(json.dumps({k: s[k] for k in ("run", "n_scored", "overall", "grains")}, indent=2))


if __name__ == "__main__":
    main()
