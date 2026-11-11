"""AR3 step 8: build the slim, committed result files from the replay and logged-run caches.

  RESEED_SAMPLE=big  uv run python scripts/autoresearch/ar3/reseed_slim.py --steps S5 D4 [--logged]
  RESEED_SAMPLE=grad uv run python scripts/autoresearch/ar3/reseed_slim.py --steps S5 D4

``--steps``: merges the shards of ``reseed_replay.py`` for each run (for ``S5`` on ``big`` the
logged run's seed steps are used when ``--steps`` has ``S5L``) into
``reseed/slim/steps_<sample>_<run>.json.gz``: one row per seed step with the voxel's truth level,
grain, distance to the nearest voxel of another true grain (in voxel widths), and the replay's
detail. The check against the stored run (evaluation count per seed step) is written to
``reseed/slim/check_<sample>.json``.

``--logged`` (big only): reads ``reseed/cache/S5_big_logged`` and writes
``reseed/slim/clusters_big.json.gz`` (one row per reseed cluster), ``members_big.json.gz`` (one row
per voxel of a cluster, with its class A / B / C / N), ``rejected_big.json.gz`` (one row per
rejected local fit of the run) and ``check_logged_big.json`` (the logged run against the stored
one).

Terms. A fit is rejected when its hit ratio / max(flood best before, hit ratio) is at most 0.9 (the
BFS's relative test, the classic mode of S5 has no retry). Classes of the last rejected fit of a
cluster voxel before its reseed, judged with the truth: A, the fit ended right (error at most the
step 7 limit) and was rejected; B, not A and the start came from an expander of the voxel's own true
grain; C, the start came from an expander of another true grain; N, no fit row (should not occur).
"""

import argparse
import glob
import json
from typing import Any, Dict, List

import numpy as np
import reseed_common as RC

import gradient_common as G  # noqa: E402

SLIM = RC.OUT / "slim"
SLIM.mkdir(parents=True, exist_ok=True)
START_BINS = (0.1, 0.2, 0.5)  # degrees, class B split by the start's angle from the truth


def geometry(T: Any) -> Dict[str, Any]:
    """Distance of each voxel to the nearest voxel of another true grain, in voxel widths."""
    from scipy.spatial import cKDTree

    xy = np.asarray(T.xy, float)
    g = np.asarray(T.gmap)
    from icenine.mic_file import MicFile

    # the layer is a triangle mesh with two voxels at about each position: the unit is the side
    # length of a voxel, read from the truth file
    spacing = float(MicFile.read(str(G.TRUTH_MIC)).voxels[0].side_length)
    dist = np.full(len(xy), np.inf)
    for gid in np.unique(g):
        inside = g == gid
        if inside.all():
            continue
        tr = cKDTree(xy[~inside])
        dist[inside] = tr.query(xy[inside], k=1)[0]
    return dict(dist=dist / spacing, spacing=spacing)


def merge_shards(run: str) -> List[Dict[str, Any]]:
    rows: List[Dict[str, Any]] = []
    for p in sorted(glob.glob(str(RC.CACHE / f"replay_{RC.SAMPLE}_{run}_*of*.json.gz"))):
        if "limit" in p:
            continue
        rows += RC.read_gz(p)
    return sorted(rows, key=lambda r: (r["voxel"], r["reseed"]))


def steps(runs: List[str], T: Any, geo: Dict[str, Any]) -> None:
    chk: Dict[str, Any] = {}
    for run in runs:
        if run == "S5L":
            rows = RC.read_gz(RC.CACHE / "S5_big_logged" / "seeds.json.gz")
            stored = {(s["voxel"], s["reseed"]): s for s in RC.seed_steps("S5")}
            for r in rows:
                s = stored[(r["voxel"], r["reseed"])]
                r["stored_evals"], r["grid"] = s["stored_evals"], s["grid"]
            name = "S5L"
        else:
            rows = merge_shards(run)
            name = run
        n_eq = sum(1 for r in rows if r["evals"] == r["stored_evals"])
        n_grid = sum(1 for r in rows if r["grid"])
        n_eq_nogrid = sum(1 for r in rows if (not r["grid"]) and r["evals"] == r["stored_evals"])
        n_nogrid = sum(1 for r in rows if not r["grid"])
        stored_n = len(RC.seed_steps("S5" if run == "S5L" else run))
        chk[name] = dict(
            rows=len(rows),
            stored_steps=stored_n,
            evals_replay=int(sum(r["evals"] for r in rows)),
            evals_stored=int(sum(r["stored_evals"] for r in rows)),
            equal=n_eq,
            grid_steps=n_grid,
            non_grid_steps=n_nogrid,
            non_grid_equal=n_eq_nogrid,
        )
        out = []
        for r in rows:
            v = r["voxel"]
            out.append(
                dict(
                    voxel=v,
                    reseed=r["reseed"],
                    level=float(T.level[v]),
                    grain=int(T.gmap[v]),
                    dist=round(float(geo["dist"][v]), 3),
                    stored_evals=r["stored_evals"],
                    grid=r["grid"],
                    evals=r["evals"],
                    costed=r["costed"],
                    refined=r["refined"],
                    fin_evals=r["fin_evals"],
                    fin_costs=[round(c, 5) for c in r["fin_costs"]],
                    fin_hit=[round(h, 4) for h in r["fin_hit"]],
                    fin_err=[round(e, 4) for e in r["fin_err"]],
                    taken=r["taken"],
                    n_cand=r["n_cand"],
                )
            )
        RC.write_gz(out, SLIM / f"steps_{RC.SAMPLE}_{name}.json.gz")
    path = SLIM / f"check_{RC.SAMPLE}.json"
    old = json.loads(path.read_text()) if path.exists() else {}
    old.update(chk)
    path.write_text(json.dumps(old, indent=1) + "\n")
    print(json.dumps(chk, indent=1))


def logged(T: Any, geo: Dict[str, Any]) -> None:
    d = RC.CACHE / "S5_big_logged"
    z = np.load(d / "fits.npz")
    meta, R = z["meta"], z["R"].astype(np.float64)
    clusters = RC.read_gz(d / "clusters.json.gz")
    floods = RC.read_gz(d / "floods.json.gz")
    run = json.loads((d / "run.json").read_text())
    stored = RC.stored_run("S5")
    gm = np.asarray(T.gmap)
    n = len(meta)
    v = meta[:, 0].astype(int)
    flood = meta[:, 1].astype(int)
    expd = meta[:, 2].astype(int)
    phase = meta[:, 3].astype(int)
    hit, best_b, ev = meta[:, 4], meta[:, 5], meta[:, 6]
    best_a = np.maximum(hit, best_b)
    ratio = np.where(best_a > 0, hit / np.where(best_a > 0, best_a, 1), 0.0)
    accepted = (best_a > 0) & (ratio > 0.9)
    Rs = R[:, :9].reshape(-1, 3, 3)
    Re = R[:, 9:].reshape(-1, 3, 3)
    err_end = T.err(Re, v)
    err_start = T.err(Rs, v)
    same = np.where(expd >= 0, gm[np.clip(expd, 0, None)] == gm[v], False)
    # the hit ratio of the expander's own accepted fit (or, for the flood's seed, of its seed step)
    seeds = RC.read_gz(d / "seeds.json.gz")
    seed_hit = {s["flood"]: s["fin_hit"][s["refined"].index(s["taken"])] for s in seeds}
    seed_vox = {s["flood"]: s["voxel"] for s in seeds}
    last_hit: Dict[int, float] = {}
    exp_hit = np.full(n, np.nan)
    seen_flood = set()
    for i in range(n):
        f = int(flood[i])
        if f not in seen_flood:
            seen_flood.add(f)
            last_hit[seed_vox[f]] = seed_hit[f]
        if expd[i] >= 0:
            exp_hit[i] = last_hit.get(int(expd[i]), np.nan)
        last_hit[int(v[i])] = float(hit[i])

    # ---- the stored run against the logged run
    st_n = run["stats"]
    ref = stored["stats"]
    keys = [
        "n_seeds",
        "n_reseeds",
        "n_evals_seed",
        "n_evals_reseed",
        "n_evals_neighbor",
        "n_evals_revisit",
        "n_neighbor_fits",
        "n_revisit_attempted",
        "index_evals",
    ]
    same_keys = {k: [st_n[k], ref[k]] for k in keys}
    res_same = [
        (a["voxel"], a["cluster_size"], a["full_search_evaluations"], a["accepted"])
        == (b["voxel"], b["cluster_size"], b["full_search_evaluations"], b["accepted"])
        for a, b in zip(st_n["reseeds"], ref["reseeds"])
    ]
    chk = dict(
        counters_logged_vs_stored=same_keys,
        counters_equal=all(a == b for a, b in same_keys.values()),
        reseeds_logged=len(st_n["reseeds"]),
        reseeds_stored=len(ref["reseeds"]),
        reseeds_equal=bool(all(res_same) and len(st_n["reseeds"]) == len(ref["reseeds"])),
        evals_total_logged=run["evals_total"],
        evals_total_stored=stored["evals_total"],
        n_fit_rows=int(n),
        fit_rows_evals=int(ev.sum()),
        wall_logged_s=run["wall_total_s"],
        timing_label=run["timing_label"],
    )
    (SLIM / "check_logged_big.json").write_text(json.dumps(chk, indent=1) + "\n")
    print(json.dumps(chk, indent=1)[:1500])

    # ---- rejected fits (slim)
    rej = ~accepted
    idx = np.nonzero(rej)[0]
    rej_rows = [
        dict(
            row=int(i),
            voxel=int(v[i]),
            flood=int(flood[i]),
            phase=int(phase[i]),
            expander=int(expd[i]),
            level=float(T.level[v[i]]),
            grain=int(gm[v[i]]),
            dist=round(float(geo["dist"][v[i]]), 3),
            hit=round(float(hit[i]), 4),
            best=round(float(best_a[i]), 4),
            evals=int(ev[i]),
            err_end=round(float(err_end[i]), 4),
            err_start=round(float(err_start[i]), 4),
            same_grain=bool(same[i]),
            exp_hit=None if np.isnan(exp_hit[i]) else round(float(exp_hit[i]), 4),
        )
        for i in idx
    ]
    RC.write_gz(rej_rows, SLIM / "rejected_big.json.gz")

    # ---- clusters and their members
    reseed_floods = [f for f in floods if f["reseed"]]
    assert len(reseed_floods) == len(clusters) == len(ref["reseeds"]), (
        len(reseed_floods),
        len(clusters),
        len(ref["reseeds"]),
    )
    first_row = {}
    for i in range(n):
        first_row.setdefault(int(flood[i]), i)
    c_out: List[Dict[str, Any]] = []
    m_out: List[Dict[str, Any]] = []
    for k, (c, f) in enumerate(zip(clusters, reseed_floods)):
        r = ref["reseeds"][k]
        assert r["voxel"] == f["voxel"] and r["cluster_size"] == len(c["members"])
        fid = f["flood"]
        lim = first_row.get(fid, n)
        in_flood = np.nonzero(flood == fid)[0]
        # the flood's own accepted / rejected fits and evaluations
        c_out.append(
            dict(
                k=k,
                voxel=int(f["voxel"]),
                size=len(c["members"]),
                level=float(T.level[f["voxel"]]),
                grain=int(gm[f["voxel"]]),
                dist=round(float(geo["dist"][f["voxel"]]), 3),
                still_refit=int(r["cluster_voxels_still_refit"]),
                settled=len(c["members"]) - int(r["cluster_voxels_still_refit"]),
                seed_evals=int(r["full_search_evaluations"]),
                reseed_evals=int(r["evaluations"]),
                flood_fits=int(len(in_flood)),
                flood_accepted=int(accepted[in_flood].sum()),
                flood_evals=int(ev[in_flood].sum()),
                n_grains=int(len(set(int(gm[m]) for m in c["members"]))),
                members_own_grain=int(sum(1 for m in c["members"] if gm[m] == gm[f["voxel"]])),
            )
        )
        by_v: Dict[int, int] = {}
        for i in range(lim):  # last fit row of each voxel before this reseed
            by_v[int(v[i])] = i
        for m in c["members"]:
            i = by_v.get(int(m))
            row = dict(k=k, voxel=int(m), level=float(T.level[m]), grain=int(gm[m]))
            row["dist"] = round(float(geo["dist"][m]), 3)
            if i is None:
                row.update(cls="N")
            else:
                right = err_end[i] <= RC.RIGHT
                if right:
                    cls = "A"
                elif same[i]:
                    cls = "B"
                else:
                    cls = "C"
                sa = float(err_start[i])
                row.update(
                    cls=cls,
                    hit=round(float(hit[i]), 4),
                    best=round(float(best_a[i]), 4),
                    err_end=round(float(err_end[i]), 4),
                    err_start=round(sa, 4),
                    same_grain=bool(same[i]),
                    accepted=bool(accepted[i]),
                    phase=int(phase[i]),
                    evals_last=int(ev[i]),
                    row=int(i),
                )
            m_out.append(row)
    RC.write_gz(c_out, SLIM / "clusters_big.json.gz")
    RC.write_gz(m_out, SLIM / "members_big.json.gz")
    # all fits by level: counts only
    lv = np.array([T.level[x] for x in v])
    sm = {
        str(L): dict(
            fits=int((lv == L).sum()),
            rejected=int(((lv == L) & rej).sum()),
            evals=int(ev[lv == L].sum()),
        )
        for L in sorted(set(lv.tolist()))
    }
    (SLIM / "fits_by_level_big.json").write_text(json.dumps(sm, indent=1) + "\n")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--steps", nargs="*", default=[])
    ap.add_argument("--logged", action="store_true")
    a = ap.parse_args()
    T = G.Truth()
    geo = geometry(T)
    if a.steps:
        steps(a.steps, T, geo)
    if a.logged:
        assert RC.SAMPLE == "big"
        logged(T, geo)


if __name__ == "__main__":
    main()
