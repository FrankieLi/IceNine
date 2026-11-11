"""AR3 step 10: compare a library run (D4L, S5L, S8L) with the stored run of the same name.

  libkeys_check.py --kind D4L [--cache DIR]

Exact checks (every number equal; the first difference is reported): the counters of the BFS
before the pass (every numeric counter except wall times), the written-down numbers of the plan
(seeds, reseeds, grid fallbacks, evaluations), the per-voxel records (all fields except
``wall_s``), the map before the pass (rebuilt from ``bfs.boundary_log`` and written as a ``.mic``)
against the stored run's ``recon.mic.gz`` byte for byte, the wrapper counts (S5L) and the early
stops (S8L). Prediction, not guaranteed: the pass in the run gives the stored pass's tries and
evaluations (the script read a ``.mic`` whose angles are rounded; the library reads the map in
memory); a difference is reported with the voxels it touches. Writes
``benchmarks/autoresearch/ar3/keys/check_<kind>.json``.
"""

import argparse
import json
from pathlib import Path
from typing import Any, Dict, List

import libkeys_common as K
import numpy as np


def first_diff_lines(a: bytes, b: bytes) -> Dict[str, Any]:
    la, lb = a.split(b"\n"), b.split(b"\n")
    if len(la) != len(lb):
        return dict(equal=False, lines=(len(la), len(lb)))
    for i, (x, y) in enumerate(zip(la, lb)):
        if x != y:
            return dict(equal=False, first_line=i, new=x.decode()[:160], stored=y.decode()[:160])
    return dict(equal=True, bytes=len(a))


def compare_counters(new: Dict[str, Any], old: Dict[str, Any]) -> Dict[str, Any]:
    a, b = K.counter_stats(new), K.counter_stats(old)
    diff = {k: (a.get(k), b.get(k)) for k in sorted(set(a) | set(b)) if k in b and a.get(k) != b[k]}
    return dict(
        n_stored_counters=len(b),
        differ=diff,
        only_in_stored=sorted(set(b) - set(a)),
        only_in_library=sorted(set(a) - set(b)),
    )


def compare_records(new: Dict[str, Any], old: Dict[str, Any]) -> Dict[str, Any]:
    if set(new) != set(old):
        return dict(equal=False, voxels=(len(new), len(old)))
    for k in sorted(new, key=int):
        a = {f: v for f, v in new[k].items() if f != "wall_s"}
        b = {f: v for f, v in old[k].items() if f != "wall_s"}
        if a != b:
            f = next(f for f in a if a[f] != b.get(f))
            return dict(equal=False, first_voxel=int(k), field=f, new=a[f], stored=b.get(f))
    return dict(equal=True, n=len(new))


def compare_pass(new_npz: Path, old_npz: Path) -> Dict[str, Any]:
    a, b = np.load(new_npz), np.load(old_npz)
    out: Dict[str, Any] = dict(
        tries=(int(len(a["voxel"])), int(len(b["voxel"]))),
        evals=(int(a["evals"].sum()), int(b["evals"].sum())),
        boundary_voxels=(int(len(a["boundary"])) if "boundary" in a else None,
                         int(len(np.unique(b["voxel"])))),
    )  # fmt: skip
    if len(a["voxel"]) == len(b["voxel"]):
        n = len(a["voxel"])
        same_key = (a["voxel"] == b["voxel"]) & (a["cand"] == b["cand"])
        same_end = (a["end"] == b["end"]).reshape(n, -1).all(1)
        end_diff = np.abs(a["end"].astype(np.float64) - b["end"].astype(np.float64)).reshape(n, -1)
        same = same_key & same_end & (a["cost"] == b["cost"]) & (a["evals"] == b["evals"])
        out["same_voxel_and_candidate"] = int(same_key.sum())
        out["same_evaluations"] = int((a["evals"] == b["evals"]).sum())
        out["same_cost_exact"] = int((a["cost"] == b["cost"]).sum())
        out["same_hit_exact"] = int((a["hit"] == b["hit"]).sum())
        out["same_end_exact"] = int(same_end.sum())
        out["max_abs_end_difference"] = float(end_diff.max())
        out["ends_further_apart_than_1e-5"] = int((end_diff.max(1) > 1e-5).sum())
        out["max_abs_cost_difference"] = float(np.abs(a["cost"] - b["cost"]).max())
        out["tries_equal_in_every_field"] = int(same.sum())
        bad = np.where(~same)[0]
        out["tries_differ"] = int(len(bad))
        out["voxels_touched_by_differing_tries"] = sorted({int(v) for v in a["voxel"][bad]})[:30]
        if len(bad):
            i = int(bad[0])
            out["first_differing_try"] = dict(
                voxel=int(a["voxel"][i]), cand=int(a["cand"][i]),
                evals=(int(a["evals"][i]), int(b["evals"][i])),
                cost=(float(a["cost"][i]), float(b["cost"][i])),
                max_abs_end_diff=float(end_diff[i].max()),
            )  # fmt: skip
    else:
        sa = set(zip(a["voxel"].tolist(), a["cand"].tolist()))
        sb = set(zip(b["voxel"].tolist(), b["cand"].tolist()))
        out["pairs_only_in_library"] = sorted(sa - sb)[:50]
        out["pairs_only_in_stored"] = sorted(sb - sa)[:50]
        out["n_pairs_only_in_library"] = len(sa - sb)
        out["n_pairs_only_in_stored"] = len(sb - sa)
        common = sorted(sa & sb)
        ia = {p: i for i, p in enumerate(zip(a["voxel"].tolist(), a["cand"].tolist()))}
        ib = {p: i for i, p in enumerate(zip(b["voxel"].tolist(), b["cand"].tolist()))}
        ja = np.array([ia[p] for p in common])
        jb = np.array([ib[p] for p in common])
        e_a, e_b = a["end"][ja].astype(np.float64), b["end"][jb].astype(np.float64)
        out["same_voxel_and_candidate"] = len(common)
        out["same_end_exact"] = int((e_a == e_b).reshape(len(common), -1).all(1).sum())
        out["same_evaluations"] = int((a["evals"][ja] == b["evals"][jb]).sum())
        out["voxels_with_a_pair_in_one_list_only"] = sorted({v for v, _ in (sa ^ sb)})
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--kind", required=True, choices=sorted(K.KEYS))
    ap.add_argument("--cache", default=str(K.CACHE))
    a = ap.parse_args()
    kind = a.kind
    run = Path(a.cache) / kind
    new = json.loads((run / "run.json").read_text())
    stored_dir = K.STORED[kind]
    old = json.loads((stored_dir / "run.json").read_text())
    exp = K.EXPECT[kind]
    st = new["stats_before_pass"]
    # ``stored_run`` is the "name" field of the stored run.json. The stored S8 run.json carries
    # "S5": reseed_s8.py calls run_seed.run("S5", ...) with the early stop patched in. Its
    # s8.json reads "S8 = S5 + E1" and its totals are S8's; check_S8L.json is not a comparison
    # with S5.
    res: Dict[str, Any] = dict(kind=kind, stored_run=old.get("name"), keys=new["keys"])
    res["plan_numbers"] = dict(
        seeds=(st["n_seeds"], exp["seeds"]),
        reseeds=(st["n_reseeds"], exp["reseeds"]),
        grid=(
            st.get("index_fallback_seed_test", 0) + st.get("index_fallback_no_candidate", 0),
            exp["grid"],
        ),
        evals_before_pass=(new["evals_total_before_pass"], exp["evals"]),
    )
    # the final stats: the pass adds its own keys (listed as "only in library") and changes no other
    res["counters"] = compare_counters(new["stats"], old["stats"])
    res["stored_evals_total"] = old.get("evals_total")
    res["records"] = compare_records(
        K.read_gz(run / "records.json.gz"), K.read_gz(stored_dir / "records.json.gz")
    )
    res["map_before_pass_vs_stored_recon_mic"] = first_diff_lines(
        K.gz_bytes(run / "recon_before_pass.mic.gz"), K.gz_bytes(stored_dir / "recon.mic.gz")
    )
    if kind == "S5L":
        w = old["wrapper"]
        fits = st["n_neighbor_fits"] + st["n_revisit_attempted"] + st["n_restart_attempted"]
        res["wrapper_counts"] = dict(
            fits=(fits, w["wrapper_calls"]),
            mc_fallbacks=(st["cma_first_fallback_runs"], w["fallback_runs"]),
            fallback_won=(st["cma_first_fallback_won"], w["fallback_won"]),
        )
    if kind == "D4L":
        w = old["wrapper"]
        fits = st["n_neighbor_fits"] + st["n_revisit_attempted"] + st["n_restart_attempted"]
        res["wrapper_counts"] = dict(
            fits=(fits, w["wrapper_calls"]),
            mc_fallbacks=(st["cma_first_fallback_runs"], w["fallback_runs"]),
            fallback_won=(st["cma_first_fallback_won"], w["fallback_won"]),
        )
    if kind == "S8L":
        s8 = json.loads((stored_dir / "s8.json").read_text())
        res["early_stops"] = dict(
            stopped_early=(st["index_stopped_early"], s8["stopped_early"]),
            seed_steps=(st["index_runs"], s8["seed_steps"]),
            seed_step_evals=(st["index_evals"], old["stats"]["index_evals"]),
        )
    pdir = K.STORED_PASS[kind][0]
    res["pass_in_run"] = compare_pass(run / "pass_tries.npz", pdir / "tries.npz")
    res["pass_counts"] = dict(
        tries=(new["stats"]["boundary_tries"], exp["pass_tries"]),
        evals=(new["stats"]["n_evals_boundary_pass"], exp["pass_evals"]),
        boundary_voxels=new["stats"]["boundary_voxels"],
        changed=new["stats"]["boundary_changed"],
    )
    res["timing_label"] = new["timing_label"]
    res["wall_total_s"] = new["wall_total_s"]
    out = K.OUT / f"check_{kind}.json"
    out.write_text(json.dumps(res, indent=1, default=str) + "\n")
    print(json.dumps(res, indent=1, default=str)[:6000])


if __name__ == "__main__":
    main()
