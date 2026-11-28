"""AR5 step 5: compare a library run (``keys_run.py``) with the stored script run of the same name.

  AR3_SAMPLE=big uv run python scripts/autoresearch/ar5/keys_check.py --kind S11L
  uv run python scripts/autoresearch/ar5/keys_check.py --kind S11L          # first gradient sample

Exact checks (the first difference is reported): the numbers written down in the plan (seeds,
reseeds, evaluations before the pass), every numeric counter of the stored run (wall times left
out), the per-voxel records (all fields except ``wall_s``), the map before the pass (rebuilt from
``bfs.boundary_log`` and written as a ``.mic``) against the stored run's ``recon.mic.gz`` byte
for byte, the counts of the stored wrapper (CMA fits, ``cma_not_fell``), the seed steps (the
stored ``e1_log``: evaluations of each refinement, the early stop) and the mixed refinement (the
stored ``mix_log``: CMA hit ratio and whether MC ran, per candidate, in order). The pass in the run
is compared with the stored pass; a difference is reported with the voxels it touches, not
repaired (step 10 found the pass in a run reads the unrounded map). The truth is read here only
(the library never sees it): wrong voxels (error above 1 degree) before and after the pass, library
against stored, paired. Writes ``benchmarks/autoresearch/ar5/keys/check_<kind>_<sample>.json``.
"""

import argparse
import json
from pathlib import Path
from typing import Any, Dict, List

import keys_common as K
import numpy as np

import gradient_common as G  # noqa: E402


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
    )
    if len(a["voxel"]) == len(b["voxel"]):
        n = len(a["voxel"])
        same_key = (a["voxel"] == b["voxel"]) & (a["cand"] == b["cand"])
        same_end = (a["end"] == b["end"]).reshape(n, -1).all(1)
        end_diff = np.abs(a["end"].astype(np.float64) - b["end"].astype(np.float64)).reshape(n, -1)
        same = same_key & same_end & (a["cost"] == b["cost"]) & (a["evals"] == b["evals"])
        out["same_voxel_and_candidate"] = int(same_key.sum())
        out["same_evaluations"] = int((a["evals"] == b["evals"]).sum())
        out["same_end_exact"] = int(same_end.sum())
        out["max_abs_end_difference"] = float(end_diff.max())
        out["tries_differ"] = int((~same).sum())
        bad = np.where(~same)[0]
        out["voxels_touched_by_differing_tries"] = sorted({int(v) for v in a["voxel"][bad]})[:30]
        if len(bad):
            i = int(bad[0])
            out["first_differing_try"] = dict(
                voxel=int(a["voxel"][i]), cand=int(a["cand"][i]),
                evals=(int(a["evals"][i]), int(b["evals"][i])),
                max_abs_end_diff=float(end_diff[i].max()),
            )  # fmt: skip
    else:
        sa = set(zip(a["voxel"].tolist(), a["cand"].tolist()))
        sb = set(zip(b["voxel"].tolist(), b["cand"].tolist()))
        out["n_pairs_only_in_library"] = len(sa - sb)
        out["n_pairs_only_in_stored"] = len(sb - sa)
        out["voxels_with_a_pair_in_one_list_only"] = sorted({v for v, _ in (sa ^ sb)})[:30]
    return out


def compare_seed_steps(index_log: List[Dict[str, Any]], e1: List[Dict[str, Any]]) -> Dict[str, Any]:
    """Library ``index_log`` against the stored seed-step log (one row per seed step)."""
    out: Dict[str, Any] = dict(steps=(len(index_log), len(e1)))
    if len(index_log) != len(e1):
        return out
    for i, (a, b) in enumerate(zip(index_log, e1)):
        la = (a["voxel"], a["finisher_evals"], a.get("stopped_at", -1), a["n_cand"])
        lb = (b["voxel"], b["fin_evals"], b["stop"], b["n_cand"])
        if la != lb:
            out.update(equal=False, first_step=i, library=la, stored=lb)
            return out
    out["equal"] = True
    return out


def compare_mixed(index_log: List[Dict[str, Any]], mix: List[Dict[str, Any]]) -> Dict[str, Any]:
    lib = [m for e in index_log for m in e.get("mixed", [])]
    out: Dict[str, Any] = dict(candidates=(len(lib), len(mix)))
    out["mc_ran"] = (sum(m["mc_ran"] for m in lib), sum(m["mc_ran"] for m in mix))
    if len(lib) != len(mix):
        return out
    for i, (a, b) in enumerate(zip(lib, mix)):
        if a["cma_hit"] != b["c_hit"] or a["mc_ran"] != b["mc_ran"]:
            out.update(equal=False, first_candidate=i, library=a, stored=b)
            return out
    out["equal"] = True
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--kind", required=True, choices=sorted(K.KEYS))
    ap.add_argument("--cache", default=str(K.CACHE))
    a = ap.parse_args()
    kind, sample = a.kind, K.SAMPLE
    run = Path(a.cache) / f"{kind}_{sample}"
    new = json.loads((run / "run.json").read_text())
    sdir = K.STORED[(kind, sample)]
    old = json.loads((sdir / "run.json").read_text())
    fit = json.loads((sdir / "fit.json").read_text()) if (sdir / "fit.json").exists() else {}
    exp = K.EXPECT.get((kind, sample), {})
    st = new["stats_before_pass"]
    res: Dict[str, Any] = dict(kind=kind, sample=sample, keys=new["keys"], stored_dir=sdir.name)
    stored_evals = old.get("evals_total", K.evals_total(old["stats"]))
    res["plan_numbers"] = dict(
        evals_before_pass=(new["evals_total_before_pass"], exp.get("evals", stored_evals)),
        seeds=(st["n_seeds"], exp.get("seeds", old["stats"]["n_seeds"])),
        reseeds=(st["n_reseeds"], exp.get("reseeds", old["stats"]["n_reseeds"])),
        grid_fallbacks=st.get("index_fallback_seed_test", 0)
        + st.get("index_fallback_no_candidate", 0),
    )  # fmt: skip
    res["stored_evals_total"] = stored_evals
    res["counters"] = compare_counters(new["stats"], old["stats"])
    res["records"] = compare_records(
        K.read_gz(run / "records.json.gz"), K.read_gz(sdir / "records.json.gz")
    )
    res["map_before_pass_vs_stored_recon_mic"] = first_diff_lines(
        K.gz_bytes(run / "recon_before_pass.mic.gz"), K.gz_bytes(sdir / "recon.mic.gz")
    )
    w = old.get("wrapper")
    if w:
        fits = st["n_neighbor_fits"] + st["n_revisit_attempted"] + st["n_restart_attempted"]
        res["wrapper_counts"] = dict(
            fits=(fits, w["wrapper_calls"]),
            cma_not_fell=(st.get("cma_not_fell"), w.get("cma_not_fell")),
            cma_first_fallback_runs=(st.get("cma_first_fallback_runs"), w.get("fallback_runs")),
        )
    index_log = K.read_gz(run / "index_log.json.gz")
    if fit:
        res["seed_steps"] = dict(
            steps=(st["index_runs"], fit["seed_steps"]),
            stopped_early=(st["index_stopped_early"], fit["stopped_early"]),
        )
    if (sdir / "e1_log.json.gz").exists():
        res["seed_step_log"] = compare_seed_steps(index_log, K.read_gz(sdir / "e1_log.json.gz"))
    if (sdir / "mix_log.json").exists():
        res["mixed_refinement"] = compare_mixed(
            index_log, json.loads((sdir / "mix_log.json").read_text())
        )
        res["mixed_counters"] = dict(
            cands=st.get("index_mixed_cands"),
            cma_reached=st.get("index_mixed_cma_reached"),
            mc_ran=st.get("index_mixed_mc_ran"),
        )
    pdir, _before, after_npy = K.STORED_PASS[(kind, sample)]
    res["pass_in_run"] = compare_pass(run / "pass_tries.npz", pdir / "tries.npz")
    s2 = new["stats"]
    res["pass_counts"] = dict(
        tries=s2["boundary_tries"],
        evals=s2["n_evals_boundary_pass"],
        boundary_voxels=s2["boundary_voxels"],
        changed=s2["boundary_changed"],
        evals_total_with_pass=new["evals_total_with_pass"],
        stored_with_pass=exp.get("with_pass"),
    )
    # the truth: wrong voxels (error above 1 degree), library against stored, paired
    T = G.Truth()
    n = T.n
    maps = {
        "lib_before": G.read_map(run / "recon_before_pass.mic.gz", n)[0],
        "lib_after": G.read_map(run / "recon.mic.gz", n)[0],
        "stored_before": G.read_map(sdir / "recon.mic.gz", n)[0],
        "stored_after": np.load(after_npy).astype(np.float64),
    }
    wr = {k: T.err(R) > G.WRONG for k, R in maps.items()}
    res["wrong"] = {k: int(v.sum()) for k, v in wr.items()}
    for tag in ("before", "after"):
        lib, stv = wr[f"lib_{tag}"], wr[f"stored_{tag}"]
        res[f"paired_{tag}"] = dict(
            fixed_vs_stored=int((stv & ~lib).sum()),
            broken_vs_stored=int((~stv & lib).sum()),
            both_wrong=int((stv & lib).sum()),
        )
    if "wrong_before" in exp:
        res["wrong"]["expected_before"] = exp["wrong_before"]
    res["timing_label"] = new["timing_label"]
    res["wall_total_s"] = new["wall_total_s"]
    K.OUT.mkdir(parents=True, exist_ok=True)
    (K.OUT / f"check_{kind}_{sample}.json").write_text(
        json.dumps(res, indent=1, default=str) + "\n"
    )
    print(json.dumps(res, indent=1, default=str)[:7000])


if __name__ == "__main__":
    main()
