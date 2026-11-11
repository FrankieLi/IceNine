"""AR5 step 1, part 3: slim results of the vote count (H3) and of check V1 (``mid``).

  uv run python scripts/autoresearch/ar5/stop_summary_part3.py

Reads ``cache/votes_<sample>.json`` and ``cache/v1_mid_shard*of*.json`` (or ``cache/v1_mid.json``);
writes ``part3_votes.npz`` (voxel, flagged, wrong, max_votes per sample), ``part3_votes.json``
(the counts at the cut that separates best, labelled as chosen after the result) and
``part3_v1.json`` (V1 rows and the pass rule). Intervals are Wilson 95% and treat voxels as
independent although they cluster in grains.
"""

import glob
import json
from typing import Any, Dict, List

import numpy as np
import stop_common as SC
import stop_rules as SR

AR5 = SC.AR5
SAMPLES = ("big", "grad", "mid")


def votes() -> Dict[str, Any]:
    arrays: Dict[str, np.ndarray] = {}
    per: Dict[str, Any] = {}
    allw: List[np.ndarray] = []
    allr: List[np.ndarray] = []
    un_all: List[np.ndarray] = []
    for s in SAMPLES:
        d = json.loads((AR5 / "cache" / f"votes_{s}.json").read_text())
        rows = d["rows"]
        v = np.array([r["voxel"] for r in rows], np.int32)
        fl = np.array([r["flagged"] for r in rows], bool)
        wr = np.array([r["wrong"] for r in rows], bool)
        mv = np.array([r["max_votes"] for r in rows], np.int32)
        arrays[f"{s}_voxel"], arrays[f"{s}_flagged"] = v, fl
        arrays[f"{s}_wrong"], arrays[f"{s}_votes"] = wr, mv
        w, r_ = mv[fl & wr], mv[fl & ~wr]
        allw.append(w)
        allr.append(r_)
        un_all.append(mv[~fl])
        per[s] = group(w, r_, mv[~fl], wr[~fl])
    per["all"] = group(
        np.concatenate(allw), np.concatenate(allr), np.concatenate(un_all), np.zeros(0, bool)
    )
    np.savez_compressed(AR5 / "part3_votes.npz", **arrays)
    out = dict(
        groups=per,
        note="the cut is chosen after the result (it is the cut with the largest share of wrong "
        "voxels at or above it minus the share of right voxels at or above it)",
        separates_rule="Wilson 95% intervals of the share at or above the cut do not overlap",
    )
    (AR5 / "part3_votes.json").write_text(json.dumps(out, indent=1, default=float) + "\n")
    return out


def group(w: np.ndarray, r: np.ndarray, un: np.ndarray, un_wrong: np.ndarray) -> Dict[str, Any]:
    cut = SR.best_cut(w, r) if len(w) and len(r) else None
    out: Dict[str, Any] = dict(
        n_wrong_flagged=int(len(w)),
        n_right_flagged=int(len(r)),
        median_wrong=float(np.median(w)) if len(w) else None,
        median_right=float(np.median(r)) if len(r) else None,
        q25_q75_wrong=[float(x) for x in np.percentile(w, [25, 75])] if len(w) else None,
        q25_q75_right=[float(x) for x in np.percentile(r, [25, 75])] if len(r) else None,
        n_unflagged_sampled=int(len(un)),
        median_unflagged_sampled=float(np.median(un)) if len(un) else None,
        wrong_among_unflagged_sampled=int(un_wrong.sum()),
        best_cut=cut,
    )
    if cut is not None:
        out["wrong_interval"] = SR.wilson_row(cut["wrong_ge"], cut["n_wrong"])
        out["right_interval"] = SR.wilson_row(cut["right_ge"], cut["n_right"])
        out["separates"] = SR.separates(
            cut["wrong_ge"], cut["n_wrong"], cut["right_ge"], cut["n_right"]
        )
    return out


def v1() -> Dict[str, Any]:
    files = sorted(glob.glob(str(AR5 / "cache" / "v1_mid_shard*of*.json")))
    files += [f for f in [str(AR5 / "cache" / "v1_mid.json")] if glob.glob(f)]
    rows: List[Dict[str, Any]] = []
    caps = []
    walls = []
    n_shard_flagged = 0
    for f in files:
        d = json.loads(open(f).read())
        rows += d["rows"]
        caps.append(bool(d["capped"]))
        walls.append(d["wall_s"])
        n_shard_flagged += d["n_flagged_in_shard"]
    flags = json.loads((AR5 / "part3_flags.json").read_text())["at_threshold"]["mid"]
    n_flagged = int(flags["flagged"])
    n_wrong_flagged = int(flags["wrong_flagged"])
    run = json.loads((SC.CACHE / "stored" / "mid" / "S8_mid" / "run.json").read_text())
    ps = json.loads(
        (SC.CACHE / "stored" / "mid" / "passes_mid" / "D4_S8" / "pass.json").read_text()
    )
    total = int(run["evals_total"]) + int(ps["evals"])
    n = len(rows)
    ev = int(sum(r["evals"] for r in rows))
    fixed = sum(1 for r in rows if r["wrong_after_pass"] and not r["wrong_with_v1"])
    broken = sum(1 for r in rows if not r["wrong_after_pass"] and r["wrong_with_v1"])
    wrong_proc = sum(1 for r in rows if r["wrong_after_pass"])
    mean = ev / n if n else float("nan")
    est = mean * n_flagged
    out = dict(
        threshold=json.loads((AR5 / "part3_pick.json").read_text())["threshold"],
        flagged_total=n_flagged,
        wrong_flagged_total=n_wrong_flagged,
        voxels_processed=n,
        shards=len(files),
        capped_shards=int(sum(caps)),
        all_flagged_processed=bool(n == n_flagged),
        evals_processed=ev,
        mean_evals_per_voxel=mean,
        evals_whole_flag_estimate=est,
        run_total_evals=total,
        share_of_run_estimate=est / total,
        share_of_run_processed=ev / total,
        kept=sum(1 for r in rows if r["kept"]),
        wrong_among_processed=wrong_proc,
        fixed=fixed,
        broken=broken,
        wrong_before_v1=wrong_proc,
        wrong_after_v1=wrong_proc - fixed + broken,
        fixed_share_of_wrong_processed=fixed / wrong_proc if wrong_proc else None,
        fixed_interval=SR.wilson_row(fixed, wrong_proc) if wrong_proc else None,
        pass_fixes_more_than_breaks=bool(fixed > broken),
        pass_fixes_quarter=bool(wrong_proc and fixed >= 0.25 * wrong_proc),
        pass_cost_under_5pc=bool(est / total < 0.05),
        wall_s_per_shard=walls,
        timing_label="contended (another agent on the machine); not a timing result",
    )
    out["passes"] = bool(
        out["pass_fixes_more_than_breaks"]
        and out["pass_fixes_quarter"]
        and out["pass_cost_under_5pc"]
    )
    (AR5 / "part3_v1.json").write_text(
        json.dumps(dict(summary=out, rows=rows), default=float) + "\n"
    )
    return out


def main() -> None:
    vt = votes()
    print(json.dumps(vt["groups"]["all"], indent=1, default=float))
    print(json.dumps(v1(), indent=1, default=float))


if __name__ == "__main__":
    main()
