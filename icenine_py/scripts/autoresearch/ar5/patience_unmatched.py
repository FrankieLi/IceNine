"""AR5 step 11, review fix 3: the seed steps that only one arm has (unmatched steps).

  RESEED_SAMPLE=<held|mid|grad|big> [PATIENCE_CACHE=<cache>] \
      uv run python scripts/autoresearch/ar5/patience_unmatched.py

``patience_stats.py`` compares seed steps matched by (voxel, occurrence at that voxel) only. This
script reads the same stored files (read only: both arms' ``e1_log.json.gz`` and the per-voxel
error maps ``stats/part2_<arm>_<sample>_err.npz`` that ``patience_stats.py`` left in the cache)
and lists the steps of each arm that the other arm does not have. For each one it gives the taken
end (``order_logs.taken_end``: hit ratio, cost, error to the truth, whether it passes the seed test
0.80), the other arm's map error at that voxel, before and after the boundary pass, and (S only)
whether any MC refinement of the step was cut (``mc_traces.json.gz``).

The broad reading of rule (4) counts, among S's unmatched steps, a taken end that is wrong (over 1
degree) where G's map before the pass is right at that voxel. The frozen rule (4) of the plan is
applied to matched steps (``stats_<sample>.json``); this file adds the unmatched part, for the
record. Writes ``patience/unmatched_<sample>.json``; writes nothing into the cache. Truth is used
for scoring only; no p-values.
"""

import gzip
import json
from typing import Any, Dict, List

import patience_common as PC  # noqa: E402  (first)

import fit_common as FC  # noqa: E402
import numpy as np  # noqa: E402
import order_logs as OL  # noqa: E402

import gradient_common as G  # noqa: E402

SAMPLE = FC.SAMPLE
SEED_TEST = 0.80


def read_gz(p: Any) -> Any:
    with gzip.open(p, "rt") as f:
        return json.load(f)


def side(
    log_own: List[Dict[str, Any]],
    keys: List[Any],
    idx: Dict[Any, Dict[str, Any]],
    err_other: Dict[str, np.ndarray],
    T: Any,
) -> List[Dict[str, Any]]:
    """One row per unmatched step of one arm: its taken end and the other arm's map there."""
    rows: List[Dict[str, Any]] = []
    for k in keys:
        e = idx[k]
        v = int(k[0])
        row: Dict[str, Any] = dict(voxel=v, occurrence=int(k[1]), reseed=bool(e.get("reseed")))
        if not e["fits"]:
            row.update(has_fits=False)
            rows.append(row)
            continue
        t = OL.taken_end(e)
        R = np.asarray(t["R"], np.float64).reshape(1, 3, 3)
        err = float(T.err(R, [v])[0])
        row.update(
            has_fits=True,
            hit=float(t["hit"]),
            cost=float(t["cost"]),
            passes_seed_test=bool(t["hit"] >= SEED_TEST),
            err_deg=float(err),
            wrong=bool(err > G.WRONG),
            other_map_err_before_deg=float(err_other["before"][v]),
            other_map_err_after_deg=float(err_other["after"][v]),
            other_map_wrong_before=bool(err_other["before"][v] > G.WRONG),
            other_map_wrong_after=bool(err_other["after"][v] > G.WRONG),
        )
        rows.append(row)
    return rows


def main() -> None:
    gd = PC.PCACHE / f"G_{SAMPLE}"
    sdir = PC.PCACHE / f"S_{SAMPLE}"
    sd = PC.PCACHE / "stats"
    lg = read_gz(gd / "e1_log.json.gz")
    ls = read_gz(sdir / "e1_log.json.gz")
    ig, isx = OL.index(lg), OL.index(ls)
    e_g = dict(np.load(sd / f"part2_G_{SAMPLE}_err.npz"))
    e_s = dict(np.load(sd / f"part2_S_{SAMPLE}_err.npz"))
    T = G.Truth()
    only_s = sorted(set(isx) - set(ig))
    only_g = sorted(set(ig) - set(isx))
    rs = side(ls, only_s, isx, e_g, T)
    rg = side(lg, only_g, ig, e_s, T)
    # was any MC refinement of an S-only step cut by P(3000)? (trace "step" is the 1-based row)
    ts = read_gz(sdir / "mc_traces.json.gz")
    cut_rows = {int(t["step"]) - 1 for t in ts if t.get("cut")}
    row_of = {id(e): i for i, e in enumerate(ls)}
    for k, r in zip(only_s, rs):
        r["cut_in_step"] = bool(row_of[id(isx[k])] in cut_rows)
    s_wrong_g_right = [r for r in rs if r.get("wrong") and not r["other_map_wrong_before"]]
    g_wrong_s_right = [r for r in rg if r.get("wrong") and not r["other_map_wrong_before"]]
    out = dict(
        sample=SAMPLE,
        matched_steps=len(set(ig) & set(isx)),
        seed_steps_g=len(lg),
        seed_steps_s=len(ls),
        only_s_n=len(only_s),
        only_g_n=len(only_g),
        only_s_wrong=sum(1 for r in rs if r.get("wrong")),
        only_g_wrong=sum(1 for r in rg if r.get("wrong")),
        only_s_wrong_where_g_map_right=len(s_wrong_g_right),
        only_g_wrong_where_s_map_right=len(g_wrong_s_right),
        only_s_with_a_cut=sum(1 for r in rs if r.get("cut_in_step")),
        only_s_wrong_where_g_map_right_voxels=[r["voxel"] for r in s_wrong_g_right],
        only_g_wrong_where_s_map_right_voxels=[r["voxel"] for r in g_wrong_s_right],
        note=(
            "unmatched = a seed step (voxel, occurrence) that only one arm has; the other arm's "
            "map is read before the pass for the counts; the frozen rule (4) counts matched steps"
        ),
        only_s=rs,
        only_g=rg,
    )
    (PC.PAT / f"unmatched_{SAMPLE}.json").write_text(
        json.dumps(out, indent=1, default=float) + "\n"
    )
    print(json.dumps({k: v for k, v in out.items() if k not in ("only_s", "only_g")}, indent=1))


if __name__ == "__main__":
    main()
