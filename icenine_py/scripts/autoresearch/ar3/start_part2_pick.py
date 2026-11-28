"""AR3 step 9, part 2: the pick (the pick set only), written to ``start/f_pick.json``.

  uv run python scripts/autoresearch/ar3/start_part2_pick.py

The split is step 8's: true grains with an even index are the pick set, with an odd index the test
set. This script reads the fits of the pick set only (the grain index of each fit's voxel is in
the feature files) and the F2 files of the pick set (``part2_pick_<sample>_f2.json``, made by
``start_part2.py --set pick``). It writes the tables of every setting on the pick set, pooled over
the two samples with a logged run of ``S5`` (``big`` and ``grad``), and the pick.

F1 (the start test): at the start of each neighbour fit and revisit one evaluation of the start
(counted); if its hit ratio is at least the threshold (0.90, 0.95 or 0.98) the start is kept and
the fit is skipped. Per fit, the stored evaluations against 1 + (0 if kept, else the stored
evaluations). A case is *fixed* if the stored fit ended wrong (error above 1 degree) and the kept
start is right, *broken* if the stored fit ended right and the kept start is wrong. The *target*
is the evaluations of the neighbour fits and revisits.

F2 (one more fit from the end of a rejected fit): see ``start_part2.py``. The saving is the
seed-step evaluations of the reseed clusters that vanish, less F2's own evaluations on every
rejected fit; the target is the seed-step evaluations of the reseed clusters.

The pick rule (the plan): among a remedy's settings, the one that saves most evaluations on the
pick set among those that break no more cases than they fix; a setting must also save more than
zero evaluations (written here, before any test result). If none qualifies, the remedy fails.
"""

import json
from pathlib import Path
from typing import Any, Dict, List

import numpy as np

HERE = Path(__file__).resolve().parent
START = HERE.parents[2] / "benchmarks" / "autoresearch" / "ar3" / "start"
THRESH = (0.90, 0.95, 0.98)
SIGMAS = (0.05, 0.1, 0.2)
SAMPLES = ("big", "grad")
RIGHT = 1.0


def load_feat(sample: str) -> Dict[str, np.ndarray]:
    return dict(np.load(START / f"feat_{sample}.npz"))


def f1_table(f: Dict[str, np.ndarray], mask: np.ndarray) -> Dict[str, Any]:
    """F1 on the fits (neighbour and revisit) selected by ``mask``, one row per threshold."""
    sel = mask & ((f["phase"] == 0) | (f["phase"] == 1))
    ev = f["evals"][sel].astype(np.int64)
    sh = f["start_hit"][sel]
    best = f["best_before"][sel]
    r_stored = f["err_end"][sel] <= RIGHT
    r_start = f["err_start"][sel] <= RIGHT
    base = int(ev.sum())
    out: Dict[str, Any] = dict(fits=int(sel.sum()), evals_stored=base)
    for t in THRESH:
        kept = sh >= t
        ev_f1 = int((1 + np.where(kept, 0, ev)).sum())
        ba = np.maximum(sh, best)
        kept_fail_rel = kept & ~((ba > 0) & (sh / np.where(ba > 0, ba, 1) > 0.9))
        out[str(t)] = dict(
            kept=int(kept.sum()),
            evals_f1=ev_f1,
            cut=float(1 - ev_f1 / base) if base else None,
            fixed=int((kept & ~r_stored & r_start).sum()),
            broken=int((kept & r_stored & ~r_start).sum()),
            both_wrong=int((kept & ~r_stored & ~r_start).sum()),
            kept_fail_relative_test=int(kept_fail_rel.sum()),
        )
    return out


def pool_f1(tables: List[Dict[str, Any]]) -> Dict[str, Any]:
    out: Dict[str, Any] = dict(
        fits=sum(t["fits"] for t in tables), evals_stored=sum(t["evals_stored"] for t in tables)
    )
    for t in THRESH:
        rows = [x[str(t)] for x in tables]
        ev = sum(r["evals_f1"] for r in rows)
        out[str(t)] = dict(
            kept=sum(r["kept"] for r in rows),
            evals_f1=ev,
            cut=1 - ev / out["evals_stored"] if out["evals_stored"] else None,
            fixed=sum(r["fixed"] for r in rows),
            broken=sum(r["broken"] for r in rows),
            both_wrong=sum(r["both_wrong"] for r in rows),
            kept_fail_relative_test=sum(r["kept_fail_relative_test"] for r in rows),
        )
    return out


def pool_f2(files: List[Dict[str, Any]]) -> Dict[str, Any]:
    base = sum(x["seed_evals_clusters"] for x in files)
    out: Dict[str, Any] = dict(
        rejected_fits=sum(x["rejected_fits"] for x in files),
        clusters=sum(x["clusters"] for x in files),
        seed_evals_clusters=base,
    )
    for s in SIGMAS:
        rows = [x["settings"][str(s)] for x in files]
        saving = sum(r["seed_evals_vanish"] for r in rows) - sum(r["f2_evals"] for r in rows)
        out[str(s)] = dict(
            f2_evals=sum(r["f2_evals"] for r in rows),
            rejected_pass=sum(r["rejected_pass"] for r in rows),
            clusters_vanish=sum(r["clusters_vanish"] for r in rows),
            seed_evals_vanish=sum(r["seed_evals_vanish"] for r in rows),
            net_saving=saving,
            cut=saving / base if base else None,
            fixed=sum(r["fixed"] for r in rows),
            broken=sum(r["broken"] for r in rows),
            both_wrong=sum(r["both_wrong"] for r in rows),
        )
    return out


def choose(tab: Dict[str, Any], keys: List[str], gain: str) -> Any:
    """The setting with the largest ``gain`` among those with broken <= fixed and gain > 0."""
    ok = [k for k in keys if tab[k]["broken"] <= tab[k]["fixed"] and tab[k][gain] > 0]
    return max(ok, key=lambda k: tab[k][gain]) if ok else None


def main() -> None:
    f1_by_sample = {}
    for s in SAMPLES:
        f = load_feat(s)
        f1_by_sample[s] = f1_table(f, (f["grain"] % 2) == 0)  # the pick set only
    f1 = pool_f1(list(f1_by_sample.values()))
    f1_pick = choose(f1, [str(t) for t in THRESH], "cut")
    out: Dict[str, Any] = dict(
        sets="pick = even true-grain index (this file); test = odd (not read here)",
        F1=dict(
            by_sample=f1_by_sample,
            pooled=f1,
            chosen_threshold=None if f1_pick is None else float(f1_pick),
        ),
    )
    f2_files = []
    for s in SAMPLES:
        p = START / f"part2_pick_{s}_f2.json"
        if p.exists():
            f2_files.append(json.loads(p.read_text()))
    if len(f2_files) == len(SAMPLES):
        f2 = pool_f2(f2_files)
        f2_pick = choose(f2, [str(s) for s in SIGMAS], "net_saving")
        out["F2"] = dict(pooled=f2, chosen_sigma_deg=None if f2_pick is None else float(f2_pick))
    else:
        out["F2"] = dict(note="F2 pick files missing")
    (START / "f_pick.json").write_text(json.dumps(out, indent=1) + "\n")
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
