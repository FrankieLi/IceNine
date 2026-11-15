"""AR5 step 14, phase A: the failed rescues of arm ``H`` (one case each) and the stored-count check.

Reads, per sample, step 12's linked ``H_<s>/{rescue_log.json.gz, records.json.gz, recon.mic.gz}``
(through ``polish_common.load_gz`` and the guards). A case is a log row with ``attempted`` and not
``rescued``: its rescue ends (the tries of ``G1``, in ``G1``'s order, with their stored hit ratio,
cost and orientation), the grid search's evaluations and the grid end (the voxel's orientation in
the map before the pass, the hit ratio of its stored record). Truth is added by the probe, for
scoring only. The stored counts that the plan's stop rule names are checked here.
"""

from typing import Any, Dict, List

import polish_common as PC  # noqa: E402  (first)

import numpy as np  # noqa: E402

EXPECT: Dict[str, Dict[str, int]] = {
    # H total (with the pass), failed rescues, their grid evaluations (plan, "What the stored
    # record says"); held2's grid evaluations are checked at gate 2, after the freeze
    "held": dict(total=29_588_170, failed=19, grid=6_563_817),
    "mid": dict(total=11_950_397, failed=4, grid=1_334_072),
    "grad": dict(total=10_180_191, failed=5, grid=1_697_698),
    "big": dict(total=3_635_256, failed=0, grid=0),
    "held2": dict(total=31_143_763, failed=21, grid=7_252_187),
}
N_VOXELS = 24570


def cache_dir(sample: str) -> Any:
    return PC.PCACHE / f"H_{sample}"


def failed_cases(sample: str) -> List[Dict[str, Any]]:
    """The failed rescues of ``H_<sample>`` (rescue ends, grid evaluations; no truth)."""
    d = cache_dir(sample)
    log = PC.load_gz(d / "rescue_log.json.gz")
    out = []
    for r in log:
        if r["attempted"] and not r["rescued"]:
            out.append(r)
    return out


def counts(sample: str) -> Dict[str, Any]:
    """The stored counts of ``H_<sample>`` against the plan's: failed rescues, their grid
    evaluations, the rescue evaluations and the run total before the pass."""
    d = cache_dir(sample)
    cs = failed_cases(sample)
    run = PC.load_json(d / "run.json")
    exp = EXPECT[sample]
    got = dict(
        failed=len(cs),
        grid=int(sum(int(c["grid_evals"]) for c in cs)),
        rescue_evals=int(sum(int(c["rescue_evals"]) for c in cs)),
        before_pass=int(run["evals_total"]),
    )
    st = PC.load_json(
        PC.ICE / "benchmarks" / "autoresearch" / "ar5" / "stallhit" / f"stats_{sample}.json"
    )
    got["total_with_pass"] = int(st["arms"]["H"]["stages"]["total"])
    ok = (
        got["failed"] == exp["failed"]
        and got["grid"] == exp["grid"]
        and got["total_with_pass"] == exp["total"]
    )
    return dict(sample=sample, expected=exp, got=got, ok=bool(ok))


def grid_ends(sample: str, voxels: List[int]) -> Dict[int, Dict[str, Any]]:
    """The grid end of each voxel: its orientation in the map before the pass (``recon.mic.gz``;
    a later neighbour fit may have moved it by a hair, as in step 9) and the hit ratio of its
    stored record."""
    import gradient_common as G  # noqa: E402

    d = cache_dir(sample)
    R_map, _ = G.read_map(PC.PG.guarded(d / "recon.mic.gz"), N_VOXELS)
    rec = PC.load_gz(d / "records.json.gz")
    return {
        int(v): dict(R=np.asarray(R_map[int(v)], np.float64), hit=float(rec[str(v)]["hit_ratio"]))
        for v in voxels
    }
