"""AR5 step 12, review fixes 2 and 6: two checks of the finished arm H run of one sample.

  RESEED_SAMPLE=<held2|held|mid|grad|big> STALLHIT_CACHE=<cache> \
      uv run python scripts/autoresearch/ar5/stallhit_extra.py

(1) Same candidate (fix 6; prediction 1's "same candidate" clause): for every seed step that both
arms ran (same voxel, same occurrence) and that has fits in both, the step's stop position
(``stop``), the position of its taken end and its fits (position, candidate, kind) are compared
between H and G. (2) The H and G maps (fix 2): voxels whose orientation differs (any matrix entry),
the largest misorientation between the two maps, and the error of the differing voxels in both
arms (the largest among those not wrong in both, and how many are wrong in both), before and after
the boundary pass. Reads only finished runs: ``<cache>/{G,H}_<s>/e1_log.json.gz`` and
``recon.mic.gz``, ``<cache>/passes_<s>/after_{G,H}.npy``, and the sample's truth (scoring only).
Writes ``benchmarks/autoresearch/ar5/stallhit/extra_<s>.json``.

It is a separate script because ``stallhit_stats.py`` imports ``stallhit_common``, whose guard
refuses ``held`` and ``held2`` once the library differs from the frozen one (phase B added the key
``SeedIndexMCStallHit``); the guard protected the runs, which are finished. This script reads no
frozen helper with a guard and changes nothing in the cache.
"""

import gzip
import json
import os
import sys
from pathlib import Path
from typing import Any, Dict, List

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "ar3"))
sys.path.insert(0, str(HERE))

import reseed_common as RC  # noqa: E402  (first: RESEED_SAMPLE -> AR3_SAMPLE)

import gradient_common as G  # noqa: E402
import numpy as np  # noqa: E402
import order_logs as OL  # noqa: E402

ICE = HERE.parents[2]
SH = ICE / "benchmarks" / "autoresearch" / "ar5" / "stallhit"
SCACHE = Path(os.environ.get("STALLHIT_CACHE", str(SH / "cache")))
N = 24570


def read_gz(p: Path) -> Any:
    with gzip.open(p, "rt") as f:
        return json.load(f)


def same_candidate(lg: List[Dict[str, Any]], lh: List[Dict[str, Any]]) -> Dict[str, Any]:
    """Matched seed steps with fits in both arms: how many stop at the same position, take the end
    of the same position and ran the same fits (position, candidate, kind)."""
    ig, ih = OL.index(lg), OL.index(lh)
    n = same_stop = same_taken = same_fits = 0
    differ: List[List[int]] = []
    for k in sorted(set(ig) & set(ih)):
        eg, eh = ig[k], ih[k]
        if not eg["fits"] or not eh["fits"]:
            continue
        n += 1
        a = int(eg["stop"]) == int(eh["stop"])
        b = int(OL.taken_end(eg)["pos"]) == int(OL.taken_end(eh)["pos"])
        c = [tuple(f[:3]) for f in eg["fits"]] == [tuple(f[:3]) for f in eh["fits"]]
        same_stop += int(a)
        same_taken += int(b)
        same_fits += int(c)
        if not (a and b and c):
            differ.append([int(k[0]), int(k[1])])
    return dict(
        matched_steps_with_fits=n,
        same_stop=same_stop,
        same_taken_position=same_taken,
        same_fits=same_fits,
        all_same=bool(n == same_stop == same_taken == same_fits),
        differing_steps=differ,
    )


def map_diff(Rg: np.ndarray, Rh: np.ndarray, T: Any) -> Dict[str, Any]:
    """Voxels whose orientation differs between the two maps, their misorientation, and their
    error in both arms (wrong = over ``G.WRONG`` degrees)."""
    d = np.any(Rg.reshape(len(Rg), -1) != Rh.reshape(len(Rh), -1), axis=1)
    idx = np.flatnonzero(d)
    out: Dict[str, Any] = dict(voxels_differ=int(len(idx)))
    if len(idx) == 0:
        out.update(
            max_mis_deg=0.0,
            max_err_not_both_wrong_deg=None,
            both_wrong=0,
            one_wrong=0,
        )
        return out
    mis = np.asarray(T.mis(Rg[idx], Rh[idx]), dtype=float)
    eg = np.asarray(T.err(Rg[idx], idx), dtype=float)
    eh = np.asarray(T.err(Rh[idx], idx), dtype=float)
    wg, wh = eg > G.WRONG, eh > G.WRONG
    both = wg & wh
    rest = ~both
    out.update(
        max_mis_deg=float(mis.max()),
        median_mis_deg=float(np.median(mis)),
        max_err_not_both_wrong_deg=float(np.maximum(eg, eh)[rest].max()) if rest.any() else None,
        both_wrong=int(both.sum()),
        one_wrong=int((wg ^ wh).sum()),
        both_wrong_err_deg=[float(x) for x in np.maximum(eg, eh)[both]],
    )
    return out


def main() -> None:
    s = RC.SAMPLE
    gd, hd = SCACHE / f"G_{s}", SCACHE / f"H_{s}"
    lg, lh = read_gz(gd / "e1_log.json.gz"), read_gz(hd / "e1_log.json.gz")
    sc = same_candidate(lg, lh)
    T = G.Truth()
    Rg_b, _ = G.read_map(gd / "recon.mic.gz", N)
    Rh_b, _ = G.read_map(hd / "recon.mic.gz", N)
    pd = SCACHE / f"passes_{s}"
    Rg_a, Rh_a = np.load(pd / "after_G.npy"), np.load(pd / "after_H.npy")
    out = dict(
        sample=s,
        same_candidate=sc,
        maps=dict(before=map_diff(Rg_b, Rh_b, T), after=map_diff(Rg_a, Rh_a, T)),
        note=(
            "maps: H against G (aligned draws), voxels whose stored orientation differs; error "
            "against the truth, symmetry-reduced; wrong = over 1 degree"
        ),
    )
    (SH / f"extra_{s}.json").write_text(json.dumps(out, indent=1) + "\n")
    print(
        json.dumps(
            dict(
                sample=s,
                same_candidate={k: v for k, v in sc.items() if k != "differing_steps"},
                maps=out["maps"],
            ),
            indent=1,
        )
    )  # noqa: E501


if __name__ == "__main__":
    main()
