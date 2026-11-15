"""AR5 step 18 (stopped at its gate): what the stored live ``J5L`` run on ``mid`` says about the
voxels it broke against ``J5`` and about its wrong skips. Reads stored outputs only (no run).

  RESEED_SAMPLE=mid LEARNSKIP_CACHE=<dir> uv run python scripts/autoresearch/ar5/learnskip_broken.py

Broken: wrong (over 1 degree from the truth) after the pass in ``J5L`` and right in ``J5``.
For each:
the fits the live run SKIPPED at that voxel, the fits it ran, and whether any skipped fit there was
one that ``J5r`` accepted (matched by voxel, kind and occurrence; exact up to the first difference
at fit 1016 only for features, the order of fits is equal to the first wrong skip, so beyond it the
match by occurrence is an approximation, said so in the output). Wrong skips live: skipped fits
whose matched ``J5r`` fit was accepted.
"""

import gzip
import json
from collections import defaultdict
from typing import Any, Dict, List

import learnskip_common as LCM  # noqa: E402  (first)

import numpy as np  # noqa: E402

import gradient_common as G  # noqa: E402
import learnskip_gate as GT  # noqa: E402

C = LCM.LCACHE


def main() -> None:
    T = G.Truth()
    n = 24570
    out: Dict[str, Any] = {}
    after = {a: np.load(C / "passes_mid" / f"after_{a}.npy") for a in ("J5", "J5L")}
    err = {a: np.asarray(T.err(after[a]), float) for a in after}
    wrong = {a: err[a] > G.WRONG for a in err}
    broken = np.where(wrong["J5L"] & ~wrong["J5"])[0].tolist()
    fixed = np.where(~wrong["J5L"] & wrong["J5"])[0].tolist()
    rec = GT.RP.read_rows("mid")
    live: List[Dict[str, Any]] = []
    with gzip.open(C / "J5L_mid" / "fits.jsonl.gz", "rt") as f:
        for line in f:
            live.append(json.loads(line))
    acc: Dict[Any, List[bool]] = defaultdict(list)
    for r in rec:
        acc[(int(r["voxel"]), r["kind"])].append(bool(r["accepted"]))
    seen: Dict[Any, int] = defaultdict(int)
    wrong_skips = 0
    wrong_in_range = 0
    unmatched_skips = 0
    per_voxel: Dict[int, Dict[str, int]] = defaultdict(lambda: dict(skipped=0, ran=0, wrong_skip=0))
    for k, r in enumerate(live):
        key = (int(r["voxel"]), r["kind"])
        j = seen[key]
        seen[key] += 1
        v = key[0]
        if r["skipped"]:
            per_voxel[v]["skipped"] += 1
            if j < len(acc[key]):
                if acc[key][j]:
                    wrong_skips += 1
                    per_voxel[v]["wrong_skip"] += 1
                    wrong_in_range += int(k < 7232)
            else:
                unmatched_skips += 1
        else:
            per_voxel[v]["ran"] += 1
    out = dict(
        n_voxels=n,
        broken_voxels=broken,
        fixed_voxels=fixed,
        broken_detail={
            str(v): dict(
                err_J5_deg=float(err["J5"][v]), err_J5L_deg=float(err["J5L"][v]), **per_voxel[v]
            )
            for v in broken
        },
        live_skips=sum(1 for r in live if r["skipped"]),
        live_wrong_skips_matched_by_occurrence=wrong_skips,
        live_wrong_skips_in_first_7232_fits=wrong_in_range,
        live_skips_without_a_matched_recording_fit=unmatched_skips,
        note="the match by (voxel, kind, occurrence) is exact up to the first wrong skip "
        "(fit 7231) "
        "and an approximation after it, because the live run's route differs from J5's there",
    )
    (LCM.LS / "broken_mid.json").write_text(json.dumps(out, indent=1) + "\n")
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
