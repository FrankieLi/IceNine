"""AR3 step 2 checks A and B: a run's map against a stored map (orientation on every voxel).

  check_ab.py A   # run A against ``runs/paint/recon.mic`` (156,766 evaluations, 62 entries)
  check_ab.py B   # run B against the stored ``key_on_fused`` map (2,816,897 evaluations)

Compares the orientation matrices voxel by voxel (to the .mic file's own precision, then exactly
as read) and the evaluation counts. Prints a JSON verdict and writes ``<run>/check.json``.
"""

import json
import sys
from typing import Any, Dict

import ar3_common as C
import numpy as np

from icenine.mic_file import MicFile
from run_layer import LAYER_OUT

EXPECT = {
    "A": dict(
        ref=C.OUT / "runs" / "paint" / "recon.mic",
        evals=156766,
        stages=dict(seed_index=119033, score=6697, fit_rejected=31036, seed_grid=0, fit_accepted=0),
        entries=62,
    ),
    "B": dict(ref=C.STORED_KEY_ON_FUSED / "recon.mic", evals=2816897),
}


def main() -> None:
    which = sys.argv[1]
    ex = EXPECT[which]
    run = LAYER_OUT / "runs" / which
    a = MicFile.read(str(run / "recon.mic"))
    b = MicFile.read(str(ex["ref"]))
    Ra = np.array([v.orientation for v in a.voxels])
    Rb = np.array([v.orientation for v in b.voxels])
    info = json.loads((run / "run.json").read_text())
    out: Dict[str, Any] = dict(
        run=which,
        n_voxels=len(a.voxels),
        n_orientation_differ=int((np.abs(Ra - Rb).max(axis=(1, 2)) > 0).sum()),
        max_abs_diff=float(np.abs(Ra - Rb).max()),
        evals_total=info["evals_total"],
        evals_expected=ex["evals"],
        evals_match=info["evals_total"] == ex["evals"],
    )
    if which == "A":
        ev = info["evals"]
        out["stages_match"] = all(ev.get(k, 0) == v for k, v in ex["stages"].items())
        out["n_entries"] = info["n_entries"]
        out["entries_match"] = info["n_entries"] == ex["entries"]
    out["passes"] = bool(
        out["n_orientation_differ"] == 0
        and out["evals_match"]
        and out.get("stages_match", True)
        and out.get("entries_match", True)
    )
    (run / "check.json").write_text(json.dumps(out, indent=2) + "\n")
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
