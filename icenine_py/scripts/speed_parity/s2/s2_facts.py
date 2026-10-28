"""S2 Step 0 facts: reflection counts of the cost functions the Step 0 timings use.

Writes ``facts.json`` (``--out`` directory): K of the per-call arms (``prof_eval``'s Q-max 8 and
Q-max 5 cost functions) and of the local cost function of the two search sets (ThreeVoxels with
``Example2.Simulation.config``; the first E0 clean voxel).

  uv run python scripts/speed_parity/s2/s2_facts.py --out benchmarks/s2_step0
"""

import argparse
import contextlib
import io
import json
import sys
from pathlib import Path
from typing import Any, Dict

HERE = Path(__file__).resolve().parent
ICE = HERE.parents[2]
sys.path.insert(0, str(ICE))
sys.path.insert(0, str(HERE.parent))


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    out = Path(a.out).resolve()
    facts: Dict[str, Any] = {}
    with contextlib.redirect_stdout(io.StringIO()):
        from s3_step0 import E0Set, ThreeVoxels

        S3 = ThreeVoxels()
        P = S3.P
        phase = P.vox.phase
        for name, args in (("per_call_q8", (0, 8.0, False)), ("per_call_q5", (3, 5.0, False))):
            facts[f"{name}_reflections"] = int(P.make(*args)._phase_recip_vecs[phase][0].shape[0])
        rec = S3.reconstructor("Example2.Simulation.config")
        cf = rec._make_local_cost_fn()
        facts["three_search_local_reflections"] = int(cf._phase_recip_vecs[phase][0].shape[0])
        SE = E0Set()
        _, _, ph, rec0 = SE.load(0)
        cf0 = rec0._make_local_cost_fn()
        facts["e0_search_local_reflections"] = int(cf0._phase_recip_vecs[ph][0].shape[0])
    (out / "facts.json").write_text(json.dumps(facts, indent=1) + "\n")
    print(facts)


if __name__ == "__main__":
    main()
