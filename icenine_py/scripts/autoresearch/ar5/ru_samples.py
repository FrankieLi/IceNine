"""AR5 step 29: support for the new clean test sample ``held8`` in the script layer (no older script
is edited). Step 23's ``ball_samples.py`` knows ``held3`` to ``held7``; this file adds ``held8`` to
the same table at import time and uses its machinery unchanged.

``held8`` is ``held`` (AR3 step 11: 250 Voronoi grains on the 24,570-voxel gradient grid, gradients
0 to 4 degrees, clean images, Q-max 8) built by the same script with the generator seed changed to
8: a new sample for a clean test, not an extra random seed for a comparison.

CLI (the build stages of ``ru_build_held8.py`` use it)::

    AR3_SAMPLE=held8 uv run python scripts/autoresearch/ar5/ru_samples.py make|stack [args]
"""

import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import ball_samples as BSM  # noqa: E402  (held3 to held7)

DS = BSM.DS
DS.SS.SEEDS["held8"] = 8  # held7: 7 (step 23), held8: 8

SEEDS = DS.SEEDS
which_sample = DS.which_sample
wanted = DS.wanted
boot = DS.boot
patched_source = DS.patched_source
run_patched = DS.run_patched

if __name__ == "__main__":
    if len(sys.argv) < 2 or sys.argv[1] not in DS.SS.LS.SCRIPTS:
        raise SystemExit("usage: ru_samples.py make|stack [args]")
    run_patched(sys.argv[1], sys.argv[2:])
