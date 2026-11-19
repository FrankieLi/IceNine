"""AR5 step 23: support for the new clean test sample ``held7`` in the script layer (no older
script is edited). Step 22's ``defer_samples.py`` knows ``held3`` to ``held6``; this file adds
``held7`` to the same table at import time and uses its machinery unchanged (the same text
replacements, now with the generator seed 7).

``held7`` is ``held`` (AR3 step 11: 250 Voronoi grains on the 24,570-voxel gradient grid, gradients
0 to 4 degrees, clean images, Q-max 8) built by the same script with the generator seed changed to
7: a new sample for a clean test, not an extra random seed for a comparison.

CLI (the build stages of ``ball_build_held7.py`` use it)::

    AR3_SAMPLE=held7 uv run python scripts/autoresearch/ar5/ball_samples.py make|stack [args]
"""

import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import defer_samples as DS  # noqa: E402

DS.SS.SEEDS["held7"] = 7  # held3: 3, held4: 4, held5: 5, held6: 6 (step 22), held7: 7

SEEDS = DS.SEEDS
which_sample = DS.which_sample
wanted = DS.wanted
boot = DS.boot
patched_source = DS.patched_source
run_patched = DS.run_patched

if __name__ == "__main__":
    if len(sys.argv) < 2 or sys.argv[1] not in DS.SS.LS.SCRIPTS:
        raise SystemExit("usage: ball_samples.py make|stack [args]")
    run_patched(sys.argv[1], sys.argv[2:])
