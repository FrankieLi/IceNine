"""AR5 step 22: support for the new clean test sample ``held6`` in the script layer (no older script
is edited). Step 21's ``screen_samples.py`` knows ``held3`` to ``held5``; this file adds ``held6``
to the same table at import time and uses its machinery unchanged (the same text replacements, now
with the generator seed 6).

``held6`` is ``held`` (AR3 step 11: 250 Voronoi grains on the 24,570-voxel gradient grid, gradients
0 to 4 degrees, clean images, Q-max 8) built by the same script with the generator seed changed to
6: a new sample for a clean test, not an extra random seed for a comparison.

CLI (the build stages of ``defer_build_held6.py`` use it)::

    AR3_SAMPLE=held6 uv run python scripts/autoresearch/ar5/defer_samples.py make|stack [args]
"""

import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import screen_samples as SS  # noqa: E402

SS.SEEDS["held6"] = 6  # held3: 3, held4: 4, held5: 5 (step 21), held6: 6

SEEDS = SS.SEEDS
which_sample = SS.which_sample
wanted = SS.wanted
boot = SS.boot
patched_source = SS.patched_source
run_patched = SS.run_patched

if __name__ == "__main__":
    if len(sys.argv) < 2 or sys.argv[1] not in SS.LS.SCRIPTS:
        raise SystemExit("usage: defer_samples.py make|stack [args]")
    run_patched(sys.argv[1], sys.argv[2:])
