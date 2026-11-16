"""AR5 step 21: support for the new clean test sample ``held5`` in the script layer (no older script
is edited, no hash of any older freeze moves). Step 18's ``learnskip_samples.py`` knows ``held3``
and ``held4``; this file adds ``held5`` to its table at import time and uses its machinery
unchanged (the same text replacements, now with the generator seed 5).

``held5`` is ``held`` (AR3 step 11: 250 Voronoi grains on the 24,570-voxel gradient grid, gradients
0 to 4 degrees, clean images, Q-max 8) built by the same script with the generator seed changed to
5: a new sample for a clean test, not an extra random seed for a comparison.

CLI (the build stages of ``screen_build_held5.py`` use it)::

    AR3_SAMPLE=held5 uv run python scripts/autoresearch/ar5/screen_samples.py make|stack [args]
"""

import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import learnskip_samples as LS  # noqa: E402

LS.SEEDS["held5"] = 5  # held3: 3, held4: 4 (step 18), held5: 5

SEEDS = LS.SEEDS
which_sample = LS.which_sample
wanted = LS.wanted
boot = LS.boot
patched_source = LS.patched_source
run_patched = LS.run_patched

if __name__ == "__main__":
    if len(sys.argv) < 2 or sys.argv[1] not in LS.SCRIPTS:
        raise SystemExit("usage: screen_samples.py make|stack [args]")
    run_patched(sys.argv[1], sys.argv[2:])
