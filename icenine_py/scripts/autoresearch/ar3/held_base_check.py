"""AR3 step 11: the base check on ``big``, repeated before the run on ``held``.

  scripts/autoresearch/ar3/held_bg.sh big base_check ar3/held_base_check.py

``ar5/fit_base_check.py`` unchanged (it compares ``<cache>/S8x_big/`` of ``ar5/fit_run.py --name
S8x``, both changes off, with the stored S8 run of ``big``: 7,201,863 evaluations, 67 wrong before
the pass), with its output sent to ``benchmarks/autoresearch/ar3/held/base_check.json`` so that the
file of AR5 step 2 is not overwritten.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "ar5"))
import fit_base_check as B  # noqa: E402
import fit_common as FC  # noqa: E402

FC.FIT = FC.ICE / "benchmarks" / "autoresearch" / "ar3" / "held"
B.main()
