"""AR5 step 8 (task ar5h): shared set-up of the sweep-order study (script layer; no library change).

Import this module first in every ``order_*.py`` script that touches a sample, with
``RESEED_SAMPLE`` set (``big``, ``grad``, ``mid`` or ``held``). It

* imports step 3's ``mix_common`` (which installs the seed arm ``X200`` on step 2's ``fit_common``)
  and then points the cache at ``benchmarks/autoresearch/ar5/order/cache`` (untracked);
* provides ``require_frozen()``: no script reads anything of ``held`` (seed logs, replays, counts)
  before ``ar5/order/frozen.json`` is committed (order rule of the plan);
* names the stored inputs (read only, in kept sibling worktrees): the replays of step 2 and the
  ``S11`` runs and passes of step 3 for ``big``, ``grad``, ``mid``, and the ``held`` sample's
  stored ``S8`` / ``S11`` runs and passes.

Nothing under ``icenine/`` is touched.
"""

import os
import subprocess
import sys
from pathlib import Path
from typing import Dict

HERE = Path(__file__).resolve().parent
ICE = HERE.parents[2]
ORDER = ICE / "benchmarks" / "autoresearch" / "ar5" / "order"  # tracked results (small)
ORDER_CACHE = ORDER / "cache"  # untracked
FROZEN = ORDER / "frozen.json"
ORDER_CACHE.mkdir(parents=True, exist_ok=True)
(ORDER_CACHE / "logs").mkdir(parents=True, exist_ok=True)
sys.path.insert(0, str(HERE))


def frozen_commit() -> str:
    """The hash of the commit that holds ``frozen.json`` ('' when it is not committed)."""
    out = subprocess.run(
        ["git", "log", "-1", "--format=%H", "--", str(FROZEN)],
        cwd=str(ICE),
        capture_output=True,
        text=True,
    )
    return out.stdout.strip()


def require_frozen() -> str:
    """Refuse (SystemExit) unless frozen.json is committed and unchanged; returns the commit."""
    h = frozen_commit()
    if not h:
        raise SystemExit("order rule: held is not read before ar5/order/frozen.json is committed")
    st = subprocess.run(
        ["git", "status", "--porcelain", "--", str(FROZEN)],
        cwd=str(ICE),
        capture_output=True,
        text=True,
    )
    if st.stdout.strip():
        raise SystemExit("order rule: frozen.json differs from its commit")
    return h


if os.environ.get("RESEED_SAMPLE") == "held":
    require_frozen()

import mix_common as MC  # noqa: E402,F401  (installs the arms; sets the step 3 cache)

import fit_common as FC  # noqa: E402
import reseed_common as RC  # noqa: E402

# point every cache user at this step's cache (mix_common set the step 3 one)
FC.CACHE = ORDER_CACHE
RC.CACHE = ORDER_CACHE
os.environ["RESEED_CACHE"] = str(ORDER_CACHE)

SAMPLE = RC.SAMPLE
WT = RC.WORKTREES
# stored inputs in kept worktrees (read only)
MIX_STORED = WT / "agent-a5a3601b746680640/icenine_py/benchmarks/autoresearch/ar5/mix/cache"
HELD_STORED = WT / "agent-a8a1ac8e80f8913cf/icenine_py/benchmarks/autoresearch/ar5/held"
STEP7_STORED = WT / "agent-a331e7882919e74aa/icenine_py/benchmarks/autoresearch/ar5"
STORED_S11: Dict[str, Path] = {s: MIX_STORED / f"S11_{s}" for s in ("big", "grad", "mid")}
STORED_S11_PASS: Dict[str, Path] = {s: MIX_STORED / f"passes_{s}" for s in ("big", "grad", "mid")}
