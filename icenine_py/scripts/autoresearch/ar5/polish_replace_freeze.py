"""AR5 step 14: step 14's freeze check replaces step 12's for ``held`` and ``held2`` (pure).

Why: ``ar5/stallhit/frozen.json`` hashed the library before step 12's phase B (commit 919c55e, the
opt-in key ``SeedIndexMCStallHit``, merged at fe3d277) changed ``config_file.py``,
``orientation_search.py`` and ``reconstructor.py``; step 12's own check can therefore no longer
agree, and ``stallhit_common`` refuses ``held`` and ``held2`` with it. Step 12 did the same to the
checks of steps 8, 9 and 11 (``stallhit_common.py``, "Steps 8, 9 and 11 froze the library as it was
then"). Step 14's own freeze (``ar5/polish/frozen.json``, committed at 84a51a4: the library as it is
now and the scripts that define the rule) is checked FIRST, for ``held`` and ``held2``; only when
it passes is step 12's check switched off. No hash is re-recorded; no frozen file is edited.
"""

from typing import Any, Callable


def install(sample: str, polish_problems: Callable[[], str], stallhit_freeze: Any) -> None:
    """For ``held`` and ``held2``: refuse (``SystemExit``) unless step 14's freeze check passes,
    then replace step 12's ``problems`` by one that agrees. Other samples: nothing changes."""
    if sample not in ("held", "held2"):
        return
    why = polish_problems()
    if why:
        raise SystemExit(f"polish rule: held and held2 are not read before the freeze ({why})")
    stallhit_freeze.problems = lambda: ""
