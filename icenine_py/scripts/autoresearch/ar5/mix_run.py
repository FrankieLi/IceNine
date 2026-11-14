"""AR5 step 3, Part 2: whole-layer runs of S11 and S11a (script layer; library unchanged).

  scripts/autoresearch/ar5/mix_chain.sh <big|grad|mid> S11 --seed-arm X200 --no-mc-fallback
  scripts/autoresearch/ar5/mix_chain.sh <big|grad|mid> S11a --seed-arm X200

This is step 2's ``fit_run.py`` with the two extra seed arms of ``mix_common`` (``X60``,
``X200``); every option of ``fit_run.py`` works as there. ``S11`` is ``S10`` (no MC fallback after
CMA in a neighbour local optimization or revisit) plus the mixed seed refinement; ``S11a`` keeps
the MC fallback in the neighbour and revisit fits (the mixed seed refinement only). Besides the
files of ``fit_run.py`` it writes ``mix_log.json`` in the run folder: one row per candidate refined
(the CMA end hit ratio, whether MC ran, the evaluations) for the analysis of how often MC ran.
"""

import json
import sys

import mix_common as MC  # noqa: E402  (first: sets the cache, installs the arms)

import fit_common as FC  # noqa: E402
import fit_run as FRUN  # noqa: E402


def main() -> None:
    name = sys.argv[sys.argv.index("--name") + 1]
    smoke = "--smoke" in sys.argv
    FRUN.main()
    out = FC.CACHE / f"{name}_{FC.SAMPLE}{'_smoke' if smoke else ''}"
    (out / "mix_log.json").write_text(json.dumps(MC.LOG) + "\n")


if __name__ == "__main__":
    main()
