"""AR5 step 24: run an older step's run, pass or arm script, or this step's ``bf_run.py``, on any
sample of this step, including ``held7`` once the gate is judged, with this step's set-up
(``bf_common``) in front.

  RESEED_SAMPLE=<sample> SCREEN_CACHE=<dir> BALL_CACHE=<main cache> uv run python \
      scripts/autoresearch/ar5/bf_launch.py <script.py> [args ...]

Only these scripts are accepted: ``fit_run.py`` (arm ``S8``: ``--name S8``), ``fit_pass.py``,
``high500_run.py`` (arm ``J5``), ``high500_pass.py`` (``--name J5``, ``--name J5D`` or ``--name
J5DB``, ``--aligned 1``), ``defer_run.py`` (arm ``J5D``), ``screen_run.py`` (step 21's ``--arm
J5S``: the screen pass and the wide repair; run on ``J5DB``'s map from a sub-cache of links named as
``J5``'s, so that the script is imported unchanged) and ``bf_run.py`` (arm ``J5DB``). They are
run unchanged (``runpy``, as ``__main__``) after ``bf_common`` has refused ``held2`` and
``fine7``, checked this step's freeze (every script and library file pinned), made ``held7`` known
to the AR3 modules and pointed every cache variable at ``SCREEN_CACHE``.
"""

import runpy
import sys

import bf_common as BC  # noqa: E402  (first)

ALLOWED = (
    "fit_run.py",
    "fit_pass.py",
    "high500_run.py",
    "high500_pass.py",
    "screen_run.py",
    "defer_run.py",
    "bf_run.py",
)


def main() -> None:
    if len(sys.argv) < 2 or sys.argv[1] not in ALLOWED:
        raise SystemExit(f"bf_launch.py: the script must be one of {', '.join(ALLOWED)}")
    script = sys.argv[1]
    sys.argv = [script] + sys.argv[2:]
    runpy.run_path(str(BC.HERE / script), run_name="__main__")


if __name__ == "__main__":
    main()
