"""AR3 step 9, part 1: the one logged rerun of S5, on the first gradient sample (cap 90 minutes).

  RESEED_SAMPLE=grad RESEED_CACHE=<start/cache> uv run python \\
      scripts/autoresearch/ar3/start_logged_grad.py

Step 8's logged run exists for ``big`` only; the first gradient sample has the seed steps but not
the neighbour fits. This is ``reseed_logged_run.py`` unchanged (same recorders, same call of
``run_seed.run("S5", ...)``, same voxel order) with the sample set to ``grad`` and the output
folder renamed to ``<RESEED_CACHE>/S5_grad_logged``. Its counters are compared with the stored
run S5 on the first gradient sample by ``start_part1.py`` (a baseline that does not reproduce
stops the task).
"""

import os

os.environ["RESEED_SAMPLE"] = "grad"

import reseed_common as RC  # noqa: E402,I001

assert RC.SAMPLE == "grad"

import reseed_logged_run as RL  # noqa: E402

RL.OUTD = RC.CACHE / "S5_grad_logged"

if __name__ == "__main__":
    RL.main()
