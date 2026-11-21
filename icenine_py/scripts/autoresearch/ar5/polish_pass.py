"""AR5 step 14: the boundary pass of arm ``J`` (``order_pass.py`` after ``polish_common``).

  RESEED_SAMPLE=<sample> POLISH_CACHE=<dir> uv run python scripts/autoresearch/ar5/polish_pass.py \
      --name J [--aligned 1]

Step 8's ``order_pass.py`` unchanged (aligned draws, step 5's rule ``hit``) on ``<cache>/J_<sample>/
recon.mic.gz``; it only sets the polish cache and the ``held2`` freeze guard first.
"""

import polish_common as PC  # noqa: E402,F401  (first: cache, guards)
import polish_freeze as PF  # noqa: E402

PF.check_or_exit()

import order_pass as OP  # noqa: E402

if __name__ == "__main__":
    OP.main()
