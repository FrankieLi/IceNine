"""AR5 step 8, Phase A: the boundary pass of one run of this step, with aligned draws.

  RESEED_SAMPLE=<sample> uv run python scripts/autoresearch/ar5/order_pass.py --name A \
      [--aligned 1|0]

Step 2's ``fit_pass.py`` unchanged (``boundary_pass.py`` on the run's map, step 5's rule ``hit``),
with one change: each try draws from its own generator seeded by (voxel, ``pass``, the number of
tries made at that voxel so far) instead of one generator for the whole pass
(``order_draws.aligned_try_fit``). ``--aligned 0`` runs the unchanged pass (the library's single
generator), for ``Yl``. Reads ``<cache>/<name>_<sample>/recon.mic.gz``; writes
``<cache>/passes_<sample>/`` as ``fit_pass.py`` does.
"""

import argparse
import json
import sys

import order_common as OC  # noqa: E402,F401  (first: sets the cache)

import boundary_common as B  # noqa: E402
import fit_pass as FPASS  # noqa: E402
import order_draws as OD  # noqa: E402


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--name", required=True)
    ap.add_argument("--aligned", type=int, default=1)
    a, rest = ap.parse_known_args()
    book = OD.DrawBook()
    if a.aligned:
        B.try_fit = OD.aligned_try_fit(book, B.try_fit)
    sys.argv = ["fit_pass.py", "--name", a.name] + rest
    FPASS.main()
    # pass.json keeps the unchanged pass's generator label; this file says what was done here
    d = FPASS.PASSES / f"D4_{a.name}"
    (d / "aligned.json").write_text(
        json.dumps(dict(aligned_draws=bool(a.aligned), aligned_tries=book.n_issued)) + "\n"
    )


if __name__ == "__main__":
    main()
