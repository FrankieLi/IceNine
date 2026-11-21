"""AR5 step 12 (review fix 1 of the second review): refuse to write outside this step's own cache.

A run or pass writes into ``<cache>/<arm>_<sample>/`` and ``<cache>/passes_<sample>/``. If one of
those names is a link (a whole-folder link, or a per-entry link to a stored file), the write would
land in the folder the link points to, which may be another worktree's cache. This guard resolves
every output path the job will write (``Path.resolve``: links followed, missing parts kept) and
refuses when the real path is not inside the real path of the cache.

  uv run python scripts/autoresearch/ar5/stallhit_cache_guard.py --cache DIR PATH [PATH ...]

A relative ``PATH`` is taken relative to ``--cache``. Exit 0 when every path is inside, exit 2 with
one line per path outside (on stderr) otherwise. No file is opened or written. The runners
(``stallhit_phase_b.sh``, ``stallhit_run_all.sh``) call it before each run and pass, since the
frozen run and pass scripts (``stallhit_run.py``, ``order_pass.py``) are not edited; the
unfrozen run scripts (``stallhit_hl.py``, ``stallhit_lib_run.py``) also call ``require_inside``.
No import of the study's other modules.
"""

import argparse
import sys
from pathlib import Path
from typing import Iterable, List, Optional, Union

PathLike = Union[str, Path]


def outside(cache: PathLike, paths: Iterable[PathLike]) -> List[str]:
    """Return one line per path whose real path is not inside the real path of ``cache``."""
    root = Path(cache).resolve()
    bad: List[str] = []
    for p in paths:
        q = Path(p)
        if not q.is_absolute():
            q = Path(cache) / q
        real = q.resolve()
        try:
            real.relative_to(root)
        except ValueError:
            bad.append(f"{q} -> {real} is outside the cache {root}")
    return bad


def require_inside(cache: PathLike, paths: Iterable[PathLike]) -> None:
    """Raise ``SystemExit`` (code 2) when an output path's real path lies outside the cache."""
    bad = outside(cache, paths)
    if bad:
        print("refusing to write: " + "; ".join(bad), file=sys.stderr)
        raise SystemExit(2)


def main(argv: Optional[List[str]] = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--cache", required=True)
    ap.add_argument("paths", nargs="+")
    a = ap.parse_args(argv)
    bad = outside(a.cache, a.paths)
    for line in bad:
        print("refusing to write: " + line, file=sys.stderr)
    return 2 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
