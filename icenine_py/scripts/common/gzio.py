"""Read result files that may be stored gzip-compressed (``name.json`` or ``name.json.gz``).

Large result JSON under ``icenine_py/benchmarks/`` is kept in git as ``.gz`` to save space. Readers
ask for the logical name (``runs/cpprule_full_v15901.json``) and get whichever copy exists. If both
exist, the plain file wins (a fresh rerun writes plain files next to the stored ``.gz``).

    sys.path.insert(0, str(ROOT / "scripts" / "common"))
    from gzio import load_json, glob_paths, resolve
"""

import glob as _glob
import gzip
import json
from pathlib import Path
from typing import IO, Any, List, Union

PathLike = Union[str, Path]
GZ = ".gz"


def logical(path: PathLike) -> Path:
    """The name without a trailing ``.gz``."""
    p = Path(path)
    return p.with_name(p.name[: -len(GZ)]) if p.name.endswith(GZ) else p


def resolve(path: PathLike) -> Path:
    """The existing file for a logical name: ``path`` itself, else ``path + '.gz'``.

    A path that already ends in ``.gz`` is returned as is. Raises FileNotFoundError if neither
    copy exists."""
    p = Path(path)
    if p.name.endswith(GZ) or p.exists():
        if p.exists():
            return p
        raise FileNotFoundError(str(p))
    gz = p.with_name(p.name + GZ)
    if gz.exists():
        return gz
    raise FileNotFoundError(f"{p} (nor {gz.name})")


def exists(path: PathLike) -> bool:
    """True if the logical name exists plain or gzip-compressed."""
    try:
        resolve(path)
    except FileNotFoundError:
        return False
    return True


def open_text(path: PathLike) -> IO[str]:
    """Open a logical name for reading text, decompressing a ``.gz`` copy transparently."""
    p = resolve(path)
    if p.name.endswith(GZ):
        return gzip.open(p, "rt", encoding="utf-8")
    return open(p, "r", encoding="utf-8")


def read_text(path: PathLike) -> str:
    with open_text(path) as fh:
        return fh.read()


def load_json(path: PathLike) -> Any:
    """``json.load`` of a logical name (plain or ``.gz``)."""
    with open_text(path) as fh:
        return json.load(fh)


def write_json_like(path: PathLike, text: str) -> Path:
    """Write ``text`` back to the stored copy of a logical name, keeping its compression.

    Writes ``path + '.gz'`` if that is the copy that exists (and no plain file does), otherwise
    the plain file. Returns the path written. The gzip header carries no name and mtime 0, so the
    output is byte-reproducible."""
    p = Path(path)
    try:
        target = resolve(p)
    except FileNotFoundError:
        target = p
    if target.name.endswith(GZ):
        with open(target, "wb") as raw:
            with gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as gz:
                gz.write(text.encode("utf-8"))
    else:
        target.write_text(text, encoding="utf-8")
    return target


def glob_paths(pattern: PathLike) -> List[Path]:
    """Sorted logical paths matching ``pattern`` plain or with ``.gz`` appended (deduplicated).

    The returned paths are logical names (no ``.gz``); pass them to ``load_json``/``resolve``."""
    pat = str(pattern)
    hits = {logical(h) for h in _glob.glob(pat)}
    hits |= {logical(h) for h in _glob.glob(pat + GZ)}
    return sorted(hits)
