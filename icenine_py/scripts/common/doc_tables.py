"""Deterministic markdown tables for documentation, generated from summary data.

Blocks are written between ``<!-- table:NAME -->`` and ``<!-- /table:NAME -->`` markers so that
``scripts/dev/sync_doc_tables.py`` can copy them into a doc.
"""

import math
from pathlib import Path
from typing import Any, Dict, List, Mapping, Optional, Sequence

import numpy as np

NOT_EVALUATED = "not evaluated"
NA = "n/a"


def _fmt(value: Any, fmt: Optional[str]) -> str:
    if value is None:
        return NA
    if isinstance(value, bool):
        return str(value)
    if isinstance(value, float) and math.isnan(value):
        return NA
    if isinstance(value, int) or (isinstance(value, np.integer)):
        return format(int(value), fmt or "d")
    if isinstance(value, (float, np.floating)):
        return format(float(value), fmt or ".3g")
    return str(value)


def markdown_table(
    rows: Sequence[Optional[Mapping[str, Any]]],
    columns: Sequence[str],
    formats: Optional[Mapping[str, str]] = None,
    labels: Optional[Sequence[str]] = None,
) -> str:
    """Render ``rows`` (dicts) as a markdown table with the given column order.

    Numbers use 3 significant figures unless ``formats[column]`` gives a format spec; None and
    NaN print "n/a"; ints use "d" and floats ".3g" by default. A row that is None prints "not
    evaluated" in every cell but the first, which holds ``labels[i]`` when labels are given.
    """
    formats = formats or {}
    lines = ["| " + " | ".join(columns) + " |", "|" + "|".join("---" for _ in columns) + "|"]
    for i, row in enumerate(rows):
        if row is None:
            first = labels[i] if labels is not None else NOT_EVALUATED
            cells = [first] + [NOT_EVALUATED] * (len(columns) - 1)
        else:
            cells = [_fmt(row.get(c), formats.get(c)) for c in columns]
        lines.append("| " + " | ".join(cells) + " |")
    return "\n".join(lines)


def write_tables(path: Path, tables: Mapping[str, str]) -> None:
    """Write ``tables`` (name -> markdown) as marker-delimited blocks, in sorted name order."""
    parts: List[str] = []
    for name in sorted(tables):
        parts.append(f"<!-- table:{name} -->\n{tables[name].strip()}\n<!-- /table:{name} -->\n")
    Path(path).write_text("\n".join(parts))


def read_tables(path: Path) -> Dict[str, str]:
    """Parse a file written by ``write_tables`` back into {name: markdown}."""
    import re

    text = Path(path).read_text()
    pat = re.compile(r"<!-- table:([\w.-]+) -->\n(.*?)\n<!-- /table:\1 -->", re.S)
    return {m.group(1): m.group(2) for m in pat.finditer(text)}
