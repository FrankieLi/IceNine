"""Copy generated markdown tables into a doc between matching table markers.

Usage: sync_doc_tables.py --doc DOC.md --tables GENERATED.md [--check]

Both files use ``<!-- table:NAME -->`` ... ``<!-- /table:NAME -->`` blocks (see
icenine_py/scripts/common/doc_tables.py). A marker in the doc with no generated block is an error
(exit 2); generated blocks the doc does not use are listed. With --check the doc is not written
and the exit code is 1 if it differs from what would be written.
"""

import argparse
import re
import sys
from pathlib import Path
from typing import Dict, Optional, Sequence

BLOCK_RE = re.compile(r"(<!-- table:([\w.-]+) -->\n)(.*?)(<!-- /table:\2 -->)", re.S)


def parse_blocks(text: str) -> Dict[str, str]:
    return {m.group(2): m.group(3).rstrip("\n") for m in BLOCK_RE.finditer(text)}


def sync(doc_text: str, generated: Dict[str, str]) -> str:
    """Return ``doc_text`` with every marker block replaced; raise KeyError on a missing one."""
    missing = [n for n in parse_blocks(doc_text) if n not in generated]
    if missing:
        raise KeyError(", ".join(missing))

    def repl(m: "re.Match[str]") -> str:
        return m.group(1) + generated[m.group(2)].strip("\n") + "\n" + m.group(4)

    return BLOCK_RE.sub(repl, doc_text)


def main(argv: Optional[Sequence[str]] = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--doc", required=True, type=Path)
    ap.add_argument("--tables", required=True, type=Path)
    ap.add_argument("--check", action="store_true")
    args = ap.parse_args(argv)

    doc = args.doc.read_text()
    generated = parse_blocks(args.tables.read_text())
    try:
        new = sync(doc, generated)
    except KeyError as e:
        print(f"error: no generated block for marker(s): {e.args[0]}", file=sys.stderr)
        return 2
    unused = sorted(set(generated) - set(parse_blocks(doc)))
    if unused:
        print(f"unused generated blocks: {', '.join(unused)}", file=sys.stderr)
    if new == doc:
        return 0
    if args.check:
        print(f"{args.doc}: tables differ from {args.tables}", file=sys.stderr)
        return 1
    args.doc.write_text(new)
    print(f"{args.doc}: updated")
    return 0


if __name__ == "__main__":
    sys.exit(main())
