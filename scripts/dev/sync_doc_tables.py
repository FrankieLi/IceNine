"""Copy generated markdown tables into a doc between matching table markers.

Usage: sync_doc_tables.py --doc DOC.md --tables GENERATED.md [MORE.md ...] [--section "Heading"]
                          [--check]

Both files use ``<!-- table:NAME -->`` ... ``<!-- /table:NAME -->`` blocks (see
icenine_py/scripts/common/doc_tables.py). A marker in the doc with no generated block is an error
(exit 2, listing each such marker with its line and heading); generated blocks the doc does not use
are listed. With --check the doc is not written and the exit code is 1 if it differs from what
would be written. Several --tables files are merged; a block name in two of them with different
contents is exit 2.

--section "Heading" syncs (or checks) only the markers inside that section: the first heading
whose text contains "Heading", up to the next heading of the same or a higher level (the rule of
audit_numbers.py --section; heading-like lines inside fenced code blocks are not headings). Markers
outside it are left as they are, with or without a generated block. Exit 2 if the section is not
found or a marker opens in it without closing in it; a section with no markers is a note, exit 0.
"""

import argparse
import re
import sys
from pathlib import Path
from typing import Dict, List, NamedTuple, Optional, Sequence, Tuple

BLOCK_RE = re.compile(r"(<!-- table:([\w.-]+) -->\n)(.*?)(<!-- /table:\2 -->)", re.S)
OPEN_RE = re.compile(r"<!-- table:([\w.-]+) -->")
HEADING_RE = re.compile(r"^(#+)\s+(.*?)\s*$")
FENCE_RE = re.compile(r"^\s*(```|~~~)")


class Heading(NamedTuple):
    line: int  # 0-based line index
    level: int
    text: str


class Marker(NamedTuple):
    name: str
    line: int  # 1-based line of the opening marker
    heading: str  # the nearest heading above it ("" before the first heading)


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


def headings(lines: Sequence[str]) -> List[Heading]:
    """Markdown ATX headings, skipping lines inside fenced code blocks."""
    out: List[Heading] = []
    fence: Optional[str] = None
    for i, line in enumerate(lines):
        f = FENCE_RE.match(line)
        if f:
            if fence is None:
                fence = f.group(1)
            elif f.group(1) == fence:
                fence = None
            continue
        if fence is None:
            m = HEADING_RE.match(line)
            if m:
                out.append(Heading(i, len(m.group(1)), m.group(2)))
    return out


def find_section(lines: Sequence[str], heading: str) -> Tuple[int, int, List[Heading]]:
    """(start, end, matches): 0-based line range (end exclusive) of the first section whose heading
    contains ``heading``, and every heading that contains it. Raises ValueError if none does.
    """
    hs = headings(lines)
    matches = [h for h in hs if heading in h.text]
    if not matches:
        raise ValueError(f"section not found: {heading!r}")
    first = matches[0]
    end = next(
        (h.line for h in hs if h.line > first.line and h.level <= first.level),
        len(lines),
    )
    return first.line, end, matches


def markers(text: str) -> List[Marker]:
    """Every opening table marker with its 1-based line and the nearest heading above it."""
    lines = text.splitlines()
    hs = headings(lines)
    out: List[Marker] = []
    for m in OPEN_RE.finditer(text):
        line = text.count("\n", 0, m.start())
        above = [h for h in hs if h.line < line]
        title = ("#" * above[-1].level + " " + above[-1].text) if above else ""
        out.append(Marker(m.group(1), line + 1, title))
    return out


def describe(missing: Sequence[Marker]) -> str:
    """One line per run of markers under the same heading: first line, heading, names."""
    groups: List[Tuple[int, str, List[str]]] = []
    for m in missing:
        if groups and groups[-1][1] == m.heading:
            groups[-1][2].append(m.name)
        else:
            groups.append((m.line, m.heading, [m.name]))
    return "\n".join(
        f"  line {line}"
        + (f' under "{h}"' if h else " (before the first heading)")
        + f": {', '.join(names)}"
        for line, h, names in groups
    )


def _offset(lines_keepends: Sequence[str], line: int) -> int:
    return sum(len(s) for s in lines_keepends[:line])


def main(argv: Optional[Sequence[str]] = None) -> int:
    ap = argparse.ArgumentParser(description=(__doc__ or "").split("\n")[0])
    ap.add_argument("--doc", required=True, type=Path)
    ap.add_argument("--tables", required=True, type=Path, nargs="+")
    ap.add_argument("--section", help="sync only the markers under this heading")
    ap.add_argument("--check", action="store_true")
    args = ap.parse_args(argv)

    doc = args.doc.read_text()
    generated: Dict[str, str] = {}
    for path in args.tables:
        for name, block in parse_blocks(path.read_text()).items():
            if generated.get(name, block) != block:
                print(f"error: block {name} differs between --tables files", file=sys.stderr)
                return 2
            generated[name] = block
    sources = ", ".join(str(p) for p in args.tables)
    keep = doc.splitlines(keepends=True)
    prefix, body, suffix = "", doc, ""
    lo, hi = 0, len(keep)  # 0-based line range of body in the doc (hi exclusive)
    if args.section is not None:
        try:
            start, end, matches = find_section(doc.splitlines(), args.section)
        except ValueError as e:
            print(f"error: {e} in {args.doc}", file=sys.stderr)
            return 2
        if len(matches) > 1:
            others = ", ".join(str(h.line + 1) for h in matches[1:])
            print(
                f"note: {len(matches)} headings contain {args.section!r}; using the first "
                f"(line {start + 1}), not line(s) {others}",
                file=sys.stderr,
            )
        a, b = _offset(keep, start), _offset(keep, end)
        prefix, body, suffix = doc[:a], doc[a:b], doc[b:]
        lo, hi = start, end
        unclosed = sorted({m.group(1) for m in OPEN_RE.finditer(body)} - set(parse_blocks(body)))
        if unclosed:
            print(
                f"error: marker(s) open in section {args.section!r} without closing in it: "
                + ", ".join(unclosed),
                file=sys.stderr,
            )
            return 2
        if not parse_blocks(body):
            print(f"note: no table markers in section {args.section!r}", file=sys.stderr)
            return 0

    used = parse_blocks(body)
    if any(n not in generated for n in used):
        missing = [
            m
            for m in markers(doc)
            if lo <= m.line - 1 < hi and m.name in used and m.name not in generated
        ]
        where = f" in section {args.section!r}" if args.section is not None else ""
        print(
            f"error: {len(missing)} marker(s){where} of {args.doc} have no generated block in "
            f"{sources}:\n{describe(missing)}",
            file=sys.stderr,
        )
        if args.section is None:
            print(
                'hint: --section "<heading>" syncs only the markers inside one section',
                file=sys.stderr,
            )
        return 2
    new = prefix + sync(body, generated) + suffix
    unused = sorted(set(generated) - set(used))
    if unused:
        print(f"unused generated blocks: {', '.join(unused)}", file=sys.stderr)
    if new == doc:
        return 0
    if args.check:
        print(f"{args.doc}: tables differ from {sources}", file=sys.stderr)
        return 1
    args.doc.write_text(new)
    print(f"{args.doc}: updated")
    return 0


if __name__ == "__main__":
    sys.exit(main())
