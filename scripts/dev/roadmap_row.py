#!/usr/bin/env python3
"""Show or edit the status cell of one ROADMAP track-table row, safely.

Usage:
    roadmap_row.py --row S4 (--status "<text>" | --state "<bold word(s)>" | --show) [--check]
                   [--doc icenine_py/docs/ROADMAP.md] [--agent I|N]

- The row is found by its first cell (exactly the row ID) in the tables under the `# Track X:`
  headings. Exactly one row must match; zero or several exit 2 with a message.
- `--status` replaces the whole last cell; `--state` replaces only the first `**...**` span in it
  and keeps the rest (exit 2 if there is none). `--show` prints the row and changes nothing.
- A row with or without a closing `|` (the Track S row S3 has none) is written back in the form it
  had. Text with a newline or an unescaped `|` is refused.
- Prints the old row (`-`) and the new row (`+`). `--check` writes nothing and exits 1 if the file
  would change. The script never commits; the caller commits the row with the change it records.
- Editing is limited to the rows the agent owns (OWN_ROWS below); any other row exits 2. The agent
  letter comes from `--agent` or `git config --local icenine.agent`; with neither the script
  refuses (it never defaults to Ilium Agent). `--show` is read-only and needs neither.

Exit codes: 0 done / nothing to change, 1 `--check` found a change, 2 refusal or usage error.
Python standard library only.
"""

import argparse
import os
import re
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import Dict, List, Optional, Tuple

# Rows each agent may edit (ROADMAP protocol rule 4 plus owner exceptions). "tracks" are the
# letters of the `# Track X:` headings whose rows are owned whole; "rows" are single extra rows in
# another track. Ilium Agent: Track S, plus row H12 in Track H (owner exception D-20261018-I3).
# Newt Agent: Track E. Every other row is refused. Change ownership only here, after a DECISIONS
# entry.
OWN_ROWS: Dict[str, Dict[str, Tuple[str, ...]]] = {
    "I": {"tracks": ("S",), "rows": ("H12",)},
    "N": {"tracks": ("E",), "rows": ()},
}

DEFAULT_DOC = "icenine_py/docs/ROADMAP.md"
TRACK_RE = re.compile(r"^# Track ([A-Z]):")
BOLD_RE = re.compile(r"\*\*(.+?)\*\*")


class RowError(Exception):
    """A refusal; the message goes to stderr and the exit code is 2."""


def split_cells(inner: str) -> List[str]:
    """Split a table row body on `|` that is not escaped as `\\|`."""
    cells: List[str] = []
    cur: List[str] = []
    i = 0
    while i < len(inner):
        ch = inner[i]
        if ch == "\\" and i + 1 < len(inner) and inner[i + 1] == "|":
            cur.append("\\|")
            i += 2
            continue
        if ch == "|":
            cells.append("".join(cur))
            cur = []
        else:
            cur.append(ch)
        i += 1
    cells.append("".join(cur))
    return cells


def parse_row(line: str) -> Tuple[List[str], bool, str]:
    """Return (raw cells, has closing pipe, text after the closing pipe) of a table row line."""
    s = line.rstrip("\r\n")
    t = s.rstrip()
    closed = t.endswith("|") and not t.endswith("\\|") and len(t) > 1
    if closed:
        return split_cells(t[1:-1]), True, s[len(t) :]
    return split_cells(s[1:]), False, ""


def build_row(cells: List[str], closed: bool, tail: str) -> str:
    return "|" + "|".join(cells) + ("|" + tail if closed else "")


def is_separator(cells: List[str]) -> bool:
    return all(re.fullmatch(r"\s*:?-{3,}:?\s*", c) for c in cells)


def find_rows(
    lines: List[str], row_id: str
) -> Tuple[List[Tuple[int, str, int]], Dict[str, List[str]]]:
    """Find rows whose first cell is `row_id` in the tables under `# Track X:` headings.

    Returns ([(line index, track letter, header cell count)], {track letter: row IDs}).
    """
    found: List[Tuple[int, str, int]] = []
    ids: Dict[str, List[str]] = {}
    track: Optional[str] = None
    in_fence = False
    header_cells = 0
    in_table = False
    table_pos = 0
    for idx, line in enumerate(lines):
        if line.lstrip().startswith("```"):
            in_fence = not in_fence
        if in_fence:
            continue
        m = TRACK_RE.match(line)
        if m:
            track = m.group(1)
            in_table = False
            continue
        if line.startswith("# "):
            track = None
            in_table = False
            continue
        if not line.startswith("|"):
            in_table = False
            continue
        if track is None:
            continue
        cells, _, _ = parse_row(line)
        if not in_table:
            in_table = True
            table_pos = 0
            header_cells = len(cells)
        table_pos += 1
        if table_pos <= 2 and (table_pos == 1 or is_separator(cells)):
            continue  # header and separator rows
        first = cells[0].strip()
        ids.setdefault(track, []).append(first)
        if first == row_id:
            found.append((idx, track, header_cells))
    return found, ids


def owns(agent: str, track: str, row_id: str) -> bool:
    own = OWN_ROWS[agent]
    return track in own["tracks"] or row_id in own["rows"]


def check_text(text: str, what: str) -> str:
    if not text.strip():
        raise RowError(f"error: {what} is empty")
    if "\n" in text or "\r" in text:
        raise RowError(f"error: {what} contains a newline")
    if re.search(r"(?<!\\)\|", text):
        raise RowError(f"error: {what} contains an unescaped '|' (write \\| for a literal bar)")
    return text.strip()


def new_last_cell(raw: str, mode: str, text: str) -> str:
    """Return the replaced raw last cell (leading and trailing spaces kept)."""
    lead = raw[: len(raw) - len(raw.lstrip())]
    trail = raw[len(raw.rstrip()) :]
    body = raw.strip()
    if mode == "status":
        return lead + text + trail
    m = BOLD_RE.search(body)
    if not m:
        raise RowError("error: the last cell has no **bold** span to replace with --state")
    state = text.strip()
    if state.startswith("**") and state.endswith("**") and len(state) > 4:
        state = state[2:-2]
    new_body = body[: m.start()] + f"**{state}**" + body[m.end() :]
    return lead + new_body + trail


def resolve_agent(arg: Optional[str]) -> str:
    agent = arg
    if not agent:
        res = subprocess.run(
            ["git", "config", "--local", "--get", "icenine.agent"], capture_output=True, text=True
        )
        agent = res.stdout.strip() if res.returncode == 0 else ""
    if agent not in OWN_ROWS:
        raise RowError(
            "error: agent unknown: pass --agent I|N or run "
            "'git config --local icenine.agent I' (or N); there is no default"
        )
    return agent


def default_doc() -> Path:
    res = subprocess.run(["git", "rev-parse", "--show-toplevel"], capture_output=True, text=True)
    root = Path(res.stdout.strip()) if res.returncode == 0 and res.stdout.strip() else Path.cwd()
    return root / DEFAULT_DOC


def run(args: argparse.Namespace) -> int:
    doc = Path(args.doc) if args.doc else default_doc()
    if not doc.is_file():
        raise RowError(f"error: {doc} not found")
    with open(doc, newline="") as fh:
        content = fh.read()
    lines = content.splitlines(keepends=True)
    found, ids = find_rows(lines, args.row)
    if not found:
        have = "; ".join(f"Track {t}: {', '.join(v)}" for t, v in sorted(ids.items()))
        raise RowError(
            f"error: no row '{args.row}' in the Track tables of {doc.name}. Rows: {have}"
        )
    if len(found) > 1:
        where = ", ".join(f"Track {t} (line {i + 1})" for i, t, _ in found)
        raise RowError(f"error: row '{args.row}' is ambiguous: found in {where}")
    idx, track, ncols = found[0]
    old_line = lines[idx]
    cells, closed, tail = parse_row(old_line)
    if len(cells) != ncols:
        raise RowError(
            f"error: row '{args.row}' has {len(cells)} cells, the table header has {ncols}; "
            "not editing a row that does not parse"
        )
    if args.show:
        print(old_line.rstrip("\r\n"))
        return 0
    agent = resolve_agent(args.agent)
    if not owns(agent, track, args.row):
        own = OWN_ROWS[agent]
        raise RowError(
            f"error: row '{args.row}' (Track {track}) is not an agent-{agent} row; agent {agent} "
            f"edits Track {', '.join(own['tracks'])}"
            + (f" and row {', '.join(own['rows'])}" if own["rows"] else "")
            + " (ROADMAP rule 4)"
        )
    if args.status is not None:
        mode, text = "status", check_text(args.status, "--status")
    else:
        mode, text = "state", check_text(args.state, "--state")
    cells[-1] = new_last_cell(cells[-1], mode, text)
    new_row = build_row(cells, closed, tail)
    old_row = old_line.rstrip("\r\n")
    eol = old_line[len(old_row) :]
    print(f"- {old_row}")
    print(f"+ {new_row}")
    if new_row == old_row:
        print("unchanged")
        return 0
    if args.check:
        print("check: the file would change")
        return 1
    lines[idx] = new_row + eol
    tmp = tempfile.NamedTemporaryFile(
        "w", newline="", dir=doc.parent, prefix=".roadmap_row.", delete=False
    )
    try:
        with tmp:
            tmp.write("".join(lines))
        os.chmod(tmp.name, doc.stat().st_mode & 0o777)
        os.replace(tmp.name, doc)
    finally:
        if os.path.exists(tmp.name):
            os.unlink(tmp.name)
    print(f"written: {doc}")
    return 0


def main(argv: Optional[List[str]] = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--row", required=True, help="row ID, the first cell, for example S4")
    mode = ap.add_mutually_exclusive_group(required=True)
    mode.add_argument("--status", help="replace the whole last cell")
    mode.add_argument("--state", help="replace the first **bold** span of the last cell")
    mode.add_argument("--show", action="store_true", help="print the row and stop")
    ap.add_argument("--check", action="store_true", help="write nothing; exit 1 on a change")
    ap.add_argument("--doc", default=None, help=f"default {DEFAULT_DOC} of the current repository")
    ap.add_argument("--agent", default=None, help="I or N (default: git config icenine.agent)")
    try:
        return run(ap.parse_args(argv))
    except RowError as exc:
        print(exc, file=sys.stderr)
        return 2


if __name__ == "__main__":
    sys.exit(main())
