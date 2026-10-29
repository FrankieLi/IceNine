"""Tests of scripts/dev/sync_doc_tables.py: --section and the exit-2 marker listing.

The whole-file round trip and the missing-marker error are also covered in test_dev_tooling.py.
Every test works on files in a temporary directory.
"""

import subprocess
import sys
from pathlib import Path
from typing import Dict, List, Tuple

import pytest

REPO = Path(__file__).resolve().parents[2]
DEV = REPO / "scripts" / "dev"
SCRIPT = DEV / "sync_doc_tables.py"
sys.path.insert(0, str(DEV))
sys.path.insert(0, str(REPO / "icenine_py" / "scripts" / "common"))

import audit_numbers  # noqa: E402
import doc_tables  # noqa: E402
import sync_doc_tables  # noqa: E402

DOC = """\
# Report

<!-- table:top -->
stale top
<!-- /table:top -->

## Alpha (2026-10-01)

<!-- table:a1 -->
stale a1
<!-- /table:a1 -->

### Alpha detail

```bash
# a shell comment, not a heading
uv run python x.py
```

<!-- table:a2 -->
stale a2
<!-- /table:a2 -->

## Beta (2026-10-02)

<!-- table:b1 -->
stale b1
<!-- /table:b1 -->

## Gamma, no tables

text
"""


def _run(doc: Path, tables: Path, *extra: str) -> "subprocess.CompletedProcess[str]":
    cmd = [sys.executable, str(SCRIPT), "--doc", str(doc), "--tables", str(tables), *extra]
    return subprocess.run(cmd, capture_output=True, text=True)


def _files(tmp_path: Path, blocks: Dict[str, str], doc_text: str = DOC) -> Tuple[Path, Path]:
    tables = tmp_path / "tables.md"
    doc_tables.write_tables(tables, blocks)
    doc = tmp_path / "doc.md"
    doc.write_text(doc_text)
    return doc, tables


def _outside(text: str, start: str, end: str) -> Tuple[str, str]:
    i, j = text.index(start), text.index(end)
    return text[:i], text[j:]


def test_section_syncs_only_its_markers(tmp_path: Path) -> None:
    # only Alpha's blocks exist: the whole file fails, the section syncs
    doc, tables = _files(tmp_path, {"a1": "| new a1 |", "a2": "| new a2 |"})
    assert _run(doc, tables).returncode == 2
    assert doc.read_text() == DOC
    res = _run(doc, tables, "--section", "Alpha")
    assert res.returncode == 0, res.stderr
    new = doc.read_text()
    assert "| new a1 |" in new and "| new a2 |" in new  # the subsection's marker too
    assert "stale top" in new and "stale b1" in new  # outside the section: untouched
    assert _outside(new, "## Alpha", "## Beta") == _outside(DOC, "## Alpha", "## Beta")
    assert _run(doc, tables, "--section", "Alpha").returncode == 0
    assert doc.read_text() == new  # idempotent


def test_section_check_never_writes(tmp_path: Path) -> None:
    doc, tables = _files(tmp_path, {"b1": "| new b1 |"})
    res = _run(doc, tables, "--section", "Beta", "--check")
    assert res.returncode == 1 and "tables differ" in res.stderr
    assert doc.read_text() == DOC
    assert _run(doc, tables, "--section", "Beta").returncode == 0
    assert _run(doc, tables, "--section", "Beta", "--check").returncode == 0


def test_section_missing_block_lists_only_section_markers(tmp_path: Path) -> None:
    doc, tables = _files(tmp_path, {"a1": "| new a1 |"})
    res = _run(doc, tables, "--section", "Alpha")
    assert res.returncode == 2
    assert "1 marker(s) in section 'Alpha'" in res.stderr
    assert 'line 20 under "### Alpha detail": a2' in res.stderr
    assert "top" not in res.stderr and "b1" not in res.stderr and "hint:" not in res.stderr
    assert doc.read_text() == DOC


def test_whole_file_exit2_lists_markers_by_heading(tmp_path: Path) -> None:
    doc, tables = _files(tmp_path, {"a1": "x"})
    res = _run(doc, tables)
    assert res.returncode == 2
    lines = res.stderr.splitlines()
    assert lines[0].startswith("error: 3 marker(s) of ")
    assert lines[1] == '  line 3 under "# Report": top'
    assert lines[2] == '  line 20 under "### Alpha detail": a2'
    assert lines[3] == '  line 26 under "## Beta (2026-10-02)": b1'
    assert lines[4].startswith('hint: --section "<heading>"')
    assert doc.read_text() == DOC


def test_marker_before_first_heading(tmp_path: Path) -> None:
    doc, tables = _files(tmp_path, {}, "<!-- table:z -->\n<!-- /table:z -->\n# H\n")
    res = _run(doc, tables)
    assert res.returncode == 2 and "line 1 (before the first heading): z" in res.stderr


def test_section_unused_blocks_relative_to_section(tmp_path: Path) -> None:
    doc, tables = _files(tmp_path, {"b1": "| b |", "a1": "| a |"})
    res = _run(doc, tables, "--section", "Beta")
    assert res.returncode == 0 and "unused generated blocks: a1" in res.stderr


def test_section_not_found_and_no_markers(tmp_path: Path) -> None:
    doc, tables = _files(tmp_path, {"b1": "| b |"})
    res = _run(doc, tables, "--section", "Delta")
    assert res.returncode == 2 and "section not found: 'Delta'" in res.stderr
    res = _run(doc, tables, "--section", "Gamma", "--check")
    assert res.returncode == 0 and "no table markers in section 'Gamma'" in res.stderr
    assert doc.read_text() == DOC


def test_section_first_match_and_note(tmp_path: Path) -> None:
    # "Alpha" is in two headings: the first wins (as audit_numbers.py --section), with a note
    doc, tables = _files(tmp_path, {"a1": "| a |", "a2": "| a2 |"})
    res = _run(doc, tables, "--section", "Alpha", "--check")
    assert res.returncode == 1
    assert (
        "note: 2 headings contain 'Alpha'; using the first (line 7), not line(s) 13" in res.stderr
    )


def test_marker_crossing_section_end_is_an_error(tmp_path: Path) -> None:
    text = "## A\n<!-- table:x -->\nold\n## B\n<!-- /table:x -->\n"
    doc, tables = _files(tmp_path, {"x": "| x |"}, text)
    res = _run(doc, tables, "--section", "A")
    assert res.returncode == 2 and "without closing in it: x" in res.stderr
    assert doc.read_text() == text


def test_section_to_end_of_file_without_trailing_newline(tmp_path: Path) -> None:
    text = "## A\nkeep\n## B\n<!-- table:x -->\nold\n<!-- /table:x -->"
    doc, tables = _files(tmp_path, {"x": "| x |"}, text)
    assert _run(doc, tables, "--section", "B").returncode == 0
    assert doc.read_text() == "## A\nkeep\n## B\n<!-- table:x -->\n| x |\n<!-- /table:x -->"


def test_headings_skip_fenced_code() -> None:
    lines = DOC.splitlines()
    hs = sync_doc_tables.headings(lines)
    assert [h.text for h in hs] == [
        "Report",
        "Alpha (2026-10-01)",
        "Alpha detail",
        "Beta (2026-10-02)",
        "Gamma, no tables",
    ]
    start, end, _ = sync_doc_tables.find_section(lines, "Alpha detail")
    assert lines[start] == "### Alpha detail" and lines[end] == "## Beta (2026-10-02)"


@pytest.mark.parametrize("heading", ["Report", "Alpha", "Beta", "Gamma", "Alpha detail"])
def test_find_section_agrees_with_audit_numbers_without_fences(heading: str) -> None:
    # same heading rule as audit_numbers.py --section on a doc with no fenced code
    lines: List[str] = [ln for ln in DOC.splitlines() if not ln.startswith(("```", "# a shell"))]
    lines = [ln for ln in lines if ln != "uv run python x.py"]
    start, end, _ = sync_doc_tables.find_section(lines, heading)
    assert (start, end) == audit_numbers.find_section(lines, heading)


def test_several_tables_files_merge(tmp_path: Path) -> None:
    # a section whose markers come from two generated files (S2 Step 4: s2_results and taska)
    doc, first = _files(tmp_path, {"a1": "| new a1 |"})
    second = tmp_path / "second.md"
    doc_tables.write_tables(second, {"a2": "| new a2 |"})
    res = _run(doc, first, str(second), "--section", "Alpha")
    assert res.returncode == 0, res.stderr
    assert "| new a1 |" in doc.read_text() and "| new a2 |" in doc.read_text()
    assert _run(doc, first, str(second), "--section", "Alpha", "--check").returncode == 0


def test_several_tables_files_conflict(tmp_path: Path) -> None:
    doc, first = _files(tmp_path, {"a1": "| one |"})
    second = tmp_path / "second.md"
    doc_tables.write_tables(second, {"a1": "| two |"})
    res = _run(doc, first, str(second), "--section", "Alpha")
    assert res.returncode == 2 and "block a1 differs between --tables files" in res.stderr
    assert doc.read_text() == DOC
