"""Pre-commit checks on staged files and added lines (stdlib only; black run via subprocess).

Prints one line per problem, ``file:line: reason``, and exits 1; silent with exit 0 when clean.
Overrides (environment): ALLOW_CLAUDE_CONFIG=1, ALLOW_LARGE=1, ALLOW_ABS_PATHS=1.

CLAUDE.md/.claude files: refused, with one exception (owner decision D-20261109-N11). A commit
that concludes a merge may bring such a file when the staged file is identical, in content and
mode, to that file on a branch being merged. Any other change to it is still refused.
"""

import fnmatch
import os
import re
import subprocess
import sys
from pathlib import Path
from typing import Dict, List, Mapping, Optional, Tuple

MAX_FILE_BYTES = 10 * 1024 * 1024
MAX_NPZ_BYTES = 5 * 1024 * 1024
MAX_LINE = 100
FORBIDDEN_GLOBS = ["*.pt", "*.pkl", "*.joblib"]
ABS_PATH_RE = re.compile(r"/(?:Users|home)/[A-Za-z0-9_.-]+/")
CORRUPT_RE = re.compile(r"\bcorrupt(s|ed|ion|ing)?\b", re.IGNORECASE)  # noqa: realistic
CORRUPT_ALIASES = (  # noqa: realistic
    "corrupt_windows",
    "corrupt_dataset",
    "CorruptionConfig",
    "--corrupt-train",
    "--corrupt",
)
DIFF_OPTS = [
    "-c",
    "core.quotePath=false",
    "diff",
    "--no-ext-diff",
    "--src-prefix=a/",
    "--dst-prefix=b/",
]
NOQA = "noqa: realistic"
REALISTIC_MSG = 'say "realistic" instead of "corrupt*"'  # noqa: realistic

Added = Dict[str, List[Tuple[int, str]]]


def _git(repo: Path, *args: str) -> str:
    res = subprocess.run(["git", *args], cwd=repo, capture_output=True, text=True, check=True)
    return str(res.stdout)


def _strip_aliases(text: str) -> str:
    for alias in sorted(CORRUPT_ALIASES, key=len, reverse=True):
        text = text.replace(alias, " ")
    return text


def staged_files(repo: Path) -> List[Tuple[str, str]]:
    """(status, path) of added/copied/modified/renamed staged files."""
    raw = _git(repo, *DIFF_OPTS, "--cached", "--name-status", "--diff-filter=ACMR", "-z")
    parts = raw.split("\0")
    files: List[Tuple[str, str]] = []
    i = 0
    while i < len(parts) and parts[i]:
        status = parts[i][0]
        if status in "RC":
            files.append((status, parts[i + 2]))
            i += 3
        else:
            files.append((status, parts[i + 1]))
            i += 2
    return files


def added_lines(repo: Path) -> Added:
    """Added lines per staged text file as (new line number, text)."""
    raw = _git(repo, *DIFF_OPTS, "--cached", "-U0", "--diff-filter=ACMR", "--no-color")
    out: Added = {}
    path: Optional[str] = None
    lineno = 0
    for line in raw.splitlines():
        if line.startswith("+++ "):
            path = line[6:].rstrip("\t") if line.startswith("+++ b/") else None
        elif line.startswith("@@"):
            m = re.match(r"@@ -\d+(?:,\d+)? \+(\d+)", line)
            lineno = int(m.group(1)) if m else 0
        elif line.startswith("+") and path is not None:
            out.setdefault(path, []).append((lineno, line[1:]))
            lineno += 1
    return out


def merge_heads(repo: Path) -> List[str]:
    """Commits being merged into HEAD (empty when no merge is in progress)."""
    path = Path(_git(repo, "rev-parse", "--git-path", "MERGE_HEAD").strip())
    if not path.is_absolute():
        path = repo / path
    return path.read_text().split() if path.is_file() else []


def _mode_and_blob(repo: Path, path: str, head: Optional[str] = None) -> Optional[Tuple[str, str]]:
    """(mode, blob id) of ``path`` in the index, or in commit ``head``; None if it is not there.

    ``path`` is taken literally (no pathspec magic such as ``:/`` or globs), and the entry git
    returns must carry exactly that path.
    """
    args = ["ls-files", "-s", "-z"] if head is None else ["ls-tree", "-z", head]
    res = subprocess.run(
        ["git", "--literal-pathspecs", *args, "--", path],
        cwd=repo,
        capture_output=True,
        text=True,
    )
    entries = [e for e in res.stdout.split("\0") if e]
    if res.returncode != 0 or len(entries) != 1:
        return None
    meta, _, name = entries[0].partition("\t")
    if name != path:
        return None
    fields = meta.split()
    return (fields[0], fields[1] if head is None else fields[2])


def same_as_merged(repo: Path, path: str, heads: List[str]) -> bool:
    """True when the staged file has the content and mode of ``path`` on a branch being merged."""
    staged = _mode_and_blob(repo, path)
    if staged is None:
        return False
    return any(_mode_and_blob(repo, path, h) == staged for h in heads)


def check(repo: Path, env: Optional[Mapping[str, str]] = None) -> List[str]:
    env = os.environ if env is None else env
    problems: List[str] = []
    files = staged_files(repo)
    added = added_lines(repo)
    heads = merge_heads(repo)

    for status, path in files:
        name = os.path.basename(path)
        parts = path.split("/")
        for g in FORBIDDEN_GLOBS:
            if fnmatch.fnmatch(name, g):
                problems.append(f"{path}:1: forbidden file type {g}")
        if "cache" in parts[:-1]:
            problems.append(f"{path}:1: file under a cache/ directory")
        if env.get("ALLOW_CLAUDE_CONFIG") != "1" and (name == "CLAUDE.md" or parts[0] == ".claude"):
            if not same_as_merged(repo, path, heads):
                problems.append(f"{path}:1: CLAUDE.md/.claude change (set ALLOW_CLAUDE_CONFIG=1)")
        size = len(
            subprocess.run(
                ["git", "cat-file", "blob", f":{path}"], cwd=repo, capture_output=True, check=True
            ).stdout
        )
        if env.get("ALLOW_LARGE") != "1":
            if size > MAX_FILE_BYTES:
                problems.append(f"{path}:1: file over 10 MB (set ALLOW_LARGE=1)")
            elif name.endswith(".npz") and size > MAX_NPZ_BYTES:
                problems.append(f"{path}:1: .npz over 5 MB (set ALLOW_LARGE=1)")

    for path, lines in added.items():
        for n, text in lines:
            if env.get("ALLOW_ABS_PATHS") != "1" and ABS_PATH_RE.search(text):
                problems.append(f"{path}:{n}: absolute home path (set ALLOW_ABS_PATHS=1)")
            if NOQA not in text and CORRUPT_RE.search(_strip_aliases(text)):
                problems.append(f"{path}:{n}: {REALISTIC_MSG}")
            if path.startswith("icenine_py/") and path.endswith(".py") and len(text) > MAX_LINE:
                problems.append(f"{path}:{n}: line longer than {MAX_LINE} characters")

    for status, path in files:
        if status == "A" and path.startswith("icenine_py/") and path.endswith(".py"):
            content = subprocess.run(
                ["git", "cat-file", "blob", f":{path}"], cwd=repo, capture_output=True, check=True
            ).stdout
            res = subprocess.run(
                [sys.executable, "-m", "black", "--check", "-q", "-l", "100", "-"],
                input=content,
                capture_output=True,
            )
            if res.returncode != 0 and b"No module named black" not in res.stderr:
                problems.append(f"{path}:1: new file fails black --check -l 100")
    return problems


def main() -> int:
    repo = Path(_git(Path.cwd(), "rev-parse", "--show-toplevel").strip())
    problems = check(repo)
    for p in problems:
        print(p)
    return 1 if problems else 0


if __name__ == "__main__":
    sys.exit(main())
