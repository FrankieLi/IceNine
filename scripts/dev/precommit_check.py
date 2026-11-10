"""Pre-commit checks on staged files and added lines (stdlib only; black run via subprocess).

Prints one line per problem, ``file:line: reason``, and exits 1; silent with exit 0 when clean.
Overrides (environment): ALLOW_CLAUDE_CONFIG=1, ALLOW_LARGE=1, ALLOW_ABS_PATHS=1.

CLAUDE.md/.claude files (any letter case, for example ``.Claude/``): every staged change is
refused, including a deletion, a rename away and a change of type (file to symlink). There is
one exception (owner decisions D-20261109-N11 and N12). A commit that concludes a merge may bring
such a file as it is on a branch being merged, when that branch's commit is on a remote-tracking
branch (so a stash or another local-only commit does not count):

- an added or changed file passes when the staged file is identical, in content and mode, to
  that file on the merged branch;
- a deletion passes when the merged branch does not have the file and the merge base has it.

The hook ``pre-merge-commit`` runs the same checks for a merge that git commits by itself. git has
not recorded the merged branch at that point, so the exception cannot apply there: such a merge
is left staged and is concluded with ``git commit``, where it does apply.
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


def is_guarded(path: str) -> bool:
    """True for CLAUDE.md at any depth and for anything under a top-level .claude/ (any case)."""
    parts = path.split("/")
    return parts[-1].lower() == "claude.md" or parts[0].lower() == ".claude"


def guarded_changes(repo: Path) -> List[Tuple[str, str]]:
    """(status, path) of every staged change to a guarded file; status is A, M, D or T.

    Renames and copies are not paired here, so a rename shows as a deletion plus an addition.
    """
    raw = _git(
        repo, *DIFF_OPTS, "--cached", "--name-status", "--no-renames", "--diff-filter=AMDT", "-z"
    )
    parts = raw.split("\0")
    pairs = [(parts[i][:1], parts[i + 1]) for i in range(0, len(parts) - 1, 2)]
    return [(status, path) for status, path in pairs if is_guarded(path)]


def published_heads(repo: Path, heads: List[str]) -> List[str]:
    """The merged commits that are on a remote-tracking branch."""
    out = []
    for h in heads:
        res = subprocess.run(
            ["git", "for-each-ref", "--count=1", "--contains", h, "refs/remotes"],
            cwd=repo,
            capture_output=True,
            text=True,
        )
        if res.returncode == 0 and res.stdout.strip():
            out.append(h)
    return out


def same_as_merged(repo: Path, path: str, heads: List[str]) -> bool:
    """True when the staged file has the content and mode of ``path`` on a branch being merged."""
    staged = _mode_and_blob(repo, path)
    if staged is None:
        return False
    return any(_mode_and_blob(repo, path, h) == staged for h in heads)


def deleted_on_merged(repo: Path, path: str, heads: List[str]) -> bool:
    """True when a branch being merged lacks ``path`` and every merge base with it has the file."""
    for h in heads:
        if _mode_and_blob(repo, path, h) is not None:
            continue
        res = subprocess.run(
            ["git", "merge-base", "--all", "HEAD", h], cwd=repo, capture_output=True, text=True
        )
        bases = res.stdout.split() if res.returncode == 0 else []
        if bases and all(_mode_and_blob(repo, path, b) is not None for b in bases):
            return True
    return False


def guarded_problems(repo: Path, env: Mapping[str, str]) -> List[str]:
    """One line per staged change to a CLAUDE.md/.claude file that the merge exception does not cover."""
    if env.get("ALLOW_CLAUDE_CONFIG") == "1":
        return []
    changes = guarded_changes(repo)
    if not changes:
        return []
    heads = published_heads(repo, merge_heads(repo))
    problems = []
    for status, path in changes:
        if status == "D":
            if not deleted_on_merged(repo, path, heads):
                problems.append(f"{path}:1: CLAUDE.md/.claude deletion (set ALLOW_CLAUDE_CONFIG=1)")
        elif not same_as_merged(repo, path, heads):
            problems.append(f"{path}:1: CLAUDE.md/.claude change (set ALLOW_CLAUDE_CONFIG=1)")
    return problems


def check(repo: Path, env: Optional[Mapping[str, str]] = None) -> List[str]:
    env = os.environ if env is None else env
    problems: List[str] = guarded_problems(repo, env)
    files = staged_files(repo)
    added = added_lines(repo)

    for status, path in files:
        name = os.path.basename(path)
        parts = path.split("/")
        for g in FORBIDDEN_GLOBS:
            if fnmatch.fnmatch(name, g):
                problems.append(f"{path}:1: forbidden file type {g}")
        if "cache" in parts[:-1]:
            problems.append(f"{path}:1: file under a cache/ directory")
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
