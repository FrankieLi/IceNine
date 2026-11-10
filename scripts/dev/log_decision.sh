#!/bin/bash
# Usage: log_decision.sh --row "<track/row>" --title "<one-line decision>" --decision "<text>"
#                        --where "<ROADMAP row, plan or doc>" [--supersedes <ID>]
#                        [--date YYYY-MM-DD] [--agent I|N] [--recorder "<name>"]
#                        [--trailer "<Co-Authored-By line>"] [--no-push] [--no-fetch] [--dry-run]
# Exit codes: 0 recorded; 1 refused or failed with nothing recorded; 2 bad usage; 3 partial
# failure (the commit was made, but shmem or the push failed; the command to rerun is printed).
# Append one owner decision to icenine_py/docs/DECISIONS.md on the target feature branch, commit it
# there without switching the caller's checkout, save it to shmem, and push the branch.
#
# Target branch: the caller's branch, followed up its recorded parents (git config
# branch.<b>.parent) to the first branch with no recorded parent; that branch must be feature/*.
# A chain that ends at develop, master or main, at a non-feature/* branch, or a detached HEAD is
# refused (no guessing); a feature/* branch with no recorded parent is itself the target.
# ID: D-<YYYYMMDD>-<agent><n>, n = 1 + the largest n on that date and agent in the target's
# DECISIONS.md, else 1. The agent letter (I or N) comes from --agent or `git config --local
# icenine.agent`; with neither the script refuses (it never defaults to I). The recorder name
# follows the letter (I = "Ilium Agent", N = "Newt Agent") unless --recorder is given.
# --supersedes must be an ID (D-YYYYMMDD-I<n> or -N<n>) that is a heading in the target's log.
# Commit trailer: --trailer or LOG_DECISION_TRAILER; required for a real run (not --dry-run).
# Origin check: origin is fetched and must not be ahead of the target. If origin is unreachable
# the script refuses; --no-fetch skips the check (with a warning). If origin has no such branch,
# or there is no origin, the check is skipped with a note.
# Commit: in the worktree that holds the target branch (refused if its DECISIONS.md is dirty), else
# in a temporary worktree removed on exit. Only DECISIONS.md is staged; the pre-commit hook runs.
# A failed shmem or push leaves the commit in place, prints the command to rerun and exits 3.
# Overrides for tests: LOG_DECISION_SHMEM_CMD (default shmem), LOG_DECISION_DOC (default
# icenine_py/docs/DECISIONS.md, relative to the repository root).
set -euo pipefail

usage() { sed -n '2,7p' "$0" | sed 's/^# \{0,1\}//' >&2; }
die() { echo "error: $*" >&2; exit 1; }
bad() { echo "error: $*" >&2; usage; exit 2; }

ROW="" TITLE="" DECISION="" WHERE="" SUPERSEDES="none" DATE="" AGENT="" RECORDER=""
TRAILER="${LOG_DECISION_TRAILER:-}"
PUSH=1 FETCH=1 DRY=0
while [ $# -gt 0 ]; do
    case "$1" in
        --row) shift; ROW="${1:-}" ;;
        --title) shift; TITLE="${1:-}" ;;
        --decision) shift; DECISION="${1:-}" ;;
        --where) shift; WHERE="${1:-}" ;;
        --supersedes) shift; SUPERSEDES="${1:-}" ;;
        --date) shift; DATE="${1:-}" ;;
        --agent) shift; AGENT="${1:-}" ;;
        --recorder) shift; RECORDER="${1:-}" ;;
        --trailer) shift; TRAILER="${1:-}" ;;
        --no-push) PUSH=0 ;;
        --no-fetch) FETCH=0 ;;
        --dry-run) DRY=1 ;;
        -h | --help) usage; exit 0 ;;
        *) bad "unknown option: $1" ;;
    esac
    shift
done

for pair in "row:$ROW" "title:$TITLE" "decision:$DECISION" "where:$WHERE"; do
    [ -n "${pair#*:}" ] || bad "--${pair%%:*} is required"
done
for pair in "row:$ROW" "title:$TITLE" "decision:$DECISION" "where:$WHERE" "supersedes:$SUPERSEDES" \
    "recorder:$RECORDER"; do
    case "${pair#*:}" in
        *$'\n'* | *$'\r'*) bad "--${pair%%:*} must not contain a newline" ;;
    esac
done
if [ "$SUPERSEDES" != "none" ]; then
    echo "$SUPERSEDES" | grep -qE '^D-[0-9]{8}-[IN][0-9]+$' ||
        bad "--supersedes must be an ID such as D-20261018-I1 (or none)"
fi
[ -n "$DATE" ] || DATE=$(date +%Y-%m-%d)
echo "$DATE" | grep -qE '^[0-9]{4}-[0-9]{2}-[0-9]{2}$' || bad "--date must be YYYY-MM-DD"
DATE8="${DATE//-/}"

git rev-parse --show-toplevel >/dev/null 2>&1 || die "not inside a git repository"
ROOT=$(git rev-parse --show-toplevel)
cd "$ROOT"

if [ -z "$AGENT" ]; then
    AGENT=$(git config --local --get icenine.agent || true)
fi
case "$AGENT" in
    I) [ -n "$RECORDER" ] || RECORDER="Ilium Agent" ;;
    N) [ -n "$RECORDER" ] || RECORDER="Newt Agent" ;;
    *) bad "agent unknown: pass --agent I|N or run 'git config --local icenine.agent I' (or N);" \
        "there is no default" ;;
esac
if [ $DRY -ne 1 ] && [ -z "$TRAILER" ]; then
    bad "no commit trailer: pass --trailer or set LOG_DECISION_TRAILER"
fi

# 1. target branch
CUR=$(git symbolic-ref -q --short HEAD || true)
[ -n "$CUR" ] || die "detached HEAD; check out a task or feature branch"
TARGET="$CUR"
hops=0
while :; do
    P=$(git config "branch.$TARGET.parent" || true)
    [ -n "$P" ] || break
    TARGET="$P"
    hops=$((hops + 1))
    [ $hops -le 20 ] || die "parent chain of '$CUR' is too long or loops"
done
case "$TARGET" in
    develop | master | main) die "target branch '$TARGET' is not allowed" ;;
    feature/*) ;;
    *) die "'$CUR' leads to '$TARGET', not a feature/* branch (no guessing)" ;;
esac
git rev-parse --verify -q "refs/heads/$TARGET" >/dev/null ||
    die "target branch '$TARGET' does not exist"

# 2. lock
COMMON=$(git rev-parse --git-common-dir)
LOCK="$COMMON/log_decision.lock"
TMPWT=""
cleanup() {
    if [ -n "$TMPWT" ] && [ -d "$TMPWT" ]; then
        git worktree remove --force "$TMPWT" >/dev/null 2>&1 || true
        rm -rf "$TMPWT"
        git worktree prune >/dev/null 2>&1 || true
    fi
    rm -f "$LOCK/owner" 2>/dev/null || true
    rmdir "$LOCK" 2>/dev/null || true
}
if ! mkdir "$LOCK" 2>/dev/null; then
    HELD=$(cat "$LOCK/owner" 2>/dev/null || echo "owner unknown")
    die "another log_decision.sh holds the lock ($LOCK; $HELD); rerun when it ends, or remove" \
        "it with rmdir if that process is gone"
fi
trap cleanup EXIT
echo "pid $$, started $(date '+%Y-%m-%d %H:%M:%S')" >"$LOCK/owner"

# 3. fetch; origin must not be ahead
wt_of() {  # print the worktree path holding refs/heads/$1, if any
    git worktree list --porcelain |
        awk -v b="refs/heads/$1" '/^worktree /{w=substr($0,10)} $1=="branch" && $2==b {print w}'
}
HOLDER=$(wt_of "$TARGET")
if ! git remote get-url origin >/dev/null 2>&1; then
    echo "note: no remote 'origin'; skipping the fetch check" >&2
elif [ $FETCH -eq 0 ]; then
    echo "warning: --no-fetch: not checking whether origin/$TARGET is ahead" >&2
else
    RC=0
    git ls-remote --exit-code --heads origin "$TARGET" >/dev/null 2>&1 || RC=$?
    case $RC in
        0)
            git fetch -q origin "+refs/heads/$TARGET:refs/remotes/origin/$TARGET" ||
                die "git fetch origin $TARGET failed; rerun with network, or pass --no-fetch"
            BEHIND=$(git rev-list --count "refs/heads/$TARGET..refs/remotes/origin/$TARGET")
            if [ "$BEHIND" -gt 0 ]; then
                if [ -n "$HOLDER" ]; then
                    HINT="git -C $HOLDER merge origin/$TARGET"
                else
                    HINT="git fetch origin $TARGET:$TARGET (fast-forward; if it refuses, merge"
                    HINT="$HINT origin/$TARGET in a checkout of $TARGET)"
                fi
                die "origin/$TARGET has $BEHIND commit(s) the local $TARGET lacks;" \
                    "bring them in first, then rerun: $HINT"
            fi
            ;;
        2) echo "note: origin has no $TARGET; skipping the ahead check" >&2 ;;
        *) die "cannot reach origin (git ls-remote exit $RC); rerun with network," \
               "or pass --no-fetch to skip the check" ;;
    esac
fi

# 4. ID from the target tip
DOC="${LOG_DECISION_DOC:-icenine_py/docs/DECISIONS.md}"
TIP=$(git show "refs/heads/$TARGET:$DOC" 2>/dev/null) || die "$DOC not found on $TARGET"
LAST=$(printf '%s\n' "$TIP" | { grep -oE "^### D-$DATE8-$AGENT[0-9]+" || true; } |
    sed -E "s/^### D-$DATE8-$AGENT//" | sort -n | tail -1)
N=$(( ${LAST:-0} + 1 ))
ID="D-$DATE8-$AGENT$N"

# 5. validate
if [ "$SUPERSEDES" != "none" ]; then
    printf '%s\n' "$TIP" | S="$SUPERSEDES" awk '
        BEGIN { h = "### " ENVIRON["S"] }
        index($0, h) == 1 {
            r = substr($0, length(h) + 1)
            if (r == "" || substr(r, 1, 1) == " ") f = 1
        }
        END { exit !f }' ||
        die "--supersedes $SUPERSEDES is not an existing entry on $TARGET"
fi
DUP=$(printf '%s\n' "$TIP" | T="$TITLE" D="$DATE8" awk '
    index($0, "### D-" ENVIRON["D"] "-") == 1 {
        sep = " · "
        s = $0
        i = index(s, sep)
        if (i == 0) next
        s = substr(s, i + length(sep))
        j = index(s, sep)
        if (j == 0) next
        s = substr(s, j + length(sep))
        if (s == ENVIRON["T"]) { print $2; exit }
    }')
[ -z "$DUP" ] ||
    die "an entry on $DATE with this title exists already ($DUP); not recording it twice"

# 6. entry
strip_dot() { local s="$1"; s="${s%"${s##*[![:space:]]}"}"; printf '%s' "${s%.}"; }
wrap() {  # wrap <first prefix> <text>: 98 columns, two-space hanging indent
    printf '%s' "$2" | awk -v first="$1" '
        { for (i = 1; i <= NF; i++) w[++n] = $i }
        END {
            line = first
            fresh = 1
            for (i = 1; i <= n; i++) {
                if (fresh) { line = line w[i]; fresh = 0 }
                else if (length(line) + 1 + length(w[i]) <= 98) line = line " " w[i]
                else { print line; line = "  " w[i] }
            }
            print line
        }'
}
DEC=$(strip_dot "$DECISION")
WHR=$(strip_dot "$WHERE")
ENTRY="### $ID · $ROW · $TITLE
- $DATE, owner. Recorded by $RECORDER. Supersedes: $SUPERSEDES.
$(wrap "- " "$DEC. Where: $WHR.")"
SHMEM_LINE="Owner decision $ID ($DATE, IceNine $ROW): $(strip_dot "$TITLE"). $DEC."

if [ $DRY -eq 1 ]; then
    echo "target: $TARGET"
    echo "id: $ID"
    echo "entry:"
    printf '%s\n' "$ENTRY"
    echo "shmem: $SHMEM_LINE"
    exit 0
fi

# 7. commit on the target branch
if [ -n "$HOLDER" ]; then
    WT="$HOLDER"
    if [ -n "$(git -C "$WT" status --porcelain -- "$DOC")" ]; then
        die "uncommitted changes to $DOC in $WT; commit or discard them first"
    fi
else
    TMPWT=$(mktemp -d "${TMPDIR:-/tmp}/log_decision.XXXXXX")
    rmdir "$TMPWT"
    git worktree add -q "$TMPWT" "$TARGET" || die "could not create a worktree of $TARGET"
    WT="$TMPWT"
fi
FILE="$WT/$DOC"
if [ -n "$(tail -c1 "$FILE")" ]; then printf '\n' >>"$FILE"; fi
printf '\n%s\n' "$ENTRY" >>"$FILE"
MSG="DECISIONS $ID: $TITLE

$TRAILER"
if ! git -C "$WT" add -- "$DOC" || ! git -C "$WT" commit -q -m "$MSG" -- "$DOC"; then
    git -C "$WT" reset -q -- "$DOC" || true
    git -C "$WT" checkout -- "$DOC" || true
    die "commit on $TARGET failed (entry removed again)"
fi
HASH=$(git -C "$WT" rev-parse --short HEAD)
FAILED=0
echo "ok: committed $ID on $TARGET as $HASH"

# 8. shmem
SHMEM_CMD="${LOG_DECISION_SHMEM_CMD:-shmem}"
if bash -c "$SHMEM_CMD"' add --source user "$1"' _ "$SHMEM_LINE" >/dev/null; then
    echo "ok: shmem"
else
    printf 'FAILED: %s add --source user %q\n' "$SHMEM_CMD" "$SHMEM_LINE"
    FAILED=3
fi

# 9. push
if [ $PUSH -eq 1 ]; then
    if git push -q origin "refs/heads/$TARGET:refs/heads/$TARGET"; then
        echo "ok: pushed $TARGET"
    else
        echo "FAILED: git push origin $TARGET"
        FAILED=3
    fi
else
    echo "skipped: push (--no-push)"
fi

# 10. result
echo "$ID $HASH"
exit $FAILED
