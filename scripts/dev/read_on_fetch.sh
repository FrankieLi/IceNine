#!/bin/bash
# Usage: read_on_fetch.sh [--agent I|N] [--other BRANCH] [--since REF] [--mark] [--save]
# ROADMAP protocol rule 2 ("read on fetch") as a read-only report on the other agent's branch:
# fetch origin, list the commits and changed files since the last-read point, diff the other
# agent's ROADMAP mailbox, list the other agent's decisions that are not in shmem (with ready
# `shmem add` lines), and print the merge command. Nothing is merged.
#
# The agent letter (I = Ilium Agent, N = Newt Agent) comes from `git config --local icenine.agent`
# or --agent; the script refuses to run when neither is set (it never defaults to one agent).
# For I the other branch is feature/eval-reduction, the mailbox "From Newt Agent" and the decisions
# D-*-N*; for N it is feature/python-speed-parity, "From Ilium Agent" and D-*-I*. --other overrides
# the branch.
#
# Last-read point: --since REF, else `git config --local icenine.<other>.readOnFetch` (stored by
# --mark; one value per clone, so all worktrees share it), else the merge base of HEAD and
# origin/<other>. --mark stores origin/<other>'s hash as the new point once the run succeeded;
# without it the point does not move, so an interrupted session reads the same entries again.
#
# Decisions: an entry is "saved" when `shmem query "Owner decision <ID>" --source user` returns a
# line containing "Owner decision <ID>". The printed `shmem add` lines are only run with --save.
# An unreachable shmem (command missing or failing) is reported as such, never read as "not saved";
# the exit status is then 3. The shmem command can be overridden with READ_ON_FETCH_SHMEM_CMD.
set -euo pipefail
AGENT=""
OTHER=""
SINCE=""
MARK=0
SAVE=0
while [ $# -gt 0 ]; do
    case "$1" in
        --agent) shift; AGENT="${1:-}" ;;
        --other) shift; OTHER="${1:-}" ;;
        --since) shift; SINCE="${1:-}" ;;
        --mark) MARK=1 ;;
        --save) SAVE=1 ;;
        *) echo "unknown option: $1" >&2; exit 2 ;;
    esac
    shift
done
ROOT=$(git rev-parse --show-toplevel)
cd "$ROOT"
[ -n "$AGENT" ] || AGENT=$(git config --local --get icenine.agent || true)
case "$AGENT" in
    I) DEFAULT_OTHER="feature/eval-reduction"; MAILBOX="From Newt Agent"; OTHER_LETTER="N" ;;
    N) DEFAULT_OTHER="feature/python-speed-parity"; MAILBOX="From Ilium Agent"; OTHER_LETTER="I" ;;
    "")
        echo "error: no agent letter; run 'git config --local icenine.agent I|N' or pass --agent" >&2
        exit 1 ;;
    *) echo "error: agent must be I or N, got '$AGENT'" >&2; exit 1 ;;
esac
[ -n "$OTHER" ] || OTHER="$DEFAULT_OTHER"
SHMEM_CMD="${READ_ON_FETCH_SHMEM_CMD:-shmem}"
ROADMAP="icenine_py/docs/ROADMAP.md"
DECISIONS="icenine_py/docs/DECISIONS.md"
KEY="icenine.$OTHER.readOnFetch"

git fetch -q origin || { echo "error: git fetch origin failed" >&2; exit 1; }
REF="origin/$OTHER"
if ! git rev-parse --verify -q "$REF^{commit}" >/dev/null; then
    echo "error: $REF does not exist after fetch" >&2
    exit 1
fi
TIP=$(git rev-parse "$REF")
if [ -n "$SINCE" ]; then
    LAST=$(git rev-parse --verify -q "$SINCE^{commit}") || {
        echo "error: --since '$SINCE' is not a commit" >&2; exit 1; }
    FROM="--since"
elif LAST=$(git config --local --get "$KEY") && [ -n "$LAST" ]; then
    git rev-parse --verify -q "$LAST^{commit}" >/dev/null || {
        echo "error: stored read point '$LAST' ($KEY) is not a commit" >&2; exit 1; }
    FROM="stored read point"
else
    LAST=$(git merge-base HEAD "$REF" 2>/dev/null) || {
        echo "error: no merge base of HEAD and $REF; pass --since" >&2; exit 1; }
    FROM="merge base with HEAD"
fi
echo "agent $AGENT; other branch $REF at $(git rev-parse --short "$TIP")"
echo "last read: $(git rev-parse --short "$LAST") ($FROM)"

echo
echo "== commits on $REF since the last read, not on your branch =="
# Commits of your own that the other agent merged are already on HEAD; leave them out.
NEW=$(git log --oneline "$LAST..$REF" --not HEAD)
ALL_N=$(git rev-list --count "$LAST..$REF")
if [ -z "$NEW" ]; then
    echo "none"
    [ "$ALL_N" -eq 0 ] || echo "($ALL_N commit(s) since the last read are already on your branch)"
else
    N=$(printf '%s\n' "$NEW" | wc -l | tr -d ' ')
    if [ "$N" -gt 40 ]; then
        printf '%s\n' "$NEW" | head -40
        echo "... ($((N - 40)) more commits)"
    else
        printf '%s\n' "$NEW"
    fi
    OWN_N=$((ALL_N - N))
    [ "$OWN_N" -eq 0 ] || echo "($OWN_N commit(s) already on your branch left out)"
    BASE=$(git merge-base HEAD "$REF" || echo "$LAST")
    echo "-- files a merge would bring ($(git rev-parse --short "$BASE") to $(git rev-parse --short "$TIP")):"
    STAT=$(git diff --stat "$BASE" "$REF")
    NS=$(printf '%s\n' "$STAT" | wc -l | tr -d ' ')
    if [ "$NS" -gt 40 ]; then
        printf '%s\n' "$STAT" | head -39
        echo "... ($((NS - 39)) more lines; the last line of the full stat follows)"
        printf '%s\n' "$STAT" | tail -1
    else
        printf '%s\n' "$STAT"
    fi
fi

# the text of one "## <name>" subsection of ROADMAP.md at a revision (empty if absent)
mailbox_at() {  # $1 revision
    git show "$1:$ROADMAP" 2>/dev/null | awk -v h="## $MAILBOX" '
        $0 == h { on = 1; print; next }
        on && /^#/ { exit }
        on { print }' || true
}
echo
echo "== mailbox \"$MAILBOX\" in $ROADMAP (diff since the last read) =="
OLD_MB=$(mktemp); NEW_MB=$(mktemp)
trap 'rm -f "$OLD_MB" "$NEW_MB"' EXIT
mailbox_at "$LAST" > "$OLD_MB"
mailbox_at "$REF" > "$NEW_MB"
if [ ! -s "$NEW_MB" ] && [ ! -s "$OLD_MB" ]; then
    echo "no \"$MAILBOX\" subsection found"
elif cmp -s "$OLD_MB" "$NEW_MB"; then
    echo "unchanged"
else
    D=$(diff -u --label "mailbox at last read" --label "mailbox at $REF" "$OLD_MB" "$NEW_MB" || true)
    ND=$(printf '%s\n' "$D" | wc -l | tr -d ' ')
    printf '%s\n' "$D" | head -80
    [ "$ND" -le 80 ] || echo "... ($((ND - 80)) more diff lines; read $ROADMAP at $REF)"
fi

echo
echo "== decisions (D-*-$OTHER_LETTER*) on $REF and shmem =="
# one record per entry: ID, date, track/row, title, text; fields separated by \x1f
US=$'\x1f'
ENTRIES=$(git show "$REF:$DECISIONS" 2>/dev/null | awk -v L="$OTHER_LETTER" -v US="$US" '
    function flush() {
        if (id != "") print id US date US row US title US text
        id = ""
    }
    $0 ~ ("^### D-[0-9]+-" L "[0-9]+ ") {
        flush()
        line = $0; sub(/^### /, "", line)
        p1 = index(line, " · ")
        id = substr(line, 1, p1 - 1)
        rest = substr(line, p1 + length(" · "))
        p2 = index(rest, " · ")
        row = substr(rest, 1, p2 - 1)
        title = substr(rest, p2 + length(" · "))
        date = ""; text = ""; n = 0
        next
    }
    /^### / { flush(); next }
    id != "" && $0 !~ /^[ \t]*$/ {
        n++
        t = $0; gsub(/\t/, " ", t); sub(/^[ \t]*-?[ \t]*/, "", t)
        if (n == 1) { if (match(t, /[0-9][0-9][0-9][0-9]-[0-9][0-9]-[0-9][0-9]/)) date = substr(t, RSTART, RLENGTH) }
        else text = (text == "" ? t : text " " t)
    }
    END { flush() }' || true)
UNREACHABLE=0
SAVE_FAILED=0
LINES=""
COUNT=0
if [ -z "$ENTRIES" ]; then
    echo "no D-*-$OTHER_LETTER* entries on $REF"
else
    TOTAL=$(printf '%s\n' "$ENTRIES" | wc -l | tr -d ' ')
    while IFS="$US" read -r ID DATE ROW TITLE TEXT; do
        set +e
        OUT=$($SHMEM_CMD query "Owner decision $ID" --source user 2>/dev/null)
        RC=$?
        set -e
        if [ $RC -ne 0 ]; then
            UNREACHABLE=1
            echo "shmem unreachable (exit $RC from '$SHMEM_CMD query'): cannot tell which entries are saved"
            break
        fi
        if printf '%s\n' "$OUT" | grep -Eq "Owner decision ${ID}([^0-9]|\$)"; then
            continue
        fi
        FACT="Owner decision $ID ($DATE, IceNine $ROW): $TITLE. $TEXT"
        COUNT=$((COUNT + 1))
        QUOTED=$(printf '%s' "$FACT" | sed "s/'/'\\\\''/g")
        LINES="${LINES}${SHMEM_CMD} add --source user '${QUOTED}'"$'\n'
        if [ $SAVE -eq 1 ]; then
            set +e
            $SHMEM_CMD add --source user "$FACT" >/dev/null 2>&1
            RC=$?
            set -e
            if [ $RC -ne 0 ]; then
                SAVE_FAILED=1
                echo "FAILED to save $ID (exit $RC); rerun with --save"
            else
                echo "saved $ID"
            fi
        else
            echo "not in shmem: $ID"
        fi
    done <<EOF
$ENTRIES
EOF
    if [ $UNREACHABLE -eq 0 ]; then
        echo "$TOTAL entries on $REF, $COUNT not in shmem"
        if [ $COUNT -gt 0 ] && [ $SAVE -eq 0 ]; then
            echo "-- ready lines (run them, or rerun with --save):"
            printf '%s' "$LINES"
        fi
    fi
fi

echo
echo "== next step (not run) =="
echo "git merge $REF    # protocol rule 2; a conflict outside your own rows and mailbox needs the owner"

STATUS=0
if [ $UNREACHABLE -eq 1 ]; then
    echo "note: decisions not checked; read point not moved" >&2
    STATUS=3
elif [ $SAVE_FAILED -eq 1 ]; then
    echo "note: some decisions were not saved; read point not moved" >&2
    STATUS=1
elif [ $MARK -eq 1 ]; then
    git config --local "$KEY" "$TIP"
    echo "read point stored: $(git rev-parse --short "$TIP") ($KEY)"
fi
exit $STATUS
