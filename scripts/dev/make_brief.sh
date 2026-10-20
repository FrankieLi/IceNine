#!/bin/bash
# Usage: make_brief.sh --plan <plan.md> --task <task> [--section "<heading>"]
#                      [--fixes <review report>] [--trailer "<Co-Authored-By line>"]
#                      [--decisions <DECISIONS.md>]
# Print the brief for an implementer agent. The brief points at study-conventions instead of
# repeating it and takes its facts from the plan and the repository:
#   - the approval: the last DECISIONS.md entry that names the plan's file name (on a path
#     boundary) and says "approved" (not "not approved", "unapproved", ...); refuses when there is
#     none, or when a later entry names the plan and withdraws it or says
#     "Supersedes: <the approval's ID>";
#   - the process limit: `git config --local icenine.maxHeavyProcs` (refuses when unset or not a
#     positive integer); a plan paragraph starting "**Process limit exception:**" replaces the
#     standard line;
#   - the scope: a plan paragraph starting "**Library changes:**", else "no change under
#     `icenine_py/icenine/`";
#   - the stop rule: every plan paragraph starting "**Stop rule:**", else "stop and report if the
#     plan cannot be followed";
#   (a marker paragraph starts at column 0 or after a "- " bullet, runs to the next blank line,
#   and is ignored inside code fences; the whole paragraph is quoted)
#   - the commit trailer from --trailer or MAKE_BRIEF_TRAILER (refuses when neither is set).
# --plan is a path from the repository root or relative to icenine_py/ or icenine_py/docs/.
# --section must match exactly one heading outside code fences (an exact match wins, else one
# heading must contain the text); it only names the section in the brief. The approval and the
# three markers are read from the WHOLE plan, never from the section alone.
# With --fixes the brief asks for the required fixes of that review report on the task branch;
# the report may be in any worktree and is printed as an absolute, resolved path.
# MAKE_BRIEF_NOW (YYYYMMDD-HHMM) overrides the clock for the report name (for testing).
set -euo pipefail

usage() {
    echo "usage: make_brief.sh --plan <plan.md> --task <task> [--section <heading>]" \
        "[--fixes <review report>] [--trailer <line>] [--decisions <DECISIONS.md>]" >&2
    exit 2
}
die() { echo "error: $*" >&2; exit 1; }

# paragraphs <marker>: every paragraph of the plan that starts with the marker (at column 0, or
# after a "- " / "* " bullet marker, which is dropped), up to the next blank line, heading or
# (for a bullet) next bullet. Fenced code blocks are skipped. The whole plan is read, not a section.
paragraphs() {
    awk -v marker="$1" '
        /^[[:space:]]*(```|~~~)/ { fence = !fence; inrule = 0; next }
        fence { next }
        inrule && ($0 ~ /^[[:space:]]*$/ || $0 ~ /^#/ || (bullet && $0 ~ /^[-*] /)) { inrule = 0 }
        !inrule && (index($0, marker) == 1 ||
                    ($0 ~ /^[-*] / && index(substr($0, 3), marker) == 1)) {
            inrule = 1
            bullet = ($0 ~ /^[-*] /)
            if (bullet) $0 = substr($0, 3)
            print
            next
        }
        inrule { print }
    ' "$PLAN_REL"
}

PLAN=""
TASK=""
SECTION=""
FIXES=""
TRAILER="${MAKE_BRIEF_TRAILER:-}"
DECISIONS=""
while [ $# -gt 0 ]; do
    case "$1" in
        --plan) shift; PLAN="${1:-}" ;;
        --task) shift; TASK="${1:-}" ;;
        --section) shift; SECTION="${1:-}" ;;
        --fixes) shift; FIXES="${1:-}" ;;
        --trailer) shift; TRAILER="${1:-}" ;;
        --decisions) shift; DECISIONS="${1:-}" ;;
        *) echo "unknown option: $1" >&2; usage ;;
    esac
    shift || true
done
[ -n "$PLAN" ] && [ -n "$TASK" ] || usage
case "$TASK" in
    *[!A-Za-z0-9._-]* | -*)
        die "task name '$TASK' may only contain letters, digits, '.', '_', '-'" ;;
esac
[ -n "$TRAILER" ] || die "no commit trailer: pass --trailer or set MAKE_BRIEF_TRAILER"

ROOT=$(git rev-parse --show-toplevel)
ROOT=$(cd "$ROOT" && pwd -P)
cd "$ROOT"

# --- plan file, as a path relative to the repository root
PLAN_REL=""
case "$PLAN" in
    /*)
        [ -f "$PLAN" ] || die "plan not found: $PLAN"
        PLAN_ABS="$(cd "$(dirname "$PLAN")" && pwd -P)/$(basename "$PLAN")"
        case "$PLAN_ABS" in
            "$ROOT"/*) PLAN_REL="${PLAN_ABS#"$ROOT"/}" ;;
            *) die "plan is outside the repository: $PLAN" ;;
        esac
        ;;
    *)
        for cand in "$PLAN" "icenine_py/$PLAN" "icenine_py/docs/$PLAN"; do
            if [ -f "$cand" ]; then PLAN_REL="$cand"; break; fi
        done
        [ -n "$PLAN_REL" ] || die "plan not found: $PLAN" \
            "(tried the repository root, icenine_py/ and icenine_py/docs/)"
        ;;
esac
PLAN_REL="${PLAN_REL#./}"
PLAN_NAME=$(basename "$PLAN_REL")

if [ -n "$SECTION" ]; then
    # Headings outside code fences only. A heading that equals the text exactly wins; otherwise
    # exactly one heading must contain it ("1" would match both "1. One" and "11. Eleven").
    COUNTS=$(awk -v s="$SECTION" '
        /^[[:space:]]*(```|~~~)/ { fence = !fence; next }
        !fence && /^#+ / {
            h = $0; sub(/^#+ +/, "", h); sub(/[[:space:]]+$/, "", h)
            if (h == s) exact++
            if (index(h, s)) contains++
        }
        END { print exact + 0, contains + 0 }
    ' "$PLAN_REL")
    EXACT=${COUNTS% *}
    CONTAINS=${COUNTS#* }
    if [ "$EXACT" -eq 1 ] || { [ "$EXACT" -eq 0 ] && [ "$CONTAINS" -eq 1 ]; }; then
        :
    elif [ "$CONTAINS" -eq 0 ]; then
        die "no heading containing '$SECTION' in $PLAN_REL"
    else
        die "'$SECTION' matches $CONTAINS headings in $PLAN_REL; give more of the heading"
    fi
fi
if [ -n "$FIXES" ]; then
    # Review reports live in .claude/reports/ of the session that wrote them (gitignored, one per
    # worktree), so any existing file is accepted and printed as an absolute, resolved path.
    [ -f "$FIXES" ] || die "review report not found: $FIXES"
    FIXES="$(cd "$(dirname "$FIXES")" && pwd -P)/$(basename "$FIXES")"
fi

# --- approval
[ -n "$DECISIONS" ] || DECISIONS="icenine_py/docs/DECISIONS.md"
[ -f "$DECISIONS" ] || die "decision log not found: $DECISIONS"
# The plan's file name counts only on a path boundary (start of a line, "/", a backtick, a space,
# a quote or "(" before it), so `old_todo_demo.md` does not name `todo_demo.md`. An entry approves
# when its lowercased text still contains "approved" after "not approved", "not yet approved",
# "never approved", "unapproved" and "disapproved" are removed. The last approving entry is the
# approval; it is refused when a later entry names the plan and withdraws it (a negated
# approval, or "withdrawn", "revoked", "rescinded", "cancelled", "rewrite", "on hold") or says
# "Supersedes: <the approval's ID>" (the log is append-only). A later entry that only mentions
# the plan (a cost decision, a note) does not cancel the approval.
RESULT=$(awk -v name="$PLAN_NAME" -v sq="'" '
    function names_plan(s,    i, before, after, n) {
        n = length(name)
        while ((i = index(s, name)) > 0) {
            before = (i == 1) ? "\n" : substr(s, i - 1, 1)
            after = substr(s, i + n, 1)
            if (index("\n/` \"" sq "(", before) > 0 && after !~ /[A-Za-z0-9_-]/) return 1
            s = substr(s, i + 1)
        }
        return 0
    }
    function withdraws(s) {
        s = tolower(s)
        return s ~ /(not|never)( +yet)?( +been)? +approved|(un|dis)approved|withdr[ae]w|withdrawn/ \
            || s ~ /revok|rescind|retract|cancel|rewrite|on hold/
    }
    function approves(s) {
        s = tolower(s)
        gsub(/(not|never)( +yet)?( +been)? +approved/, "", s)
        gsub(/(un|dis)approved/, "", s)
        return s ~ /approved/
    }
    function flush(    id) {
        if (head == "") return
        id = head; sub(/ .*/, "", id)
        if (names_plan(body) && approves(body)) {
            approved_head = head; approved_id = id; refused = ""
        } else if (approved_head != "" && names_plan(body) && withdraws(body)) {
            refused = head
        } else if (approved_head != "" && index(body, approved_id) > 0 \
                   && tolower(body) ~ /supersede/) {
            refused = head
        }
    }
    /^### / { flush(); head = substr($0, 5); body = $0 "\n"; next }
    head != "" { body = body $0 "\n" }
    END {
        flush()
        if (refused != "") print "REFUSED\t" refused "\t" approved_head
        else print "OK\t" approved_head
    }
' "$DECISIONS")
APPROVAL=${RESULT#*$'\t'}
if [ "${RESULT%%$'\t'*}" = "REFUSED" ]; then
    LATER=${APPROVAL%%$'\t'*}
    EARLIER=${APPROVAL#*$'\t'}
    die "the approval '$EARLIER' is followed by '$LATER', which names $PLAN_NAME and withdraws it" \
        "or supersedes the approval; the plan needs a new owner approval"
fi
[ -n "$APPROVAL" ] || die "no entry in $DECISIONS names $PLAN_NAME as approved;" \
    "the plan needs an owner approval first"

# --- process limit
MAXPROCS=$(git config --local --get icenine.maxHeavyProcs || true)
[ -n "$MAXPROCS" ] || die "icenine.maxHeavyProcs is not set;" \
    "run: git config --local icenine.maxHeavyProcs <n>"
case "$MAXPROCS" in
    *[!0-9]* | 0 | 0[0-9]*)
        die "icenine.maxHeavyProcs must be a positive integer, got '$MAXPROCS'" ;;
esac
EXCEPTION=$(paragraphs '**Process limit exception:**')
if [ -n "$EXCEPTION" ]; then
    COMPUTE="Compute (plan exception to the machine limit of $MAXPROCS heavy processes): $EXCEPTION"
else
    COMPUTE="Compute: at most $MAXPROCS heavy processes at once on this machine."
fi

# --- scope
LIBRARY=$(paragraphs '**Library changes:**')
if [ -n "$LIBRARY" ]; then
    SCOPE="Scope: $LIBRARY"
else
    SCOPE='Scope: no change under `icenine_py/icenine/`.'
fi

# --- stop rule (every paragraph starting with the marker, up to the next blank line)
STOP=$(paragraphs '**Stop rule:**' | sed 's/^/> /')

# --- branch
CUR=$(git symbolic-ref --short -q HEAD || true)
[ -n "$CUR" ] || die "detached HEAD; check out a feature/* branch"
case "$CUR" in
    feature/*) ;;
    *) die "current branch '$CUR' is not a feature/* branch" ;;
esac
CUR_PARENT=$(git config "branch.$CUR.parent" || true)
if [ -n "$FIXES" ]; then
    if [ -n "$CUR_PARENT" ] && [ "${CUR%-"$TASK"}" != "$CUR" ]; then
        BRANCH_LINE="Branch: work on the existing task branch \`$CUR\` (parent \`$CUR_PARENT\`);"
        BRANCH_LINE="$BRANCH_LINE add new commits there."
    elif git rev-parse --verify -q "refs/heads/$CUR-$TASK" >/dev/null; then
        BRANCH_LINE="Branch: work on the existing task branch \`$CUR-$TASK\` (parent \`$CUR\`);"
        BRANCH_LINE="$BRANCH_LINE add new commits there."
    else
        BRANCH_LINE="Branch: the task is already merged, so run"
        BRANCH_LINE="$BRANCH_LINE \`scripts/dev/start_task.sh $TASK-fixes\` from \`$CUR\`"
        BRANCH_LINE="$BRANCH_LINE (branch \`$CUR-$TASK-fixes\`)."
    fi
else
    BRANCH_LINE="Branch: run \`scripts/dev/start_task.sh $TASK\` from \`$CUR\`"
    BRANCH_LINE="$BRANCH_LINE (branch \`$CUR-$TASK\`)."
fi

# --- report path
NOW="${MAKE_BRIEF_NOW:-$(date +%Y%m%d-%H%M)}"
case "$NOW" in
    [0-9][0-9][0-9][0-9][0-9][0-9][0-9][0-9]-[0-9][0-9][0-9][0-9]) ;;
    *) die "MAKE_BRIEF_NOW must look like YYYYMMDD-HHMM, got '$NOW'" ;;
esac
SUFFIX=""
[ -z "$FIXES" ] || SUFFIX="-fixes"
REPORT=".claude/reports/$NOW-$TASK$SUFFIX.md"

# --- print
echo "Follow study-conventions."
if [ -n "$FIXES" ]; then
    echo "Plan: $PLAN_REL${SECTION:+, section \"$SECTION\"}."
    echo "Apply the required fixes in $FIXES, follow study-conventions, re-run the affected" \
        "tests, and add the new commit hashes to your report. Name the fix numbers in your" \
        "commit subjects."
elif [ -n "$SECTION" ]; then
    echo "Implement section \"$SECTION\" of $PLAN_REL."
else
    echo "Implement the plan $PLAN_REL."
fi
echo "Approval: $APPROVAL"
echo "$BRANCH_LINE"
echo "$SCOPE"
echo 'No new dependency (`pyproject.toml`, `uv.lock` unchanged); no default flip; no golden change.'
if [ -n "$STOP" ]; then
    echo "Stop rule of the plan:"
    echo "$STOP"
else
    echo "Stop rule: stop and report if the plan cannot be followed."
fi
echo "$COMPUTE"
echo "Commits: stage by explicit path; end each commit message with:"
echo "$TRAILER"
echo "Report: write the full report to \`$REPORT\`; the final message is at most about 250 words."
