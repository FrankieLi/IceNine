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
# heading must contain the text). The approval is read from the whole plan. With --section, the
# "**Library changes:**" and "**Process limit exception:**" markers are read from the section (its
# heading up to the next heading of the same or a higher level) when it has them, else from the
# plan-wide text: the plan outside every task-scoped "## " section, where the other tasks'
# markers are. Task-scoped are the "## " sections under "# Steps" or, in a plan with no "# Steps"
# heading, every "## " section (each up to the next "#" or "## " heading). When neither has one,
# the brief says "none stated in this section" and a note goes to stderr; with no "# Steps"
# heading a stderr note also names each marker paragraph left out (in another "## " section) and
# its line. The stop rule is read from the section, else from the whole plan. Without --section
# every marker is read from the whole plan.
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
# (for a bullet) next bullet. Fenced code blocks are skipped. Lines LO..HI only when given, and
# never the lines of SKIP ("a b;c d;": line ranges a..b, c..d) when given. With a fifth argument
# "lines", print only the line number where each paragraph starts.
paragraphs() {
    awk -v marker="$1" -v lo="${2:-1}" -v hi="${3:-999999999}" -v skip="${4:-}" \
        -v mode="${5:-}" '
        BEGIN { ns = split(skip, S, ";") }
        function skipped(n,    i, p) {
            for (i = 1; i <= ns; i++) {
                if (S[i] == "") continue
                split(S[i], p, " ")
                if (p[1] <= n && n <= p[2]) return 1
            }
            return 0
        }
        /^[[:space:]]*(```|~~~)/ { fence = !fence; inrule = 0; next }
        fence { next }
        NR < lo || NR > hi || skipped(NR) { inrule = 0; next }
        inrule && ($0 ~ /^[[:space:]]*$/ || $0 ~ /^#/ || (bullet && $0 ~ /^[-*] /)) { inrule = 0 }
        !inrule && (index($0, marker) == 1 ||
                    ($0 ~ /^[-*] / && index(substr($0, 3), marker) == 1)) {
            inrule = 1
            bullet = ($0 ~ /^[-*] /)
            if (bullet) $0 = substr($0, 3)
            if (mode == "lines") print NR; else print
            next
        }
        inrule && mode != "lines" { print }
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

# --- plan file, as a path relative to the repository root (absolute if outside it)
PLAN_REL=""
case "$PLAN" in
    /*)
        [ -f "$PLAN" ] || die "plan not found: $PLAN"
        PLAN_ABS="$(cd "$(dirname "$PLAN")" && pwd -P)/$(basename "$PLAN")"
        case "$PLAN_ABS" in
            "$ROOT"/*) PLAN_REL="${PLAN_ABS#"$ROOT"/}" ;;
            *)
                # e.g. Track E's plan in ~/.claude/plans on Newt's machine: the brief names it by
                # its absolute path, which only an agent on this machine can read.
                PLAN_REL="$PLAN_ABS"
                echo "note: plan outside the repository; the brief names it by absolute path" >&2
                ;;
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
    # Line range of the section: its heading up to the next heading of the same or a higher level.
    RANGE=$(awk -v s="$SECTION" -v exact="$EXACT" '
        /^[[:space:]]*(```|~~~)/ { fence = !fence; next }
        fence || !/^#+ / { next }
        {
            h = $0; sub(/^#+ +/, "", h); sub(/[[:space:]]+$/, "", h)
            lvl = index($0, " ") - 1
        }
        start && lvl <= slvl { print start, NR - 1; done = 1; exit }
        !start && ((exact == 1 && h == s) || (exact == 0 && index(h, s))) { start = NR; slvl = lvl }
        END { if (start && !done) print start, 999999999 }
    ' "$PLAN_REL")
    SEC_LO=${RANGE% *}
    SEC_HI=${RANGE#* }
    # Task-scoped "## " sections (each up to the next "#" or "## " heading): the other tasks'
    # markers live there, so the plan-wide default is the text outside all of them. They are the
    # "## " sections under "# Steps"; a plan with no "# Steps" heading (the new_todo.sh template)
    # has no such boundary, so there every "## " section is task-scoped.
    HAS_STEPS=$(awk '
        /^[[:space:]]*(```|~~~)/ { fence = !fence; next }
        !fence && /^# / { h = $0; sub(/^# +/, "", h); sub(/[[:space:]]+$/, "", h) }
        !fence && /^# / && h == "Steps" { print 1; found = 1; exit }
        END { if (!found) print 0 }
    ' "$PLAN_REL")
    SCOPED_SECTIONS=$(awk -v all="$((1 - HAS_STEPS))" '
        /^[[:space:]]*(```|~~~)/ { fence = !fence; next }
        fence || !/^#+ / { next }
        {
            h = $0; sub(/^#+ +/, "", h); sub(/[[:space:]]+$/, "", h)
            lvl = index($0, " ") - 1
        }
        lvl <= 2 && sub_start { printf "%d %d;", sub_start, NR - 1; sub_start = 0 }
        lvl == 1 { steps = (h == "Steps"); next }
        (all || steps) && lvl == 2 { sub_start = NR }
        END { if (sub_start) printf "%d %d;", sub_start, 999999999 }
    ' "$PLAN_REL")
fi
# heading_before <line>: the last "## " heading above the line, outside code fences
heading_before() {
    awk -v n="$1" '
        NR >= n { exit }
        /^[[:space:]]*(```|~~~)/ { fence = !fence; next }
        !fence && /^## / { h = $0 }
        END { print h }
    ' "$PLAN_REL"
}
# marker <text>: the paragraphs of the whole plan; with --section, those inside the section, else
# those outside every task-scoped "## " section (the plan-wide default), else nothing. In a plan
# with no "# Steps" heading, a stderr note names each paragraph left out and its line.
marker() {
    if [ -z "$SECTION" ]; then
        paragraphs "$1"
        return
    fi
    local found="" all="" kept="" ln=""
    found=$(paragraphs "$1" "$SEC_LO" "$SEC_HI")
    if [ -z "$found" ]; then
        found=$(paragraphs "$1" 1 999999999 "$SCOPED_SECTIONS")
        if [ "$HAS_STEPS" -eq 0 ]; then
            all=$(paragraphs "$1" 1 999999999 "" lines)
            kept=$(paragraphs "$1" 1 999999999 "$SCOPED_SECTIONS" lines)
            for ln in $all; do
                if ! printf '%s\n' "$kept" | grep -qx "$ln"; then
                    echo "note: skipped the $1 paragraph at line $ln of $PLAN_REL: it is in" \
                        "\"$(heading_before "$ln")\" and the plan has no \"# Steps\" heading," \
                        "so every \"## \" section is its own task's" >&2
                fi
            done
        fi
    fi
    printf '%s' "$found"
}
# marker_or_whole <text>: the paragraphs inside --section if it has any, else the whole plan's
# (the stop rule: a rule written in another section still applies)
marker_or_whole() {
    local found=""
    [ -z "$SECTION" ] || found=$(paragraphs "$1" "$SEC_LO" "$SEC_HI")
    [ -n "$found" ] || found=$(paragraphs "$1")
    printf '%s' "$found"
}
# none_note <marker>: the stderr note for a --section brief that states the marker nowhere
none_note() {
    local where='the "## " sections under "# Steps"'
    [ "$HAS_STEPS" -eq 1 ] || where='every "## " section (the plan has no "# Steps" heading)'
    echo "note: section \"$SECTION\" has no $1 paragraph and the plan has none outside" \
        "$where; the brief says \"none stated in this section\"" >&2
}
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
EXCEPTION=$(marker '**Process limit exception:**')
if [ -n "$EXCEPTION" ]; then
    COMPUTE="Compute (plan exception to the machine limit of $MAXPROCS heavy processes): $EXCEPTION"
elif [ -n "$SECTION" ]; then
    none_note '**Process limit exception:**'
    COMPUTE="Compute: process limit exception none stated in this section; at most $MAXPROCS heavy"
    COMPUTE="$COMPUTE processes at once on this machine."
else
    COMPUTE="Compute: at most $MAXPROCS heavy processes at once on this machine."
fi

# --- scope
LIBRARY=$(marker '**Library changes:**')
if [ -n "$LIBRARY" ]; then
    SCOPE="Scope: $LIBRARY"
elif [ -n "$SECTION" ]; then
    none_note '**Library changes:**'
    SCOPE='Scope: library changes none stated in this section, so no change under'
    SCOPE="$SCOPE \`icenine_py/icenine/\`."
else
    SCOPE='Scope: no change under `icenine_py/icenine/`.'
fi

# --- stop rule (every paragraph starting with the marker, up to the next blank line)
STOP=$(marker_or_whole '**Stop rule:**' | sed 's/^/> /')

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
        BRANCH_LINE="Branch: in your own worktree, run \`git checkout $CUR\` (parent"
        BRANCH_LINE="$BRANCH_LINE \`$CUR_PARENT\`; the main session frees it first) and add new"
        BRANCH_LINE="$BRANCH_LINE commits there."
    elif git rev-parse --verify -q "refs/heads/$CUR-$TASK" >/dev/null; then
        BRANCH_LINE="Branch: in your own worktree, run \`git checkout $CUR-$TASK\` (parent"
        BRANCH_LINE="$BRANCH_LINE \`$CUR\`; the main session frees it first) and add new"
        BRANCH_LINE="$BRANCH_LINE commits there."
    else
        NEW="$CUR-$TASK-fixes"
        BRANCH_LINE="Branch: the task is already merged, so in your own worktree run"
        BRANCH_LINE="$BRANCH_LINE \`git checkout -b $NEW $CUR\` and"
        BRANCH_LINE="$BRANCH_LINE \`git config branch.$NEW.parent $CUR\`."
    fi
else
    NEW="$CUR-$TASK"
    BRANCH_LINE="Branch: in your own worktree run \`git checkout -b $NEW $CUR\` and"
    BRANCH_LINE="$BRANCH_LINE \`git config branch.$NEW.parent $CUR\`."
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
echo "Refusals: if a tool or hook refuses an action (Edit or Write refused for a path, a hook" \
    "denial, a declined permission), stop and report it; never redo the action another way," \
    "such as writing files through Bash."
echo "Interruptions: commit each finished piece as you go (sleep or an app restart can stop you" \
    "at any point); when resumed, start from \`git status\` and \`git log\`."
echo "$COMPUTE"
echo "Commits: stage by explicit path; end each commit message with:"
echo "$TRAILER"
echo "Report: write the full report to \`$REPORT\`; the final message is at most about 250 words."
