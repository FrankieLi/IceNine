#!/bin/bash
# Usage: verify_handback.sh <parent> <report> [--branch B] [--fixes REVIEW_REPORT] [--no-tests]
#                           [--agent I|N]
# Mechanical gate on an implementer hand-back, run before /review-task (it is not a review).
# Each check prints `ok`, `FLAG` (the main session decides) or `FAIL` (exit 1):
#   commits of <parent>..<branch> (none = FAIL; hashes named in the report and on the branch must
#   agree, FLAG; a named hash the branch contains, such as the base or parent tip it started from,
#   is a note, not a FLAG); uncommitted changes in the branch's worktree outside CLAUDE.md and
#   .claude/**;
#   changed files by class (shared Track E/S files, other icenine_py/icenine/**, Src/**,
#   pyproject.toml and uv.lock, golden files: FLAG; DECISIONS.md, CLAUDE.md, .claude/** committed:
#   FAIL; ROADMAP.md lines outside the caller's own rows and mailbox: FLAG); the eight report
#   fields of .claude/agents/implementer.md (missing = FAIL; over 80 lines = FLAG); timing labels
#   (a Results line with a time unit and none of single-worker / contended / quiet, on the line or
#   in the Results heading line: FLAG, a heuristic); with --fixes, every number of the review's
#   "Required fixes" list must appear as a fix number in the report (FLAG per missing number);
#   the tests/...py[::name] paths (or bare test_*.py[::name] names, taken as icenine_py/tests/<name>;
#   FLAG if absent) named on the report's "Tests" field, run in the branch's worktree
#   (FAIL on a failure, FLAG if none named or no worktree holds the branch; --no-tests skips); the
#   line gives the passed / failed / skipped counts (from a JUnit XML file pytest writes through
#   PYTEST_ADDOPTS, else pytest's count line) and pytest's exit code (FLAG when exit 0 comes
#   with neither, since then no count shows that a test ran).
# The caller's own ROADMAP territory comes from `git config --local icenine.agent` or --agent:
# I = Track S, row H12, "From Ilium Agent"; N = Track E, "From Newt Agent". With neither set the
# script refuses to run (exit 2; it never defaults to one agent). The pytest command can be
# overridden with VERIFY_HANDBACK_PYTEST_CMD (default `uv run pytest`; the test paths are
# appended; it runs in <worktree>/icenine_py). <report> and --fixes may be relative paths.
set -euo pipefail
PARENT=""
REPORT=""
BRANCH=""
FIXES=""
NO_TESTS=0
AGENT=""
POS=0
while [ $# -gt 0 ]; do
    case "$1" in
        --branch) shift; BRANCH="${1:-}" ;;
        --fixes) shift; FIXES="${1:-}" ;;
        --no-tests) NO_TESTS=1 ;;
        --agent) shift; AGENT="${1:-}" ;;
        -*) echo "unknown option: $1" >&2; exit 2 ;;
        *)
            POS=$((POS + 1))
            if [ $POS -eq 1 ]; then PARENT="$1"; elif [ $POS -eq 2 ]; then REPORT="$1"
            else echo "unexpected argument: $1" >&2; exit 2; fi ;;
    esac
    shift
done
if [ -z "$PARENT" ] || [ -z "$REPORT" ]; then
    echo "usage: verify_handback.sh <parent> <report> [--branch B] [--fixes REVIEW] [--no-tests]" \
        "[--agent I|N]" >&2
    exit 2
fi
# make the file arguments absolute before the cd, so a relative path means what the caller meant
case "$REPORT" in /*) ;; *) REPORT="$PWD/$REPORT" ;; esac
if [ -n "$FIXES" ]; then case "$FIXES" in /*) ;; *) FIXES="$PWD/$FIXES" ;; esac; fi
ROOT=$(git rev-parse --show-toplevel)
cd "$ROOT"
[ -n "$BRANCH" ] || BRANCH=$(git rev-parse --abbrev-ref HEAD)
git rev-parse --verify -q "$PARENT^{commit}" >/dev/null || { echo "error: no such parent '$PARENT'" >&2; exit 2; }
git rev-parse --verify -q "$BRANCH^{commit}" >/dev/null || { echo "error: no such branch '$BRANCH'" >&2; exit 2; }
[ -f "$REPORT" ] || { echo "error: report '$REPORT' not found" >&2; exit 2; }
if [ -n "$FIXES" ] && [ ! -f "$FIXES" ]; then echo "error: review report '$FIXES' not found" >&2; exit 2; fi
[ -n "$AGENT" ] || AGENT=$(git config --local --get icenine.agent || true)
case "$AGENT" in
    I) TRACK_HDR="# Track S:"; OWN_ROW="| H12 |"; OWN_MB="## From Ilium Agent" ;;
    N) TRACK_HDR="# Track E:"; OWN_ROW=""; OWN_MB="## From Newt Agent" ;;
    "")
        echo "error: no agent letter; run 'git config --local icenine.agent I|N' or pass --agent" >&2
        exit 2 ;;
    *) echo "error: agent must be I or N, got '$AGENT'" >&2; exit 2 ;;
esac

NFAIL=0
NFLAG=0
ok() { echo "ok    $*"; }
flag() { echo "FLAG  $*"; NFLAG=$((NFLAG + 1)); }
fail() { echo "FAIL  $*"; NFAIL=$((NFAIL + 1)); }
wt_of() {  # print the worktree path holding refs/heads/$1, if any
    git worktree list --porcelain |
        awk -v b="refs/heads/$1" '/^worktree /{w=substr($0,10)} $1=="branch" && $2==b {print w}'
}
short_list() {  # first 4 lines of stdin joined by ", ", with a count of the rest
    awk 'NR<=4 { s = (s == "" ? $0 : s ", " $0) } END { if (NR > 4) s = s " (+" NR-4 " more)"; print s }'
}

# Field lines: `Field:`, `**Field:**`, `**Field** (note):`, `Field (note)`, `## Field (note)`,
# `4. Field:`, `Branch/commits:`; case-insensitive; `Fieldfoo:` and `Result:` are not.
FIELDS="branch / commits|files|tests|results|criteria|deviations from the plan|suspicious / unexpected|left undone"
field_block() {  # $1 field name (lowercase): its lines in the report, field line included
    awk -v want="$1" -v fields="$FIELDS" '
        BEGIN { n = split(fields, F, "|") }
        {
            l = tolower($0); head = (l ~ /^[ \t]*#/); start = ""
            sub(/^[#>*_ \t-]+/, "", l); sub(/^[0-9]+[.)][ \t]*/, "", l); sub(/^[#>*_ \t-]+/, "", l)
            gsub(/[ \t]*\/[ \t]*/, " / ", l)
            for (i = 1; i <= n; i++) {
                f = F[i]
                if (index(l, f) == 1) {
                    r = substr(l, length(f) + 1)
                    # the field name, then optionally a parenthetical, then a colon or the line end;
                    # a heading line may continue with any text after the name
                    if (r ~ /^[*_]*[ \t]*([(][^)]*[)])?[*_]*[ \t]*(:|$)/) { start = f; break }
                    if (head && r ~ /^[^a-z0-9]/) { start = f; break }
                }
            }
            if (start != "") cur = start
            if (cur == want) print
        }' "$REPORT"
}

echo "verify_handback: $BRANCH against $PARENT, report $(basename "$REPORT")"

# --- commits ---------------------------------------------------------------------------------
BR_COMMITS=$(git rev-list "$PARENT..$BRANCH")
if [ -z "$BR_COMMITS" ]; then
    fail "no commits on $BRANCH since $PARENT"
else
    NC=$(printf '%s\n' "$BR_COMMITS" | wc -l | tr -d ' ')
    ok "$NC commit(s) since $PARENT"
    git log --oneline "$PARENT..$BRANCH" | head -8 | sed 's/^/        /'
    TOKENS=$(field_block "branch / commits" | tr -c '0-9a-zA-Z' '\n' | grep -xE '[0-9a-f]{7,40}' | sort -u || true)
    if [ -z "$TOKENS" ]; then
        flag "report names no commit hash under \"Branch / commits\""
    else
        MISSING=""
        for H in $BR_COMMITS; do
            hit=0
            for T in $TOKENS; do case "$H" in "$T"*) hit=1; break ;; esac; done
            [ $hit -eq 1 ] || MISSING="$MISSING$(printf '%s' "$H" | cut -c1-9)"$'\n'
        done
        EXTRA=""
        ANCESTORS=""
        for T in $TOKENS; do
            on=0
            for H in $BR_COMMITS; do case "$H" in "$T"*) on=1; break ;; esac; done
            [ $on -eq 0 ] || continue
            # a commit the branch contains (the base or parent tip it started from) is on the
            # branch's history, so only hashes the branch does not contain are flagged
            if git rev-parse --verify -q "$T^{commit}" >/dev/null; then
                if git merge-base --is-ancestor "$T" "$BRANCH"; then
                    ANCESTORS="$ANCESTORS$T"$'\n'
                else
                    EXTRA="$EXTRA$T"$'\n'
                fi
            elif printf '%s' "$T" | grep -q '[a-f]' && printf '%s' "$T" | grep -q '[0-9]'; then
                EXTRA="$EXTRA$T"$'\n'
            fi
        done
        [ -z "$MISSING" ] || flag "branch commits not named in the report: $(printf '%s' "$MISSING" | short_list)"
        [ -z "$EXTRA" ] || flag "report names hashes not on the branch: $(printf '%s' "$EXTRA" | short_list)"
        if [ -z "$MISSING" ] && [ -z "$EXTRA" ]; then ok "commit hashes agree between report and branch"; fi
        [ -z "$ANCESTORS" ] ||
            echo "note  report names commits from before $PARENT..$BRANCH (base or parent):" \
                "$(printf '%s' "$ANCESTORS" | short_list)"
    fi
fi

# --- working tree ----------------------------------------------------------------------------
WT=$(wt_of "$BRANCH")
if [ -z "$WT" ]; then
    flag "no worktree holds $BRANCH: working tree not checked, tests cannot run"
else
    ST=$(git -C "$WT" status --porcelain | awk '{ p = substr($0, 4); gsub(/"/, "", p)
        if (p == "CLAUDE.md" || p ~ /^\.claude\//) print "C " p; else print "O " p }')
    DIRTY=$(printf '%s\n' "$ST" | sed -n 's/^O //p')
    CLAUDE_FILES=$(printf '%s\n' "$ST" | sed -n 's/^C //p')
    if [ -n "$DIRTY" ]; then
        flag "uncommitted changes in the worktree: $(printf '%s\n' "$DIRTY" | short_list)"
    else
        ok "worktree clean outside CLAUDE.md and .claude/**"
    fi
    [ -z "$CLAUDE_FILES" ] || echo "note  left uncommitted for the owner: $(printf '%s\n' "$CLAUDE_FILES" | short_list)"
fi

# --- changed files ---------------------------------------------------------------------------
CHANGED=$(git diff --name-only "$PARENT...$BRANCH")
class() {  # $1 label, $2 ERE on the path, $3 flag|fail, $4 note, $5 optional ERE to leave out
    local hits
    hits=$(printf '%s\n' "$CHANGED" | grep -E "$2" || true)
    if [ -n "${5:-}" ]; then hits=$(printf '%s\n' "$hits" | grep -vE "$5" || true); fi
    [ -n "$hits" ] || return 0
    local msg="$1: $(printf '%s\n' "$hits" | short_list)"
    [ -z "$4" ] || msg="$msg ($4)"
    if [ "$3" = fail ]; then fail "$msg"; else flag "$msg"; fi
}
SHARED='^icenine_py/icenine/(cost_functions\.py|orientation_search\.py|batched_screen\.py|reconstructor\.py|_rasterize\.c)$'
class "shared Track E/S file changed" "$SHARED" flag "a heads-up is due, protocol rule 3(a)"
class "library change" '^icenine_py/icenine/' flag "the plan must allow it" "$SHARED"
class "C++ change" '^Src/' flag ""
class "dependency files" '(^|/)(pyproject\.toml|uv\.lock)$' flag "new dependency?"
class "golden file" '^icenine_py/tests/data/[^/]*golden' flag "no golden regeneration inside a task"
class "DECISIONS.md committed" '^icenine_py/docs/DECISIONS\.md$' fail "the main session records decisions"
class "CLAUDE.md committed" '^CLAUDE\.md$' fail ""
class ".claude/** committed" '^\.claude/' fail "left for the owner"
if printf '%s\n' "$CHANGED" | grep -qx 'icenine_py/docs/ROADMAP.md'; then
    RM=icenine_py/docs/ROADMAP.md
    BASE=$(git merge-base "$PARENT" "$BRANCH")
    ranges() {  # stdin ROADMAP text; prints "start end;" for the caller's own territory
        awk -v track="$TRACK_HDR" -v mb="$OWN_MB" -v row="$OWN_ROW" '
            /^# / { if (ton) { printf "%d %d;", ts, NR - 1; ton = 0 } }
            /^#/ { if (mon) { printf "%d %d;", ms, NR - 1; mon = 0 } }
            track != "" && index($0, track) == 1 { ton = 1; ts = NR }
            mb != "" && $0 == mb { mon = 1; ms = NR }
            row != "" && index($0, row) == 1 { printf "%d %d;", NR, NR }
            END { if (ton) printf "%d %d;", ts, NR; if (mon) printf "%d %d;", ms, NR }'
    }
    OLDR=$(git show "$BASE:$RM" | ranges)
    NEWR=$(git show "$BRANCH:$RM" | ranges)
    BAD=$(git diff -U0 "$BASE" "$BRANCH" -- "$RM" | awk -v oldr="$OLDR" -v newr="$NEWR" '
        function inr(rs, lo, hi,   n, a, i, p) {
            n = split(rs, a, ";")
            for (i = 1; i <= n; i++) { if (a[i] == "") continue; split(a[i], p, " "); if (p[1] <= lo && hi <= p[2]) return 1 }
            return 0
        }
        /^@@ / {
            split($2, o, ","); split($3, w, ",")
            oa = substr(o[1], 2) + 0; ob = (o[2] == "" ? 1 : o[2] + 0)
            wa = substr(w[1], 2) + 0; wb = (w[2] == "" ? 1 : w[2] + 0)
            bad = 0
            if (ob > 0 && !inr(oldr, oa, oa + ob - 1)) bad = 1
            if (wb > 0 && !inr(newr, wa, wa + wb - 1)) bad = 1
            if (bad) print "lines " (wb > 0 ? wa "-" wa + wb - 1 : "removed at " oa)
        }')
    if [ -n "$BAD" ]; then
        flag "ROADMAP.md changed outside the caller's own rows and mailbox: $(printf '%s\n' "$BAD" | short_list)"
    else
        ok "ROADMAP.md changes are inside the caller's own rows and mailbox"
    fi
fi
NCH=$(printf '%s\n' "$CHANGED" | grep -c . || true)
echo "note  $NCH file(s) changed since $PARENT"

# --- report ----------------------------------------------------------------------------------
MISSING_F=""
for F in "branch / commits" "files" "tests" "results" "criteria" "deviations from the plan" "suspicious / unexpected" "left undone"; do
    [ -n "$(field_block "$F")" ] || MISSING_F="$MISSING_F$F"$'\n'
done
if [ -n "$MISSING_F" ]; then
    fail "report fields missing: $(printf '%s' "$MISSING_F" | awk '{ s = (s == "" ? $0 : s "; " $0) } END { print s }')"
else
    ok "all eight report fields present"
fi
RL=$(wc -l < "$REPORT" | tr -d ' ')
if [ "$RL" -gt 80 ]; then flag "report has $RL lines (about 60 expected)"; else ok "report length $RL lines"; fi

# timing labels: a Results line with a time unit needs a label on the line or in the heading line
RES=$(field_block results)
if [ -n "$RES" ]; then
    HEAD1=$(printf '%s\n' "$RES" | head -1)
    UNIT='[0-9][[:space:]]?(µs|us|ms|s)([^A-Za-z]|$)'
    LABEL='single-worker|contended|quiet'
    if printf '%s\n' "$HEAD1" | grep -Eiq "$LABEL"; then
        NOLAB=0
    else
        NOLAB=$(printf '%s\n' "$RES" | grep -E "$UNIT" | grep -Eivc "$LABEL" || true)
    fi
    if [ "${NOLAB:-0}" -gt 0 ]; then
        flag "$NOLAB Results line(s) with a time unit and no single-worker / contended / quiet label"
    else
        ok "no unlabelled timing line in Results (heuristic)"
    fi
fi

# --- required fixes --------------------------------------------------------------------------
if [ -n "$FIXES" ]; then
    REQ=$(awk '
        on && /^#/ { exit }
        tolower($0) ~ /^[#*[:space:]]*required fixes/ { on = 1; next }
        on { print }' "$FIXES")
    NUMS=$(printf '%s\n' "$REQ" | sed -nE 's/^[[:space:]]*([0-9]+)[.)].*/\1/p')
    if [ -z "$NUMS" ]; then
        flag "no numbered items under \"Required fixes\" in $(basename "$FIXES")"
    else
        SEEN=$(awk '
            /^[[:space:]]*[-*]?[[:space:]]*[0-9]+[.):]/ { t = $0; sub(/^[^0-9]*/, "", t); sub(/[^0-9].*/, "", t); print t }
            tolower($0) ~ /fix/ {
                s = $0
                gsub(/[0-9a-f][0-9a-f][0-9a-f][0-9a-f][0-9a-f][0-9a-f][0-9a-f]+/, "", s)
                while (match(s, /[0-9]+([ ]?(-|to)[ ]?[0-9]+)?/)) {
                    tok = substr(s, RSTART, RLENGTH); s = substr(s, RSTART + RLENGTH)
                    if (tok ~ /[^0-9]/) { n = split(tok, p, /[^0-9]+/); for (k = p[1]; k <= p[n]; k++) print k }
                    else print tok
                }
            }' "$REPORT" | sort -u)
        MISS=""
        for N in $NUMS; do
            printf '%s\n' "$SEEN" | grep -qx "$N" || MISS="$MISS$N"$'\n'
        done
        if [ -n "$MISS" ]; then
            flag "required fixes not named in the report: $(printf '%s' "$MISS" | awk '{ s = (s == "" ? $0 : s ", " $0) } END { print s }')"
        else
            ok "all $(printf '%s\n' "$NUMS" | wc -l | tr -d ' ') required fixes named in the report"
        fi
    fi
fi

# --- tests -----------------------------------------------------------------------------------
if [ $NO_TESTS -eq 1 ]; then
    echo "note  tests skipped (--no-tests)"
elif [ -n "$WT" ]; then
    # `tests/...py[::name]` as written, or a bare `test_*.py[::name]`, which is taken to be
    # icenine_py/tests/<name>; a bare name that is not there is flagged and not run
    TWORDS=$(field_block tests | awk '{
        gsub(/[`,;()]/, " ")
        for (i = 1; i <= NF; i++) {
            w = $i; sub(/[.:]+$/, "", w)
            if (match(w, /tests\/[A-Za-z0-9_\/.-]+\.py(::.*)?$/)) print "P " substr(w, RSTART)
            else if (w ~ /^[A-Za-z0-9_-]*test_[A-Za-z0-9_.-]*\.py(::.*)?$/) print "B tests/" w
        }}' | sort -u || true)
    TPATHS=""
    NOTFOUND=""
    while read -r kind tp; do
        [ -n "${tp:-}" ] || continue
        if [ "$kind" = B ] && [ ! -f "$WT/icenine_py/${tp%%::*}" ]; then
            NOTFOUND="$NOTFOUND${tp#tests/}"$'\n'
        else
            TPATHS="$TPATHS$tp"$'\n'
        fi
    done <<EOF2
$TWORDS
EOF2
    TPATHS=$(printf '%s' "$TPATHS" | sort -u)
    [ -z "$NOTFOUND" ] ||
        flag "test file(s) named on the Tests field not found in icenine_py/tests: $(printf '%s' "$NOTFOUND" | short_list)"
    if [ -z "$TPATHS" ] && [ -n "$NOTFOUND" ]; then
        echo "note  no existing test path left: tests not run"
    elif [ -z "$TPATHS" ]; then
        flag "report names no tests/...py or test_*.py path on its Tests field: tests not run"
    else
        ARGS=()
        while IFS= read -r line; do ARGS+=("$line"); done <<EOF
$TPATHS
EOF
        # No -q here: icenine_py's addopts already has -q, and a second -q (verbosity -2) drops
        # pytest's final count line, so the last output line was a warning-summary line. The
        # counts come from a JUnit XML file requested through PYTEST_ADDOPTS (so an overriding
        # command gets the same arguments), else from pytest's count line in the output.
        CMD="${VERIFY_HANDBACK_PYTEST_CMD:-uv run pytest}"
        TDIR=$(mktemp -d)
        TLOG="$TDIR/out.log"
        JUNIT="$TDIR/junit.xml"
        set +e
        (cd "$WT/icenine_py" &&
            PYTEST_ADDOPTS="${PYTEST_ADDOPTS:-} --junitxml='$JUNIT'" \
                bash -c "$CMD \"\$@\"" verify "${ARGS[@]}") > "$TLOG" 2>&1
        RC=$?
        set -e
        OUT=$(tail -20 "$TLOG")
        COUNTS=""
        if [ -s "$JUNIT" ]; then
            # the first <testsuite ...> tag: tests includes skipped, failures and errors
            COUNTS=$(tr '\n' ' ' < "$JUNIT" | awk '{
                if (!match($0, /<testsuite[ \t][^>]*>/)) exit
                tag = substr($0, RSTART, RLENGTH)
                n = split("tests failures errors skipped", K, " ")
                for (i = 1; i <= n; i++) {
                    v[K[i]] = 0
                    if (match(tag, "[ \t]" K[i] "=\"[0-9]+\"")) {
                        s = substr(tag, RSTART, RLENGTH); gsub(/[^0-9]/, "", s); v[K[i]] = s + 0
                    }
                }
                p = v["tests"] - v["failures"] - v["errors"] - v["skipped"]
                out = p " passed, " v["failures"] " failed, " v["skipped"] " skipped"
                if (v["errors"] > 0) out = out ", " v["errors"] " errors"
                print out
            }')
        fi
        if [ -z "$COUNTS" ]; then
            # pytest's count line, e.g. "==== 1 failed, 2 passed, 3 warnings in 0.12s ====";
            # never a warning or "-- Docs:" line
            COUNTS=$(awk '{
                l = $0; gsub(/^[= \t]+|[= \t]+$/, "", l); sub(/ in [0-9.]+s( \([^)]*\))?$/, "", l)
                if (l !~ /^[0-9]+ [a-z]+(, [0-9]+ [a-z]+)*$/) next
                p = 0; f = 0; s = 0; e = 0; n = split(l, it, ", ")
                for (i = 1; i <= n; i++) {
                    split(it[i], w, " ")
                    if (w[2] == "passed") p = w[1]; else if (w[2] == "failed") f = w[1]
                    else if (w[2] == "skipped") s = w[1]
                    else if (w[2] == "error" || w[2] == "errors") e = w[1]
                }
                c = p " passed, " f " failed, " s " skipped"
                if (e > 0) c = c ", " e " errors"
            } END { print c }' "$TLOG")
        fi
        rm -rf "$TDIR"
        [ -n "$COUNTS" ] || COUNTS="no count line in the output"
        BROKEN=$(printf '%s\n' "$COUNTS" | awk -F', ' '{ for (i = 1; i <= NF; i++) {
            split($i, w, " "); if ((w[2] == "failed" || w[2] == "errors") && w[1] > 0) b = 1 } }
            END { print b + 0 }')
        if [ "$RC" -ne 0 ] || [ "$BROKEN" -ne 0 ]; then
            fail "named tests failed (${#ARGS[@]} path(s)): $COUNTS; pytest exit $RC; last lines:"
            printf '%s\n' "$OUT" | tail -6 | sed 's/^/        | /'
        elif [ "$COUNTS" = "no count line in the output" ]; then
            # exit 0 but neither a JUnit file nor a count line (an overriding command that drops
            # PYTEST_ADDOPTS and the count line): the gate cannot show that any test ran
            flag "named tests (${#ARGS[@]} path(s)): no counts (no JUnit file, no count line in" \
                "the output); pytest exit 0"
        else
            ok "named tests passed (${#ARGS[@]} path(s)): $COUNTS; pytest exit 0"
        fi
    fi
fi

echo "summary: $NFAIL FAIL, $NFLAG FLAG"
[ $NFAIL -eq 0 ] || exit 1
exit 0
