#!/bin/bash
# PreToolUse hook (Bash): deny bare python/pytest/pip, blanket git staging/committing, git commit
# --no-verify and ALLOW_* override assignments. Reads the hook JSON on stdin; prints a
# permissionDecision "deny" object to block, else nothing. (Constrains agents only; the owner runs
# git outside Claude Code.)
INPUT=$(cat)
CMD=$(printf '%s' "$INPUT" | jq -r '.tool_input.command // empty')
[ -z "$CMD" ] && exit 0

block() {
    jq -cn --arg r "$1" '{hookSpecificOutput: {hookEventName: "PreToolUse",
        permissionDecision: "deny", permissionDecisionReason: $r}}'
    exit 0
}

# 1. join backslash continuations; 2. drop heredoc bodies; 3. drop quoted strings (commit
# messages, echo text) so their content is never inspected.
STRIPPED=$(printf '%s' "$CMD" | perl -0pe '
    s/\\\n//g;
    s/(?<!<)<<(?!<)-?\s*["\x27]?(\w+)["\x27]?([^\n]*)\n.*?^\s*\1$/<<$2/gms;
    s/"(?:[^"\\]|\\.)*"|\x27[^\x27]*\x27/""/gs;
')
# Split compound commands on && || ; | & newlines, subshell parens, braces and backticks.
SEGMENTS=$(printf '%s' "$STRIPPED" |
    sed -e 's/&&/\n/g' -e 's/||/\n/g' -e 's/\$(/\n/g' -e 's/[;|&(){}`]/\n/g')

BARE_RE='^(python[0-9.]*|pip[0-9.]*|pytest)$'
ASSIGN_RE='^[A-Za-z_][A-Za-z_0-9]*='
UV_MSG='Use uv: "uv run pytest" instead of bare "pytest", "uv run python" instead of "python3", '
UV_MSG+='"uv pip" instead of "pip". This project uses uv for Python environment management.'

while IFS= read -r seg; do
    set -f
    # shellcheck disable=SC2206
    tok=($seg)
    set +f
    [ ${#tok[@]} -eq 0 ] && continue
    for t in "${tok[@]}"; do
        if [[ "$t" =~ ^(export)?ALLOW_(CLAUDE_CONFIG|LARGE|ABS_PATHS)= ]]; then
            block 'Do not set ALLOW_* pre-commit overrides from an agent; ask the owner.'
        fi
    done
    # skip leading assignments and wrappers (nohup, time, sudo, nice [-n N], env [-u X] ...)
    i=0
    while [ $i -lt ${#tok[@]} ]; do
        t="${tok[$i]}"
        if [[ "$t" =~ $ASSIGN_RE ]]; then
            i=$((i + 1))
        elif [[ "$t" =~ ^(nohup|time|exec|command)$ ]]; then
            i=$((i + 1))
        elif [ "$t" = "sudo" ]; then
            i=$((i + 1))
            while [[ "${tok[$i]:-}" == -* ]]; do i=$((i + 1)); done
        elif [ "$t" = "nice" ]; then
            i=$((i + 1))
            nxt="${tok[$i]:-}"
            if [ "$nxt" = "-n" ]; then
                i=$((i + 2))
            elif [[ "$nxt" =~ ^-[0-9]+$ ]]; then
                i=$((i + 1))
            fi
        elif [ "$t" = "env" ]; then
            i=$((i + 1))
            while [ $i -lt ${#tok[@]} ]; do
                if [ "${tok[$i]}" = "-u" ]; then i=$((i + 2))
                elif [ "${tok[$i]}" = "-i" ] || [[ "${tok[$i]}" =~ $ASSIGN_RE ]]; then i=$((i + 1))
                else break; fi
            done
        else
            break
        fi
    done
    [ $i -lt ${#tok[@]} ] || continue
    first="${tok[$i]}"
    # (a) bare python / pytest / pip as the command
    if [[ "$first" =~ $BARE_RE ]]; then
        block "$UV_MSG"
    fi
    # (b) git add / git commit
    [ "$first" = "git" ] || continue
    j=$((i + 1))
    while [ $j -lt ${#tok[@]} ]; do # skip global options such as -C <dir>, -c k=v
        case "${tok[$j]}" in
            -C | -c) j=$((j + 2)) ;;
            -*) j=$((j + 1)) ;;
            *) break ;;
        esac
    done
    sub="${tok[$j]:-}"
    args=("${tok[@]:$((j + 1))}")
    if [ "$sub" = "add" ]; then
        paths=0
        update=0
        for a in "${args[@]}"; do
            case "$a" in
                --all) block 'Stage files explicitly by path: "git add --all" is not allowed.' ;;
                --update) update=1 ;;
                --*) ;;
                -*A*) block 'Stage files explicitly by path: "git add -A" is not allowed.' ;;
                -*u*) update=1 ;;
                -*) ;;
                . | ./ | :/ | :/. | '*' | :/\*)
                    block 'Stage files explicitly by path: "git add ." is not allowed.' ;;
                *) paths=$((paths + 1)) ;;
            esac
        done
        if [ $update -eq 1 ] && [ $paths -eq 0 ]; then
            block 'Stage files explicitly by path: "git add -u" without paths is not allowed.'
        fi
    elif [ "$sub" = "commit" ]; then
        for a in "${args[@]}"; do
            if [ "$a" = "--no-verify" ] || [[ "$a" =~ ^-[^-mFCcS]*n ]]; then
                block 'Do not skip the git hooks (--no-verify / -n); fix the reported problem.'
            fi
            if [ "$a" = "--all" ] || [[ "$a" =~ ^-[^-mFCcS]*a ]]; then
                block 'Stage files explicitly and commit without "-a"/"--all" (also "-am").'
            fi
        done
    fi
done <<<"$SEGMENTS"
exit 0
