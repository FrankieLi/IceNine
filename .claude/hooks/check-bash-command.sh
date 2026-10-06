#!/bin/bash
# PreToolUse hook (Bash): block bare python/pytest/pip and blanket git staging/committing.
# Reads the hook JSON on stdin; prints {"decision":"block","reason":...} to block, else nothing.
INPUT=$(cat)
CMD=$(printf '%s' "$INPUT" | jq -r '.tool_input.command // empty')
[ -z "$CMD" ] && exit 0

block() {
    jq -cn --arg r "$1" '{decision: "block", reason: $r}'
    exit 0
}

# Drop quoted strings (commit messages, echo text) so their content is never inspected, then
# split compound commands on && || ; | & and newlines.
STRIPPED=$(printf '%s' "$CMD" | perl -0pe "s/\"(?:[^\"\\\\]|\\\\.)*\"|'[^']*'/\"\"/gs")
SEGMENTS=$(printf '%s' "$STRIPPED" | sed -e 's/&&/\n/g' -e 's/||/\n/g' -e 's/[;|&]/\n/g')

while IFS= read -r seg; do
    # shellcheck disable=SC2206
    tok=($seg)
    [ ${#tok[@]} -eq 0 ] && continue
    # (a) bare python / pytest / pip at the start of a command segment
    if [[ "${tok[0]}" =~ ^(python3?|pytest|pip)$ ]]; then
        block 'Use uv: "uv run pytest" instead of bare "pytest", "uv run python" instead of "python3", "uv pip" instead of "pip". This project uses uv for Python environment management.'
    fi
    # (b) git add / git commit blanket flags
    [ "${tok[0]}" = "git" ] || continue
    i=1
    while [ $i -lt ${#tok[@]} ]; do # skip global options such as -C <dir>, -c k=v
        case "${tok[$i]}" in
            -C | -c) i=$((i + 2)) ;;
            -*) i=$((i + 1)) ;;
            *) break ;;
        esac
    done
    sub="${tok[$i]:-}"
    args=("${tok[@]:$((i + 1))}")
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
                .) block 'Stage files explicitly by path: "git add ." is not allowed.' ;;
                *) paths=$((paths + 1)) ;;
            esac
        done
        if [ $update -eq 1 ] && [ $paths -eq 0 ]; then
            block 'Stage files explicitly by path: "git add -u" without paths is not allowed.'
        fi
    elif [ "$sub" = "commit" ]; then
        for a in "${args[@]}"; do
            if [ "$a" = "--all" ] || [[ "$a" =~ ^-[^-mFCcS]*a ]]; then
                block 'Stage files explicitly and commit without "-a"/"--all" (also "-am").'
            fi
        done
    fi
done <<<"$SEGMENTS"
exit 0
