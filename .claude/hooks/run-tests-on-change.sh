#!/bin/bash
# PostToolUse hook: auto-run pytest when Python files in icenine_py are edited
INPUT=$(cat)
FILE_PATH=$(echo "$INPUT" | jq -r '.tool_input.file_path // empty')

# Only trigger for Python files in the icenine_py directory
if [[ "$FILE_PATH" == *"icenine_py/"* && "$FILE_PATH" == *.py ]]; then
    # Determine which test file to run based on the changed file
    BASENAME=$(basename "$FILE_PATH" .py)
    PROJECT_DIR="$CLAUDE_PROJECT_DIR/icenine_py"

    # If editing a test file, run that specific test
    if [[ "$BASENAME" == test_* ]]; then
        cd "$PROJECT_DIR" && pytest "tests/$BASENAME.py" -q --tb=short 2>&1
        exit $?
    fi

    # If editing a source file, look for a matching test
    TEST_FILE="$PROJECT_DIR/tests/test_${BASENAME}.py"
    if [[ -f "$TEST_FILE" ]]; then
        cd "$PROJECT_DIR" && pytest "tests/test_${BASENAME}.py" -q --tb=short 2>&1
        exit $?
    fi
fi

exit 0
