#!/bin/bash
# Activate the repository git hooks in .githooks (applies to every worktree of this clone).
set -e
cd "$(git rev-parse --show-toplevel)"
chmod +x .githooks/*
git config core.hooksPath .githooks
echo "Hooks active: core.hooksPath=.githooks"
echo "To undo: git config --unset core.hooksPath"
echo "To bypass once: git commit --no-verify (or set ALLOW_CLAUDE_CONFIG/ALLOW_LARGE/ALLOW_ABS_PATHS=1)"
