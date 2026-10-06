#!/bin/bash
# Usage: new_todo.sh <slug> "<title>" ["<subtitle>"]  -- scaffold icenine_py/docs/todo_<slug>.md
set -e
SLUG="$1"
TITLE="$2"
SUB="${3:-}"
[ -n "$SLUG" ] && [ -n "$TITLE" ] || { echo 'usage: new_todo.sh <slug> "<title>" ["<subtitle>"]' >&2; exit 2; }
ROOT=$(git rev-parse --show-toplevel)
OUT="$ROOT/icenine_py/docs/todo_$SLUG.md"
mkdir -p "$(dirname "$OUT")"
[ ! -e "$OUT" ] || { echo "error: $OUT already exists" >&2; exit 1; }
cat >"$OUT" <<DOC
---
title: "TODO: $TITLE"
subtitle: "$SUB"
date: "$(date +%Y-%m-%d)"
geometry: margin=1in
fontsize: 11pt
---

# Status

**Idea, not started.** (Who raised it, when, and what it waits for.)

# What

(One paragraph: what would be built or measured.)

# Why it might help

(The problem it addresses and the evidence for it.)

# Framing

(How to state the problem; prefer symmetry-agnostic framing.)

# Experiments (sketch)

(The smallest experiments that would test the idea, with success criteria.)

# Risks

(What could make it not work, or mislead.)

# Relation to other work

(Links to related TODOs, reports and MIGRATION_HISTORY sections.)
DOC
echo "$OUT"
