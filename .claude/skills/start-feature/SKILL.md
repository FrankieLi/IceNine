---
name: start-feature
description: Start a new feature branch off develop (gitflow). Use when beginning new work.
user-invocable: true
allowed-tools: Bash(git *)
argument-hint: "<feature-name>"
---

Create a new gitflow feature branch for: $ARGUMENTS

Follow these steps exactly:

1. Fetch latest from origin: `git fetch origin`
2. Checkout develop and pull: `git checkout develop && git pull origin develop`
3. Create the feature branch: `git checkout -b feature/$ARGUMENTS`
4. Push and set upstream: `git push -u origin feature/$ARGUMENTS`
5. Confirm the branch was created and show `git log --oneline -3`

Branch naming rules:
- Always prefix with `feature/`
- Use kebab-case (e.g., `feature/add-cost-function`)
- If `$ARGUMENTS` already starts with `feature/`, don't double-prefix
