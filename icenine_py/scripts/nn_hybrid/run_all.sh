#!/bin/sh
# Full Task-1 run, one pipeline at a time (10 workers each) so timings stay comparable.
cd "$(dirname "$0")/../.." || exit 1
for P in H3 H0 H1 H3c H3m HG; do
  uv run python scripts/nn_hybrid/run.py run --pipeline $P --workers 10 || exit 1
done
uv run python scripts/nn_hybrid/run.py run --pipeline H3 --model realistic_s1 --variants all --workers 10
