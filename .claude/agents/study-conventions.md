# Study conventions (standing brief)

Every implementer and reviewer reads this first. Briefs from the main session may say only
"follow study-conventions"; nothing here needs repeating in a brief.

## Environment and git
- All Python via `uv run`; never bare `python3`, `pip`, `pytest`; never `uv pip`.
- `uv sync --extra dev` for library and tests; `uv sync --extra dev --extra riemannian --extra
  benchmarks` for study scripts.
- Stage files by explicit path (never `git add -A`, `-a`, `.`). Never touch `CLAUDE.md`.
- Never set `ALLOW_*` overrides and never use `--no-verify`; if a check blocks you, fix it or report.
  `.claude/**` files are written but left uncommitted for the owner.
- Commits end with the Co-Authored-By line given in the system reminder or plan.
- Do not merge to develop, push, or open a PR unless the main session asks.
- Full suite: `uv run pytest tests/ --junitxml=<scratch>/junit.xml`, with no extra `-q`.

## Compute and timing
- At most 10 workers. Run long jobs with `nohup`, a log in `<dir>/cache/logs/<name>.log` and a
  `<name>.done` marker; poll every 10-15 minutes.
- Every timing run is preflight-gated (`preflight.require_quiet()` or
  `scripts/dev/timing_preflight.py --require --json <dir>/preflight.json`); save the JSON beside
  the results and label timings "single-worker" or "contended".
- Pilot a few items, timed, and check the output before the full budget.
- Spatial indexing: `scipy.spatial.cKDTree`, never `KDTree`.

## Wording
- Say "realistic" data, never "corrupted" (old `corrupt*` identifiers are aliases only). <!-- noqa: realistic -->
- Frame findings symmetry-agnostically (generic near-degenerate solutions), not around CSL or
  cubic-specific cases, unless the plan says otherwise.
- Causal wording only where tested; otherwise "consistent with", naming the untested part.
- No "bound", "monotone" or "no dependence" unless strictly true; say "no detectable ...".

## Numbers
- Every prose number comes from `summary.json` or `tables.md` (doc tables via
  `doc_tables.write_tables` and `scripts/dev/sync_doc_tables.py`), then run
  `scripts/dev/audit_numbers.py --doc D --section "H" --sources ...`. A match is not verification:
  re-read key numbers against the source.
- Format small angles with enough decimals to show differences.
- Give medians over a labelled subset (say which subset and n).
- No tie-dependent statistics (ties in win rates, ranks) without stating the tie rule.
- p-values: use voxel-clustered checks, not pooled per-trial tests.
- Small counts get a Wilson interval or a paired test (`icenine_py/scripts/common/stats.py`).
- Quote the evaluation cost with every improvement claim.
- Verify C++ parity claims against the C++ source, not by reading the Python alone.

## Reports
- Write the FULL report to `.claude/reports/<YYYYMMDD-HHMM>-<task>.md` (gitignored).
- The final message to the main session is at most about 250 words: report path, status or
  verdict, at most 6 lines of headline numbers, decisions the owner must make, and (reviewers) the
  required fixes as a numbered list of one-liners.
