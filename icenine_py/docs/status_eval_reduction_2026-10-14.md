# Where the evaluation-reduction work stands (2026-10-14)

Written for the owner by Newt Agent (Track E). A snapshot after Task C was merged into
`feature/eval-reduction` (e609c74). Details live in the linked documents; this page says where we are.

## The goal

Cut the number of cost-function evaluations the orientation search makes, without losing accuracy,
and, where the physics allows it, make the cost function smoother. The search runs in levels: at each
level it scores a grid of candidate orientations (the "screen"), runs a short Monte Carlo polish on each
survivor (the "quick MC", 12 evaluations per candidate), keeps the best quarter, and refines the grid
around them for the next level.

We test on two fixed voxel sets: six (and in Task C, thirty) interior voxels of the full 500-grain
sample, and "E0", 200 single voxels each with its own simulated image. Images are clean (no noise)
until the noise model is redone.

## What is done

| Task | Question | Answer in one line | Code change |
|---|---|---|---|
| A | Why do C++ and Python make different numbers of evaluations? | With the C++ overlap rule, Python makes 0.964 times C++'s evaluations; C++'s extra is a second, redundant re-cost of every returned candidate. The default (lenient) rule makes about 1.46 times C++'s, which is a method difference, not a bug. | `BatchedScreen` became the default (owner); `CppEdgeRules` and `CppOmegaInterp` added, off |
| A-exact | Can Python light a detector exactly as C++ does? | Yes: with `CppEdgeRules` on, the per-peak counts of the complete cost agree with C++ on 2,230 of 2,230 peaks. | inside `CppEdgeRules` (still off) |
| B | Why does the cost function have odd edges? | Near the truth (under 0.05 degree), 73.7% of the cost changes come from a spot crossing an omega-frame boundary (physics); farther out, pixel boundaries dominate (an artefact of pixels). | none (study only) |
| C | Do the candidates overlap, and can removing duplicates save evaluations? | They overlap only modestly, and every filter that saves evaluations also loses some voxels. Nothing goes forward yet. | `CandidateDedup` added, off |

### Task C in a little more detail

- **How much overlap.** On the six full-sample voxels, about 1 in 5 of the candidates handed to the
  quick MC has a better candidate within half a grid step. After the quick MC almost none end at the same
  place (0.06% within 0.1 degree), so the duplicates are near neighbours, not copies.
- **What the filters do.** Keeping only the best candidate of each C++ group cuts evaluations to about 76%
  on the full-sample voxels (0 of 6 wrong; 2 of 30 wrong against the baseline's 1). Removing candidates
  within half a grid step of a better one cuts them to about 91% (0 of 30 wrong).
- **Why that is not enough.** On the 200 E0 voxels the same filters save under 1% of evaluations, because
  92% of an E0 search is the screen and the re-cost, which these filters do not touch. They also change
  which voxels come out wrong: the best-per-group filter breaks 24 voxels and fixes 10 (55 wrong against
  41).
- **Choice for Task G.** The plan's rule takes a filter forward only if it breaks no more voxels than the
  baseline. None does, so Task G starts from the baseline (key off).

Full report: `.claude/reports/20261014-1507-eval-reduction-C.md` (on Newt Agent's machine). Results:
MIGRATION_HISTORY "Task C" and `benchmarks/eval_reduction/C/`.

## What the work has taught us so far

- The evaluations are spent in two places, and they need different fixes. On full-sample voxels the
  quick MC is about 55% of the work, so candidate filters can save up to about a quarter. On E0 voxels
  the screen and re-cost dominate, so only cheaper or fewer screen evaluations will help (Track S, S4).
- The truth is mostly lost at the "keep the best quarter" cut, not later. Any filter that thins the
  candidates before that cut risks removing the right one, which is what Task C's E0 numbers show.
- E0 counts move by themselves between code versions: today's baseline differs from the stored earlier
  run in 30 of the 200 voxels (41 wrong now, 39 before), cause not tested. Small broken/fixed counts
  on E0 should be read against that.

## A process note

Decision D-20261014-N1 says the code review runs on the task branch before the merge. Task C was
merged before its review, at the owner's instruction today (D-20261014-N3). The review now runs on the
merged commits (e5a2cbb..e609c74).

## What is next

1. Code review of Task C, with fixes from a fresh implementer.
2. Task D: which spots on the two detectors belong to one peak, and which peaks belong to a voxel.
3. Task E: smoother cost modes and ranking keys built from the geometry.
4. Task F: searching outward from a start orientation to 2 to 5 degrees.
5. Task G: putting it together, with a 2,000-voxel BFS pilot.

## Waiting on the owner

- Whether Task G's BFS pilot should also try the two Task C filters that look harmless on full-sample
  voxels (`clique=cpp:1`, `pre=0.5`), even though they lose voxels on E0.
- Whether to delete the untracked `benchmarks/eval_reduction/A/memory/cache/` folder (it holds logs).
- No defaults to flip and no golden tests to regenerate from Track E right now.

## The other track

Ilium Agent's Track S (cheaper evaluations) has merged S3, a batched quick MC behind the opt-in key
`BatchedQuickMC` (off). Its mailbox and rows are in `docs/ROADMAP.md`; owner decisions for both
tracks are in `docs/DECISIONS.md`.
