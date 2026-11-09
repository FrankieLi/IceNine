---
title: "TODO: Later stage: autoresearch"
subtitle: "Ideas the owner set aside for an automated research stage"
date: "2026-11-01"
geometry: margin=1in
fontsize: 11pt
---

# Status

**Ideas, not started.** This file collects ideas the owner set aside for a later automated
research stage. Each section is one idea; none is planned. They wait for that stage and the
owner's pick.

# A1. Order of the layer-wide screen and the BFS

Raised by the owner on 2026-11-01 on the review page (S4 thread), word for word: "Add to a todo
for a lager stage for autoresearch.... What order of operation, so that we can preserve the
advantage of the scan while still take advanget of the BFS? Perhaps scan across in batch, and for
anyone that takes, perform BFS. Then alternate back to scan? Do experiments" (D-20261101-I2).

## What

Find an order of operations that keeps the advantage of the layer-wide level-0 screen (one
table-driven pass over many voxels at once) and the advantage of the BFS (only seeds pay for a
full search; neighbours inherit an orientation cheaply). The owner's sketch: screen the layer in
batches; every voxel whose candidates "take" (a candidate passes the acceptance test after its
re-cost and refinement) becomes a seed, and the BFS grows from it; when the BFS stalls, go back to
screening the voxels still unassigned; alternate until the layer is done.

## Why it might help

* S4 (MIGRATION_HISTORY "S4", table `s4_projection`, Phase D, one process, projected): today's BFS
  takes 7.65-10.8 h, of which the search after the seeds' screen and re-cost is 82-84%. A full
  search for every voxel after one global screen takes 447 h (41-58x slower). The hybrid with a
  stored layer-wide local-cost table, BFS as today, saves about 10%.
* The screen of the whole layer takes 6.58 s on 10 threads (measured), against 3.94 h projected for
  `BatchedScreen` on one worker, so screening many voxels is cheap; what costs is the per-voxel
  search after it.
* Track E Task G: every wrong BFS voxel on the pilot is a boundary voxel carrying a touching
  grain's orientation. A pilot of the Phase D BFS left 28 of 60 grains with no accepted seed
  (MIGRATION_HISTORY "Phase D pilot results"). A screen candidate list for every voxel could show
  both: a boundary voxel whose inherited orientation is not among its own candidates, and a grain
  that no seed reached.

## Framing

A scheduling problem with three operations: a cheap batched screen (per voxel, table-driven), a
costly per-voxel re-cost and search (seeds and failed inheritances), and a cheap neighbour fit that
succeeds inside grains and can carry the wrong grain across a boundary. Choose the order and the
batch size. Measures: wall time on one process (quiet machine), wrong voxels (more than 1 degree
from the truth), and grains with no seed.

## Experiments (sketch)

On the 2,000-voxel BFS pilot region, clean images, one process:

1. Baseline: today's BFS (Track E's key-off rerun, D-20261029-N2).
2. One scan: screen the whole region in one pass; seed the BFS from the voxels with the best
   screened candidates in turn; reuse the screen results for every reseed.
3. Alternating (the owner's sketch): screen a batch of k unassigned voxels, accept those whose
   candidates take, grow the BFS from them until it stalls, screen the rest; sweep k.
4. As 3, plus a check of each boundary voxel's inherited orientation against its own screen
   candidates.

Success: no more wrong voxels than the baseline, judged against the seed-order spread (14 fixed and
14 broken between two code versions, findings doc "Limits"); no more grains with no seed; and a
lower wall time.

## Risks

* Under the Python default qualified-peak rule the level-0 screen passes 78.5% of the layer's
  orientation-voxel pairs (26.5% under `CppOverlapRule`), so it barely filters. "Takes" has to be
  decided by the re-cost and refinement, not by the screen.
* The layer-wide local-cost table is about 0.62 GB and takes 20.1 s to build (S4); S4's library
  path is parked.
* The acceptance test can accept a neighbouring grain at a boundary (Task F; D-20261029-N2).
* One seed order moves results by about 14 voxels each way, so small differences need several
  orders.

## Relation to other work

S4 (MIGRATION_HISTORY "S4: translation correctness and loop ordering", `todo_translation_screen_prototype.md`);
review page item 21 (port of the screen to the fused C path); Track E Task G (BFS pilot,
`findings_eval_reduction_2026-10.md`); Track X item X1 and `todo_500grain_bfs_end_to_end.md`.
