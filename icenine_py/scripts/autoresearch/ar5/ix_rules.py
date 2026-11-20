"""AR5 step 26 (task ar5z): rule ``IX``, the gates and the predictions, as constants (pure).

``ix_freeze`` writes these into ``ar5/ix/frozen.json`` before any exactness or timing run.
"""

from typing import Any, Dict

RULE = (
    "Rule IX: a new indexer returns, for every voxel, exactly the old index_voxel output "
    "(rotations, explained-ray and vote counts, attached rays, in the same order). Allowed: "
    "vectorizing, reordering, reusing arrays, skipping work that provably cannot change the "
    "output. Not allowed: any change of IndexParams (tolerances, n_seeds 150, max_pairs 600,000, "
    "the pair subsample generator), of top_k 20, or of anything after indexing. No truth, no "
    "threshold, nothing crystal-specific; any symmetry group works."
)
SEED_GENERATOR = 26  # 300 seed steps per sample for held3 to held6, grad, big
N_PER_SAMPLE = 300
FULL_SAMPLES = ("mid", "held")  # every stored seed voxel
SAMPLED = ("held3", "held4", "held5", "held6", "grad", "big")
TIMING_SAMPLES = ("mid", "held")
BLOCK = 100  # voxels per alternating block of the index-only timing

GATES: Dict[str, Any] = dict(
    z0=(
        "quiet: the old indexer on the 200 mid voxels (every 8th stored seed step), scaled to the "
        "1,626 stored steps, is at least 30% of a quiet J5SD estimate on mid (stored evaluations "
        "by stage at the session's quiet ms per evaluation)"
    ),
    e1=(
        "exactness, index only: old and new on every exactness voxel (5,311 of mid and held, 300 "
        "each of held3 to held6, grad, big). Tier 1: bit-identical outputs. Tier 2, only if tier "
        "1 fails on at most 0.1% of voxels: same order, counts and attached rays, rotation "
        "entries within 1e-12. Anything else: stop"
    ),
    e2=(
        "exactness, live: the J5SD script arm with the new indexer on mid and held. Tier 1: maps "
        "and counts equal stored J5SD exactly (7,383,368 and 17,114,589). Tier 2: evaluations "
        "within 0.5% on each; at most 0.1% of voxels move over 0.001 degree; against stored J5SD "
        "on mid and held pooled, broken minus fixed at most 3 per sample and median error at "
        "most 0.001 degree higher"
    ),
    w1=(
        "seconds, quiet, index only: summed indexing seconds over the 5,311 voxels at least 50% "
        "lower than the old indexer's, same session, old and new in alternating blocks of 100 "
        "voxels"
    ),
    w2=(
        "seconds, quiet, whole run: J5SD old and new, each run alone, on mid and held: pooled "
        "wall seconds (sum of both) at least 25% lower, and at least 15% lower on each sample"
    ),
    tier2_only="held7 house rule (not built or read under tier 1)",
)
PREDICTIONS = dict(
    e1="tier 1 holds",
    w1="indexing 2 to 5 times faster",
    w2="held 35% to 55% less wall, mid 25% to 45%",
    evals="evaluation counts unchanged everywhere",
)
TIMING_RULE = (
    "a timed run needs a passing timing_preflight --require first; if it refuses, wait and retry "
    "up to 3 times 20 minutes apart; a contended second never enters a gate; a gate that could "
    "not be timed quietly is 'not judged: machine not quiet'"
)
