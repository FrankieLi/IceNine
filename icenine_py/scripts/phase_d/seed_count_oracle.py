"""Estimate how many no-start seed searches a Phase D BFS makes, with an oracle in place of the
fits: a seed succeeds with probability 1 - p_fail; an expansion accepts a neighbour iff it is in
the same grain (else REFIT). Uses the BFS rules of BFSReconstruction (random voxel order, seeds
only from NOT_VISITED voxels, neighbours within 2 side lengths via a cKDTree, REFIT voxels never
become seeds). Optionally the revisit option (REFIT neighbours of an expansion are retried, so
voxels that are REFIT only because another grain reached them first are recovered; this does not
change the seed count in the oracle, it changes the unresolved count).

Usage (from icenine_py/): uv run python scripts/phase_d/seed_count_oracle.py
"""

import json
from pathlib import Path

import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
from scipy.spatial import cKDTree

HERE = Path(__file__).resolve().parent
EX = HERE.parents[2] / "Examples" / "Example2.ManyGrains"
SIDE = 0.009375


def main() -> None:
    raw = np.loadtxt(EX / "SimInput" / "rand_500grains_1mm_neworient_s0.mic", skiprows=1)
    pos = raw[:, :2]
    grain = np.load(EX / "SimInput" / "rand_500grains_1mm_neworient_s0_grainmap.npy")
    n = len(pos)
    tree = cKDTree(pos)
    nbrs = [
        np.array([j for j in nb if j != i])
        for i, nb in enumerate(tree.query_ball_point(pos, 2.0 * SIDE))
    ]
    out = {"n_voxels": n, "n_grains": int(grain.max() + 1)}
    for name, rad in (("edge_vertex_1.01side", 1.01 * SIDE), ("bfs_2side", 2.0 * SIDE)):
        pairs = tree.query_pairs(rad, output_type="ndarray")
        same = grain[pairs[:, 0]] == grain[pairs[:, 1]]
        p = pairs[same]
        m = coo_matrix((np.ones(len(p)), (p[:, 0], p[:, 1])), shape=(n, n))
        out[f"pieces_{name}"] = int(connected_components(m, directed=False)[0])
    NV, VI, FI, RF = 0, 1, 2, 3
    res = {}
    for p_fail in (0.0, 0.05, 0.1, 0.2, 0.3):
        seeds = []
        for rep in range(20):
            rng = np.random.default_rng(100 + rep)
            state = np.zeros(n, np.int8)
            n_seeds = 0
            for s in rng.permutation(n):
                if state[s] != NV:
                    continue
                n_seeds += 1
                if rng.random() < p_fail:
                    state[s] = RF
                    continue
                state[s] = FI
                q = [s]
                head = 0

                # insert_seed(s)
                def insert(i: int) -> None:
                    for j in nbrs[i]:
                        if state[j] == NV:
                            state[j] = VI
                            q.append(j)

                insert(s)
                head = 1
                while head < len(q):
                    j = q[head]
                    head += 1
                    if state[j] == FI:
                        continue
                    # an expansion centre: the grain of the voxel that queued j is not tracked;
                    # oracle: accepted iff j is in the same grain as the seed
                    if grain[j] == grain[s]:
                        state[j] = FI
                        insert(j)
                    else:
                        state[j] = RF
            seeds.append(n_seeds)
        res[str(p_fail)] = dict(
            mean=float(np.mean(seeds)), min=int(min(seeds)), max=int(max(seeds))
        )
    out["oracle_seeds_by_p_fail"] = res
    (HERE.parents[1] / "benchmarks" / "phase_d_seed_diag" / "seed_count_oracle.json").write_text(
        json.dumps(out, indent=1) + "\n"
    )
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
