"""Task C: CandidateDedup in the level loop of ``AdaptiveVoxelReconstructor.reconstruct_voxel``:
the key parses, off is bit-identical to the default (no new recorder keys), the "discrete" record is
unchanged with it on, the filters act where the plan says (before both quick-MC paths, and before
the ranking cut) and ``SearchCandidate.clique`` is the FZ point."""

import math
from pathlib import Path
from typing import Any, Dict, List, Tuple

import numpy as np
import pytest
from scipy.spatial.transform import Rotation

from icenine import candidate_dedup as CD
from icenine import cost_functions as CF
from icenine.orientation_search import SearchParameters, run_discrete_search_spaced
from icenine.reconstructor import _get_voxel_vertices

from test_discrete_recorder import _args  # noqa: F401
from test_discrete_recorder import _rz

needs_c = pytest.mark.skipif(not CF._HAS_C_RASTERIZE, reason="C extension missing")


# ------------------------------------------------------------------------------ config and params
def test_config_key_default_off_and_parse(tmp_path: Any) -> None:
    from icenine.config_file import ConfigFile

    base = Path(__file__).resolve().parents[1] / "scripts" / "phase_d" / "configs"
    text = (base / "ReconstructPhaseD_mc_clean.config").read_text()
    sp0 = SearchParameters.from_config(
        ConfigFile.from_file(str(base / "ReconstructPhaseD_mc_clean.config"))
    )
    assert sp0.candidate_dedup == "off"
    f = tmp_path / "on.config"
    f.write_text(text + "\nCandidateDedup clique=cpp:2,pre=0.5,post=0.25\n")
    sp = SearchParameters.from_config(ConfigFile.from_file(str(f)))
    spec = CD.DedupSpec.parse(sp.candidate_dedup)
    assert (spec.clique, spec.clique_k, spec.pre, spec.post) == ("cpp", 2, 0.5, 0.25)
    f.write_text(text + "\nCandidateDedup clique=banana:1\n")
    with pytest.raises(ValueError, match="CandidateDedup"):
        ConfigFile.from_file(str(f))
    with pytest.raises(ValueError, match="CandidateDedup"):
        SearchParameters(candidate_dedup="pre=-1")


def test_clique_is_the_fz_point_index() -> None:
    out = run_discrete_search_spaced(**_args())
    assert len(out) > 0
    assert {c.clique for c in out} <= set(range(6))
    # each candidate came from local_grid[j] @ fz[clique] (angles 0.002 j + 0.01 clique, j < 5)
    for c in out:
        ang = math.atan2(c.orientation[1, 0], c.orientation[0, 0])
        j = (ang - 0.01 * c.clique) / 0.002
        assert abs(j - round(j)) < 1e-6 and 0 <= round(j) < 5


# ------------------------------------------------------------------------------ the level loop
def _recon(**kw: Any) -> Tuple[Any, Any, np.ndarray]:
    """The ThreeVoxels toy problem (one level) with eight FZ points within 2 deg of the truth, so
    that the candidates of neighbouring cliques overlap."""
    from test_findoptimal_refactor import _build

    rec, voxel, R_true = _build()
    rng = np.random.default_rng(3)
    jit = Rotation.from_rotvec(rng.normal(size=(8, 3)) * np.radians(1.2)).as_matrix()
    rec.setup.fz_orientations = np.concatenate([jit @ R_true, rec.setup.fz_orientations])
    for k, v in kw.items():
        setattr(rec.setup.search_params, k, v)
    return rec, voxel, R_true


def _run(**kw: Any) -> Tuple[Any, List[Tuple[str, Dict[str, Any]]], Any]:
    rec, voxel, R_true = _recon(**kw)
    events: List[Tuple[str, Dict[str, Any]]] = []
    rec.recorder = lambda name, data: events.append((name, data))
    res = rec.reconstruct_voxel(
        _get_voxel_vertices(voxel), voxel.phase, rng=np.random.default_rng(7)
    )
    return res, events, rec


def test_off_is_bit_identical_to_the_default_and_adds_no_keys() -> None:
    res0, ev0, _ = _run()
    res1, ev1, _ = _run(candidate_dedup="off")
    assert np.array_equal(res0.orientation, res1.orientation) and res0.cost == res1.cost
    assert [n for n, _ in ev0] == [n for n, _ in ev1]
    for (_, a), (_, b) in zip(ev0, ev1):
        assert a.keys() == b.keys()
        assert not any(k.startswith("dd_") for k in a)


def test_c0_clique_filter_runs_the_quick_mc_on_fewer_candidates() -> None:
    _, base, _ = _run()
    res, ev, rec = _run(candidate_dedup="clique=cpp:1")
    disc0, disc1 = base[0][1], ev[0][1]
    assert np.array_equal(disc0["R"], disc1["R"])  # the "discrete" record is unchanged
    quick = ev[1][1]
    n = len(disc1["R"])
    assert quick["dd_n_discrete"] == n and quick["dd_n_mc"] < n
    assert len(quick["perm"]) == quick["dd_n_mc"] and set(quick["perm"]) <= set(range(n))
    assert quick["dd_n_distinct"] == -1
    # at most one candidate per FZ point (8 jittered + 3 original), and the quarter is counted
    # from the discrete list, capped by what ran
    assert quick["dd_n_mc"] <= 11
    assert quick["n_keep"] == min(max(1, int(n * rec.keep_fraction)), quick["dd_n_mc"])
    assert np.isfinite(res.cost)


def test_c1_dedup_keeps_the_greedy_leaders_of_the_discrete_list() -> None:
    _, base, rec0 = _run()
    sym = np.asarray(
        __import__("icenine.orientation_search", fromlist=["x"]).get_symmetry_quaternions(
            rec0.setup.exp_setup.get_sample_symmetry()
        )
    )
    disc = base[0][1]
    t = 0.5
    kept, _ = CD.greedy_nms(disc["R"], t * disc["diameter"], sym)  # list is in score order
    _, ev, _ = _run(candidate_dedup=f"pre={t}")
    quick = ev[1][1]
    assert quick["dd_n_mc"] == len(kept) < len(disc["R"])
    assert sorted(quick["perm"].tolist()) == sorted(kept.tolist())


def test_c2_orders_distinct_candidates_first_and_never_keeps_a_passed_over_one() -> None:
    _, ev, rec = _run(candidate_dedup="post=1.0")
    quick = ev[1][1]
    assert quick["dd_n_mc"] == quick["dd_n_discrete"]  # nothing is dropped before the quick MC
    nk, nd = quick["n_keep"], quick["dd_n_distinct"]
    assert 1 <= nk <= nd
    from icenine.orientation_search import get_symmetry_quaternions

    sym = get_symmetry_quaternions(rec.setup.exp_setup.get_sample_symmetry())
    disc_d = ev[0][1]["diameter"]
    kept, _ = CD.greedy_nms(quick["R"][:nd], 1.0 * disc_d, sym)
    assert len(kept) == nd  # the first nd are mutually distinct at the tolerance
    assert len(CD.greedy_nms(quick["R"][:nk], 1.0 * disc_d, sym)[0]) == nk


@needs_c
def test_the_batched_quick_mc_sees_the_filtered_list() -> None:
    res, ev, rec = _run(candidate_dedup="clique=cpp:1", batched_quick_mc=True)
    quick = ev[1][1]
    (evaluator,) = rec._quick_mc_evaluators.values()
    assert evaluator.counts["orient_compact"] == 12 * quick["dd_n_mc"]
    assert quick["dd_n_mc"] < quick["dd_n_discrete"]
    assert np.isfinite(res.cost)


def _quick(spec: str) -> Tuple[Dict[str, Any], float]:
    _, ev, rec = _run(candidate_dedup=spec)
    (quick,) = [d for n, d in ev if n == "quick_mc"]  # the toy problem has one level
    return quick, float(rec.keep_fraction)


def test_quarter_filtered_counts_the_quarter_from_the_list_the_filter_left() -> None:
    # pre: quarter=orig counts from the discrete list (capped by what ran), quarter=filtered
    # from the shorter quick-MC list
    q_orig, kf = _quick("pre=0.5")
    q_filt, _ = _quick("pre=0.5,quarter=filtered")
    n, m = q_orig["dd_n_discrete"], q_orig["dd_n_mc"]
    assert (q_filt["dd_n_discrete"], q_filt["dd_n_mc"]) == (n, m) and m < n
    assert q_orig["n_keep"] == min(max(1, int(n * kf)), m)
    assert q_filt["n_keep"] == max(1, int(m * kf))
    assert q_filt["n_keep"] < q_orig["n_keep"]  # the toy case tells the two rules apart


def test_quarter_filtered_with_post_counts_the_quarter_from_the_distinct_candidates() -> None:
    q_orig, kf = _quick("post=0.25")
    q_filt, _ = _quick("post=0.25,quarter=filtered")
    n, nd = q_orig["dd_n_discrete"], q_orig["dd_n_distinct"]
    assert q_orig["dd_n_mc"] == n and q_filt["dd_n_distinct"] == nd and 1 <= nd < n
    assert q_orig["n_keep"] == min(max(1, int(n * kf)), nd)
    assert q_filt["n_keep"] == min(max(1, int(nd * kf)), nd)
    assert q_filt["n_keep"] < q_orig["n_keep"]
