"""AR5 step 25 phase 0: the small models and guards of ``es_sizing.py`` (no sample is opened)."""

import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "scripts" / "autoresearch" / "ar5"))
import es_sizing as E  # noqa: E402


def test_auc_known_values() -> None:
    y = np.array([0, 0, 1, 1])
    assert E.auc(y, np.array([0.1, 0.2, 0.3, 0.4])) == 1.0
    assert E.auc(y, np.array([0.4, 0.3, 0.2, 0.1])) == 0.0
    assert E.auc(y, np.full(4, 0.5)) == 0.5
    assert np.isnan(E.auc(np.zeros(3, int), np.arange(3.0)))


def test_logistic_separates_and_tree_finds_split() -> None:
    rng = np.random.default_rng(0)
    X = rng.normal(size=(400, 3))
    y = (X[:, 0] + 0.2 * rng.normal(size=400) > 0).astype(int)
    w = E.fit_logistic(X, y)
    assert E.auc(y, E.predict_logistic(w, X)) > 0.95
    tree = E.fit_tree(X, y)
    assert E.auc(y, E.predict_tree(tree, X)) > 0.9
    assert isinstance(tree, tuple) and tree[0] == 0


def test_guard_refuses_held_out_names() -> None:
    for name in ("held2", "held7", "fine7"):
        with pytest.raises(SystemExit):
            E._guard(name)
    E._guard("held3")


def _write_cache(root: Path, s: str, e1: list, mc: list) -> None:
    import gzip
    import json

    d, r = root / f"J5D_{s}", root / f"stored_J5SD_{s}"
    d.mkdir(parents=True)
    r.mkdir(parents=True)
    with gzip.open(d / "e1_log.json.gz", "wt") as fh:
        json.dump(e1, fh)
    with gzip.open(d / "mc_traces.json.gz", "wt") as fh:
        json.dump(mc, fh)
    stats = dict(n_evals_seed=10, n_evals_reseed=0, n_evals_neighbor=0, n_evals_revisit=0)
    (d / "run.json").write_text(json.dumps(dict(stats=stats)))
    (r / "repair.json").write_text(json.dumps(dict(pass_evals=0, n_evals=0)))


def test_extract_pins_D_label_and_check_evaluation(tmp_path: Path) -> None:
    """Step 0 reaches D (never), step 1 reaches D (rescued), step 2 has a CMA end at 0.95 (no D)."""

    def step(i: int, cma: list, stop: int) -> dict:
        fits = [[p, -1, "C", h, 1.0 - h, 100] for p, h in enumerate(cma)]
        fits.append([4, 0, "M", 0.7, 0.3, 700])
        return dict(voxel=i, costed=10, fits=fits, reseed=False, stop=stop, evals=2000)

    e1 = [
        step(0, [0.5, 0.6, 0.7, 0.4], -1),
        step(1, [0.5, 0.6, 0.7, 0.4], 3),
        step(2, [0.5, 0.96, 0.7, 0.4], -1),
    ]
    mc = [
        dict(
            step=i + 1,
            j=0,
            voxel=i,
            evals=700,
            checks=[[500, 0.6]],
            trace=[[1, 2.0], [300, 1.5], [600, 1.0]],
        )
        for i in range(3)
    ]
    _write_cache(tmp_path, "tmp", e1, mc)
    rows, info = E.extract(tmp_path, "tmp")
    assert info["steps"] == 3 and info["bad_map"] == 0
    assert [r["y"] for r in rows] == [0, 1]  # step 2 has a 0.95 end before D
    # costed 10 + 4 CMA fits x 100 + the paused refinement: 500 evaluations + the check (501)
    assert rows[0]["evD"] == 10 + 400 + 501
    assert rows[0]["after"] == 2000 - 911
    # the trace is cut at the check (evaluation 500): the cost there is 1.5, the later gain
    # (evaluation 600, cost 1.0) is not used
    assert rows[0]["gains_best"] == 2
    assert rows[0]["ng2_best"] == 1  # one gain after evaluation 208 (the one at 300)
    assert rows[0]["cost1"] == pytest.approx(0.3)  # best CMA end (hit 0.7) beats the MC cost 1.5
