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
