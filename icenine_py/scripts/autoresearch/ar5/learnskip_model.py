"""AR5 step 18 (task ar5r): the boosted-tree model of the learned skip: training (scikit-learn,
lazily imported), export to JSON and a plain-Python evaluator (no numpy, no scikit-learn at run
time).

``train`` is the frozen model: ``HistGradientBoostingClassifier(max_iter=300, learning_rate=0.1,
random_state=0)``, all else at the scikit-learn 1.6.1 defaults. ``export`` writes the trees
(feature,
threshold, children, leaf value, the baseline). ``PlainModel.score`` is the predicted chance that a
fit is accepted. A model call is not a cost-function evaluation: it is counted as zero evaluations.
"""

import json
import math
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

import numpy as np

import learnskip_rules as LR


def train(X: np.ndarray, y: np.ndarray) -> Any:
    from sklearn.ensemble import HistGradientBoostingClassifier

    m = HistGradientBoostingClassifier(
        max_iter=LR.MODEL_SPEC["max_iter"],
        learning_rate=LR.MODEL_SPEC["learning_rate"],
        random_state=LR.MODEL_SPEC["random_state"],
    )
    m.fit(X, np.asarray(y, bool))
    return m


def export(
    m: Any, threshold: Optional[float], extra: Optional[Dict[str, Any]] = None
) -> Dict[str, Any]:
    """The trees of a fitted classifier as a JSON-able dict (``threshold`` None: no skip)."""
    trees = []
    for (p,) in m._predictors:
        n = p.nodes
        trees.append(
            dict(
                f=n["feature_idx"].tolist(),
                t=n["num_threshold"].tolist(),
                l=n["left"].tolist(),
                r=n["right"].tolist(),
                leaf=n["is_leaf"].astype(int).tolist(),
                v=n["value"].tolist(),
            )
        )
    out: Dict[str, Any] = dict(
        baseline=float(np.ravel(m._baseline_prediction)[0]),
        n_rounds=int(m.n_iter_),
        threshold=threshold,
        trees=trees,
    )
    out.update(extra or {})
    return out


def dump(model: Dict[str, Any], path: Path) -> None:
    path.write_text(json.dumps(model, separators=(",", ":")))


class PlainModel:
    """The exported trees, evaluated in plain Python."""

    def __init__(self, d: Dict[str, Any]) -> None:
        self.d = d
        self.baseline = float(d["baseline"])
        self.threshold: Optional[float] = d["threshold"]
        self.trees: List[Tuple[Any, ...]] = [
            (t["f"], t["t"], t["l"], t["r"], t["leaf"], t["v"]) for t in d["trees"]
        ]

    @classmethod
    def load(cls, path: Path) -> "PlainModel":
        return cls(json.loads(Path(path).read_text()))

    def raw(self, x: Sequence[float]) -> float:
        total = self.baseline
        for f, th, left, right, leaf, v in self.trees:
            n = 0
            while not leaf[n]:
                n = left[n] if x[f[n]] <= th[n] else right[n]
            total += v[n]
        return total

    def score(self, x: Sequence[float]) -> float:
        z = self.raw(x)
        if z >= 0:
            return 1.0 / (1.0 + math.exp(-z))
        e = math.exp(z)
        return e / (1.0 + e)

    def skips(self, score: float) -> bool:
        """Skip a fit whose score is under the threshold (no threshold: never)."""
        return self.threshold is not None and score < self.threshold

    def score_many(self, X: np.ndarray) -> np.ndarray:
        return np.array([self.score(row) for row in np.asarray(X, np.float64).tolist()])
