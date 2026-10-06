"""
Stage timer and cost-evaluation counters, installed from scripts by patching module attributes
(nothing in icenine/ is modified).

* StageTimer.stage(name): context manager; stages nest, and both the inclusive time and the
  exclusive time (inclusive minus nested stages) are accumulated per name.
* StageTimer.patch(owner, attr, name): replace owner.attr (a module function or a class method) by a
  wrapper that runs it inside stage(name). Wrappers return the original result untouched, so a
  patched run is bit-identical to an unpatched one.
* StageTimer.count_evaluations(cls): patch cls.evaluate (VoxelCostFunction) so every call is counted
  per INSTANCE (attribute `_st_evals`) and its time accumulated under "evaluate". Label an instance
  with timer.label(obj, "global" | "local") to read counts by role afterwards.

Shared by scripts/nn_hybrid/run.py (Task 1) and the run-time profiling of Task 3.
"""

import contextlib
import functools
import time
from collections import defaultdict
from typing import Any, Callable, Dict, Iterator, List, Optional, Tuple


class StageTimer:
    def __init__(self) -> None:
        self.inclusive: Dict[str, float] = defaultdict(float)
        self.exclusive: Dict[str, float] = defaultdict(float)
        self.calls: Dict[str, int] = defaultdict(int)
        self._stack: List[List[Any]] = []  # [name, start, nested seconds]
        self._patches: List[Tuple[Any, str, Any]] = []
        self.labels: Dict[int, str] = {}
        self._instances: Dict[int, Any] = {}

    # -- timing ------------------------------------------------------------------------------

    @contextlib.contextmanager
    def stage(self, name: str) -> Iterator[None]:
        frame: List[Any] = [name, time.perf_counter(), 0.0]
        self._stack.append(frame)
        try:
            yield
        finally:
            dt = time.perf_counter() - frame[1]
            self._stack.pop()
            self.inclusive[name] += dt
            self.exclusive[name] += dt - frame[2]
            self.calls[name] += 1
            if self._stack:
                self._stack[-1][2] += dt

    def reset(self) -> None:
        self.inclusive.clear()
        self.exclusive.clear()
        self.calls.clear()

    def snapshot(self) -> Dict[str, Dict[str, float]]:
        return dict(
            inclusive=dict(self.inclusive),
            exclusive=dict(self.exclusive),
            calls={k: float(v) for k, v in self.calls.items()},
        )

    def delta_since(
        self, before: Dict[str, Dict[str, float]], kind: str = "exclusive"
    ) -> Dict[str, float]:
        """Seconds ("exclusive" or "inclusive") accumulated since `before` (a snapshot())."""
        now = self.exclusive if kind == "exclusive" else self.inclusive
        b = before[kind]
        return {k: v - b.get(k, 0.0) for k, v in now.items() if v - b.get(k, 0.0) != 0}

    # -- patching ----------------------------------------------------------------------------

    def patch(self, owner: Any, attr: str, name: str) -> None:
        orig = getattr(owner, attr)

        @functools.wraps(orig)
        def wrapper(*a: Any, **k: Any) -> Any:
            with self.stage(name):
                return orig(*a, **k)

        self._patches.append((owner, attr, orig))
        setattr(owner, attr, wrapper)

    def count_evaluations(self, cls: Any, name: str = "evaluate") -> None:
        orig = cls.evaluate
        timer = self

        @functools.wraps(orig)
        def wrapper(self: Any, *a: Any, **k: Any) -> Any:
            self._st_evals = getattr(self, "_st_evals", 0) + 1
            timer._instances[id(self)] = self
            with timer.stage(name):
                return orig(self, *a, **k)

        self._patches.append((cls, "evaluate", orig))
        cls.evaluate = wrapper

    def label(self, obj: Any, role: str) -> None:
        self.labels[id(obj)] = role
        self._instances[id(obj)] = obj

    def evaluations(self, obj: Any) -> int:
        return int(getattr(obj, "_st_evals", 0))

    def evaluations_by_role(self) -> Dict[str, int]:
        out: Dict[str, int] = defaultdict(int)
        for i, obj in self._instances.items():
            out[self.labels.get(i, "unlabelled")] += int(getattr(obj, "_st_evals", 0))
        return dict(out)

    def uninstall(self) -> None:
        for owner, attr, orig in reversed(self._patches):
            setattr(owner, attr, orig)
        self._patches.clear()
