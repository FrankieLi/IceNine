"""
Stage timer and cost-evaluation counters, installed from scripts by patching module attributes
(nothing in icenine/ is modified).

* StageTimer.stage(name): context manager; stages nest, and both the inclusive time and the
  exclusive time (inclusive minus nested stages) are accumulated per name.
* StageTimer.patch(owner, attr, name): replace owner.attr (a module function or a class method) by a
  wrapper that runs it inside stage(name). Wrappers return the original result untouched, so a
  patched run is bit-identical to an unpatched one. `name` may be a callable (args, kwargs) -> str
  to pick the stage per call (e.g. quick MC vs FindOptimal by max_mc_steps). Patching the same
  (owner, attr) twice raises. staticmethod / classmethod descriptors are preserved, and uninstall
  restores the exact original object (deleting the attribute again when it was inherited).
* StageTimer.count_evaluations(cls): patch cls.evaluate (VoxelCostFunction) so every call is counted
  per INSTANCE (weakly referenced) and its time accumulated under "evaluate". Label an instance with
  timer.label(obj, "global" | "local") to read counts by role afterwards (evaluations_by_role()).
* `with timer.installed():` uninstalls every patch on exit.

IMPORTANT for by-name imports: `from m import f` binds f in the importing module, so patching m.f
does not affect it. Patch the name where it is looked up, e.g.
`timer.patch(icenine.reconstructor, "run_discrete_search_spaced", "discrete")` for the
reconstructor.

Shared by scripts/nn_hybrid/run.py (Task 1) and the run-time profiling of Task 3.
"""

import contextlib
import functools
import inspect
import time
import weakref
from collections import defaultdict
from typing import Any, Callable, Dict, Iterator, List, Tuple, Union

StageName = Union[str, Callable[[tuple, dict], str]]
_MISSING = object()


class StageTimer:
    def __init__(self) -> None:
        self.inclusive: Dict[str, float] = defaultdict(float)
        self.exclusive: Dict[str, float] = defaultdict(float)
        self.calls: Dict[str, int] = defaultdict(int)
        self._stack: List[List[Any]] = []  # [name, start, nested seconds]
        self._patches: List[Tuple[Any, str, Any]] = []  # (owner, attr, original or _MISSING)
        self._roles: "weakref.WeakKeyDictionary[Any, str]" = weakref.WeakKeyDictionary()
        self._evals: "weakref.WeakKeyDictionary[Any, int]" = weakref.WeakKeyDictionary()
        self._evals_base: Dict[str, int] = {}

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
        """Clear the times and calls, and make the evaluation counts start from zero again."""
        self.inclusive.clear()
        self.exclusive.clear()
        self.calls.clear()
        self._evals_base = dict(self._by_role())

    def snapshot(self) -> Dict[str, Any]:
        return dict(
            inclusive=dict(self.inclusive),
            exclusive=dict(self.exclusive),
            calls={k: float(v) for k, v in self.calls.items()},
            evaluations=self.evaluations_by_role(),
        )

    def delta_since(self, before: Dict[str, Any], kind: str = "exclusive") -> Dict[str, float]:
        """Seconds ("exclusive" or "inclusive") accumulated since `before` (a snapshot())."""
        now = self.exclusive if kind == "exclusive" else self.inclusive
        b = before[kind]
        return {k: v - b.get(k, 0.0) for k, v in now.items() if v - b.get(k, 0.0) != 0}

    def evaluations_since(self, before: Dict[str, Any]) -> Dict[str, int]:
        now, b = self.evaluations_by_role(), before["evaluations"]
        return {k: v - b.get(k, 0) for k, v in now.items()}

    # -- patching ----------------------------------------------------------------------------

    def _install(self, owner: Any, attr: str, make: Callable[[Any], Any]) -> None:
        if any(o is owner and a == attr for o, a, _ in self._patches):
            raise RuntimeError(f"{owner!r}.{attr} is already patched by this timer")
        static = inspect.getattr_static(owner, attr)
        orig_fn = static.__func__ if isinstance(static, (staticmethod, classmethod)) else static
        wrapper = make(orig_fn)
        if isinstance(static, staticmethod):
            wrapper = staticmethod(wrapper)
        elif isinstance(static, classmethod):
            wrapper = classmethod(wrapper)
        own = (
            vars(owner).get(attr, _MISSING) if hasattr(owner, "__dict__") else getattr(owner, attr)
        )
        self._patches.append((owner, attr, own))
        setattr(owner, attr, wrapper)

    def patch(self, owner: Any, attr: str, name: StageName) -> None:
        def make(orig: Callable[..., Any]) -> Callable[..., Any]:
            @functools.wraps(orig)
            def wrapper(*a: Any, **k: Any) -> Any:
                with self.stage(name(a, k) if callable(name) else name):
                    return orig(*a, **k)

            return wrapper

        self._install(owner, attr, make)

    def count_evaluations(self, cls: Any, name: str = "evaluate") -> None:
        timer = self

        def make(orig: Callable[..., Any]) -> Callable[..., Any]:
            @functools.wraps(orig)
            def wrapper(self: Any, *a: Any, **k: Any) -> Any:
                timer._evals[self] = timer._evals.get(self, 0) + 1
                with timer.stage(name):
                    return orig(self, *a, **k)

            return wrapper

        self._install(cls, "evaluate", make)

    def uninstall(self) -> None:
        for owner, attr, own in reversed(self._patches):
            if own is _MISSING:
                delattr(owner, attr)  # the attribute was inherited
            else:
                setattr(owner, attr, own)
        self._patches.clear()

    @contextlib.contextmanager
    def installed(self) -> Iterator["StageTimer"]:
        try:
            yield self
        finally:
            self.uninstall()

    # -- evaluation counters ---------------------------------------------------------------------

    def label(self, obj: Any, role: str) -> None:
        self._roles[obj] = role

    def evaluations(self, obj: Any) -> int:
        return int(self._evals.get(obj, 0))

    def _by_role(self) -> Dict[str, int]:
        out: Dict[str, int] = defaultdict(int)
        for obj, n in list(self._evals.items()):
            out[self._roles.get(obj, "unlabelled")] += n
        return dict(out)

    def evaluations_by_role(self) -> Dict[str, int]:
        """Evaluation counts by label since the last reset() (live instances only)."""
        return {k: v - self._evals_base.get(k, 0) for k, v in self._by_role().items()}
