"""AR5 step 8 (task ar5h): aligned random draws in the script layer (library untouched).

A test device only: the library keeps its own draws.

In a whole-layer run one generator is shared by every local optimization, grid fallback and pass
try, so a change at one voxel shifts the draws of every later call and moves other voxels (the
lesson of step 7). Aligned draws give every such call its own generator, seeded from three things
only: the voxel, the kind of call and the number of calls of that kind made at that voxel so far.
The same (voxel, kind, count) then gets the same draws whatever ran before, in any order.

* ``voxel_key(vertices)``: a 64-bit hash of the voxel's vertex array (the layer's calls carry the
  vertices, not the voxel number); a plain integer is used as is (the pass knows the number).
* ``DrawBook.rng(key, kind)``: the next generator for that voxel and kind (counts per book).
* ``aligned_local_optimization``: wraps a ``local_optimization``-like callable (keyword call as the
  reconstructor makes it); the kind is the name of the calling function.
* ``aligned_reconstruct_voxel``: the same for the grid fallback (kind ``grid``).
* ``aligned_try_fit``: the same for the boundary pass's ``try_fit`` (kind ``pass``).

The seed step's refinements already use a fresh ``default_rng(0)`` each (library code), so they are
aligned as they stand. Nothing here changes a value other than the generator passed in.
"""

import hashlib
import sys
import zlib
from typing import Any, Callable, Dict, Tuple

import numpy as np


def voxel_key(vertices: Any) -> int:
    """A 64-bit integer key of a voxel from its vertices (array or tensor) or its number."""
    if isinstance(vertices, (int, np.integer)):
        return int(vertices)
    arr = vertices.detach().cpu().numpy() if hasattr(vertices, "detach") else np.asarray(vertices)
    a = np.ascontiguousarray(arr, dtype=np.float64)
    return int.from_bytes(hashlib.blake2b(a.tobytes(), digest_size=8).digest(), "little")


def kind_id(kind: str) -> int:
    return zlib.crc32(kind.encode("utf-8"))


class DrawBook:
    """Counts the calls per (voxel, kind) and hands out the generators."""

    def __init__(self) -> None:
        self.counts: Dict[Tuple[int, str], int] = {}
        self.n_issued = 0

    def rng(self, key: int, kind: str) -> np.random.Generator:
        n = self.counts.get((key, kind), 0)
        self.counts[(key, kind)] = n + 1
        self.n_issued += 1
        return make_rng(key, kind, n)


def make_rng(key: int, kind: str, count: int) -> np.random.Generator:
    """The generator of the ``count``-th (from 0) call of ``kind`` at the voxel with ``key``."""
    return np.random.default_rng(np.random.SeedSequence([int(key), kind_id(kind), int(count)]))


def aligned_local_optimization(book: DrawBook, inner: Callable[..., Any]) -> Callable[..., Any]:
    """``inner(voxel_vertices=, phase_index=, initial_orientation=, rng=, **kw)`` with the ``rng``
    replaced by the aligned generator of (voxel, calling function, count)."""

    def wrapped(
        voxel_vertices: Any,
        phase_index: int,
        initial_orientation: Any,
        rng: Any = None,
        **kw: Any,
    ) -> Any:
        kind = sys._getframe(1).f_code.co_name
        g = book.rng(voxel_key(voxel_vertices), kind)
        return inner(
            voxel_vertices=voxel_vertices,
            phase_index=phase_index,
            initial_orientation=initial_orientation,
            rng=g,
            **kw,
        )

    return wrapped


def aligned_reconstruct_voxel(book: DrawBook, inner: Callable[..., Any]) -> Callable[..., Any]:
    """The grid fallback (``reconstruct_voxel(voxel_vertices=, phase_index=, rng=)``)."""

    def wrapped(voxel_vertices: Any, phase_index: int, rng: Any = None, **kw: Any) -> Any:
        g = book.rng(voxel_key(voxel_vertices), "grid")
        return inner(voxel_vertices=voxel_vertices, phase_index=phase_index, rng=g, **kw)

    return wrapped


def aligned_try_fit(book: DrawBook, inner: Callable[..., Any]) -> Callable[..., Any]:
    """The boundary pass's ``try_fit(hooks, bfs, v, R0, rng)``; ``v`` is the voxel number."""

    def wrapped(hooks: Any, bfs: Any, v: int, R0: Any, rng: Any = None) -> Any:
        return inner(hooks, bfs, v, R0, book.rng(int(v), "pass"))

    return wrapped
