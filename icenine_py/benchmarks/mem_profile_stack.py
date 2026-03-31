#!/usr/bin/env python3
"""
Memory profiling script for MultiScaleImageStack construction.

Loads progressively more omega frames (10, 30, 60, 120, 180) and measures
RSS memory at each stage:
  1. After loading SparseImageStack (baseline)
  2. After _downsample_stack (scale 1: 4×)
  3. After _downsample_stack (scale 2: 8×)
  4. After _omega_blend (scale 1 + scale 2)

Also checks torch allocator cache separately from true RSS.
"""

import gc
import os
import sys
import tracemalloc
from pathlib import Path

import psutil
import torch
import torch.nn.functional as F

project_root = Path(__file__).parent.parent.parent
sys.path.insert(0, str(project_root / "icenine_py"))

EXAMPLE_DIR = project_root / "Examples" / "Example2.ThreeVoxels"
BASENAME = "3Grains.sim"
DATA_DIR = EXAMPLE_DIR / "ScatteringData_Python"


def rss_mb() -> float:
    proc = psutil.Process(os.getpid())
    return proc.memory_info().rss / 1024 / 1024


def torch_mb() -> str:
    alloc = torch.cuda.memory_allocated() / 1024 / 1024 if torch.cuda.is_available() else 0
    reserved = torch.cuda.memory_reserved() / 1024 / 1024 if torch.cuda.is_available() else 0
    # For CPU, report allocator cache if available
    try:
        cache = torch.cuda.memory_reserved() / 1024**2
    except Exception:
        cache = 0
    return f"cuda_alloc={alloc:.0f}MB" if torch.cuda.is_available() else "cpu"


def banner(msg: str):
    print(f"\n{'─'*60}")
    print(f"  {msg}")
    print(f"  RSS: {rss_mb():.0f} MB")
    print(f"{'─'*60}", flush=True)


def measure_stage(tag: str, before_mb: float):
    after = rss_mb()
    delta = after - before_mb
    print(f"  [{tag}]  RSS: {after:.0f} MB  (Δ {delta:+.0f} MB)", flush=True)
    return after


def run_profile(n_omega: int):
    from icenine.differentiable_cost import (
        ExperimentalImageStack,
        MultiScaleImageStack,
        SparseImageStack,
    )

    print(f"\n{'='*60}")
    print(f"  n_omega = {n_omega}  ({n_omega * 2} image files)")
    print(f"{'='*60}", flush=True)

    gc.collect()
    t0 = rss_mb()
    print(f"  Baseline RSS: {t0:.0f} MB")

    # ---- Stage 1: Load SparseImageStack ----
    sparse = SparseImageStack.from_image_directory(
        directory=str(DATA_DIR),
        basename=BASENAME,
        ext="d",
        serial_length=5,
        n_omega=n_omega,
        n_detectors=2,
        num_rows=2048,
        num_cols=2048,
        binary=True,
    )
    t1 = measure_stage(f"SparseImageStack (n_omega={n_omega})", t0)
    print(f"         sparse.memory_bytes = {sparse.memory_bytes/1024:.1f} KB")
    total_nnz = sum(c.shape[0] for c in sparse._pixel_coords)
    print(f"         total nnz pixels    = {total_nnz}")

    # ---- Stage 2: _downsample_stack factor=4 ----
    # Create a minimal MultiScaleImageStack-like object to test one step at a time
    ms = MultiScaleImageStack.__new__(MultiScaleImageStack)
    ms.downsample_factors = [1, 4, 8]
    ms.n_scales = 3
    ms.omega_window = 0
    ms.scales = []

    # Scale 0: pass-through (just reference the sparse stack)
    ms.scales.append(sparse)
    gc.collect()
    t_after_s0 = measure_stage("scale 0 (pass-through, no densify)", t1)

    # Scale 1: factor=4
    ds1 = ms._downsample_stack(sparse, factor=4)
    gc.collect()
    t_after_s1 = measure_stage("scale 1 (_downsample factor=4)", t_after_s0)
    sz1 = ds1.images.numel() * 4 / 1024 / 1024
    print(f"         scale1.images shape = {tuple(ds1.images.shape)}  "
          f"({sz1:.1f} MB tensor)")
    ms.scales.append(ds1)

    # Scale 2: factor=8
    ds2 = ms._downsample_stack(sparse, factor=8)
    gc.collect()
    t_after_s2 = measure_stage("scale 2 (_downsample factor=8)", t_after_s1)
    sz2 = ds2.images.numel() * 4 / 1024 / 1024
    print(f"         scale2.images shape = {tuple(ds2.images.shape)}  "
          f"({sz2:.1f} MB tensor)")
    ms.scales.append(ds2)

    # ---- Stage 3: omega_blend on scale 1 ----
    blended1 = ms._omega_blend(ds1, window=1)
    gc.collect()
    t_after_b1 = measure_stage("_omega_blend scale 1 (window=1)", t_after_s2)

    # ---- Stage 4: omega_blend on scale 2 ----
    blended2 = ms._omega_blend(ds2, window=1)
    gc.collect()
    t_after_b2 = measure_stage("_omega_blend scale 2 (window=1)", t_after_b1)

    # Summary
    print(f"\n  SUMMARY  (n_omega={n_omega})")
    print(f"    Baseline             : {t0:.0f} MB")
    print(f"    After SparseStack    : {t1:.0f} MB  Δ={t1-t0:+.0f}")
    print(f"    After scale1 (4×)   : {t_after_s1:.0f} MB  Δ={t_after_s1-t1:+.0f}")
    print(f"    After scale2 (8×)   : {t_after_s2:.0f} MB  Δ={t_after_s2-t_after_s1:+.0f}")
    print(f"    After blend1 (±1)   : {t_after_b1:.0f} MB  Δ={t_after_b1-t_after_s2:+.0f}")
    print(f"    After blend2 (±1)   : {t_after_b2:.0f} MB  Δ={t_after_b2-t_after_b1:+.0f}")
    print(f"    Total growth        : {t_after_b2-t0:+.0f} MB")

    # Keep tensors in scope to prevent GC
    del blended1, blended2, ds1, ds2, ms, sparse
    gc.collect()
    t_end = rss_mb()
    print(f"    After del+gc        : {t_end:.0f} MB  Δ={t_end-t0:+.0f}")

    return t_after_b2 - t0  # total growth


def run_tracemalloc_profile(n_omega: int = 30):
    """Python-level allocation breakdown using tracemalloc."""
    from icenine.differentiable_cost import MultiScaleImageStack, SparseImageStack

    print(f"\n{'='*60}")
    print(f"  tracemalloc profile  (n_omega={n_omega})")
    print(f"{'='*60}", flush=True)

    tracemalloc.start()
    snap0 = tracemalloc.take_snapshot()

    sparse = SparseImageStack.from_image_directory(
        directory=str(DATA_DIR),
        basename=BASENAME,
        ext="d",
        serial_length=5,
        n_omega=n_omega,
        n_detectors=2,
        num_rows=2048,
        num_cols=2048,
        binary=True,
    )
    snap1 = tracemalloc.take_snapshot()

    ms = MultiScaleImageStack.__new__(MultiScaleImageStack)
    ms.downsample_factors = [1, 4, 8]
    ms.n_scales = 3
    ms.omega_window = 0
    ms.scales = []
    ms.scales.append(sparse)
    ds1 = ms._downsample_stack(sparse, factor=4)
    ms.scales.append(ds1)
    ds2 = ms._downsample_stack(sparse, factor=8)
    ms.scales.append(ds2)

    snap2 = tracemalloc.take_snapshot()
    tracemalloc.stop()

    print("\n  Top 10 allocations (SparseStack load):")
    stats = snap1.compare_to(snap0, "lineno")
    for s in stats[:10]:
        if s.size_diff > 0:
            print(f"    {s.size_diff/1024:.0f} KB  {s.traceback[0]}")

    print("\n  Top 10 allocations (_downsample_stack x2):")
    stats = snap2.compare_to(snap1, "lineno")
    for s in stats[:10]:
        if s.size_diff > 0:
            print(f"    {s.size_diff/1024:.0f} KB  {s.traceback[0]}")

    del ds1, ds2, ms, sparse
    gc.collect()


if __name__ == "__main__":
    os.chdir(EXAMPLE_DIR)

    print("Memory profiling: MultiScaleImageStack construction")
    print(f"Initial RSS: {rss_mb():.0f} MB")

    # Progressive scan: 10, 30, 60, 120, 180 omega frames
    for n_omega in [10, 30, 60, 120, 180]:
        run_profile(n_omega)

    # Python-level breakdown at n_omega=30
    run_tracemalloc_profile(n_omega=30)
