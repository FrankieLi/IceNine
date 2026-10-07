#!/usr/bin/env python3
"""Profilers (Task 3): cProfile of 5 representative runs per pipeline, torch.profiler (CPU) of the
network stage with record_function around its sub-stages, and the forward pass on MPS vs CPU.

  uv run python scripts/profiling/prof_profile.py cprofile-u0       # reconstruct pipelines + F1
  uv run python scripts/profiling/prof_profile.py cprofile-seeded   # H0 H1 H3 MCr HG
  uv run python scripts/profiling/prof_profile.py torch-profiler    # net sub-stages
  uv run python scripts/profiling/prof_profile.py mps               # forward pass MPS vs CPU

cProfile adds a per-call overhead, so its shares overstate Python-level time next to the stage
timer's; both are reported. Files are written to benchmarks/profiling/ with directories stripped
from the function names (no absolute paths).
"""

import argparse
import contextlib
import cProfile
import io
import json
import pstats
import sys
import time
from pathlib import Path
from typing import Any, Dict, List, Tuple

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import prof_common as PC  # noqa: E402

OUT = PC.OUT
TOP = 30
N_RUNS = 5


# ---------------------------------------------------------------------------
# cProfile helpers
# ---------------------------------------------------------------------------


def _category(key: Tuple[str, int, str]) -> str:
    fn, _line, name = key
    if fn == "~" or fn.startswith("<") or fn.startswith("{"):
        return "C builtins and methods (numpy / torch / scipy kernels, len, ...)"
    if "/icenine/" in fn:
        return "icenine/ python code"
    if "/site-packages/" in fn:
        return "third-party python code (numpy / torch / scipy / pymatgen wrappers)"
    if "/scripts/" in fn or "/benchmarks/" in fn:
        return "scripts (harness and study code)"
    return "other python (stdlib)"


def shares(prof: cProfile.Profile) -> Dict[str, float]:
    st = pstats.Stats(prof)
    tot = sum(v[2] for v in st.stats.values())  # type: ignore[attr-defined]
    out: Dict[str, float] = {}
    for key, v in st.stats.items():  # type: ignore[attr-defined]
        c = _category(key)
        out[c] = out.get(c, 0.0) + v[2]
    return {k: v / tot for k, v in sorted(out.items(), key=lambda kv: -kv[1])} | {
        "_total_tottime_s": tot
    }


def dump(prof: cProfile.Profile, name: str, header: str) -> Dict[str, Any]:
    OUT.mkdir(parents=True, exist_ok=True)
    s = io.StringIO()
    st = pstats.Stats(prof, stream=s)
    st.strip_dirs().sort_stats("cumulative").print_stats(TOP)
    (OUT / f"cprofile_{name}.txt").write_text(header + "\n" + s.getvalue())
    return shares(prof)


# ---------------------------------------------------------------------------
# cProfile: U0 reconstruct pipelines (+ F1 post hocs) and the seeded pipelines
# ---------------------------------------------------------------------------


def cprofile_u0() -> None:
    import prof_u0 as U

    C = U.C
    U.init_worker(C.worker_args())
    from optimizer_sweep import voxel_context

    cases = []
    for vpos, v in enumerate(U.B.voxels()):
        for var in C.VARIANTS:
            e0 = C.CACHE_DIR / "e0" / f"v{v}_{var}.npz"
            if e0.exists() and "unbuildable" not in np.load(e0).files:
                cases.append((v, vpos, var))
    pick = cases[:: max(1, len(cases) // N_RUNS)][:N_RUNS]
    profs = {p: cProfile.Profile() for p in U.PIPES}
    W = C.get_worker()
    walls: Dict[str, List[float]] = {p: [] for p in U.PIPES}
    for v, vpos, var in pick:
        keys = np.load(C.CACHE_DIR / "images" / f"v{v}_{var}.npz")["keys"]
        C.attach(keys)
        vctx = voxel_context(W.ctx, v)
        e0 = np.load(C.CACHE_DIR / "e0" / f"v{v}_{var}.npz")
        phase = vctx.voxel.phase
        ctxd = dict(
            W=W,
            vctx=vctx,
            vpos=vpos,
            keys=keys,
            fold=int(U.M.fold_of(np.array([vpos]))[0]),
            fe_full=U._fe(None, phase),
            fe_low=U._fe(U.PROXY_Q, phase),
            diameter=float(e0["s0_L3_disc_diameter"]),
        )
        U._models(ctxd["fold"])
        for pipe in U.RECON:
            profs[pipe].enable()
            r = U.run_recon(pipe, ctxd)
            profs[pipe].disable()
            walls[pipe].append(r["wall"])
            if pipe in ("baseline", "proxy"):
                name = "F1" if pipe == "baseline" else "proxy+F1"
                profs[name].enable()
                f = U.run_f1(r, ctxd)
                profs[name].disable()
                walls[name].append(f["wall"])
        print("profiled", v, var, flush=True)
    res = {}
    for p in U.PIPES:
        hdr = (
            f"cProfile of {p}: {len(pick)} runs (voxel, variant): "
            f"{[(v, var) for v, _, var in pick]}; wall under the profiler, s: "
            f"{[round(x, 2) for x in walls[p]]}; top {TOP} by cumulative time"
        )
        res[p] = dump(profs[p], f"u0_{p.replace('+', '_')}", hdr)
        res[p]["wall_under_profiler_s"] = walls[p]
    PC.write_json(OUT / "cprofile_u0_shares.json", res)


def cprofile_seeded() -> None:
    import prof_seeded as SD

    NH, fs = SD.NH, SD.fs
    wargs, sweep = NH.worker_args([SD.MODEL])
    SD.init_worker(wargs)
    voxels = [int(v) for v in sweep["voxel_indices"]]
    ris_all = [1, 3, 6, 7, 8]
    ris_hg = [6, 7, 8, 9, 6]
    profs = {p: cProfile.Profile() for p in SD.PIPES}
    walls: Dict[str, List[float]] = {p: [] for p in SD.PIPES}
    cases = []
    for k in range(N_RUNS):
        for p in SD.PIPES:
            ri = ris_hg[k] if p == "HG" else ris_all[k]
            cases.append((p, ri, k % 2))
    W = fs._W
    for pipe, ri, vi in cases:
        vpos = 0
        batches = list(
            NH.variant_batches(
                W.ctx,
                voxels[vpos],
                vpos,
                ri,
                sweep["n_roi"][vpos, ri],
                sweep["fail_pass1"][vpos, ri],
                [NH.VARIANTS[vi]],
            )
        )
        vb = batches[0]
        keys = NH.case_keys(W.ctx, vb, 0)
        fs.attach_images(keys)
        cs = dict(
            vctx=vb.vctx,
            vb=vb,
            j=0,
            ri=ri,
            vidx=voxels[vpos],
            keys=keys,
            draws=vb.draws,
            seed_b=NH.b_seed(W.ctx.args, vpos, ri, 0, vb.vi),
            render1_s=vb.render1_s,
            n_batch=max(int(vb.have.sum()), 1),
        )
        profs[pipe].enable()
        r = SD.run_pipeline(pipe, cs)
        profs[pipe].disable()
        walls[pipe].append(r["wall"])
        print("profiled", pipe, ri, flush=True)
    res = {}
    for p in SD.PIPES:
        used = [ri for q, ri, _ in cases if q == p]
        hdr = (
            f"cProfile of {p}: {N_RUNS} cases, radius indices {used} (variants alternate clean / "
            f"realistic, voxel 0, direction 0); wall under the profiler, s: "
            f"{[round(x, 2) for x in walls[p]]}; top {TOP} by cumulative time"
        )
        res[p] = dump(profs[p], f"seeded_{p}", hdr)
        res[p]["wall_under_profiler_s"] = walls[p]
    PC.write_json(OUT / "cprofile_seeded_shares.json", res)


# ---------------------------------------------------------------------------
# torch.profiler (CPU) with record_function around the net sub-stages
# ---------------------------------------------------------------------------


class RecordTimer(PC.ST.StageTimer):
    """A StageTimer whose stages are also torch.profiler record_function ranges."""

    @contextlib.contextmanager
    def stage(self, name: str) -> Any:  # type: ignore[override]
        from torch.profiler import record_function

        with record_function(f"net::{name}"):
            with super().stage(name):
                yield


def torch_profiler() -> None:
    import torch
    from torch.profiler import ProfilerActivity, profile

    import prof_seeded as SD

    NH, fs = SD.NH, SD.fs
    wargs, sweep = NH.worker_args([SD.MODEL])
    SD.init_worker(wargs)
    W = fs._W
    T = RecordTimer()
    I = PC.Instrument(T, W.rec.params.max_mc_steps, net=True)  # noqa: E741
    I.install()
    I.mark(W.local_fn, "local")
    SD._S.update(T=T, I=I)
    voxels = [int(v) for v in sweep["voxel_indices"]]
    ri = 7  # 2 deg
    cases = []
    for vi in (0, 1):
        vb = list(
            NH.variant_batches(
                W.ctx,
                voxels[0],
                0,
                ri,
                sweep["n_roi"][0, ri],
                sweep["fail_pass1"][0, ri],
                [NH.VARIANTS[vi]],
            )
        )[0]
        for j in range(3):
            keys = NH.case_keys(W.ctx, vb, j)
            cases.append((vb, j, keys))
    snaps_before = T.snapshot()
    with profile(activities=[ProfilerActivity.CPU], record_shapes=False) as prof:
        for vb, j, keys in cases[: N_RUNS + 1]:
            fs.attach_images(keys)
            cs = dict(
                vctx=vb.vctx,
                vb=vb,
                j=j,
                ri=ri,
                vidx=voxels[0],
                keys=keys,
                draws=vb.draws,
                seed_b=NH.b_seed(W.ctx.args, 0, ri, j, vb.vi),
                render1_s=vb.render1_s,
                n_batch=max(int(vb.have.sum()), 1),
            )
            # the network stage only (H3's passes), not the finisher
            net = SD._S["net"]
            ctx = W.ctx
            with T.stage("net"):
                p1, _ = SD.ps.prepare_nominal(ctx, voxels[0], vb.R_nom0[j])
                NH.net_stage(
                    ctx,
                    net,
                    voxels[0],
                    [p1],
                    SD._rows(vb.b1, j),
                    np.array([p1 is not None and bool(vb.ok1[j])]),
                    vb.delta0[j : j + 1],
                    vb.R_nom0[j][None],
                    vb.vctx.R_true,
                    [vb.draws[j]],
                    vb.vctx.sources,
                    vb.variant,
                    vb.seed,
                    timer=T,
                )
    ev = prof.key_averages()
    lines = [
        f"torch.profiler (CPU), net x3 stage of {len(cases[:N_RUNS + 1])} cases (radius 2 deg, "
        f"voxel 0, clean and realistic), threads={torch.get_num_threads()}",
        "record_function ranges (net::<stage>): total self+children CPU time in ms and calls",
    ]
    rf = [e for e in ev if e.key.startswith("net::")]
    for e in sorted(rf, key=lambda e: -e.cpu_time_total):
        lines.append(
            f"  {e.key:28s} total {e.cpu_time_total / 1e3:9.2f} ms  calls {e.count:4d}  "
            f"per call {e.cpu_time_total / 1e3 / max(e.count, 1):8.3f} ms"
        )
    lines.append("")
    lines.append("top operators by self CPU time:")
    lines.append(ev.table(sort_by="self_cpu_time_total", row_limit=25, max_name_column_width=60))
    (OUT / "torch_profiler_net.txt").write_text("\n".join(lines))
    excl = {k: round(v, 4) for k, v in T.delta_since(snaps_before, "exclusive").items()}
    PC.write_json(OUT / "torch_profiler_net_stage_exclusive.json", excl)
    print("\n".join(lines[:12]))


# ---------------------------------------------------------------------------
# Forward pass: MPS vs CPU
# ---------------------------------------------------------------------------


def mps_vs_cpu(reps: int = 200) -> None:
    import torch

    import prof_seeded as SD

    NH, fs, ps = SD.NH, SD.fs, SD.ps
    from icenine.orientation_eval import decode_windows

    wargs, sweep = NH.worker_args([SD.MODEL])
    SD.init_worker(wargs)
    W = fs._W
    ctx = W.ctx
    voxels = [int(v) for v in sweep["voxel_indices"]]
    vb = list(
        NH.variant_batches(
            ctx, voxels[0], 0, 7, sweep["n_roi"][0, 7], sweep["fail_pass1"][0, 7], [NH.VARIANTS[1]]
        )
    )[0]
    rows = np.nonzero(vb.ok1)[0]
    K = ctx.args.frame_half_width
    net = SD._S["net"]
    with torch.no_grad():
        x_all = decode_windows(vb.b1["windows"][rows], K)
    ctx_all, nom_all = vb.b1["context"][rows], vb.b1["nom_off"][rows]
    res: Dict[str, Any] = dict(
        torch=torch.__version__,
        mps_available=bool(torch.backends.mps.is_available()),
        n_rows=int(len(rows)),
        x_shape=list(x_all.shape),
        reps=reps,
    )

    def bench(fn: Any, sync: Any) -> Dict[str, float]:
        for _ in range(20):
            fn()
        sync()
        t = []
        for _ in range(reps):
            t0 = time.perf_counter()
            fn()
            sync()
            t.append(time.perf_counter() - t0)
        return PC.summarize_times([x * 1e3 for x in t])  # ms

    nop = lambda: None  # noqa: E731
    for bs in (1, 20):
        n = min(bs, len(rows))
        x, c, o = x_all[:n], ctx_all[:n], nom_all[:n]
        out: Dict[str, Any] = {"batch": n}
        for threads in (1, 8):
            torch.set_num_threads(threads)
            with torch.no_grad():
                out[f"cpu_{threads}thread"] = bench(lambda: net(x, c, {"nom_off": o}), nop)
        torch.set_num_threads(1)
        if res["mps_available"]:
            dev = torch.device("mps")
            netm = SD.ps.load_model(dict(ctx.args.models)[SD.MODEL])[0].to(dev)
            xm, cm, om = x.to(dev), c.to(dev), o.to(dev)
            with torch.no_grad():
                out["mps_resident"] = bench(
                    lambda: netm(xm, cm, {"nom_off": om}), torch.mps.synchronize
                )
                out["mps_with_transfer"] = bench(
                    lambda: [t.cpu() for t in netm(x.to(dev), c.to(dev), {"nom_off": o.to(dev)})],
                    torch.mps.synchronize,
                )
                m_c, L_c = net(x, c, {"nom_off": o})
                m_m, L_m = netm(xm, cm, {"nom_off": om})
            out["max_abs_diff_mean_deg"] = float((m_c - m_m.cpu()).abs().max())
            out["max_abs_diff_chol"] = float((L_c - L_m.cpu()).abs().max())
        res[f"batch{bs}"] = out
    PC.write_json(OUT / "mps_vs_cpu.json", res)
    print(json.dumps(res, indent=1))


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    ap.add_argument("cmd", choices=["cprofile-u0", "cprofile-seeded", "torch-profiler", "mps"])
    a = ap.parse_args()
    PC.write_json(OUT / f"isolation_{a.cmd}.json", PC.isolation_record(a.cmd))
    {
        "cprofile-u0": cprofile_u0,
        "cprofile-seeded": cprofile_seeded,
        "torch-profiler": torch_profiler,
        "mps": mps_vs_cpu,
    }[a.cmd]()


if __name__ == "__main__":
    main()
