"""S3 Step 0 (go/no-go): what a batched quick MC could gain. Script only, no library change.

Two measurements (docs/todo_batched_quick_mc.md, Step 0):

``trace``   Instrumented sequential ``reconstruct_voxel`` runs. Per level: the number of quick-MC
            candidates (the batch size), which candidates' single block failed (no trial strictly
            below the start cost), the random draws each ``optimize`` call consumed (checked: 30 per
            candidate, 33 when its block failed), and the rounds option A (exact, speculative) would
            need for speculation depth k = 1, 2, 4 (computed from the sequential outcomes). Also the
            share of evaluations and wall time the quick MC has in the search.
``batched`` A throwaway batched local cost (stages A-C in torch for B orientations, host hand-off,
            stage D in C once per orientation) on CPU and on MPS, timed against
            ``VoxelCostFunction.evaluate`` on the same cost function, with a per-orientation
            equality check against ``evaluate`` and a count of the calls that would force a host
            sync on a GPU.

Usage (from icenine_py/; each --set in its own process, they chdir into different examples):
  uv run python scripts/speed_parity/s3_step0.py trace   --set three --out DIR [--seeds 3]
  uv run python scripts/speed_parity/s3_step0.py trace   --set e0    --out DIR [--n-voxels 10]
  uv run python scripts/speed_parity/s3_step0.py batched --set three --out DIR
  uv run python scripts/speed_parity/s3_step0.py batched --set e0    --out DIR
  uv run python scripts/speed_parity/s3_step0.py summary --out DIR

Timings call ``preflight.require_quiet`` first (``--allow-busy`` records the preflight and
continues, for a machine that is not quiet; the JSON then says so).
"""

import argparse
import contextlib
import copy
import io
import json
import math
import sys
import time
from collections import defaultdict
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
import torch

HERE = Path(__file__).resolve().parent
ICENINE_PY = HERE.parents[1]
sys.path.insert(0, str(ICENINE_PY))
sys.path.insert(0, str(ICENINE_PY / "scripts" / "common"))
sys.path.insert(0, str(HERE))
torch.set_num_threads(1)  # single-worker timings

import preflight  # noqa: E402
import stats  # noqa: E402

from icenine import cost_functions as CF  # noqa: E402
from icenine.cost_functions import OverlapInfo, VoxelCostFunction  # noqa: E402
from icenine.diffraction_core import get_scattering_omegas_torch  # noqa: E402
from icenine.orientation_search import MCOptimizer  # noqa: E402

KS = (1, 2, 4)  # speculation depths
BATCH_SIZES = (1, 4, 16, 64, 256, 1024)
N_TRIAL_STEPS = 10  # max_mc_steps of the quick MC
N_START_EVALS = 2  # start orientation as a matrix, and again after the quaternion round trip


# ---------------------------------------------------------------------------------------------
# set-up of the two example sets
# ---------------------------------------------------------------------------------------------


class ThreeVoxels:
    """Example2.ThreeVoxels voxel 0 (images of the three-grain simulation), two configs."""

    name = "three"

    def __init__(self) -> None:
        import prof_eval as P  # chdirs into the example and builds the physics

        self.P = P

    def reconstructor(self, cfgname: str, cpp_rule: bool = False) -> Any:
        from icenine.config_file import ConfigFile
        from icenine.reconstructor import AdaptiveVoxelReconstructor, setup_reconstruction
        from icenine.sampling import load_fundamental_zone_file

        P = self.P
        cfg = ConfigFile.from_file(str(P.ex / "ConfigFiles" / cfgname))
        fz = load_fundamental_zone_file(cfg.fundamental_zone_filename)
        setup = setup_reconstruction(cfg, exp_data=P.exp_data, fz_orientations=fz)
        if cpp_rule:
            setup.search_params.cpp_overlap_rule = True
        return AdaptiveVoxelReconstructor(setup)

    def case(self) -> Tuple[Any, np.ndarray, int]:
        P = self.P
        return P.verts, np.asarray(P.vox.orientation, dtype=np.float64), P.vox.phase


class E0Set:
    """Voxels of the E0 clean set (FindOptimal robustness study), ReconstructQ8.config physics."""

    name = "e0"

    def __init__(self) -> None:
        sys.path.insert(0, str(ICENINE_PY / "scripts" / "findoptimal_robustness"))
        import common as C

        self.C = C
        C.init_worker(C.worker_args())
        self.info = dict(np.load(C.OUT_DIR / "voxels.npz"))
        self._case: Optional[Tuple[Any, ...]] = None

    def load(self, vpos: int) -> Tuple[Any, np.ndarray, int, Any]:
        """Build the clean case of voxel position ``vpos``: (vertices, R_true, phase, rec)."""
        C = self.C
        vidx = int(self.info["voxel_indices"][vpos])
        case = C.build_case(vidx, vpos, "clean")
        if case is None:
            raise RuntimeError(f"voxel {vidx} (position {vpos}) is unbuildable")
        keys, R_true, vctx = case
        C.attach(keys)
        W = C.get_worker()
        self._case = (vctx.vertices, np.asarray(R_true, dtype=np.float64), vctx.voxel.phase)
        return (*self._case, W.rec)

    def run_seed(self, vpos: int, seed: int) -> np.random.Generator:
        return self.C.run_seed(vpos, seed)

    def quiet(self) -> Any:
        return self.C.quiet()


# ---------------------------------------------------------------------------------------------
# trace: instrumented sequential reconstruct_voxel
# ---------------------------------------------------------------------------------------------


def draws_consumed(st0: Dict[str, Any], st1: Dict[str, Any]) -> Optional[int]:
    """Number of doubles drawn between two PCG64 states (searches 0..40), or None."""
    g = np.random.Generator(np.random.PCG64())
    for n in range(0, 41):
        g.bit_generator.state = copy.deepcopy(st0)
        if n:
            g.random(n)
        if g.bit_generator.state == st1:
            return n
    return None


def trace_reconstruct(
    rec: Any,
    vertices: Any,
    phase: int,
    rng: np.random.Generator,
    quiet: Any = None,
    sim_ks: Tuple[int, ...] = (),
) -> Dict[str, Any]:
    """One reconstruct_voxel with the quick-MC calls recorded (the result is unchanged)."""
    state: Dict[str, Any] = dict(level=None, phase="other")
    # uninstrumented run first (same generator state): wall time and the reference result
    rng_plain = copy.deepcopy(rng)
    t0 = time.perf_counter()
    with contextlib.redirect_stdout(io.StringIO()):
        ref = rec.reconstruct_voxel(vertices, phase, rng=rng_plain)
    wall_plain = time.perf_counter() - t0
    pending: List[Tuple[List[Optional[int]], int, Dict[str, Any], Dict[str, Any]]] = []
    captured: Dict[int, Dict[str, Any]] = {}  # level -> state at level start, one entry per call
    levels: Dict[int, Dict[str, Any]] = {}
    ev: Dict[str, List[float]] = defaultdict(lambda: [0, 0.0])
    opt_t: Dict[str, float] = defaultdict(float)
    orig_opt, orig_eval = MCOptimizer.optimize, VoxelCostFunction.evaluate

    def eval_wrap(self: Any, *a: Any, **k: Any) -> Any:
        t0 = time.perf_counter()
        r = orig_eval(self, *a, **k)
        e = ev[state["phase"]]
        e[0] += 1
        e[1] += time.perf_counter() - t0
        return r

    def opt_wrap(self: Any, *a: Any, **k: Any) -> Any:
        quick = state["level"] is not None and k.get("max_mc_steps") == 10
        prev = state["phase"]
        state["phase"] = "quick_mc" if quick else "refine"
        st0 = copy.deepcopy(self.rng.bit_generator.state) if quick else None
        t0 = time.perf_counter()
        res = orig_opt(self, *a, **k)
        opt_t[state["phase"]] += time.perf_counter() - t0
        if quick:
            lr = self.last_run
            assert lr["n_blocks"] == 1 and lr["steps_run"] == N_TRIAL_STEPS, lr
            lv = levels[state["level"]]
            lv["failed"].append(bool(lr["n_restarts"] == 1))
            lv["draws"].append(None)
            pending.append(
                (
                    lv["draws"],
                    len(lv["draws"]) - 1,
                    st0,
                    copy.deepcopy(self.rng.bit_generator.state),
                )
            )
            lv["start_cost"].append(float(lr["cost_start"]))
            lv["end_cost"].append(float(lr["cost_end"]))
            cap = captured.setdefault(state["level"], dict(st_level0=st0, calls=[]))
            cap["calls"].append(
                dict(
                    R0=np.array(k["initial_orientation"], copy=True),
                    step=float(k["angular_step"]),
                    R=np.array(res.orientation, copy=True),
                    cost=float(res.cost),
                    opt=self,
                )
            )
        state["phase"] = prev
        return res

    def recorder(name: str, d: Dict[str, Any]) -> None:
        if name == "discrete":
            state["level"] = int(d["level"])
            levels[state["level"]] = dict(
                level=state["level"],
                n_q_max=float(d["n_q_max"]),
                diameter_deg=math.degrees(float(d["diameter"])),
                failed=[],
                draws=[],
                start_cost=[],
                end_cost=[],
            )
        elif name == "quick_mc":
            state["level"] = None

    MCOptimizer.optimize = opt_wrap  # type: ignore[method-assign]
    VoxelCostFunction.evaluate = eval_wrap  # type: ignore[method-assign]
    rec.recorder = recorder
    try:
        t0 = time.perf_counter()
        with quiet or contextlib.redirect_stdout(io.StringIO()):
            res = rec.reconstruct_voxel(vertices, phase, rng=rng)
        wall = time.perf_counter() - t0
    finally:
        MCOptimizer.optimize = orig_opt  # type: ignore[method-assign]
        VoxelCostFunction.evaluate = orig_eval  # type: ignore[method-assign]
        rec.recorder = None
    for draws, i, st0_, st1_ in pending:
        draws[i] = draws_consumed(st0_, st1_)
    if sim_ks:
        for lv_no, cap in captured.items():
            levels[lv_no]["sim"] = {
                str(kk): simulate_option_a(cap["calls"], cap["st_level0"], kk) for kk in sim_ks
            }
    g, loc, tot = rec.last_eval_counts
    same = bool(np.array_equal(res.orientation, ref.orientation) and res.cost == ref.cost)
    return dict(
        wall_s=wall,
        wall_plain_s=wall_plain,
        instrumented_equals_plain=same,
        cost=float(res.cost),
        eval_counts=dict(global_=g, local=loc, total=tot),
        evaluate_calls={k: dict(n=int(v[0]), seconds=v[1]) for k, v in ev.items()},
        optimize_seconds=dict(opt_t),
        levels=[levels[k] for k in sorted(levels)],
    )


def option_a_rounds(failed: List[bool], k: int) -> Dict[str, Any]:
    """Rounds of option A for a level whose candidates have the given sequential block outcomes.

    Round r starts at the first unresolved candidate p. Every unresolved candidate j >= p runs k
    chains (assumed offsets F_known + 0..k-1, F_known = failures among 0..p-1). Walking in order,
    candidate j is resolved while the failures among p..j-1 are at most k-1 (its true offset is
    then one of the offsets it ran); the first candidate with k failures before it, and all after
    it, go to the next round."""
    n = len(failed)
    p, rounds = 0, 0
    batches: List[int] = []
    while p < n:
        rounds += 1
        batches.append((n - p) * k)
        extra, j = 0, p
        while j < n and extra <= k - 1:
            extra += int(failed[j])
            j += 1
        p = j
    chains = sum(batches)
    return dict(
        rounds=rounds,
        chains=chains,
        batches=batches,
        lockstep_steps=rounds * N_TRIAL_STEPS + 1,  # +1 for the start evaluations
        computed_evals=N_START_EVALS * n + N_TRIAL_STEPS * chains,
        counted_evals=(N_START_EVALS + N_TRIAL_STEPS) * n,
    )


def simulate_option_a(
    calls: List[Dict[str, Any]], st_level0: Dict[str, Any], k: int
) -> Dict[str, Any]:
    """Option A (exact) run for real with ``evaluate`` as the cost: the speculative rounds with k
    chains per unresolved candidate, on the generator state at the start of the level. Reports the
    rounds, the chain-trial evaluations of the plan as written (10 per chain), the distinct
    evaluations when a chain's trial evaluation is reused (memoised) across offsets, and whether
    every resolved chain equals the sequential result (orientation, cost). Logic and counting only,
    no batching. Memoisation key: (candidate, triple index m, accepted triples so far); sound
    because trial t of a chain at offset off uses the triple m = 10 j + off + t of the level's
    stream and its proposal depends only on the candidate's start and the accepted triples."""
    from icenine.sampling import QuaternionGrid, _quat_multiply, matrix_to_quaternion
    from icenine.sampling import quaternion_to_matrix

    n = len(calls)
    opt = calls[0]["opt"]
    cost_fn, verts, phase = opt.cost_fn, opt.voxel_vertices, opt.phase_index
    step = calls[0]["step"]
    radius = math.tan(step) / math.sqrt(12.0) if step > 0 else 0.01
    g = np.random.Generator(np.random.PCG64())
    g.bit_generator.state = copy.deepcopy(st_level0)
    # scalar draws, as today; the buffer is 3 k + 30 draws longer than the 33 n the true chains
    # can use, for the speculative chains of the last candidates (those chains are discarded)
    n_buf = 33 * n + 3 * k + 3 * N_TRIAL_STEPS
    u = np.array([g.uniform(-radius, radius) for _ in range(n_buf)])
    grid = QuaternionGrid()
    starts = []
    for c in calls:  # the two start evaluations of optimize (matrix, and after the quaternion trip)
        iq = matrix_to_quaternion(c["R0"])
        gm = cost_fn.evaluate(c["R0"], verts, phase).cost
        oc = cost_fn.evaluate(quaternion_to_matrix(iq), verts, phase).cost
        starts.append((iq, gm, oc))
    memo: Dict[Tuple[Any, ...], float] = {}
    n_chains = 0

    def chain(j: int, off: int) -> Tuple[bool, np.ndarray, float]:
        nonlocal n_chains
        n_chains += 1
        iq, gm, oc = starts[j]
        opt_q, opt_cost, acc = iq.copy(), oc, ()
        for t in range(N_TRIAL_STEPS):
            m = N_TRIAL_STEPS * j + off + t
            delta = grid.get_near_identity_point(u[3 * m], u[3 * m + 1], u[3 * m + 2])
            trial_q = _quat_multiply(delta, opt_q)
            key = (j, m, acc)
            cost = memo.get(key)
            if cost is None:
                cost = memo[key] = cost_fn.evaluate(
                    quaternion_to_matrix(trial_q), verts, phase
                ).cost
            if cost < opt_cost:
                opt_q, opt_cost, acc = trial_q, cost, acc + (m,)
        failed = opt_cost >= gm
        return failed, quaternion_to_matrix(iq if failed else opt_q), (gm if failed else opt_cost)

    p, f_known, rounds, ok = 0, 0, 0, True
    new_per_round: List[int] = []
    while p < n:
        rounds += 1
        before = len(memo)
        got = {(j, f_known + e): chain(j, f_known + e) for j in range(p, n) for e in range(k)}
        new_per_round.append(len(memo) - before)
        off, j = f_known, p
        while j < n and off - f_known <= k - 1:
            failed, R, cost = got[(j, off)]
            ok &= bool(np.array_equal(R, calls[j]["R"]) and cost == calls[j]["cost"])
            off += int(failed)
            j += 1
        f_known, p = off, j
    return dict(
        k=k,
        rounds=rounds,
        chains=n_chains,
        plan_trial_evals=N_TRIAL_STEPS * n_chains,
        memo_trial_evals=len(memo),
        start_evals=N_START_EVALS * n,
        counted_evals=(N_START_EVALS + N_TRIAL_STEPS) * n,
        new_trial_evals_per_round=new_per_round,
        matches_sequential=ok,
    )


def cmd_trace(a: argparse.Namespace) -> None:
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    runs: List[Dict[str, Any]] = []
    sim_ks = tuple(int(x) for x in a.sim_ks.split(",") if x)
    if a.set == "three":
        S = ThreeVoxels()
        verts, _, phase = S.case()
        arms = [
            ("Example2.Simulation.config", False),
            ("ReconstructQ8.config", False),
            ("ReconstructQ8.config", True),
        ]
        for ai, (cfgname, cpp) in enumerate(arms):
            if a.arms and ai not in [int(x) for x in a.arms.split(",")]:
                continue
            rec = S.reconstructor(cfgname, cpp)
            for seed in range(a.seeds):
                r = trace_reconstruct(rec, verts, phase, np.random.default_rng(seed), sim_ks=sim_ks)
                r.update(
                    set="three",
                    config=cfgname,
                    cpp_overlap_rule=cpp,
                    voxel=0,
                    seed=seed,
                    levels_config=dict(
                        max_local_resolution=rec.params.max_local_resolution,
                        mc_radius_scale_factor=rec.params.mc_radius_scale_factor,
                    ),
                )
                runs.append(r)
                print(_run_line(r), flush=True)
    else:
        S = E0Set()
        for vpos in range(a.n_voxels):
            verts, _, phase, rec = S.load(vpos)
            for seed in range(a.seeds):
                r = trace_reconstruct(rec, verts, phase, S.run_seed(vpos, seed), sim_ks=sim_ks)
                r.update(
                    set="e0",
                    config="ReconstructQ8.config (E0 clean)",
                    cpp_overlap_rule=False,
                    voxel=vpos,
                    seed=seed,
                    levels_config=dict(
                        max_local_resolution=rec.params.max_local_resolution,
                        mc_radius_scale_factor=rec.params.mc_radius_scale_factor,
                    ),
                )
                runs.append(r)
                print(_run_line(r), flush=True)
    (out / f"trace_{a.set}{a.tag}.json").write_text(json.dumps(runs, indent=1, default=float))


def _run_line(r: Dict[str, Any]) -> str:
    nc = [len(lv["failed"]) for lv in r["levels"]]
    nf = [int(sum(lv["failed"])) for lv in r["levels"]]
    q = r["evaluate_calls"].get("quick_mc", dict(n=0, seconds=0.0))
    return (
        f"{r['set']} {r['config']} cpp={r['cpp_overlap_rule']} voxel {r['voxel']} seed {r['seed']}:"
        f" wall {r['wall_s']:.1f}s, cand/level {nc}, failed {nf},"
        f" quick-MC evaluate calls {q['n']} ({q['seconds']:.2f}s)"
    )


# ---------------------------------------------------------------------------------------------
# batched: throwaway batched local cost
# ---------------------------------------------------------------------------------------------


class SyncCounter(torch.overrides.TorchFunctionMode):
    """Counts the tensor calls that would force a host sync on a GPU: device-to-host copies,
    scalar conversions and data-dependent output shapes (nonzero, boolean-mask indexing)."""

    NAMES = {
        "item",
        "cpu",
        "numpy",
        "tolist",
        "__bool__",
        "__int__",
        "__float__",
        "nonzero",
        "masked_select",
        "unique",
        "argwhere",
    }

    def __init__(self, assume_gpu: bool) -> None:
        """assume_gpu: count the calls on CPU tensors as well (what a GPU run would pay);
        else only calls on tensors that are on a non-CPU device."""
        super().__init__()
        self.assume_gpu = assume_gpu
        self.counts: Dict[str, int] = defaultdict(int)

    def __torch_function__(self, func: Any, types: Any, args: Any = (), kwargs: Any = None) -> Any:
        name = getattr(func, "__name__", str(func))
        hit = name in self.NAMES
        t0 = args[0] if args and isinstance(args[0], torch.Tensor) else None
        if hit and not self.assume_gpu and (t0 is None or t0.device.type == "cpu"):
            hit = False
        if hit and self.assume_gpu and name == "cpu":
            hit = False  # a no-op on CPU tensors; the .numpy() that follows is the counted sync
        if name == "__getitem__" and len(args) > 1:
            idx = args[1] if isinstance(args[1], tuple) else (args[1],)
            hit = any(isinstance(i, torch.Tensor) and i.dtype == torch.bool for i in idx)
            hit = hit and (self.assume_gpu or (t0 is not None and t0.device.type != "cpu"))
            name = "__getitem__(bool mask)"
        if name == "__setitem__" and len(args) > 1:
            idx = args[1] if isinstance(args[1], tuple) else (args[1],)
            hit = any(isinstance(i, torch.Tensor) and i.dtype == torch.bool for i in idx)
            hit = hit and (self.assume_gpu or (t0 is not None and t0.device.type != "cpu"))
            name = "__setitem__(bool mask)"
        if hit:
            self.counts[name] += 1
        return func(*args, **(kwargs or {}))


class ThrowawayBatchedCost:
    """Batched local cost of one VoxelCostFunction (pixel_radius 0, hard mode, no ray rule).

    ``variant='fixed'``: stages A-C on all 2K omega rows of every orientation (validity as a mask,
    no data-dependent shape until the single host transfer); the frame lookup and the compaction
    are done on the host from the transferred omegas. ``variant='compact'`` (CPU only): rows of
    the observable reflections only, boolean-mask compaction after the eta filter and after the
    frame lookup, as ``ScreenEvaluator`` does. Stage D is the existing C call, once per
    orientation, against a prebuilt table of all frames."""

    def __init__(self, cf: VoxelCostFunction, phase: int, device: str, variant: str) -> None:
        assert cf.pixel_radius == 0 and cf.ray_rule is None and cf.mode == "hard"
        assert not cf.cpp_edge_rules and cf.eta_limit is not None
        assert variant in ("fixed", "compact") and (variant == "fixed" or device == "cpu")
        self.cf, self.phase, self.variant = cf, phase, variant
        self.dev = torch.device(device)
        dev = self.dev
        f32 = dict(dtype=torch.float32, device=dev)
        g_hkl, g_mag = cf._phase_recip_vecs[phase]
        self.K = g_hkl.shape[0]
        self.g_hkl_view = g_hkl.T  # the non-contiguous view evaluate multiplies with
        self.g_hkl_T = g_hkl.T.contiguous().to(dev)
        self.g_mag = g_mag.to(dev)
        sim = cf.simulator
        self.energy, self.chi = sim.beam_energy, sim.beam_deflection_chi
        self.beam = sim.beam_direction.to(**f32)
        m = cf.sample.sample_to_lab_matrix
        self.base_rot = m[:3, :3].to(**f32)
        self.base_t = m[:3, 3].to(**f32)
        self.dets = []
        for d in cf.detector_list:
            pl = d._detector_plane
            self.dets.append(
                dict(
                    n=pl.normal.to(**f32),
                    d=torch.as_tensor(pl.d, **f32),
                    pos=d._position.to(**f32),
                    org=d._lab_frame_coord_origin.to(**f32),
                    bj=d._lab_frame_basis_j.to(**f32),
                    bk=d._lab_frame_basis_k.to(**f32),
                    hw=float(d.pixel_half_width),
                    pw=float(d.pixel_width),
                    hh=float(d.pixel_half_height),
                    ph=float(d.pixel_height),
                )
            )
        self.n_det = len(self.dets)
        rng = cf.range_map
        self.low, self.width = rng.low, rng.width
        self.index = np.array([-1 if i is None else int(i) for i in rng.index_list], dtype=np.int64)
        self.n_bins = len(self.index)
        n_w = int(self.index.max()) + 1
        # prebuilt table of all frames (a stand-in for S1's image table)
        self.table = [
            cf.exp_data.get_image(w, d).get_binary_numpy()
            for w in range(n_w)
            for d in range(self.n_det)
        ]
        first = cf.exp_data.get_image(int(self.index.max()), 0)
        self.rows, self.cols = first.num_rows, first.num_cols
        self.rule_bits = (1 if cf.cpp_overlap_rule else 0) | (2 if cf.cpp_edge_rules else 0)
        self.vert4 = None  # set per call (voxel vertices)

    # --- device stages ----------------------------------------------------------------------

    def _bragg(self, R: torch.Tensor) -> Tuple[torch.Tensor, ...]:
        """g = R g_hkl^T for every orientation, then the Bragg solve. On the CPU this reproduces
        ``evaluate`` bit for bit (found in Step 0): one 3x3 @ 3xK gemm per orientation against the
        transposed view of g_hkl (a batched matmul, or the same gemm on a contiguous copy of
        g_hkl^T, rounds 45% of the g rows differently), and g kept component-major (3, B K) so the
        elementwise kernels see contiguous columns as in ``evaluate``. On a GPU: one batched
        matmul (a GPU's transcendental functions differ from the CPU's in any case)."""
        B, K = R.shape[0], self.K
        if self.dev.type == "cpu":
            gT = torch.cat([R[b] @ self.g_hkl_view for b in range(B)], dim=1)  # (3, B K)
            g = gT.T
        else:
            g = torch.matmul(R, self.g_hkl_T).transpose(1, 2).reshape(B * K, 3)
        br = get_scattering_omegas_torch(g, self.g_mag.repeat(B), self.energy, self.chi)
        return g.contiguous(), br.observable, br.omega1, br.omega2

    def _rows_eta(self, omegas: torch.Tensor, normals: torch.Tensor) -> Tuple[torch.Tensor, ...]:
        """Rz(omega) @ base, lab normals, reflected beam and the eta mask; any leading shape."""
        cf = self.cf
        c, s = torch.cos(omegas), torch.sin(omegas)
        z, o = torch.zeros_like(c), torch.ones_like(c)
        Rz = torch.stack([c, -s, z, s, c, z, z, z, o], dim=-1).reshape(*c.shape, 3, 3)
        full_rot = Rz @ self.base_rot
        lab_n = torch.matmul(full_rot, normals.unsqueeze(-1)).squeeze(-1)
        dot = (self.beam * lab_n).sum(dim=-1, keepdim=True)
        refl = self.beam - 2.0 * dot * lab_n
        rd = torch.norm(refl, dim=-1)
        safe = torch.where(rd > 0, rd, torch.ones_like(rd))
        eta = torch.atan2(torch.abs(refl[..., 1]) / safe, torch.abs(refl[..., 2]) / safe)
        valid = (eta < cf.eta_limit) & (rd > 0)
        if cf.min_sin_eta > 0.0:
            valid = valid & (torch.sin(eta) >= cf.min_sin_eta)
        return full_rot, refl, valid

    def _stage_bc(
        self, full_rot: torch.Tensor, refl: torch.Tensor, verts: torch.Tensor
    ) -> Tuple[torch.Tensor, torch.Tensor]:
        """Stages B and C: lab vertices and detector pixel coordinates; any leading shape."""
        lead = full_rot.shape[:-2]
        f4 = torch.zeros(*lead, 4, 4, dtype=torch.float32, device=self.dev)
        f4[..., :3, :3] = full_rot
        f4[..., :3, 3] = self.base_t
        f4[..., 3, 3] = 1.0
        v4 = torch.cat([verts, torch.ones(3, 1, device=self.dev)], dim=1)  # (3, 4)
        lab = torch.einsum("...ij,vj->...vi", f4, v4)[..., :3]  # (..., 3, 3)
        dirs = refl.unsqueeze(-2).expand(*lead, 3, 3)
        hit_all, pix = [], []
        for d in self.dets:
            denom = (dirs * d["n"]).sum(dim=-1)
            numer = -((lab * d["n"]).sum(dim=-1) + d["d"])
            parallel = torch.abs(denom) < 1e-8
            safe_d = torch.where(parallel, torch.ones_like(denom), denom)
            t = torch.where(parallel, torch.zeros_like(denom), numer / safe_d)
            hits = (~parallel) & (t > 0)
            rel = lab + t.unsqueeze(-1) * dirs - d["pos"]
            loc = rel - d["org"]
            j = (loc * d["bj"]).sum(dim=-1)
            k = (loc * d["bk"]).sum(dim=-1)
            col = (j + d["hw"]) / d["pw"]
            row = (k + d["hh"]) / d["ph"]
            hit_all.append(hits.all(dim=-1))
            pix.append(torch.stack([col, row], dim=-1))
        return torch.stack(hit_all, dim=-1), torch.stack(
            pix, dim=-3
        )  # (..., n_det), (..., n_det, 3, 2)

    # --- the two variants -------------------------------------------------------------------

    def _forward_fixed(self, R: torch.Tensor, verts: torch.Tensor) -> torch.Tensor:
        B, K = R.shape[0], self.K
        g, obs, om1, om2 = self._bragg(R)
        normals = (g / torch.norm(g, dim=1, keepdim=True)).reshape(B, K, 3)
        omegas = torch.cat([om1.reshape(B, K), om2.reshape(B, K)], dim=1)  # (B, 2K)
        normals = torch.cat([normals, normals], dim=1)
        obs2 = torch.cat([obs.reshape(B, K)] * 2, dim=1)
        full_rot, refl, valid = self._rows_eta(omegas, normals)
        valid = valid & obs2
        det_hit, pix = self._stage_bc(full_rot, refl, verts)
        return torch.cat(
            [
                pix.reshape(B, 2 * K, -1),
                det_hit.float(),
                omegas.unsqueeze(-1),
                valid.float().unsqueeze(-1),
            ],
            dim=-1,
        )  # (B, 2K, 6 n_det + n_det + 2) float32: the single device-to-host transfer

    def evaluate_many(
        self, R_np: np.ndarray, verts: torch.Tensor, timing: Optional[Dict[str, float]] = None
    ) -> List[OverlapInfo]:
        B, K = len(R_np), self.K
        n_det = self.n_det
        t0 = time.perf_counter()
        R = torch.from_numpy(np.asarray(R_np, dtype=np.float64)).float().to(self.dev)
        verts = verts.to(self.dev)
        if self.variant == "fixed":
            packed = self._forward_fixed(R, verts)
            arr = packed.cpu().numpy()  # the one host sync of the call
            t1 = time.perf_counter()
            omega = arr[..., 6 * n_det + n_det]
            ok = arr[..., -1] > 0.5
            wedge = self._wedge(omega)
            keep = ok & (wedge >= 0)
            counts = keep.sum(axis=1)
            sel = arr[keep]
            pix = np.ascontiguousarray(sel[:, : 6 * n_det].reshape(-1, n_det, 3, 2))
            hit = np.ascontiguousarray(sel[:, 6 * n_det : 7 * n_det] > 0.5).astype(np.uint8)
            wedge_k = np.ascontiguousarray(wedge[keep].astype(np.int32))
        else:
            g, obs, om1, om2 = self._bragg(R)
            og = g[obs]
            normals = og / torch.norm(og, dim=1, keepdim=True)
            seg1 = torch.arange(B).repeat_interleave(K)[obs]
            omegas = torch.cat([om1[obs], om2[obs]])
            normals = torch.cat([normals, normals])
            seg = torch.cat([seg1, seg1])
            full_rot, refl, valid = self._rows_eta(omegas, normals)
            vidx = torch.nonzero(valid).squeeze(1)
            om_np = omegas[vidx].numpy()
            wedge_all = self._wedge(om_np)
            ok = wedge_all >= 0
            sel_i = vidx[torch.from_numpy(np.nonzero(ok)[0])]
            seg_np = seg[sel_i].numpy()
            order = np.argsort(seg_np, kind="stable")
            sel_i = sel_i[torch.from_numpy(order)]
            seg_np = seg_np[order]
            wedge_k = np.ascontiguousarray(wedge_all[ok][order].astype(np.int32))
            det_hit_t, pix_t = self._stage_bc(full_rot[sel_i], refl[sel_i], verts)
            t1 = time.perf_counter()
            pix = np.ascontiguousarray(pix_t.numpy())
            hit = np.ascontiguousarray(det_hit_t.numpy().astype(np.uint8))
            counts = np.bincount(seg_np, minlength=B)
        t2 = time.perf_counter()
        offs = np.concatenate([[0], np.cumsum(counts)])
        rb = (self.rule_bits,) if self.rule_bits else ()
        out: List[OverlapInfo] = []
        for b in range(B):
            s, e = int(offs[b]), int(offs[b + 1])
            if e == s:
                out.append(OverlapInfo())
                continue
            r = CF._c_stage_d_overlap(
                self.table,
                wedge_k[s:e],
                hit[s:e],
                None,
                pix[s:e],
                n_det,
                0,
                e - s,
                self.rows,
                self.cols,
                *rb,
            )
            out.append(
                OverlapInfo(
                    pixel_overlap=r[0],
                    pixel_on_detector=r[1],
                    peak_overlap=r[2],
                    peak_on_detector=r[3],
                    quality=r[4],
                    n_quality_points=r[5],
                    detectors_overlap=r[6],
                )
            )
        t3 = time.perf_counter()
        self.cf.eval_count += B
        if timing is not None:
            for k, v in (("torch", t1 - t0), ("host_prep", t2 - t1), ("stage_d", t3 - t2)):
                timing[k] = timing.get(k, 0.0) + v
        return out

    def _wedge(self, omega: np.ndarray) -> np.ndarray:
        """Frame index of each omega exactly as calculate_diffraction_overlap_batched (-1: none)."""
        bins = CF._omega_bins(omega, self.low, self.width, self.cf.cpp_edge_rules)
        okb = (bins >= 0) & (bins < self.n_bins)
        return np.where(okb, self.index[np.clip(bins, 0, self.n_bins - 1)], -1)


def info_tuple(oi: OverlapInfo) -> Tuple[Any, ...]:
    return (
        oi.pixel_overlap,
        oi.pixel_on_detector,
        oi.peak_overlap,
        oi.peak_on_detector,
        oi.detectors_overlap,
        oi.quality,
        oi.n_quality_points,
    )


def rot(axis: np.ndarray, ang: float) -> np.ndarray:
    axis = axis / np.linalg.norm(axis)
    K = np.array([[0, -axis[2], axis[1]], [axis[2], 0, -axis[0]], [-axis[1], axis[0], 0]])
    return np.eye(3) + np.sin(ang) * K + (1 - np.cos(ang)) * K @ K


def orientation_pool(truth: np.ndarray, n: int, max_deg: float = 3.0, seed: int = 0) -> np.ndarray:
    rng = np.random.default_rng(seed)
    return np.stack(
        [rot(rng.normal(size=3), np.deg2rad(rng.uniform(0, max_deg))) @ truth for _ in range(n)]
    )


def median_us(fn: Any, n_rep: int, per: int) -> Tuple[float, float, List[float]]:
    """(median, min, all) of fn() wall time in microseconds per ``per`` items."""
    ts = []
    for _ in range(n_rep):
        t0 = time.perf_counter()
        fn()
        ts.append((time.perf_counter() - t0) / per * 1e6)
    return float(np.median(ts)), float(min(ts)), ts


def time_arm(
    cf: VoxelCostFunction, verts: Any, phase: int, pool: np.ndarray, devices: List[str]
) -> Dict[str, Any]:
    n_eval = 400
    for o in pool[:30]:
        cf.evaluate(o, verts, phase)
    ev_med, ev_min, _ = median_us(
        lambda: [cf.evaluate(o, verts, phase) for o in pool[:n_eval]], 5, n_eval
    )
    ref = [cf.evaluate(o, verts, phase) for o in pool]
    # sync-forcing calls in one evaluate (CPU tensors: these calls are free here, not on a GPU)
    sc = SyncCounter(assume_gpu=True)
    with sc:
        cf.evaluate(pool[0], verts, phase)
    res: Dict[str, Any] = dict(
        K=int(cf._phase_recip_vecs[phase][0].shape[0]),
        evaluate_us_median=ev_med,
        evaluate_us_min=ev_min,
        evaluate_sync_calls=dict(sc.counts),
        evaluate_sync_total=int(sum(sc.counts.values())),
        n_valid_peaks_median=int(np.median([r.n_quality_points for r in ref])),
        batched=[],
    )
    for device in devices:
        for variant in ("fixed", "compact"):
            if variant == "compact" and device != "cpu":
                continue
            bc = ThrowawayBatchedCost(cf, phase, device, variant)
            sc = SyncCounter(assume_gpu=device == "cpu")
            with sc:
                bc.evaluate_many(pool[:16], verts)
            n_sync = int(sum(sc.counts.values()))
            # equality with evaluate on the whole pool, in chunks of 256
            n_diff, n_diff_counts, n_peak, n_cost = 0, 0, 0, 0
            dcost: List[float] = []
            for s in range(0, len(pool), 256):
                got = bc.evaluate_many(pool[s : s + 256], verts)
                for g_, r_ in zip(got, ref[s : s + 256]):
                    ga, ra = info_tuple(g_), info_tuple(r_)
                    if ga != ra:
                        n_diff += 1
                        n_diff_counts += int(ga[:5] != ra[:5] or ga[6] != ra[6])
                        n_peak += int(g_.peak_overlap != r_.peak_overlap)
                        if g_.cost != r_.cost:
                            n_cost += 1
                            dcost.append(abs(g_.cost - r_.cost))
            diff_stats = dict(
                n_peak_overlap_differs=n_peak,
                n_cost_differs=n_cost,
                max_abs_dcost=float(max(dcost)) if dcost else 0.0,
                median_abs_dcost=float(np.median(dcost)) if dcost else 0.0,
                median_cost_of_pool=float(np.median([r_.cost for r_ in ref])),
            )
            for B in BATCH_SIZES:
                batch = pool[:B]
                for _ in range(3):  # warm-up (shader compile on MPS)
                    bc.evaluate_many(batch, verts)
                if device == "mps":
                    torch.mps.synchronize()
                n_rep = 25 if B <= 16 else (9 if B <= 256 else 5)
                med, mn, _ = median_us(lambda: bc.evaluate_many(batch, verts), n_rep, B)
                tm: Dict[str, float] = {}
                for _ in range(n_rep):
                    bc.evaluate_many(batch, verts, timing=tm)
                res["batched"].append(
                    dict(
                        device=device,
                        variant=variant,
                        B=B,
                        us_median=med,
                        us_min=mn,
                        speedup_median=ev_med / med,
                        parts_us={k: v / n_rep / B * 1e6 for k, v in tm.items()},
                        sync_calls_B16=n_sync,
                        pool_n=len(pool),
                        n_info_differs=n_diff,
                        n_counts_differ=n_diff_counts,
                        **diff_stats,
                    )
                )
                print(
                    f"    {device:3s} {variant:7s} B={B:5d}: {med:7.1f} us/orient (min {mn:7.1f}),"
                    f" {ev_med / med:5.2f}x vs evaluate {ev_med:.0f} us; differs from evaluate in"
                    f" {n_diff}/{len(pool)}",
                    flush=True,
                )
    return res


def cmd_batched(a: argparse.Namespace) -> None:
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    info = preflight.preflight()
    try:
        preflight.require_quiet(info=info)
        info["quiet"] = True
    except preflight.MachineBusyError as exc:
        if not a.allow_busy:
            raise
        info["quiet"] = False
        info["not_quiet_reason"] = str(exc)
        print("WARNING: machine not quiet, timings are 'contended':", exc, flush=True)
    devices = ["cpu"] + (["mps"] if torch.backends.mps.is_available() and not a.no_mps else [])
    arms: List[Tuple[str, VoxelCostFunction, Any, int, np.ndarray]] = []
    if a.set == "three":
        S = ThreeVoxels()
        P = S.P
        verts, truth, phase = S.case()
        for name, args in (
            ("ThreeVoxels Q-max 8, default rule", (0, 8.0, False)),
            ("ThreeVoxels Q-max 8, C++ rule", (0, 8.0, True)),
            ("ThreeVoxels all reflections (Simulation.config), default rule", (0, 0.0, False)),
        ):
            arms.append((name, P.make(*args), verts, phase, truth))
    else:
        S = E0Set()
        verts, truth, phase, rec = S.load(a.vpos)
        arms.append(
            (
                f"E0 clean voxel {a.vpos}, default rule",
                rec._make_local_cost_fn(),
                verts,
                phase,
                truth,
            )
        )
        rec.setup.search_params.cpp_overlap_rule = True
        arms.append(
            (f"E0 clean voxel {a.vpos}, C++ rule", rec._make_local_cost_fn(), verts, phase, truth)
        )
    results: Dict[str, Any] = dict(
        preflight=info,
        torch=torch.__version__,
        torch_threads=torch.get_num_threads(),
        label="single-worker" if info.get("quiet") else "contended",
        arms={},
    )
    for name, cf, verts, phase, truth in arms:
        print(f"[{name}]", flush=True)
        pool = orientation_pool(truth, max(BATCH_SIZES))
        results["arms"][name] = time_arm(cf, verts, phase, pool, devices)
    (out / f"batched_{a.set}.json").write_text(json.dumps(results, indent=1, default=float))
    (out / f"preflight_{a.set}.json").write_text(json.dumps(info, indent=1, default=float))


# ---------------------------------------------------------------------------------------------
# summary
# ---------------------------------------------------------------------------------------------


ARM_OF_TRACE = {  # trace arm -> batched-timing arm with the same cost function
    "Example2.Simulation.config, cpp rule False": (
        "ThreeVoxels all reflections (Simulation.config), default rule"
    ),
    "ReconstructQ8.config, cpp rule False": "ThreeVoxels Q-max 8, default rule",
    "ReconstructQ8.config, cpp rule True": "ThreeVoxels Q-max 8, C++ rule",
    "ReconstructQ8.config (E0 clean), cpp rule False": "E0 clean voxel 0, default rule",
}


def interp_total_us(points: Dict[int, float], b: float) -> float:
    """Wall microseconds of one batched call of b orientations, from the measured us/orientation
    at the batch sizes of ``points`` (linear in b between them, and beyond the largest)."""
    xs = sorted(points)
    tot = [x * points[x] for x in xs]
    if b <= xs[0]:
        return tot[0] * b / xs[0]
    for lo, hi, t_lo, t_hi in zip(xs, xs[1:], tot, tot[1:]):
        if b <= hi:
            return t_lo + (t_hi - t_lo) * (b - lo) / (hi - lo)
    slope = (tot[-1] - tot[-2]) / (xs[-1] - xs[-2])
    return tot[-1] + slope * (b - xs[-1])


def model_run_us(
    levels: List[Dict[str, Any]], points: Dict[int, float], c_rot: float, scheme: str, k: int = 1
) -> float:
    """Modelled quick-MC wall time (us) of one run with the lockstep function. ``scheme``: 'plan'
    (option A as written: every chain runs all its 10 trials, k chains per unresolved
    candidate), 'memo' (distinct trial evaluations only, simulated), 'single' (one round, no
    speculation: draws that do not depend on other candidates). Each step costs one batched call
    plus c_rot per orientation for building the trial rotations (a scalar Python loop today)."""
    total = 0.0
    for lv in levels:
        n = len(lv["failed"])
        total += interp_total_us(points, N_START_EVALS * n)
        if scheme == "single":
            batches = [float(n)]
        elif scheme == "plan":
            batches = [float(b) for b in option_a_rounds(lv["failed"], k)["batches"]]
        else:
            sim = lv["sim"][str(k)]
            batches = [e / N_TRIAL_STEPS for e in sim["new_trial_evals_per_round"]]
        for b in batches:
            total += N_TRIAL_STEPS * (interp_total_us(points, b) + c_rot * b)
    return total


def cmd_summary(a: argparse.Namespace) -> None:
    sys.path.insert(0, str(ICENINE_PY / "scripts" / "common"))
    import doc_tables

    out = Path(a.out)
    summary: Dict[str, Any] = {}
    tables: Dict[str, str] = {}
    batched: Dict[str, Any] = {}
    for tag in ("three", "e0"):
        p = out / f"batched_{tag}.json"
        if p.exists():
            batched.update(json.loads(p.read_text())["arms"])
    c_rot = a.c_rot

    for tag, fname in (
        ("three", "trace_three.json"),
        ("three_sim", "trace_three_sim.json"),
        ("e0", "trace_e0.json"),
    ):
        p = out / fname
        if not p.exists():
            continue
        runs = json.loads(p.read_text())
        groups: Dict[str, List[Dict[str, Any]]] = defaultdict(list)
        for r in runs:
            groups[f"{r['config']}, cpp rule {r['cpp_overlap_rule']}"].append(r)
        for gname, rs in groups.items():
            key = f"{tag}: {gname}"
            g: Dict[str, Any] = dict(n_runs=len(rs), levels=[])
            rows = []
            for L in range(max(len(r["levels"]) for r in rs)):
                lv = [r["levels"][L] for r in rs if L < len(r["levels"])]
                nc = np.array([len(x["failed"]) for x in lv])
                nf = np.array([sum(x["failed"]) for x in lv])
                lo, hi = stats.wilson(int(nf.sum()), int(nc.sum()))
                row: Dict[str, Any] = dict(
                    level=L,
                    runs=len(lv),
                    cand_median=float(np.median(nc)),
                    cand_min=int(nc.min()),
                    cand_max=int(nc.max()),
                    failed=int(nf.sum()),
                    cand=int(nc.sum()),
                    rate=float(nf.sum() / nc.sum()),
                    lo=lo,
                    hi=hi,
                    draws_ok=all(d in (30, 33) for x in lv for d in x["draws"])
                    and all(
                        d == 30 + 3 * int(f) for x in lv for d, f in zip(x["draws"], x["failed"])
                    ),
                )
                barm0 = batched.get(ARM_OF_TRACE.get(gname, ""))
                if barm0 is not None:
                    for dev, variant in (("cpu", "compact"), ("cpu", "fixed"), ("mps", "fixed")):
                        pts0 = {
                            e["B"]: e["us_median"]
                            for e in barm0["batched"]
                            if e["device"] == dev and e["variant"] == variant
                        }
                        row[f"speedup_{dev}_{variant}"] = float(
                            np.median(
                                [
                                    barm0["evaluate_us_median"]
                                    * len(x["failed"])
                                    / interp_total_us(pts0, len(x["failed"]))
                                    for x in lv
                                ]
                            )
                        )
                for kk in KS:
                    oa = [option_a_rounds(x["failed"], kk) for x in lv]
                    row[f"rounds_k{kk}"] = float(np.mean([o["rounds"] for o in oa]))
                    row[f"comp_over_counted_k{kk}"] = float(
                        sum(o["computed_evals"] for o in oa) / sum(o["counted_evals"] for o in oa)
                    )
                rows.append(row)
                g["levels"].append(row)
            tables[f"s3_trace_{tag}_{len(tables)}"] = (
                f"**{key}** ({len(rs)} runs)\n\n"
                + doc_tables.markdown_table(
                    rows,
                    [
                        "level",
                        "runs",
                        "cand_median",
                        "cand_min",
                        "cand_max",
                        "failed",
                        "cand",
                        "rate",
                        "lo",
                        "hi",
                        "rounds_k1",
                        "rounds_k2",
                        "rounds_k4",
                        "comp_over_counted_k1",
                        "comp_over_counted_k2",
                        "comp_over_counted_k4",
                        "draws_ok",
                        "speedup_cpu_compact",
                        "speedup_cpu_fixed",
                        "speedup_mps_fixed",
                    ],
                    formats=dict(
                        rate=".3f",
                        lo=".3f",
                        hi=".3f",
                        cand_median=".0f",
                        rounds_k1=".1f",
                        rounds_k2=".1f",
                        rounds_k4=".1f",
                        comp_over_counted_k1=".1f",
                        comp_over_counted_k2=".1f",
                        comp_over_counted_k4=".1f",
                        speedup_cpu_compact=".1f",
                        speedup_cpu_fixed=".1f",
                        speedup_mps_fixed=".1f",
                    ),
                )
            )
            # shares of the search
            wall = [r["wall_plain_s"] for r in rs]
            q_s = [r["optimize_seconds"].get("quick_mc", 0.0) for r in rs]
            q_n = [r["evaluate_calls"].get("quick_mc", dict(n=0))["n"] for r in rs]
            all_n = [sum(v["n"] for v in r["evaluate_calls"].values()) for r in rs]
            all_s = [sum(v["seconds"] for v in r["evaluate_calls"].values()) for r in rs]
            inst = [r["wall_s"] for r in rs]
            g.update(
                wall_plain_s_median=float(np.median(wall)),
                quick_mc_seconds_median=float(np.median(q_s)),
                quick_mc_share_of_wall=float(np.median([q / w for q, w in zip(q_s, inst)])),
                quick_mc_evaluate_calls_median=float(np.median(q_n)),
                evaluate_calls_median=float(np.median(all_n)),
                quick_mc_share_of_evaluate_calls=float(
                    np.median([q / t for q, t in zip(q_n, all_n)])
                ),
                evaluate_seconds_share_of_wall=float(
                    np.median([s_ / w for s_, w in zip(all_s, inst)])
                ),
                quick_mc_us_per_evaluation=float(
                    np.median([1e6 * q / n for q, n in zip(q_s, q_n)])
                ),
                instrumented_equals_plain=all(r["instrumented_equals_plain"] for r in rs),
                global_evals_median=float(np.median([r["eval_counts"]["global_"] for r in rs])),
                local_evals_median=float(np.median([r["eval_counts"]["local"] for r in rs])),
            )
            # simulation (option A for real, with and without reuse of trial evaluations)
            if all("sim" in lv for r in rs for lv in r["levels"]):
                sim_rows = []
                for kk in sorted({int(x) for x in rs[0]["levels"][0]["sim"]}):
                    per_run = []
                    for r in rs:
                        sl = [lv["sim"][str(kk)] for lv in r["levels"]]
                        per_run.append(
                            dict(
                                rounds=sum(x["rounds"] for x in sl),
                                plan=sum(x["plan_trial_evals"] + x["start_evals"] for x in sl),
                                memo=sum(x["memo_trial_evals"] + x["start_evals"] for x in sl),
                                counted=sum(x["counted_evals"] for x in sl),
                                ok=all(x["matches_sequential"] for x in sl),
                            )
                        )
                    row = dict(
                        k=kk,
                        rounds_per_run=float(np.mean([x["rounds"] for x in per_run])),
                        counted=float(np.mean([x["counted"] for x in per_run])),
                        plan=float(np.mean([x["plan"] for x in per_run])),
                        memo=float(np.mean([x["memo"] for x in per_run])),
                        plan_over_counted=float(
                            sum(x["plan"] for x in per_run) / sum(x["counted"] for x in per_run)
                        ),
                        memo_over_counted=float(
                            sum(x["memo"] for x in per_run) / sum(x["counted"] for x in per_run)
                        ),
                        all_match=all(x["ok"] for x in per_run),
                    )
                    sim_rows.append(row)
                g["simulation"] = sim_rows
                tables[f"s3_sim_{tag}"] = (
                    f"**{key}**, option A run for real ({len(rs)} runs; "
                    "sums over the levels of a run)\n\n"
                    + doc_tables.markdown_table(
                        sim_rows,
                        [
                            "k",
                            "rounds_per_run",
                            "counted",
                            "plan",
                            "memo",
                            "plan_over_counted",
                            "memo_over_counted",
                            "all_match",
                        ],
                        formats=dict(
                            rounds_per_run=".1f",
                            counted=".0f",
                            plan=".0f",
                            memo=".0f",
                            plan_over_counted=".1f",
                            memo_over_counted=".1f",
                        ),
                    )
                )
            # modelled quick-MC time with the measured batched cost
            barm = batched.get(ARM_OF_TRACE.get(gname, ""))
            if barm is not None:
                model_rows = []
                seq_us = float(np.median([1e6 * q for q in q_s]))
                for dev, variant in (("cpu", "compact"), ("cpu", "fixed"), ("mps", "fixed")):
                    pts = {
                        e["B"]: e["us_median"]
                        for e in barm["batched"]
                        if e["device"] == dev and e["variant"] == variant
                    }
                    if not pts:
                        continue
                    for scheme, kk in (
                        ("single", 1),
                        ("plan", 1),
                        ("plan", 2),
                        ("plan", 4),
                        ("memo", 4),
                        ("memo", 16),
                    ):
                        if scheme == "memo" and "sim" not in rs[0]["levels"][0]:
                            continue
                        for crot in (c_rot, 0.0):
                            m = float(
                                np.median(
                                    [model_run_us(r["levels"], pts, crot, scheme, kk) for r in rs]
                                )
                            )
                            model_rows.append(
                                dict(
                                    device=dev,
                                    variant=variant,
                                    scheme=scheme,
                                    k=kk,
                                    c_rot=crot,
                                    modelled_ms=m / 1e3,
                                    sequential_ms=seq_us / 1e3,
                                    ratio=seq_us / m,
                                    whole_search=1.0
                                    / (
                                        1.0
                                        - g["quick_mc_share_of_wall"]
                                        + g["quick_mc_share_of_wall"] * m / seq_us
                                    ),
                                )
                            )
                g["model"] = model_rows
                tables[f"s3_model_{tag}_{len(tables)}"] = (
                    f"**{key}**: modelled quick-MC time per run vs the measured sequential "
                    f"quick MC (trial-rotation build c_rot us per orientation); whole_search = "
                    f"Amdahl estimate from the quick-MC share of the wall time\n\n"
                    + doc_tables.markdown_table(
                        model_rows,
                        [
                            "device",
                            "variant",
                            "scheme",
                            "k",
                            "c_rot",
                            "modelled_ms",
                            "sequential_ms",
                            "ratio",
                            "whole_search",
                        ],
                        formats=dict(
                            c_rot=".1f",
                            modelled_ms=".0f",
                            sequential_ms=".0f",
                            ratio=".2f",
                            whole_search=".2f",
                        ),
                    )
                )
            summary[key] = g

    # batched timings
    for name, arm in batched.items():
        rows = []
        for e in arm["batched"]:
            rows.append(
                dict(
                    device=e["device"],
                    variant=e["variant"],
                    B=e["B"],
                    us=e["us_median"],
                    us_min=e["us_min"],
                    speedup=e["speedup_median"],
                    torch_us=e["parts_us"].get("torch"),
                    host_us=e["parts_us"].get("host_prep"),
                    stage_d_us=e["parts_us"].get("stage_d"),
                    differs=e["n_info_differs"],
                    counts_differ=e["n_counts_differ"],
                    peak_differs=e.get("n_peak_overlap_differs"),
                    cost_differs=e.get("n_cost_differs"),
                    max_dcost=e.get("max_abs_dcost"),
                    med_dcost=e.get("median_abs_dcost"),
                    pool=e["pool_n"],
                    syncs=e["sync_calls_B16"],
                )
            )
        summary.setdefault("batched", {})[name] = dict(
            K=arm["K"],
            evaluate_us_median=arm["evaluate_us_median"],
            evaluate_us_min=arm["evaluate_us_min"],
            evaluate_sync_total=arm["evaluate_sync_total"],
            evaluate_sync_calls=arm["evaluate_sync_calls"],
            n_valid_peaks_median=arm["n_valid_peaks_median"],
            rows=rows,
        )
        tables[f"s3_batched_{len(tables)}"] = (
            f"**{name}**: evaluate {arm['evaluate_us_median']:.0f} us "
            f"(min {arm['evaluate_us_min']:.0f}), K = {arm['K']} reflections, "
            f"{arm['evaluate_sync_total']} sync-forcing calls per evaluate\n\n"
            + doc_tables.markdown_table(
                rows,
                [
                    "device",
                    "variant",
                    "B",
                    "us",
                    "us_min",
                    "speedup",
                    "torch_us",
                    "host_us",
                    "stage_d_us",
                    "differs",
                    "counts_differ",
                    "peak_differs",
                    "cost_differs",
                    "max_dcost",
                    "med_dcost",
                    "pool",
                    "syncs",
                ],
                formats=dict(
                    us=".1f",
                    us_min=".1f",
                    speedup=".2f",
                    torch_us=".1f",
                    host_us=".1f",
                    stage_d_us=".1f",
                    max_dcost=".4f",
                    med_dcost=".5f",
                ),
            )
        )
    doc_tables.write_tables(out / "tables.md", tables)
    (out / "summary.json").write_text(json.dumps(summary, indent=1, default=float))
    print("wrote", out / "tables.md", "and", out / "summary.json")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    ap.add_argument("cmd", choices=["trace", "batched", "summary"])
    ap.add_argument("--set", choices=["three", "e0"], default="three")
    ap.add_argument("--out", required=True)
    ap.add_argument("--seeds", type=int, default=3)
    ap.add_argument("--n-voxels", type=int, default=10)
    ap.add_argument("--vpos", type=int, default=0)
    ap.add_argument("--arms", default="", help="three: arm indices to trace (0,1,2), default all")
    ap.add_argument("--tag", default="", help="suffix of the trace json name")
    ap.add_argument("--sim-ks", default="", help="speculation depths to run for real (trace)")
    ap.add_argument("--allow-busy", action="store_true")
    ap.add_argument("--no-mps", action="store_true")
    ap.add_argument("--c-rot", type=float, default=14.0, help="summary: us per trial rotation")
    a = ap.parse_args()
    a.out = str(Path(a.out).resolve())  # the example sets chdir
    {"trace": cmd_trace, "batched": cmd_batched, "summary": cmd_summary}[a.cmd](a)


if __name__ == "__main__":
    main()
