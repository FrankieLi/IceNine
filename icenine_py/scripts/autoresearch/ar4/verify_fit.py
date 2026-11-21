"""AR4 step 1: train B1 and N1 on the train samples and write the file of thresholds.

  uv run python scripts/autoresearch/ar4/verify_fit.py

Train samples: ``big`` and ``grad`` (the first gradient sample), all their pairs
(``verify_inputs.py``). One training run per model, generator seed 0, no other seed; nothing here
reads or touches ``mid``. Writes, tracked: ``benchmarks/autoresearch/ar4/models.npz`` (B1 weights
and N1 weights, the standardization of the summary numbers) and
``benchmarks/autoresearch/ar4/thresholds.json`` (the use-1 cut of each score). A score is in
[0, 1] (higher: more likely right); a voxel or seed is flagged when its score is at or under
the cut.

Use-1 cut of a score: the cut that flags the fewest voxels while flagging at least 95% of the wrong
voxels (error above 1 degree) of the S8 maps after the boundary pass, the two train samples pooled
(pairs of kind ``after``). Voxels tied with the cut are all flagged.
"""

import json
import time
from pathlib import Path
from typing import Any, Dict, Optional

import numpy as np
import torch
import verify_common as VC
import verify_lib as VL
import verify_models as VM

KIND_AFTER = 0
SCORES = ("B0", "B1", "N1", "lit")


def load_inputs(sample: str, suffix: str = "") -> Dict[str, np.ndarray]:
    z = np.load(VC.CACHE / f"inputs_{sample}{suffix}.npz")
    return {k: z[k] for k in z.files}


def concat(parts: Any) -> Dict[str, np.ndarray]:
    out: Dict[str, np.ndarray] = {}
    for k in parts[0]:
        if k == "offsets":
            continue
        out[k] = np.concatenate([p[k] for p in parts])
    base = 0
    offs = [np.zeros(1, np.int64)]
    for p in parts:
        offs.append(p["offsets"][1:] + base)
        base += int(p["offsets"][-1])
    out["offsets"] = np.concatenate(offs)
    return out


def standardizer(summary: np.ndarray) -> Dict[str, np.ndarray]:
    s = summary[:, list(VM.SUMMARY_USED)]
    mu, sd = s.mean(axis=0), s.std(axis=0)
    sd[sd < 1e-12] = 1.0
    return dict(mu=mu, sd=sd)


def _summ(inp: Dict[str, np.ndarray], st: Dict[str, np.ndarray]) -> np.ndarray:
    s = inp["summary"][:, list(VM.SUMMARY_USED)]
    return ((s - st["mu"]) / st["sd"]).astype(np.float32)


def train_n1(inp: Dict[str, np.ndarray], st: Dict[str, np.ndarray]) -> Dict[str, Any]:
    """One training run of N1 (generator seed 0)."""
    torch.manual_seed(VM.SEED)
    gen = torch.Generator().manual_seed(VM.SEED)
    net = VM.SetNet()
    y = VL.right_label(inp["err"]).astype(np.float32)
    n_pos, n_neg = float(y.sum()), float(len(y) - y.sum())
    w_pos, w_neg = 0.5 / max(n_pos, 1.0), 0.5 / max(n_neg, 1.0)
    summ = _summ(inp, st)
    opt = torch.optim.Adam(net.parameters(), lr=VM.LR, weight_decay=VM.WEIGHT_DECAY)
    P = len(y)
    t0 = time.time()
    losses = []
    for ep in range(VM.EPOCHS):
        perm = torch.randperm(P, generator=gen).numpy()
        tot = 0.0
        for lo in range(0, P, VM.BATCH):
            idx = perm[lo : lo + VM.BATCH]
            x, m = VM.gather_batch(inp["peaks"], inp["offsets"], idx)
            yb = torch.from_numpy(y[idx])
            wb = torch.where(yb > 0.5, torch.tensor(w_pos), torch.tensor(w_neg)) * len(y)
            logit = net(torch.from_numpy(x), torch.from_numpy(m), torch.from_numpy(summ[idx]))
            loss = (
                torch.nn.functional.binary_cross_entropy_with_logits(logit, yb, reduction="none")
                * wb
            ).mean()
            opt.zero_grad()
            loss.backward()
            opt.step()
            tot += float(loss) * len(idx)
        losses.append(tot / P)
        print(
            f"N1 epoch {ep + 1}/{VM.EPOCHS} loss {losses[-1]:.4f} {time.time() - t0:.0f} s",
            flush=True,
        )
        assert time.time() - t0 < 1200, "N1 training passed twice its 10-minute budget"
    return dict(net=net, losses=losses, wall_s=time.time() - t0)


def n1_scores(net: VM.SetNet, inp: Dict[str, np.ndarray], st: Dict[str, np.ndarray]) -> np.ndarray:
    summ = _summ(inp, st)
    P = len(summ)
    out = np.zeros(P)
    net.eval()
    with torch.no_grad():
        for lo in range(0, P, 2048):
            idx = np.arange(lo, min(lo + 2048, P))
            x, m = VM.gather_batch(inp["peaks"], inp["offsets"], idx)
            out[idx] = torch.sigmoid(
                net(torch.from_numpy(x), torch.from_numpy(m), torch.from_numpy(summ[idx]))
            ).numpy()
    return out


def all_scores(
    inp: Dict[str, np.ndarray], models: Dict[str, Any], net: Optional[VM.SetNet]
) -> Dict[str, np.ndarray]:
    """The four scores of every pair (higher: more likely right)."""
    summ = inp["summary"]
    out = {
        "B0": summ[:, VL.SUMMARY_NAMES.index("hit")],
        "lit": summ[:, VL.SUMMARY_NAMES.index("lit_share")],
        "B1": VL.predict_logistic(
            summ[:, list(VM.SUMMARY_USED)],
            models["b1_w"],
            float(models["b1_b"]),
            models["b1_mu"],
            models["b1_sd"],
        ),
    }
    assert net is not None
    out["N1"] = n1_scores(net, inp, dict(mu=models["s_mu"], sd=models["s_sd"]))
    return out


def load_models(path: Path) -> Dict[str, Any]:
    z = np.load(path)
    return {k: z[k] for k in z.files}


def net_from(models: Dict[str, Any]) -> VM.SetNet:
    net = VM.SetNet()
    sd = {k[3:]: torch.from_numpy(v) for k, v in models.items() if k.startswith("n1_")}
    net.load_state_dict(sd)
    return net


def main() -> None:
    parts = [load_inputs(s) for s in VC.TRAIN]
    inp = concat(parts)
    P = len(inp["err"])
    y = VL.right_label(inp["err"])
    print("train pairs", P, "right", int(y.sum()), flush=True)
    st = standardizer(inp["summary"])
    # B1
    t0 = time.time()
    w, b, mu, sd = VL.fit_logistic(
        inp["summary"][:, list(VM.SUMMARY_USED)],
        y.astype(np.float64),
        epochs=VM.LOGISTIC_EPOCHS,
        lr=VM.LOGISTIC_LR,
    )
    wall_b1 = time.time() - t0
    # N1
    tr = train_n1(inp, st)
    net: VM.SetNet = tr["net"]
    models: Dict[str, Any] = dict(
        b1_w=w, b1_b=np.array(b), b1_mu=mu, b1_sd=sd, s_mu=st["mu"], s_sd=st["sd"]
    )
    for k, v in net.state_dict().items():
        models["n1_" + k] = v.numpy()
    np.savez(VC.AR4 / "models.npz", **models)
    sc = all_scores(inp, models, net)
    use1 = inp["kind"] == KIND_AFTER
    wrong = ~y
    out: Dict[str, Any] = dict(
        train_samples=list(VC.TRAIN),
        rule="cut = flags the fewest voxels while flagging >= 95% of the wrong voxels (error > 1 "
        "deg) of the S8 maps after the pass, train samples pooled; flagged: score <= cut; ties "
        "at the cut are all flagged",
        n_pairs=int(P),
        n_use1_voxels=int(use1.sum()),
        n_use1_wrong=int((use1 & wrong).sum()),
        n1_parameters=VM.n_parameters(net),
        n1_epochs=VM.EPOCHS,
        n1_final_loss=tr["losses"][-1],
        n1_wall_s=tr["wall_s"],
        b1_wall_s=wall_b1,
        b1_weights={n: float(x) for n, x in zip([VL.SUMMARY_NAMES[i] for i in VM.SUMMARY_USED], w)},
        generator_seed=VM.SEED,
        thresholds={},
        timing_label="contended (another agent on the machine); not a timing result",
    )
    for name in SCORES:
        s = sc[name][use1]
        cut = VL.pick_threshold(s, wrong[use1], 0.95)
        fl, wf = VL.flag_counts(s, wrong[use1], cut)
        out["thresholds"][name] = dict(
            cut=cut,
            train_voxels_flagged=fl,
            train_wrong_flagged=wf,
            train_flagged_share=fl / int(use1.sum()),
        )
    (VC.AR4 / "thresholds.json").write_text(json.dumps(out, indent=1) + "\n")
    print(json.dumps(out["thresholds"]), flush=True)


if __name__ == "__main__":
    main()
