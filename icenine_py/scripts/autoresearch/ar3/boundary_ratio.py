"""AR3 step 5, reading 6 (candidate feature): stored acceptance ratio of P4's fit-settled voxels.

  uv run python scripts/autoresearch/ar3/boundary_ratio.py

Step 4's review found that every wrong voxel among the fits the caller accepted in P4 has an
acceptance ratio (hit ratio of the fit over the entry's best hit ratio) below 0.95. This is a
candidate for the skip rule: skip the tries of a fit-settled voxel whose stored ratio is at least
theta. The ratio exists only for P4's fit-settled voxels (``paint.json.gz``: the fit and its entry).
D4's and E4's records do not store the entry's best value (the BFS ``_expand`` keeps it in a local
variable); there the skip rule of ``boundary_summary.py`` uses, in its place, the voxel's hit ratio
over the highest hit ratio among its neighbours in the map. Here the rule is only listed and
measured on P4R4, picked on nothing: a fixed grid, shown on A and B. Writes
``boundary/ratio.json`` and the table ``ar3e_ratio`` in ``boundary/tables.md``.
"""

import json

import boundary_common as B
import boundary_lib as BL
import numpy as np
import run_gradient as RG
from fine_summary import wil, wil_str

import doc_tables as DT  # noqa: E402

GRID = [0.90, 0.92, 0.95, 0.98]


def main() -> None:
    env = B.Env("grad")
    R0 = B.load_map("P4R4", env.n)
    e0 = env.err(R0)
    t = BL.Tries(B.RUNS / "P4R4")
    after, _took = BL.apply_rule(t, R0, "hit")
    e1 = env.err(after)
    res = RG.load_paint(B.FINE_RUNS / "P4")
    ratio = np.full(env.n, np.nan)
    for v, f in res.fits.items():
        e = res.entries[f.from_grain]
        ratio[v] = f.hit / max(e.best_hit, f.hit)
    have = ~np.isnan(ratio)
    wb, wa = e0 > B.WRONG, e1 > B.WRONG
    fixed, broken = wb & ~wa, ~wb & wa
    tv = BL.tries_per_voxel(t, env.n)
    rows = []
    out = {"n_fit_settled": int(have.sum()), "wrong_fit_settled": int((wb & have).sum())}
    out["wrong_fit_settled_max_ratio"] = float(np.nanmax(np.where(wb & have, ratio, np.nan)))
    for part, sel in (("A", env.isA), ("B", ~env.isA)):
        for th in GRID:
            sk = have & (ratio >= th) & sel
            rows.append(
                dict(
                    set=part,
                    theta=th,
                    fit_settled=int((have & sel).sum()),
                    skipped=int(sk.sum()),
                    wrong_skipped=int((wb & sk).sum()),
                    tries_skipped=int(tv[sk].sum()),
                    fixed_lost=int((fixed & sk).sum()),
                    fixed_in_set=int((fixed & sel).sum()),
                    broken_avoided=int((broken & sk).sum()),
                    lost_ci=wil_str(wil(int((fixed & sk).sum()), max(int((fixed & sel).sum()), 1))),
                )
            )
    out["rows"] = rows
    (B.BND / "ratio.json").write_text(json.dumps(out, indent=1) + "\n")
    tabs = DT.read_tables(B.BND / "tables.md")
    tabs["ar3e_ratio"] = DT.markdown_table(
        rows,
        [
            "set",
            "theta",
            "fit_settled",
            "skipped",
            "wrong_skipped",
            "tries_skipped",
            "fixed_lost",
            "fixed_in_set",
            "broken_avoided",
        ],
        {"theta": ".2f"},
    )
    DT.write_tables(B.BND / "tables.md", tabs)
    print(json.dumps({k: v for k, v in out.items() if k != "rows"}), flush=True)


if __name__ == "__main__":
    main()
