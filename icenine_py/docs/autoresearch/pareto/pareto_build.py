"""Build the cost-against-accuracy section (inline SVG panels + table) from pareto_rows.json.

Usage: pareto_build.py -> writes pareto_parts.json and preview SVGs (pareto_prev_*.svg, not
committed), prints the frontiers. Rows: pareto_rows.json, every number read from repo files.
Label text and offsets per point: LABELS. Look at the previews before publishing.
"""

import html
import json
import os
import math

S = os.path.dirname(os.path.abspath(__file__)) + "/"
d = json.load(open(S + "pareto_rows.json"))

PANELS = [
    dict(
        id="gradient_24570",
        big=True,
        title="Whole layer, grains with orientation gradients (24,570 voxels)",
    ),
    dict(
        id="layer_24570", big=True, title="Whole layer, one orientation per grain (24,570 voxels)"
    ),
    dict(id="pilot_2000", big=False, title="Pilot region (2,000 voxels), whole run"),
    dict(
        id="pilot_full_searches_69",
        big=False,
        title="Pilot: the 69 voxels that need a global search",
    ),
    dict(id="e0_200", big=False, title="200 hard voxels, one global search each"),
    dict(id="holdout_300", big=False, title="Holdout (300 voxels), no grid-search arm"),
]

# label text and offset for chosen points: method -> (text, dx, dy, anchor)
LABELS = {
    "gradient_24570": {
        "BFS indexed seeds, CMA local optimization": ("BFS, indexed seeds, CMA", 12, -8, "start"),
        "BFS indexed, CMA + boundary pass": ("+ boundary pass", -14, 4, "end"),
        "BFS indexed, CMA + 2 boundary passes": ("+ 2 passes", 0, 25, "middle"),
        "Paint + stock + 2 boundary passes": ("paint + stock + 2 passes", 0, 25, "middle"),
        "BFS indexed seeds, stock local optimization": (
            "BFS, indexed seeds, stock (the base run)",
            -11,
            -8,
            "end",
        ),
        "Paint grains (P)": ("paint grains", 6, -12, "end"),
        "Paint + stock local optimization (P plus R)": ("paint + stock", -6, -12, "start"),
        "Paint + CMA local optimization (P4 plus R4)": ("paint + CMA", 12, -8, "start"),
        "BFS indexed, stock + boundary pass": ("BFS stock + pass", 0, -13, "middle"),
        "Paint + stock + boundary pass": ("paint + stock + pass", -12, 4, "end"),
        "Paint + CMA + boundary pass": ("paint + CMA + pass", 12, -5, "start"),
    },
    "layer_24570": {
        "Paint grains (run C)": ("paint grains", 12, -10, "start"),
        "BFS indexed, lib. CMA + boundary pass": (
            "BFS, library CMA + boundary pass",
            22,
            -24,
            "end",
        ),
        "BFS indexed, CMA + boundary pass": ("BFS CMA + pass", -12, -24, "start"),
        "BFS indexed seeds, stock local optimization": (
            "BFS, indexed seeds, stock (the base run)",
            -11,
            -8,
            "end",
        ),
        "BFS indexed seeds, library CMA keys": ("BFS, library CMA", 11, 4, "start"),
        "BFS indexed seeds, CMA local optimization": ("BFS, CMA", 11, 4, "start"),
        "BFS indexed, stock + boundary pass": ("BFS stock + pass", 0, -13, "middle"),
    },
    "pilot_2000": {
        "Paint grains (FusedCost)": ("paint grains", 14, -12, "start"),
        "BFS, grid search for every seed (key off)": (
            "BFS, a grid search for every seed",
            -11,
            4,
            "end",
        ),
        "BFS with indexed seeds (key on)": ("BFS, indexed seeds", 0, -14, "middle"),
        "BFS grid-search seeds, FusedCost (key off)": ("the same, FusedCost on", -11, 4, "end"),
    },
    "pilot_full_searches_69": {
        "Ray indexing, best of first 1 candidate": (
            "ray indexing, first candidate only",
            -4,
            -14,
            "start",
        ),
        "Ray indexing, best of first 5 candidates": ("best of 5", 12, -12, "start"),
        "Ray indexing + local optimization, best 1": (
            "+ local optimization, best 1",
            4,
            -12,
            "end",
        ),
        "Ray indexing + local optimization, best 2": (
            "+ local optimization, best 2",
            12,
            -14,
            "start",
        ),
        "Grid search (stored searches)": ("grid search", 0, -13, "middle"),
    },
    "e0_200": {
        "Ray indexing, best of first 1 candidate": (
            "ray indexing, first candidate only",
            -4,
            -15,
            "start",
        ),
        "Ray indexing + local optimization, best 1": (
            "+ local optimization on the best 1, 2, 3",
            -4,
            -15,
            "start",
        ),
        "Grid search, default overlap rule": ("grid search", -11, 4, "end"),
        "Grid search, C++ overlap rule": ("grid search, C++ overlap rule", -11, 4, "end"),
    },
    "holdout_300": {
        "Ray indexing + local optimization, best 1": (
            "ray indexing + local optimization on the best 1",
            12,
            -8,
            "start",
        ),
        "Ray indexing + local optimization, best 2": ("on the best 2", 14, 21, "start"),
        "Ray indexing + local optimization, best 3": ("on the best 3", 0, -14, "middle"),
    },
}

LIGHT = dict(
    surface="#ffffff",
    ink="#142033",
    ink2="#4a5568",
    ink3="#74808f",
    rule="#d9dfe6",
    ray="#2a78d6",
    bfs="#eb6834",
    paint="#1baf7a",
)


def wilson(k, n):
    z = 1.959964
    p = k / n
    den = 1 + z * z / n
    c = (p + z * z / (2 * n)) / den
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / den
    return max(0.0, c - h) * 100, min(1.0, c + h) * 100


def fam(r):
    f = r["family"]
    if f == "grid search":
        return "grid"
    if f in ("BFS with indexed seeds",):
        return "bfs"
    if f in ("paint grains", "paint + local optimization"):
        return "paint"
    assert f in ("ray indexing + local optimization", "other"), f
    assert r["method"].startswith("Ray indexing"), r["method"]
    return "ray"


def fe(x):
    if x >= 100:
        return f"{x:,.0f}"
    return f"{x:,.0f}" if abs(x - round(x)) < 1e-9 else f"{x:,.1f}"


def fp(p):
    if p == 0:
        return "0"
    return f"{p:.2f}" if p < 1 else f"{p:.1f}"


def nice_step(span, target=4):
    raw = span / target
    mag = 10 ** math.floor(math.log10(raw))
    for m in (1, 2, 2.5, 5, 10):
        if raw <= m * mag:
            return m * mag
    return 10 * mag


def prep(sid):
    rows = [r for r in d["rows"] if r["sample"] == sid]
    for r in rows:
        r["pct"] = 100.0 * r["wrong"] / r["n"]
        r["lo"], r["hi"] = wilson(r["wrong"], r["n"])
        st = r.get("interval_pct")
        if st:
            assert abs(st[0] - r["lo"]) < 0.06 and abs(st[1] - r["hi"]) < 0.06, (
                r["method"],
                st,
                r["lo"],
                r["hi"],
            )
        r["off"] = "keep rule" in r["method"]
        r["fam"] = fam(r)
    drawn = [r for r in rows if not r["off"]]
    best = float("inf")
    for r in rows:
        r["front"] = False
    for r in sorted(drawn, key=lambda r: (r["evals_per_voxel"], r["pct"])):
        if r["pct"] < best:
            r["front"] = True
            best = r["pct"]
    return rows, drawn


def panel(p, preview=False):
    rows, drawn = prep(p["id"])
    W, H = (860, 400) if p["big"] else (860, 300)
    ml, mr, mt, mb = 52, 34, 40, 46
    xs = [r["evals_per_voxel"] for r in drawn]
    lx0, lx1 = math.log10(min(xs)), math.log10(max(xs))
    pad = max(0.12, 0.09 * (lx1 - lx0))
    lx0, lx1 = lx0 - pad, lx1 + pad
    ymaxd = max(r["pct"] for r in drawn)
    step = nice_step(ymaxd * 1.12 if ymaxd > 0 else 1.0)
    ymax = step * math.ceil(ymaxd * 1.12 / step - 1e-9) if ymaxd > 0 else 1.0

    def X(v):
        return ml + (math.log10(v) - lx0) / (lx1 - lx0) * (W - ml - mr)

    def Y(v):
        return H - mb - v / ymax * (H - mt - mb)

    cands = []
    for k in range(-1, 7):
        for m in (1, 2, 5):
            v = m * 10.0**k
            if lx0 <= math.log10(v) <= lx1:
                cands.append(v)
    if len(cands) > 7:
        cands = [v for v in cands if abs(math.log10(v) - round(math.log10(v))) < 1e-9]
    if lx1 - lx0 < 1.2:
        cands = [
            m * 10.0**k
            for k in range(-1, 7)
            for m in (1, 1.5, 2, 3, 4, 5, 6, 8)
            if lx0 <= math.log10(m * 10.0**k) <= lx1
        ]
    o = []
    aria = (
        p["title"]
        + ". "
        + "; ".join(
            f"{r['method']}: {fe(r['evals_per_voxel'])} evaluations per voxel, {r['wrong']} of "
            f"{r['n']} wrong"
            for r in drawn
        )
    )
    o.append(f'<svg viewBox="0 0 {W} {H}" role="img" aria-label="{html.escape(aria, quote=True)}">')
    if preview:
        c = LIGHT
        o.append(
            f'<style>text{{font-family:Helvetica,Arial;fill:{c["ink2"]}}}.tick{{font-size:11.5px;'
            f'fill:{c["ink3"]}}}'
            f'.grid{{stroke:{c["rule"]};stroke-width:1}}.ptitle{{font-size:13.5px;font-weight:600;'
            f'fill:{c["ink"]}}}'
            f'.plab{{font-size:12px;fill:{c["ink"]}}}.front{{fill:none;stroke:{c["ink"]};'
            f"stroke-width:1.5}}"
            f'.pt{{stroke:{c["surface"]};stroke-width:2}}.pt-grid{{fill:{c["ink3"]}}}.'
            f'pt-ray{{fill:{c["ray"]}}}'
            f'.pt-bfs{{fill:{c["bfs"]}}}.pt-paint{{fill:{c["paint"]}}}.base{{stroke:{c["ink3"]};'
            f'stroke-width:1}}.fring{{fill:none;stroke:{c["ink"]};stroke-width:1.5}}</style>'
            f'<rect width="{W}" height="{H}" fill="{c["surface"]}"/>'
        )
    o.append(f'<text class="ptitle" x="{ml - 40}" y="18">{html.escape(p["title"])}</text>')
    o.append("<g>")
    yv = 0.0
    while yv <= ymax + 1e-9:
        y = Y(yv)
        cls = "base" if yv == 0 else "grid"
        o.append(f'<line class="{cls}" x1="{ml}" y1="{y:.1f}" x2="{W - mr}" y2="{y:.1f}"/>')
        lab = ("0" if yv == 0 else (f"{yv:g}")) + "%"
        o.append(f'<text class="tick" x="{ml - 8}" y="{y + 4:.1f}" text-anchor="end">{lab}</text>')
        yv += step
    for v in cands:
        x = X(v)
        o.append(f'<line class="grid" x1="{x:.1f}" y1="{mt}" x2="{x:.1f}" y2="{H - mb}"/>')
        o.append(
            f'<text class="tick" x="{x:.1f}" y="{H - mb + 17}" text-anchor="middle">{fe(v)}</text>'
        )
    o.append(
        f'<text class="tick" x="{(ml + W - mr) / 2:.1f}" y="{H - 8}" text-anchor="middle">'
        "evaluations per voxel (log scale); left is cheaper</text>"
    )
    o.append(
        f'<text class="tick" x="{ml - 40}" y="{mt - 7}">% of voxels wrong; down is better</text>'
    )
    o.append("</g>")
    fr = sorted([r for r in drawn if r["front"]], key=lambda r: r["evals_per_voxel"])
    path = f'M{X(fr[0]["evals_per_voxel"]):.1f},{Y(fr[0]["pct"]):.1f}'
    for r in fr[1:]:
        path += f' H{X(r["evals_per_voxel"]):.1f} V{Y(r["pct"]):.1f}'
    path += f" H{W - mr}"
    o.append(f'<path class="front" d="{path}"/>')
    o.append("<g>")
    for r in sorted(drawn, key=lambda r: (r["front"], -r["pct"])):
        x, y = X(r["evals_per_voxel"]), Y(r["pct"])
        tip = (
            f'{r["method"]} ({r["step"]}): {fe(r["evals_per_voxel"])} evaluations per voxel; '
            f'{r["wrong"]:,} of {r["n"]:,} wrong, {fp(r["pct"])}% (95% interval {fp(r["lo"])} to '
            f'{fp(r["hi"])}%)' + ("; on the frontier" if r["front"] else "")
        )
        t = f"<title>{html.escape(tip)}</title>"
        if r["front"]:
            o.append(f'<circle class="fring" cx="{x:.1f}" cy="{y:.1f}" r="11"/>')
        if r["boundary_pass"]:
            s = 7.2
            o.append(
                f'<path class="pt pt-{r["fam"]}" d="M{x:.1f},{y - s:.1f} L{x + s:.1f},{y:.1f} '
                f"L{x:.1f},{y + s:.1f} "
                f'L{x - s:.1f},{y:.1f} Z">{t}</path>'
            )
        else:
            o.append(
                f'<circle class="pt pt-{r["fam"]}" cx="{x:.1f}" cy="{y:.1f}" r="5.8">{t}</circle>'
            )
    o.append("</g><g>")
    labs = LABELS.get(p["id"], {})
    for m in labs:
        assert any(r["method"] == m for r in drawn), (p["id"], m)
    for r in drawn:
        if r["method"] in labs:
            text, dx, dy, anc = labs[r["method"]]
            x, y = X(r["evals_per_voxel"]) + dx, Y(r["pct"]) + dy
            o.append(
                f'<text class="plab" x="{x:.1f}" y="{y:.1f}" text-anchor="{anc}">'
                f"{html.escape(text)}</text>"
            )
    o.append("</g></svg>")
    return "\n".join(o), rows, fr


def note(r):
    n = []
    f = r.get("flag") or ""
    if r["off"]:
        n.append("another keep rule of the pass; not drawn (same cost, more wrong)")
    if "derived" in f:
        n.append("cost is a sum of two stated numbers")
    if r["fam"] == "ray":
        n.append("the indexing's own vote is not counted")
    if "chosen after" in f:
        n.append("arm chosen after the best-1 result was seen")
    if "unresolved" in f:
        n.append(f.split(";")[0])
    if (
        "contested voxels only" in f
        or "contested voxels" in r["method"]
        or "contested" in r["method"]
    ):
        n.append("local optimization on the contested voxels only")
    if r["method"] == "Grid search, C++ overlap rule":
        n.append("evaluation total from the result file only")
    return "; ".join(n)


def build():
    svgs, table, fronts = {}, [], {}
    for p in PANELS:
        svg, rows, fr = panel(p)
        svgs[p["id"]] = svg
        fronts[p["id"]] = fr
        open(S + f"pareto_prev_{p['id']}.svg", "w").write(
            panel(p, preview=True)[0].replace(
                "<svg ", '<svg xmlns="http://www.w3.org/2000/svg" ', 1
            )
        )
        table.append(f'<tr><th colspan="6">{html.escape(p["title"])}</th></tr>')
        for r in sorted(rows, key=lambda r: (r["evals_per_voxel"], r["pct"])):
            table.append(
                "<tr><td>"
                + html.escape(r["method"])
                + ("" if not r["front"] else " <strong>(frontier)</strong>")
                + f'</td><td>{html.escape(r["step"])}</td><td class="num">'
                f'{fe(r["evals_per_voxel"])}</td>'
                f'<td class="num">{r["wrong"]:,} of {r["n"]:,}</td>'
                f'<td class="num">{fp(r["pct"])}% ({fp(r["lo"])} to {fp(r["hi"])})</td>'
                f"<td>{html.escape(note(r))}</td></tr>"
            )
    return svgs, "\n".join(table), fronts


if __name__ == "__main__":
    svgs, table, fronts = build()
    for k, fr in fronts.items():
        print(k)
        for r in fr:
            print(
                "   F",
                fe(r["evals_per_voxel"]),
                r["wrong"],
                "/",
                r["n"],
                fp(r["pct"]),
                "|",
                r["method"],
            )
    json.dump({"svgs": svgs, "table": table}, open(S + "pareto_parts.json", "w"))
    print("sizes", {k: len(v) for k, v in svgs.items()}, len(table))
