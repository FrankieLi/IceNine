"""Build the AR6 tracking page from docs/grain_opt/IDEAS.md, STATE.md ("## Now" section) and
benchmarks/grain_opt/progress.json. Writes benchmarks/grain_opt/page/index.html (self-contained).

Usage (from icenine_py/): uv run python scripts/grain_opt/build_page.py
"""

import html
import json
import math
import re
from pathlib import Path
from typing import Any, Dict, List

HERE = Path(__file__).resolve().parent
ICE = HERE.parents[1]
DOCS = ICE / "docs" / "grain_opt"
OUT = ICE / "benchmarks" / "grain_opt" / "page" / "index.html"

GOALS = {
    "G1a": "Pareto frontier of evaluations against quality, fixed objective",
    "G1b": "Robustness: linear gradients to random voxel-to-voxel misorientation (0-5 deg)",
    "G1c": "Sensitivity to elastic strain (lattice deformation)",
    "G2a": "Ways to use intensity",
    "G2b": "Richer diffraction simulation",
    "G2c": "Learned and neural models",
}
STATUS_ORDER = ["running", "proposed", "done", "dropped"]


def esc(s: Any) -> str:
    return html.escape(str(s))


def wilson(k: int, n: int, z: float = 1.96) -> List[float]:
    p = k / n
    den = 1 + z * z / n
    c = (p + z * z / (2 * n)) / den
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / den
    return [100 * p, 100 * max(c - h, 0), 100 * min(c + h, 1)]


def ideas() -> List[Dict[str, str]]:
    rows = []
    for line in (DOCS / "IDEAS.md").read_text().splitlines():
        if not re.match(r"^\| I\d+ ", line):
            continue
        cells = [c.strip() for c in line.strip().strip("|").split("|")]
        cells += [""] * (6 - len(cells))
        rows.append(dict(zip(["id", "goal", "idea", "status", "result", "findings"], cells)))
    return rows


def state_now() -> List[str]:
    text = (DOCS / "STATE.md").read_text()
    m = re.search(r"^## Now\n(.*?)(?=^## |\Z)", text, re.S | re.M)
    return [ln for ln in (m.group(1) if m else "").strip().splitlines() if ln.strip()]


def md_inline(s: str) -> str:
    s = esc(s)
    s = re.sub(r"`([^`]+)`", r"<code>\1</code>", s)
    s = re.sub(r"\*\*([^*]+)\*\*", r"<strong>\1</strong>", s)
    return s


# arm groups: (key, legend name, marker shape, colour token)
GROUPS = [
    ("ref", "Reference, all 19,013 painted grad voxels", "circle", "--ref"),
    ("start", "Painted start (G20 voxels)", "square", "--muted"),
    ("c", "AR3 per-voxel CMA, C (G20 voxels)", "diamond", "--run"),
    ("field", "Grain field, strict cost, N evaluations", "circle", "--accent"),
    ("quad", "Quadratic grain field, strict cost", "square", "--accent"),
    ("polish", "Grain field 100 + per-voxel polish", "triangle", "--done"),
    ("voxel", "Per-voxel CMA only (this program)", "diamond", "--prop"),
    ("smooth", "Sinkhorn / LS first, then strict", "cross", "--bad"),
]


def group_of(p: Dict[str, Any]) -> str:
    if p.get("group"):
        return str(p["group"])
    lab = p["label"]
    if p.get("family") == "reference":
        return "ref"
    if lab.startswith("start"):
        return "start"
    if lab.startswith("C "):
        return "c"
    if lab.startswith("quadratic"):
        return "quad"
    if "+ voxel" in lab:
        return "polish"
    if lab.startswith("voxel"):
        return "voxel"
    if lab.startswith("field"):
        return "field"
    return "smooth"


def marker(shape: str, x: float, y: float, r: float, color: str, extra: str = "") -> str:
    st = f'fill="var({color})" stroke="var({color})"'
    if shape == "circle":
        return f'<circle cx="{x:.1f}" cy="{y:.1f}" r="{r}" {st} {extra}/>'
    if shape == "square":
        return f'<rect x="{x-r:.1f}" y="{y-r:.1f}" width="{2*r}" height="{2*r}" {st} {extra}/>'
    if shape == "diamond":
        pts = f"{x:.1f},{y-r-1:.1f} {x+r+1:.1f},{y:.1f} {x:.1f},{y+r+1:.1f} {x-r-1:.1f},{y:.1f}"
        return f'<polygon points="{pts}" {st} {extra}/>'
    if shape == "triangle":
        pts = f"{x:.1f},{y-r-1:.1f} {x+r+1:.1f},{y+r:.1f} {x-r-1:.1f},{y+r:.1f}"
        return f'<polygon points="{pts}" {st} {extra}/>'
    d = r + 0.5  # cross
    return (f'<path d="M{x-d:.1f},{y-d:.1f}L{x+d:.1f},{y+d:.1f}M{x-d:.1f},{y+d:.1f}L{x+d:.1f},'
            f'{y-d:.1f}" stroke="var({color})" stroke-width="2.4" fill="none" {extra}/>')


def frontier(pts: List[Any]) -> List[Any]:
    """Lower-left frontier: points no other point beats on both renders and share (ties kept)."""
    pts = sorted(pts, key=lambda t: (t[0], t[1]))
    out, best = [], float("inf")
    for t in pts:
        if t[1] < best:
            out.append(t)
            best = t[1]
    return out


def pareto_svg(points: List[Dict[str, Any]]) -> str:
    W, H, L, R, T, B = 660, 320, 64, 18, 16, 50
    xs = [p["renders_per_voxel"] for p in points]
    xmax = max(max(xs), 100)
    ys = [wilson(p["k"], p["n"]) for p in points]
    ymax = max(5.0, max(y[2] for y in ys) * 1.12)
    gmap = {g[0]: g for g in GROUPS}

    def X(v: float) -> float:
        lv = math.log10(max(v, 1.0))
        return L + (W - L - R) * (0.07 + 0.93 * lv / math.log10(xmax * 1.25)) if v > 0 else L + 8

    def Y(v: float) -> float:
        return T + (H - T - B) * (1 - v / ymax)

    out = [f'<svg viewBox="0 0 {W} {H}" role="img" aria-label="Cost against quality" class="chart">']
    step = 1 if ymax <= 6 else (2 if ymax <= 12 else 5)
    v = 0.0
    while v <= ymax:
        out.append(f'<line x1="{L}" x2="{W-R}" y1="{Y(v):.1f}" y2="{Y(v):.1f}" class="grid"/>')
        out.append(f'<text x="{L-8}" y="{Y(v)+4:.1f}" class="tick" text-anchor="end">{v:g}%</text>')
        v += step
    for t in [0, 1, 10, 100, 1000]:
        if t > xmax * 1.25:
            break
        out.append(f'<text x="{X(t):.1f}" y="{H-B+18}" class="tick" text-anchor="middle">{t}</text>')
    out.append(f'<text x="{(L+W-R)/2:.0f}" y="{H-8}" class="axis" text-anchor="middle">'
               "num objective eval / voxel (log scale; 0 = no fit)</text>")
    out.append(f'<text x="16" y="{(T+H-B)/2:.0f}" class="axis" text-anchor="middle" '
               f'transform="rotate(-90 16 {(T+H-B)/2:.0f})">share of voxels above 0.1 deg (%)</text>')
    # frontier over the G20 points (excluding the reference set)
    g20 = [(p["renders_per_voxel"], c) for p, (c, _l, _h) in zip(points, ys) if group_of(p) != "ref"]
    fr = frontier(g20)
    if len(fr) > 1:
        path = " ".join(f"{X(x):.1f},{Y(y):.1f}" for x, y in fr)
        out.append(f'<polyline points="{path}" class="front"/>')
    for idx, (p, (c, lo, hi)) in enumerate(zip(points, ys), 1):
        g = gmap.get(group_of(p), GROUPS[-1])
        x, y = X(p["renders_per_voxel"]), Y(c)
        out.append(f'<line x1="{x:.1f}" x2="{x:.1f}" y1="{Y(lo):.1f}" y2="{Y(hi):.1f}" '
                   f'stroke="var({g[3]})" stroke-width="1.6" opacity="0.6"/>')
        tip = (f'<title>{esc(p["label"])}: {c:.2f}% ({lo:.2f} to {hi:.2f}), '
               f'{p["renders_per_voxel"]} objective eval / voxel, {esc(p.get("set", ""))}</title>')
        out.append(f"<g>{marker(g[2], x, y, 4.5, g[3])}{tip}</g>")
        out.append(f'<text x="{x+8:.1f}" y="{y-7:.1f}" class="lab">{idx}</text>')
    out.append("</svg>")
    used = []
    for p in points:
        k = group_of(p)
        if k not in used:
            used.append(k)
    items = []
    for k in [g[0] for g in GROUPS if g[0] in used]:
        g = gmap[k]
        items.append(f'<li><svg viewBox="0 0 18 18" width="18" height="18" aria-hidden="true">'
                     f'{marker(g[2], 9, 9, 4.5, g[3])}</svg>{esc(g[1])}</li>')
    items.append('<li><svg viewBox="0 0 26 18" width="26" height="18" aria-hidden="true">'
                 '<line x1="1" y1="9" x2="25" y2="9" class="front"/></svg>Frontier (G20 arms)</li>')
    legend = f'<ul class="legend">{"".join(items)}</ul>'
    howto = ('<p class="note">How to read: lower and further left is better (fewer voxels off by more '
             'than 0.1 deg, for fewer objective evaluations). Each marker is one arm; the vertical line through it is '
             'the Wilson 95% interval; the number next to it matches the table below. The dashed line '
             'joins the arms that no other arm beats on both axes.</p>')
    return "".join(out) + legend + howto


def build() -> str:
    prog = json.loads((ICE / "benchmarks" / "grain_opt" / "progress.json").read_text())
    rows = ideas()
    counts = {s: sum(1 for r in rows if r["status"] == s) for s in STATUS_ORDER}
    now = state_now()
    css = """
:root{--bg:#f6f7f9;--panel:#ffffff;--fg:#1c2330;--muted:#5d6879;--line:#dde2ea;
--accent:#2f6f9f;--run:#b26b00;--done:#2e7d4f;--drop:#8a8f98;--prop:#4a5a73;
--ref:#8a8f98;--bad:#b23b3b;--sans:"IBM Plex Sans",system-ui,sans-serif;--mono:"IBM Plex Mono",ui-monospace,monospace}
@media (prefers-color-scheme:dark){:root:not([data-theme="light"]){--bg:#11151b;--panel:#171c24;
--fg:#e3e8ef;--muted:#9aa5b5;--line:#2a313c;--accent:#6fb0de;--run:#e3a33f;--done:#5cbf86;
--drop:#7c838d;--prop:#a9b6cb;--ref:#7c838d;--bad:#e06a6a;color-scheme:dark}}
:root[data-theme="dark"]{--bg:#11151b;--panel:#171c24;--fg:#e3e8ef;--muted:#9aa5b5;--line:#2a313c;
--accent:#6fb0de;--run:#e3a33f;--done:#5cbf86;--drop:#7c838d;--prop:#a9b6cb;--ref:#7c838d;
--bad:#e06a6a;color-scheme:dark}
body{background:var(--bg);color:var(--fg);font-family:var(--sans);font-size:15px;line-height:1.5}
.wrap{max-width:1040px;margin:0 auto;padding-block:28px 48px;
padding-left:max(16px,env(safe-area-inset-left,0px));padding-right:max(16px,env(safe-area-inset-right,0px))}
h1{font-size:26px;margin:0 0 4px;text-wrap:balance}h2{font-size:18px;margin:32px 0 10px}
.sub{color:var(--muted);margin:0 0 18px}
.tiles{display:flex;flex-wrap:wrap;gap:10px}
.tile{background:var(--panel);border:1px solid var(--line);border-radius:8px;padding:10px 14px;min-width:110px}
.tile b{display:block;font-size:22px;font-variant-numeric:tabular-nums}
.tile span{color:var(--muted);font-size:12px;letter-spacing:.04em;text-transform:uppercase}
.now{background:var(--panel);border:1px solid var(--line);border-left:4px solid var(--accent);
border-radius:6px;padding:12px 16px}.now p{margin:4px 0}
.goals{display:grid;grid-template-columns:repeat(auto-fit,minmax(280px,1fr));gap:8px;margin:0;padding:0;list-style:none}
.goals li{background:var(--panel);border:1px solid var(--line);border-radius:6px;padding:8px 12px;font-size:14px}
.goals code{font-weight:600}
.tbl{overflow-x:auto;background:var(--panel);border:1px solid var(--line);border-radius:8px}
table{border-collapse:collapse;width:100%;font-size:14px}
th,td{text-align:left;padding:7px 10px;border-bottom:1px solid var(--line);vertical-align:top}
th{font-size:12px;color:var(--muted);text-transform:uppercase;letter-spacing:.04em}
td.id{font-family:var(--mono);white-space:nowrap}td.res{color:var(--muted);min-width:180px}
.pill{display:inline-block;font-size:12px;padding:1px 8px;border-radius:999px;border:1px solid currentColor;white-space:nowrap}
.running{color:var(--run)}.done{color:var(--done)}.dropped{color:var(--drop)}.proposed{color:var(--prop)}
code{font-family:var(--mono);font-size:.92em}
.chart{width:100%;height:auto;max-width:720px;background:var(--panel);border:1px solid var(--line);border-radius:8px}
.chart .grid{stroke:var(--line);stroke-width:1}.chart .tick{fill:var(--muted);font-size:11px;font-family:var(--mono)}
.chart .axis{fill:var(--muted);font-size:12px}.chart .lab{fill:var(--fg);font-size:11px}
.chart .dot{fill:var(--accent)}.chart .dot.reference{fill:var(--ref)}.chart .dot.g20ref{fill:var(--run)}
.chart .whisk.g20ref{stroke:var(--run)}
.chart .whisk{stroke:var(--accent);stroke-width:2}.chart .whisk.reference{stroke:var(--ref)}
.chart .front{fill:none;stroke:var(--fg);stroke-width:1.2;stroke-dasharray:5 4;opacity:.55}
.legend{list-style:none;padding:0;margin:10px 0 0;display:flex;flex-wrap:wrap;gap:6px 18px;font-size:13px;max-width:720px}
.legend li{display:flex;align-items:center;gap:6px}
.note{color:var(--muted);font-size:13px;max-width:70ch}
.filters{display:flex;gap:6px;flex-wrap:wrap;margin:0 0 8px}
.filters button{font:inherit;font-size:13px;background:var(--panel);color:var(--fg);border:1px solid var(--line);
border-radius:6px;padding:3px 10px;cursor:pointer}.filters button[aria-pressed="true"]{border-color:var(--accent);color:var(--accent)}
.filters button:focus-visible{outline:2px solid var(--accent);outline-offset:1px}
"""
    parts = [
        "<title>Grain Fit Research</title>",
        '<meta name="claude:viewport-fit-sides" content="cover">',
        '<link rel="preconnect" href="https://fonts.googleapis.com">',
        '<link rel="stylesheet" href="https://fonts.googleapis.com/css2?family=IBM+Plex+Mono:wght@400;600'
        '&family=IBM+Plex+Sans:wght@400;600&display=swap">',
        f"<style>{css}</style>",
        '<div class="wrap">',
        "<h1>Grain fit research (AR6)</h1>",
        f'<p class="sub">Fitting every voxel of a grain together against the detector images. '
        f'Branch <code>feature/grain-opt</code>. Updated {esc(prog["updated"])}.</p>',
        '<div class="tiles">'
        + "".join(f'<div class="tile"><b class="{s}">{counts[s]}</b><span>{s}</span></div>' for s in STATUS_ORDER)
        + f'<div class="tile"><b>{sum(1 for s in prog["steps"] if s["status"] == "done")}'
        f'/{len(prog["steps"])}</b><span>steps done</span></div></div>',
        "<h2>Now</h2>",
        '<div class="now">' + "".join(f"<p>{md_inline(ln.lstrip('- '))}</p>" for ln in now) + "</div>",
        "<h2>Goals</h2>",
        '<ul class="goals">' + "".join(f"<li><code>{k}</code> {esc(v)}</li>" for k, v in GOALS.items()) + "</ul>",
        "<h2>Steps</h2>",
        '<div class="tbl"><table><thead><tr><th>step</th><th>what</th><th>status</th><th>result</th>'
        "</tr></thead><tbody>",
    ]
    for s in reversed(prog["steps"]):
        parts.append(
            f'<tr><td class="id">{esc(s["id"])}</td><td>{esc(s["title"])}<br><span class="note">'
            f'{esc(s["date"])} {esc(s.get("commit", ""))}</span></td>'
            f'<td><span class="pill {esc(s["status"])}">{esc(s["status"])}</span></td>'
            f'<td class="res">{esc(s["result"])}</td></tr>'
        )
    parts.append("</tbody></table></div>")
    parts.append("<h2>Cost against quality</h2>")
    parts.append(f'<p class="note">{esc(prog["pareto_note"])}</p>')
    parts.append('<p class="note"><strong>Unit.</strong> One objective evaluation is one call of the '
                 'objective being fitted. In a grain fit that call scores every voxel of the grain '
                 'once, so N grain evaluations count as N per voxel. In a per-voxel fit it is one '
                 'strict-cost call for one voxel. Chained stages add up. Log scale; 0 = no fit.</p>')
    parts.append(pareto_svg(prog["pareto"]))
    key_rows = "".join(
        f'<tr><td class="id">{i}</td><td>{esc(p["label"])}</td><td class="id">{esc(p.get("set", ""))}</td>'
        f'<td class="id">{p["renders_per_voxel"]}</td><td class="id">{p["k"]} of {p["n"]} '
        f'({wilson(p["k"], p["n"])[1]:.1f} to {wilson(p["k"], p["n"])[2]:.1f}%)</td>'
        f'<td class="id">{(f"{p[chr(109)+chr(101)+chr(100)+chr(105)+chr(97)+chr(110)]:.4f}" if "median" in p else "")}</td></tr>'
        for i, p in enumerate(prog["pareto"], 1)
    )
    parts.append('<div class="tbl" style="margin-top:10px"><table><thead><tr><th>#</th><th>arm</th>'
                 '<th>voxels</th><th>num objective eval / voxel</th><th>above 0.1 deg</th><th>median deg</th>'
                 f'</tr></thead><tbody>{key_rows}</tbody></table></div>')
    parts.append("<h2>Ideas</h2>")
    parts.append('<div class="filters" role="group" aria-label="Filter ideas by status">'
                 '<button id="f-all" data-f="all" aria-pressed="true">all</button>'
                 + "".join(f'<button id="f-{s}" data-f="{s}" aria-pressed="false">{s}</button>' for s in STATUS_ORDER)
                 + "</div>")
    parts.append('<div class="tbl"><table id="ideas"><thead><tr><th>id</th><th>goal</th><th>idea</th>'
                 "<th>status</th><th>result</th></tr></thead><tbody>")
    key = {s: i for i, s in enumerate(STATUS_ORDER)}
    for r in sorted(rows, key=lambda r: (key.get(r["status"], 9), int(r["id"][1:]))):
        parts.append(
            f'<tr data-s="{esc(r["status"])}"><td class="id">{esc(r["id"])}</td><td class="id">{esc(r["goal"])}</td>'
            f'<td>{md_inline(r["idea"])}</td><td><span class="pill {esc(r["status"])}">{esc(r["status"])}</span></td>'
            f'<td class="res">{md_inline(r["result"])}</td></tr>'
        )
    parts.append("</tbody></table></div>")
    parts.append(
        '<p class="note">Source files: <code>icenine_py/docs/grain_opt/IDEAS.md</code>, <code>STATE.md</code>, '
        "<code>findings.md</code>, <code>benchmarks/grain_opt/progress.json</code>.</p></div>"
    )
    parts.append(
        "<script>document.querySelectorAll('.filters button').forEach(b=>b.addEventListener('click',()=>{"
        "const f=b.dataset.f;document.querySelectorAll('.filters button').forEach(x=>x.setAttribute('aria-pressed',x===b));"
        "document.querySelectorAll('#ideas tbody tr').forEach(r=>{r.hidden=!(f==='all'||r.dataset.s===f)})}));</script>"
    )
    return "\n".join(parts)


def main() -> None:
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(build())
    print(OUT)


if __name__ == "__main__":
    main()
