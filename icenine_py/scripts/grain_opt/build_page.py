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


def pareto_svg(points: List[Dict[str, Any]]) -> str:
    W, H, L, R, T, B = 640, 300, 56, 16, 14, 44
    xs = [p["renders_per_voxel"] for p in points]
    xmax = max(max(xs), 100)
    ys = [wilson(p["k"], p["n"]) for p in points]
    ymax = max(5.0, max(y[2] for y in ys) * 1.15)

    def X(v: float) -> float:
        # log scale with 0 drawn at the left edge ("0" tick)
        lv = math.log10(max(v, 1.0))
        return L + (W - L - R) * (0.06 + 0.94 * lv / math.log10(xmax * 1.2)) if v > 0 else L + 4

    def Y(v: float) -> float:
        return T + (H - T - B) * (1 - v / ymax)

    out = [f'<svg viewBox="0 0 {W} {H}" role="img" aria-label="Pareto chart" class="chart">']
    step = 1 if ymax <= 6 else (2 if ymax <= 12 else 5)
    v = 0.0
    while v <= ymax:
        out.append(f'<line x1="{L}" x2="{W-R}" y1="{Y(v):.1f}" y2="{Y(v):.1f}" class="grid"/>')
        out.append(f'<text x="{L-8}" y="{Y(v)+4:.1f}" class="tick" text-anchor="end">{v:g}%</text>')
        v += step
    for t in [0, 1, 10, 100, 1000]:
        if t > xmax * 1.2:
            break
        out.append(f'<text x="{X(t):.1f}" y="{H-B+18}" class="tick" text-anchor="middle">{t}</text>')
    out.append(f'<text x="{(L+W-R)/2:.0f}" y="{H-6}" class="axis" text-anchor="middle">'
               "renders per voxel (log)</text>")
    out.append(f'<text x="14" y="{(T+H-B)/2:.0f}" class="axis" text-anchor="middle" '
               f'transform="rotate(-90 14 {(T+H-B)/2:.0f})">voxels above 0.1 deg</text>')
    for p, (c, lo, hi) in zip(points, ys):
        fam = esc(p.get("family", "g20"))
        x = X(p["renders_per_voxel"])
        out.append(f'<line x1="{x:.1f}" x2="{x:.1f}" y1="{Y(lo):.1f}" y2="{Y(hi):.1f}" '
                   f'class="whisk {fam}"/>')
        out.append(f'<circle cx="{x:.1f}" cy="{Y(c):.1f}" r="4.5" class="dot {fam}"><title>'
                   f'{esc(p["label"])}: {c:.2f}% ({lo:.2f} to {hi:.2f}), '
                   f'{p["renders_per_voxel"]} renders/voxel, {esc(p.get("set", ""))}</title></circle>')
        out.append(f'<text x="{x+8:.1f}" y="{Y(c)-6:.1f}" class="lab">{esc(p["label"])}</text>')
    out.append("</svg>")
    return "".join(out)


def build() -> str:
    prog = json.loads((ICE / "benchmarks" / "grain_opt" / "progress.json").read_text())
    rows = ideas()
    counts = {s: sum(1 for r in rows if r["status"] == s) for s in STATUS_ORDER}
    now = state_now()
    css = """
:root{--bg:#f6f7f9;--panel:#ffffff;--fg:#1c2330;--muted:#5d6879;--line:#dde2ea;
--accent:#2f6f9f;--run:#b26b00;--done:#2e7d4f;--drop:#8a8f98;--prop:#4a5a73;
--ref:#8a8f98;--sans:"IBM Plex Sans",system-ui,sans-serif;--mono:"IBM Plex Mono",ui-monospace,monospace}
@media (prefers-color-scheme:dark){:root:not([data-theme="light"]){--bg:#11151b;--panel:#171c24;
--fg:#e3e8ef;--muted:#9aa5b5;--line:#2a313c;--accent:#6fb0de;--run:#e3a33f;--done:#5cbf86;
--drop:#7c838d;--prop:#a9b6cb;--ref:#7c838d;color-scheme:dark}}
:root[data-theme="dark"]{--bg:#11151b;--panel:#171c24;--fg:#e3e8ef;--muted:#9aa5b5;--line:#2a313c;
--accent:#6fb0de;--run:#e3a33f;--done:#5cbf86;--drop:#7c838d;--prop:#a9b6cb;--ref:#7c838d;
color-scheme:dark}
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
.chart .dot{fill:var(--accent)}.chart .dot.reference{fill:var(--ref)}
.chart .whisk{stroke:var(--accent);stroke-width:2}.chart .whisk.reference{stroke:var(--ref)}
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
    parts.append(f'<p class="note">{esc(prog["pareto_note"])} Whiskers: Wilson 95% intervals.</p>')
    parts.append(pareto_svg(prog["pareto"]))
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
