"""Put the cost-against-accuracy section into the review page.

usage: pareto_insert.py <live html in> <html out>

Run pareto_build.py first (it writes pareto_parts.json next to this file). On a page that
already has the section, only the section is replaced. The texts and captions are in
pareto_section.html, the styles in pareto_page.css: after new rows, check every number in the
texts against the frontiers that pareto_build.py prints.
"""

import json
import sys
from pathlib import Path

HERE = Path(__file__).parent
PARTS = json.loads((HERE / "pareto_parts.json").read_text())
CSS = (HERE / "pareto_page.css").read_text()


def swatch(cls: str, diamond: bool = False) -> str:
    mark = (
        f'<path class="pt {cls}" d="M9 1.8 L16.2 9 L9 16.2 L1.8 9 Z"/>'
        if diamond
        else f'<circle class="pt {cls}" cx="9" cy="9" r="5.8"/>'
    )
    return f'<svg viewBox="0 0 18 18" aria-hidden="true">{mark}</svg>'


LEGEND = (
    '<div class="plegend" aria-label="Legend for the six plots">'
    f'<span>{swatch("pt-grid")}grid search</span>'
    f'<span>{swatch("pt-ray")}ray indexing, one voxel at a time</span>'
    f'<span>{swatch("pt-bfs")}BFS with indexed seeds</span>'
    f'<span>{swatch("pt-paint")}paint grains</span>'
    f'<span>{swatch("pt-bfs", True)}diamond: with a boundary pass</span>'
    '<span><svg viewBox="0 0 26 18" aria-hidden="true" style="width:26px">'
    '<path class="front" d="M1 4 H9 V13 H25"/><circle class="fring" cx="9" cy="13" r="4"/></svg>'
    "line and ring: on the frontier</span>"
    "</div>"
)


def section() -> str:
    s = (HERE / "pareto_section.html").read_text()
    s = s.replace("{{LEGEND}}", LEGEND).replace("{{TABLE}}", PARTS["table"])
    for sid, svg in PARTS["svgs"].items():
        assert s.count("{{SVG:%s}}" % sid) == 1, sid
        s = s.replace("{{SVG:%s}}" % sid, svg)
    assert "{{" not in s
    return s


def rep(s: str, old: str, new: str, n: int = 1) -> str:
    assert s.count(old) == n, (old[:70], s.count(old))
    return s.replace(old, new)


SUMMARY_BULLET = (
    "    <li><strong>Cost against accuracy (new, 11 Nov):</strong> every method tried so far is "
    "on one plot per sample, with the frontier marked. On the gradient sample the best trade is "
    "the BFS with indexed seeds, the CMA local optimization and one boundary pass: about 3 times "
    "fewer evaluations than the base run and fewer wrong voxels. See "
    '<a href="#frontier">Cost against accuracy</a>.</li>\n'
)
UPDATED = (
    'Last updated 2026-11-11: a new section, "Cost against accuracy, all methods", puts every '
    "method tried so far on a plot of cost against accuracy, one per sample (your request). "
    "Before that, 2026-11-10 at 23:50: "
)


def apply(html: str) -> str:
    found = '<section id="found">\n'
    if 'id="frontier"' in html:  # an update: swap the section, keep the one-time edits
        a, b = html.index('<section id="frontier">'), html.index(found)
        return html[:a] + section() + html[b:]
    html = rep(
        html, "  --series: #2a78d6;\n", "  --series: #2a78d6; --s2: #eb6834; --s3: #1baf7a;\n"
    )
    html = rep(html, "--series: #3987e5;", "--series: #3987e5; --s2: #d95926; --s3: #199e70;", 2)
    html = rep(html, "footer { font-size: 13.5px;", CSS + "footer { font-size: 13.5px;")
    nav = '    <a href="#found">What was found</a>\n'
    html = rep(html, nav, '    <a href="#frontier">Cost against accuracy</a>\n' + nav)
    pr42 = "    <li><strong>PR #42</strong> is merged on your word.</li>\n"
    html = rep(html, pr42, SUMMARY_BULLET + pr42)
    html = rep(html, "as of 2026-11-10</div>", "as of 2026-11-11</div>")
    html = rep(html, "Last updated 2026-11-10 at 23:50: ", UPDATED)
    return rep(html, found, section() + found)


if __name__ == "__main__":
    out = apply(Path(sys.argv[1]).read_text())
    Path(sys.argv[2]).write_text(out)
    print("bytes", len(out.encode()), "figures", out.count("<figure>"), "svg", out.count("<svg"))
    for w in ("finisher", "local fit"):
        print(w, out.count(w))
