"""Build an artifact-ready HTML page of orientation and misorientation maps from a JSON spec.

  uv run python scripts/mic_maps_report.py spec.json --out-dir DIR

Writes DIR/index.html, DIR/maps/<figure>.png and DIR/stats.json. Publish DIR/index.html as an
Artifact with every DIR/maps/*.png in `files` (published path "maps/<name>.png").

Spec (paths relative to the spec file unless absolute):
{
  "title": "BFS Reseed Maps",                      # 2-4 word page name
  "summary": "One paragraph: what is compared and the headline.",
  "reference": {"mic": "truth.mic", "index": "region_index.npy", "label": "Original (truth)"},
  "z": null,                                         # layer z, required for multi-layer files
  "figures": [
    {"name": "mc_clean", "title": "MC, clean images", "caption": "...",
     "reference": {...optional per-figure override...},
     "mics": [{"path": "a/recon.mic", "label": "No reseed"}, ...]}
  ],
  "notes": ["Optional caveats, one sentence each."]
}
Figures are drawn by plot_mic.plot_maps (same colour scales in every figure).
"""

import argparse
import html
import json
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))

import plot_mic  # noqa: E402

CSS = """
:root {
  --bg: #f5f7f9; --surface: #ffffff; --fg: #1a2028; --muted: #5a6470; --rule: #d8dde3;
  --accent: #2c5d8a; --accent-soft: #e6eef6;
  --sans: "IBM Plex Sans", -apple-system, "Segoe UI", Roboto, sans-serif;
  --serif: "IBM Plex Serif", Georgia, serif;
  --mono: "IBM Plex Mono", ui-monospace, Menlo, Consolas, monospace;
}
@media (prefers-color-scheme: dark) {
  :root:not([data-theme="light"]) {
    --bg: #12161b; --surface: #1a2027; --fg: #e3e7ec; --muted: #9aa4af; --rule: #2e3742;
    --accent: #8db8e0; --accent-soft: #1d2a37; color-scheme: dark;
  }
}
:root[data-theme="dark"] {
  --bg: #12161b; --surface: #1a2027; --fg: #e3e7ec; --muted: #9aa4af; --rule: #2e3742;
  --accent: #8db8e0; --accent-soft: #1d2a37; color-scheme: dark;
}
body { background: var(--bg); color: var(--fg); font: 16px/1.6 var(--sans); padding: 0 16px; }
main { max-width: 72rem; margin: 0 auto; padding-block: 40px 72px; display: flex;
  flex-direction: column; gap: 10px; }
.prose { max-width: 46rem; }
h1, h2 { font-family: var(--serif); text-wrap: balance; line-height: 1.25; margin: 0; }
h1 { font-size: 1.9rem; font-weight: 600; }
h2 { font-size: 1.3rem; font-weight: 600; margin-top: 28px; padding-top: 14px;
  border-top: 1px solid var(--rule); }
p, ul { margin: 4px 0; }
.meta { color: var(--muted); font: 0.85rem var(--mono); }
code { font-family: var(--mono); font-size: 0.86em; background: var(--accent-soft);
  padding: 1px 5px; border-radius: 3px; overflow-wrap: anywhere; }
.table-wrap { overflow-x: auto; background: var(--surface); border: 1px solid var(--rule);
  border-radius: 4px; margin: 8px 0; }
table { border-collapse: collapse; width: 100%; font-size: 0.9rem; }
th, td { padding: 7px 12px; text-align: left; border-bottom: 1px solid var(--rule); }
th { color: var(--muted); font-size: 0.76rem; letter-spacing: 0.04em; text-transform: uppercase;
  font-weight: 600; }
tr:last-child td { border-bottom: 0; }
td.n, th.n { text-align: right; font-family: var(--mono); font-variant-numeric: tabular-nums;
  white-space: nowrap; }
figure { margin: 8px 0 0; overflow-x: auto; }
figure img { display: block; max-width: 100%; height: auto; background: #fff; border-radius: 4px;
  border: 1px solid var(--rule); }
figcaption { color: var(--muted); font-size: 0.88rem; margin-top: 6px; max-width: 46rem; }
"""

FONTS = (
    "https://fonts.googleapis.com/css2?family=IBM+Plex+Mono:wght@400"
    "&family=IBM+Plex+Sans:wght@400;600&family=IBM+Plex+Serif:wght@600&display=swap"
)

HOW_TO_READ = (
    "Each figure has two rows. Top row: orientation maps, the reference first, then each "
    "reconstruction. Colour is the orientation reduced to the cubic fundamental zone, shown as its "
    "Rodrigues vector in red, green and blue on a fixed scale, so one orientation has the same "
    "colour in every panel. Bottom row: misorientation of each reconstruction against the "
    "reference, voxel by voxel, after cubic symmetry, from blue (0°) to red (62.8°, the largest "
    "misorientation cubic symmetry allows). Grey voxels have no matching voxel in that "
    "reconstruction. A wrong voxel is one more than 1° from the reference."
)


def _resolve(base: Path, p: Optional[str]) -> Optional[str]:
    if p is None:
        return None
    q = Path(p).expanduser()
    return str(q if q.is_absolute() else (base / q).resolve())


def _fmt(x: Optional[float], nd: int = 3) -> str:
    return "–" if x is None else f"{x:.{nd}f}"


def build(spec_path: Path, out_dir: Path) -> Dict[str, Any]:
    spec = json.loads(spec_path.read_text())
    base = spec_path.parent
    (out_dir / "maps").mkdir(parents=True, exist_ok=True)
    z = spec.get("z")
    all_stats: List[Dict[str, Any]] = []
    for fig in spec["figures"]:
        ref = {**spec.get("reference", {}), **fig.get("reference", {})}
        idx = np.load(_resolve(base, ref["index"])) if ref.get("index") else None
        png = out_dir / "maps" / f"{fig['name']}.png"
        stats = plot_mic.plot_maps(
            _resolve(base, ref["mic"]),
            [_resolve(base, m["path"]) for m in fig["mics"]],
            [m["label"] for m in fig["mics"]],
            str(png),
            idx,
            ref.get("label", "Reference"),
            fig.get("title", fig["name"]),
            z,
        )
        all_stats.append({"figure": fig["name"], "title": fig.get("title", ""), "panels": stats})
    (out_dir / "stats.json").write_text(json.dumps(all_stats, indent=2))

    e = html.escape
    rows = []
    for fs in all_stats:
        for st in fs["panels"]:
            rows.append(
                "<tr>"
                f"<td>{e(fs['title'] or fs['figure'])}</td><td>{e(st['label'])}</td>"
                f"<td class='n'>{st['n_matched']:,} of {st['n_reference']:,}</td>"
                f"<td class='n'>{st['n_over_1deg']:,}</td>"
                f"<td class='n'>{_fmt(st['median_right_deg'], 4)}</td>"
                f"<td class='n'>{_fmt(st['p90_deg'], 2)}</td>"
                f"<td class='n'>{_fmt(st['max_deg'], 1)}</td>"
                "</tr>"
            )
    figs = []
    for fig in spec["figures"]:
        cap = fig.get("caption", "")
        figs.append(
            f"<h2>{e(fig.get('title', fig['name']))}</h2>"
            f"<figure><img src='maps/{e(fig['name'])}.png' "
            f"alt='Orientation and misorientation maps: {e(fig.get('title', fig['name']))}'>"
            + (f"<figcaption>{e(cap)}</figcaption>" if cap else "")
            + "</figure>"
        )
    notes = spec.get("notes") or []
    ref0 = spec.get("reference", {})
    ref_name = Path(str(ref0.get("mic", ""))).name
    notes_html = (
        '<h2>Notes</h2><ul class="prose">' + "".join(f"<li>{e(n)}</li>" for n in notes) + "</ul>"
        if notes
        else ""
    )
    page = f"""<title>{e(spec['title'])}</title>
<link rel="preconnect" href="https://fonts.googleapis.com">
<link rel="stylesheet" href="{FONTS}">
<style>{CSS}</style>
<main>
  <h1>{e(spec['title'])}</h1>
  <p class="meta">Reference: {e(ref0.get('label', 'Reference'))} · {e(ref_name)}</p>
  <div class="prose">
    <p>{e(spec.get('summary', ''))}</p>
    <h2>How to read the maps</h2>
    <p>{HOW_TO_READ}</p>
  </div>
  <h2>Numbers</h2>
  <div class="table-wrap"><table>
    <thead><tr><th>Figure</th><th>Reconstruction</th><th class="n">Voxels compared</th>
    <th class="n">Wrong (over 1°)</th><th class="n">Median error of right voxels (°)</th>
    <th class="n">90th percentile (°)</th><th class="n">Max (°)</th></tr></thead>
    <tbody>{''.join(rows)}</tbody>
  </table></div>
  {''.join(figs)}
  {notes_html}
  <p class="meta">Drawn with icenine_py/scripts/mic_maps_report.py and plot_mic.py.</p>
</main>
"""
    (out_dir / "index.html").write_text(page, encoding="utf-8")
    return {"index": str(out_dir / "index.html"), "stats": all_stats}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    ap.add_argument("spec", help="JSON spec (see module docstring)")
    ap.add_argument("--out-dir", required=True)
    a = ap.parse_args()
    res = build(Path(a.spec).resolve(), Path(a.out_dir).resolve())
    print(res["index"])
    for fs in res["stats"]:
        for st in fs["panels"]:
            print(
                f"{fs['figure']} | {st['label']}: wrong {st['n_over_1deg']} of {st['n_matched']}"
                f" matched ({st['n_reference']} reference)"
            )


if __name__ == "__main__":
    main()
