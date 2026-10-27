"""3D superposition viewer: a small HTML page that uses the 3Dmol.js library.

By default the browser loads 3Dmol.js from the jsDelivr CDN. If a local copy exists
in the cache (``vendor_3dmol()`` fetches it from the npm registry) it is inlined
instead, which makes the viewer work without CDN access (e.g. for screenshots).
"""

from __future__ import annotations

import io
import json
import tarfile
import urllib.request
from pathlib import Path

import numpy as np

from proteinlens.align import AlignmentResult
from proteinlens.data import Structure, cache_dir

THREEDMOL_VERSION = "2.5.5"
CDN_URL = f"https://cdn.jsdelivr.net/npm/3dmol@{THREEDMOL_VERSION}/build/3Dmol-min.js"
NPM_URL = f"https://registry.npmjs.org/3dmol/-/3dmol-{THREEDMOL_VERSION}.tgz"
COLOR_A, COLOR_B = "#2a78d6", "#eb6834"  # blue = A (moved onto B), orange = B


def _local_js() -> Path:
    return cache_dir() / "3dmol" / f"3Dmol-{THREEDMOL_VERSION}-min.js"


def vendor_3dmol() -> Path:
    """Download 3Dmol-min.js from the npm registry into the cache (once)."""
    dest = _local_js()
    if dest.exists():
        return dest
    with urllib.request.urlopen(NPM_URL, timeout=120) as resp:
        data = resp.read()
    with tarfile.open(fileobj=io.BytesIO(data), mode="r:gz") as tar:
        member = tar.getmember("package/build/3Dmol-min.js")
        js = tar.extractfile(member).read()  # type: ignore[union-attr]
    dest.parent.mkdir(parents=True, exist_ok=True)
    dest.write_bytes(js)
    return dest


def _script_tag() -> str:
    local = _local_js()
    if local.exists():
        return f"<script>{local.read_text()}</script>"
    return f'<script src="{CDN_URL}"></script>'


def superposition_html(
    a: Structure, b: Structure, result: AlignmentResult, height: int = 460
) -> str:
    """HTML showing A (blue, moved onto B with the TM-align superposition) and B (orange).
    Residues of A in the top deviating regions are drawn thicker in red."""
    pdb_a = a.to_pdb(np.asarray(result.rotation), np.asarray(result.translation), chain_id="A")
    pdb_b = b.to_pdb(chain_id="B")
    highlight = []
    for g in result.top_regions:
        lo, hi = g.residues_a.split("-") if "-" in g.residues_a[1:] else (g.residues_a,) * 2
        try:
            highlight.append([int(lo), int(hi)])
        except ValueError:  # insertion codes; skip highlighting
            continue
    return f"""<div id="pl-viewer" style="width:100%;height:{height}px;position:relative;"></div>
{_script_tag()}
<script>
const viewer = $3Dmol.createViewer("pl-viewer", {{backgroundColor: "white"}});
viewer.addModel({json.dumps(pdb_a)}, "pdb");
viewer.addModel({json.dumps(pdb_b)}, "pdb");
viewer.setStyle({{model: 0}}, {{cartoon: {{color: "{COLOR_A}"}}}});
viewer.setStyle({{model: 1}}, {{cartoon: {{color: "{COLOR_B}", opacity: 0.85}}}});
for (const [lo, hi] of {json.dumps(highlight)}) {{
  viewer.setStyle({{model: 0, resi: lo + "-" + hi}},
                  {{cartoon: {{color: "#c0392b", thickness: 1.0}}}});
}}
viewer.zoomTo();
viewer.render();
</script>"""


def single_html(s: Structure, height: int = 380) -> str:
    """HTML showing one structure colored from N- to C-terminus."""
    return f"""<div id="pl-single" style="width:100%;height:{height}px;position:relative;"></div>
{_script_tag()}
<script>
const v = $3Dmol.createViewer("pl-single", {{backgroundColor: "white"}});
v.addModel({json.dumps(s.to_pdb())}, "pdb");
v.setStyle({{}}, {{cartoon: {{color: "spectrum"}}}});
v.zoomTo();
v.render();
</script>"""
