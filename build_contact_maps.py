"""
Extract contact maps from scatterplot JSON (*_common structures).

Input:
  P01732,P10966,P01730,1RHH,1RHH_scatterplot.json

Each structure*_common has:
  - nodes1: residues (often annotated with IgStrand, e.g. R25.A.P01732=>A1549)
  - links:  contacts (source/target residue pairs)

Outputs:
  - contact_maps.json   (matrices + edge lists)
  - contact_maps.html   (interactive grid visualization)
"""

from __future__ import annotations

import argparse
import json
import re
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parent
DEFAULT_JSON = ROOT / "P01732,P10966,P01730,1RHH,1RHH_scatterplot.json"
OUT_JSON = ROOT / "contact_maps.json"
OUT_HTML = ROOT / "contact_maps.html"

# Residues may include insertion codes: L82C, N35A, T82A
ID_PAT = re.compile(
    r"^([A-Z])(\d+[A-Za-z]?)\.([A-Za-z0-9]+)\.([A-Za-z0-9]+)(?:={1,2}>(.+))?$"
)
BLOCK_PAT = re.compile(r'(?:^|[\n,])\s*"(structure\d+(?:_common|_diff)?)"\s*:\s*\{')


def extract_object(text: str, start: int) -> str:
    depth = 0
    i = start
    in_str = False
    esc = False
    while i < len(text):
        ch = text[i]
        if in_str:
            if esc:
                esc = False
            elif ch == "\\":
                esc = True
            elif ch == '"':
                in_str = False
        else:
            if ch == '"':
                in_str = True
            elif ch == "{":
                depth += 1
            elif ch == "}":
                depth -= 1
                if depth == 0:
                    return text[start : i + 1]
        i += 1
    raise ValueError("unbalanced braces")


def fix_json_blob(s: str) -> str:
    # trailing commas
    s = re.sub(r",\s*([}\]])", r"\1", s)
    return s


def load_structures(path: Path) -> dict[str, dict]:
    raw = path.read_text(encoding="utf-8", errors="replace")
    # Fix glued objects: }}"structure0_common" → }, "structure0_common"
    raw = re.sub(r'\}\s*"structure', '}, "structure', raw)

    structs: dict[str, dict] = {}
    matches = list(BLOCK_PAT.finditer(raw))
    # also catch structure0_common if glued at start of line without newline pattern
    if '"structure0_common"' in raw and not any(m.group(1) == "structure0_common" for m in matches):
        # re-scan after fix
        matches = list(BLOCK_PAT.finditer(raw))

    for m in matches:
        name = m.group(1)
        brace = raw.find("{", m.end() - 1)
        blob = extract_object(raw, brace)
        structs[name] = json.loads(fix_json_blob(blob))
    return structs


def parse_node_id(s: str) -> dict | None:
    m = ID_PAT.match(s)
    if not m:
        return None
    aa, pos_raw, chain, protein, ig = m.groups()
    # numeric part of position (ignore insertion letter for matrix index)
    num = int(re.match(r"\d+", pos_raw).group(0))
    return {
        "aa": aa,
        "pos_label": pos_raw,
        "pos": num,
        "chain": chain,
        "protein": protein,
        "igstrand": ig,
        "raw": s,
    }


def residue_sort_key(label: str) -> tuple:
    m = re.match(r"(\d+)([A-Za-z]?)$", label)
    if not m:
        return (10**9, label)
    return (int(m.group(1)), m.group(2) or "")


def build_contact_map(struct: dict, name: str) -> dict:
    protein = struct.get("id", name)
    nodes = []
    by_id = {}
    for n in struct.get("nodes1", []):
        p = parse_node_id(n["id"])
        if not p:
            continue
        p["color"] = n.get("c")
        p["x"] = n.get("x")
        p["y"] = n.get("y")
        nodes.append(p)
        by_id[n["id"]] = p

    # ordered unique residue labels on this structure
    labels = sorted({p["pos_label"] for p in nodes}, key=residue_sort_key)
    index = {lab: i for i, lab in enumerate(labels)}
    n = len(labels)
    matrix = [[0 for _ in range(n)] for _ in range(n)]

    contacts = []
    ig_contacts = []
    seen = set()

    for link in struct.get("links", []):
        a = parse_node_id(link["source"])
        b = parse_node_id(link["target"])
        if not a or not b:
            continue
        weight = link.get("n") or link.get("v") or 1
        ia, ib = index.get(a["pos_label"]), index.get(b["pos_label"])
        if ia is None or ib is None:
            # node may be missing from nodes1; extend labels dynamically
            for p in (a, b):
                if p["pos_label"] not in index:
                    index[p["pos_label"]] = len(labels)
                    labels.append(p["pos_label"])
            # rebuild sizes if extended
            if len(labels) != n:
                n = len(labels)
                new_mat = [[0 for _ in range(n)] for _ in range(n)]
                for r in range(len(matrix)):
                    for c in range(len(matrix[r])):
                        new_mat[r][c] = matrix[r][c]
                matrix = new_mat
                ia, ib = index[a["pos_label"]], index[b["pos_label"]]

        matrix[ia][ib] = max(matrix[ia][ib], int(weight))
        matrix[ib][ia] = max(matrix[ib][ia], int(weight))

        key = tuple(sorted([a["pos_label"], b["pos_label"]]))
        if key not in seen and a["pos_label"] != b["pos_label"]:
            seen.add(key)
            contacts.append(
                {
                    "a": {
                        "aa": a["aa"],
                        "pos": a["pos_label"],
                        "igstrand": a["igstrand"],
                    },
                    "b": {
                        "aa": b["aa"],
                        "pos": b["pos_label"],
                        "igstrand": b["igstrand"],
                    },
                    "weight": int(weight),
                }
            )
            if a["igstrand"] and b["igstrand"]:
                ig_contacts.append(
                    {
                        "a": a["igstrand"],
                        "b": b["igstrand"],
                        "weight": int(weight),
                        "res_a": f"{a['aa']}{a['pos_label']}",
                        "res_b": f"{b['aa']}{b['pos_label']}",
                    }
                )

    return {
        "structure": name,
        "protein": protein,
        "n_residues": len(labels),
        "n_contacts": len(contacts),
        "residues": labels,
        "residue_meta": {
            p["pos_label"]: {
                "aa": p["aa"],
                "igstrand": p["igstrand"],
                "chain": p["chain"],
            }
            for p in nodes
        },
        "matrix": matrix,
        "contacts": contacts,
        "igstrand_contacts": ig_contacts,
    }


def consensus_igstrand_map(maps: list[dict]) -> dict:
    """Contacts shared in IgStrand space across common structures."""
    counts: dict[tuple[str, str], int] = defaultdict(int)
    examples: dict[tuple[str, str], list] = defaultdict(list)
    for m in maps:
        seen = set()
        for c in m["igstrand_contacts"]:
            key = tuple(sorted([c["a"], c["b"]]))
            if key in seen:
                continue
            seen.add(key)
            counts[key] += 1
            examples[key].append(
                {
                    "protein": m["protein"],
                    "res_a": c["res_a"],
                    "res_b": c["res_b"],
                    "weight": c["weight"],
                }
            )

    ig_labels = sorted({x for pair in counts for x in pair})
    index = {lab: i for i, lab in enumerate(ig_labels)}
    n = len(ig_labels)
    matrix = [[0 for _ in range(n)] for _ in range(n)]
    edges = []
    for (a, b), cnt in sorted(counts.items(), key=lambda kv: (-kv[1], kv[0])):
        ia, ib = index[a], index[b]
        matrix[ia][ib] = cnt
        matrix[ib][ia] = cnt
        edges.append(
            {
                "a": a,
                "b": b,
                "n_proteins": cnt,
                "examples": examples[(a, b)],
            }
        )

    return {
        "coordinate_system": "IgStrand",
        "n_positions": n,
        "positions": ig_labels,
        "matrix": matrix,  # cell = how many proteins share this contact
        "contacts": edges,
        "note": "matrix[i][j] = number of *_common structures that contain this IgStrand contact",
    }


HTML_TEMPLATE = r"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>Common contact maps</title>
<style>
:root{--ink:#162033;--muted:#64748b;--line:#dbe3ee;--bg:#eef2f7}
*{box-sizing:border-box}
body{margin:0;font:14px/1.45 Helvetica,Arial,sans-serif;color:var(--ink);background:var(--bg)}
header{padding:16px 20px;background:#fff;border-bottom:1px solid var(--line)}
header h1{margin:0;font-size:18px}
header p{margin:4px 0 0;color:var(--muted);font-size:13px}
main{display:grid;grid-template-columns:240px 1fr 300px;gap:12px;padding:12px;height:calc(100vh - 78px)}
.card{background:#fff;border:1px solid var(--line);border-radius:10px;padding:12px;overflow:auto}
nav button{display:block;width:100%;text-align:left;margin:0 0 6px;padding:8px 10px;border:1px solid #cbd5e1;border-radius:7px;background:#fff;cursor:pointer;font-size:12px}
nav button.active{background:#172554;color:#fff;border-color:#172554}
#wrap{overflow:auto}
canvas{image-rendering:pixelated;cursor:crosshair;background:#111}
#legend{font-size:12px;margin-top:8px}
#info h2{margin:0 0 8px;font-size:15px}
#info pre{white-space:pre-wrap;font-size:11px;background:#f8fafc;padding:8px;border-radius:6px;border:1px solid var(--line)}
@media(max-width:1000px){main{grid-template-columns:1fr;height:auto}}
</style>
</head>
<body>
<header>
  <h1>Contact maps from common interactions</h1>
  <p>Residue×residue grids from <code>*_common</code> links · dark cell = contact · consensus uses IgStrand IDs</p>
</header>
<main>
  <nav class="card" id="nav"></nav>
  <section class="card">
    <div id="title" style="font-weight:700;margin-bottom:8px"></div>
    <div id="wrap"><canvas id="cv"></canvas></div>
    <div id="legend"></div>
  </section>
  <aside class="card" id="info"><h2>Inspector</h2><p style="color:#888">hover a cell</p></aside>
</main>
<script>
const DATA = __DATA__;
const maps = DATA.maps;
const cons = DATA.consensus_igstrand;

function ramp(t){
  // black → green (contact strength)
  const r = Math.round(10+20*(1-t));
  const g = Math.round(30+200*t);
  const b = Math.round(20+40*(1-t));
  return `rgb(${r},${g},${b})`;
}

function drawMap(map, mode){
  const labels = map.residues || map.positions;
  const M = map.matrix;
  const n = labels.length;
  const cell = Math.max(4, Math.min(14, Math.floor(700/Math.max(n,1))));
  const pad = 48;
  const W = pad + n*cell + 8;
  const H = pad + n*cell + 8;
  const cv = document.getElementById('cv');
  cv.width = W; cv.height = H;
  const ctx = cv.getContext('2d');
  ctx.fillStyle = '#0b1220';
  ctx.fillRect(0,0,W,H);

  let maxv = 1;
  for(let i=0;i<n;i++) for(let j=0;j<n;j++) maxv = Math.max(maxv, M[i][j]||0);

  for(let i=0;i<n;i++){
    for(let j=0;j<n;j++){
      const v = M[i][j]||0;
      if(v<=0){
        ctx.fillStyle = '#1e293b';
      } else {
        ctx.fillStyle = ramp(v/maxv);
      }
      ctx.fillRect(pad+j*cell, pad+i*cell, cell-1, cell-1);
    }
  }

  // sparse axis labels
  ctx.fillStyle = '#94a3b8';
  ctx.font = '9px monospace';
  const step = Math.max(1, Math.ceil(n/25));
  for(let i=0;i<n;i+=step){
    const lab = String(labels[i]);
    ctx.save();
    ctx.translate(pad+i*cell+cell/2, pad-4);
    ctx.rotate(-Math.PI/2);
    ctx.fillText(lab, 0, 0);
    ctx.restore();
    ctx.fillText(lab, 2, pad+i*cell+cell*0.75);
  }

  document.getElementById('title').textContent =
    (map.protein? map.protein+' · ' : '') + (map.structure||'consensus') +
    ` · ${n}×${n} · ${map.n_contacts|| (map.contacts&&map.contacts.length) || 0} contacts`;

  document.getElementById('legend').innerHTML =
    mode==='consensus'
      ? `IgStrand consensus: cell = # of proteins sharing contact (max ${maxv}). Green = more conserved contact.`
      : `Residue contact map: cell = link weight (max ${maxv}). Green = contact present.`;

  cv.onmousemove = (ev)=>{
    const rect = cv.getBoundingClientRect();
    const scaleX = cv.width/rect.width, scaleY = cv.height/rect.height;
    const x = (ev.clientX-rect.left)*scaleX;
    const y = (ev.clientY-rect.top)*scaleY;
    const j = Math.floor((x-pad)/cell);
    const i = Math.floor((y-pad)/cell);
    const info = document.getElementById('info');
    if(i<0||j<0||i>=n||j>=n){ info.innerHTML='<h2>Inspector</h2><p style="color:#888">hover a cell</p>'; return; }
    const v = M[i][j]||0;
    const ri = labels[i], rj = labels[j];
    let extra='';
    if(mode==='protein'){
      const meta = map.residue_meta||{};
      const A = meta[ri]||{}, B = meta[rj]||{};
      extra = `\n${A.aa||''}${ri}  (IgStrand ${A.igstrand||'—'})`+
              `\n${B.aa||''}${rj}  (IgStrand ${B.igstrand||'—'})`;
      const hit = (map.contacts||[]).find(c =>
        (c.a.pos===ri && c.b.pos===rj) || (c.a.pos===rj && c.b.pos===ri));
      if(hit) extra += `\nweight ${hit.weight}`;
    } else {
      const hit = (map.contacts||[]).find(c =>
        (c.a===ri && c.b===rj) || (c.a===rj && c.b===ri));
      if(hit){
        extra = `\nshared by ${hit.n_proteins} protein(s)\n`+
          hit.examples.map(e=>`${e.protein}: ${e.res_a}–${e.res_b}`).join('\n');
      }
    }
    info.innerHTML = `<h2>Inspector</h2><pre>${ri} × ${rj}\nvalue: ${v}${extra}</pre>`;
  };
}

const nav = document.getElementById('nav');
const items = [
  {id:'consensus', label:'Consensus (IgStrand)', mode:'consensus', map:cons},
  ...maps.map((m,i)=>({id:'m'+i, label:`${m.protein} (${m.structure})`, mode:'protein', map:m}))
];
items.forEach((it,idx)=>{
  const b=document.createElement('button');
  b.textContent=it.label;
  if(idx===0) b.classList.add('active');
  b.onclick=()=>{
    [...nav.children].forEach(x=>x.classList.remove('active'));
    b.classList.add('active');
    drawMap(it.map, it.mode);
  };
  nav.appendChild(b);
});
drawMap(items[0].map, items[0].mode);
</script>
</body>
</html>
"""


def main() -> None:
    ap = argparse.ArgumentParser(description="Build contact maps from *_common structures")
    ap.add_argument("json_path", nargs="?", type=Path, default=DEFAULT_JSON)
    ap.add_argument("-o", "--output", type=Path, default=OUT_JSON)
    ap.add_argument("--html", type=Path, default=OUT_HTML)
    args = ap.parse_args()

    structs = load_structures(args.json_path.expanduser().resolve())

    # Prefer *_common; if structure0_common missing, still use 1-4
    common_names = sorted(
        [k for k in structs if k.endswith("_common")],
        key=lambda s: int(re.search(r"\d+", s).group(0)),
    )
    if not common_names:
        raise SystemExit("No *_common structures found in JSON")

    maps = [build_contact_map(structs[name], name) for name in common_names]
    consensus = consensus_igstrand_map(maps)

    doc = {
        "schema_version": "0.1",
        "type": "contact_maps",
        "source": str(args.json_path),
        "structures_used": common_names,
        "proteins": [m["protein"] for m in maps],
        "maps": maps,
        "consensus_igstrand": consensus,
        "note": (
            "Each map is built from links inside structure*_common. "
            "consensus_igstrand counts how many proteins share the same IgStrand–IgStrand contact."
        ),
    }

    args.output.write_text(json.dumps(doc, indent=2))
    html = HTML_TEMPLATE.replace("__DATA__", json.dumps(doc, separators=(",", ":")))
    args.html.write_text(html)

    print(f"wrote {args.output}")
    print(f"wrote {args.html}")
    for m in maps:
        print(f"  {m['structure']} ({m['protein']}): {m['n_residues']} residues, {m['n_contacts']} contacts")
    print(
        f"  consensus IgStrand: {consensus['n_positions']} positions, "
        f"{len(consensus['contacts'])} unique contacts"
    )


if __name__ == "__main__":
    main()
