# IgFolds2D Residue Interactions

2D Ig-fold maps and interaction layers for **CD8α** (UniProt [P01732](https://www.uniprot.org/uniprotkb/P01732), PDB **1CD8**), plus comparative contact maps across related Ig domains.

This repo turns spreadsheet geometry and iCn3D interaction tables into **layered JSON** and **static HTML** viewers you can open in a browser.

---

## Quick start

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install openpyxl
```

| What to look at | Open in browser |
|-----------------|-----------------|
| Layer 1 + H-bond density map | [`layer1_viz.html`](layer1_viz.html) |
| IgStrand strand colors only | [`layer1_viz_igstrand_colors.html`](layer1_viz_igstrand_colors.html) |
| Common residue–residue contacts | [`contact_maps.html`](contact_maps.html) |

No server required — double-click the HTML or use **File → Open**.

---

## What’s in the pipeline

```
Excel Ig diagram ──► Layer 1 identity ──┐
                                        ├──► CD8a-P01732.layers.json ──► layer1_viz.html
iCn3D H-bond table ──► Layer 4 density ─┘

Scatterplot *_common links ──► contact_maps.json ──► contact_maps.html
```

| Layer / product | Question it answers |
|-----------------|---------------------|
| **Layer 1 — static map** | Where is each residue on the IgV lace? (segment, IgStrand, PDB #) |
| **H-bonds** | How many H-bonds does each residue have? Who are the partners? |
| **Common contacts** | Which residue–residue interactions are shared across related Ig structures? |

---

## 1. Static representation (Layer 1)

**Input:** [`1CD8_ig_diagram.xlsx`](1CD8_ig_diagram.xlsx) — sheet `1. IgV`

| Sheet content | Becomes |
|---------------|---------|
| Cells like `C22` | amino acid + PDB position |
| Twin cells (+21 cols) | IgStrand codes (e.g. `2550`) |
| Row 6 headers | strand / segment labels (`A`, `B`, `A'`, …) |

**Build:**

```bash
.venv/bin/python build_layer1_cd8a.py
```

**Output fields** (per residue in `CD8a-P01732.layers.json` → `layers["1_identity"]`):

- `residue`, `pdb.pos`, `foldnum` (IgStrand), `segment`
- `spreadsheet_cell` / template row·col
- `seq_num`, `chain`

**View:** open [`layer1_viz.html`](layer1_viz.html) or [`layer1_viz_igstrand_colors.html`](layer1_viz_igstrand_colors.html)

---

## 2. Hydrogen bonds

**Input:** [`P01732-allinteraction.html`](P01732-allinteraction.html) — iCn3D “All interactions” H-bond table for P01732.

The same builder attaches H-bond density to Layer 4:

```bash
.venv/bin/python build_layer1_cd8a.py
# Layer 1 only:
.venv/bin/python build_layer1_cd8a.py --skip-hbonds
```

**Per residue** (`layers["4_interaction"]`):

- `hbond_density` / `hbond_count`
- backbone / sidechain / mixed counts
- partner list (atoms, distance, weight)

UniProt positions are filled using an auto-detected offset (`uniprot = pdb + 21` for 1CD8 / P01732).

**View** ([`layer1_viz.html`](layer1_viz.html)):

| Couleur mode | Meaning |
|--------------|---------|
| `none` | white cells |
| `segment` | IgStrand strand colors |
| `hbond_density` | pale → saturated teal (more H-bonds = darker) |

Toggle **H-bonds** to draw edges; the inspector lists partners, distances, and weights.

> **Note:** H-bond density is about hydrogen bonds on one protein. It is **not** the same as the common residue–residue contact maps below.

---

## 3. Common residue–residue interactions

Shared contacts across Ig-related structures (from scatterplot `*_common` blocks), shown as **contact-map grids** — not pairwise line drawings.

**Proteins:** P01732 · P10966 · P01730 · 1RHH · 1RHH2

**Input:** `P01732,P10966,P01730,1RHH,1RHH_scatterplot.json`  
(place next to the script if regenerating; uses `structure*_common` → `links`)

**Build:**

```bash
.venv/bin/python build_contact_maps.py
```

**Outputs:**

| File | Role |
|------|------|
| [`contact_maps.json`](contact_maps.json) | matrices + contact lists |
| [`contact_maps.html`](contact_maps.html) | interactive grids |

**How to read the HTML**

- **Per-protein map** — axes = PDB residue numbers; green cell = common contact in that structure  
- **Consensus (IgStrand)** — axes = shared fold seats (`B2550`, `E7548`, …); green = that fold-position contact is in the common set  
- Dark = no contact; green = contact (brighter ≈ stronger / more shared)

**What this solves:** residue numbers differ between proteins (`C43` vs `C23`). IgStrand puts them on one floor plan so the same structural contacts line up.

> These are **common interactions**, not MSA sequence conservation.

---

## Layered JSON (CD8α)

[`CD8a-P01732.layers.json`](CD8a-P01732.layers.json)

| Layer key | Contents |
|-----------|----------|
| `1_identity` | addressing from Excel |
| `2_intrinsic` | stub |
| `3_structural` | stub |
| `4_interaction` | H-bond density & partners |
| `5_comparative` | stub (MSA / mutations later) |

---

## MSA (related, separate track)

This repo also includes FASTA / MSA artifacts for multi-sequence alignment work (e.g. `P01732_P10966_P01730_1RHH_1RHH2_MSA.fasta`, `msa_script.py`).  

MSA answers “which letters align in sequence?” — complementary to contact maps (“which residues touch in the fold?”).

---

## Repo layout

```
1CD8_ig_diagram.xlsx          # Layer 1 geometry
P01732-allinteraction.html    # iCn3D H-bonds
build_layer1_cd8a.py          # build Layer 1 + H-bonds JSON
CD8a-P01732.layers.json
layer1_viz.html               # Layer 1 + H-bond viz
layer1_viz_igstrand_colors.html

build_contact_maps.py         # common contacts → maps
contact_maps.json
contact_maps.html
```

---

## Tips

- Prefer opening the committed `*.html` files directly; regenerate JSON when inputs change.
- IgStrand is the shared coordinate system for comparative Ig views.
- For collaboration: Layer 1 shows the map, H-bonds color chemistry on one chain, common contacts compare fold wiring across homologs.
