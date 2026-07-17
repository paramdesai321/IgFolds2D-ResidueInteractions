# Feature vectors (layered) — Excel-only Layer 1

## Input (required)

Only this file:

```
1CD8_ig_diagram.xlsx
```

Sheet `1. IgV` must have:

- left map: residues like `C22`
- right map (same places, +21 columns): IgStrand numbers like `2550`
- row 6 strand headers (`A`, `B`, `A'`, …)

## Build

```bash
.venv/bin/python feature_vectors/build_layer1_cd8a.py ~/Downloads/1CD8_ig_diagram.xlsx
```

Or with defaults:

```bash
.venv/bin/python feature_vectors/build_layer1_cd8a.py
```

## Output

`feature_vectors/CD8a-P01732.layers.json`

## What Layer 1 contains

| Field | Source in Excel |
|-------|-----------------|
| residue | left cell letter |
| pdb.pos | left cell number |
| foldnum (IgStrand) | right twin cell |
| spreadsheet_cell | Excel address |
| template_position | same as spreadsheet row/col |
| segment | row-6 header (else inferred from IgStrand) |
| domain | sheet title `IgV` |
| seq_num | order by IgStrand |
| chain | default `A` (1CD8) |

`uniprot.pos` is `null` for now — that index is not written on the sheet.

Layers 2–5 are empty stubs for later (H-bonds, MSA, etc.).
