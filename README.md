# IgFolds2D Residue Interactions

Branch: `Hbonds`

Layered feature vectors for an IgV domain (CD8α / PDB **1CD8**).

## Inputs

| File | Role |
|------|------|
| `1CD8_ig_diagram.xlsx` | Layer 1 identity (sheet `1. IgV`) |
| `P01732-allinteraction.html` | Layer 4 H-bond density (iCn3D table) |

## Build data

```bash
.venv/bin/python build_layer1_cd8a.py
```

**Output:** `CD8a-P01732.layers.json`

## View

Open this file in a browser (double-click or File → Open):

```
layer1_viz.html
```

It is self-contained (data embedded). No extra build step.

**Couleur modes** (like GPCRfold2D):
- `none` — white cells
- `segment` — strand colors
- `hbond_density` — pale→saturated teal ramp (same as `helix4_density`)

Toggle **H-bonds** for edges; inspector lists partners, distances, weights.
