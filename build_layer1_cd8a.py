"""
Build Layer-1 feature vectors from the Excel Ig diagram ONLY.

Input:
  1CD8_ig_diagram.xlsx  (sheet "1. IgV")

From the sheet we read:
  - left cells  : residue letter + PDB-style number (e.g. C22)
  - right cells : IgStrand code (same geometry, col + 21)
  - row 6       : strand / segment headers (A, B, A', …)
  - title       : IgV domain

No assignment.json / template.json required.
"""

from __future__ import annotations

import json
import re
from pathlib import Path

from openpyxl import load_workbook
from openpyxl.utils import get_column_letter

ROOT = Path(__file__).resolve().parents[1]
EXCEL_DEFAULT = Path("/Users/nicolemathias/Downloads/1CD8_ig_diagram.xlsx")
OUT = ROOT / "feature_vectors" / "CD8a-P01732.layers.json"

# Sensible defaults for THIS spreadsheet (not read as external files)
CHAIN = "A"
PDB_ID = "1CD8"
UNIPROT_ID = "P01732"  # known for 1CD8; positions left null unless labeled in sheet
DOMAIN = "IgV"

RES_PAT = re.compile(r"^([A-Z])(\d+)$")

STRAND_FROM_IG = {
    1: "A",
    2: "B",
    3: "C",
    4: "C'",
    5: "C''",
    6: "D",
    7: "E",
    8: "F",
    9: "G",
}


def empty_layers() -> dict:
    return {
        "2_intrinsic": {},
        "3_structural": {},
        "4_interaction": {},
        "5_comparative": {},
    }


def segment_from_igstrand(ig: int) -> str:
    thousands = ig // 1000
    hundreds = (ig // 100) % 10
    if thousands == 1 and hundreds == 8:
        return "A'"
    return STRAND_FROM_IG.get(thousands, f"strand_{thousands}")


def is_anchor(ig: int) -> bool:
    return ig % 100 == 50


def parse_excel(path: Path) -> list[dict]:
    wb = load_workbook(path, data_only=True)
    ws = wb["1. IgV"]

    # Strand headers on the LEFT map (row 6)
    strand_by_col: dict[int, str] = {}
    for c in range(4, 17):
        v = ws.cell(6, c).value
        if v is not None and str(v).strip():
            strand_by_col[c] = str(v).strip()

    domain = DOMAIN
    title = ws.cell(1, 9).value  # "IgV" on the sheet
    if title and str(title).strip():
        domain = str(title).strip()

    records: list[dict] = []
    for row in ws.iter_rows(min_row=1, max_row=40, max_col=20):
        for cell in row:
            if cell.value is None:
                continue
            m = RES_PAT.match(str(cell.value).strip())
            if not m:
                continue

            aa = m.group(1)
            pdb_pos = int(m.group(2))
            ig_raw = ws.cell(cell.row, cell.column + 21).value
            if not isinstance(ig_raw, int):
                # skip residue cells without a twin IgStrand number
                continue
            ig = int(ig_raw)

            segment = strand_by_col.get(cell.column) or segment_from_igstrand(ig)
            address = f"{get_column_letter(cell.column)}{cell.row}"

            records.append(
                {
                    "aa": aa,
                    "pdb_pos": pdb_pos,
                    "igstrand": ig,
                    "segment": segment,
                    "domain": domain,
                    "spreadsheet_cell": address,
                    "spreadsheet_row": cell.row,
                    "spreadsheet_col": cell.column,
                }
            )

    # Chain order ≈ IgStrand numerical order (canonical for this numbering)
    records.sort(key=lambda r: r["igstrand"])
    return records


def build(excel_path: Path) -> dict:
    records = parse_excel(excel_path)

    residues = []
    for seq_num, rec in enumerate(records, start=1):
        ig = rec["igstrand"]
        layer1 = {
            "seq_num": seq_num,
            "residue": rec["aa"],
            "chain": CHAIN,
            "domain": rec["domain"],
            "segment": rec["segment"],
            # Template position = spreadsheet geometry (Excel is the template)
            "template_position": {
                "row": rec["spreadsheet_row"],
                "col": rec["spreadsheet_col"],
            },
            "spreadsheet_cell": rec["spreadsheet_cell"],
            "spreadsheet_position": {
                "row": rec["spreadsheet_row"],
                "col": rec["spreadsheet_col"],
                "address": rec["spreadsheet_cell"],
            },
            "foldnum": {
                "scheme": "IgStrand",
                "code": ig,
                "anchor": is_anchor(ig),
            },
            # UniProt absolute index is NOT on the sheet → leave null for now
            "uniprot": {"id": UNIPROT_ID, "pos": None},
            "pdb": {"id": PDB_ID, "chain": CHAIN, "pos": rec["pdb_pos"]},
        }

        residues.append(
            {
                "id": f"{PDB_ID}:{CHAIN}:{rec['pdb_pos']}",
                "layers": {
                    "1_identity": layer1,
                    **empty_layers(),
                },
            }
        )

    return {
        "schema_version": "0.1",
        "fold": "IgV",
        "protein": {
            "id": UNIPROT_ID,
            "pdb": PDB_ID,
            "name": "CD8A",
            "domain": DOMAIN,
        },
        "template": {
            "id": "excel:1CD8_ig_diagram",
            "numbering_scheme": "IgStrand",
            "source": str(excel_path),
        },
        "sources": {
            "excel": str(excel_path),
            "pdb": {"id": PDB_ID, "chain": CHAIN},
            "uniprot": {"id": UNIPROT_ID, "note": "id known; residue positions not on sheet"},
        },
        "layer_defs": {
            "1_identity": {
                "description": "Addressing & identity from Excel only",
                "fields": [
                    "seq_num",
                    "residue",
                    "chain",
                    "domain",
                    "segment",
                    "template_position",
                    "spreadsheet_cell",
                    "foldnum",
                    "uniprot",
                    "pdb",
                ],
            },
            "2_intrinsic": {"description": "From amino-acid type alone"},
            "3_structural": {"description": "From 3D model"},
            "4_interaction": {"description": "H-bonds, contacts, weights"},
            "5_comparative": {"description": "MSA / mutations"},
        },
        "counts": {
            "residues": len(residues),
            "anchors": sum(
                1 for r in residues if r["layers"]["1_identity"]["foldnum"]["anchor"]
            ),
            "segments": sorted(
                {
                    r["layers"]["1_identity"]["segment"]
                    for r in residues
                    if r["layers"]["1_identity"]["segment"]
                }
            ),
        },
        "residues": residues,
    }


def main() -> None:
    import argparse

    p = argparse.ArgumentParser(description="Build Layer-1 JSON from Ig Excel diagram only")
    p.add_argument(
        "excel",
        nargs="?",
        type=Path,
        default=EXCEL_DEFAULT,
        help=f"Path to 1CD8_ig_diagram.xlsx (default: {EXCEL_DEFAULT})",
    )
    p.add_argument(
        "-o",
        "--output",
        type=Path,
        default=OUT,
        help=f"Output JSON path (default: {OUT})",
    )
    args = p.parse_args()

    doc = build(args.excel.expanduser().resolve())
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(doc, indent=2))
    print(f"wrote {args.output}")
    print(
        f"residues={doc['counts']['residues']} "
        f"anchors={doc['counts']['anchors']} "
        f"segments={doc['counts']['segments']}"
    )


if __name__ == "__main__":
    main()
