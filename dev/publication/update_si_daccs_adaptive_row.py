from __future__ import annotations

from pathlib import Path

from docx import Document


DOCX = Path(
    "/Users/romain/Library/CloudStorage/OneDrive-PaulScherrerInstitut/"
    "trails/manuscript/Supplementary Information.docx"
)

ROW_VALUES = [
    "Adaptive routing (0.01% cutoff)",
    "49,905",
    "118,050",
    "92.70",
    "-5.3",
    "-5.4",
    "89.6",
    "72.8",
]


def _cell_text(cell) -> str:
    return " ".join(cell.text.split())


def main() -> None:
    doc = Document(str(DOCX))
    target_table = None
    for table in doc.tables:
        if not table.rows:
            continue
        header = [_cell_text(cell) for cell in table.rows[0].cells]
        if header[:8] == [
            "Depth",
            "Routed nodes",
            "Routed edges",
            "PM score",
            "Deviation from static",
            "Deviation (%)",
            "Routing time (s)",
            "LCA time (s)",
        ]:
            target_table = table
            break

    if target_table is None:
        raise RuntimeError("Could not find the DACCS PM depth-sweep table.")

    target_row = None
    for row in target_table.rows:
        if _cell_text(row.cells[0]).startswith("Adaptive routing"):
            target_row = row
            break

    if target_row is None:
        target_row = target_table.add_row()

    if len(target_row.cells) < len(ROW_VALUES):
        raise RuntimeError("Target table has fewer cells than expected.")

    for cell, value in zip(target_row.cells, ROW_VALUES, strict=True):
        cell.text = value

    doc.save(str(DOCX))


if __name__ == "__main__":
    main()
