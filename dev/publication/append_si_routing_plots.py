from __future__ import annotations

import csv
import shutil
from pathlib import Path

from docx import Document
from docx.enum.section import WD_ORIENT, WD_SECTION
from docx.enum.table import WD_CELL_VERTICAL_ALIGNMENT
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.shared import Emu, Inches, Pt


REPO_ROOT = Path(__file__).resolve().parents[2]
RUN_DIR = (
    REPO_ROOT
    / "dev"
    / "publication"
    / "notebook_runs"
    / "temporal_lci_routing_comparison_clean"
)
SCORES_CSV = RUN_DIR / "routing_comparison_scores.csv"

SOURCE_DOCX = Path(
    "/Users/romain/Library/CloudStorage/OneDrive-PaulScherrerInstitut/"
    "trails/manuscript/Supplementary Information.docx"
)
OUTPUT_DOCX = Path("/tmp/Supplementary Information.with_routing_plots.docx")

CASE_ORDER = ["bev", "polyol", "marine", "daccs"]
CASE_LABELS = {
    "bev": "Battery electric vehicle",
    "polyol": "CCU polyol",
    "marine": "Marine fuel switch",
    "daccs": "Direct air capture and carbon storage",
}


def slug(value: str, max_length: int = 120) -> str:
    text = "".join(char.lower() if char.isalnum() else "_" for char in str(value))
    text = "_".join(part for part in text.split("_") if part)
    return (text or "item")[:max_length]


def method_label(method: str) -> str:
    prefix = "EF v3.1 - "
    return method.removeprefix(prefix) if method.startswith(prefix) else method


def methods_by_case() -> dict[str, list[str]]:
    out: dict[str, list[str]] = {}
    with SCORES_CSV.open(newline="") as handle:
        for row in csv.DictReader(handle):
            case = row["case"]
            out.setdefault(case, []).append(row["method"])
    return out


def image_paths(case: str, method: str) -> tuple[Path, Path]:
    stem = f"{slug(method, 100)}_routing_comparison"
    case_dir = RUN_DIR / case
    unstacked = case_dir / f"{stem}_unstacked.png"
    stacked = case_dir / f"{stem}_stacked.png"
    missing = [path for path in [unstacked, stacked] if not path.exists()]
    if missing:
        raise FileNotFoundError("Missing routing plot(s):\n" + "\n".join(map(str, missing)))
    return unstacked, stacked


def clear_cell(cell) -> None:
    for paragraph in cell.paragraphs:
        paragraph.clear()


def set_cell_text(cell, text: str, *, bold: bool = False, size: int = 8) -> None:
    clear_cell(cell)
    paragraph = cell.paragraphs[0]
    paragraph.alignment = WD_ALIGN_PARAGRAPH.LEFT
    run = paragraph.add_run(text)
    run.bold = bold
    run.font.size = Pt(size)


def add_picture(cell, path: Path, width: Emu) -> None:
    clear_cell(cell)
    paragraph = cell.paragraphs[0]
    paragraph.alignment = WD_ALIGN_PARAGRAPH.CENTER
    paragraph.paragraph_format.space_after = Pt(0)
    paragraph.add_run().add_picture(str(path), width=width)
    cell.vertical_alignment = WD_CELL_VERTICAL_ALIGNMENT.CENTER


def set_table_column_widths(table, width: Emu) -> None:
    for row in table.rows:
        for cell in row.cells:
            cell.width = width


def append_case_table(
    document: Document,
    *,
    case: str,
    methods: list[str],
    picture_width: Emu,
    column_width: Emu,
) -> None:
    document.add_heading(CASE_LABELS.get(case, case), level=2)

    table = document.add_table(rows=1, cols=2)
    table.style = "Table Grid"
    table.autofit = False
    set_table_column_widths(table, column_width)
    set_cell_text(table.rows[0].cells[0], "Unstacked annual impacts", bold=True, size=8)
    set_cell_text(table.rows[0].cells[1], "Stacked annual impacts", bold=True, size=8)

    for index, method in enumerate(methods, start=1):
        label_row = table.add_row()
        label_cell = label_row.cells[0].merge(label_row.cells[1])
        set_cell_text(
            label_cell,
            f"Figure S-routing-{case}-{index:02d}. {method_label(method)}.",
            bold=True,
            size=8,
        )

        image_row = table.add_row()
        unstacked, stacked = image_paths(case, method)
        add_picture(image_row.cells[0], unstacked, picture_width)
        add_picture(image_row.cells[1], stacked, picture_width)


def append_routing_plot_appendix(document: Document) -> None:
    section = document.add_section(WD_SECTION.NEW_PAGE)
    section.orientation = WD_ORIENT.LANDSCAPE
    section.page_width, section.page_height = section.page_height, section.page_width
    section.left_margin = Inches(0.45)
    section.right_margin = Inches(0.45)
    section.top_margin = Inches(0.45)
    section.bottom_margin = Inches(0.45)

    usable_width = (
        section.page_width - section.left_margin - section.right_margin
    )
    column_width = Emu(int(usable_width / 2))
    picture_width = Emu(int(usable_width / 2 * 0.94))

    document.add_heading("Routing-comparison plots for all case-study indicators", level=1)
    paragraph = document.add_paragraph()
    paragraph.add_run(
        "This appendix reports the routing-comparison plots generated for all "
        "European Footprint v3.1 indicators in each case study. Each row pairs "
        "the unstacked annual-impact view with the corresponding stacked "
        "annual-impact view. The cumulative foreground-only, adaptive-routing, "
        "and static reference scores are shown as in the main case-study figures."
    )

    methods = methods_by_case()
    for position, case in enumerate(CASE_ORDER):
        if case not in methods:
            raise ValueError(f"Missing case in {SCORES_CSV}: {case}")
        if position:
            document.add_page_break()
        append_case_table(
            document,
            case=case,
            methods=methods[case],
            picture_width=picture_width,
            column_width=column_width,
        )


def main() -> None:
    if not SOURCE_DOCX.exists():
        raise FileNotFoundError(SOURCE_DOCX)
    if not SCORES_CSV.exists():
        raise FileNotFoundError(SCORES_CSV)

    working_copy = OUTPUT_DOCX.with_suffix(".work.docx")
    shutil.copy2(SOURCE_DOCX, working_copy)
    document = Document(str(working_copy))
    append_routing_plot_appendix(document)
    document.save(OUTPUT_DOCX)
    print(OUTPUT_DOCX)


if __name__ == "__main__":
    main()
