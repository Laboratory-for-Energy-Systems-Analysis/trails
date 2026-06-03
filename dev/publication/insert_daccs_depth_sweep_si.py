from __future__ import annotations

import csv
import shutil
import tempfile
from datetime import datetime
from io import BytesIO
from pathlib import Path
from zipfile import ZIP_DEFLATED, ZipFile
import xml.etree.ElementTree as ET

REPO_ROOT = Path(__file__).resolve().parents[2]
DOCX_PATH = Path(
    "/Users/romain/Library/CloudStorage/OneDrive-PaulScherrerInstitut/"
    "trails/manuscript/Supplementary Information.docx"
)
CSV_PATH = (
    REPO_ROOT
    / "dev"
    / "notebook_runs"
    / "daccs_pm_depth_sweep"
    / "daccs_pm_depth_sweep.csv"
)

W_NS = "http://schemas.openxmlformats.org/wordprocessingml/2006/main"
XML_NS = "http://www.w3.org/XML/1998/namespace"


def qn(name: str) -> str:
    prefix, local = name.split(":")
    if prefix == "w":
        return f"{{{W_NS}}}{local}"
    if prefix == "xml":
        return f"{{{XML_NS}}}{local}"
    raise ValueError(name)


def collect_namespaces(xml_bytes: bytes) -> None:
    seen: set[tuple[str, str]] = set()
    for _event, namespace in ET.iterparse(BytesIO(xml_bytes), events=("start-ns",)):
        if namespace in seen:
            continue
        seen.add(namespace)
        prefix, uri = namespace
        if prefix == "xml":
            continue
        ET.register_namespace(prefix, uri)


def text_of(element: ET.Element) -> str:
    return "".join(t.text or "" for t in element.findall(f".//{qn('w:t')}")).strip()


def p_style(element: ET.Element) -> str:
    pstyle = element.find(f"./{qn('w:pPr')}/{qn('w:pStyle')}")
    return "" if pstyle is None else pstyle.get(qn("w:val"), "")


def paragraph(text: str, style: str | None = "para") -> ET.Element:
    p = ET.Element(qn("w:p"))
    if style:
        ppr = ET.SubElement(p, qn("w:pPr"))
        ET.SubElement(ppr, qn("w:pStyle"), {qn("w:val"): style})
    r = ET.SubElement(p, qn("w:r"))
    t = ET.SubElement(r, qn("w:t"))
    if text.startswith(" ") or text.endswith(" "):
        t.set(qn("xml:space"), "preserve")
    t.text = text
    return p


def cell(text: str, *, bold: bool = False) -> ET.Element:
    tc = ET.Element(qn("w:tc"))
    tcpr = ET.SubElement(tc, qn("w:tcPr"))
    ET.SubElement(tcpr, qn("w:tcW"), {qn("w:w"): "0", qn("w:type"): "auto"})
    p = ET.SubElement(tc, qn("w:p"))
    r = ET.SubElement(p, qn("w:r"))
    if bold:
        rpr = ET.SubElement(r, qn("w:rPr"))
        ET.SubElement(rpr, qn("w:b"))
    t = ET.SubElement(r, qn("w:t"))
    if text.startswith(" ") or text.endswith(" "):
        t.set(qn("xml:space"), "preserve")
    t.text = text
    return tc


def table(rows: list[list[str]]) -> ET.Element:
    tbl = ET.Element(qn("w:tbl"))
    tblpr = ET.SubElement(tbl, qn("w:tblPr"))
    ET.SubElement(tblpr, qn("w:tblStyle"), {qn("w:val"): "TableGrid"})
    ET.SubElement(tblpr, qn("w:tblW"), {qn("w:w"): "0", qn("w:type"): "auto"})
    tblgrid = ET.SubElement(tbl, qn("w:tblGrid"))
    for _ in rows[0]:
        ET.SubElement(tblgrid, qn("w:gridCol"), {qn("w:w"): "1200"})
    for row_index, values in enumerate(rows):
        tr = ET.SubElement(tbl, qn("w:tr"))
        for value in values:
            tr.append(cell(value, bold=row_index == 0))
    return tbl


def as_float(row: dict[str, str], key: str) -> float:
    value = str(row.get(key, "")).strip()
    return float(value) if value else float("nan")


def format_seconds(value: float) -> str:
    if value != value:
        return "-"
    return f"{value:.1f}"


def format_int(value: float) -> str:
    if value != value:
        return "-"
    return f"{int(value):,}"


def load_table_rows() -> tuple[list[list[str]], dict[str, float]]:
    with CSV_PATH.open("r", newline="", encoding="utf-8") as handle:
        data = list(csv.DictReader(handle))

    static = next(row for row in data if row["mode"] == "static")
    temporal = [row for row in data if row["mode"] == "temporal"]
    temporal.sort(key=lambda row: int(float(row["depth"])))

    static_score = as_float(static, "static_score")
    rows = [
        [
            "Depth",
            "Routed nodes",
            "Routed edges",
            "PM score",
            "Deviation from static",
            "Deviation (%)",
            "Routing time (s)",
            "LCA time (s)",
        ],
        [
            "Static",
            "-",
            "-",
            f"{static_score:.3f}",
            "0.000",
            "0.00",
            "-",
            format_seconds(as_float(static, "static_lca_seconds")),
        ],
    ]

    for row in temporal:
        score = as_float(row, "score")
        deviation = as_float(row, "score_deviation_from_static")
        relative = 100.0 * as_float(row, "relative_deviation_from_static")
        rows.append(
            [
                str(int(float(row["depth"]))),
                format_int(as_float(row, "graph_nodes")),
                format_int(as_float(row, "graph_edges")),
                f"{score:.3f}",
                f"{deviation:.3f}",
                f"{relative:.2f}",
                format_seconds(as_float(row, "routing_seconds")),
                format_seconds(as_float(row, "temporal_lca_seconds")),
            ]
        )

    depth7 = next(row for row in temporal if int(float(row["depth"])) == 7)
    return rows, {
        "static_score": static_score,
        "depth1_score": as_float(temporal[0], "score"),
        "depth1_relative": 100.0
        * as_float(temporal[0], "relative_deviation_from_static"),
        "depth6_score": as_float(
            next(r for r in temporal if r["depth"] == "6"), "score"
        ),
        "depth6_relative": 100.0
        * as_float(
            next(r for r in temporal if r["depth"] == "6"),
            "relative_deviation_from_static",
        ),
        "depth7_score": as_float(depth7, "score"),
        "depth7_relative": 100.0 * as_float(depth7, "relative_deviation_from_static"),
        "depth7_nodes": as_float(depth7, "graph_nodes"),
        "depth7_edges": as_float(depth7, "graph_edges"),
        "depth7_routing": as_float(depth7, "routing_seconds"),
        "depth7_lca": as_float(depth7, "temporal_lca_seconds"),
    }


def section_elements() -> list[ET.Element]:
    table_rows, stats = load_table_rows()
    minutes = stats["depth7_routing"] / 60.0

    return [
        paragraph(
            "Depth sensitivity and computational cost for a direct-air-capture case",
            "Heading2",
        ),
        paragraph(
            "We additionally tested the effect of temporalisation depth on both "
            "runtime and score for the direct-air-capture case study. The functional "
            "unit was 20 billion kg of captured CO2 in 2025 for the activity "
            '"carbon dioxide, captured, with a solvent-based direct air capture '
            'system, 1MtCO2". The calculation used the SSP2-PkBudg1000 Trails data '
            "package, the same foreground inventories as the depth-sweep notebook, "
            "a routing cutoff of 1e-3, and the EF v3.1 particulate matter formation "
            "indicator. The static matrix LCA score was "
            f"{stats['static_score']:.3f} disease incidence.",
        ),
        paragraph(
            "The experiment was designed as a computational sensitivity test rather "
            "than as a recommended production setting. Depth 1 routes only the first "
            "temporalised layer below the functional unit, whereas larger depths "
            "allow temporal distributions to be followed further into the "
            "background supply chain before remaining frontier demands are solved "
            "with year-specific matrices.",
        ),
        table(table_rows),
        paragraph(
            "The routed graph grew strongly with temporalisation depth. The number "
            "of routed nodes increased from 197 at depth 1 to "
            f"{int(stats['depth7_nodes']):,} at depth 7, while the number of edges "
            f"increased from 196 to {int(stats['depth7_edges']):,}. Routing time "
            f"therefore dominated the high-depth calculations: at depth 7, routing "
            f"took {stats['depth7_routing']:.1f} s ({minutes:.1f} min), whereas "
            f"the subsequent LCA solve and scoring step took "
            f"{stats['depth7_lca']:.1f} s. This behaviour is consistent with an "
            "exponential worst-case expansion of the process-year routing graph, "
            "although the realised growth is moderated by repeated-node aggregation "
            "and the routing cutoff.",
        ),
        paragraph(
            "The score deviation from the static result decreased substantially "
            "between shallow and intermediate depths, from "
            f"{stats['depth1_relative']:.2f}% at depth 1 to "
            f"{stats['depth6_relative']:.2f}% at depth 6. The depth-7 score was "
            f"{stats['depth7_score']:.3f} disease incidence, corresponding to a "
            f"{stats['depth7_relative']:.2f}% deviation from the static score. The "
            "small non-monotonic change between depths 6 and 7 shows that increasing "
            "depth is not a simple numerical convergence parameter: additional "
            "routed paths can shift contributions across years and can therefore "
            "change both the timing and the target-year background coefficients of "
            "the assessed inventory.",
        ),
    ]


def remove_existing_section(body: ET.Element) -> None:
    children = list(body)
    start = None
    for index, element in enumerate(children):
        if element.tag == qn("w:p") and text_of(element) == (
            "Depth sensitivity and computational cost for a direct-air-capture case"
        ):
            start = index
            break
    if start is None:
        return

    end = len(children)
    for index in range(start + 1, len(children)):
        element = children[index]
        if element.tag == qn("w:p") and p_style(element) in {"Heading1", "Heading2"}:
            end = index
            break
    for element in children[start:end]:
        body.remove(element)


def insertion_index(body: ET.Element) -> int:
    for index, element in enumerate(list(body)):
        if element.tag == qn("w:p") and text_of(element) == (
            "Coupling Trails inventories with FaIR"
        ):
            return index
    raise ValueError("Could not locate the FaIR coupling heading.")


def update_docx() -> Path:
    if not DOCX_PATH.exists():
        raise FileNotFoundError(DOCX_PATH)
    if not CSV_PATH.exists():
        raise FileNotFoundError(CSV_PATH)

    timestamp = datetime.now().strftime("%Y%m%d-%H%M%S")
    backup = DOCX_PATH.with_name(
        f"{DOCX_PATH.stem}.backup-before-daccs-depth-sweep-{timestamp}.docx"
    )
    shutil.copy2(DOCX_PATH, backup)

    with ZipFile(DOCX_PATH, "r") as zin:
        document_xml = zin.read("word/document.xml")
        collect_namespaces(document_xml)
        root = ET.fromstring(document_xml)
        body = root.find(qn("w:body"))
        if body is None:
            raise RuntimeError("Document body not found.")

        remove_existing_section(body)
        idx = insertion_index(body)
        for offset, element in enumerate(section_elements()):
            body.insert(idx + offset, element)

        updated_xml = ET.tostring(root, encoding="utf-8", xml_declaration=True)

        with tempfile.NamedTemporaryFile(
            "wb", suffix=".docx", dir=DOCX_PATH.parent, delete=False
        ) as handle:
            temp_path = Path(handle.name)

        with ZipFile(temp_path, "w", ZIP_DEFLATED) as zout:
            for item in zin.infolist():
                payload = (
                    updated_xml
                    if item.filename == "word/document.xml"
                    else zin.read(item.filename)
                )
                zout.writestr(item, payload)

    shutil.move(str(temp_path), DOCX_PATH)
    with ZipFile(DOCX_PATH, "r") as check:
        bad = check.testzip()
    if bad is not None:
        raise RuntimeError(f"Updated DOCX failed ZIP validation at {bad}.")
    return backup


if __name__ == "__main__":
    backup_path = update_docx()
    print(f"Updated: {DOCX_PATH}")
    print(f"Backup: {backup_path}")
