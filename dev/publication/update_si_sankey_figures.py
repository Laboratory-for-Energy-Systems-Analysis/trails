from __future__ import annotations

import os
import shutil
import time
import zipfile
from io import BytesIO
from pathlib import Path

from lxml import etree
from PIL import Image

REPO_ROOT = Path(__file__).resolve().parents[2]
DOCX_PATH = Path(
    "/Users/romain/Library/CloudStorage/OneDrive-PaulScherrerInstitut/trails/"
    "manuscript/Supplementary Information.docx"
)
OUTPUT_ROOT = REPO_ROOT / "dev/notebook_runs/temporal_lci_selected_indicator_sankeys"

REPLACEMENTS = {
    "word/media/image21.png": OUTPUT_ROOT
    / "polyol_42300/polyol_adaptive_cutoff_1e-4_sankey.png",
    "word/media/image22.png": OUTPUT_ROOT
    / "bev_42305/bev_adaptive_cutoff_1e-4_sankey.png",
    "word/media/image23.png": OUTPUT_ROOT
    / "daccs_42302/daccs_adaptive_cutoff_1e-4_sankey.png",
    "word/media/image24.png": OUTPUT_ROOT
    / "marine_42303/marine_adaptive_cutoff_1e-4_sankey.png",
}

NS = {
    "w": "http://schemas.openxmlformats.org/wordprocessingml/2006/main",
    "a": "http://schemas.openxmlformats.org/drawingml/2006/main",
    "r": "http://schemas.openxmlformats.org/officeDocument/2006/relationships",
    "wp": "http://schemas.openxmlformats.org/drawingml/2006/wordprocessingDrawing",
}


def _png_size(path: Path) -> tuple[int, int]:
    with Image.open(path) as image:
        return int(image.width), int(image.height)


def _read_rels(package: zipfile.ZipFile) -> dict[str, str]:
    rels = etree.fromstring(package.read("word/_rels/document.xml.rels"))
    return {rel.get("Id"): f"word/{rel.get('Target')}" for rel in rels}


def _update_document_xml(package: zipfile.ZipFile) -> tuple[bytes, dict[str, int]]:
    document = etree.fromstring(package.read("word/document.xml"))
    rel_targets = _read_rels(package)
    target_sizes = {
        target: _png_size(source) for target, source in REPLACEMENTS.items()
    }
    update_counts = {target: 0 for target in REPLACEMENTS}

    for paragraph in document.xpath(".//w:body/w:p", namespaces=NS):
        blips = paragraph.xpath(".//a:blip/@r:embed", namespaces=NS)
        targets = [rel_targets.get(blip) for blip in blips]
        matching_targets = [target for target in targets if target in REPLACEMENTS]
        if not matching_targets:
            continue

        target = matching_targets[0]
        image_width, image_height = target_sizes[target]
        extents = paragraph.xpath(".//wp:extent", namespaces=NS)
        if not extents:
            raise RuntimeError(f"Missing wp:extent for {target}")
        width_emu = int(extents[0].get("cx"))
        height_emu = round(width_emu * image_height / image_width)

        for extent in extents:
            extent.set("cx", str(width_emu))
            extent.set("cy", str(height_emu))
            update_counts[target] += 1
        for extent in paragraph.xpath(".//a:ext", namespaces=NS):
            extent.set("cx", str(width_emu))
            extent.set("cy", str(height_emu))

    missing = [target for target, count in update_counts.items() if count == 0]
    if missing:
        raise RuntimeError("Did not find document drawings for: " + ", ".join(missing))

    return (
        etree.tostring(
            document,
            xml_declaration=True,
            encoding="UTF-8",
            standalone="yes",
        ),
        update_counts,
    )


def _replace_docx() -> Path:
    missing = [
        path for path in [DOCX_PATH, *REPLACEMENTS.values()] if not path.exists()
    ]
    if missing:
        raise FileNotFoundError("Missing files:\n" + "\n".join(map(str, missing)))

    timestamp = time.strftime("%Y%m%d-%H%M")
    backup = DOCX_PATH.with_name(f"{DOCX_PATH.name}.bak-codex-{timestamp}")
    shutil.copy2(DOCX_PATH, backup)

    tmp_path = DOCX_PATH.with_name(f"{DOCX_PATH.name}.tmp-codex")
    try:
        with zipfile.ZipFile(DOCX_PATH, "r") as source:
            document_xml, update_counts = _update_document_xml(source)
            with zipfile.ZipFile(tmp_path, "w") as target:
                for info in source.infolist():
                    data = source.read(info.filename)
                    if info.filename == "word/document.xml":
                        data = document_xml
                    elif info.filename in REPLACEMENTS:
                        data = REPLACEMENTS[info.filename].read_bytes()
                    target.writestr(info, data)

        os.replace(tmp_path, DOCX_PATH)
    finally:
        if tmp_path.exists():
            tmp_path.unlink()

    with zipfile.ZipFile(DOCX_PATH, "r") as package:
        bad_file = package.testzip()
        if bad_file is not None:
            raise RuntimeError(f"Updated DOCX zip failed at {bad_file}")
        for target, source in REPLACEMENTS.items():
            if package.read(target) != source.read_bytes():
                raise RuntimeError(f"Replacement mismatch for {target}")
        document = etree.fromstring(package.read("word/document.xml"))
        for target in REPLACEMENTS:
            assert update_counts[target] > 0

    return backup


if __name__ == "__main__":
    backup_path = _replace_docx()
    print(f"Updated {DOCX_PATH}")
    print(f"Backup: {backup_path}")
    for target, source in REPLACEMENTS.items():
        print(f"{target} <- {source}")
