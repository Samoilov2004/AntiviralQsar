from __future__ import annotations

import csv
from dataclasses import dataclass
from io import BytesIO, StringIO
from pathlib import PurePosixPath
import re
import zipfile

from .sources import clean_text, xml_parse


@dataclass
class Chunk:
    locator: str
    text: str
    context: str = ""


def table_grid(table):
    rows, occupied = [], {}
    for ri, tr in enumerate(table.xpath(".//tr")):
        cells, col = [], 0
        for cell in tr.xpath("./th|./td"):
            while (ri, col) in occupied:
                cells.append(occupied[(ri, col)])
                col += 1
            content = clean_text(cell)
            try:
                colspan, rowspan = int(cell.get("colspan", 1)), int(cell.get("rowspan", 1))
            except ValueError:
                colspan = rowspan = 1
            if not 1 <= colspan <= 100 or not 1 <= rowspan <= 1000:
                raise ValueError("unsupported_table_spans")
            for dr in range(rowspan):
                for dc in range(colspan):
                    occupied[(ri + dr, col + dc)] = content
            cells.extend([content] * colspan)
            col += colspan
        while (ri, col) in occupied:
            cells.append(occupied[(ri, col)])
            col += 1
        rows.append(cells)
    return rows


def row_chunks(rows, locator, caption="", context="", limit=12000, header_rows=1):
    header = "\n".join(" | ".join(map(str, r)) for r in rows[:header_rows])
    prefix = f"{caption}\nColumn headers (context only):\n{header}\n"
    current, first = [], 0
    for i, row in enumerate(rows):
        text = f"ROW {i+1}: " + " | ".join("" if x is None else str(x) for x in row)
        if current and (len(current) >= 8 or len(prefix) + sum(map(len, current)) + len(text) > limit):
            yield Chunk(f"{locator}:rows-{first+1}-{i}", prefix + "\n".join(current), context)
            current, first = [], i
        current.append(text)
    if current:
        yield Chunk(f"{locator}:rows-{first+1}-{len(rows)}", prefix + "\n".join(current), context)


def text_chunks(text, locator, context="", limit=12000):
    paragraphs = re.split(r"\n\s*\n", text)
    current, size, part = [], 0, 1
    for paragraph in paragraphs:
        # Long unbroken paragraphs are explicitly segmented, without silent truncation.
        pieces = [paragraph[i:i+limit] for i in range(0, len(paragraph), limit)] or [""]
        for piece in pieces:
            if current and size + len(piece) > limit:
                yield Chunk(f"{locator}:part-{part}", "\n\n".join(current), context)
                current, size, part = [], 0, part + 1
            current.append(piece)
            size += len(piece)
    if current:
        yield Chunk(f"{locator}:part-{part}", "\n\n".join(current), context)


def jats_chunks(data):
    root = xml_parse(data)
    title = " ".join(root.xpath(".//article-title//text()"))
    paragraphs = root.xpath(".//abstract//p | .//body//p[not(ancestor::table-wrap) and not(ancestor::fig)]")
    context_parts = [clean_text(p) for p in paragraphs if re.search(
        r"\b(cells?|cellular|cytotoxicity|antiviral assay|infection|infected)\b", clean_text(p), re.I)]
    context = (title + "\n" + "\n".join(context_parts))[:6000]
    chunks = []
    for i, table in enumerate(root.xpath(".//table-wrap")):
        label = table.get("id") or f"table-{i+1}"
        caption = " ".join(clean_text(x) for x in table.xpath("./label|./caption|./table-wrap-foot"))
        rows = table_grid(table)
        if rows:
            header_rows = 0
            for tr in table.xpath(".//tr"):
                if tr.xpath("./th"):
                    header_rows += 1
                else:
                    break
            chunks.extend(row_chunks(rows, label, caption, context, header_rows=header_rows or 1))
    text = "\n\n".join(clean_text(p) for p in paragraphs)
    chunks.extend(text_chunks(text, "body", context))
    attachments = []
    for element in root.xpath(".//supplementary-material//*[@*[local-name()='href']] | .//supplementary-material[@*[local-name()='href']]"):
        for key, href in element.attrib.items():
            if key.endswith("href"):
                attachments.append(href)
    return chunks, list(dict.fromkeys(attachments))


def attachment_chunks(name, data, max_bytes=50_000_000):
    if len(data) > max_bytes:
        return [], [f"attachment_too_large:{name}"]
    suffix = PurePosixPath(name.lower()).suffix
    chunks, warnings = [], []
    if suffix == ".pdf":
        import pymupdf
        with pymupdf.open(stream=data, filetype="pdf") as pdf:
            for i, page in enumerate(pdf):
                text = page.get_text(sort=True)
                if len(text.strip()) < 40:
                    warnings.append(f"image_only_page:{name}:page-{i+1}")
                    continue
                chunks.extend(text_chunks(text, f"{name}:page-{i+1}"))
                try:
                    for j, table in enumerate(page.find_tables().tables):
                        chunks.extend(row_chunks(table.extract(), f"{name}:page-{i+1}:table-{j+1}"))
                except (ValueError, RuntimeError):
                    warnings.append(f"pdf_table_parse_failed:{name}:page-{i+1}")
    elif suffix in (".csv", ".tsv", ".txt"):
        text = data.decode("utf-8-sig", errors="replace")
        if suffix == ".txt":
            chunks.extend(text_chunks(text, name))
        else:
            rows = list(csv.reader(StringIO(text), delimiter="\t" if suffix == ".tsv" else ","))
            chunks.extend(row_chunks(rows, name))
    elif suffix == ".xlsx":
        from openpyxl import load_workbook
        book = load_workbook(BytesIO(data), read_only=True, data_only=True)
        try:
            for sheet in book:
                chunks.extend(row_chunks(list(sheet.iter_rows(values_only=True)), f"{name}:{sheet.title}"))
        finally:
            book.close()
    elif suffix in (".xml", ".nxml"):
        chunks, _ = jats_chunks(data)
        for chunk in chunks:
            chunk.locator = name + ":" + chunk.locator
    else:
        warnings.append(f"unsupported_attachment:{name}")
    return chunks, warnings


def zip_chunks(data, max_bytes=50_000_000):
    chunks, warnings = [], []
    total = 0
    with zipfile.ZipFile(BytesIO(data)) as archive:
        for item in archive.infolist():
            if item.is_dir():
                continue
            total += item.file_size
            if item.file_size > max_bytes or total > max_bytes * 4:
                warnings.append(f"attachment_archive_size_limit:{item.filename}")
                continue
            # Read in memory; never extract untrusted paths or executables to disk.
            try:
                parsed, errors = attachment_chunks(item.filename, archive.read(item), max_bytes)
                chunks.extend(parsed)
                warnings.extend(errors)
            except Exception as exc:
                warnings.append(f"attachment_parse_failed:{item.filename}:{type(exc).__name__}")
    return chunks, warnings


def relevant_chunk(chunk: Chunk):
    return bool(re.search(r"\b(?:IC|EC|CC)\s*[_\-]?\s*50\b|(?:inhibitory|effective|cytotoxic) concentration",
                          chunk.text, re.I))
