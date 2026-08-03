#!/usr/bin/env python3
"""Replace manual glossary-symbol text in pubs/**/*.tex with \gls{...}.

The script parses pubs/shared/symbol-glossary-shared.tex, builds a mapping from
ensuremath bodies to glossary keys, then scans tex files under pubs/ and replaces
matching symbol text only inside LaTeX math spans with the corresponding
glossary reference.
"""

from __future__ import annotations

import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
GLOSSARY_FILE = ROOT / "pubs" / "shared" / "symbol-glossary-shared.tex"
PUBS_DIR = ROOT / "pubs"
TARGET_FILES = [
    PUBS_DIR / "applied-ocean-research-model" / "sections" / "module-details.tex",
    PUBS_DIR / "applied-ocean-research-model" / "sections" / "meem-appendix.tex",
    PUBS_DIR / "applied-ocean-research-model" / "sections" / "other-appendices.tex",
]

ENTRY_RE = re.compile(
    r"\\glsxtrnewsymbol\[(?P<attrs>[^\]]*)\]\{(?P<key>[^}]+)\}\{\\ensuremath\{(?P<body>.*)\}\}"
)
GLS_COMMAND_RE = re.compile(r"\\gls\{[^{}]*\}")


def parse_entries() -> list[tuple[str, str]]:
    entries: list[tuple[str, str]] = []
    for line in GLOSSARY_FILE.read_text().splitlines():
        stripped = line.strip()
        if not stripped.startswith("\\glsxtrnewsymbol["):
            continue
        match = ENTRY_RE.fullmatch(stripped)
        if not match:
            continue
        key = match.group("key")
        body = match.group("body")
        entries.append((key, body))
    return entries


def target_files() -> list[Path]:
    return [path for path in TARGET_FILES if path.exists()]


def find_math_spans(text: str) -> list[tuple[int, int]]:
    spans: list[tuple[int, int]] = []
    index = 0
    length = len(text)

    while index < length:
        char = text[index]
        if char == "\\" and index + 1 < length and text[index + 1] in "([":
            end_token = "\\)" if text[index + 1] == "(" else "\\]"
            start = index + 2
            end = text.find(end_token, start)
            if end != -1:
                spans.append((start, end))
                index = end + 2
                continue
        if char == "$":
            if index + 1 < length and text[index + 1] == "$":
                start = index + 2
                end = text.find("$$", start)
                if end != -1:
                    spans.append((start, end))
                    index = end + 2
                    continue
            else:
                start = index + 1
                end = start
                while end < length:
                    if text[end] == "$" and text[end - 1] != "\\":
                        spans.append((start, end))
                        index = end + 1
                        break
                    end += 1
                else:
                    break
                continue
        index += 1

    return spans


def replace_in_math_spans(text: str, mapping: list[tuple[str, str]]) -> tuple[str, int]:
    spans = find_math_spans(text)
    if not spans:
        return text, 0

    pieces: list[str] = []
    last_index = 0
    total = 0

    for start, end in spans:
        pieces.append(text[last_index:start])
        span_text, count = replace_in_math_span(text[start:end], mapping)
        total += count
        pieces.append(span_text)
        last_index = end

    pieces.append(text[last_index:])
    return "".join(pieces), total


def replace_in_math_span(span_text: str, mapping: list[tuple[str, str]]) -> tuple[str, int]:
    pieces: list[str] = []
    index = 0
    total = 0

    while index < len(span_text):
        if span_text.startswith("\\gls{", index):
            match = GLS_COMMAND_RE.match(span_text, index)
            if match:
                pieces.append(match.group(0))
                index = match.end()
                continue

        for key, body in mapping:
            if span_text.startswith(body, index) and body_matches_boundaries(span_text, index, body):
                pieces.append(f"\\gls{{{key}}}")
                index += len(body)
                total += 1
                break
        else:
            pieces.append(span_text[index])
            index += 1

    return "".join(pieces), total


def body_matches_boundaries(span_text: str, index: int, body: str) -> bool:
    if body.startswith("\\"):
        return True

    previous_char = span_text[index - 1] if index > 0 else ""
    next_index = index + len(body)
    next_char = span_text[next_index] if next_index < len(span_text) else ""
    if previous_char and (previous_char.isalnum() or previous_char == "\\"):
        return False
    if next_char and (next_char.isalnum() or next_char == "\\"):
        return False
    return True


def main() -> None:
    entries = parse_entries()
    # Replace longest bodies first to avoid partial overlaps such as \vec{F}
    # inside \vec{F}_{rad}.
    entries.sort(key=lambda item: len(item[1]), reverse=True)

    files = target_files()
    changed_files = 0
    total_replacements = 0

    for path in files:
        original = path.read_text()
        updated, count = replace_in_math_spans(original, entries)
        updated = "\n".join(line.rstrip() for line in updated.splitlines())
        if original.endswith("\n"):
            updated += "\n"
        if count and updated != original:
            path.write_text(updated)
            changed_files += 1
            total_replacements += count
            print(f"{path}: {count} replacements")

    print(f"Processed {len(files)} files")
    print(f"Updated {changed_files} files")
    print(f"Total replacements: {total_replacements}")


if __name__ == "__main__":
    main()
