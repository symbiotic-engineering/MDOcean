#!/usr/bin/env python3
"""Replace math-mode symbols with glossary references in AOR TeX files."""

from __future__ import annotations

import re
import subprocess
from pathlib import Path

from math_symbol_parser import canonicalize_symbol, find_math_symbol_matches, load_glossary_symbol_map, symbol_key

ROOT = Path(__file__).resolve().parents[2]
PUBS_DIR = ROOT / "pubs" / "applied-ocean-research-model"
GLOSSARY_FILE = ROOT / "pubs" / "shared" / "symbol-glossary-aor.tex"
GLOSSARY_SEED_FILE = ROOT / "pubs" / "shared" / "symbol-glossary-aor-seed.tex"
PAPERLINT_FILE = ROOT / "dev" / "latex" / "Paper-Linter" / "paperlint.py"
PAPERLINT_SETTINGS = ROOT / "dev" / "latex" / "Paper-Linter" / "settings.txt"
IGNORE_FILE = PUBS_DIR / "numeric-results.tex"
MAIN_FILE = PUBS_DIR / "main.tex"


def target_files() -> list[Path]:
    files = []
    for path in sorted(PUBS_DIR.rglob("*.tex")):
        if path == IGNORE_FILE:
            continue
        files.append(path)
    return files


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


def replace_symbols_in_span(span_text: str, glossary_symbol_map: dict[str, str]):
    replacements: list[tuple[int, int, str]] = []
    missing_symbols: set[str] = set()
    for start, end, symbol in find_math_symbol_matches(span_text):
        label = glossary_symbol_map.get(symbol_key(canonicalize_symbol(symbol)))
        if label is None:
            missing_symbols.add(symbol)
            continue
        replacements.append((start, end, label))

    if not replacements:
        return span_text, 0, missing_symbols

    updated = span_text
    for start, end, label in sorted(replacements, key=lambda item: item[0], reverse=True):
        updated = updated[:start] + f"\\gls{{{label}}}" + updated[end:]
    return updated, len(replacements), missing_symbols


def replace_symbols_in_text(text: str, glossary_symbol_map: dict[str, str]):
    spans = find_math_spans(text)
    if not spans:
        return text, 0, set()

    pieces = []
    last = 0
    total = 0
    missing_symbols: set[str] = set()
    for start, end in spans:
        pieces.append(text[last:start])
        updated, count, missing = replace_symbols_in_span(text[start:end], glossary_symbol_map)
        pieces.append(updated)
        total += count
        missing_symbols.update(missing)
        last = end
    pieces.append(text[last:])
    return "".join(pieces), total, missing_symbols


def update_glossary_with_paperlint() -> None:
    subprocess.run(
        [
            "python3",
            str(PAPERLINT_FILE),
            str(MAIN_FILE),
            "--ignore",
            str(IGNORE_FILE),
            "--settings",
            str(PAPERLINT_SETTINGS),
            "--symbol-glossary-output",
            str(GLOSSARY_FILE),
            "--symbol-glossary-seed",
            str(GLOSSARY_SEED_FILE),
        ],
        cwd=ROOT,
        check=True,
    )


def process_files(files: list[Path], glossary_symbol_map: dict[str, str]):
    changed_files = 0
    total_replacements = 0
    missing_symbols: set[str] = set()
    for path in files:
        original = path.read_text()
        updated, count, missing = replace_symbols_in_text(original, glossary_symbol_map)
        missing_symbols.update(missing)
        if count and updated != original:
            path.write_text(updated)
            changed_files += 1
            total_replacements += count
            print(f"{path}: {count} replacements")
    return changed_files, total_replacements, missing_symbols


def main() -> None:
    files = target_files()
    glossary_symbol_map = load_glossary_symbol_map(GLOSSARY_FILE)
    changed_files, total_replacements, missing_symbols = process_files(files, glossary_symbol_map)

    if missing_symbols:
        print(f"Found {len(missing_symbols)} missing symbols; updating glossary with paperlint.")
        update_glossary_with_paperlint()
        glossary_symbol_map = load_glossary_symbol_map(GLOSSARY_FILE)
        changed_files_2, total_replacements_2, missing_symbols = process_files(files, glossary_symbol_map)
        changed_files += changed_files_2
        total_replacements += total_replacements_2

    if missing_symbols:
        print(f"Unresolved symbols after glossary update: {len(missing_symbols)}")
        for symbol in sorted(missing_symbols):
            print(f"  - {symbol}")

    print(f"Processed {len(files)} files")
    print(f"Updated {changed_files} files")
    print(f"Total replacements: {total_replacements}")


if __name__ == "__main__":
    main()
