#!/usr/bin/env python3
"""Merge the four source symbol glossaries and handle duplicate keys."""

import re
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
SOURCE_FILES = [
    ROOT / "pubs" / "applied-ocean-research-model" / "symbol-glossary-aor.tex",
    ROOT / "pubs" / "renewable-energy-mdo" / "symbol-glossary-re.tex",
    ROOT / "mdocean" / "simulation" / "modules" / "OpenFLASH" / "pubs" / "JFM" / "symbol-glossary-jfm.tex",
    # ROOT / "pubs" / "UMERC-2025-grid-value" / "symbol-glossary-umerc.tex",
]
OUTPUT_FILE = ROOT / "pubs" / "dissertation" / "symbol-glossary-dissertation.tex"


def is_content_line(line: str) -> bool:
    stripped = line.strip()
    return bool(stripped) and not stripped.startswith("%")


def get_symbol_key(line: str) -> str | None:
    match = re.search(r"\]\{([^}]*)\}\{", line)
    if match is None:
        return None
    return match.group(1)


def main() -> None:
    merged_lines: list[str] = []
    seen: set[str] = set()
    seen_keys: set[str] = set()

    for source_file in SOURCE_FILES:
        for line in source_file.read_text(encoding="utf-8").splitlines():
            if not is_content_line(line):
                continue
            if line in seen:
                continue
            key = get_symbol_key(line)
            if key is not None and key in seen_keys:
                merged_lines.append(f"% {line}")
                seen.add(line)
                continue
            seen.add(line)
            if key is not None:
                seen_keys.add(key)
            merged_lines.append(line)

    OUTPUT_FILE.write_text("\n".join(merged_lines) + "\n", encoding="utf-8")
    print(f"Wrote {OUTPUT_FILE.relative_to(ROOT)} from {len(SOURCE_FILES)} sources.")


if __name__ == "__main__":
    main()