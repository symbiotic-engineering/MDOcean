#!/usr/bin/env python3
"""Merge the four source symbol glossaries and deduplicate exact lines."""

from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
SOURCE_FILES = [
    ROOT / "pubs" / "shared" / "symbol-glossary-aor.tex",
    ROOT / "pubs" / "shared" / "symbol-glossary-re.tex",
    ROOT / "pubs" / "shared" / "symbol-glossary-jfm.tex",
    ROOT / "pubs" / "shared" / "symbol-glossary-umerc.tex",
]
OUTPUT_FILE = ROOT / "pubs" / "shared" / "symbol-glossary-shared.tex"


def is_content_line(line: str) -> bool:
    stripped = line.strip()
    return bool(stripped) and not stripped.startswith("%")


def main() -> None:
    merged_lines: list[str] = []
    seen: set[str] = set()

    for source_file in SOURCE_FILES:
        for line in source_file.read_text(encoding="utf-8").splitlines():
            if not is_content_line(line):
                continue
            if line in seen:
                continue
            seen.add(line)
            merged_lines.append(line)

    OUTPUT_FILE.write_text("\n".join(merged_lines) + "\n", encoding="utf-8")
    print(f"Wrote {OUTPUT_FILE.relative_to(ROOT)} from {len(SOURCE_FILES)} sources.")


if __name__ == "__main__":
    main()