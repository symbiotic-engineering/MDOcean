#!/usr/bin/env python3
"""Prefix refs and labels in the OpenFLASH appendix copy to avoid collisions."""

import re
from pathlib import Path

FILE = Path("/home/becca/Documents/git/MDOcean-worktrees/review/mdocean/simulation/modules/OpenFLASH/pubs/JFM/jfm-appendix.tex")
PREFIX = "jfm:"
REF_CMDS = ["ref", "eqref", "cref", "Cref", "autoref", "pageref", "nameref"]


def main() -> None:
    text = FILE.read_text(encoding="utf-8")
    labels = sorted(set(re.findall(r"\\label\{([^}]+)\}", text)), key=len, reverse=True)
    labels = [label for label in labels if not label.startswith(PREFIX)]
    if not labels:
        print("No labels to prefix.")
        return

    mapping = {label: f"{PREFIX}{label}" for label in labels}
    pattern = re.compile(r"\\(" + "|".join(REF_CMDS) + r")\{([^}]*)\}")

    def repl(match: re.Match[str]) -> str:
        command = match.group(1)
        body = match.group(2)
        parts = [part.strip() for part in body.split(",")]
        updated = [mapping.get(part, part) for part in parts]
        return f"\\{command}{{{','.join(updated)}}}"

    updated = pattern.sub(repl, text)
    for old, new in mapping.items():
        updated = updated.replace(f"\\label{{{old}}}", f"\\label{{{new}}}")

    if updated != text:
        FILE.write_text(updated, encoding="utf-8")
        print(f"Updated {FILE} ({len(mapping)} labels)")


if __name__ == "__main__":
    main()
