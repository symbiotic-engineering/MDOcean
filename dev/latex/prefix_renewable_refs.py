#!/usr/bin/env python3
"""Prefix references in the renewable-energy chapter to match remdo labels."""

import re
from pathlib import Path

ROOT = Path("/home/becca/Documents/git/MDOcean")
TARGET_DIR = ROOT / "pubs" / "renewable-energy-mdo"

REF_CMDS = ["ref", "eqref", "cref", "Cref", "autoref", "pageref", "nameref"]


def main() -> None:
    labels = set()
    tex_files = []
    for path in TARGET_DIR.rglob("*.tex"):
        if "/aux/" in str(path) or path.name.endswith(".bak"):
            continue
        tex_files.append(path)
        text = path.read_text(encoding="utf-8")
        for label in re.findall(r"\\label\{remdo:([^}]+)\}", text):
            labels.add(label)

    if not labels:
        print("No remdo labels found.")
        return

    pattern = re.compile(r"\\(" + "|".join(REF_CMDS) + r")\{([^}]*)\}")
    changed = 0

    for path in tex_files:
        original = path.read_text(encoding="utf-8")

        def repl(match: re.Match[str]) -> str:
            command = match.group(1)
            body = match.group(2)
            parts = [part.strip() for part in body.split(",")]
            updated_parts = []
            for part in parts:
                if part.startswith("remdo:"):
                    updated_parts.append(part)
                elif part in labels:
                    updated_parts.append(f"remdo:{part}")
                else:
                    updated_parts.append(part)
            return f"\\{command}{{{','.join(updated_parts)}}}"

        updated = pattern.sub(repl, original)
        if updated != original:
            path.write_text(updated, encoding="utf-8")
            changed += 1
            print(f"Updated {path}")

    print(f"Done. Files changed: {changed}")


if __name__ == "__main__":
    main()
