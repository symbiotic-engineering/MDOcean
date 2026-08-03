#!/usr/bin/env python3
"""Update AOR symbol glossary and replace math symbols with \\gls references."""

from __future__ import annotations

import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
PUBS_DIR = ROOT / "pubs" / "applied-ocean-research-model"
GLOSSARY_FILE = ROOT / "pubs" / "shared" / "symbol-glossary-aor.tex"
GLOSSARY_SEED_FILE = ROOT / "pubs" / "shared" / "symbol-glossary-aor-seed.tex"
PAPERLINT_FILE = ROOT / "dev" / "latex" / "Paper-Linter" / "paperlint.py"
PAPERLINT_SETTINGS = ROOT / "dev" / "latex" / "Paper-Linter" / "settings.txt"
IGNORE_FILE = PUBS_DIR / "numeric-results.tex"
MAIN_FILE = PUBS_DIR / "main.tex"


def run_paperlint(extra_args: list[str]) -> None:
    subprocess.run(
        [
            "python3",
            str(PAPERLINT_FILE),
            str(MAIN_FILE),
            "--ignore",
            str(IGNORE_FILE),
            "--settings",
            str(PAPERLINT_SETTINGS),
            "--output",
            str(ROOT / "dev" / "latex" / "Paper-Linter" / "lint-aor.txt"),
            "--symbol-glossary-output",
            str(GLOSSARY_FILE),
            "--symbol-glossary-seed",
            str(GLOSSARY_SEED_FILE),
            *extra_args,
        ],
        cwd=ROOT,
        check=True,
    )


def main() -> None:
    run_paperlint([])
    run_paperlint(["--replace-glossary-refs"])


if __name__ == "__main__":
    main()
