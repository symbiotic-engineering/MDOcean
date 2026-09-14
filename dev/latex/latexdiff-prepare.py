#!/usr/bin/env python3

import json
import os
import re
import sys
from pathlib import Path


REPO_ROOT = Path(os.environ["DVC_REPO_ROOT"]).resolve()
CACHE_DIR = Path(os.environ["DVC_CACHE_DIR"]).resolve()
MAPPING_FILE = Path(os.environ["DVC_DEPS_JSON"])


# ----------------------------------------------------------------------
# DVC cache
# ----------------------------------------------------------------------

def cache_object(md5):
    for path in (
        CACHE_DIR / md5[:2] / md5[2:],
        CACHE_DIR / "files" / "md5" / md5[:2] / md5[2:],
        CACHE_DIR / "files" / "md5" / md5[:2] / md5[2:].with_suffix(".dir"),
    ):
        if path.is_file():
            return path

    return None


# ----------------------------------------------------------------------
# Directory -> .dir MD5 mapping
# ----------------------------------------------------------------------

with MAPPING_FILE.open() as f:
    DVC_DIRS = json.load(f)


# ----------------------------------------------------------------------
# Expand .dir manifests into:
#
#     repo-relative file -> file MD5
# ----------------------------------------------------------------------

FILE_HASHES = {}

for directory, dir_md5 in DVC_DIRS.items():

    # dvc.lock stores directory hashes as "...dir"
    if dir_md5.endswith(".dir"):
        dir_md5 = dir_md5[:-4]

    manifest = cache_object(dir_md5)

    if manifest is None:
        print(
            f"WARNING: DVC .dir object not found: "
            f"{directory}: {dir_md5}",
            file=sys.stderr,
        )
        continue

    try:
        entries = json.loads(manifest.read_text())
    except Exception as exc:
        print(
            f"WARNING: could not read DVC .dir object "
            f"{manifest}: {exc}",
            file=sys.stderr,
        )
        continue

    for entry in entries:
        relpath = entry.get("relpath")
        md5 = entry.get("md5")

        if relpath and md5:
            FILE_HASHES[
                (Path(directory) / relpath).as_posix()
            ] = md5


# ----------------------------------------------------------------------
# Resolve an includegraphics path
# ----------------------------------------------------------------------

def repo_relative_path(filename):
    path = Path(filename)

    if path.is_absolute():
        try:
            return path.resolve().relative_to(REPO_ROOT).as_posix()
        except ValueError:
            return None

    try:
        return (
            REPO_ROOT / path
        ).resolve().relative_to(REPO_ROOT).as_posix()

    except ValueError:
        return None


# ----------------------------------------------------------------------
# \includegraphics
# ----------------------------------------------------------------------

INCLUDEGRAPHICS = re.compile(
    r"""
    (\\includegraphics
        \s*
        (?:\[[^\]]*\])?
        \s*
        \{)
    ([^{}]+)
    (\})
    """,
    re.VERBOSE,
)


def replace_graphic(match):
    filename = match.group(2).strip()

    # Don't touch dynamically constructed paths.
    if filename.startswith("\\"):
        return match.group(0)

    repo_path = repo_relative_path(filename)

    if repo_path is None:
        return match.group(0)

    md5 = FILE_HASHES.get(repo_path)

    if md5 is None:
        return match.group(0)

    cached = cache_object(md5)

    if cached is None:
        print(
            f"WARNING: DVC cache object missing: "
            f"{repo_path}: {md5}",
            file=sys.stderr,
        )
        return match.group(0)

    return (
        match.group(1)
        + str(cached)
        + match.group(3)
    )


# ----------------------------------------------------------------------
# stdin -> stdout
# ----------------------------------------------------------------------

text = sys.stdin.read()

text = INCLUDEGRAPHICS.sub(
    replace_graphic,
    text,
)

sys.stdout.write(text)