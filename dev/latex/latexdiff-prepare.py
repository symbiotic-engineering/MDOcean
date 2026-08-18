#!/usr/bin/env python3

import json
import os
import re
import subprocess
import sys
from pathlib import Path


# ======================================================================
# Configuration
# ======================================================================

# The temporary git-latexdiff tree in which this script is running.
WORK_ROOT = Path.cwd().resolve()

# The real repository containing dvc.lock and .dvc/cache.
DVC_REPO_ROOT = Path(
    os.environ["DVC_REPO_ROOT"]
).resolve()

CACHE_DIR = Path(
    os.environ.get(
        "DVC_CACHE_DIR",
        DVC_REPO_ROOT / ".dvc" / "cache",
    )
).resolve()


# ======================================================================
# DVC cache helpers
# ======================================================================

def cache_object(md5):
    """
    Return the path to a DVC cache object for an MD5 hash.

    Supports both:

        .dvc/cache/ab/cdef...
        .dvc/cache/files/md5/ab/cdef...
    """

    candidates = [
        CACHE_DIR / md5[:2] / md5[2:],
        CACHE_DIR / "files" / "md5" / md5[:2] / md5[2:],
    ]

    for candidate in candidates:
        if candidate.is_file():
            return candidate

    return None


# ======================================================================
# Read dvc.lock using yq
# ======================================================================

def build_file_hash_map():
    """
    Build:

        repository-relative file path -> DVC MD5

    from the dvc.lock belonging to this working tree.

    yq is used to parse YAML; no PyYAML dependency is required.
    """

    lockfile = WORK_ROOT / "dvc.lock"

    if not lockfile.is_file():
        print(
            f"WARNING: {lockfile} not found",
            file=sys.stderr,
        )
        return {}

    # Extract all directory outputs as:
    #
    #   <directory>\t<directory-md5.dir>
    #
    result = subprocess.run(
        [
            "yq",
            "-r",
            """
            .stages[]
            | .outs[]?
            | select(.path != null and .md5 != null)
            | select(.md5 | endswith(".dir"))
            | [.path, .md5]
            | @tsv
            """,
            str(lockfile),
        ],
        cwd=WORK_ROOT,
        capture_output=True,
        text=True,
    )

    if result.returncode != 0:
        print(
            f"WARNING: yq failed to read {lockfile}: "
            f"{result.stderr.strip()}",
            file=sys.stderr,
        )
        return {}

    file_hashes = {}

    for line in result.stdout.splitlines():

        if not line.strip():
            continue

        directory, directory_md5 = line.split("\t", 1)

        # DVC stores directory hashes as "<md5>.dir".
        dir_md5 = directory_md5.removesuffix(".dir")

        manifest_path = cache_object(dir_md5)

        if manifest_path is None:
            print(
                f"WARNING: DVC .dir object not found: "
                f"{dir_md5}",
                file=sys.stderr,
            )
            continue

        try:
            manifest = json.loads(
                manifest_path.read_text()
            )
        except Exception as exc:
            print(
                f"WARNING: could not read DVC .dir object "
                f"{manifest_path}: {exc}",
                file=sys.stderr,
            )
            continue

        for entry in manifest:

            relpath = entry.get("relpath")
            md5 = entry.get("md5")

            if not relpath or not md5:
                continue

            repo_path = (
                Path(directory) / relpath
            ).as_posix()

            file_hashes[repo_path] = md5

    return file_hashes


HASHES = build_file_hash_map()


# ======================================================================
# Figure-path normalization
# ======================================================================

def repo_relative_figure_path(filename):
    """
    Convert an \\includegraphics filename into the repository-relative
    path used by the DVC .dir manifest.
    """

    path = Path(filename)

    if path.is_absolute():
        try:
            return path.resolve().relative_to(
                DVC_REPO_ROOT
            ).as_posix()
        except ValueError:
            return None

    # The LaTeX source is in the temporary tree, but the repository
    # structure is the same as the real repository.
    #
    # Find the path relative to the temporary tree.
    try:
        return (
            WORK_ROOT / path
        ).resolve().relative_to(
            WORK_ROOT
        ).as_posix()
    except ValueError:
        return None


# ======================================================================
# \includegraphics rewriting
# ======================================================================

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
    prefix = match.group(1)
    filename = match.group(2).strip()
    suffix = match.group(3)

    # Don't try to resolve dynamically generated paths.
    if filename.startswith("\\"):
        return match.group(0)

    repo_path = repo_relative_figure_path(filename)

    if repo_path is None:
        return match.group(0)

    md5 = HASHES.get(repo_path)

    if md5 is None:
        # Not a DVC-tracked file.
        return match.group(0)

    cached = cache_object(md5)

    if cached is None:
        print(
            f"WARNING: DVC cache object missing for "
            f"{repo_path}: {md5}",
            file=sys.stderr,
        )
        return match.group(0)

    print(
        f"DVC: {repo_path} -> {md5}",
        file=sys.stderr,
    )

    return (
        prefix
        + str(cached)
        + suffix
    )


# ======================================================================
# Process TeX files in the temporary tree
# ======================================================================

for texfile in WORK_ROOT.rglob("*.tex"):

    if ".git" in texfile.parts or ".dvc" in texfile.parts:
        continue

    try:
        text = texfile.read_text()
    except UnicodeDecodeError:
        continue

    new_text = INCLUDEGRAPHICS.sub(
        replace_graphic,
        text,
    )

    if new_text != text:
        texfile.write_text(new_text)

        print(
            f"Updated {texfile}",
            file=sys.stderr,
        )