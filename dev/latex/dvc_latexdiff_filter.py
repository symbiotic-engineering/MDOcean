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

REPO_ROOT = Path(
    os.environ["DVC_REPO_ROOT"]
).resolve()

OLD_REV = os.environ["DVC_OLD_REV"]

CACHE_DIR = Path(
    os.environ["DVC_CACHE_DIR"]
).resolve()


# ======================================================================
# DVC cache helpers
# ======================================================================

def cache_object(md5):
    """
    Return the path to the DVC cache object for an MD5 hash.

    Supports both:
        .dvc/cache/ab/cdef...
    and:
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
# Read dvc.lock
# ======================================================================

def load_yaml(text):
    """
    Load YAML using PyYAML.

    dvc.lock is YAML, so this avoids requiring DVC itself to parse
    the file.
    """

    import yaml

    return yaml.safe_load(text) or {}


def current_dvc_lock():
    """
    Read the current working-tree dvc.lock.
    """

    lockfile = REPO_ROOT / "dvc.lock"

    if not lockfile.is_file():
        return {}

    return load_yaml(lockfile.read_text())


def old_dvc_lock():
    """
    Read dvc.lock from the requested historical Git revision.
    """

    result = subprocess.run(
        [
            "git",
            "show",
            f"{OLD_REV}:dvc.lock",
        ],
        cwd=REPO_ROOT,
        capture_output=True,
        text=True,
    )

    if result.returncode != 0:
        print(
            f"WARNING: could not read dvc.lock from {OLD_REV}",
            file=sys.stderr,
        )
        return {}

    return load_yaml(result.stdout)


# ======================================================================
# .dir manifest handling
# ======================================================================

def find_dir_cache_object(dir_md5):
    """
    Locate the DVC .dir cache object.

    The .dir suffix is stripped before this function is called.
    """

    path = cache_object(dir_md5)

    if path is None:
        print(
            f"WARNING: DVC .dir object not found: {dir_md5}",
            file=sys.stderr,
        )

    return path


def load_dir_manifest(dir_md5):
    """
    Load a DVC .dir manifest.

    Returns:
        [
            {
                "md5": "...",
                "relpath": "foo.pdf"
            },
            ...
        ]
    """

    path = find_dir_cache_object(dir_md5)

    if path is None:
        return []

    try:
        return json.loads(path.read_text())
    except Exception as exc:
        print(
            f"WARNING: could not read DVC .dir object "
            f"{path}: {exc}",
            file=sys.stderr,
        )
        return []


# ======================================================================
# Build a mapping:
#
#     repo-relative path -> MD5
#
# for one dvc.lock.
# ======================================================================

def build_file_hash_map(lock):
    result = {}

    for stage in lock.get("stages", {}).values():

        for output in stage.get("outs", []):

            if not isinstance(output, dict):
                continue

            directory = output.get("path")
            directory_md5 = output.get("md5")

            if not directory or not directory_md5:
                continue

            # We only care about directory outputs.
            if not directory_md5.endswith(".dir"):
                continue

            dir_md5 = directory_md5[:-4]

            manifest = load_dir_manifest(dir_md5)

            for entry in manifest:

                relpath = entry.get("relpath")
                md5 = entry.get("md5")

                if not relpath or not md5:
                    continue

                repo_path = (
                    Path(directory) / relpath
                ).as_posix()

                result[repo_path] = md5

    return result


# Build both maps once when the filter process starts.
OLD_HASHES = build_file_hash_map(
    old_dvc_lock()
)

NEW_HASHES = build_file_hash_map(
    current_dvc_lock()
)


# ======================================================================
# Figure-path normalization
# ======================================================================

def repo_relative_figure_path(filename):
    """
    Convert an \\includegraphics filename to the repository-relative
    path used in dvc.lock.

    For ordinary LaTeX such as:

        \\includegraphics{pubs/.../foo.pdf}

    this simply normalizes the path.

    For absolute paths, strip the repository root if possible.
    """

    path = Path(filename)

    if path.is_absolute():

        try:
            return path.resolve().relative_to(
                REPO_ROOT
            ).as_posix()

        except ValueError:
            return None

    # Paths in the TeX source are normally relative to the
    # repository/project directory.
    try:
        absolute = (
            REPO_ROOT / path
        ).resolve()

        return absolute.relative_to(
            REPO_ROOT
        ).as_posix()

    except ValueError:
        return None


# ======================================================================
# \\includegraphics filter
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

    # Don't try to resolve dynamically generated LaTeX paths.
    if filename.startswith("\\"):
        return match.group(0)

    repo_path = repo_relative_figure_path(filename)

    if repo_path is None:
        return match.group(0)

    # --------------------------------------------------------------
    # Determine the MD5.
    #
    # If old and new hashes are equal, either is fine.
    #
    # If they differ, we need to determine which revision's source
    # is currently being filtered.
    # --------------------------------------------------------------

    old_md5 = OLD_HASHES.get(repo_path)
    new_md5 = NEW_HASHES.get(repo_path)

    if old_md5 is None and new_md5 is None:
        # Not a DVC-tracked file.
        return match.group(0)

    if old_md5 == new_md5:
        md5 = old_md5

    else:
        # ----------------------------------------------------------
        # IMPORTANT:
        #
        # latexdiff's filter-script interface does not explicitly
        # identify the old/new input. With --flatten, however,
        # latexdiff-vc processes the old and new retrieved files from
        # separate temporary directories.
        #
        # We therefore use the file's actual existence when possible.
        # ----------------------------------------------------------

        actual = Path(filename)

        old_exists = False
        new_exists = False

        if actual.is_absolute():
            old_exists = actual.is_file()

        else:
            new_exists = (
                REPO_ROOT / actual
            ).is_file()

        # If this is an ordinary current-worktree path, use new.
        if new_exists:
            md5 = new_md5

        elif old_exists:
            #
            # This case can occur when latexdiff-vc/flatten supplies
            # an absolute path into its temporary old checkout.
            #
            # The repo-relative path has already been extracted above.
            #
            md5 = old_md5

        else:
            #
            # We cannot safely identify the revision.
            #
            # Do not alter the figure rather than generating an
            # incorrect diff.
            #
            return match.group(0)

    if md5 is None:
        return match.group(0)

    # --------------------------------------------------------------
    # Replace the figure with the DVC cache object.
    # --------------------------------------------------------------

    cached = cache_object(md5)

    if cached is None:
        print(
            f"WARNING: DVC cache object missing for "
            f"{repo_path}: {md5}",
            file=sys.stderr,
        )

        return match.group(0)

    return (
        prefix
        + str(cached)
        + suffix
    )


# ======================================================================
# Process stdin -> stdout
# ======================================================================

text = sys.stdin.read()

text = INCLUDEGRAPHICS.sub(
    replace_graphic,
    text,
)

sys.stdout.write(text)