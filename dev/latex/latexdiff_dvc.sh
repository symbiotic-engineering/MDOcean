#!/usr/bin/env bash
set -euo pipefail

# Usage:
#
#   ./dvc_latexdiff.sh <old_git_revision> <main_latex_file.tex>
#
# Example:
#
#   ./dvc_latexdiff.sh HEAD~1 pubs/applied-ocean-research-model/main.tex

OLD_REV="$1"
TEX_FILE="$2"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
FILTER="$SCRIPT_DIR/dvc_latexdiff_filter.py"

# Absolute repository root.
REPO_ROOT="$(git rev-parse --show-toplevel)"

# ----------------------------------------------------------------------
# Fetch the DVC objects for the current revision.
# ----------------------------------------------------------------------

echo "Fetching current DVC objects..."

if ! dvc fetch; then
    echo "WARNING: current dvc fetch failed; continuing." >&2
fi

# ----------------------------------------------------------------------
# Fetch the DVC objects referenced by the old Git revision.
#
# No Git worktree and no DVC checkout are required.
# ----------------------------------------------------------------------

echo "Fetching DVC objects for revision $OLD_REV..."

if ! dvc fetch --rev "$OLD_REV"; then
    echo "WARNING: dvc fetch --rev $OLD_REV failed; continuing." >&2
fi

# ----------------------------------------------------------------------
# DVC cache directory.
# ----------------------------------------------------------------------

DVC_CACHE_DIR="$(dvc cache dir --show)"

# Make it absolute in case DVC returns a relative path.
if [[ "$DVC_CACHE_DIR" != /* ]]; then
    DVC_CACHE_DIR="$REPO_ROOT/$DVC_CACHE_DIR"
fi

DVC_CACHE_DIR="$(cd "$DVC_CACHE_DIR" && pwd)"

# ----------------------------------------------------------------------
# Tell the filter:
#
#   - repository root
#   - old Git revision
#   - DVC cache location
#
# The filter obtains the old dvc.lock with:
#
#   git show OLD_REV:dvc.lock
#
# and reads the current dvc.lock directly.
# ----------------------------------------------------------------------

export DVC_REPO_ROOT="$REPO_ROOT"
export DVC_OLD_REV="$OLD_REV"
export DVC_CACHE_DIR

# ----------------------------------------------------------------------
# Run latexdiff-vc.
#
# --flatten makes latexdiff-vc retrieve the old revision into a
# temporary directory and passes --flatten to latexdiff.
#
# The filter is therefore applied to the old and new source files
# independently.
# ----------------------------------------------------------------------

latexdiff-vc \
    --git \
    --flatten \
    --force \
    --graphics-markup=both \
    # --filter-script="$FILTER" \
    # --ignore-filter-stderr \
    -r "$OLD_REV" \
    "$TEX_FILE"