#!/usr/bin/env python3
"""
update_calkit.py -- Merge generated analysis/postpro stages into calkit.yaml.

Workflow
--------
1. In MATLAB, run ``analysis_to_calkit.m``.
   This writes ``calkit_stages.yaml`` containing the auto-generated
   analysis-* and postpro-* stages for every analysis class.

2. Run this script from the repository root::

       python mdocean/analysis/update_calkit.py

   The script merges the generated stages into ``calkit.yaml`` while:

   * Preserving all non-analysis/postpro stages (e.g. mermaid-*, make-drag-integral).
   * Preserving special inputs such as ``from_stage_outputs`` that calkit
     adds manually and that are not emitted by ``write_calkit_stage``.
   * Prompting the user before discarding any inputs or outputs that exist in
     the current ``calkit.yaml`` but are absent from the generated stages.
   * Inserting brand-new stages in alphabetical order among the existing
     analysis-*/postpro-* stages.

Usage
-----
    python mdocean/analysis/update_calkit.py [--calkit FILE] [--stages FILE] [--yes]

Options
-------
    --calkit FILE   Path to calkit.yaml          (default: calkit.yaml)
    --stages FILE   Path to calkit_stages.yaml   (default: calkit_stages.yaml)
    --yes           Keep all manual additions without prompting
"""

import sys
import argparse
from pathlib import Path

try:
    from ruamel.yaml import YAML
    from ruamel.yaml.comments import CommentedMap
except ImportError:
    sys.exit(
        "ERROR: ruamel.yaml is required.\n"
        "Install with:  pip install ruamel.yaml"
    )


# ---------------------------------------------------------------------------
# Helpers for YAML item inspection
# ---------------------------------------------------------------------------

def get_path(item):
    """Return the canonical file path from a YAML list item.

    Items are either plain strings (``- ./path/to/file``) or dicts with a
    ``path`` key (``- path: ./file\n  storage: git``).
    Returns ``None`` for other dict types (e.g. ``from_stage_outputs``).
    """
    if isinstance(item, str):
        return item
    if isinstance(item, dict) and "path" in item:
        return item["path"]
    return None


def is_special(item):
    """True for dict items that are NOT plain file paths.

    Examples: ``{from_stage_outputs: analysis-X}``.
    These are calkit-specific directives that ``write_calkit_stage`` does not
    emit and must always be preserved.
    """
    return isinstance(item, dict) and "path" not in item


# ---------------------------------------------------------------------------
# Merge logic
# ---------------------------------------------------------------------------

def ask_keep(kind, path, yes):
    """Prompt the user whether to keep a manual addition.  Default: keep."""
    if yes:
        print(f"    [keep] manual {kind}: {path!r}")
        return True
    ans = input(
        f"  Keep manual {kind} not in generated: {path!r}  [Y/n] "
    ).strip().lower()
    return ans != "n"


def merge_inputs(old_inputs, new_inputs, stage_key, yes):
    """Return merged inputs list for *stage_key*.

    Order: special items (``from_stage_outputs`` etc.) first, then
    generated inputs, then any manually-added path items the user wants
    to keep.
    """
    old_inputs = list(old_inputs or [])
    new_inputs = list(new_inputs or [])

    # Special dict items (e.g. from_stage_outputs) are always preserved.
    special = [x for x in old_inputs if is_special(x)]

    # Paths that the generated stages already cover.
    new_paths = {get_path(x) for x in new_inputs if get_path(x)}

    # Plain-path items in the old list that are absent from the generated list.
    manual = [
        x for x in old_inputs
        if not is_special(x) and get_path(x) not in new_paths
    ]

    kept = []
    if manual:
        print(
            f"  Stage {stage_key!r}: "
            f"{len(manual)} input(s) in calkit.yaml not in generated:"
        )
        for item in manual:
            if ask_keep("input", get_path(item), yes):
                kept.append(item)

    return special + new_inputs + kept


def merge_outputs(old_outputs, new_outputs, stage_key, yes):
    """Return merged outputs list for *stage_key*.

    Generated outputs come first; any manually-added outputs the user
    wishes to keep are appended at the end.
    """
    old_outputs = list(old_outputs or [])
    new_outputs = list(new_outputs or [])

    new_paths = {get_path(x) for x in new_outputs if get_path(x)}

    manual = [
        x for x in old_outputs
        if get_path(x) not in new_paths
    ]

    kept = []
    if manual:
        print(
            f"  Stage {stage_key!r}: "
            f"{len(manual)} output(s) in calkit.yaml not in generated:"
        )
        for item in manual:
            if ask_keep("output", get_path(item), yes):
                kept.append(item)

    return new_outputs + kept


# ---------------------------------------------------------------------------
# Insertion helpers
# ---------------------------------------------------------------------------

def analysis_name(stage_key):
    """Extract the analysis class name from a stage key.

    E.g. ``'analysis-AllFigCompare'`` → ``'AllFigCompare'``.
    """
    for prefix in ("analysis-", "postpro-"):
        if stage_key.startswith(prefix):
            return stage_key[len(prefix):]
    return stage_key


def find_insert_after(stages_map, new_key):
    """Return the key after which *new_key* should be inserted.

    New stages are placed in alphabetical order among the existing
    analysis-* / postpro-* stages.  For a given analysis class name,
    ``analysis-X`` always comes before ``postpro-X``.

    Returns ``None`` if the new stage should be inserted before all
    existing analysis/postpro stages.
    """
    new_name = analysis_name(new_key)
    new_is_postpro = new_key.startswith("postpro-")

    insert_after = None
    for key in stages_map:
        if not (key.startswith("analysis-") or key.startswith("postpro-")):
            continue
        name = analysis_name(key)
        is_postpro = key.startswith("postpro-")

        if name < new_name:
            insert_after = key
        elif name == new_name and (not is_postpro) and new_is_postpro:
            # analysis-X is the immediate predecessor of postpro-X
            insert_after = key

    return insert_after


def insert_stage(stages_map, insert_after_key, new_key, new_value):
    """Return a new CommentedMap with *new_key*/*new_value* inserted.

    The new entry is placed immediately after *insert_after_key* (or at
    the position of the first analysis/postpro stage if *insert_after_key*
    is ``None``).
    """
    items = list(stages_map.items())

    if insert_after_key is None:
        # Insert before the first analysis/postpro stage.
        idx = next(
            (
                i for i, (k, _) in enumerate(items)
                if k.startswith("analysis-") or k.startswith("postpro-")
            ),
            len(items),
        )
    else:
        idx = next(
            (i for i, (k, _) in enumerate(items) if k == insert_after_key),
            len(items) - 1,
        ) + 1

    new_map = CommentedMap()
    for i, (k, v) in enumerate(items):
        if i == idx:
            new_map[new_key] = new_value
        new_map[k] = v
    if idx >= len(items):
        new_map[new_key] = new_value

    return new_map


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--calkit",
        default="calkit.yaml",
        metavar="FILE",
        help="Path to calkit.yaml  (default: calkit.yaml)",
    )
    parser.add_argument(
        "--stages",
        default="calkit_stages.yaml",
        metavar="FILE",
        help="Path to generated stages file  (default: calkit_stages.yaml)",
    )
    parser.add_argument(
        "--yes",
        action="store_true",
        help="Keep all manual additions without prompting",
    )
    args = parser.parse_args()

    calkit_path = Path(args.calkit)
    stages_path = Path(args.stages)

    if not calkit_path.exists():
        sys.exit(f"ERROR: {calkit_path} not found.")
    if not stages_path.exists():
        sys.exit(
            f"ERROR: {stages_path} not found.\n"
            "Run analysis_to_calkit.m in MATLAB first."
        )

    # ------------------------------------------------------------------
    # Load files
    # ------------------------------------------------------------------

    # Round-trip loader preserves comments and original formatting for
    # parts of calkit.yaml that we do NOT touch.
    yaml_rt = YAML()
    yaml_rt.preserve_quotes = True
    yaml_rt.width = 4096  # prevent unwanted line wrapping

    print(f"Loading {calkit_path} ...")
    with open(calkit_path) as fh:
        calkit = yaml_rt.load(fh)

    # Safe (plain dict) load is fine for the generated stages file.
    yaml_safe = YAML(typ="safe")
    print(f"Loading {stages_path} ...")
    with open(stages_path) as fh:
        gen_stages = yaml_safe.load(fh)

    if not gen_stages:
        sys.exit(f"ERROR: {stages_path} is empty or could not be parsed.")

    # ------------------------------------------------------------------
    # Merge
    # ------------------------------------------------------------------

    stages = calkit["pipeline"]["stages"]

    for stage_key, gen_stage in gen_stages.items():
        if not (
            stage_key.startswith("analysis-")
            or stage_key.startswith("postpro-")
        ):
            print(f"Skipping non-analysis entry in generated file: {stage_key!r}")
            continue

        if stage_key in stages:
            print(f"Updating  {stage_key}")
            old = stages[stage_key]

            merged_in = merge_inputs(
                old.get("inputs"), gen_stage.get("inputs"), stage_key, args.yes
            )
            merged_out = merge_outputs(
                old.get("outputs"), gen_stage.get("outputs"), stage_key, args.yes
            )

            # Only replace inputs and outputs; preserve kind, command, etc.
            stages[stage_key]["inputs"] = merged_in
            stages[stage_key]["outputs"] = merged_out

        else:
            print(f"Inserting {stage_key}")
            after = find_insert_after(stages, stage_key)
            new_stages = insert_stage(stages, after, stage_key, gen_stage)
            calkit["pipeline"]["stages"] = new_stages
            stages = calkit["pipeline"]["stages"]

    # ------------------------------------------------------------------
    # Write back
    # ------------------------------------------------------------------

    print(f"\nWriting {calkit_path} ...")
    with open(calkit_path, "w") as fh:
        yaml_rt.dump(calkit, fh)

    print("Done.")


if __name__ == "__main__":
    main()
