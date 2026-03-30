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
   * Overwriting all plain-path inputs and outputs with the generated values
     so that stale removed dependencies do not remain.
   * Inserting brand-new stages in alphabetical order among the existing
     analysis-*/postpro-* stages.

Usage
-----
    python mdocean/analysis/update_calkit.py [--calkit FILE] [--stages FILE]

Options
-------
    --calkit FILE   Path to calkit.yaml          (default: calkit.yaml)
    --stages FILE   Path to calkit_stages.yaml   (default: calkit_stages.yaml)
"""

import sys
import argparse
from pathlib import Path

try:
    from ruamel.yaml import YAML
    from ruamel.yaml.comments import CommentedMap, CommentedSeq
except ImportError:
    sys.exit(
        "ERROR: ruamel.yaml is required.\n"
        "Install with:  pip install ruamel.yaml"
    )


# ---------------------------------------------------------------------------
# Merge logic
# ---------------------------------------------------------------------------

def _to_commented_map(d):
    cm = CommentedMap()
    for k, v in d.items():
        cm[k] = v
    return cm


def _to_commented_seq(seq):
    cs = CommentedSeq()
    for x in seq:
        if isinstance(x, dict):
            cs.append(_to_commented_map(x))
        else:
            cs.append(x)
    return cs


def merge_inputs(old_inputs, new_inputs):
    """Return merged inputs list.

    Generated inputs completely replace the old inputs (overwrite behaviour)
    so that stale removed dependencies do not remain.
    """
    new_inputs = list(new_inputs or [])

    seq = CommentedSeq()
    for x in new_inputs:
        if isinstance(x, dict):
            seq.append(_to_commented_map(x))
        else:
            seq.append(x)

    return seq


def merge_outputs(old_outputs, new_outputs):
    """Return merged outputs list.

    Generated outputs completely replace the old outputs (overwrite
    behaviour) so that stale removed outputs do not remain.
    """
    new_outputs = list(new_outputs or [])

    seq = CommentedSeq()
    for x in new_outputs:
        if isinstance(x, dict):
            seq.append(_to_commented_map(x))
        else:
            seq.append(x)

    return seq


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
    yaml_rt.width = 72  # keep reasonable line width to avoid large diffs
    yaml_rt.indent(mapping=2, sequence=4, offset=2)
    # Represent Python None as explicit 'null' in the YAML output
    try:
        yaml_rt.representer.add_representer(type(None),
                                           lambda dumper, value: dumper.represent_scalar('tag:yaml.org,2002:null', 'null'))
    except Exception:
        # Fallback: some ruamel versions expose representer differently; ignore if not available
        pass

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
                old.get("inputs"), gen_stage.get("inputs")
            )
            merged_out = merge_outputs(
                old.get("outputs"), gen_stage.get("outputs")
            )

            # Only replace inputs and outputs; preserve kind, command, etc.
            stages[stage_key]["inputs"] = merged_in
            stages[stage_key]["outputs"] = merged_out

        else:
            print(f"Inserting {stage_key}")
            after = find_insert_after(stages, stage_key)
            # Convert generated stage (plain dict) into CommentedMap
            new_value = CommentedMap()
            for k, v in gen_stage.items():
                if k in ("inputs", "outputs"):
                    new_value[k] = _to_commented_seq(v or [])
                else:
                    new_value[k] = v
            new_stages = insert_stage(stages, after, stage_key, new_value)
            calkit["pipeline"]["stages"] = new_stages
            stages = calkit["pipeline"]["stages"]

    # ------------------------------------------------------------------
    # Write back
    # ------------------------------------------------------------------

    # Normalize storage values: ensure empty storage entries are explicit nulls
    def _normalize_storage(node):
        # Handle mappings
        if isinstance(node, dict):
            # If this mapping contains a storage key with empty string, set to None
            if "storage" in node and (node.get("storage") == "" or node.get("storage") is None):
                node["storage"] = None
            for k, v in list(node.items()):
                _normalize_storage(v)
        # Handle sequences
        elif isinstance(node, list):
            for i, v in enumerate(node):
                # If sequence contains a mapping with a storage key
                if isinstance(v, dict) and "storage" in v and (v.get("storage") == "" or v.get("storage") is None):
                    v["storage"] = None
                _normalize_storage(v)

    _normalize_storage(calkit)

    print(f"\nWriting {calkit_path} ...")

    with open(calkit_path, "w") as fh:
        yaml_rt.dump(calkit, fh)

    print("Done.")


if __name__ == "__main__":
    main()
