#!/usr/bin/env python3

import sys
import yaml
import hashlib
import subprocess
from copy import deepcopy
from pathlib import Path

BASE_LOCK, OURS_LOCK, THEIRS_LOCK = sys.argv[1:4]
REPO_ROOT = Path.cwd()


# -------------------------------
# Utility
# -------------------------------

def load_yaml(path):
    with open(path) as f:
        return yaml.safe_load(f) or {}


def dump_yaml(data, path):
    with open(path, "w") as f:
        yaml.safe_dump(data, f, sort_keys=False)


def repo_file_modified(filepath):
    """
    Returns True if git sees changes/conflicts in the given file.
    Conservative: abort merge if modified.
    """
    try:
        result = subprocess.run(
            ["git", "status", "--porcelain", filepath],
            capture_output=True,
            text=True,
            check=False,
        )
        return bool(result.stdout.strip())
    except Exception:
        return True  # safest failure mode


def check_pipeline_definition_safety():
    """
    Abort merge if dvc.yaml is modified or conflicted.
    """
    if repo_file_modified("dvc.yaml"):
        print("Pipeline definition file dvc.yaml is modified or conflicted.")
        print("Resolve dvc.yaml merge first, then retry.")
        sys.exit(1)


# -------------------------------
# Stage-aware merge logic
# -------------------------------

def merge_lockfiles(base, ours, theirs):
    base_stages = base.get("stages", {})
    our_stages = ours.get("stages", {})
    their_stages = theirs.get("stages", {})

    merged_stages = {}
    conflict_stage_names = []
    conflict_blocks = []

    all_stages = set(base_stages) | set(our_stages) | set(their_stages)

    for stage in sorted(all_stages):
        b = base_stages.get(stage)
        o = our_stages.get(stage)
        t = their_stages.get(stage)

        # Identical
        if o == t:
            merged_stages[stage] = o
            print(f"[MERGED] {stage} (identical)")
            continue

        # Only ours changed
        if b == t and o != b:
            merged_stages[stage] = o
            print(f"[MERGED] {stage} (chose OURS)")
            continue

        # Only theirs changed
        if b == o and t != b:
            merged_stages[stage] = t
            print(f"[MERGED] {stage} (chose THEIRS)")
            continue

        # Newly added stage
        if b is None:
            if o and not t:
                merged_stages[stage] = o
                print(f"[MERGED] {stage} (new, OURS)")
                continue
            if t and not o:
                merged_stages[stage] = t
                print(f"[MERGED] {stage} (new, THEIRS)")
                continue

        # TRUE CONFLICT
        print(f"[CONFLICT] {stage}")
        conflict_stage_names.append(stage)

        ours_yaml = yaml.safe_dump(o, sort_keys=False)
        theirs_yaml = yaml.safe_dump(t, sort_keys=False)

        block = (
            f"{stage}:\n"
            f"<<<<<<< OURS\n"
            f"{ours_yaml}"
            f"=======\n"
            f"{theirs_yaml}"
            f">>>>>>> THEIRS\n"
        )

        conflict_blocks.append(block)

    return merged_stages, conflict_stage_names, conflict_blocks


# -------------------------------
# Main
# -------------------------------

def main():
    def has_conflict_markers(path):
        with open(path, "r", encoding="utf-8") as f:
            text = f.read()
        return "<<<<<<<" in text or "=======" in text or ">>>>>>>" in text
    
    for name, path in [
        ("BASE", BASE_LOCK),
        ("OURS", OURS_LOCK),
        ("THEIRS", THEIRS_LOCK),
    ]:
        if has_conflict_markers(path):
            print(f"{name} file already contains conflict markers. Aborting.")
            sys.exit(0)
    
    check_pipeline_definition_safety()

    base = load_yaml(BASE_LOCK)
    ours = load_yaml(OURS_LOCK)
    theirs = load_yaml(THEIRS_LOCK)

    merged_stages, conflict_stage_names, conflict_blocks = merge_lockfiles(
        base, ours, theirs
    )
    
    with open(OURS_LOCK, "w") as f:
        f.write("stages:\n")
    
        # Write clean stages
        for stage, data in merged_stages.items():
            yaml_text = yaml.safe_dump({stage: data}, sort_keys=False)
            indented = "\n".join("  " + line for line in yaml_text.splitlines())
            f.write(indented + "\n")
    
        # Write conflict stages
        for block in conflict_blocks:
            indented = "\n".join("  " + line for line in block.splitlines())
            f.write(indented + "\n")
    
    if conflict_stage_names:
        print("\nConflicts detected in these stages:")
        for stage in conflict_stage_names:
            print(f"  - {stage}")
    
        print("\nResolve these conflicts manually in dvc.lock.")
        sys.exit(1)
    # ---------------------------------------
    # No conflicts: safe to run global dry repro
    # ---------------------------------------

    print("Merged dvc.lock successfully.")
    print("Running `dvc repro --dry` to detect stages that need recomputation...\n")

    result = subprocess.run(["dvc", "repro", "--dry"])

    if result.returncode != 0:
        print("\nSome stages need to be recomputed. Run `dvc repro` to update them.")
    else:
        print("\nAll stages are up-to-date after merge.")

    sys.exit(0)


if __name__ == "__main__":
    main()
