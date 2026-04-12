#!/usr/bin/env python3

import sys
import yaml
import hashlib
import subprocess
from copy import deepcopy
from pathlib import Path

BASE_LOCK, OURS_LOCK, THEIRS_LOCK = sys.argv[1:4]
REPO_ROOT = Path.cwd()
CALKIT_YAML = REPO_ROOT / "calkit.yaml"


# -------------------------------
# Utility
# -------------------------------

def load_yaml(path):
    with open(path) as f:
        return yaml.safe_load(f) or {}


def dump_yaml(data, path):
    with open(path, "w") as f:
        yaml.safe_dump(data, f, sort_keys=False)


# Helper: normalize path strings for comparison
def _norm_path(p):
    if not isinstance(p, str):
        return ""
    if p.startswith("./"):
        return p[2:]
    return p


# Helper: extract 'path' entries from a stage's deps list
def _extract_dep_paths(stage_dict):
    if not stage_dict:
        return []
    deps = stage_dict.get("deps") or []
    paths = []
    for dep in deps:
        if isinstance(dep, dict):
            p = dep.get("path")
        else:
            p = dep
        if p:
            paths.append(_norm_path(p))
    return paths


# Helper: discover .dvc files in repo and normalize names (without .dvc suffix)
def _repo_dot_dvc_paths():
    try:
        result = subprocess.check_output(["find", ".", "-name", "*.dvc", "-print"], text=True)
        lines = [line.strip() for line in result.splitlines() if line.strip()]
    except Exception:
        lines = []

    normalized = set()
    for p in lines:
        if p.startswith("./"):
            p = p[2:]
        if p.endswith(".dvc"):
            p = p[:-4]
        normalized.add(p)
    return normalized


# -------------------------------
# calkit.yaml helpers
# -------------------------------

def get_git_tracked_outputs(calkit_path=CALKIT_YAML):
    """Return a list of output paths from calkit.yaml stages that have storage: git."""
    try:
        with open(calkit_path) as f:
            calkit = yaml.safe_load(f) or {}
    except FileNotFoundError:
        print(f"Warning: {calkit_path} not found; no git-tracked outputs to check out.")
        return []

    paths = []
    stages = (calkit.get("pipeline") or {}).get("stages", {})
    for stage_data in stages.values():
        if not isinstance(stage_data, dict):
            continue
        for out in stage_data.get("outputs", []):
            if isinstance(out, dict) and out.get("storage") == "git":
                p = out.get("path")
                if p:
                    paths.append(_norm_path(p))
    return paths


def checkout_git_tracked_outputs(prefer):
    """Write ours/theirs blob content for each calkit.yaml output with storage: git.

    prefer must be "O" (ours) or "T" (theirs).  Any other value is a no-op.

    Uses `git cat-file blob :2:<path>` (ours) or `:3:<path>` (theirs) to read
    content directly from git's object store without touching the index.  This
    avoids the index.lock conflict that occurs when a git merge driver tries to
    run `git checkout --ours/--theirs` while git already holds the lock.

    Because the index is not updated here, the caller should print a reminder
    for the user to `git add` the resolved files.
    """
    if prefer not in {"O", "T"}:
        return

    # git index stages: 1=base, 2=ours, 3=theirs
    stage = "2" if prefer == "O" else "3"
    strategy = "ours" if prefer == "O" else "theirs"
    paths = get_git_tracked_outputs()

    if not paths:
        return

    print(f"\nWriting {strategy} version of git-tracked stage outputs:")
    resolved = []
    for p in paths:
        result = subprocess.run(
            ["git", "cat-file", "blob", f":{stage}:{p}"],
            capture_output=True,
        )
        if result.returncode == 0:
            try:
                dest = Path(p)
                dest.parent.mkdir(parents=True, exist_ok=True)
                dest.write_bytes(result.stdout)
                print(f"  [OK] {p}")
                resolved.append(p)
            except OSError as exc:
                print(f"  [ERROR] {p}: {exc}")
        else:
            # File not staged at this conflict stage (not in conflict, or
            # missing on one side) — nothing to do.
            msg = result.stderr.decode("utf-8", errors="replace").strip()
            print(f"  [SKIP] {p}: {msg}")

    if resolved:
        print(
            "\nThese files were written to the working tree but not staged.\n"
            "After the merge completes, run:\n"
            f"  git add {' '.join(resolved)}"
        )


# -------------------------------
# Stage-aware merge logic
# -------------------------------

def merge_lockfiles(base, ours, theirs):
    base_stages = base.get("stages", {}) or {}
    our_stages = ours.get("stages", {}) or {}
    their_stages = theirs.get("stages", {}) or {}

    merged_stages = {}
    conflict_stage_names = []
    conflict_blocks = []

    all_stages = set(base_stages) | set(our_stages) | set(their_stages)

    prefer = input("Prefer (O)URS, (T)HEIRS, or (C)ONFLICT for stages with non-.dvc conflicts? [O/T/C]: ").strip().upper()
    if prefer not in {"O", "T", "C"}:
        print("Invalid choice. Defaulting to CONFLICT.")
        prefer = "C"

    # precompute repo .dvc targets
    dot_dvc_targets = _repo_dot_dvc_paths()

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

        # TRUE CONFLICT: decide whether the conflict involves .dvc deps
        dep_paths = set()
        for st in (b, o, t):
            dep_paths.update(_extract_dep_paths(st))

        # If none of the deps reference a .dvc target in the repo, treat as non-.dvc conflict
        has_dot_dvc_dep = any(p in dot_dvc_targets for p in dep_paths)

        if not has_dot_dvc_dep:
            if prefer == "O":
                merged_stages[stage] = o
                print(f"[MERGED] {stage} (non-.dvc conflict, chose OURS)")
                continue
            elif prefer == "T":
                merged_stages[stage] = t
                print(f"[MERGED] {stage} (non-.dvc conflict, chose THEIRS)")
                continue

        # Real conflict that needs manual resolution
        print(f"[CONFLICT] {stage} (conflicting .dvc deps or explicit conflict)")
        print("Resolve this conflict manually in dvc.lock.")

        conflict_stage_names.append(stage)

        # Prepare YAML fragments for conflict block. We want the inner content
        # of the stage (no leading "stage:"), and we will emit these lines
        # indented appropriately when writing the final file.
        ours_yaml = yaml.safe_dump(o or {}, sort_keys=False).rstrip()
        theirs_yaml = yaml.safe_dump(t or {}, sort_keys=False).rstrip()

        # Build block where the first line is the stage key and following lines
        # are prefixed with two spaces. When the entire block is later prefixed
        # by two spaces (for top-level 'stages:' indentation), the inner lines
        # will be indented by four spaces (correct nesting).
        block_lines = [f"{stage}:"]
        block_lines.append("  <<<<<<< OURS")
        if ours_yaml:
            for line in ours_yaml.splitlines():
                block_lines.append("  " + line)
        block_lines.append("  =======")
        if theirs_yaml:
            for line in theirs_yaml.splitlines():
                block_lines.append("  " + line)
        block_lines.append("  >>>>>>> THEIRS")

        conflict_blocks.append("\n".join(block_lines))

    return merged_stages, conflict_stage_names, conflict_blocks, prefer


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

    base = load_yaml(BASE_LOCK)
    ours = load_yaml(OURS_LOCK)
    theirs = load_yaml(THEIRS_LOCK)

    merged_stages, conflict_stage_names, conflict_blocks, prefer = merge_lockfiles(
        base, ours, theirs
    )

    # If there are no conflicts, write a full merged lock preserving other metadata.
    if not conflict_stage_names:
        merged_lock = deepcopy(ours or {})
        merged_lock["stages"] = merged_stages
        dump_yaml(merged_lock, OURS_LOCK)
    else:
        # Preserve non-stage top-level metadata (prefer OURS for those keys),
        # then emit a stages: block with both merged and conflict entries.
        meta = {k: v for k, v in (ours or {}).items() if k != "stages"}
        # Write metadata first using safe_dump so ordering is standard YAML
        with open(OURS_LOCK, "w", encoding="utf-8") as f:
            if meta:
                f.write(yaml.safe_dump(meta, sort_keys=False))
            f.write("stages:\n")

            # Write clean stages (these are dicts -> emit YAML and indent two spaces)
            for stage, data in merged_stages.items():
                yaml_text = yaml.safe_dump({stage: data}, sort_keys=False)
                # indent every line by two spaces
                indented = "\n".join("  " + line for line in yaml_text.splitlines())
                f.write(indented + "\n")

            # Write conflict blocks (already prepared with stage header
            # and inner lines prefixed by two spaces). Prefix each line with two spaces.
            for block in conflict_blocks:
                indented = "\n".join("  " + line for line in block.splitlines())
                f.write(indented + "\n")

    if conflict_stage_names:
        print("\nConflicts detected in these stages:")
        for stage in conflict_stage_names:
            print(f"  - {stage}")

        print("\nResolve these conflicts manually in dvc.lock.")
        checkout_git_tracked_outputs(prefer)
        sys.exit(1)
    # ---------------------------------------
    # No conflicts: safe to run global dry repro
    # ---------------------------------------

    print("Merged dvc.lock successfully.")
    checkout_git_tracked_outputs(prefer)
    print("Running `dvc repro --dry` to detect stages that need recomputation...\n")

    result = subprocess.run(["dvc", "repro", "--dry"])

    if result.returncode != 0:
        print("\nSome stages need to be recomputed. Run `dvc repro` to update them.")
    else:
        print("\nAll stages are up-to-date after merge.")

    sys.exit(0)


if __name__ == "__main__":
    main()
