#!/usr/bin/env bash
set -e

repo_root="$(git rev-parse --show-toplevel)"
deps_file="$repo_root/aor-latexdiff-dvc-deps.json"

for tree in old new; do
    while IFS= read -r subfolder; do
        subfolder="${subfolder%/}"
        source="$repo_root/$subfolder"
        target="../$tree/$subfolder"

        if [[ ! -e "$target" ]]; then
            mkdir -p "$(dirname "$target")"
            ln -s "$source" "$target"
            echo "Linked $source to $target"
        fi
    done < <(python3 -c 'import json, sys; print("\n".join(json.load(open(sys.argv[1]))))' "$deps_file")
done