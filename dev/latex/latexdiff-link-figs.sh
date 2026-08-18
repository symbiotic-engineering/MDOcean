#!/usr/bin/env bash
set -e

repo_root="$(git rev-parse --show-toplevel)"
deps_file="$repo_root/aor-latexdiff-dvc-deps.txt"

for tree in old new; do
    while IFS= read -r subfolder; do
        subfolder="${subfolder%/}"
        source="$repo_root/$subfolder"
        target="../$tree/$subfolder"

        if [[ ! -e "$target" ]]; then
            mkdir -p "$(dirname "$target")"
            ln -s "$source" "$target"
        fi
    done < "$deps_file"
done