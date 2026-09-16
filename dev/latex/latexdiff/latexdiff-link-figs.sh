#!/usr/bin/env bash
set -e

repo_root="${DVC_REPO_ROOT:-$(git rev-parse --show-toplevel)}"
cache_dir="${DVC_CACHE_DIR:-$repo_root/.dvc/cache}"
deps_file="${DVC_DEPS_JSON:-$repo_root/dev/latex/latexdiff/aor-latexdiff-dvc-deps.json}"

# resolve a DVC cache object path for a given file/directory-manifest md5
resolve_cache_object() {
    local md5="$1"
    local candidate
    for candidate in \
        "$cache_dir/${md5:0:2}/${md5:2}" \
        "$cache_dir/files/md5/${md5:0:2}/${md5:2}" \
        "$cache_dir/files/md5/${md5:0:2}/${md5:2}.dir"
    do
        if [[ -e "$candidate" ]]; then
            printf '%s\n' "$candidate"
            return 0
        fi
    done
    return 1
}

# materialize a DVC-tracked directory under $target from its .dir manifest
link_directory() {
    local md5="${1%.dir}"
    local target="$2"
    local manifest
    manifest="$(resolve_cache_object "$md5")" || return 1

    mkdir -p "$target"
    while IFS=$'\t' read -r relpath file_md5; do
        local file_source
        file_source="$(resolve_cache_object "$file_md5")" || continue
        mkdir -p "$(dirname "$target/$relpath")"
        ln -s "$file_source" "$target/$relpath"
    done < <(python3 -c '
import json, sys
with open(sys.argv[1]) as f:
    for entry in json.load(f):
        relpath = entry["relpath"]
        md5 = entry["md5"]
        print(relpath + "\t" + md5)
' "$manifest")
}

for tree in old new; do
    while IFS=$'\t' read -r subfolder md5; do
        subfolder="${subfolder%/}"
        target="../$tree/$subfolder"

        if [[ -e "$target" ]]; then
            continue
        fi

        if [[ "$md5" == *.dir ]]; then
            link_directory "$md5" "$target" || { echo "Skipping $subfolder ($tree, $md5): cache object not found" >&2; continue; }
        else
            source="$(resolve_cache_object "$md5")" || { echo "Skipping $subfolder ($tree, $md5): cache object not found" >&2; continue; }
            mkdir -p "$(dirname "$target")"
            ln -s "$source" "$target"
        fi

        echo "Linked $subfolder ($tree, $md5) to $target"
    done < <(python3 -c '
import json, sys
with open(sys.argv[1]) as f:
    deps = json.load(f).get(sys.argv[2], {})
for path, md5 in deps.items():
    print(f"{path}\t{md5}")
' "$deps_file" "$tree")
done