#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 2 ]]; then
  echo "usage: $0 <old-hash> <new-hash>" >&2
  exit 1
fi

old_hash="$1"
new_hash="$2"

# build a {path: md5} json object for the DVC-managed deps of build-AOR-paper
# as recorded at a given git revision
deps_json_for_hash() {
  local hash="$1"

  yq -r '.stages["build-AOR-paper"].deps[]' <(git show "$hash:dvc.yaml") |
    while read -r dep; do
      # skip deps that are tracked by git (ie not DVC-managed) at this revision
      if git ls-tree -r --name-only "$hash" -- "$dep" | grep -q .; then
        continue
      fi

      md5=$(yq -r --arg dep "$dep" '
          # a dep may be a plain DVC-tracked path (md5 recorded directly on
          # build-AOR-paper itself) or a pipeline output produced by another stage
          (.stages["build-AOR-paper"].deps[]? | select(.path == $dep) | .md5)
          //
          ([.stages[].outs[]? | select(.path == $dep) | .md5] | .[0])
          // ""
        ' <(git show "$hash:dvc.lock"))

      # if not a pipeline output; check for a standalone .dvc sidecar file
      if [[ -z "$md5" || "$md5" == "null" ]]; then
        dvc_file="${dep%/}.dvc"
        if git cat-file -e "$hash:$dvc_file" 2>/dev/null; then
          md5=$(git show "$hash:$dvc_file" | yq -r '.outs[0].md5')
        fi
      fi

      # pipeline output: use md5 from dvc.lock
      if [[ -n "$md5" && "$md5" != "null" ]]; then
        printf '%s\t%s\n' "$dep" "$md5"
      fi
    done |
    jq -Rn '[inputs | split("\t") | {(.[0]): .[1]}] | add // {}'
}

old_json=$(deps_json_for_hash "$old_hash")
new_json=$(deps_json_for_hash "$new_hash")

jq -n --argjson old "$old_json" --argjson new "$new_json" '{old: $old, new: $new}' \
  > aor-latexdiff-dvc-deps.json
