#!/usr/bin/env bash
set -euo pipefail

yq -r '.stages["build-AOR-paper"].deps[]' dvc.yaml |
  while read -r dep; do
    git status --porcelain --ignored "$dep" |
    awk '$1 == "??" || $1 == "!!" {print $2}' |
    while read -r dep; do
      md5=$(yq -r --arg dep "$dep" '
        .stages[].outs[]?
          | select(.path == $dep)
          | .md5
        ' dvc.lock | head -n1)

      # not a pipeline output; check for a standalone .dvc sidecar file
      if [[ -z "$md5" || "$md5" == "null" ]]; then
        dvc_file="${dep%/}.dvc"
        if [[ -f "$dvc_file" ]]; then
          md5=$(yq -r '.outs[0].md5' "$dvc_file")
        fi
      fi

      if [[ -n "$md5" && "$md5" != "null" ]]; then
        printf '%s\t%s\n' "$dep" "$md5"
      fi
    done
  done |
  jq -Rn '[inputs | split("\t") | {(.[0]): .[1]}] | add' \
    > aor-latexdiff-dvc-deps.json
