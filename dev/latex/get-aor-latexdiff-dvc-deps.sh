#!/usr/bin/env bash
set -euo pipefail

yq -r '.stages["build-AOR-paper"].deps[]' dvc.yaml |
  while read -r dep; do
    git status --porcelain --ignored "$dep" | # get dependencies that are not git-tracked
    awk '$1 == "??" || $1 == "!!" {print $2}' |
    while read -r dep; do
      md5=$(yq -r --arg dep "$dep" '
        ([.stages[].outs[]? | select(.path == $dep)] | length) as $found
          | ($found | debug) as $_
          | if $found > 0 then
              # if this dep is a pipeline output, use the md5 from build-AOR-paper deps
              (.stages["build-AOR-paper"].deps[]? | select(.path == $dep) | .md5)
              //
              # if this dep is a pipeline output but not found in build-AOR-paper deps, ie because it is a subfolder, use the md5 from the first other pipeline stage
              ([.stages[].outs[]? | select(.path == $dep) | .md5] | .[0])
              // ""
            else
              # if this dep is not a pipeline output, return an empty string
              ""
            end
        ' dvc.lock )

      # if not a pipeline output; check for a standalone .dvc sidecar file
      if [[ -z "$md5" || "$md5" == "null" ]]; then
        dvc_file="${dep%/}.dvc"
        if [[ -f "$dvc_file" ]]; then
          md5=$(yq -r '.outs[0].md5' "$dvc_file")
        fi
      fi

      # pipeline output: use md5 from dvc.lock
      if [[ -n "$md5" && "$md5" != "null" ]]; then
        printf '%s\t%s\n' "$dep" "$md5"
      fi
    done
  done |
  jq -Rn '[inputs | split("\t") | {(.[0]): .[1]}] | add' \
    > aor-latexdiff-dvc-deps.json
