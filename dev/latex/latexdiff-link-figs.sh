#!/usr/bin/env bash
set -e

repo_root="$(git rev-parse --show-toplevel)"


for tree in old new; do
    for subfolder in from-matlab manual ; do
        figs="$repo_root/pubs/applied-ocean-research-model/figs/$subfolder"
        target="../$tree/pubs/applied-ocean-research-model/figs/$subfolder"

        if [[ ! -e "$target" ]]; then
            ln -s "$figs" "$target"
        fi
    done
done