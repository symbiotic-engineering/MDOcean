#!/bin/sh
# take stdin, write it to a temp .tex file, run de-macro, delete the temp file, return as stdout
set -eu

tmp=$(mktemp --tmpdir=. .de-macro-filter-XXXXXX.tex)
clean="${tmp%.tex}-clean.tex"

trap 'rm -f "$tmp" "$clean"' EXIT

cat > "$tmp"
de-macro "$tmp"
cat "$clean"