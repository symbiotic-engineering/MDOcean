#!/bin/bash

PAPER_DIR="paper"
FIGS_DIR="$PAPER_DIR/figs"

# Your macro expansion pattern: replace #1 with the macro argument
MATLAB_PATTERN="figs/2025-04-11_10.21.55/Figure_#1.pdf"

echo "Scanning .tex files under $PAPER_DIR for included graphics..."

# Extract the argument inside { } for all \includegraphics
mapfile -t included_args < <(
  grep -r --include="*.tex" '\\includegraphics' "$PAPER_DIR" | \
  perl -nle '
    while (/\\includegraphics(?:\[[^\]]*\])?{((?:[^{}]+|{[^{}]*})*)}/g) {
      print $1;
    }
  ' | sort -u
)

echo "Found ${#included_args[@]} unique includegraphics arguments."

declare -A used_figs

echo "Expanding macros and recording used figure paths..."

for arg in "${included_args[@]}"; do
  arg="$(echo "$arg" | xargs)"  # Trim whitespace

  if [[ "$arg" =~ matlabFilepath\{([^}]+)\} ]]; then
    num="${BASH_REMATCH[1]}"
    expanded="${MATLAB_PATTERN//#1/$num}"
    expanded_lower=$(echo "$expanded" | tr '[:upper:]' '[:lower:]')
    used_figs["$expanded_lower"]=1

  elif [[ "$arg" == figs/* ]]; then
    arg_lower=$(echo "$arg" | tr '[:upper:]' '[:lower:]')
    used_figs["$arg_lower"]=1

  else
    arg_lower=$(echo "$arg" | tr '[:upper:]' '[:lower:]')
    used_figs["$arg_lower"]=1
  fi
done

echo "After macro expansion, ${#used_figs[@]} figures marked as used."
echo "Checking files in $FIGS_DIR ..."

find "$FIGS_DIR" -type f -print0 | while IFS= read -r -d '' figfile; do
  relpath="${figfile#$PAPER_DIR/}"
  relpath_lower=$(echo "$relpath" | tr '[:upper:]' '[:lower:]')

  if [[ -n "${used_figs[$relpath_lower]}" ]]; then
    # Used — keep
    :
  else
    echo "Deleting unused figure: $figfile"
    rm "$figfile"
  fi
done

echo "Done."
