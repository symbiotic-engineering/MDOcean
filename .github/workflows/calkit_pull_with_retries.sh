#!/usr/bin/env bash
set -u

max_tries=10
try=1

while (( try <= max_tries )); do
  echo "=== Attempt $try / $max_tries ==="

  git rev-parse HEAD
  git submodule status --recursive
  git ls-files --stage | grep 160000 || true
  git config --get-regexp submodule || true

  # Run pull and capture ALL output (stdout+stderr)
  out="$(calkit pull 2>&1)"
  echo "$out"

  # Check for the specific retry-able error pattern
  if echo "$out" | grep -Eq 'ERROR: failed to pull data from the cloud - [0-9]+ files failed to download|ERROR: unexpected error'; then
    echo "Retryable error detected. Retrying..."
    ((try++))
    sleep 2
  else
    echo "No matching error found. Done."
    exit 0
  fi
done

echo "ERROR: exceeded maximum retries ($max_tries)."
exit 1