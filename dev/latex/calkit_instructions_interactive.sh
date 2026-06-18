#!/usr/bin/env bash

set -u

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "$REPO_ROOT"
SAFE_NAME_PATTERN='^[A-Za-z0-9._-]+$'

prompt_yes_no() {
  local message="$1"
  local answer
  while true; do
    read -r -p "$message [y/n]: " answer
    case "$answer" in
      y|Y) return 0 ;;
      n|N) return 1 ;;
      *) echo "Please answer y or n." ;;
    esac
  done
}

prompt_nonempty() {
  local message="$1"
  local value
  while true; do
    read -r -p "$message" value
    if [ -n "$value" ]; then
      printf '%s\n' "$value"
      return 0
    fi
    echo "Input cannot be empty."
  done
}

prompt_valid_branch_name() {
  local message="$1"
  local value
  while true; do
    value="$(prompt_nonempty "$message")"
    if git check-ref-format --branch "$value" >/dev/null 2>&1; then
      printf '%s\n' "$value"
      return 0
    fi
    echo "Invalid branch name. Use a valid git branch name."
  done
}

prompt_valid_folder_name() {
  local message="$1"
  local value
  while true; do
    value="$(prompt_nonempty "$message")"
    if [[ "$value" =~ $SAFE_NAME_PATTERN ]]; then
      printf '%s\n' "$value"
      return 0
    fi
    echo "Invalid folder name. Use letters, numbers, dot, underscore, or hyphen only (no slashes)."
  done
}

prompt_valid_filename() {
  local message="$1"
  local value
  while true; do
    value="$(prompt_nonempty "$message")"
    if [[ "$value" =~ $SAFE_NAME_PATTERN ]]; then
      printf '%s\n' "$value"
      return 0
    fi
    echo "Invalid filename. Use letters, numbers, dot, underscore, or hyphen only (no slashes)."
  done
}

run_step() {
  local description="$1"
  shift
  local -a cmd=( "$@" )
  local answer
  local command_preview

  while true; do
    echo
    echo "Step: $description"
    printf -v command_preview '%q ' "${cmd[@]}"
    echo "Command: ${command_preview% }"
    read -r -p "Press Enter to run, 's' to skip, or 'q' to quit: " answer

    case "$answer" in
      q|Q)
        echo "Stopped by user."
        exit 1
        ;;
      s|S)
        echo "Skipped."
        return 0
        ;;
      "")
        if "${cmd[@]}"; then
          return 0
        fi

        echo "Command failed."
        read -r -p "Enter to retry, 's' to skip, or 'q' to quit: " answer
        case "$answer" in
          q|Q)
            echo "Stopped by user."
            exit 1
            ;;
          s|S)
            echo "Skipped after failure."
            return 0
            ;;
          *)
            ;;
        esac
        ;;
      *)
        echo "Invalid choice."
        ;;
    esac
  done
}

run_calkit_pull_with_retry() {
  local -a pull_cmd=( "$@" )
  local answer
  local command_preview

  while true; do
    echo
    echo "Step: Pull data from Calkit"
    printf -v command_preview '%q ' "${pull_cmd[@]}"
    echo "Command: ${command_preview% }"
    read -r -p "Press Enter to run, 's' to skip, or 'q' to quit: " answer

    case "$answer" in
      q|Q)
        echo "Stopped by user."
        exit 1
        ;;
      s|S)
        echo "Skipped."
        return 0
        ;;
      "")
        if "${pull_cmd[@]}"; then
          return 0
        fi

        echo "calkit pull failed. The maintainer docs suggest retrying until success."
        if prompt_yes_no "Retry calkit pull now?"; then
          continue
        fi
        if prompt_yes_no "Skip this pull step?"; then
          return 0
        fi
        echo "Stopped by user."
        exit 1
        ;;
      *)
        echo "Invalid choice."
        ;;
    esac
  done
}

run_up_to_date_check() {
  local stage="$1"
  local output_file
  local answer

  output_file="$(mktemp)" || {
    echo "Failed to create temporary file for status check output."
    exit 1
  }

  while true; do
    echo
    echo "Check: Data and pipelines must be up to date"
    if calkit dvc status -d "$stage" | tee "$output_file"; then
      if grep -q "Data and pipelines are up to date." "$output_file"; then
        rm -f "$output_file"
        return 0
      fi

      echo "Expected status message not found in output."
      if prompt_yes_no "Continue anyway?"; then
        rm -f "$output_file"
        return 0
      fi
    else
      echo "Status command failed."
      if prompt_yes_no "Retry status check?"; then
        continue
      fi
      if prompt_yes_no "Skip this status check?"; then
        rm -f "$output_file"
        return 0
      fi
      echo "Stopped by user."
      rm -f "$output_file"
      exit 1
    fi

    read -r -p "Enter to retry, 's' to skip, or 'q' to quit: " answer
    case "$answer" in
      q|Q)
        echo "Stopped by user."
        rm -f "$output_file"
        exit 1
        ;;
      s|S)
        rm -f "$output_file"
        return 0
        ;;
      *)
        ;;
    esac
  done
}

stage_paper_tex_updates() {
  local tex_found=0
  local tex_file

  while IFS= read -r -d '' tex_file; do
    git add "$tex_file" || {
      echo "Failed to stage file: $tex_file. Check file permissions and git status."
      return 1
    }
    tex_found=1
  done < <(find "pubs/$PAPER_FOLDER" -maxdepth 1 -type f -name "*.tex" -print0)

  if [ "$tex_found" -eq 0 ]; then
    echo "No top-level .tex files found in pubs/$PAPER_FOLDER."
  fi

  if [ -d "pubs/$PAPER_FOLDER/sections" ]; then
    git add "pubs/$PAPER_FOLDER/sections/" || {
      echo "Failed to stage directory: pubs/$PAPER_FOLDER/sections/. Check directory permissions and git status."
      return 1
    }
  else
    echo "No sections directory in pubs/$PAPER_FOLDER; skipping sections staging."
  fi
}

choose_stage_and_folder() {
  echo
  echo "Select paper target:"
  echo "  1) build-AOR-paper / applied-ocean-research-model"
  echo "  2) build-RE-paper / renewable-energy-mdo"
  echo "  3) Custom"

  local choice
  while true; do
    read -r -p "Choice [1/2/3]: " choice
    case "$choice" in
      1)
        PAPER_STAGE="build-AOR-paper"
        PAPER_FOLDER="applied-ocean-research-model"
        return 0
        ;;
      2)
        PAPER_STAGE="build-RE-paper"
        PAPER_FOLDER="renewable-energy-mdo"
        return 0
        ;;
      3)
        PAPER_STAGE="$(prompt_nonempty 'Enter <paper-stage>: ')"
        PAPER_FOLDER="$(prompt_valid_folder_name 'Enter <paper-folder> (no slashes): ')"
        return 0
        ;;
      *)
        echo "Please enter 1, 2, or 3."
        ;;
    esac
  done
}

workflow_github_to_overleaf() {
  local branch_name

  choose_stage_and_folder
  branch_name="$(prompt_valid_branch_name 'Enter branch name to create (e.g., overleaf-my-update): ')"

  run_step "Checkout main" git checkout main
  run_step "Pull latest main" git pull
  run_step "Create new branch" git checkout -b "$branch_name"
  run_step "Push branch and set upstream" git push --set-upstream origin "$branch_name"
  run_calkit_pull_with_retry calkit pull -f
  run_step "Check pipeline" calkit check pipeline -c
  run_up_to_date_check "$PAPER_STAGE"
  run_step "Sync GitHub changes to Overleaf" calkit overleaf sync --push-only "pubs/$PAPER_FOLDER"

  echo
  echo "Done: GitHub changes were synced toward Overleaf flow."
}

workflow_overleaf_to_github() {
  local branch_name
  local shared_file

  choose_stage_and_folder

  run_step "Checkout main" git checkout main
  run_calkit_pull_with_retry calkit pull -f
  run_step "Check pipeline" calkit check pipeline -c
  run_up_to_date_check "$PAPER_STAGE"
  run_step "Show Overleaf status" calkit overleaf status "pubs/$PAPER_FOLDER"

  if ! prompt_yes_no "Do the listed Overleaf changes look correct?"; then
    echo "Stopped so you can review differences."
    exit 1
  fi

  run_step "Sync Overleaf changes locally without commit" calkit overleaf sync --no-commit "pubs/$PAPER_FOLDER"

  echo
  echo "If merge conflicts occurred, resolve them now before continuing."
  if ! prompt_yes_no "Have conflicts been resolved (or no conflicts occurred)?"; then
    echo "Stopped by user."
    exit 1
  fi

  branch_name="$(prompt_valid_branch_name 'Enter branch to work on (new or existing): ')"
  if prompt_yes_no "Create this branch with -b now?"; then
    run_step "Create and checkout branch" git checkout -b "$branch_name"
    run_step "Push branch with upstream" git push -u origin HEAD
  else
    run_step "Checkout existing branch" git checkout "$branch_name"
  fi

  # These files are shared outputs mapped across papers in the maintainer workflow docs.
  if prompt_yes_no "Did Overleaf change references.bib/shared-pkg.tex/elsarticle-num-names.bst?"; then
    while true; do
      shared_file="$(prompt_valid_filename 'Enter one changed filename (e.g., references.bib): ')"
      if [ ! -f "pubs/$PAPER_FOLDER/$shared_file" ]; then
        echo "File not found: pubs/$PAPER_FOLDER/$shared_file"
        if prompt_yes_no "Try a different filename?"; then
          continue
        fi
        echo "Stopped by user."
        exit 1
      fi
      run_step "Copy shared file" cp "pubs/$PAPER_FOLDER/$shared_file" "pubs/shared/$shared_file"
      run_step "Stage shared file" git add "pubs/shared/$shared_file"
      run_step "Commit shared file update" git commit -m "update shared $shared_file from overleaf"
      if ! prompt_yes_no "Add another shared file update?"; then
        break
      fi
    done
  fi

  if prompt_yes_no "Did you change figure order in Overleaf?"; then
    echo
    echo "Please update mdocean/plots/fig_tab_pub_mapping.m before continuing."
    if ! prompt_yes_no "Continue after updating that file?"; then
      echo "Stopped by user."
      exit 1
    fi
    run_step "Stage figure mapping" git add mdocean/plots/fig_tab_pub_mapping.m
    run_step "Commit figure order update" git commit -m "update figure order"

    if prompt_yes_no "Run figure-order pipeline steps locally now?"; then
      run_step "Generate calkit stages" calkit run make-calkit-stages
      run_step "Update calkit config" python mdocean/analysis/update_calkit.py
      run_step "Run paper stage" calkit run --log "$PAPER_STAGE"
      run_step "Save pipeline outputs" calkit save dvc.lock dvc.yaml .calkit/ calkit.yaml calkit_stages.yaml "pubs/$PAPER_FOLDER/numeric-results.tex" results/**/end.json results/**/*.tex -m "Run pipeline with updated fig order"
    else
      echo "Skipping local pipeline run; remember to verify CI for your branch."
    fi
  fi

  run_step "Stage paper tex updates" stage_paper_tex_updates
  run_step "Commit paper sync" git commit -m "Update paper from Overleaf sync"
  run_step "Push branch" git push

  echo
  echo "Done: Overleaf-to-GitHub flow completed."
}

main() {
  if [[ ! -t 0 ]]; then
    echo "This script requires an interactive terminal."
    exit 1
  fi

  echo "Interactive Calkit maintainer helper"
  echo "Source: docs/for_maintainers.rst (Writing section)"
  echo
  echo "Choose workflow:"
  echo "  1) Getting GitHub changes into Overleaf"
  echo "  2) Making changes on Overleaf and getting them onto GitHub"
  echo "  q) Quit"

  local workflow
  while true; do
    read -r -p "Choice [1/2/q]: " workflow
    case "$workflow" in
      1)
        workflow_github_to_overleaf
        return 0
        ;;
      2)
        workflow_overleaf_to_github
        return 0
        ;;
      q|Q)
        echo "Stopped by user."
        return 0
        ;;
      *)
        echo "Please enter 1, 2, or q."
        ;;
    esac
  done
}

main "$@"
