# This script iterates through recent runs of a specific workflow and searches for a string
workflow_id="calkit-run.yml"
grep_for="ld.so"
grep_for_2="Running stage 'analysis"
gh run list -w $workflow_id --limit 50 --json databaseId,headBranch,headSha,createdAt | jq -r '.[] | "\(.databaseId)|\(.headBranch)|\(.headSha[0:7])|\(.createdAt)"' | while IFS='|' read -r run_id branch commit timestamp; do
  log=$(gh run view $run_id --log)
  analysis_ran=false
  ld_so_errored=false
  echo "$log" | grep -q "$grep_for_2" && analysis_ran=true
  echo "$log" | grep -q "$grep_for" && ld_so_errored=true
  if [ "$analysis_ran" = true ] && [ "$ld_so_errored" = true ] ; then
    echo "ANALYSIS RAN WITH LD.SO ERROR: $run_id (Branch: $branch, Commit: $commit, Date: $timestamp)"
  elif [ "$analysis_ran" = true ] ; then
    echo "ANALYSIS RAN WITHOUT LD.SO ERROR: $run_id (Branch: $branch, Commit: $commit, Date: $timestamp)"
  #else
    #echo "DID NOT RUN ANALYSIS: $run_id (Branch: $branch, Commit: $commit, Date: $timestamp)"
  fi
done
# earliest found: STRING WAS FOUND IN RUN: 20153627609 (Branch: figures-class-2, Commit: 985a6b3, Date: 2025-12-12T01:53:52Z)