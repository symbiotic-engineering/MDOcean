#!/bin/bash

# Path to the target MATLAB script
MATLAB_FILE="/opt/hostedtoolcache/MATLAB/2024.2.999/x64/toolbox/globaloptim/+globaloptim/+internal/@FcnEvaluator/evaluate_constr_obj_on_feasible_parallel.m"

# Check if the file exists
if [ ! -f "$MATLAB_FILE" ]; then
  echo "File does not exist: $MATLAB_FILE"
  exit 1
fi

# Print current contents of the file
echo "Current contents of the file:"
cat "$MATLAB_FILE"

# The lines of code to insert
CODE_TO_INSERT="
if isRowOriented
    nonlinCineq = nonlinCineq.';
    nonlinCeq = nonlinCeq.';
end
"

# Insert the code at line 17
sed -i "17i$CODE_TO_INSERT" "$MATLAB_FILE"

# Print new contents of the file to verify the changes
echo "New contents of the file:"
cat "$MATLAB_FILE"
