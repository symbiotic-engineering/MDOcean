toolboxes = {ver().Name};

req_toolboxes = {'Global Optimization Toolbox', 'Statistics and Machine Learning Toolbox','Optimization Toolbox'};
req_for_test_toolboxes = {'Simscape','Simscape Multibody','Simulink'};
if ~isMATLABReleaseOlderThan('R2023a')
    req_for_test_toolboxes(end+1) = {'Report Generation Toolbox'};
end
optional_toolboxes = {'Symbolic Math Toolbox','Parallel Computing Toolbox'};

req = contains(req_toolboxes, toolboxes);
assert(all(req),['Missing required toolboxes: ' req_toolboxes{~req}])

test = contains(req_for_test_toolboxes, toolboxes);
assert(all(test),['Missing required for testing toolboxes: ' req_for_test_toolboxes{~test}])

opt = contains(optional_toolboxes, toolboxes);
if ~all(opt)
    warning(['Missing optional toolboxes: ' req_toolboxes{~opt}])
end