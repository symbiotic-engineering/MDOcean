function intermed_result_struct = analysis_fcn(p,b)
    % Run runtime comparison analysis
    [f1, f2, f3] = module_runtime_compare(p,b);

    % Store results for post-processing
    intermed_result_struct.dynamics_runtime_fig = f1;
    intermed_result_struct.hydro_runtime_fig = f2;
    intermed_result_struct.sim_runtime_fig = f3;
end
