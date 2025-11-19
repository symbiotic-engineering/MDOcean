function intermed_result_struct = analysis_fcn(p,b)
    % Run single analysis with drag convergence plotting
    p.operational_dynamics_debug_plots_on = true;
    figs = run_single(p,b);

    % Store figure numbers for post-processing
    intermed_result_struct.created_figs = figs;
end
