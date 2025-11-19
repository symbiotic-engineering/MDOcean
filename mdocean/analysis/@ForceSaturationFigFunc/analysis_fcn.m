function intermed_result_struct = analysis_fcn(p,b)
    % Run force saturation analysis
    [fig1, fig2] = force_sat_results(p, b);
    
    % Store figures for post-processing
    intermed_result_struct.fig1 = fig1;
    intermed_result_struct.fig2 = fig2;
end
