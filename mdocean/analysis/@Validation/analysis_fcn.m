function intermed_result_struct = analysis_fcn(~,~)
    % Run validation analysis
    [~, ~, ~, ~, tab1a, fig_cost_vs_N_WEC] = validate_nominal_RM3('report');
    [~,~,~,~,tab1b] = validate_nominal_RM3('wecsim');
    
    % Store results for post-processing
    intermed_result_struct.tab1a = tab1a;
    intermed_result_struct.tab1b = tab1b;
    intermed_result_struct.fig_cost_vs_N_WEC = fig_cost_vs_N_WEC;
end
