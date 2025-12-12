function intermed_result_struct = analysis_fcn(~,~)
    % Run validation analysis
    [~, ~, ~, ~, tab_report, fig_cost_vs_N_WEC] = validate_nominal_RM3('report');
    [~,~,~,~,tab_wecsim] = validate_nominal_RM3('wecsim');
    
    % Store results for post-processing
    intermed_result_struct.tab_report = tab_report;
    intermed_result_struct.tab_wecsim = tab_wecsim;
    intermed_result_struct.fig_cost_vs_N_WEC = fig_cost_vs_N_WEC;
end
