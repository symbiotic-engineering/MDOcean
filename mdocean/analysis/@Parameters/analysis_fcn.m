function intermed_result_struct = analysis_fcn(~,~)
    % Get parameters information
    [~, params_table] = parameters();
    
    % Store results for post-processing
    intermed_result_struct.params_table = params_table;
end
