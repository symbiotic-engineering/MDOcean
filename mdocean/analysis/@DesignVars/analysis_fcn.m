function intermed_result_struct = analysis_fcn(~,b)
    % Get design variables information
    % Store results for post-processing
    intermed_result_struct.X_mins = b.X_mins;
    intermed_result_struct.X_noms = b.X_noms;
    intermed_result_struct.X_maxs = b.X_maxs;
    intermed_result_struct.var_names = b.var_names;
end
