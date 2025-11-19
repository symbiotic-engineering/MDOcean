function intermed_result_struct = analysis_fcn(~,~)
    % Run WEC-Sim dynamics validation
    [~, ~, ~, tab, fig_singlebody, fig_multibody] = validate_dynamics();
    
    % Store results for post-processing
    intermed_result_struct.validation_table = tab;
    intermed_result_struct.fig_singlebody = fig_singlebody;
    intermed_result_struct.fig_multibody = fig_multibody;
    intermed_result_struct.all_sea_states_fig = gcf();
end
