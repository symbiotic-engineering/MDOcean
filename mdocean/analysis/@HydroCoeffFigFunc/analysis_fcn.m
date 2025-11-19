function intermed_result_struct = analysis_fcn(~,~)
    % Run hydrodynamic coefficient error analysis
    [~, ~, fig] = hydro_coeff_err();
    
    % Store figure for post-processing
    intermed_result_struct.created_figs = fig;
end
